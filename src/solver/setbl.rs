//! SETBL (xbl.f): sets up the BL Newton system coefficients for the current BL variables and
//! the edge velocities, and incorporates the local BL system coefficients into the global
//! Newton system (VA/VB/VM/VDEL/VZ) that BLSOLV consumes. Line-for-line on `SolverState`.
//!
//! Also MRCL (xfoil.f) — the Mach/Re dependence on CL selected by MATYP/RETYP — which SETBL
//! calls first. `IDAMPV = IDAMP` is pinned at 0 (the envelope e^N `DAMPL` model).

use crate::bl::blsolv::NewtonSystem;
use crate::bl::blsys::{assemble_interval_system, assemble_te_system, IntervalFlags};
use crate::bl::mrchdu::march_prescribed_dstar;
use crate::bl::mrchue::march_direct;
use crate::bl::system::{
    check_transition, FlowParameters, FlowRegime, IntervalSystem, MachClDependence, ReClDependence, TransitionCheck,
};
use crate::geometry::spline_value;
use crate::solver::blstate::SolverState;
use crate::solver::pointers::xi_trip;
use crate::solver::velocity::set_ue_with_sources;

/// Everything SETBL produces besides the arrays it writes into `SolverState`.
#[derive(Debug, Clone)]
pub struct AssembledSystem {
    /// The global BL Newton system as BLSOLV sees it (VA, VB, VDEL, VM, VZ, NSYS, IVTE1, IVZ).
    pub newton: NewtonSystem,
    /// The BL parameters SETBL derived (REYBL, HSTINV, ...): the marches and UPDATE share them.
    pub flow: FlowParameters,
    /// RE_CLMR, MSQ_CLMR: d(Re)/d(CL) and d(M²)/d(CL) for the fixed-CL sensitivity column.
    pub re_d_cl: f64,
    pub machsqd_d_cl: f64,
    /// MA_CLMR: d(M)/d(CL) (MRCL's M_CLS); VISCAL's MINF_CL is the same quantity
    pub mach_d_cl: f64,
    /// DULE1, DULE2: the LE Ue mismatch between UEDG and USAV = UINV + DIJ·MASS, per side.
    pub ue_le_mismatch: [f64; 3],
}

/// MRCL: sets the actual Mach and Reynolds numbers from the unit-CL values and the specified
/// CLS according to MATYP/RETYP. Returns (M_CLS, R_CLS).
pub fn set_mach_re_from_cl(st: &mut SolverState, cls: f64) -> (f64, f64) {
    let cla = cls.max(0.000001);
    // XFOIL's 'MRCL: Illegal Re(CL) / Mach(CL) dependence trigger. Setting fixed ...' branches
    // are unreachable here: the dependence enums cannot hold an illegal index.
    let mut m_cls;
    match st.mach_cl_dependence {
        MachClDependence::Fixed => {
            st.mach = st.mach_cl1;
            m_cls = 0.0;
        }
        MachClDependence::InverseSqrtCl => {
            st.mach = st.mach_cl1 / cla.sqrt();
            m_cls = -0.5 * st.mach / cla;
        }
    }
    let mut r_cls;
    match st.re_cl_dependence {
        ReClDependence::Fixed => {
            st.re = st.re_cl1;
            r_cls = 0.0;
        }
        ReClDependence::InverseSqrtCl => {
            st.re = st.re_cl1 / cla.sqrt();
            r_cls = -0.5 * st.re / cla;
        }
        ReClDependence::InverseCl => {
            st.re = st.re_cl1 / cla;
            r_cls = -st.re / cla;
        }
    }
    if st.mach >= 0.99 {
        // 'MRCL: CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99'
        st.mach = 0.99;
        m_cls = 0.0;
    }
    let mut rrat = 1.0;
    if st.re_cl1 > 0.0 {
        rrat = st.re / st.re_cl1;
    }
    if rrat > 100.0 {
        // 'MRCL: CL too low for chosen Re(CL) dependence — artificially limiting Re'
        st.re = st.re_cl1 * 100.0;
        r_cls = 0.0;
    }
    (m_cls, r_cls)
}

/// SETBL. Requires the pointer layer, DIJ, UINV/UINV_A, the BL arrays (initialised by MRCHUE
/// on the first call when `st.lblini` is false), SST_GO/SST_GP of the current stagnation point
/// and the control fields (`lalfa`, `cl`/`clspec`, `matyp`/`retyp`, `minf1`/`reinf1`, `acrit`,
/// `vaccel`). Writes TAU/DIS/CTQ/DELT/USLP, ITRAN/XSSITR/TFORCE, XOCTR/YOCTR/TINDEX and leaves
/// UEDG holding the *marched* Ue (USAV = UINV + DIJ·MASS is the mismatch reference).
pub fn assemble_newton_system(st: &mut SolverState) -> AssembledSystem {
    // set the CL used to define Mach, Reynolds numbers
    let clmr = if st.alpha_specified { st.cl } else { st.cl_specified };

    // set current MINF(CL)
    let (ma_clmr, re_clmr) = set_mach_re_from_cl(st, clmr);
    let msq_clmr = 2.0 * st.mach * ma_clmr;

    // set compressibility parameter TKLAM and derivative TK_MSQ (COMSET), gas constant, the
    // parameters for compressibility correction, stagnation density and 1/enthalpy, and the
    // Reynolds number based on freestream density, velocity, viscosity
    let mut params = FlowParameters::new(st.mach, st.re, st.gamma_gas);

    // IDAMPV = IDAMP
    params.amplification_model = st.amplification_model;
    // save TE thickness
    let _dwte = st.wake_gap[1];

    if !st.bl_initialised {
        // initialize BL by marching with Ue (fudge at separation)
        let acrit = st.ncrit;
        march_direct(st, &params, acrit, None);
        st.bl_initialised = true;
    }

    // march BL with current Ue and Ds to establish transition
    let acrit = st.ncrit;
    march_prescribed_dstar(st, &params, acrit, None);

    let mut usav: [Vec<f64>; 3] = [Vec::new(), st.ue[1].clone(), st.ue[2].clone()];
    set_ue_with_sources(st);
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            std::mem::swap(&mut usav[is][ibl], &mut st.ue[is][ibl]);
        }
    }

    let ile1 = st.i_node[1][2];
    let ile2 = st.i_node[2][2];
    let ite1 = st.i_node[1][st.i_te_station[1]];
    let ite2 = st.i_node[2][st.i_te_station[2]];

    let jvte1 = st.i_row[1][st.i_te_station[1]];
    let jvte2 = st.i_row[2][st.i_te_station[2]];

    let dule1 = st.ue[1][2] - usav[1][2];
    let dule2 = st.ue[2][2] - usav[2][2];

    let nsys = st.n_rows;

    // set LE and TE Ue sensitivities wrt all m values
    let mut ule1_m = vec![0.0; nsys + 1];
    let mut ule2_m = vec![0.0; nsys + 1];
    let mut ute1_m = vec![0.0; nsys + 1];
    let mut ute2_m = vec![0.0; nsys + 1];
    for js in 1..=2 {
        for jbl in 2..=st.n_stations[js] {
            let j = st.i_node[js][jbl];
            let jv = st.i_row[js][jbl];
            ule1_m[jv] = -st.velocity_sign[1][2] * st.velocity_sign[js][jbl] * st.dij[ile1][j];
            ule2_m[jv] = -st.velocity_sign[2][2] * st.velocity_sign[js][jbl] * st.dij[ile2][j];
            ute1_m[jv] = -st.velocity_sign[1][st.i_te_station[1]] * st.velocity_sign[js][jbl] * st.dij[ite1][j];
            ute2_m[jv] = -st.velocity_sign[2][st.i_te_station[2]] * st.velocity_sign[js][jbl] * st.dij[ite2][j];
        }
    }

    let ule1_a = st.ue_inviscid_d_alpha[1][2];
    let ule2_a = st.ue_inviscid_d_alpha[2][2];

    st.transition_node_fraction[1] = 0.0;
    st.transition_node_fraction[2] = 0.0;

    // the global system
    let mut va = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vb = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vdel = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vm = vec![vec![[0.0f64; 3]; nsys]; nsys];
    let mut vz = [[0.0f64; 2]; 3];

    // COM1/COM2 and the XT sensitivities are COMMON
    let mut s1 = std::mem::take(&mut st.station1);
    let mut s2 = std::mem::take(&mut st.station2);
    let mut trloc = std::mem::take(&mut st.transition);
    let mut sys = IntervalSystem::default();
    let (mut ami, mut cti) = (0.0, 0.0);
    let mut trforc = false;

    // Go over each boundary layer/wake
    for is in 1..=2 {
        // there is no station "1" at similarity, so zero everything out
        let mut u1_m = vec![0.0; nsys + 1];
        let mut d1_m = vec![0.0; nsys + 1];
        let mut u2_m = vec![0.0; nsys + 1];
        let mut d2_m = vec![0.0; nsys + 1];
        let mut u1_a = 0.0;
        let mut d1_a = 0.0;
        let mut due1 = 0.0;
        let mut dds1 = 0.0;

        // similarity station pressure gradient parameter  x/u du/dx: BULE = 1.0 (BLDIF's constant)
        let amcrit = st.ncrit[is];

        // set forced transition arc length position
        let xiforc = xi_trip(st, is);

        let mut tran;
        let mut turb;

        // Sweep downstream setting up BL equation linearizations
        for ibl in 2..=st.n_stations[is] {
            let iv = st.i_row[is][ibl];
            let simi = ibl == 2;
            let wake = ibl > st.i_te_station[is];
            tran = ibl == st.i_transition_station[is];
            turb = ibl > st.i_transition_station[is];
            let i = st.i_node[is][ibl];

            // set primary variables for current station
            let xsi = st.xi[is][ibl];
            if ibl < st.i_transition_station[is] {
                ami = st.sqrtctau[is][ibl];
            }
            if ibl >= st.i_transition_station[is] {
                cti = st.sqrtctau[is][ibl];
            }
            let uei = st.ue[is][ibl];
            let thi = st.theta[is][ibl];
            let mdi = st.mass_defect[is][ibl];
            let dsi = mdi / uei;
            let dswaki = if wake {
                st.wake_gap[ibl - st.i_te_station[is]]
            } else {
                0.0
            };

            // set derivatives of DSI (= D2)
            let d2_m2 = 1.0 / uei;
            let d2_u2 = -dsi / uei;
            for js in 1..=2 {
                for jbl in 2..=st.n_stations[js] {
                    let j = st.i_node[js][jbl];
                    let jv = st.i_row[js][jbl];
                    u2_m[jv] = -st.velocity_sign[is][ibl] * st.velocity_sign[js][jbl] * st.dij[i][j];
                    d2_m[jv] = d2_u2 * u2_m[jv];
                }
            }
            d2_m[iv] += d2_m2;

            let u2_a = st.ue_inviscid_d_alpha[is][ibl];
            let d2_a = d2_u2 * u2_a;

            // "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
            let due2 = st.ue[is][ibl] - usav[is][ibl];
            let dds2 = d2_u2 * due2;

            s2.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, &params);
            s2.set_kinematic_variables(&params);

            // check for transition and set TRAN, XT, etc. if found
            if tran {
                match check_transition(&s1, &s2, s1.ampl, amcrit, xiforc, &params) {
                    TransitionCheck::None { ampl2 } => {
                        ami = ampl2;
                        tran = false;
                        trloc.xi_transition = s2.xi;
                    }
                    TransitionCheck::Free {
                        transition: location,
                        ampl2,
                    } => {
                        ami = ampl2;
                        trforc = false;
                        trloc = location;
                    }
                    TransitionCheck::Forced { transition: location } => {
                        trforc = true;
                        trloc = location;
                    }
                }
                s2.ampl = ami;
            }
            // if ibl == st.itran[is] && !tran: 'SETBL: Xtr???  n1 n2: ' AMPL1, AMPL2

            // assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
            // at the previous "1" station and the current "2" station
            let (mut tte_tte1, mut tte_tte2) = (0.0, 0.0);
            let (mut cte_cte1, mut cte_cte2, mut cte_tte1, mut cte_tte2) = (0.0, 0.0, 0.0, 0.0);
            if ibl == st.i_te_station[is] + 1 {
                // define quantities at start of wake, adding TE base thickness to Dstar
                let (it1, it2) = (st.i_te_station[1], st.i_te_station[2]);
                let tte = st.theta[1][it1] + st.theta[2][it2];
                let dte = st.dstar[1][it1] + st.dstar[2][it2] + st.te_thickness_normal;
                let cte = (st.sqrtctau[1][it1] * st.theta[1][it1] + st.sqrtctau[2][it2] * st.theta[2][it2]) / tte;
                assemble_te_system(&mut sys, &mut s2, cte, tte, dte, &params);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                let dte_mte1 = 1.0 / st.ue[1][it1];
                let dte_ute1 = -st.dstar[1][it1] / st.ue[1][it1];
                let dte_mte2 = 1.0 / st.ue[2][it2];
                let dte_ute2 = -st.dstar[2][it2] / st.ue[2][it2];
                cte_cte1 = st.theta[1][it1] / tte;
                cte_cte2 = st.theta[2][it2] / tte;
                cte_tte1 = (st.sqrtctau[1][it1] - cte) / tte;
                cte_tte2 = (st.sqrtctau[2][it2] - cte) / tte;

                // re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
                for js in 1..=2 {
                    for jbl in 2..=st.n_stations[js] {
                        let jv = st.i_row[js][jbl];
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                    }
                }
                d1_m[jvte1] += dte_mte1;
                d1_m[jvte2] += dte_mte2;

                // "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
                due1 = 0.0;
                dds1 = dte_ute1 * (st.ue[1][it1] - usav[1][it1]) + dte_ute2 * (st.ue[2][it2] - usav[2][it2]);
            } else {
                let flags = IntervalFlags {
                    similarity: simi,
                    transition: tran,
                    turbulent: turb,
                    wake,
                };
                assemble_interval_system(&mut sys, &mut s1, &mut s2, flags, Some(&trloc), amcrit, &params);
            }

            // Save wall shear and equil. max shear coefficient for plotting output
            st.tau[is][ibl] = 0.5 * s2.rho * s2.ue * s2.ue * s2.cf;
            st.dissipation[is][ibl] = s2.rho * s2.ue * s2.ue * s2.ue * s2.cdiss * s2.hstar * 0.5;
            st.sqrtctaueq[is][ibl] = s2.sqrtctaueq;
            st.delta[is][ibl] = s2.delta;
            st.us_plot_scale[is][ibl] = 1.60 / (1.0 + s2.us);

            // set XI sensitivities wrt LE Ue changes
            let (xi_ule1, xi_ule2) = if is == 1 {
                (st.s_stagnation_d_gamma_node0, -st.s_stagnation_d_gamma_node1)
            } else {
                (-st.s_stagnation_d_gamma_node0, st.s_stagnation_d_gamma_node1)
            };

            // stuff BL system coefficients into main Jacobian matrix
            for k in 0..3 {
                let (vs1, vs2) = (&sys.jacobian_station1[k], &sys.jacobian_station2[k]);
                let vsx = sys.residual_d_xi[k];
                for jv in 1..=nsys {
                    vm[iv - 1][jv - 1][k] = vs1[2] * d1_m[jv]
                        + vs1[3] * u1_m[jv]
                        + vs2[2] * d2_m[jv]
                        + vs2[3] * u2_m[jv]
                        + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * ule1_m[jv] + xi_ule2 * ule2_m[jv]);
                }
                vb[iv - 1][k][0] = vs1[0];
                vb[iv - 1][k][1] = vs1[1];
                va[iv - 1][k][0] = vs2[0];
                va[iv - 1][k][1] = vs2[1];
                if st.alpha_specified {
                    vdel[iv - 1][k][1] = sys.residual_d_re[k] * re_clmr + sys.residual_d_machsqd[k] * msq_clmr;
                } else {
                    vdel[iv - 1][k][1] = (vs1[3] * u1_a + vs1[2] * d1_a)
                        + (vs2[3] * u2_a + vs2[2] * d2_a)
                        + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
                }
                vdel[iv - 1][k][0] = sys.residual[k]
                    + (vs1[3] * due1 + vs1[2] * dds1)
                    + (vs2[3] * due2 + vs2[2] * dds2)
                    + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * dule1 + xi_ule2 * dule2);
            }

            if ibl == st.i_te_station[is] + 1 {
                // redefine coefficients for TTE, DTE, etc
                for k in 0..3 {
                    let vs1 = &sys.jacobian_station1[k];
                    vz[k][0] = vs1[0] * cte_cte1;
                    vz[k][1] = vs1[0] * cte_tte1 + vs1[1] * tte_tte1;
                    vb[iv - 1][k][0] = vs1[0] * cte_cte2;
                    vb[iv - 1][k][1] = vs1[0] * cte_tte2 + vs1[1] * tte_tte2;
                }
            }

            // turbulent intervals will follow if currently at transition interval
            if tran {
                turb = true;

                // save transition location
                st.i_transition_station[is] = ibl;
                st.transition_forced[is] = trforc;
                st.xi_transition[is] = trloc.xi_transition;

                // interpolate airfoil geometry to find transition x/c (for user output)
                let str = if is == 1 {
                    st.s_stagnation - trloc.xi_transition
                } else {
                    st.s_stagnation + trloc.xi_transition
                };
                let chx = st.x_te - st.x_le;
                let chy = st.y_te - st.y_le;
                let chsq = chx * chx + chy * chy;
                let n = st.n_foil_nodes;
                let xtr = spline_value(str, &st.x[1..=n], &st.dxds[1..=n], &st.s[1..=n]);
                let ytr = spline_value(str, &st.y[1..=n], &st.dyds[1..=n], &st.s[1..=n]);
                st.x_transition[is] = ((xtr - st.x_le) * chx + (ytr - st.y_le) * chy) / chsq;
                st.y_transition[is] = ((ytr - st.y_le) * chx - (xtr - st.x_le) * chy) / chsq;
            }
            // (TRAN = .FALSE.; both flags are recomputed from ITRAN at the top of the sweep)
            let _ = turb;

            if ibl == st.i_te_station[is] {
                // set "2" variables at TE to wake correlations for next station
                // (BLVAR(3); BLMID(3) only sets the interval CFM, which the next BLSYS recomputes)
                s2.set_closure_variables(FlowRegime::Wake, &params);
            }

            for js in 1..=2 {
                for jbl in 2..=st.n_stations[js] {
                    let jv = st.i_row[js][jbl];
                    u1_m[jv] = u2_m[jv];
                    d1_m[jv] = d2_m[jv];
                }
            }
            u1_a = u2_a;
            d1_a = d2_a;
            due1 = due2;
            dds1 = dds2;

            if ibl == st.i_transition_station[is] && s2.xi > s1.xi {
                let frac = (trloc.xi_transition - s1.xi) / (s2.xi - s1.xi);
                st.transition_node_fraction[is] = if is == 1 {
                    (st.i_stagnation_node as i64 - st.i_transition_station[is] as i64 + 3) as f64 - frac
                } else {
                    (st.i_stagnation_node as i64 + st.i_transition_station[is] as i64 - 2) as f64 + frac
                };
            }

            // set BL variables for next station
            s1 = s2.clone();
        }
        // 'Side IS forced/free transition at x/c = XOCTR(IS) ITRAN(IS)'
    }

    st.station1 = s1;
    st.station2 = s2;
    st.transition = trloc;

    let sys = NewtonSystem {
        n_rows: nsys,
        diagonal: va,
        subdiagonal: vb,
        rhs: vdel,
        mass_influence: vm,
        te_block: vz,
        i_te_row_upper: Some(st.i_row[1][st.i_te_station[1]] - 1),
        i_wake_row: Some(st.i_row[2][st.i_te_station[2] + 1] - 1),
        elimination_threshold: st.elimination_threshold,
        s_total: Some(st.s[st.n_foil_nodes] - st.s[1]),
    };
    AssembledSystem {
        newton: sys,
        flow: params,
        re_d_cl: re_clmr,
        machsqd_d_cl: msq_clmr,
        mach_d_cl: ma_clmr,
        ue_le_mismatch: [0.0, dule1, dule2],
    }
}
