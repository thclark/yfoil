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
#[doc(alias = "MRCL")]
pub fn set_mach_re_from_cl(state: &mut SolverState, cls: f64) -> (f64, f64) {
    let cla = cls.max(0.000001);
    // XFOIL's 'MRCL: Illegal Re(CL) / Mach(CL) dependence trigger. Setting fixed ...' branches
    // are unreachable here: the dependence enums cannot hold an illegal index.
    let mut m_cls;
    match state.mach_cl_dependence {
        MachClDependence::Fixed => {
            state.mach = state.mach_cl1;
            m_cls = 0.0;
        }
        MachClDependence::InverseSqrtCl => {
            state.mach = state.mach_cl1 / cla.sqrt();
            m_cls = -0.5 * state.mach / cla;
        }
    }
    let mut r_cls;
    match state.re_cl_dependence {
        ReClDependence::Fixed => {
            state.re = state.re_cl1;
            r_cls = 0.0;
        }
        ReClDependence::InverseSqrtCl => {
            state.re = state.re_cl1 / cla.sqrt();
            r_cls = -0.5 * state.re / cla;
        }
        ReClDependence::InverseCl => {
            state.re = state.re_cl1 / cla;
            r_cls = -state.re / cla;
        }
    }
    if state.mach >= 0.99 {
        // 'MRCL: CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99'
        state.mach = 0.99;
        m_cls = 0.0;
    }
    let mut rrat = 1.0;
    if state.re_cl1 > 0.0 {
        rrat = state.re / state.re_cl1;
    }
    if rrat > 100.0 {
        // 'MRCL: CL too low for chosen Re(CL) dependence — artificially limiting Re'
        state.re = state.re_cl1 * 100.0;
        r_cls = 0.0;
    }
    (m_cls, r_cls)
}

/// SETBL. Requires the pointer layer, DIJ, UINV/UINV_A, the BL arrays (initialised by MRCHUE
/// on the first call when `st.lblini` is false), SST_GO/SST_GP of the current stagnation point
/// and the control fields (`lalfa`, `cl`/`clspec`, `matyp`/`retyp`, `minf1`/`reinf1`, `acrit`,
/// `vaccel`). Writes TAU/DIS/CTQ/DELT/USLP, ITRAN/XSSITR/TFORCE, XOCTR/YOCTR/TINDEX and leaves
/// UEDG holding the *marched* Ue (USAV = UINV + DIJ·MASS is the mismatch reference).
#[doc(alias = "SETBL")]
pub fn assemble_newton_system(state: &mut SolverState) -> AssembledSystem {
    // set the CL used to define Mach, Reynolds numbers
    let clmr = if state.alpha_specified {
        state.cl
    } else {
        state.cl_specified
    };

    // set current MINF(CL)
    let (ma_clmr, re_clmr) = set_mach_re_from_cl(state, clmr);
    let msq_clmr = 2.0 * state.mach * ma_clmr;

    // set compressibility parameter TKLAM and derivative TK_MSQ (COMSET), gas constant, the
    // parameters for compressibility correction, stagnation density and 1/enthalpy, and the
    // Reynolds number based on freestream density, velocity, viscosity
    let mut params = FlowParameters::new(state.mach, state.re, state.gamma_gas);

    // IDAMPV = IDAMP
    params.amplification_model = state.amplification_model;
    // save TE thickness
    let _dwte = state.wake_gap[1];

    if !state.bl_initialised {
        // initialize BL by marching with Ue (fudge at separation)
        let acrit = state.ncrit;
        march_direct(state, &params, acrit, None);
        state.bl_initialised = true;
    }

    // march BL with current Ue and Ds to establish transition
    let acrit = state.ncrit;
    march_prescribed_dstar(state, &params, acrit, None);

    let mut usav: [Vec<f64>; 3] = [Vec::new(), state.ue[1].clone(), state.ue[2].clone()];
    set_ue_with_sources(state);
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            std::mem::swap(&mut usav[side][i_station], &mut state.ue[side][i_station]);
        }
    }

    let ile1 = state.i_node[1][2];
    let ile2 = state.i_node[2][2];
    let ite1 = state.i_node[1][state.i_te_station[1]];
    let ite2 = state.i_node[2][state.i_te_station[2]];

    let jvte1 = state.i_row[1][state.i_te_station[1]];
    let jvte2 = state.i_row[2][state.i_te_station[2]];

    let dule1 = state.ue[1][2] - usav[1][2];
    let dule2 = state.ue[2][2] - usav[2][2];

    let nsys = state.n_rows;

    // set LE and TE Ue sensitivities wrt all m values
    let mut ule1_m = vec![0.0; nsys + 1];
    let mut ule2_m = vec![0.0; nsys + 1];
    let mut ute1_m = vec![0.0; nsys + 1];
    let mut ute2_m = vec![0.0; nsys + 1];
    for j_side in 1..=2 {
        for j_station in 2..=state.n_stations[j_side] {
            let j = state.i_node[j_side][j_station];
            let j_row = state.i_row[j_side][j_station];
            ule1_m[j_row] = -state.velocity_sign[1][2] * state.velocity_sign[j_side][j_station] * state.dij[ile1][j];
            ule2_m[j_row] = -state.velocity_sign[2][2] * state.velocity_sign[j_side][j_station] * state.dij[ile2][j];
            ute1_m[j_row] = -state.velocity_sign[1][state.i_te_station[1]]
                * state.velocity_sign[j_side][j_station]
                * state.dij[ite1][j];
            ute2_m[j_row] = -state.velocity_sign[2][state.i_te_station[2]]
                * state.velocity_sign[j_side][j_station]
                * state.dij[ite2][j];
        }
    }

    let ule1_a = state.ue_inviscid_d_alpha[1][2];
    let ule2_a = state.ue_inviscid_d_alpha[2][2];

    state.transition_node_fraction[1] = 0.0;
    state.transition_node_fraction[2] = 0.0;

    // the global system
    let mut va = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vb = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vdel = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vm = vec![vec![[0.0f64; 3]; nsys]; nsys];
    let mut vz = [[0.0f64; 2]; 3];

    // COM1/COM2 and the XT sensitivities are COMMON
    let mut s1 = std::mem::take(&mut state.station1);
    let mut s2 = std::mem::take(&mut state.station2);
    let mut trloc = std::mem::take(&mut state.transition);
    let mut sys = IntervalSystem::default();
    let (mut ami, mut cti) = (0.0, 0.0);
    let mut trforc = false;

    // Go over each boundary layer/wake
    for side in 1..=2 {
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
        let amcrit = state.ncrit[side];

        // set forced transition arc length position
        let xiforc = xi_trip(state, side);

        let mut tran;
        let mut turb;

        // Sweep downstream setting up BL equation linearizations
        for i_station in 2..=state.n_stations[side] {
            let i_row = state.i_row[side][i_station];
            let simi = i_station == 2;
            let wake = i_station > state.i_te_station[side];
            tran = i_station == state.i_transition_station[side];
            turb = i_station > state.i_transition_station[side];
            let i = state.i_node[side][i_station];

            // set primary variables for current station
            let xsi = state.xi[side][i_station];
            if i_station < state.i_transition_station[side] {
                ami = state.sqrtctau[side][i_station];
            }
            if i_station >= state.i_transition_station[side] {
                cti = state.sqrtctau[side][i_station];
            }
            let uei = state.ue[side][i_station];
            let thi = state.theta[side][i_station];
            let mdi = state.mass_defect[side][i_station];
            let dsi = mdi / uei;
            let dswaki = if wake {
                state.wake_gap[i_station - state.i_te_station[side]]
            } else {
                0.0
            };

            // set derivatives of DSI (= D2)
            let d2_m2 = 1.0 / uei;
            let d2_u2 = -dsi / uei;
            for j_side in 1..=2 {
                for j_station in 2..=state.n_stations[j_side] {
                    let j = state.i_node[j_side][j_station];
                    let j_row = state.i_row[j_side][j_station];
                    u2_m[j_row] = -state.velocity_sign[side][i_station]
                        * state.velocity_sign[j_side][j_station]
                        * state.dij[i][j];
                    d2_m[j_row] = d2_u2 * u2_m[j_row];
                }
            }
            d2_m[i_row] += d2_m2;

            let u2_a = state.ue_inviscid_d_alpha[side][i_station];
            let d2_a = d2_u2 * u2_a;

            // "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
            let due2 = state.ue[side][i_station] - usav[side][i_station];
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
            if i_station == state.i_te_station[side] + 1 {
                // define quantities at start of wake, adding TE base thickness to Dstar
                let (it1, it2) = (state.i_te_station[1], state.i_te_station[2]);
                let tte = state.theta[1][it1] + state.theta[2][it2];
                let dte = state.dstar[1][it1] + state.dstar[2][it2] + state.te_thickness_normal;
                let cte =
                    (state.sqrtctau[1][it1] * state.theta[1][it1] + state.sqrtctau[2][it2] * state.theta[2][it2]) / tte;
                assemble_te_system(&mut sys, &mut s2, cte, tte, dte, &params);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                let dte_mte1 = 1.0 / state.ue[1][it1];
                let dte_ute1 = -state.dstar[1][it1] / state.ue[1][it1];
                let dte_mte2 = 1.0 / state.ue[2][it2];
                let dte_ute2 = -state.dstar[2][it2] / state.ue[2][it2];
                cte_cte1 = state.theta[1][it1] / tte;
                cte_cte2 = state.theta[2][it2] / tte;
                cte_tte1 = (state.sqrtctau[1][it1] - cte) / tte;
                cte_tte2 = (state.sqrtctau[2][it2] - cte) / tte;

                // re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
                for j_side in 1..=2 {
                    for j_station in 2..=state.n_stations[j_side] {
                        let j_row = state.i_row[j_side][j_station];
                        d1_m[j_row] = dte_ute1 * ute1_m[j_row] + dte_ute2 * ute2_m[j_row];
                    }
                }
                d1_m[jvte1] += dte_mte1;
                d1_m[jvte2] += dte_mte2;

                // "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
                due1 = 0.0;
                dds1 = dte_ute1 * (state.ue[1][it1] - usav[1][it1]) + dte_ute2 * (state.ue[2][it2] - usav[2][it2]);
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
            state.tau[side][i_station] = 0.5 * s2.rho * s2.ue * s2.ue * s2.cf;
            state.dissipation[side][i_station] = s2.rho * s2.ue * s2.ue * s2.ue * s2.cdiss * s2.hstar * 0.5;
            state.sqrtctaueq[side][i_station] = s2.sqrtctaueq;
            state.delta[side][i_station] = s2.delta;
            state.us_plot_scale[side][i_station] = 1.60 / (1.0 + s2.us);

            // set XI sensitivities wrt LE Ue changes
            let (xi_ule1, xi_ule2) = if side == 1 {
                (state.s_stagnation_d_gamma_node0, -state.s_stagnation_d_gamma_node1)
            } else {
                (-state.s_stagnation_d_gamma_node0, state.s_stagnation_d_gamma_node1)
            };

            // stuff BL system coefficients into main Jacobian matrix
            for k in 0..3 {
                let (vs1, vs2) = (&sys.jacobian_station1[k], &sys.jacobian_station2[k]);
                let vsx = sys.residual_d_xi[k];
                for j_row in 1..=nsys {
                    vm[i_row - 1][j_row - 1][k] = vs1[2] * d1_m[j_row]
                        + vs1[3] * u1_m[j_row]
                        + vs2[2] * d2_m[j_row]
                        + vs2[3] * u2_m[j_row]
                        + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * ule1_m[j_row] + xi_ule2 * ule2_m[j_row]);
                }
                vb[i_row - 1][k][0] = vs1[0];
                vb[i_row - 1][k][1] = vs1[1];
                va[i_row - 1][k][0] = vs2[0];
                va[i_row - 1][k][1] = vs2[1];
                if state.alpha_specified {
                    vdel[i_row - 1][k][1] = sys.residual_d_re[k] * re_clmr + sys.residual_d_machsqd[k] * msq_clmr;
                } else {
                    vdel[i_row - 1][k][1] = (vs1[3] * u1_a + vs1[2] * d1_a)
                        + (vs2[3] * u2_a + vs2[2] * d2_a)
                        + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * ule1_a + xi_ule2 * ule2_a);
                }
                vdel[i_row - 1][k][0] = sys.residual[k]
                    + (vs1[3] * due1 + vs1[2] * dds1)
                    + (vs2[3] * due2 + vs2[2] * dds2)
                    + (vs1[4] + vs2[4] + vsx) * (xi_ule1 * dule1 + xi_ule2 * dule2);
            }

            if i_station == state.i_te_station[side] + 1 {
                // redefine coefficients for TTE, DTE, etc
                for k in 0..3 {
                    let vs1 = &sys.jacobian_station1[k];
                    vz[k][0] = vs1[0] * cte_cte1;
                    vz[k][1] = vs1[0] * cte_tte1 + vs1[1] * tte_tte1;
                    vb[i_row - 1][k][0] = vs1[0] * cte_cte2;
                    vb[i_row - 1][k][1] = vs1[0] * cte_tte2 + vs1[1] * tte_tte2;
                }
            }

            // turbulent intervals will follow if currently at transition interval
            if tran {
                turb = true;

                // save transition location
                state.i_transition_station[side] = i_station;
                state.transition_forced[side] = trforc;
                state.xi_transition[side] = trloc.xi_transition;

                // interpolate airfoil geometry to find transition x/c (for user output)
                let str = if side == 1 {
                    state.s_stagnation - trloc.xi_transition
                } else {
                    state.s_stagnation + trloc.xi_transition
                };
                let chx = state.x_te - state.x_le;
                let chy = state.y_te - state.y_le;
                let chsq = chx * chx + chy * chy;
                let n = state.n_foil_nodes;
                let xtr = spline_value(str, &state.x[1..=n], &state.dxds[1..=n], &state.s[1..=n]);
                let ytr = spline_value(str, &state.y[1..=n], &state.dyds[1..=n], &state.s[1..=n]);
                state.x_transition[side] = ((xtr - state.x_le) * chx + (ytr - state.y_le) * chy) / chsq;
                state.y_transition[side] = ((ytr - state.y_le) * chx - (xtr - state.x_le) * chy) / chsq;
            }
            // (TRAN = .FALSE.; both flags are recomputed from ITRAN at the top of the sweep)
            let _ = turb;

            if i_station == state.i_te_station[side] {
                // set "2" variables at TE to wake correlations for next station
                // (BLVAR(3); BLMID(3) only sets the interval CFM, which the next BLSYS recomputes)
                s2.set_closure_variables(FlowRegime::Wake, &params);
            }

            for j_side in 1..=2 {
                for j_station in 2..=state.n_stations[j_side] {
                    let j_row = state.i_row[j_side][j_station];
                    u1_m[j_row] = u2_m[j_row];
                    d1_m[j_row] = d2_m[j_row];
                }
            }
            u1_a = u2_a;
            d1_a = d2_a;
            due1 = due2;
            dds1 = dds2;

            if i_station == state.i_transition_station[side] && s2.xi > s1.xi {
                let frac = (trloc.xi_transition - s1.xi) / (s2.xi - s1.xi);
                state.transition_node_fraction[side] = if side == 1 {
                    (state.i_stagnation_node as i64 - state.i_transition_station[side] as i64 + 3) as f64 - frac
                } else {
                    (state.i_stagnation_node as i64 + state.i_transition_station[side] as i64 - 2) as f64 + frac
                };
            }

            // set BL variables for next station
            s1 = s2.clone();
        }
        // 'Side IS forced/free transition at x/c = XOCTR(IS) ITRAN(IS)'
    }

    state.station1 = s1;
    state.station2 = s2;
    state.transition = trloc;

    let sys = NewtonSystem {
        n_rows: nsys,
        diagonal: va,
        subdiagonal: vb,
        rhs: vdel,
        mass_influence: vm,
        te_block: vz,
        i_te_row_upper: Some(state.i_row[1][state.i_te_station[1]] - 1),
        i_wake_row: Some(state.i_row[2][state.i_te_station[2] + 1] - 1),
        elimination_threshold: state.elimination_threshold,
        s_total: Some(state.s[state.n_foil_nodes] - state.s[1]),
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
