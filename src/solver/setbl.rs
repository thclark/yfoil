//! SETBL (xbl.f): sets up the BL Newton system coefficients for the current BL variables and
//! the edge velocities, and incorporates the local BL system coefficients into the global
//! Newton system (VA/VB/VM/VDEL/VZ) that BLSOLV consumes. Line-for-line on `BlState`.
//!
//! Also MRCL (xfoil.f) — the Mach/Re dependence on CL selected by MATYP/RETYP — which SETBL
//! calls first. `IDAMPV = IDAMP` is pinned at 0 (the envelope e^N `DAMPL` model).

use crate::bl::blsolv::BlsolvInput;
use crate::bl::blsys::{assemble_interval_system, assemble_te_system, IntervalFlags};
use crate::bl::mrchdu::mrchdu;
use crate::bl::mrchue::mrchue;
use crate::bl::system::{check_transition, FlowParameters, FlowRegime, IntervalSystem, TransitionCheck};
use crate::geometry::seval;
use crate::solver::blstate::BlState;
use crate::solver::pointers::xifset;
use crate::solver::velocity::ueset;

/// Everything SETBL produces besides the arrays it writes into `BlState`.
#[derive(Debug, Clone)]
pub struct SetblResult {
    /// The global BL Newton system as BLSOLV sees it (VA, VB, VDEL, VM, VZ, NSYS, IVTE1, IVZ).
    pub sys: BlsolvInput,
    /// The BL parameters SETBL derived (REYBL, HSTINV, ...): the marches and UPDATE share them.
    pub params: FlowParameters,
    /// RE_CLMR, MSQ_CLMR: d(Re)/d(CL) and d(M²)/d(CL) for the fixed-CL sensitivity column.
    pub re_clmr: f64,
    pub msq_clmr: f64,
    /// MA_CLMR: d(M)/d(CL) (MRCL's M_CLS); VISCAL's MINF_CL is the same quantity
    pub ma_clmr: f64,
    /// DULE1, DULE2: the LE Ue mismatch between UEDG and USAV = UINV + DIJ·MASS, per side.
    pub dule: [f64; 3],
}

/// MRCL: sets the actual Mach and Reynolds numbers from the unit-CL values and the specified
/// CLS according to MATYP/RETYP. Returns (M_CLS, R_CLS).
pub fn mrcl(st: &mut BlState, cls: f64) -> (f64, f64) {
    let cla = cls.max(0.000001);
    if st.retyp < 1 || st.retyp > 3 {
        // 'MRCL:  Illegal Re(CL) dependence trigger. Setting fixed Re.'
        st.retyp = 1;
    }
    if st.matyp < 1 || st.matyp > 3 {
        // 'MRCL:  Illegal Mach(CL) dependence trigger. Setting fixed Mach.'
        st.matyp = 1;
    }
    let mut m_cls;
    match st.matyp {
        1 => {
            st.minf = st.minf1;
            m_cls = 0.0;
        }
        2 => {
            st.minf = st.minf1 / cla.sqrt();
            m_cls = -0.5 * st.minf / cla;
        }
        _ => {
            st.minf = st.minf1;
            m_cls = 0.0;
        }
    }
    let mut r_cls;
    match st.retyp {
        1 => {
            st.reinf = st.reinf1;
            r_cls = 0.0;
        }
        2 => {
            st.reinf = st.reinf1 / cla.sqrt();
            r_cls = -0.5 * st.reinf / cla;
        }
        _ => {
            st.reinf = st.reinf1 / cla;
            r_cls = -st.reinf / cla;
        }
    }
    if st.minf >= 0.99 {
        // 'MRCL: CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99'
        st.minf = 0.99;
        m_cls = 0.0;
    }
    let mut rrat = 1.0;
    if st.reinf1 > 0.0 {
        rrat = st.reinf / st.reinf1;
    }
    if rrat > 100.0 {
        // 'MRCL: CL too low for chosen Re(CL) dependence — artificially limiting Re'
        st.reinf = st.reinf1 * 100.0;
        r_cls = 0.0;
    }
    (m_cls, r_cls)
}

/// SETBL. Requires the pointer layer, DIJ, UINV/UINV_A, the BL arrays (initialised by MRCHUE
/// on the first call when `st.lblini` is false), SST_GO/SST_GP of the current stagnation point
/// and the control fields (`lalfa`, `cl`/`clspec`, `matyp`/`retyp`, `minf1`/`reinf1`, `acrit`,
/// `vaccel`). Writes TAU/DIS/CTQ/DELT/USLP, ITRAN/XSSITR/TFORCE, XOCTR/YOCTR/TINDEX and leaves
/// UEDG holding the *marched* Ue (USAV = UINV + DIJ·MASS is the mismatch reference).
pub fn setbl(st: &mut BlState) -> SetblResult {
    // set the CL used to define Mach, Reynolds numbers
    let clmr = if st.lalfa { st.cl } else { st.clspec };

    // set current MINF(CL)
    let (ma_clmr, re_clmr) = mrcl(st, clmr);
    let msq_clmr = 2.0 * st.minf * ma_clmr;

    // set compressibility parameter TKLAM and derivative TK_MSQ (COMSET), gas constant, the
    // parameters for compressibility correction, stagnation density and 1/enthalpy, and the
    // Reynolds number based on freestream density, velocity, viscosity
    let mut params = FlowParameters::new(st.minf, st.reinf, st.gamma);

    // IDAMPV = IDAMP
    params.idampv = st.idamp;
    // save TE thickness
    let _dwte = st.wgap[1];

    if !st.lblini {
        // initialize BL by marching with Ue (fudge at separation)
        let acrit = st.acrit;
        mrchue(st, &params, acrit, None);
        st.lblini = true;
    }

    // march BL with current Ue and Ds to establish transition
    let acrit = st.acrit;
    mrchdu(st, &params, acrit, None);

    let mut usav: [Vec<f64>; 3] = [Vec::new(), st.uedg[1].clone(), st.uedg[2].clone()];
    ueset(st);
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            std::mem::swap(&mut usav[is][ibl], &mut st.uedg[is][ibl]);
        }
    }

    let ile1 = st.ipan[1][2];
    let ile2 = st.ipan[2][2];
    let ite1 = st.ipan[1][st.iblte[1]];
    let ite2 = st.ipan[2][st.iblte[2]];

    let jvte1 = st.isys[1][st.iblte[1]];
    let jvte2 = st.isys[2][st.iblte[2]];

    let dule1 = st.uedg[1][2] - usav[1][2];
    let dule2 = st.uedg[2][2] - usav[2][2];

    let nsys = st.nsys;

    // set LE and TE Ue sensitivities wrt all m values
    let mut ule1_m = vec![0.0; nsys + 1];
    let mut ule2_m = vec![0.0; nsys + 1];
    let mut ute1_m = vec![0.0; nsys + 1];
    let mut ute2_m = vec![0.0; nsys + 1];
    for js in 1..=2 {
        for jbl in 2..=st.nbl[js] {
            let j = st.ipan[js][jbl];
            let jv = st.isys[js][jbl];
            ule1_m[jv] = -st.vti[1][2] * st.vti[js][jbl] * st.dij[ile1][j];
            ule2_m[jv] = -st.vti[2][2] * st.vti[js][jbl] * st.dij[ile2][j];
            ute1_m[jv] = -st.vti[1][st.iblte[1]] * st.vti[js][jbl] * st.dij[ite1][j];
            ute2_m[jv] = -st.vti[2][st.iblte[2]] * st.vti[js][jbl] * st.dij[ite2][j];
        }
    }

    let ule1_a = st.uinv_a[1][2];
    let ule2_a = st.uinv_a[2][2];

    st.tindex[1] = 0.0;
    st.tindex[2] = 0.0;

    // the global system
    let mut va = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vb = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vdel = vec![[[0.0f64; 2]; 3]; nsys];
    let mut vm = vec![vec![[0.0f64; 3]; nsys]; nsys];
    let mut vz = [[0.0f64; 2]; 3];

    // COM1/COM2 and the XT sensitivities are COMMON
    let mut s1 = std::mem::take(&mut st.com1);
    let mut s2 = std::mem::take(&mut st.com2);
    let mut trloc = std::mem::take(&mut st.trloc);
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
        let amcrit = st.acrit[is];

        // set forced transition arc length position
        let xiforc = xifset(st, is);

        let mut tran;
        let mut turb;

        // Sweep downstream setting up BL equation linearizations
        for ibl in 2..=st.nbl[is] {
            let iv = st.isys[is][ibl];
            let simi = ibl == 2;
            let wake = ibl > st.iblte[is];
            tran = ibl == st.itran[is];
            turb = ibl > st.itran[is];
            let i = st.ipan[is][ibl];

            // set primary variables for current station
            let xsi = st.xssi[is][ibl];
            if ibl < st.itran[is] {
                ami = st.ctau[is][ibl];
            }
            if ibl >= st.itran[is] {
                cti = st.ctau[is][ibl];
            }
            let uei = st.uedg[is][ibl];
            let thi = st.thet[is][ibl];
            let mdi = st.mass[is][ibl];
            let dsi = mdi / uei;
            let dswaki = if wake { st.wgap[ibl - st.iblte[is]] } else { 0.0 };

            // set derivatives of DSI (= D2)
            let d2_m2 = 1.0 / uei;
            let d2_u2 = -dsi / uei;
            for js in 1..=2 {
                for jbl in 2..=st.nbl[js] {
                    let j = st.ipan[js][jbl];
                    let jv = st.isys[js][jbl];
                    u2_m[jv] = -st.vti[is][ibl] * st.vti[js][jbl] * st.dij[i][j];
                    d2_m[jv] = d2_u2 * u2_m[jv];
                }
            }
            d2_m[iv] += d2_m2;

            let u2_a = st.uinv_a[is][ibl];
            let d2_a = d2_u2 * u2_a;

            // "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
            let due2 = st.uedg[is][ibl] - usav[is][ibl];
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
            if ibl == st.iblte[is] + 1 {
                // define quantities at start of wake, adding TE base thickness to Dstar
                let (it1, it2) = (st.iblte[1], st.iblte[2]);
                let tte = st.thet[1][it1] + st.thet[2][it2];
                let dte = st.dstr[1][it1] + st.dstr[2][it2] + st.ante;
                let cte = (st.ctau[1][it1] * st.thet[1][it1] + st.ctau[2][it2] * st.thet[2][it2]) / tte;
                assemble_te_system(&mut sys, &mut s2, cte, tte, dte, &params);

                tte_tte1 = 1.0;
                tte_tte2 = 1.0;
                let dte_mte1 = 1.0 / st.uedg[1][it1];
                let dte_ute1 = -st.dstr[1][it1] / st.uedg[1][it1];
                let dte_mte2 = 1.0 / st.uedg[2][it2];
                let dte_ute2 = -st.dstr[2][it2] / st.uedg[2][it2];
                cte_cte1 = st.thet[1][it1] / tte;
                cte_cte2 = st.thet[2][it2] / tte;
                cte_tte1 = (st.ctau[1][it1] - cte) / tte;
                cte_tte2 = (st.ctau[2][it2] - cte) / tte;

                // re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
                for js in 1..=2 {
                    for jbl in 2..=st.nbl[js] {
                        let jv = st.isys[js][jbl];
                        d1_m[jv] = dte_ute1 * ute1_m[jv] + dte_ute2 * ute2_m[jv];
                    }
                }
                d1_m[jvte1] += dte_mte1;
                d1_m[jvte2] += dte_mte2;

                // "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
                due1 = 0.0;
                dds1 = dte_ute1 * (st.uedg[1][it1] - usav[1][it1]) + dte_ute2 * (st.uedg[2][it2] - usav[2][it2]);
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
            st.dis[is][ibl] = s2.rho * s2.ue * s2.ue * s2.ue * s2.cdiss * s2.hstar * 0.5;
            st.ctq[is][ibl] = s2.sqrtctaueq;
            st.delt[is][ibl] = s2.delta;
            st.uslp[is][ibl] = 1.60 / (1.0 + s2.us);

            // set XI sensitivities wrt LE Ue changes
            let (xi_ule1, xi_ule2) = if is == 1 {
                (st.sst_go, -st.sst_gp)
            } else {
                (-st.sst_go, st.sst_gp)
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
                if st.lalfa {
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

            if ibl == st.iblte[is] + 1 {
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
                st.itran[is] = ibl;
                st.tforce[is] = trforc;
                st.xssitr[is] = trloc.xi_transition;

                // interpolate airfoil geometry to find transition x/c (for user output)
                let str = if is == 1 {
                    st.sst - trloc.xi_transition
                } else {
                    st.sst + trloc.xi_transition
                };
                let chx = st.xte - st.xle;
                let chy = st.yte - st.yle;
                let chsq = chx * chx + chy * chy;
                let n = st.n;
                let xtr = seval(str, &st.x[1..=n], &st.xp[1..=n], &st.s[1..=n]);
                let ytr = seval(str, &st.y[1..=n], &st.yp[1..=n], &st.s[1..=n]);
                st.xoctr[is] = ((xtr - st.xle) * chx + (ytr - st.yle) * chy) / chsq;
                st.yoctr[is] = ((ytr - st.yle) * chx - (xtr - st.xle) * chy) / chsq;
            }
            // (TRAN = .FALSE.; both flags are recomputed from ITRAN at the top of the sweep)
            let _ = turb;

            if ibl == st.iblte[is] {
                // set "2" variables at TE to wake correlations for next station
                // (BLVAR(3); BLMID(3) only sets the interval CFM, which the next BLSYS recomputes)
                s2.set_closure_variables(FlowRegime::Wake, &params);
            }

            for js in 1..=2 {
                for jbl in 2..=st.nbl[js] {
                    let jv = st.isys[js][jbl];
                    u1_m[jv] = u2_m[jv];
                    d1_m[jv] = d2_m[jv];
                }
            }
            u1_a = u2_a;
            d1_a = d2_a;
            due1 = due2;
            dds1 = dds2;

            if ibl == st.itran[is] && s2.xi > s1.xi {
                let frac = (trloc.xi_transition - s1.xi) / (s2.xi - s1.xi);
                st.tindex[is] = if is == 1 {
                    (st.ist as i64 - st.itran[is] as i64 + 3) as f64 - frac
                } else {
                    (st.ist as i64 + st.itran[is] as i64 - 2) as f64 + frac
                };
            }

            // set BL variables for next station
            s1 = s2.clone();
        }
        // 'Side IS forced/free transition at x/c = XOCTR(IS) ITRAN(IS)'
    }

    st.com1 = s1;
    st.com2 = s2;
    st.trloc = trloc;

    let sys = BlsolvInput {
        nsys,
        va,
        vb,
        vdel,
        vm,
        vz,
        ivte1: Some(st.isys[1][st.iblte[1]] - 1),
        ivz: Some(st.isys[2][st.iblte[2] + 1] - 1),
        vaccel: st.vaccel,
        arc_length: Some(st.s[st.n] - st.s[1]),
    };
    SetblResult {
        sys,
        params,
        re_clmr,
        msq_clmr,
        ma_clmr,
        dule: [0.0, dule1, dule2],
    }
}
