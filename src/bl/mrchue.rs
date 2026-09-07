//! MRCHUE (xbl.f): march the BLs and wake in direct mode using the UEDG array, switching to an
//! inverse (prescribed-Hk) step where the direct step would separate. Line-for-line on
//! `BlState`, with an optional trace at the same points as the reference instrumentation.

use crate::bl::blsys::{blsys, tesys, IntervalFlags};
use crate::bl::gauss::gauss_solve_4x4;
use crate::bl::system::{limit_dstar, trchek, BLLocalSystem, FlowParameters, TransitionResult};
use crate::solver::blstate::BlState;
use crate::solver::pointers::xifset;

/// One Newton iteration of one station, as the reference trace records it.
#[derive(Debug, Clone, Default)]
pub struct StationIter {
    pub is: usize,
    pub ibl: usize,
    pub itbl: usize,
    /// AMPL1 AMPL2 XT AMCRIT, TRAN, ITRAN(IS) — logged right after TRCHEK
    pub ampl: [f64; 4],
    pub tran: bool,
    pub itran: usize,
    /// X2 U2 T2 D2 S2
    pub primary: [f64; 5],
    /// M2 H2 HK2 RT2 V2
    pub kinematic: [f64; 5],
    /// HS2 US2 CQ2 CF2 DI2
    pub closure: [f64; 5],
    /// VSREZ(1..4) before GAUSS
    pub residual: [f64; 4],
    /// VS2 rows 1..3, 5 columns, before GAUSS
    pub vs2: [[f64; 5]; 3],
    /// VSREZ(1..4) after GAUSS
    pub solution: [f64; 4],
    pub dmax: f64,
    pub rlx: f64,
    /// CTI THI DSI UEI after the update and DSLIM
    pub updated: [f64; 4],
    /// false on the direct→inverse "try again" path, which performs no update
    pub has_update: bool,
    pub converged: bool,
}

#[derive(Debug, Clone, Default)]
pub struct MrchueTrace {
    pub iters: Vec<StationIter>,
}

/// MRCHUE. Requires the pointer layer (XSSI, IPAN, IBLTE, NBL, WGAP), UEDG initialised
/// (UINV on the first call), ANTE, XSTRIP and the transition thresholds.
pub fn mrchue(st: &mut BlState, params: &FlowParameters, acrit: [f64; 3], mut trace: Option<&mut MrchueTrace>) {
    // shape parameters for separation criteria
    let hlmax = 3.8;
    let htmax = 2.5;

    // COM1/COM2 and XT are COMMON in XFOIL: they persist across sides (the side-2 similarity
    // station's BLSYS copies COM2 into COM1 before anything reads it, so this is trace-faithful
    // rather than algorithmic).
    let mut s1 = std::mem::take(&mut st.com1);
    let mut s2 = std::mem::take(&mut st.com2);
    let mut trloc = std::mem::take(&mut st.trloc);
    for is in 1..=2 {
        let amcrit = acrit[is];

        // set forced transition arc length position
        let xiforc = xifset(st, is);

        // initialize similarity station with Thwaites' formula
        let ibl0 = 2;
        let xsi0 = st.xssi[is][ibl0];
        let uei0 = st.uedg[is][ibl0];
        let bule = 1.0_f64;
        let ucon = uei0 / xsi0.powf(bule);
        let tsq = 0.45 / (ucon * (5.0 * bule + 1.0) * params.re) * xsi0.powf(1.0 - bule);
        let mut thi = tsq.sqrt();
        let mut dsi = 2.2 * thi;
        let mut ami = 0.0;

        // initialize Ctau for first turbulent station
        let mut cti = 0.03;

        let mut tran = false;
        let mut turb = false;
        st.itran[is] = st.iblte[is];

        let mut sys = BLLocalSystem::default();
        let mut trforc = false;
        // TE quantities defined at the first wake station (used again by the fallback)
        let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);

        // march downstream
        for ibl in 2..=st.nbl[is] {
            let ibm = ibl - 1;
            let simi = ibl == 2;
            let wake = ibl > st.iblte[is];

            // prescribed quantities
            let xsi = st.xssi[is][ibl];
            let mut uei = st.uedg[is][ibl];
            let dswaki = if wake { st.wgap[ibl - st.iblte[is]] } else { 0.0 };

            let mut direct = true;
            let mut htarg = 0.0;
            let mut dmax = 0.0;
            let mut converged = false;
            let mut hk2_snapshot = 0.0;

            // Newton iteration loop for current station
            for itbl in 1..=25 {
                // assemble 10x3 linearized system at the "1" and "2" stations
                s2.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, params);
                s2.set_kinematic_variables(params);
                // the reference trace records AMPL1/AMPL2/XT/TRAN/ITRAN here, before TRCHEK
                let pre = ([s1.ampl, s2.ampl, trloc.xt, amcrit], tran, st.itran[is]);

                // check for transition and set appropriate flags and things
                if !simi && !turb {
                    match trchek(&s1, &s2, s1.ampl, amcrit, xiforc, params) {
                        TransitionResult::NoTransition { ampl2 } => {
                            ami = ampl2;
                            tran = false;
                            trloc.xt = s2.xi; // TRCHEK2 leaves XT = X2 (and the XT_* as they were)
                            st.itran[is] = ibl + 2;
                        }
                        TransitionResult::FreeTransition { location, ampl2 } => {
                            ami = ampl2;
                            tran = true;
                            trforc = false;
                            trloc = location;
                            st.itran[is] = ibl;
                            if cti <= 0.0 {
                                cti = 0.03;
                                s2.sqrtctau = cti;
                            }
                        }
                        TransitionResult::ForcedTransition { location } => {
                            tran = true;
                            trforc = true;
                            trloc = location;
                            st.itran[is] = ibl;
                            if cti <= 0.0 {
                                cti = 0.03;
                                s2.sqrtctau = cti;
                            }
                        }
                    }
                    s2.ampl = ami;
                }

                let flags = IntervalFlags { simi, tran, turb, wake };
                if ibl == st.iblte[is] + 1 {
                    tte = st.thet[1][st.iblte[1]] + st.thet[2][st.iblte[2]];
                    dte = st.dstr[1][st.iblte[1]] + st.dstr[2][st.iblte[2]] + st.ante;
                    cte = (st.ctau[1][st.iblte[1]] * st.thet[1][st.iblte[1]]
                        + st.ctau[2][st.iblte[2]] * st.thet[2][st.iblte[2]])
                        / tte;
                    tesys(&mut sys, &mut s2, cte, tte, dte, params);
                } else {
                    blsys(&mut sys, &mut s1, &mut s2, flags, Some(&trloc), amcrit, params);
                }
                hk2_snapshot = s2.hk;

                if direct {
                    // try direct mode (set dUe = 0 in currently empty 4th line)
                    sys.vs2[3][0] = 0.0;
                    sys.vs2[3][1] = 0.0;
                    sys.vs2[3][2] = 0.0;
                    sys.vs2[3][3] = 1.0;
                    sys.vsrez[3] = 0.0;
                } else {
                    // inverse mode (force Hk to prescribed value HTARG)
                    sys.vs2[3][0] = 0.0;
                    sys.vs2[3][1] = s2.hk_d_theta;
                    sys.vs2[3][2] = s2.hk_d_dstar;
                    sys.vs2[3][3] = s2.hk_d_ue;
                    sys.vsrez[3] = htarg - s2.hk;
                }

                let mut rec = StationIter {
                    is,
                    ibl,
                    itbl,
                    ampl: pre.0,
                    tran: pre.1,
                    itran: pre.2,
                    primary: [s2.xi, s2.ue, s2.theta, s2.dstar, s2.sqrtctau],
                    kinematic: [s2.machsqd_edge, s2.h, s2.hk, s2.retheta, s2.nu],
                    closure: [s2.hstar, s2.us, s2.sqrtctaueq, s2.cf, s2.cdiss],
                    residual: sys.vsrez,
                    vs2: [sys.vs2[0], sys.vs2[1], sys.vs2[2]],
                    has_update: true,
                    ..Default::default()
                };

                // solve Newton system for current "2" station
                let mut z = [[0.0f64; 4]; 4];
                let mut r = [0.0f64; 4];
                for k in 0..4 {
                    r[k] = sys.vsrez[k];
                    for l in 0..4 {
                        z[k][l] = sys.vs2[k][l];
                    }
                }
                gauss_solve_4x4(&mut z, &mut r);
                rec.solution = r;

                let rlx;
                if direct {
                    // determine max changes and underrelax if necessary
                    dmax = (r[1] / thi).abs().max((r[2] / dsi).abs());
                    if ibl < st.itran[is] {
                        dmax = dmax.max((r[0] / 10.0).abs());
                    }
                    if ibl >= st.itran[is] {
                        dmax = dmax.max((r[0] / cti).abs());
                    }
                    rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

                    // see if direct mode is not applicable
                    if ibl != st.iblte[is] + 1 {
                        // calculate resulting kinematic shape parameter Hk
                        let msq = uei * uei * params.h_stagnation_inv
                            / (params.gamma_gas_m1 * (1.0 - 0.5 * uei * uei * params.h_stagnation_inv));
                        let htest = (dsi + rlx * r[2]) / (thi + rlx * r[1]);
                        let (hktest, _, _) = crate::bl::hkin(htest, msq);

                        // decide whether to do direct or inverse problem based on Hk
                        let hmax = if ibl < st.itran[is] { hlmax } else { htmax };
                        direct = hktest < hmax;
                    }

                    if direct {
                        // update as usual
                        if ibl >= st.itran[is] {
                            cti += rlx * r[0];
                        }
                        thi += rlx * r[1];
                        dsi += rlx * r[2];
                    } else {
                        // set prescribed Hk for inverse calculation at the current station
                        let hmax = if ibl < st.itran[is] { hlmax } else { htmax };
                        htarg = if ibl < st.itran[is] {
                            // laminar case: relatively slow increase in Hk downstream
                            s1.hk + 0.03 * (s2.xi - s1.xi) / s1.theta
                        } else if ibl == st.itran[is] {
                            // transition interval: weighted laminar and turbulent case
                            s1.hk + (0.03 * (trloc.xt - s1.xi) - 0.15 * (s2.xi - trloc.xt)) / s1.theta
                        } else if wake {
                            // turbulent wake case: asymptotic wake behavior with approximate Backward Euler
                            let cnst = 0.03 * (s2.xi - s1.xi) / s1.theta;
                            let hk1 = s1.hk;
                            let mut hk2 = hk1;
                            for _ in 0..3 {
                                hk2 -=
                                    (hk2 + cnst * (hk2 - 1.0).powi(3) - hk1) / (1.0 + 3.0 * cnst * (hk2 - 1.0).powi(2));
                            }
                            hk2
                        } else {
                            // turbulent case: relatively fast decrease in Hk downstream
                            s1.hk - 0.15 * (s2.xi - s1.xi) / s1.theta
                        };
                        // limit specified Hk to something reasonable
                        htarg = if wake { htarg.max(1.01) } else { htarg.max(hmax) };
                        // try again with prescribed Hk
                        rec.dmax = dmax;
                        rec.rlx = rlx;
                        rec.has_update = false;
                        if let Some(t) = trace.as_mut() {
                            t.iters.push(rec);
                        }
                        continue;
                    }
                } else {
                    // added Ue clamp   MD  3 Apr 03
                    dmax = (r[1] / thi).abs().max((r[2] / dsi).abs()).max((r[3] / uei).abs());
                    if ibl >= st.itran[is] {
                        dmax = dmax.max((r[0] / cti).abs());
                    }
                    rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };
                    // update variables
                    if ibl >= st.itran[is] {
                        cti += rlx * r[0];
                    }
                    thi += rlx * r[1];
                    dsi += rlx * r[2];
                    uei += rlx * r[3];
                }

                // eliminate absurd transients
                if ibl >= st.itran[is] {
                    cti = cti.min(0.30);
                    cti = cti.max(0.0000001);
                }
                let hklim = if ibl <= st.iblte[is] { 1.02 } else { 1.00005 };
                let msq = uei * uei * params.h_stagnation_inv
                    / (params.gamma_gas_m1 * (1.0 - 0.5 * uei * uei * params.h_stagnation_inv));
                let mut dsw = dsi - dswaki;
                limit_dstar(&mut dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                rec.dmax = dmax;
                rec.rlx = rlx;
                rec.updated = [cti, thi, dsi, uei];
                if dmax <= 1.0e-5 {
                    converged = true;
                    rec.converged = true;
                }
                if let Some(t) = trace.as_mut() {
                    t.iters.push(rec);
                }
                if converged {
                    break;
                }
            }

            if !converged {
                // 'MRCHUE: Convergence failed at IBL side IS Res = DMAX'
                // the current unconverged solution might still be reasonable...
                if dmax > 0.1 {
                    // the current solution is garbage --> extrapolate values instead
                    if ibl > 3 {
                        if ibl <= st.iblte[is] {
                            thi = st.thet[is][ibm] * (st.xssi[is][ibl] / st.xssi[is][ibm]).powf(0.5);
                            dsi = st.dstr[is][ibm] * (st.xssi[is][ibl] / st.xssi[is][ibm]).powf(0.5);
                        } else if ibl == st.iblte[is] + 1 {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                        } else {
                            thi = st.thet[is][ibm];
                            let ratlen = (st.xssi[is][ibl] - st.xssi[is][ibm]) / (10.0 * st.dstr[is][ibm]);
                            dsi = (st.dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                        }
                        if ibl == st.itran[is] {
                            cti = 0.05;
                        }
                        if ibl > st.itran[is] {
                            cti = st.ctau[is][ibm];
                        }
                        uei = st.uedg[is][ibl];
                        if ibl > 2 && ibl < st.nbl[is] {
                            uei = 0.5 * (st.uedg[is][ibl - 1] + st.uedg[is][ibl + 1]);
                        }
                    }
                }
                // label 109
                s2.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, params);
                s2.set_kinematic_variables(params);
                // check for transition and set appropriate flags and things
                if !simi && !turb {
                    match trchek(&s1, &s2, s1.ampl, amcrit, xiforc, params) {
                        TransitionResult::NoTransition { ampl2 } => {
                            ami = ampl2;
                            tran = false;
                            trloc.xt = s2.xi;
                            st.itran[is] = ibl + 2;
                        }
                        TransitionResult::FreeTransition { location, ampl2 } => {
                            ami = ampl2;
                            tran = true;
                            trforc = false;
                            trloc = location;
                            st.itran[is] = ibl;
                        }
                        TransitionResult::ForcedTransition { location } => {
                            tran = true;
                            trforc = true;
                            trloc = location;
                            st.itran[is] = ibl;
                        }
                    }
                    s2.ampl = ami;
                }
                // set all other extrapolated values for current station (BLVAR/BLMID)
                let ityp = if wake {
                    crate::bl::system::FlowRegime::Wake
                } else if ibl < st.itran[is] {
                    crate::bl::system::FlowRegime::Laminar
                } else {
                    crate::bl::system::FlowRegime::Turbulent
                };
                s2.set_closure_variables(ityp, params);
                hk2_snapshot = s2.hk;
            }
            let _ = hk2_snapshot;

            // label 110: store primary variables
            st.ctau[is][ibl] = if ibl < st.itran[is] { ami } else { cti };
            st.thet[is][ibl] = thi;
            st.dstr[is][ibl] = dsi;
            st.uedg[is][ibl] = uei;
            st.mass[is][ibl] = dsi * uei;
            st.tau[is][ibl] = 0.5 * s2.rho * s2.ue * s2.ue * s2.cf;
            st.dis[is][ibl] = s2.rho * s2.ue * s2.ue * s2.ue * s2.cdiss * s2.hstar * 0.5;
            st.ctq[is][ibl] = s2.sqrtctaueq;
            st.delt[is][ibl] = s2.delta;
            st.tstr[is][ibl] = s2.hstar * s2.theta;

            // set "1" variables to "2" variables for next streamwise station
            s2.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, params);
            s2.set_kinematic_variables(params);
            s1 = s2.clone();

            // turbulent intervals will follow transition interval or TE
            if tran || ibl == st.iblte[is] {
                turb = true;
                // save transition location (TFORCE, XSSITR)
                let _ = trforc;
            }
            tran = false;

            if ibl == st.iblte[is] {
                thi = st.thet[1][st.iblte[1]] + st.thet[2][st.iblte[2]];
                dsi = st.dstr[1][st.iblte[1]] + st.dstr[2][st.iblte[2]] + st.ante;
            }
        }
    }
    st.com1 = s1;
    st.com2 = s2;
    st.trloc = trloc;
}
