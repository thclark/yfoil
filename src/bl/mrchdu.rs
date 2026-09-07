//! MRCHDU (xbl.f): march the BLs and wake in mixed mode using the current Ue and Hk. The
//! calculated Ue and Hk lie along a line quasi-normal to the natural Ue-Hk characteristic line
//! of the current BL so that the Goldstein or Levy-Lees singularity is never encountered.
//! Continuous checking of transition onset is performed. Line-for-line on `BlState`, with an
//! optional trace at the same points as the reference instrumentation (`xfoil_mrchdu_trace.dat`).

use crate::bl::blsys::{blsys, tesys, IntervalFlags};
use crate::bl::gauss::gauss_solve_4x4;
use crate::bl::hkin;
use crate::bl::system::{dslim, trchek, BLLocalSystem, FlowParameters, FlowRegime, TransitionResult};
use crate::solver::blstate::BlState;
use crate::solver::pointers::xifset;

/// One Newton iteration of one station, as the reference trace records it.
#[derive(Debug, Clone, Default)]
pub struct MrchduIter {
    pub is: usize,
    pub ibl: usize,
    pub itbl: usize,
    /// AMPL1 AMPL2 XT AMCRIT, TRAN, ITRAN(IS) — logged after BLKIN, before TRCHEK
    pub ampl: [f64; 4],
    pub tran: bool,
    pub itran: usize,
    /// X2 U2 T2 D2 S2
    pub primary: [f64; 5],
    /// M2 H2 HK2 RT2 V2
    pub kinematic: [f64; 5],
    /// HS2 US2 CQ2 CF2 DI2
    pub closure: [f64; 5],
    /// UEREF HKREF
    pub ueref: f64,
    pub hkref: f64,
    /// SENNEW SENS — only where the Ue-Hk characteristic slope is formed
    pub sens: Option<[f64; 2]>,
    /// VSREZ(1..4) before GAUSS
    pub residual: [f64; 4],
    /// VS2 rows 1..4, 5 columns, before GAUSS
    pub vs2: [[f64; 5]; 4],
    /// VSREZ(1..4) after GAUSS
    pub solution: [f64; 4],
    pub dmax: f64,
    pub rlx: f64,
    /// CTI THI DSI UEI AMI after the update and DSLIM
    pub updated: [f64; 5],
    pub converged: bool,
}

#[derive(Debug, Clone, Default)]
pub struct MrchduTrace {
    pub iters: Vec<MrchduIter>,
    /// (IBL, IS, DMAX) for every station whose 25 Newton iterations did not converge
    pub failed: Vec<(usize, usize, f64)>,
}

/// MRCHDU. Requires the pointer layer (XSSI, IBLTE, NBL, WGAP), the current
/// UEDG/THET/DSTR/CTAU, ITRAN from the previous march, ANTE, XSTRIP and the transition
/// thresholds. Updates THET/DSTR/CTAU/UEDG/MASS/TAU/DIS/CTQ/DELT/TSTR, ITRAN, XSSITR and TFORCE.
pub fn mrchdu(st: &mut BlState, params: &FlowParameters, acrit: [f64; 3], mut trace: Option<&mut MrchduTrace>) {
    const DEPS: f64 = 5.0e-6;

    // constant controlling how far Hk is allowed to deviate from the specified value
    let senswt = 1000.0;

    // COM1/COM2 and XT are COMMON, and AMI, SENS/SENNEW, UEREF/HKREF and CTE/TTE/DTE are
    // MRCHDU locals never re-initialised per side: all persist across sides and stations.
    let mut s1 = std::mem::take(&mut st.com1);
    let mut s2 = std::mem::take(&mut st.com2);
    let mut trloc = std::mem::take(&mut st.trloc);
    let mut ami = 0.0;
    let (mut sens, mut sennew) = (0.0, 0.0);
    let (mut ueref, mut hkref) = (0.0, 0.0);
    let (mut cte, mut tte, mut dte) = (0.0, 0.0, 0.0);
    let mut sys = BLLocalSystem::default();
    let mut trforc = false;

    for is in 1..=2 {
        let amcrit = acrit[is];

        // set forced transition arc length position
        let xiforc = xifset(st, is);

        // (leading edge pressure gradient parameter BULE = 1.0 is the constant BLDIF uses)

        // old transition station
        let itrold = st.itran[is];

        let mut tran = false;
        let mut turb = false;
        st.itran[is] = st.iblte[is];

        // march downstream
        for ibl in 2..=st.nbl[is] {
            let ibm = ibl - 1;
            let simi = ibl == 2;
            let wake = ibl > st.iblte[is];

            // initialize current station to existing variables
            let xsi = st.xssi[is][ibl];
            let mut uei = st.uedg[is][ibl];
            let mut thi = st.thet[is][ibl];
            let mut dsi = st.dstr[is][ibl];

            // fixed BUG   MD 7 June 99
            let mut cti;
            if ibl < itrold {
                ami = st.ctau[is][ibl];
                cti = 0.03;
            } else {
                cti = st.ctau[is][ibl];
                if cti <= 0.0 {
                    cti = 0.03;
                }
            }

            let dswaki = if wake { st.wgap[ibl - st.iblte[is]] } else { 0.0 };
            if ibl <= st.iblte[is] {
                dsi = (dsi - dswaki).max(1.02000 * thi) + dswaki;
            }
            if ibl > st.iblte[is] {
                dsi = (dsi - dswaki).max(1.00005 * thi) + dswaki;
            }

            let mut dmax = 0.0;
            let mut converged = false;

            // Newton iteration loop for current station
            for itbl in 1..=25 {
                // assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
                // at the previous "1" station and the current "2" station
                // (the "1" station coefficients will be ignored)
                s2.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, params);
                s2.blkin(params);
                let pre = ([s1.ampl, s2.ampl, trloc.xt, amcrit], tran, st.itran[is]);

                // check for transition and set appropriate flags and things
                if !simi && !turb {
                    match trchek(&s1, &s2, s1.ampl, amcrit, xiforc, params) {
                        TransitionResult::NoTransition { ampl2 } => {
                            ami = ampl2;
                            tran = false;
                            trloc.xt = s2.x;
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

                let mut rec = MrchduIter {
                    is,
                    ibl,
                    itbl,
                    ampl: pre.0,
                    tran: pre.1,
                    itran: pre.2,
                    primary: [s2.x, s2.u, s2.theta, s2.dstar, s2.ctau],
                    kinematic: [s2.msq, s2.h, s2.hk, s2.rt, s2.v],
                    closure: [s2.hs, s2.us, s2.cq, s2.cf, s2.di],
                    ..Default::default()
                };

                // set stuff at first iteration...
                if itbl == 1 {
                    // set "baseline" Ue and Hk for forming  Ue(Hk)  relation
                    ueref = s2.u;
                    hkref = s2.hk;

                    // if current point IBL was turbulent and is now laminar, then...
                    if ibl < st.itran[is] && ibl >= itrold {
                        // extrapolate baseline Hk
                        let uem = st.uedg[is][ibl - 1];
                        let dsm = st.dstr[is][ibl - 1];
                        let thm = st.thet[is][ibl - 1];
                        let msq = uem * uem * params.h_stagnation_inv
                            / (params.gamma_gas_m1 * (1.0 - 0.5 * uem * uem * params.h_stagnation_inv));
                        let (hk, _, _) = hkin(dsm / thm, msq);
                        hkref = hk;
                    }

                    // if current point IBL was laminar, then...
                    if ibl < itrold {
                        // reinitialize or extrapolate Ctau if it's now turbulent
                        if tran {
                            st.ctau[is][ibl] = 0.03;
                        }
                        if turb {
                            st.ctau[is][ibl] = st.ctau[is][ibl - 1];
                        }
                        if tran || turb {
                            cti = st.ctau[is][ibl];
                            s2.ctau = cti;
                        }
                    }
                }

                if simi || ibl == st.iblte[is] + 1 {
                    // for similarity station or first wake point, prescribe Ue
                    sys.vs2[3][0] = 0.0;
                    sys.vs2[3][1] = 0.0;
                    sys.vs2[3][2] = 0.0;
                    sys.vs2[3][3] = s2.u_uei;
                    sys.vsrez[3] = ueref - s2.u;
                } else {
                    // calculate Ue-Hk characteristic slope
                    let mut vtmp = [[0.0f64; 4]; 4];
                    let mut vztmp = [0.0f64; 4];
                    for k in 0..4 {
                        vztmp[k] = sys.vsrez[k];
                        for l in 0..4 {
                            vtmp[k][l] = sys.vs2[k][l];
                        }
                    }
                    // set unit dHk
                    vtmp[3][0] = 0.0;
                    vtmp[3][1] = s2.hk_t;
                    vtmp[3][2] = s2.hk_d;
                    vtmp[3][3] = s2.hk_u * s2.u_uei;
                    vztmp[3] = 1.0;

                    // calculate dUe response
                    gauss_solve_4x4(&mut vtmp, &mut vztmp);

                    // set  SENSWT * (normalized dUe/dHk)
                    sennew = senswt * vztmp[3] * hkref / ueref;
                    if itbl <= 5 {
                        sens = sennew;
                    } else if itbl <= 15 {
                        sens = 0.5 * (sens + sennew);
                    }
                    rec.sens = Some([sennew, sens]);

                    // set prescribed Ue-Hk combination
                    sys.vs2[3][0] = 0.0;
                    sys.vs2[3][1] = s2.hk_t * hkref;
                    sys.vs2[3][2] = s2.hk_d * hkref;
                    sys.vs2[3][3] = (s2.hk_u * hkref + sens / ueref) * s2.u_uei;
                    sys.vsrez[3] = -(hkref * hkref) * (s2.hk / hkref - 1.0) - sens * (s2.u / ueref - 1.0);
                }

                rec.ueref = ueref;
                rec.hkref = hkref;
                rec.residual = sys.vsrez;
                rec.vs2 = sys.vs2;

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

                // determine max changes and underrelax if necessary
                // (added Ue clamp   MD  3 Apr 03)
                dmax = (r[1] / thi).abs().max((r[2] / dsi).abs()).max((r[3] / uei).abs());
                if ibl >= st.itran[is] {
                    dmax = dmax.max((r[0] / (10.0 * cti)).abs());
                }
                let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

                // update as usual
                if ibl < st.itran[is] {
                    ami += rlx * r[0];
                }
                if ibl >= st.itran[is] {
                    cti += rlx * r[0];
                }
                thi += rlx * r[1];
                dsi += rlx * r[2];
                uei += rlx * r[3];

                // eliminate absurd transients
                if ibl >= st.itran[is] {
                    cti = cti.min(0.30);
                    cti = cti.max(0.0000001);
                }
                let hklim = if ibl <= st.iblte[is] { 1.02 } else { 1.00005 };
                let msq = uei * uei * params.h_stagnation_inv
                    / (params.gamma_gas_m1 * (1.0 - 0.5 * uei * uei * params.h_stagnation_inv));
                let mut dsw = dsi - dswaki;
                dslim(&mut dsw, thi, uei, msq, hklim);
                dsi = dsw + dswaki;

                rec.dmax = dmax;
                rec.rlx = rlx;
                rec.updated = [cti, thi, dsi, uei, ami];
                if dmax <= DEPS {
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
                // 'MRCHDU: Convergence failed at IBL side IS Res = DMAX'
                if let Some(t) = trace.as_mut() {
                    t.failed.push((ibl, is, dmax));
                }
                // the current unconverged solution might still be reasonable...
                if dmax > 0.1 {
                    // the current solution is garbage --> extrapolate values instead
                    if ibl > 3 {
                        if ibl <= st.iblte[is] {
                            thi = st.thet[is][ibm] * (st.xssi[is][ibl] / st.xssi[is][ibm]).powf(0.5);
                            dsi = st.dstr[is][ibm] * (st.xssi[is][ibl] / st.xssi[is][ibm]).powf(0.5);
                            uei = st.uedg[is][ibm];
                        } else if ibl == st.iblte[is] + 1 {
                            cti = cte;
                            thi = tte;
                            dsi = dte;
                            uei = st.uedg[is][ibm];
                        } else {
                            thi = st.thet[is][ibm];
                            let ratlen = (st.xssi[is][ibl] - st.xssi[is][ibm]) / (10.0 * st.dstr[is][ibm]);
                            dsi = (st.dstr[is][ibm] + thi * ratlen) / (1.0 + ratlen);
                            uei = st.uedg[is][ibm];
                        }
                        if ibl == st.itran[is] {
                            cti = 0.05;
                        }
                        if ibl > st.itran[is] {
                            cti = st.ctau[is][ibm];
                        }
                    }
                }
                // label 109
                s2.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, params);
                s2.blkin(params);
                // check for transition and set appropriate flags and things
                if !simi && !turb {
                    match trchek(&s1, &s2, s1.ampl, amcrit, xiforc, params) {
                        TransitionResult::NoTransition { ampl2 } => {
                            ami = ampl2;
                            tran = false;
                            trloc.xt = s2.x;
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
                // set all other extrapolated values for current station — XFOIL calls BLVAR in
                // this order, each call clamping HK2 in place, so the sequence is kept.
                // (BLMID only sets the interval CFM, which nothing reads after this point.)
                if ibl < st.itran[is] {
                    s2.blvar(FlowRegime::Laminar, params);
                }
                if ibl >= st.itran[is] {
                    s2.blvar(FlowRegime::Turbulent, params);
                }
                if wake {
                    s2.blvar(FlowRegime::Wake, params);
                }
            }

            // label 110: pick up here after the Newton iterations
            sens = sennew;

            // store primary variables
            st.ctau[is][ibl] = if ibl < st.itran[is] { ami } else { cti };
            st.thet[is][ibl] = thi;
            st.dstr[is][ibl] = dsi;
            st.uedg[is][ibl] = uei;
            st.mass[is][ibl] = dsi * uei;
            st.tau[is][ibl] = 0.5 * s2.r * s2.u * s2.u * s2.cf;
            st.dis[is][ibl] = s2.r * s2.u * s2.u * s2.u * s2.di * s2.hs * 0.5;
            st.ctq[is][ibl] = s2.cq;
            st.delt[is][ibl] = s2.de;
            st.tstr[is][ibl] = s2.hs * s2.theta;

            // set "1" variables to "2" variables for next streamwise station
            s2.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, params);
            s2.blkin(params);
            s1 = s2.clone();

            // turbulent intervals will follow transition interval or TE
            if tran || ibl == st.iblte[is] {
                turb = true;
                // save transition location
                st.tforce[is] = trforc;
                st.xssitr[is] = trloc.xt;
            }
            tran = false;
        }
    }
    st.com1 = s1;
    st.com2 = s2;
    st.trloc = trloc;
}
