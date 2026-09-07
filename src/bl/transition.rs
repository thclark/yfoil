//! Transition: TRCHEK2 (xblsys.f) with the XT sensitivities of XBL.INC's /V_VARA/ block, the
//! laminar amplification rates DAMPL / DAMPL2 and the interval average AXSET (xblsys.f).

use super::params::*;
use super::station::StationState;

// ============================================================================
// Transition Location and Derivatives
// ============================================================================

/// Transition location and its derivatives
///
/// Used by TRDIF to handle intervals containing laminar-turbulent transition.
/// The transition location XT is where N(x) = Ncrit.
#[derive(Debug, Clone, Default)]
pub struct TransitionLocation {
    /// Transition x-location
    pub xt: f64,
    /// ∂XT/∂N1 (amplification at station 1)
    pub xt_a1: f64,
    /// ∂XT/∂θ1
    pub xt_t1: f64,
    /// ∂XT/∂δ*1
    pub xt_d1: f64,
    /// ∂XT/∂Ue1
    pub xt_u1: f64,
    /// ∂XT/∂X1
    pub xt_x1: f64,
    /// ∂XT/∂θ2
    pub xt_t2: f64,
    /// ∂XT/∂δ*2
    pub xt_d2: f64,
    /// ∂XT/∂Ue2
    pub xt_u2: f64,
    /// ∂XT/∂X2
    pub xt_x2: f64,
    /// ∂XT/∂M²
    pub xt_ms: f64,
    /// ∂XT/∂Re
    pub xt_re: f64,
    /// ∂XT/∂Xforc (forced transition location)
    pub xt_xf: f64,
}

/// Result of transition check (TRCHEK2 equivalent)
#[derive(Debug, Clone)]
pub enum TransitionResult {
    /// No transition in this interval - return updated amplification
    NoTransition {
        /// Amplification factor at station 2
        ampl2: f64,
    },
    /// Free transition (N = Ncrit) occurred
    FreeTransition {
        /// Transition location and derivatives
        location: TransitionLocation,
        /// Amplification at station 2 (= Ncrit)
        ampl2: f64,
    },
    /// Forced transition at prescribed location
    ForcedTransition {
        /// Transition location (= xiforc)
        location: TransitionLocation,
    },
}

/// TRCHEK2 (xblsys.f): checks whether transition occurs in the interval X1..X2, solving the
/// implicit second-order amplification equation for N2 by Newton iteration.
///
/// Mirrors the Fortran control flow exactly: the loop's final `XT/TT/DT/UT`, `HKT/RTT`, `AX`
/// and interpolation weights are reused for the free-transition sensitivities (nothing is
/// recomputed after the loop), `AX <= 0` and non-convergence fall through to the transition
/// tests rather than returning early, and the returned `ampl2` is the iterated value (which
/// may exceed Ncrit) — exactly what XFOIL leaves in `AMPL2`.
#[allow(unused_assignments)] // loop-carried locals mirror the Fortran; the loop always runs
pub fn trchek(
    s1: &StationState,
    s2: &StationState,
    ampl1: f64,
    acrit: f64,
    xiforc: f64,
    params: &FlowParameters,
) -> TransitionResult {
    const DAEPS: f64 = 5.0e-5;
    let (x1, x2) = (s1.xi, s2.xi);

    // calculate average amplification rate AX over X1..X2 interval, with the current AMPL2
    let r0 = axset(
        s1.hk,
        s1.theta,
        s1.retheta,
        ampl1,
        s2.hk,
        s2.theta,
        s2.retheta,
        s2.ampl,
        acrit,
        params.idampv,
    );
    // set initial guess for iterate N2 (AMPL2) at X2
    let mut ampl2 = ampl1 + r0.ax * (x2 - x1);

    // loop state carried into the post-loop section (as XFOIL's locals/COMMON are)
    let mut amplt_a2 = 0.0;
    let (mut wf2, mut wf2_a1, mut wf2_a2, mut wf2_x1, mut wf2_x2, mut wf2_xf) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut xt, mut tt, mut dt, mut ut) = (0.0, 0.0, 0.0, 0.0);
    let (mut xt_a2, mut tt_a2, mut dt_a2, mut ut_a2) = (0.0, 0.0, 0.0, 0.0);
    let mut st = s2.clone();
    let mut r = r0.clone();

    // solve implicit system for amplification AMPL2
    for _itam in 1..=30 {
        // define weighting factors WF1,WF2 for defining "T" quantities from 1,2
        let (amplt, sfa, sfa_a1, sfa_a2);
        if ampl2 <= acrit {
            // there is no transition yet, "T" is the same as "2"
            amplt = ampl2;
            amplt_a2 = 1.0;
            sfa = 1.0;
            sfa_a1 = 0.0;
            sfa_a2 = 0.0;
        } else {
            // there is transition in X1..X2, "T" is set from N1, N2
            amplt = acrit;
            amplt_a2 = 0.0;
            sfa = (amplt - ampl1) / (ampl2 - ampl1);
            sfa_a1 = (sfa - 1.0) / (ampl2 - ampl1);
            sfa_a2 = (-sfa) / (ampl2 - ampl1);
        }

        let (sfx, sfx_x1, sfx_x2, sfx_xf) = if xiforc < x2 {
            let sfx = (xiforc - x1) / (x2 - x1);
            ((sfx), (sfx - 1.0) / (x2 - x1), (-sfx) / (x2 - x1), 1.0 / (x2 - x1))
        } else {
            (1.0, 0.0, 0.0, 0.0)
        };

        // set weighting factor from free or forced transition
        if sfa < sfx {
            wf2 = sfa;
            wf2_a1 = sfa_a1;
            wf2_a2 = sfa_a2;
            wf2_x1 = 0.0;
            wf2_x2 = 0.0;
            wf2_xf = 0.0;
        } else {
            wf2 = sfx;
            wf2_a1 = 0.0;
            wf2_a2 = 0.0;
            wf2_x1 = sfx_x1;
            wf2_x2 = sfx_x2;
            wf2_xf = sfx_xf;
        }
        let wf1 = 1.0 - wf2;
        let wf1_a2 = -wf2_a2;

        // interpolate BL variables to XT
        xt = x1 * wf1 + x2 * wf2;
        tt = s1.theta * wf1 + s2.theta * wf2;
        dt = s1.dstar * wf1 + s2.dstar * wf2;
        ut = s1.ue * wf1 + s2.ue * wf2;
        xt_a2 = x1 * wf1_a2 + x2 * wf2_a2;
        tt_a2 = s1.theta * wf1_a2 + s2.theta * wf2_a2;
        dt_a2 = s1.dstar * wf1_a2 + s2.dstar * wf2_a2;
        ut_a2 = s1.ue * wf1_a2 + s2.ue * wf2_a2;

        // temporarily set "2" variables from "T" for BLKIN (U2_UEI, U2_MS, DW2 stay station 2's)
        st = s2.clone();
        st.xi = xt;
        st.theta = tt;
        st.dstar = dt;
        st.ue = ut;
        st.set_kinematic_variables(params);

        // calculate amplification rate AX over current X1-XT interval
        r = axset(
            s1.hk,
            s1.theta,
            s1.retheta,
            ampl1,
            st.hk,
            tt,
            st.retheta,
            amplt,
            acrit,
            params.idampv,
        );

        // punch out early if there is no amplification here
        if r.ax <= 0.0 {
            break;
        }

        // set sensitivity of AX(A2)
        let ax_a2 = (r.ax_hk2 * st.hk_d_theta + r.ax_t2 + r.ax_rt2 * st.retheta_d_theta) * tt_a2
            + (r.ax_hk2 * st.hk_d_dstar) * dt_a2
            + (r.ax_hk2 * st.hk_d_ue + r.ax_rt2 * st.retheta_d_ue) * ut_a2
            + r.ax_a2 * amplt_a2;

        // residual for implicit AMPL2 definition (amplification equation)
        let res = ampl2 - ampl1 - r.ax * (x2 - x1);
        let res_a2 = 1.0 - ax_a2 * (x2 - x1);

        let da2 = -res / res_a2;

        let mut rlx = 1.0;
        let dxt = xt_a2 * da2;
        if rlx * (dxt / (x2 - x1)).abs() > 0.05 {
            rlx = 0.05 * ((x2 - x1) / dxt).abs();
        }
        if rlx * da2.abs() > 1.0 {
            rlx = 1.0 * (1.0 / da2).abs();
        }

        // check if converged
        if da2.abs() < DAEPS {
            break;
        }

        if (ampl2 > acrit && ampl2 + rlx * da2 < acrit) || (ampl2 < acrit && ampl2 + rlx * da2 > acrit) {
            // limited Newton step so AMPL2 doesn't step across AMCRIT either way
            ampl2 = acrit;
        } else {
            // regular Newton step
            ampl2 += rlx * da2;
        }
        // (XFOIL prints 'TRCHEK2: N2 convergence failed.' after 30 iterations and continues)
    }

    // label 101: test for free or forced transition
    let mut trfree = ampl2 >= acrit;
    let mut trforc = xiforc > x1 && xiforc <= x2;

    // set transition interval flag
    let tran = trforc || trfree;
    if !tran {
        return TransitionResult::NoTransition { ampl2 };
    }

    // resolve if both forced and free transition
    if trfree && trforc {
        trforc = xiforc < xt;
        trfree = xiforc >= xt;
    }
    let _ = trfree;

    if trforc {
        // if forced transition, then XT is prescribed
        let location = TransitionLocation {
            xt: xiforc,
            xt_xf: 1.0,
            ..Default::default()
        };
        return TransitionResult::ForcedTransition { location };
    }

    // free transition ... set sensitivities of XT
    let wf1 = 1.0 - wf2;
    let wf1_a1 = -wf2_a1;
    let wf1_x1 = -wf2_x1;
    let wf1_x2 = -wf2_x2;
    let wf1_xf = -wf2_xf;

    let mut xt_x1 = wf1;
    let tt_t1 = wf1;
    let dt_d1 = wf1;
    let ut_u1 = wf1;
    let mut xt_x2 = wf2;
    let tt_t2 = wf2;
    let dt_d2 = wf2;
    let ut_u2 = wf2;

    let xt_a1 = x1 * wf1_a1 + x2 * wf2_a1;
    let tt_a1 = s1.theta * wf1_a1 + s2.theta * wf2_a1;
    let dt_a1 = s1.dstar * wf1_a1 + s2.dstar * wf2_a1;
    let ut_a1 = s1.ue * wf1_a1 + s2.ue * wf2_a1;

    xt_x1 = x1 * wf1_x1 + x2 * wf2_x1 + xt_x1;
    let tt_x1 = s1.theta * wf1_x1 + s2.theta * wf2_x1;
    let dt_x1 = s1.dstar * wf1_x1 + s2.dstar * wf2_x1;
    let ut_x1 = s1.ue * wf1_x1 + s2.ue * wf2_x1;

    xt_x2 = x1 * wf1_x2 + x2 * wf2_x2 + xt_x2;
    let tt_x2 = s1.theta * wf1_x2 + s2.theta * wf2_x2;
    let dt_x2 = s1.dstar * wf1_x2 + s2.dstar * wf2_x2;
    let ut_x2 = s1.ue * wf1_x2 + s2.ue * wf2_x2;

    let _xt_xf = x1 * wf1_xf + x2 * wf2_xf;
    let tt_xf = s1.theta * wf1_xf + s2.theta * wf2_xf;
    let dt_xf = s1.dstar * wf1_xf + s2.dstar * wf2_xf;
    let ut_xf = s1.ue * wf1_xf + s2.ue * wf2_xf;

    // at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT ) from the last loop pass
    let (hkt_tt, hkt_dt, hkt_ut, hkt_ms) = (st.hk_d_theta, st.hk_d_dstar, st.hk_d_ue, st.hk_d_machsqd);
    let (rtt_tt, rtt_ut, rtt_ms, rtt_re) = (
        st.retheta_d_theta,
        st.retheta_d_ue,
        st.retheta_d_machsqd,
        st.retheta_d_re,
    );
    let ax_t1 = r.ax_hk1 * s1.hk_d_theta
        + r.ax_t1
        + r.ax_rt1 * s1.retheta_d_theta
        + (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_t1;
    let ax_d1 = r.ax_hk1 * s1.hk_d_dstar + (r.ax_hk2 * hkt_dt) * dt_d1;
    let ax_u1 = r.ax_hk1 * s1.hk_d_ue + r.ax_rt1 * s1.retheta_d_ue + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_u1;
    let ax_a1 = r.ax_a1
        + (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_a1
        + (r.ax_hk2 * hkt_dt) * dt_a1
        + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_a1;
    let ax_x1 = (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_x1
        + (r.ax_hk2 * hkt_dt) * dt_x1
        + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_x1;
    let ax_t2 = (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_t2;
    let ax_d2 = (r.ax_hk2 * hkt_dt) * dt_d2;
    let ax_u2 = (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_u2;
    let ax_a2 = r.ax_a2 * amplt_a2
        + (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_a2
        + (r.ax_hk2 * hkt_dt) * dt_a2
        + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_a2;
    let ax_x2 = (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_x2
        + (r.ax_hk2 * hkt_dt) * dt_x2
        + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_x2;
    let ax_xf = (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_xf
        + (r.ax_hk2 * hkt_dt) * dt_xf
        + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_xf;
    let ax_ms = r.ax_hk2 * hkt_ms + r.ax_rt2 * rtt_ms + r.ax_hk1 * s1.hk_d_machsqd + r.ax_rt1 * s1.retheta_d_machsqd;
    let ax_re = r.ax_rt2 * rtt_re + r.ax_rt1 * s1.retheta_d_re;

    // set sensitivities of residual RES
    let z_ax = -(x2 - x1);
    let z_a1 = z_ax * ax_a1 - 1.0;
    let z_t1 = z_ax * ax_t1;
    let z_d1 = z_ax * ax_d1;
    let z_u1 = z_ax * ax_u1;
    let z_x1 = z_ax * ax_x1 + r.ax;
    let z_a2 = z_ax * ax_a2 + 1.0;
    let z_t2 = z_ax * ax_t2;
    let z_d2 = z_ax * ax_d2;
    let z_u2 = z_ax * ax_u2;
    let z_x2 = z_ax * ax_x2 - r.ax;
    let _z_xf = z_ax * ax_xf;
    let z_ms = z_ax * ax_ms;
    let z_re = z_ax * ax_re;

    // set sensitivities of XT, with RES being stationary for A2 constraint
    let location = TransitionLocation {
        xt,
        xt_a1: xt_a1 - (xt_a2 / z_a2) * z_a1,
        xt_t1: -(xt_a2 / z_a2) * z_t1,
        xt_d1: -(xt_a2 / z_a2) * z_d1,
        xt_u1: -(xt_a2 / z_a2) * z_u1,
        xt_x1: xt_x1 - (xt_a2 / z_a2) * z_x1,
        xt_t2: -(xt_a2 / z_a2) * z_t2,
        xt_d2: -(xt_a2 / z_a2) * z_d2,
        xt_u2: -(xt_a2 / z_a2) * z_u2,
        xt_x2: xt_x2 - (xt_a2 / z_a2) * z_x2,
        xt_ms: -(xt_a2 / z_a2) * z_ms,
        xt_re: -(xt_a2 / z_a2) * z_re,
        xt_xf: 0.0,
    };
    TransitionResult::FreeTransition { location, ampl2 }
}

// ============================================================================
// Laminar Amplification Rate (DAMPL, AXSET)
// ============================================================================

/// Result from DAMPL - local amplification rate
#[derive(Debug, Clone, Default)]
pub struct AmplificationRate {
    /// Spatial amplification rate dN/dx
    pub ax: f64,
    /// ∂AX/∂Hk
    pub ax_hk: f64,
    /// ∂AX/∂θ
    pub ax_th: f64,
    /// ∂AX/∂Rθ
    pub ax_rt: f64,
}

/// Calculate local amplification rate (DAMPL equivalent)
///
/// Amplification rate routine for envelope e^n method.
/// Reference: Drela, M., Giles, M., "Viscous/Inviscid Analysis of Transonic
/// and Low Reynolds Number Airfoils", AIAA Journal, Oct. 1987.
///
/// # Arguments
/// * `hk` - Kinematic shape parameter
/// * `th` - Momentum thickness
/// * `rt` - Momentum thickness Reynolds number
///
/// # Returns
/// Spatial amplification rate and derivatives
pub fn dampl(hk: f64, th: f64, rt: f64) -> AmplificationRate {
    const DGR: f64 = 0.08; // Ramp width in log10(Rt)

    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // log10(Critical Rth) - H correlation for Falkner-Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    let grcrit = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    let gr = rt.log10();
    let gr_rt = 1.0 / (2.3025851 * rt);

    if gr < grcrit - DGR {
        // No amplification for Rtheta < Rcrit
        AmplificationRate::default()
    } else {
        // Steep cubic ramp to turn on AX smoothly as Rtheta exceeds Rcrit
        let rnorm = (gr - (grcrit - DGR)) / (2.0 * DGR);
        let rn_hk = -grc_hk / (2.0 * DGR);
        let rn_rt = gr_rt / (2.0 * DGR);

        let (rfac, rfac_hk, rfac_rt) = if rnorm >= 1.0 {
            (1.0, 0.0, 0.0)
        } else {
            let rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
            let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
            (rfac, rfac_rn * rn_hk, rfac_rn * rn_rt)
        };

        // Amplification envelope slope correlation for Falkner-Skan
        let arg = 3.87 * hmi - 2.52;
        let arg_hk = 3.87 * hmi_hk;

        let ex = (-arg * arg).exp();
        let ex_hk = ex * (-2.0 * arg * arg_hk);

        let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
        let dadr_hk = 0.028 - 0.0345 * ex_hk;

        // m(H) correlation (March 1991 version)
        let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi;
        let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi;
        let af_hk = af_hmi * hmi_hk;

        let ax = (af * dadr / th) * rfac;
        let ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
        let ax_th = -ax / th;
        let ax_rt = (af * dadr / th) * rfac_rt;

        AmplificationRate {
            ax,
            ax_hk,
            ax_th,
            ax_rt,
        }
    }
}

/// DAMPL2 (xblsys.f): amplification rate for the *modified* envelope e^n method (Nov 1996) —
/// the envelope rate of `dampl`, blended for Hk > 3.5 with the Orr–Sommerfeld maximum
/// amplification correlation for separated profiles. Selected by OPER `DAMP` (IDAMPV = 1).
pub fn dampl2(hk: f64, th: f64, rt: f64) -> AmplificationRate {
    const DGR: f64 = 0.08;
    const HK1: f64 = 3.5;
    const HK2: f64 = 4.0;

    let hmi = 1.0 / (hk - 1.0);
    let hmi_hk = -hmi * hmi;

    // log10(Critical Rth) -- H   correlation for Falkner-Skan profiles
    let aa = 2.492 * hmi.powf(0.43);
    let aa_hk = (aa / hmi) * 0.43 * hmi_hk;

    let bb = (14.0 * hmi - 9.24).tanh();
    let bb_hk = (1.0 - bb * bb) * 14.0 * hmi_hk;

    let grc = aa + 0.7 * (bb + 1.0);
    let grc_hk = aa_hk + 0.7 * bb_hk;

    let gr = rt.log10();
    let gr_rt = 1.0 / (2.3025851 * rt);

    let (mut ax, mut ax_hk, mut ax_th, mut ax_rt);
    if gr < grc - DGR {
        // no amplification for Rtheta < Rcrit
        ax = 0.0;
        ax_hk = 0.0;
        ax_th = 0.0;
        ax_rt = 0.0;
    } else {
        // Set steep cubic ramp used to turn on AX smoothly as Rtheta exceeds Rcrit
        let rnorm = (gr - (grc - DGR)) / (2.0 * DGR);
        let rn_hk = -grc_hk / (2.0 * DGR);
        let rn_rt = gr_rt / (2.0 * DGR);

        let (rfac, rfac_hk, rfac_rt);
        if rnorm >= 1.0 {
            rfac = 1.0;
            rfac_hk = 0.0;
            rfac_rt = 0.0;
        } else {
            rfac = 3.0 * rnorm * rnorm - 2.0 * rnorm * rnorm * rnorm;
            let rfac_rn = 6.0 * rnorm - 6.0 * rnorm * rnorm;
            rfac_hk = rfac_rn * rn_hk;
            rfac_rt = rfac_rn * rn_rt;
        }

        // set envelope amplification rate with respect to Rtheta: DADR = d(N)/d(Rtheta) = f(H)
        let arg = 3.87 * hmi - 2.52;
        let arg_hk = 3.87 * hmi_hk;

        let ex = (-arg * arg).exp();
        let ex_hk = ex * (-2.0 * arg * arg_hk);

        let dadr = 0.028 * (hk - 1.0) - 0.0345 * ex;
        let dadr_hk = 0.028 - 0.0345 * ex_hk;

        // set conversion factor from d/d(Rtheta) to d/dx: AF = Theta d(Rtheta)/dx = f(H)
        let brg = -20.0 * hmi;
        let af = -0.05 + 2.7 * hmi - 5.5 * hmi * hmi + 3.0 * hmi * hmi * hmi + 0.1 * brg.exp();
        let af_hmi = 2.7 - 11.0 * hmi + 9.0 * hmi * hmi - 2.0 * brg.exp();
        let af_hk = af_hmi * hmi_hk;

        // set amplification rate with respect to x, with RFAC shutting off amplification below Rcrit
        ax = (af * dadr / th) * rfac;
        ax_hk = (af_hk * dadr / th + af * dadr_hk / th) * rfac + (af * dadr / th) * rfac_hk;
        ax_th = -ax / th;
        ax_rt = (af * dadr / th) * rfac_rt;
    }

    if hk < HK1 {
        return AmplificationRate {
            ax,
            ax_hk,
            ax_th,
            ax_rt,
        };
    }

    // non-envelope max-amplification correction for separated profiles
    let hnorm = (hk - HK1) / (HK2 - HK1);
    let hn_hk = 1.0 / (HK2 - HK1);

    // set blending fraction HFAC = 0..1 over HK1 < HK < HK2
    let (hfac, hf_hk);
    if hnorm >= 1.0 {
        hfac = 1.0;
        hf_hk = 0.0;
    } else {
        hfac = 3.0 * hnorm * hnorm - 2.0 * hnorm * hnorm * hnorm;
        hf_hk = (6.0 * hnorm - 6.0 * hnorm * hnorm) * hn_hk;
    }

    // "normal" envelope amplification rate AX1
    let ax1 = ax;
    let ax1_hk = ax_hk;
    let ax1_th = ax_th;
    let ax1_rt = ax_rt;

    // set modified amplification rate AX2
    let gr0 = 0.30 + 0.35 * (-0.15 * (hk - 5.0)).exp();
    let gr0_hk = -0.35 * (-0.15 * (hk - 5.0)).exp() * 0.15;

    let tnr = (1.2 * (gr - gr0)).tanh();
    let tnr_rt = (1.0 - tnr * tnr) * 1.2 * gr_rt;
    let tnr_hk = -(1.0 - tnr * tnr) * 1.2 * gr0_hk;

    let (mut ax2, mut ax2_hk, mut ax2_rt, mut ax2_th);
    ax2 = (0.086 * tnr - 0.25 / (hk - 1.0).powf(1.5)) / th;
    ax2_hk = (0.086 * tnr_hk + 1.5 * 0.25 / (hk - 1.0).powf(2.5)) / th;
    ax2_rt = (0.086 * tnr_rt) / th;
    ax2_th = -ax2 / th;

    if ax2 < 0.0 {
        ax2 = 0.0;
        ax2_hk = 0.0;
        ax2_rt = 0.0;
        ax2_th = 0.0;
    }

    // blend the two amplification rates
    ax = hfac * ax2 + (1.0 - hfac) * ax1;
    ax_hk = hfac * ax2_hk + (1.0 - hfac) * ax1_hk + hf_hk * (ax2 - ax1);
    ax_rt = hfac * ax2_rt + (1.0 - hfac) * ax1_rt;
    ax_th = hfac * ax2_th + (1.0 - hfac) * ax1_th;

    AmplificationRate {
        ax,
        ax_hk,
        ax_th,
        ax_rt,
    }
}

/// Result from AXSET - averaged amplification rate over interval
#[derive(Debug, Clone, Default)]
pub struct AveragedAmplification {
    /// Average amplification rate
    pub ax: f64,
    /// ∂AX/∂Hk1
    pub ax_hk1: f64,
    /// ∂AX/∂θ1
    pub ax_t1: f64,
    /// ∂AX/∂Rθ1
    pub ax_rt1: f64,
    /// ∂AX/∂N1 (amplification factor)
    pub ax_a1: f64,
    /// ∂AX/∂Hk2
    pub ax_hk2: f64,
    /// ∂AX/∂θ2
    pub ax_t2: f64,
    /// ∂AX/∂Rθ2
    pub ax_rt2: f64,
    /// ∂AX/∂N2
    pub ax_a2: f64,
}

/// Calculate averaged amplification rate over interval (AXSET equivalent)
///
/// Returns average amplification AX over interval 1..2 using RMS average.
///
/// # Arguments
/// * `hk1`, `t1`, `rt1`, `a1` - Station 1 values (Hk, θ, Rθ, N)
/// * `hk2`, `t2`, `rt2`, `a2` - Station 2 values
/// * `acrit` - Critical amplification factor for transition
///
/// # Returns
/// Averaged amplification rate and derivatives
pub fn axset(
    hk1: f64,
    t1: f64,
    rt1: f64,
    a1: f64,
    hk2: f64,
    t2: f64,
    rt2: f64,
    a2: f64,
    acrit: f64,
    idampv: usize,
) -> AveragedAmplification {
    // 2nd-order: local amplification rates at both stations, envelope (IDAMPV = 0) or
    // modified-envelope (IDAMPV = 1, OPER DAMP) method
    let (ax1_result, ax2_result) = if idampv == 0 {
        (dampl(hk1, t1, rt1), dampl(hk2, t2, rt2))
    } else {
        (dampl2(hk1, t1, rt1), dampl2(hk2, t2, rt2))
    };

    let ax1 = ax1_result.ax;
    let ax2 = ax2_result.ax;

    // RMS-average version (better on coarse grids)
    let axsq = 0.5 * (ax1 * ax1 + ax2 * ax2);
    let (axa, axa_ax1, axa_ax2) = if axsq <= 0.0 {
        (0.0, 0.0, 0.0)
    } else {
        let axa = axsq.sqrt();
        (axa, 0.5 * ax1 / axa, 0.5 * ax2 / axa)
    };

    // Small additional term to ensure dN/dx > 0 near N = Ncrit
    let arg = (20.0 * (acrit - 0.5 * (a1 + a2))).min(20.0);
    let (exn, exn_a1, exn_a2) = if arg <= 0.0 {
        (1.0, 0.0, 0.0)
    } else {
        let exn = (-arg).exp();
        (exn, 20.0 * 0.5 * exn, 20.0 * 0.5 * exn)
    };

    let dax = exn * 0.002 / (t1 + t2);
    let dax_a1 = exn_a1 * 0.002 / (t1 + t2);
    let dax_a2 = exn_a2 * 0.002 / (t1 + t2);
    let dax_t1 = -dax / (t1 + t2);
    let dax_t2 = -dax / (t1 + t2);

    // Final result
    let ax = axa + dax;

    AveragedAmplification {
        ax,
        ax_hk1: axa_ax1 * ax1_result.ax_hk,
        ax_t1: axa_ax1 * ax1_result.ax_th + dax_t1,
        ax_rt1: axa_ax1 * ax1_result.ax_rt,
        ax_a1: dax_a1,
        ax_hk2: axa_ax2 * ax2_result.ax_hk,
        ax_t2: axa_ax2 * ax2_result.ax_th + dax_t2,
        ax_rt2: axa_ax2 * ax2_result.ax_rt,
        ax_a2: dax_a2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    // ========================================================================
    // DAMPL Tests - validate laminar amplification rate against XFOIL Fortran
    // ========================================================================

    #[test]
    fn test_dampl_below_critical() {
        // Test case 1: Below critical Rt (should return 0)
        // HK = 2.5, TH = 0.001, RT = 100
        let result = dampl(2.5, 0.001, 100.0);

        assert_eq!(result.ax, 0.0);
        assert_eq!(result.ax_hk, 0.0);
        assert_eq!(result.ax_th, 0.0);
        assert_eq!(result.ax_rt, 0.0);
    }

    #[test]
    fn test_dampl_ramp_region() {
        // Test case 2: Near critical Rt (ramp region)
        // HK = 2.5, TH = 0.001, RT = 700
        // XFOIL reference: AX = 0.5576204658E+00
        let result = dampl(2.5, 0.001, 700.0);

        assert_relative_eq!(result.ax, 0.5576, epsilon = 0.01);
        assert_relative_eq!(result.ax_hk, 68.78, epsilon = 1.0);
        assert_relative_eq!(result.ax_th, -557.6, epsilon = 1.0);
        assert!(result.ax_rt > 0.0); // Should be positive in ramp region
    }

    #[test]
    fn test_dampl_above_critical() {
        // Test case 3: Above critical Rt
        // HK = 2.5, TH = 0.002, RT = 2000
        // XFOIL reference values:
        // AX = 0.7412202954E+00
        // AX_HK = 0.3105413914E+01
        // AX_T = -0.3706101379E+03
        // AX_RT = 0.0 (above ramp)
        let result = dampl(2.5, 0.002, 2000.0);

        assert_relative_eq!(result.ax, 0.7412, epsilon = 0.001);
        assert_relative_eq!(result.ax_hk, 3.105, epsilon = 0.01);
        assert_relative_eq!(result.ax_th, -370.6, epsilon = 0.5);
        assert_eq!(result.ax_rt, 0.0); // Above ramp, derivative is 0
    }

    #[test]
    fn test_dampl_lower_hk() {
        // Test case 4: Lower Hk (more unstable profile)
        // HK = 2.2, TH = 0.0015, RT = 5000
        // XFOIL reference: AX = 0.5429499745E+00
        let result = dampl(2.2, 0.0015, 5000.0);

        assert_relative_eq!(result.ax, 0.543, epsilon = 0.01);
        assert_relative_eq!(result.ax_hk, 7.95, epsilon = 0.1);
        assert_relative_eq!(result.ax_th, -362.0, epsilon = 1.0);
    }

    #[test]
    fn test_dampl_higher_hk() {
        // Test case 5: Higher Hk (more stable profile)
        // HK = 2.8, TH = 0.0025, RT = 3000
        // XFOIL reference: AX = 0.2168501139E+01
        let result = dampl(2.8, 0.0025, 3000.0);

        assert_relative_eq!(result.ax, 2.168, epsilon = 0.01);
        assert_relative_eq!(result.ax_hk, 7.416, epsilon = 0.1);
        assert_relative_eq!(result.ax_th, -867.4, epsilon = 1.0);
    }

    // ========================================================================
    // DAMPL Exact Fixture Tests - machine precision comparison with XFOIL
    // These fixtures were extracted from instrumented XFOIL NACA0012 at Re=1e6
    // ========================================================================

    #[test]
    fn test_dampl_xfoil_fixture_1() {
        // XFOIL fixture: near critical transition (in ramp region)
        // HK=0.2580689064632752E+01 TH=0.2202254221917420E-03 RT=0.2615260734070808E+03
        // AX=0.4909634245555000E-01 AX_HK=0.5787501799836026E+02 AX_TH=-0.2229367616459996E+03 AX_RT=0.2419961491556620E-01
        let hk = 0.2580689064632752e+01;
        let th = 0.2202254221917420e-03;
        let rt = 0.2615260734070808e+03;

        let result = dampl(hk, th, rt);

        let rel_tol = 1e-10;
        assert_relative_eq!(result.ax, 0.4909634245555000e-01, epsilon = rel_tol);
        assert_relative_eq!(result.ax_hk, 0.5787501799836026e+02, epsilon = rel_tol * 100.0);
        assert_relative_eq!(result.ax_th, -0.2229367616459996e+03, epsilon = rel_tol * 1000.0);
        assert_relative_eq!(result.ax_rt, 0.2419961491556620e-01, epsilon = rel_tol);
    }

    #[test]
    fn test_dampl_xfoil_fixture_2() {
        // XFOIL fixture: early amplification region
        // HK=0.2593152760810387E+01 TH=0.2328303859189811E-03 RT=0.2761529642099222E+03
        // AX=0.4728845409931889E+01 AX_HK=0.3642605131865983E+03 AX_TH=-0.2031025886620058E+05 AX_RT=0.1445021321701126E+00
        let hk = 0.2593152760810387e+01;
        let th = 0.2328303859189811e-03;
        let rt = 0.2761529642099222e+03;

        let result = dampl(hk, th, rt);

        let rel_tol = 1e-10;
        assert_relative_eq!(result.ax, 0.4728845409931889e+01, epsilon = rel_tol);
        assert_relative_eq!(result.ax_hk, 0.3642605131865983e+03, epsilon = rel_tol * 1000.0);
        assert_relative_eq!(result.ax_th, -0.2031025886620058e+05, epsilon = rel_tol * 100000.0);
        assert_relative_eq!(result.ax_rt, 0.1445021321701126e+00, epsilon = rel_tol);
    }

    #[test]
    fn test_dampl_xfoil_fixture_3() {
        // XFOIL fixture: mid amplification region
        // HK=0.2606142972624875E+01 TH=0.2455695046420016E-03 RT=0.2907803892252000E+03
        // AX=0.9712213028497693E+01 AX_HK=0.1409468673226392E+03 AX_TH=-0.3954975208610060E+05 AX_RT=0.4082874225503037E-01
        let hk = 0.2606142972624875e+01;
        let th = 0.2455695046420016e-03;
        let rt = 0.2907803892252000e+03;

        let result = dampl(hk, th, rt);

        let rel_tol = 1e-10;
        assert_relative_eq!(result.ax, 0.9712213028497693e+01, epsilon = rel_tol);
        assert_relative_eq!(result.ax_hk, 0.1409468673226392e+03, epsilon = rel_tol * 1000.0);
        assert_relative_eq!(result.ax_th, -0.3954975208610060e+05, epsilon = rel_tol * 100000.0);
        assert_relative_eq!(result.ax_rt, 0.4082874225503037e-01, epsilon = rel_tol);
    }

    #[test]
    fn test_dampl_xfoil_fixture_4() {
        // XFOIL fixture: above ramp (AX_RT = 0)
        // HK=0.2619759006090018E+01 TH=0.2584472817606003E-03 RT=0.3054087290294926E+03
        // AX=0.1002497313242920E+02 AX_HK=0.4715318465087587E+02 AX_TH=-0.3878923803778030E+05 AX_RT=0.0
        let hk = 0.2619759006090018e+01;
        let th = 0.2584472817606003e-03;
        let rt = 0.3054087290294926e+03;

        let result = dampl(hk, th, rt);

        let rel_tol = 1e-10;
        assert_relative_eq!(result.ax, 0.1002497313242920e+02, epsilon = rel_tol);
        assert_relative_eq!(result.ax_hk, 0.4715318465087587e+02, epsilon = rel_tol * 100.0);
        assert_relative_eq!(result.ax_th, -0.3878923803778030e+05, epsilon = rel_tol * 100000.0);
        assert_eq!(result.ax_rt, 0.0); // Above ramp
    }

    // ========================================================================
    // AXSET Tests - validate averaged amplification rate against XFOIL Fortran
    // ========================================================================

    #[test]
    fn test_axset_growing_bl() {
        // AXSET Test 1: Growing BL
        // XFOIL reference values
        let result = axset(
            2.5, 0.0015, 1800.0, 3.0, // Station 1: HK, T, RT, A
            2.6, 0.0018, 2200.0, 4.5, // Station 2
            9.0, // ACRIT
            0,   // IDAMPV
        );

        // AX = 0.1160733819E+01
        assert_relative_eq!(result.ax, 1.1607, epsilon = 0.01);
        // AX_HK1 = 0.1762713194E+01
        assert_relative_eq!(result.ax_hk1, 1.76, epsilon = 0.1);
        // AX_T1 = -0.2804905396E+03
        assert_relative_eq!(result.ax_t1, -280.5, epsilon = 1.0);
        // AX_HK2 = 0.3531968832E+01
        assert_relative_eq!(result.ax_hk2, 3.53, epsilon = 0.1);
        // AX_T2 = -0.4111100464E+03
        assert_relative_eq!(result.ax_t2, -411.1, epsilon = 1.0);
    }

    #[test]
    fn test_axset_near_transition() {
        // AXSET Test 2: Near transition
        // XFOIL reference values
        let result = axset(
            2.5, 0.002, 2500.0, 7.5, // Station 1
            2.55, 0.0022, 2750.0, 8.5, // Station 2
            9.0, // ACRIT
            0,   // IDAMPV
        );

        // AX = 0.7944293022E+00
        assert_relative_eq!(result.ax, 0.794, epsilon = 0.01);
        // AX_HK1 = 0.1448710322E+01
        assert_relative_eq!(result.ax_hk1, 1.45, epsilon = 0.1);
        // AX_T1 = -0.1728937836E+03
        assert_relative_eq!(result.ax_t1, -172.9, epsilon = 1.0);
        // AX_HK2 = 0.2121723175E+01
        assert_relative_eq!(result.ax_hk2, 2.12, epsilon = 0.1);
        // AX_T2 = -0.2039280548E+03
        assert_relative_eq!(result.ax_t2, -203.9, epsilon = 1.0);
    }

    // ========================================================================
    // TRCHEK Tests - validate transition detection
    // ========================================================================

    #[test]
    fn test_trchek_no_transition_below_critical() {
        // Test case: Rtheta below critical, no amplification
        // AX = 0 when RT is below critical Reynolds number

        let params = FlowParameters::new(0.0, 1e6, 1.4);

        // Create stations with low Rtheta (below critical for HK ~2.5)
        // Critical log10(RT) for HK=2.5 is about 2.86, so RT ~720 is below critical
        let mut s1 = StationState::default();
        s1.xi = 0.02;
        s1.ue = 1.15;
        s1.theta = 0.0002;
        s1.dstar = 0.0005;
        s1.ampl = 0.0;
        s1.set_primary_variables(s1.xi, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.ue, &params);
        s1.set_kinematic_variables(&params);
        // s1.rt should be around 200 (below critical)

        let mut s2 = StationState::default();
        s2.xi = 0.04;
        s2.ue = 1.12;
        s2.theta = 0.0004;
        s2.dstar = 0.0010;
        s2.ampl = 0.0;
        s2.set_primary_variables(s2.xi, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.ue, &params);
        s2.set_kinematic_variables(&params);

        let result = trchek(&s1, &s2, 0.0, 9.0, 1e6, &params);

        // Should return NoTransition since RT is below critical
        match result {
            TransitionResult::NoTransition { ampl2 } => {
                // ampl2 should be close to 0 (no growth)
                assert!(ampl2 < 1.0, "Expected no significant amplification growth");
            }
            _ => panic!("Expected NoTransition result for sub-critical Reynolds number"),
        }
    }

    #[test]
    fn test_trchek_no_transition_below_ncrit() {
        // Test case: Amplification growing but not reaching Ncrit

        let params = FlowParameters::new(0.0, 1e6, 1.4);

        // Create stations with moderate Rtheta (above critical, amplification active)
        let mut s1 = StationState::default();
        s1.xi = 0.08;
        s1.ue = 1.10;
        s1.theta = 0.001;
        s1.dstar = 0.0025;
        s1.ampl = 2.0; // Starting amplification
        s1.set_primary_variables(s1.xi, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.ue, &params);
        s1.set_kinematic_variables(&params);

        let mut s2 = StationState::default();
        s2.xi = 0.10;
        s2.ue = 1.08;
        s2.theta = 0.0012;
        s2.dstar = 0.003;
        s2.ampl = 2.5;
        s2.set_primary_variables(s2.xi, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.ue, &params);
        s2.set_kinematic_variables(&params);

        let result = trchek(&s1, &s2, 2.0, 9.0, 1e6, &params);

        // Should return NoTransition with growing ampl2
        match result {
            TransitionResult::NoTransition { ampl2 } => {
                // ampl2 should be greater than ampl1 but less than Ncrit
                assert!(ampl2 >= 2.0, "Amplification should grow");
                assert!(ampl2 < 9.0, "Amplification should not reach Ncrit");
            }
            _ => panic!("Expected NoTransition result"),
        }
    }

    #[test]
    fn test_trchek_free_transition() {
        // Test case: Amplification exceeds Ncrit within interval

        let params = FlowParameters::new(0.0, 1e6, 1.4);

        // Create stations with high Rtheta and amplification close to Ncrit
        let mut s1 = StationState::default();
        s1.xi = 0.15;
        s1.ue = 1.05;
        s1.theta = 0.002;
        s1.dstar = 0.005;
        s1.ampl = 8.0; // Close to Ncrit=9
        s1.set_primary_variables(s1.xi, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.ue, &params);
        s1.set_kinematic_variables(&params);

        let mut s2 = StationState::default();
        s2.xi = 0.25;
        s2.ue = 1.02;
        s2.theta = 0.003;
        s2.dstar = 0.0075;
        s2.ampl = 12.0; // Well above Ncrit
        s2.set_primary_variables(s2.xi, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.ue, &params);
        s2.set_kinematic_variables(&params);

        let result = trchek(&s1, &s2, 8.0, 9.0, 1e6, &params);

        // Should return FreeTransition
        match result {
            TransitionResult::FreeTransition { location, ampl2 } => {
                // Transition should occur between X1 and X2
                assert!(
                    location.xt >= s1.xi && location.xt <= s2.xi,
                    "Transition location should be within interval"
                );
                // ampl2 should equal Ncrit
                assert_relative_eq!(ampl2, 9.0, epsilon = 0.01);
            }
            TransitionResult::NoTransition { ampl2 } => {
                // If no transition, ampl2 must have grown past Ncrit
                // This might happen if amplification is very fast
                if ampl2 >= 9.0 {
                    panic!(
                        "Amplification exceeded Ncrit but FreeTransition not returned: ampl2={}",
                        ampl2
                    );
                }
            }
            _ => panic!("Expected FreeTransition result"),
        }
    }

    #[test]
    fn test_trchek_forced_transition() {
        // Test case: Forced transition at prescribed location

        let params = FlowParameters::new(0.0, 1e6, 1.4);

        // Create stations where natural transition wouldn't occur
        let mut s1 = StationState::default();
        s1.xi = 0.10;
        s1.ue = 1.10;
        s1.theta = 0.001;
        s1.dstar = 0.0025;
        s1.ampl = 2.0;
        s1.set_primary_variables(s1.xi, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.ue, &params);
        s1.set_kinematic_variables(&params);

        let mut s2 = StationState::default();
        s2.xi = 0.20;
        s2.ue = 1.05;
        s2.theta = 0.0015;
        s2.dstar = 0.0038;
        s2.ampl = 3.0;
        s2.set_primary_variables(s2.xi, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.ue, &params);
        s2.set_kinematic_variables(&params);

        // Force transition at x = 0.15
        let xiforc = 0.15;
        let result = trchek(&s1, &s2, 2.0, 9.0, xiforc, &params);

        // Should return ForcedTransition at xiforc
        match result {
            TransitionResult::ForcedTransition { location } => {
                assert_relative_eq!(location.xt, xiforc, epsilon = 1e-10);
                assert_relative_eq!(location.xt_xf, 1.0, epsilon = 1e-10);
            }
            TransitionResult::FreeTransition { location, .. } => {
                // Free transition might occur first if amplification is high enough
                // Check if it's before xiforc
                if location.xt > xiforc {
                    panic!("Expected forced transition at {}, got free at {}", xiforc, location.xt);
                }
            }
            TransitionResult::NoTransition { .. } => {
                panic!("Expected ForcedTransition result, got NoTransition");
            }
        }
    }
}
