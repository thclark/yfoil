//! BL Newton system data structures and solver
//!
//! This module implements XFOIL's coupled Newton system for the
//! viscous-inviscid iteration. The key insight is that XFOIL solves
//! ALL BL equations simultaneously with the DIJ coupling embedded
//! in the Jacobian matrix.
//!
//! ## Newton System Structure
//!
//! The system has 3 equations per BL station:
//! 1. Momentum integral equation
//! 2. Shape parameter (energy) equation
//! 3. Lag/shear stress equation
//!
//! The unknowns at each station are:
//! 1. dCtau (or dAmpl for laminar) - shear/amplification change
//! 2. dTheta - momentum thickness change
//! 3. dMass - mass defect change (where Mass = Ue * δ*)
//!
//! The Jacobian has special structure:
//! - Block tridiagonal (VA, VB) for local station coupling
//! - Dense VM matrix for global DIJ coupling
//!
//! ## Reference
//!
//! XFOIL source files:
//! - xbl.f: SETBL (Newton system assembly)
//! - xblsys.f: BLSYS, BLDIF (local BL equations)
//! - xsolve.f: BLSOLV (custom block solver)

use super::closure::{cf_lam, cf_turb, di_lam, dilw, hc_turb, hkin, hs_lam, hs_turb};
use nalgebra::{DMatrix, DVector};

// ============================================================================
// BL Closure Constants (from XFOIL's BLPAR.INC)
// ============================================================================

/// Shear coefficient lag constant
pub const SCCON: f64 = 5.6;

/// G-beta locus constant (G-beta relation)
pub const GACON: f64 = 6.70;

/// G-beta locus constant
pub const GBCON: f64 = 0.75;

/// Wall term constant for G-beta
pub const GCCON: f64 = 18.0;

/// Wall/wake dissipation length ratio Lo/L
pub const DLCON: f64 = 0.9;

/// Ctau weighting coefficient (derived from G-beta constants)
/// CTCON = 0.5 / (GACON² * GBCON)
pub const CTCON: f64 = 0.5 / (GACON * GACON * GBCON);

/// Skin friction factor (usually 1.0)
pub const CFFAC: f64 = 1.0;

/// Shear lag UxEQ weight
pub const DUXCON: f64 = 1.0;

/// Similarity station pressure gradient parameter (x/U dU/dx)
/// Set to 1.0 for stagnation point
pub const BULE: f64 = 1.0;

/// Initial turbulent Ctau coefficient at transition
/// CTR = CTRCON * exp(-CTRCEX/(Hk-1))
pub const CTRCON: f64 = 1.8;

/// Initial turbulent Ctau exponent at transition
pub const CTRCEX: f64 = 3.3;

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
    s1: &BLStationState,
    s2: &BLStationState,
    ampl1: f64,
    acrit: f64,
    xiforc: f64,
    params: &BLGlobalParams,
) -> TransitionResult {
    const DAEPS: f64 = 5.0e-5;
    let (x1, x2) = (s1.x, s2.x);

    // calculate average amplification rate AX over X1..X2 interval, with the current AMPL2
    let r0 = axset(
        s1.hk,
        s1.theta,
        s1.rt,
        ampl1,
        s2.hk,
        s2.theta,
        s2.rt,
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
        ut = s1.u * wf1 + s2.u * wf2;
        xt_a2 = x1 * wf1_a2 + x2 * wf2_a2;
        tt_a2 = s1.theta * wf1_a2 + s2.theta * wf2_a2;
        dt_a2 = s1.dstar * wf1_a2 + s2.dstar * wf2_a2;
        ut_a2 = s1.u * wf1_a2 + s2.u * wf2_a2;

        // temporarily set "2" variables from "T" for BLKIN (U2_UEI, U2_MS, DW2 stay station 2's)
        st = s2.clone();
        st.x = xt;
        st.theta = tt;
        st.dstar = dt;
        st.u = ut;
        st.blkin(params);

        // calculate amplification rate AX over current X1-XT interval
        r = axset(
            s1.hk,
            s1.theta,
            s1.rt,
            ampl1,
            st.hk,
            tt,
            st.rt,
            amplt,
            acrit,
            params.idampv,
        );

        // punch out early if there is no amplification here
        if r.ax <= 0.0 {
            break;
        }

        // set sensitivity of AX(A2)
        let ax_a2 = (r.ax_hk2 * st.hk_t + r.ax_t2 + r.ax_rt2 * st.rt_t) * tt_a2
            + (r.ax_hk2 * st.hk_d) * dt_a2
            + (r.ax_hk2 * st.hk_u + r.ax_rt2 * st.rt_u) * ut_a2
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
    let ut_a1 = s1.u * wf1_a1 + s2.u * wf2_a1;

    xt_x1 = x1 * wf1_x1 + x2 * wf2_x1 + xt_x1;
    let tt_x1 = s1.theta * wf1_x1 + s2.theta * wf2_x1;
    let dt_x1 = s1.dstar * wf1_x1 + s2.dstar * wf2_x1;
    let ut_x1 = s1.u * wf1_x1 + s2.u * wf2_x1;

    xt_x2 = x1 * wf1_x2 + x2 * wf2_x2 + xt_x2;
    let tt_x2 = s1.theta * wf1_x2 + s2.theta * wf2_x2;
    let dt_x2 = s1.dstar * wf1_x2 + s2.dstar * wf2_x2;
    let ut_x2 = s1.u * wf1_x2 + s2.u * wf2_x2;

    let _xt_xf = x1 * wf1_xf + x2 * wf2_xf;
    let tt_xf = s1.theta * wf1_xf + s2.theta * wf2_xf;
    let dt_xf = s1.dstar * wf1_xf + s2.dstar * wf2_xf;
    let ut_xf = s1.u * wf1_xf + s2.u * wf2_xf;

    // at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT ) from the last loop pass
    let (hkt_tt, hkt_dt, hkt_ut, hkt_ms) = (st.hk_t, st.hk_d, st.hk_u, st.hk_ms);
    let (rtt_tt, rtt_ut, rtt_ms, rtt_re) = (st.rt_t, st.rt_u, st.rt_ms, st.rt_re);
    let ax_t1 =
        r.ax_hk1 * s1.hk_t + r.ax_t1 + r.ax_rt1 * s1.rt_t + (r.ax_hk2 * hkt_tt + r.ax_t2 + r.ax_rt2 * rtt_tt) * tt_t1;
    let ax_d1 = r.ax_hk1 * s1.hk_d + (r.ax_hk2 * hkt_dt) * dt_d1;
    let ax_u1 = r.ax_hk1 * s1.hk_u + r.ax_rt1 * s1.rt_u + (r.ax_hk2 * hkt_ut + r.ax_rt2 * rtt_ut) * ut_u1;
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
    let ax_ms = r.ax_hk2 * hkt_ms + r.ax_rt2 * rtt_ms + r.ax_hk1 * s1.hk_ms + r.ax_rt1 * s1.rt_ms;
    let ax_re = r.ax_rt2 * rtt_re + r.ax_rt1 * s1.rt_re;

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

// ============================================================================
// BL Flow Type Enum
// ============================================================================

/// Type of BL flow at a station
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BLFlowType {
    /// Laminar flow
    Laminar = 1,
    /// Turbulent flow (attached)
    Turbulent = 2,
    /// Wake
    Wake = 3,
}

// ============================================================================
// Global BL Parameters
// ============================================================================

/// Global BL parameters (from XFOIL's V_VAR common block)
///
/// These parameters are constant throughout the BL calculation for a given
/// flow condition (Mach, Reynolds number).
#[derive(Debug, Clone)]
pub struct BLGlobalParams {
    /// IDAMPV: amplification model selected in SETBL from IDAMP (0 = DAMPL, 1 = DAMPL2)
    pub idampv: usize,
    /// Freestream velocity qinf
    pub qinf: f64,

    /// Karman-Tsien parameter TKBL = (1 - M²)^(-1/2) - 1 for subsonic
    pub tk: f64,
    /// d(TKBL)/d(M²)
    pub tk_ms: f64,

    /// Stagnation density ratio ρ_stag/ρ_∞
    pub rst: f64,
    /// d(RST)/d(M²)
    pub rst_ms: f64,

    /// 1 / stagnation enthalpy
    pub hstinv: f64,
    /// d(HSTINV)/d(M²)
    pub hstinv_ms: f64,

    /// Reynolds number based on freestream
    pub reybl: f64,
    /// d(REYBL)/d(M²)
    pub reybl_ms: f64,
    /// d(REYBL)/d(Re) = normalized sensitivity
    pub reybl_re: f64,

    /// Gas constants
    pub gamma: f64, // Cp/Cv (typically 1.4)
    pub gm1: f64, // gamma - 1

    /// Viscosity ratio (Hvrat in XFOIL)
    pub hvrat: f64,
}

impl BLGlobalParams {
    /// Create global BL parameters from flow conditions
    ///
    /// This is equivalent to the parameter setup in SETBL/COMSET.
    ///
    /// # Arguments
    /// * `mach` - Freestream Mach number
    /// * `reynolds` - Reynolds number based on chord
    /// * `gamma` - Ratio of specific heats (1.4 for air)
    pub fn new(mach: f64, reynolds: f64, gamma: f64) -> Self {
        let gm1 = gamma - 1.0;
        let msq = mach * mach;

        // Karman-Tsien parameter TKLAM and its M² derivative, as COMSET forms them
        let beta = (1.0 - msq).sqrt();
        let beta_msq = -0.5 / beta;
        let tk = msq / ((1.0 + beta) * (1.0 + beta));
        let tk_ms = 1.0 / ((1.0 + beta) * (1.0 + beta)) - 2.0 * tk / (1.0 + beta) * beta_msq;

        // Stagnation density ratio (isentropic)
        // RST = (1 + (γ-1)/2 M²)^(1/(γ-1))
        let tr = 1.0 + 0.5 * gm1 * msq;
        let rst = tr.powf(1.0 / gm1);
        let rst_ms = 0.5 * rst / tr;

        // 1/stagnation enthalpy
        // HSTINV = (γ-1) M² / (q²∞ (1 + (γ-1)/2 M²))
        // For normalized q∞ = 1:
        let hstinv = gm1 * msq / tr;
        let hstinv_ms = gm1 / tr - 0.5 * gm1 * hstinv / tr;

        // Reynolds number adjustment for compressibility
        // Based on edge temperature/viscosity
        let qinf = 1.0; // Normalized
                        // HVRAT (Sutherland's constant ratio) is never assigned on XFOIL's analysis path: only
                        // the plotting routines (BLPLOT/DPLOT) set it to 0.35. In a non-plotting run it keeps the
                        // static zero of uninitialised COMMON, and the viscosity law reduces to HERAT**1.5.
                        // Reproduced here (0.35 gave REYBL = 1e6 - 1 ULP on the reference case; XFOIL: 1e6).
        let hvrat = 0.0;

        let herat = 1.0 - 0.5 * qinf * qinf * hstinv;
        let herat_ms = -0.5 * qinf * qinf * hstinv_ms;

        let reybl = reynolds * (herat * herat * herat).sqrt() * (1.0 + hvrat) / (herat + hvrat);
        let reybl_re = (herat * herat * herat).sqrt() * (1.0 + hvrat) / (herat + hvrat);
        let reybl_ms = reybl * (1.5 / herat - 1.0 / (herat + hvrat)) * herat_ms;

        Self {
            idampv: 0,
            qinf,
            tk,
            tk_ms,
            rst,
            rst_ms,
            hstinv,
            hstinv_ms,
            reybl,
            reybl_ms,
            reybl_re,
            gamma,
            gm1,
            hvrat,
        }
    }

    /// Create incompressible parameters (M = 0)
    pub fn incompressible(reynolds: f64) -> Self {
        Self::new(0.0, reynolds, 1.4)
    }
}

/// BL station state variables (from XFOIL's V_VAR1/V_VAR2)
///
/// This holds all primary and derived variables at a single BL station,
/// along with their derivatives with respect to the primary variables.
///
/// The naming follows XFOIL conventions where derivatives are denoted by
/// suffixes like _u (w.r.t. U), _t (w.r.t. θ), _d (w.r.t. δ*), etc.
#[derive(Debug, Clone, Default)]
pub struct BLStationState {
    // ========================================================================
    // Primary variables (set by blprv)
    // ========================================================================
    /// Arc length position (X2 in XFOIL)
    pub x: f64,
    /// Edge velocity (compressible) (U2 in XFOIL)
    pub u: f64,
    /// Momentum thickness θ (T2 in XFOIL)
    pub theta: f64,
    /// Displacement thickness δ* without wake gap (D2 in XFOIL)
    pub dstar: f64,
    /// Shear stress coefficient Ctau (turbulent) (S2 in XFOIL)
    pub ctau: f64,
    /// Amplification factor (laminar) (AMPL2 in XFOIL)
    pub ampl: f64,
    /// Wake gap contribution to δ* (DW2 in XFOIL)
    pub dw: f64,

    // Velocity transformation (compressible to incompressible)
    /// ∂U/∂Uei (compressible w.r.t. incompressible)
    pub u_uei: f64,
    /// ∂U/∂M²
    pub u_ms: f64,

    // ========================================================================
    // Kinematic secondary variables (set by blkin)
    // ========================================================================
    /// Shape factor H = δ*/θ
    pub h: f64,
    pub h_t: f64, // ∂H/∂θ
    pub h_d: f64, // ∂H/∂δ*

    /// Edge Mach number squared M²
    pub msq: f64,
    pub msq_u: f64,  // ∂M²/∂U
    pub msq_ms: f64, // ∂M²/∂(M∞²)

    /// Density ratio ρ/ρ∞ (static to freestream)
    pub r: f64,
    pub r_u: f64,
    pub r_ms: f64,

    /// Kinematic viscosity ν ratio
    pub v: f64,
    pub v_u: f64,
    pub v_ms: f64,
    pub v_re: f64,

    /// Kinematic shape factor Hk
    pub hk: f64,
    pub hk_u: f64,  // ∂Hk/∂U
    pub hk_t: f64,  // ∂Hk/∂θ
    pub hk_d: f64,  // ∂Hk/∂δ*
    pub hk_ms: f64, // ∂Hk/∂M²

    /// Momentum Reynolds number Rθ = ρ·Ue·θ/μ
    pub rt: f64,
    pub rt_u: f64,
    pub rt_t: f64,
    pub rt_ms: f64,
    pub rt_re: f64,

    // ========================================================================
    // Turbulence-dependent secondary variables (set by blvar)
    // ========================================================================
    /// Density thickness shape factor H**
    pub hc: f64,
    pub hc_u: f64,
    pub hc_t: f64,
    pub hc_d: f64,
    pub hc_ms: f64,

    /// Energy shape factor H*
    pub hs: f64,
    pub hs_u: f64,
    pub hs_t: f64,
    pub hs_d: f64,
    pub hs_ms: f64,
    pub hs_re: f64,

    /// Normalized slip velocity Us
    pub us: f64,
    pub us_u: f64,
    pub us_t: f64,
    pub us_d: f64,
    pub us_ms: f64,
    pub us_re: f64,

    /// Equilibrium shear stress coefficient CQ (Ctau^(1/2)_eq)
    pub cq: f64,
    pub cq_u: f64,
    pub cq_t: f64,
    pub cq_d: f64,
    pub cq_ms: f64,
    pub cq_re: f64,

    /// Skin friction coefficient Cf
    pub cf: f64,
    pub cf_u: f64,
    pub cf_t: f64,
    pub cf_d: f64,
    pub cf_ms: f64,
    pub cf_re: f64,

    /// Dissipation coefficient 2*CD/H*
    pub di: f64,
    pub di_u: f64,
    pub di_t: f64,
    pub di_d: f64,
    pub di_s: f64, // ∂Di/∂Ctau
    pub di_ms: f64,
    pub di_re: f64,

    /// BL thickness δ (from Green's correlation)
    pub de: f64,
    pub de_u: f64,
    pub de_t: f64,
    pub de_d: f64,
    pub de_ms: f64,

    // ========================================================================
    // For convenience
    // ========================================================================
    /// Mass defect M = Ue·δ* (stored, not computed with derivatives here)
    pub mass: f64,
}

impl BLStationState {
    /// Set primary BL variables from input parameters (BLPRV equivalent)
    ///
    /// This converts incompressible edge velocity Uei to compressible U
    /// using the Karman-Tsien transformation.
    ///
    /// # Arguments
    /// * `xsi` - Arc length position
    /// * `ami` - Amplification factor (laminar)
    /// * `cti` - Shear stress coefficient (turbulent)
    /// * `thi` - Momentum thickness θ
    /// * `dsi` - Total displacement thickness δ* (including wake gap)
    /// * `dswaki` - Wake gap contribution
    /// * `uei` - Edge velocity (incompressible)
    /// * `params` - Global BL parameters
    pub fn blprv(
        &mut self,
        xsi: f64,
        ami: f64,
        cti: f64,
        thi: f64,
        dsi: f64,
        dswaki: f64,
        uei: f64,
        params: &BLGlobalParams,
    ) {
        self.x = xsi;
        self.ampl = ami;
        self.ctau = cti;
        self.theta = thi;
        self.dstar = dsi - dswaki; // D2 is delta* without wake gap
        self.dw = dswaki;

        // Karman-Tsien velocity transformation: Ue_compressible from Ue_incompressible
        // U2 = UEI*(1-TKBL) / (1 - TKBL*(UEI/QINFBL)^2)
        let tk = params.tk;
        let qinf = params.qinf;
        let uei_q = uei / qinf;
        let uei_q2 = uei_q * uei_q;
        let denom = 1.0 - tk * uei_q2;

        self.u = uei * (1.0 - tk) / denom;

        // Derivative: ∂U/∂Uei
        // U2_UEI = (1 + TKBL*(2*U2*UEI/QINFBL^2 - 1)) / (1 - TKBL*(UEI/QINFBL)^2)
        self.u_uei = (1.0 + tk * (2.0 * self.u * uei / (qinf * qinf) - 1.0)) / denom;

        // Derivative: ∂U/∂M² (through TK)
        // U2_MS = (U2*(UEI/QINFBL)^2 - UEI) * TKBL_MS / denom
        self.u_ms = (self.u * uei_q2 - uei) * params.tk_ms / denom;
    }

    /// Calculate turbulence-independent secondary variables (BLKIN equivalent)
    ///
    /// This computes the kinematic BL quantities that don't depend on whether
    /// the flow is laminar or turbulent: M², density ratio, shape factor H,
    /// kinematic shape factor Hk, and momentum Reynolds number Rθ.
    ///
    /// # Arguments
    /// * `params` - Global BL parameters
    pub fn blkin(&mut self, params: &BLGlobalParams) {
        let u = self.u;
        let t = self.theta;
        let d = self.dstar;
        let gm1 = params.gm1;
        let hstinv = params.hstinv;
        let hstinv_ms = params.hstinv_ms;
        let hvrat = params.hvrat;
        let rst = params.rst;
        let rst_ms = params.rst_ms;
        let reybl = params.reybl;
        let reybl_ms = params.reybl_ms;
        let reybl_re = params.reybl_re;

        // Edge Mach number squared
        // M² = U²·HSTINV / (γ-1)·(1 - 0.5·U²·HSTINV))
        let u2_hstinv = u * u * hstinv;
        let denom = gm1 * (1.0 - 0.5 * u2_hstinv);
        self.msq = u2_hstinv / denom;

        let tr = 1.0 + 0.5 * gm1 * self.msq;
        self.msq_u = 2.0 * self.msq * tr / u;
        self.msq_ms = u * u * tr / denom * hstinv_ms;

        // Edge static density (isentropic)
        // R = RST * TR^(-1/(γ-1))
        self.r = rst * tr.powf(-1.0 / gm1);
        self.r_u = -self.r / tr * 0.5 * self.msq_u;
        self.r_ms = -self.r / tr * 0.5 * self.msq_ms + rst_ms * tr.powf(-1.0 / gm1);

        // Shape factor H = δ*/θ
        self.h = d / t;
        self.h_d = 1.0 / t;
        self.h_t = -self.h / t;

        // Edge static/stagnation enthalpy ratio
        let herat = 1.0 - 0.5 * u * u * hstinv;
        let he_u = -u * hstinv;
        let he_ms = -0.5 * u * u * hstinv_ms;

        // Molecular viscosity ratio (Sutherland-type)
        // V = sqrt(HERAT^3) * (1+HVRAT)/(HERAT+HVRAT) / REYBL
        let herat_32 = (herat * herat * herat).sqrt(); // SQRT(HERAT**3): cube first, as the Fortran does
        self.v = herat_32 * (1.0 + hvrat) / (herat + hvrat) / reybl;
        let v_he = self.v * (1.5 / herat - 1.0 / (herat + hvrat));

        self.v_u = v_he * he_u;
        self.v_ms = -self.v / reybl * reybl_ms + v_he * he_ms;
        self.v_re = -self.v / reybl * reybl_re;

        // Kinematic shape factor Hk (compressibility correction)
        let (hk, hk_h, hk_msq) = hkin(self.h, self.msq);
        self.hk = hk;

        self.hk_u = hk_msq * self.msq_u;
        self.hk_t = hk_h * self.h_t;
        self.hk_d = hk_h * self.h_d;
        self.hk_ms = hk_msq * self.msq_ms;

        // Momentum thickness Reynolds number Rθ = ρ·U·θ/μ = R·U·T/V
        self.rt = self.r * u * t / self.v;
        self.rt_u = self.rt * (1.0 / u + self.r_u / self.r - self.v_u / self.v);
        self.rt_t = self.rt / t;
        self.rt_ms = self.rt * (self.r_ms / self.r - self.v_ms / self.v);
        self.rt_re = self.rt * (-self.v_re / self.v);
    }

    /// Copy all values from another station state
    ///
    /// Used for the COM1 = COM2 operation in XFOIL
    pub fn copy_from(&mut self, other: &BLStationState) {
        *self = other.clone();
    }

    /// Calculate all secondary BL variables (BLVAR equivalent)
    ///
    /// This calculates the turbulence-dependent secondary variables
    /// based on the flow type (laminar, turbulent, wake).
    ///
    /// # Arguments
    /// * `flow_type` - Type of BL flow (Laminar, Turbulent, or Wake)
    /// * `params` - Global BL parameters
    pub fn blvar(&mut self, flow_type: BLFlowType, _params: &BLGlobalParams) {
        let hk = self.hk;
        let rt = self.rt;
        let msq = self.msq;
        let h = self.h;
        let t = self.theta;
        let d = self.dstar;
        let s = self.ctau; // Ctau (or amplification for laminar)

        // XFOIL clamps HK2 itself in COMMON (BLVAR: `HK2 = MAX(HK2, ...)`) and leaves the Hk
        // derivatives as BLKIN set them. The clamped value persists: BLMID, the next station's
        // HK1 (after COM1 = COM2) and MRCHUE's HTARG all read it, so it is written back here.
        let hk = match flow_type {
            BLFlowType::Wake => hk.max(1.00005),
            _ => hk.max(1.05),
        };
        self.hk = hk;

        // ====================================================================
        // Density thickness shape parameter H** (from HCT)
        // ====================================================================
        let (hc, hc_hk, hc_msq) = hc_turb(hk, msq);

        self.hc = hc;
        self.hc_u = hc_hk * self.hk_u + hc_msq * self.msq_u;
        self.hc_t = hc_hk * self.hk_t;
        self.hc_d = hc_hk * self.hk_d;
        self.hc_ms = hc_hk * self.hk_ms + hc_msq * self.msq_ms;

        // ====================================================================
        // Energy shape factor H* (from HSL or HST)
        // ====================================================================
        let hs_result = match flow_type {
            BLFlowType::Laminar => hs_lam(hk, rt, msq),
            BLFlowType::Turbulent | BLFlowType::Wake => hs_turb(hk, rt, msq),
        };

        self.hs = hs_result.val;
        self.hs_u = hs_result.val_hk * self.hk_u + hs_result.val_rt * self.rt_u + hs_result.val_msq * self.msq_u;
        self.hs_t = hs_result.val_hk * self.hk_t + hs_result.val_rt * self.rt_t;
        self.hs_d = hs_result.val_hk * self.hk_d;
        self.hs_ms = hs_result.val_hk * self.hk_ms + hs_result.val_rt * self.rt_ms + hs_result.val_msq * self.msq_ms;
        self.hs_re = hs_result.val_rt * self.rt_re;

        // ---- normalized slip velocity  Us
        let us = 0.5 * self.hs * (1.0 - (hk - 1.0) / (GBCON * h));
        let us_hs = 0.5 * (1.0 - (hk - 1.0) / (GBCON * h));
        let us_hk = 0.5 * self.hs * (-1.0 / (GBCON * h));
        let us_h = 0.5 * self.hs * (hk - 1.0) / (GBCON * (h * h));
        self.us = us;
        self.us_u = us_hs * self.hs_u + us_hk * self.hk_u;
        self.us_t = us_hs * self.hs_t + us_hk * self.hk_t + us_h * self.h_t;
        self.us_d = us_hs * self.hs_d + us_hk * self.hk_d + us_h * self.h_d;
        self.us_ms = us_hs * self.hs_ms + us_hk * self.hk_ms;
        self.us_re = us_hs * self.hs_re;
        if flow_type != BLFlowType::Wake && self.us > 0.95 {
            self.us = 0.98;
            self.us_u = 0.0;
            self.us_t = 0.0;
            self.us_d = 0.0;
            self.us_ms = 0.0;
            self.us_re = 0.0;
        }
        if flow_type == BLFlowType::Wake && self.us > 0.99995 {
            self.us = 0.99995;
            self.us_u = 0.0;
            self.us_t = 0.0;
            self.us_d = 0.0;
            self.us_ms = 0.0;
            self.us_re = 0.0;
        }
        let us = self.us;

        // ---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
        let mut hkc = hk - 1.0;
        let mut hkc_hk = 1.0;
        let mut hkc_rt = 0.0;
        if flow_type == BLFlowType::Turbulent {
            let gcc = GCCON;
            hkc = hk - 1.0 - gcc / rt;
            hkc_hk = 1.0;
            hkc_rt = gcc / (rt * rt);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hk = 0.0;
                hkc_rt = 0.0;
            }
        }
        let hkb = hk - 1.0;
        let usb = 1.0 - us;
        let cq = (CTCON * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk))).sqrt();
        let cq_hs = CTCON * hkb * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq;
        let cq_us = CTCON * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / usb * 0.5 / cq;
        let cq_hk = CTCON * self.hs * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq
            - CTCON * self.hs * hkb * (hkc * hkc) / (usb * h * ((hk * hk) * hk)) * 2.0 * 0.5 / cq
            + CTCON * self.hs * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_hk;
        let cq_rt = CTCON * self.hs * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_rt;
        let cq_h = -(CTCON * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / h * 0.5 / cq);
        self.cq = cq;
        self.cq_u = cq_hs * self.hs_u + cq_us * self.us_u + cq_hk * self.hk_u;
        self.cq_t = cq_hs * self.hs_t + cq_us * self.us_t + cq_hk * self.hk_t;
        self.cq_d = cq_hs * self.hs_d + cq_us * self.us_d + cq_hk * self.hk_d;
        self.cq_ms = cq_hs * self.hs_ms + cq_us * self.us_ms + cq_hk * self.hk_ms;
        self.cq_re = cq_hs * self.hs_re + cq_us * self.us_re;
        self.cq_u = self.cq_u + cq_rt * self.rt_u;
        self.cq_t = self.cq_t + cq_h * self.h_t + cq_rt * self.rt_t;
        self.cq_d = self.cq_d + cq_h * self.h_d;
        self.cq_ms = self.cq_ms + cq_rt * self.rt_ms;
        self.cq_re = self.cq_re + cq_rt * self.rt_re;

        // ---- set skin friction coefficient
        let (cf, cf_hk, cf_rt, cf_m) = match flow_type {
            // wake
            BLFlowType::Wake => (0.0, 0.0, 0.0, 0.0),
            // laminar
            BLFlowType::Laminar => {
                let r = cf_lam(hk, rt, msq);
                (r.val, r.val_hk, r.val_rt, r.val_msq)
            }
            // turbulent
            BLFlowType::Turbulent => {
                let r = cf_turb(hk, rt, msq, CFFAC);
                let l = cf_lam(hk, rt, msq);
                if l.val > r.val {
                    // laminar Cf is greater than turbulent Cf -- use laminar
                    // (this will only occur for unreasonably small Rtheta)
                    (l.val, l.val_hk, l.val_rt, l.val_msq)
                } else {
                    (r.val, r.val_hk, r.val_rt, r.val_msq)
                }
            }
        };
        self.cf = cf;
        self.cf_u = cf_hk * self.hk_u + cf_rt * self.rt_u + cf_m * self.msq_u;
        self.cf_t = cf_hk * self.hk_t + cf_rt * self.rt_t;
        self.cf_d = cf_hk * self.hk_d;
        self.cf_ms = cf_hk * self.hk_ms + cf_rt * self.rt_ms + cf_m * self.msq_ms;
        self.cf_re = cf_rt * self.rt_re;

        // ---- dissipation function    2 CD / H*
        match flow_type {
            BLFlowType::Laminar => {
                // laminar
                let r = di_lam(hk, rt);
                self.di = r.val;
                self.di_u = r.val_hk * self.hk_u + r.val_rt * self.rt_u;
                self.di_t = r.val_hk * self.hk_t + r.val_rt * self.rt_t;
                self.di_d = r.val_hk * self.hk_d;
                self.di_s = 0.0;
                self.di_ms = r.val_hk * self.hk_ms + r.val_rt * self.rt_ms;
                self.di_re = r.val_rt * self.rt_re;
            }
            BLFlowType::Turbulent => {
                // turbulent wall contribution
                let c = cf_turb(hk, rt, msq, CFFAC);
                let cf2t = c.val;
                let cf2t_u = c.val_hk * self.hk_u + c.val_rt * self.rt_u + c.val_msq * self.msq_u;
                let cf2t_t = c.val_hk * self.hk_t + c.val_rt * self.rt_t;
                let cf2t_d = c.val_hk * self.hk_d;
                let cf2t_ms = c.val_hk * self.hk_ms + c.val_rt * self.rt_ms + c.val_msq * self.msq_ms;
                let cf2t_re = c.val_rt * self.rt_re;
                let mut di = (0.5 * cf2t * us) * 2.0 / self.hs;
                let di_hs = -((0.5 * cf2t * us) * 2.0 / (self.hs * self.hs));
                let di_us = (0.5 * cf2t) * 2.0 / self.hs;
                let di_cf2t = (0.5 * us) * 2.0 / self.hs;
                let mut di_s = 0.0;
                let mut di_u = di_hs * self.hs_u + di_us * self.us_u + di_cf2t * cf2t_u;
                let mut di_t = di_hs * self.hs_t + di_us * self.us_t + di_cf2t * cf2t_t;
                let mut di_d = di_hs * self.hs_d + di_us * self.us_d + di_cf2t * cf2t_d;
                let mut di_ms = di_hs * self.hs_ms + di_us * self.us_ms + di_cf2t * cf2t_ms;
                let mut di_re = di_hs * self.hs_re + di_us * self.us_re + di_cf2t * cf2t_re;

                // set minimum Hk for wake layer to still exist
                let grt = rt.ln();
                let hmin = 1.0 + 2.1 / grt;
                let hm_rt = -(2.1 / (grt * grt)) / rt;

                // set factor DFAC for correcting wall dissipation for very low Hk
                let fl = (hk - 1.0) / (hmin - 1.0);
                let fl_hk = 1.0 / (hmin - 1.0);
                let fl_rt = (-fl / (hmin - 1.0)) * hm_rt;
                let tfl = fl.tanh();
                let dfac = 0.5 + 0.5 * tfl;
                let df_fl = 0.5 * (1.0 - tfl * tfl);
                let df_hk = df_fl * fl_hk;
                let df_rt = df_fl * fl_rt;

                di_s *= dfac;
                di_u = di_u * dfac + di * (df_hk * self.hk_u + df_rt * self.rt_u);
                di_t = di_t * dfac + di * (df_hk * self.hk_t + df_rt * self.rt_t);
                di_d = di_d * dfac + di * (df_hk * self.hk_d);
                di_ms = di_ms * dfac + di * (df_hk * self.hk_ms + df_rt * self.rt_ms);
                di_re = di_re * dfac + di * (df_rt * self.rt_re);
                di *= dfac;

                self.di = di;
                self.di_s = di_s;
                self.di_u = di_u;
                self.di_t = di_t;
                self.di_d = di_d;
                self.di_ms = di_ms;
                self.di_re = di_re;
            }
            BLFlowType::Wake => {
                // zero wall contribution for wake
                self.di = 0.0;
                self.di_s = 0.0;
                self.di_u = 0.0;
                self.di_t = 0.0;
                self.di_d = 0.0;
                self.di_ms = 0.0;
                self.di_re = 0.0;
            }
        }

        // ---- Add on turbulent outer layer contribution
        if flow_type != BLFlowType::Laminar {
            let dd = (s * s) * (0.995 - us) * 2.0 / self.hs;
            let dd_hs = -((s * s) * (0.995 - us) * 2.0 / (self.hs * self.hs));
            let dd_us = -((s * s) * 2.0 / self.hs);
            let dd_s = s * 2.0 * (0.995 - us) * 2.0 / self.hs;
            self.di = self.di + dd;
            self.di_s = dd_s;
            self.di_u = self.di_u + dd_hs * self.hs_u + dd_us * self.us_u;
            self.di_t = self.di_t + dd_hs * self.hs_t + dd_us * self.us_t;
            self.di_d = self.di_d + dd_hs * self.hs_d + dd_us * self.us_d;
            self.di_ms = self.di_ms + dd_hs * self.hs_ms + dd_us * self.us_ms;
            self.di_re = self.di_re + dd_hs * self.hs_re + dd_us * self.us_re;

            // add laminar stress contribution to outer layer CD
            let dd = 0.15 * ((0.995 - us) * (0.995 - us)) / rt * 2.0 / self.hs;
            let dd_us = -0.15 * (0.995 - us) * 2.0 / rt * 2.0 / self.hs;
            let dd_hs = -dd / self.hs;
            let dd_rt = -dd / rt;
            self.di = self.di + dd;
            self.di_u = self.di_u + dd_hs * self.hs_u + dd_us * self.us_u + dd_rt * self.rt_u;
            self.di_t = self.di_t + dd_hs * self.hs_t + dd_us * self.us_t + dd_rt * self.rt_t;
            self.di_d = self.di_d + dd_hs * self.hs_d + dd_us * self.us_d;
            self.di_ms = self.di_ms + dd_hs * self.hs_ms + dd_us * self.us_ms + dd_rt * self.rt_ms;
            self.di_re = self.di_re + dd_hs * self.hs_re + dd_us * self.us_re + dd_rt * self.rt_re;
        }

        if flow_type == BLFlowType::Turbulent {
            let l = di_lam(hk, rt);
            if l.val > self.di {
                // laminar CD is greater than turbulent CD -- use laminar
                // (this will only occur for unreasonably small Rtheta)
                self.di = l.val;
                self.di_s = 0.0;
                self.di_u = l.val_hk * self.hk_u + l.val_rt * self.rt_u;
                self.di_t = l.val_hk * self.hk_t + l.val_rt * self.rt_t;
                self.di_d = l.val_hk * self.hk_d;
                self.di_ms = l.val_hk * self.hk_ms + l.val_rt * self.rt_ms;
                self.di_re = l.val_rt * self.rt_re;
            }
        }

        if flow_type == BLFlowType::Wake {
            // laminar wake CD
            let l = dilw(hk, rt);
            if l.val > self.di {
                // laminar wake CD is greater than turbulent CD -- use laminar
                self.di = l.val;
                self.di_s = 0.0;
                self.di_u = l.val_hk * self.hk_u + l.val_rt * self.rt_u;
                self.di_t = l.val_hk * self.hk_t + l.val_rt * self.rt_t;
                self.di_d = l.val_hk * self.hk_d;
                self.di_ms = l.val_hk * self.hk_ms + l.val_rt * self.rt_ms;
                self.di_re = l.val_rt * self.rt_re;
            }
        }

        if flow_type == BLFlowType::Wake {
            // double dissipation for the wake (two wake halves)
            self.di *= 2.0;
            self.di_s *= 2.0;
            self.di_u *= 2.0;
            self.di_t *= 2.0;
            self.di_d *= 2.0;
            self.di_ms *= 2.0;
            self.di_re *= 2.0;
        }

        // ====================================================================
        // BL thickness Delta (Green's correlation)
        // DE = (3.15 + 1.72/(HK-1)) * T + D
        // ====================================================================
        let de = (3.15 + 1.72 / (hk - 1.0)) * t + d;
        let de_hk = -1.72 / ((hk - 1.0) * (hk - 1.0)) * t;

        self.de = de;
        self.de_u = de_hk * self.hk_u;
        self.de_t = de_hk * self.hk_t + 3.15 + 1.72 / (hk - 1.0);
        self.de_d = de_hk * self.hk_d + 1.0;
        self.de_ms = de_hk * self.hk_ms;

        // Clamp DE to reasonable values
        let hdmax = 12.0;
        if self.de > hdmax * t {
            self.de = hdmax * t;
            self.de_u = 0.0;
            self.de_t = hdmax;
            self.de_d = 0.0;
            self.de_ms = 0.0;
        }
    }
}

// ============================================================================
// Midpoint Skin Friction (BLMID)
// ============================================================================

/// Midpoint skin friction result (from XFOIL's BLMID)
///
/// This holds the skin friction coefficient at the midpoint between two
/// stations, along with its derivatives.
#[derive(Debug, Clone, Default)]
pub struct MidpointCf {
    /// Midpoint skin friction coefficient
    pub cfm: f64,

    /// Derivatives w.r.t. station 1
    pub cfm_u1: f64,
    pub cfm_t1: f64,
    pub cfm_d1: f64,

    /// Derivatives w.r.t. station 2
    pub cfm_u2: f64,
    pub cfm_t2: f64,
    pub cfm_d2: f64,

    /// Derivatives w.r.t. global parameters
    pub cfm_ms: f64,
    pub cfm_re: f64,
}

impl MidpointCf {
    /// Calculate midpoint skin friction (BLMID equivalent)
    ///
    /// Calculates the skin friction coefficient at the midpoint between
    /// two BL stations. For turbulent flow, uses the maximum of turbulent
    /// and laminar Cf.
    ///
    /// # Arguments
    /// * `s1` - Station 1 state
    /// * `s2` - Station 2 state
    /// * `flow_type` - Type of BL flow
    /// * `is_similarity` - True if this is a similarity station (copy s2→s1)
    pub fn compute(s1: &BLStationState, s2: &BLStationState, flow_type: BLFlowType, is_similarity: bool) -> Self {
        let mut result = Self::default();

        // For similarity station, station 1 equals station 2
        let (hk1, rt1, m1) = if is_similarity {
            (s2.hk, s2.rt, s2.msq)
        } else {
            (s1.hk, s1.rt, s1.msq)
        };

        let (hk1_u1, hk1_t1, hk1_d1, hk1_ms) = if is_similarity {
            (s2.hk_u, s2.hk_t, s2.hk_d, s2.hk_ms)
        } else {
            (s1.hk_u, s1.hk_t, s1.hk_d, s1.hk_ms)
        };

        let (rt1_u1, rt1_t1, rt1_ms, rt1_re) = if is_similarity {
            (s2.rt_u, s2.rt_t, s2.rt_ms, s2.rt_re)
        } else {
            (s1.rt_u, s1.rt_t, s1.rt_ms, s1.rt_re)
        };

        let (m1_u1, m1_ms) = if is_similarity {
            (s2.msq_u, s2.msq_ms)
        } else {
            (s1.msq_u, s1.msq_ms)
        };

        // Midpoint averages
        let hka = 0.5 * (hk1 + s2.hk);
        let rta = 0.5 * (rt1 + s2.rt);
        let ma = 0.5 * (m1 + s2.msq);

        // Midpoint skin friction coefficient
        let (cfm, cfm_hka, cfm_rta, cfm_ma) = match flow_type {
            BLFlowType::Wake => {
                // Zero skin friction in wake
                (0.0, 0.0, 0.0, 0.0)
            }
            BLFlowType::Laminar => {
                // Laminar Cf
                let cf_result = cf_lam(hka, rta, ma);
                (cf_result.val, cf_result.val_hk, cf_result.val_rt, 0.0)
            }
            BLFlowType::Turbulent => {
                // Turbulent Cf
                let cf_turb_result = cf_turb(hka, rta, ma, CFFAC);
                // Check if laminar is higher
                let cf_lam_result = cf_lam(hka, rta, ma);

                if cf_lam_result.val > cf_turb_result.val {
                    (cf_lam_result.val, cf_lam_result.val_hk, cf_lam_result.val_rt, 0.0)
                } else {
                    (
                        cf_turb_result.val,
                        cf_turb_result.val_hk,
                        cf_turb_result.val_rt,
                        cf_turb_result.val_msq,
                    )
                }
            }
        };

        result.cfm = cfm;

        // Derivatives w.r.t. station 1 variables (factor of 0.5 from averaging)
        result.cfm_u1 = 0.5 * (cfm_hka * hk1_u1 + cfm_ma * m1_u1 + cfm_rta * rt1_u1);
        result.cfm_t1 = 0.5 * (cfm_hka * hk1_t1 + cfm_rta * rt1_t1);
        result.cfm_d1 = 0.5 * cfm_hka * hk1_d1;

        // Derivatives w.r.t. station 2 variables
        result.cfm_u2 = 0.5 * (cfm_hka * s2.hk_u + cfm_ma * s2.msq_u + cfm_rta * s2.rt_u);
        result.cfm_t2 = 0.5 * (cfm_hka * s2.hk_t + cfm_rta * s2.rt_t);
        result.cfm_d2 = 0.5 * cfm_hka * s2.hk_d;

        // Derivatives w.r.t. global parameters
        result.cfm_ms = 0.5
            * (cfm_hka * hk1_ms
                + cfm_ma * m1_ms
                + cfm_rta * rt1_ms
                + cfm_hka * s2.hk_ms
                + cfm_ma * s2.msq_ms
                + cfm_rta * s2.rt_ms);
        result.cfm_re = 0.5 * (cfm_rta * rt1_re + cfm_rta * s2.rt_re);

        result
    }
}

// ============================================================================
// Local BL Equation Coefficients
// ============================================================================

/// Upwinding parameter calculation
///
/// Returns (upw, upw_u1, upw_t1, upw_d1, upw_u2, upw_t2, upw_d2, upw_ms)
fn compute_upwinding(s1: &BLStationState, s2: &BLStationState, is_wake: bool) -> UpwindParams {
    let hk1 = s1.hk;
    let hk2 = s2.hk;

    // Upwinding constant (less in wake)
    let hupwt = 1.0;
    let hdcon = if is_wake {
        hupwt / (hk2 * hk2)
    } else {
        5.0 * hupwt / (hk2 * hk2)
    };
    let hd_hk1 = 0.0;
    let hd_hk2 = -hdcon * 2.0 / hk2;

    // Local upwinding based on log(Hk-1) change
    let arg = ((hk2 - 1.0) / (hk1 - 1.0)).abs();
    let hl = arg.ln();
    let hl_hk1 = -1.0 / (hk1 - 1.0);
    let hl_hk2 = 1.0 / (hk2 - 1.0);

    // Upwinding parameter: 0.5 = trapezoidal, 1.0 = backward Euler
    let hlsq = (hl * hl).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;
    let upw_hl = ehh * hl * hdcon;
    let upw_hd = 0.5 * ehh * hlsq;

    let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    UpwindParams {
        upw,
        upw_u1: upw_hk1 * s1.hk_u,
        upw_t1: upw_hk1 * s1.hk_t,
        upw_d1: upw_hk1 * s1.hk_d,
        upw_u2: upw_hk2 * s2.hk_u,
        upw_t2: upw_hk2 * s2.hk_t,
        upw_d2: upw_hk2 * s2.hk_d,
        upw_ms: upw_hk1 * s1.hk_ms + upw_hk2 * s2.hk_ms,
    }
}

/// Upwinding parameters
#[derive(Debug, Clone, Default)]
struct UpwindParams {
    upw: f64,
    upw_u1: f64,
    upw_t1: f64,
    upw_d1: f64,
    upw_u2: f64,
    upw_t2: f64,
    upw_d2: f64,
    upw_ms: f64,
}

/// Local BL equation coefficients (from XFOIL's V_SYS)
///
/// These are the Jacobian entries for a single station pair (1→2).
#[derive(Debug, Clone, Default)]
pub struct BLLocalSystem {
    /// Jacobian w.r.t. previous station: VS1(4,5)
    /// Rows: 4 equations (momentum, shape, lag, auxiliary)
    /// Cols: 5 variables (Ctau, Theta, Dstar, Ue, X)
    pub vs1: [[f64; 5]; 4],

    /// Jacobian w.r.t. current station: VS2(4,5)
    pub vs2: [[f64; 5]; 4],

    /// Residual vector: VSREZ(4)
    pub vsrez: [f64; 4],

    /// Sensitivity to Reynolds number: VSR(4)
    pub vsr: [f64; 4],

    /// Sensitivity to Mach squared: VSM(4)
    pub vsm: [f64; 4],

    /// Sensitivity to arc length: VSX(4)
    pub vsx: [f64; 4],
}

impl BLLocalSystem {
    /// Set up the Newton system for a BL interval (BLDIF equivalent)
    ///
    /// This sets up the Jacobian and residual for the BL equations between
    /// two stations. The equations are:
    /// - Row 1: Amplification (laminar) or Shear lag (turbulent/wake)
    /// - Row 2: Momentum integral
    /// - Row 3: Shape parameter (energy)
    ///
    /// # Arguments
    /// * `s1` - Station 1 state (upstream)
    /// * `s2` - Station 2 state (downstream)
    /// * `cfm` - Midpoint skin friction
    /// * `flow_type` - Type of BL flow
    /// * `is_similarity` - True if station 2 is a similarity station (LE)
    pub fn bldif(
        &mut self,
        s1: &BLStationState,
        s2: &BLStationState,
        cfm: &MidpointCf,
        flow_type: BLFlowType,
        is_similarity: bool,
        acrit: f64,
        idampv: usize,
    ) {
        // Initialize to zero
        for k in 0..4 {
            self.vsrez[k] = 0.0;
            self.vsm[k] = 0.0;
            self.vsr[k] = 0.0;
            self.vsx[k] = 0.0;
            for l in 0..5 {
                self.vs1[k][l] = 0.0;
                self.vs2[k][l] = 0.0;
            }
        }

        // Logarithmic differences
        let (xlog, ulog, tlog, hlog, ddlog) = if is_similarity {
            // Similarity station: prescribed differences
            (1.0, BULE, 0.5 * (1.0 - BULE), 0.0, 0.0)
        } else {
            // Normal station: compute from values
            let xlog = (s2.x / s1.x).ln();
            let ulog = (s2.u / s1.u).ln();
            let tlog = (s2.theta / s1.theta).ln();
            let hlog = (s2.hs / s1.hs).ln();
            (xlog, ulog, tlog, hlog, 1.0)
        };

        // Compute upwinding parameters
        let is_wake = flow_type == BLFlowType::Wake;
        let upw = compute_upwinding(s1, s2, is_wake);

        // Equation 1: Amplification (laminar) or Shear lag (turbulent/wake)
        match flow_type {
            BLFlowType::Laminar if is_similarity => {
                // LE point: set zero amplification factor (XFOIL: VS2(1,1) = 1.0)
                // This ensures the pivot in BLSOLV is well-conditioned
                self.vs2[0][0] = 1.0;
                self.vsrez[0] = -s2.ampl;
            }
            BLFlowType::Laminar => {
                // laminar part --> set amplification equation (BLDIF ITYP=1), verbatim:
                // set average amplification AX over interval X1..X2
                let r = axset(
                    s1.hk, s1.theta, s1.rt, s1.ampl, s2.hk, s2.theta, s2.rt, s2.ampl, acrit, idampv,
                );
                let ax = r.ax;
                let rezc = s2.ampl - s1.ampl - ax * (s2.x - s1.x);
                let z_ax = -(s2.x - s1.x);

                self.vs1[0][0] = z_ax * r.ax_a1 - 1.0;
                self.vs1[0][1] = z_ax * (r.ax_hk1 * s1.hk_t + r.ax_t1 + r.ax_rt1 * s1.rt_t);
                self.vs1[0][2] = z_ax * (r.ax_hk1 * s1.hk_d);
                self.vs1[0][3] = z_ax * (r.ax_hk1 * s1.hk_u + r.ax_rt1 * s1.rt_u);
                self.vs1[0][4] = ax;
                self.vs2[0][0] = z_ax * r.ax_a2 + 1.0;
                self.vs2[0][1] = z_ax * (r.ax_hk2 * s2.hk_t + r.ax_t2 + r.ax_rt2 * s2.rt_t);
                self.vs2[0][2] = z_ax * (r.ax_hk2 * s2.hk_d);
                self.vs2[0][3] = z_ax * (r.ax_hk2 * s2.hk_u + r.ax_rt2 * s2.rt_u);
                self.vs2[0][4] = -ax;
                self.vsm[0] =
                    z_ax * (r.ax_hk1 * s1.hk_ms + r.ax_rt1 * s1.rt_ms + r.ax_hk2 * s2.hk_ms + r.ax_rt2 * s2.rt_ms);
                self.vsr[0] = z_ax * (r.ax_rt1 * s1.rt_re + r.ax_rt2 * s2.rt_re);
                self.vsx[0] = 0.0;
                self.vsrez[0] = -rezc;
            }
            BLFlowType::Turbulent | BLFlowType::Wake => {
                // Shear lag equation
                self.setup_shear_lag_equation(s1, s2, &upw, flow_type);
            }
        }

        // Equation 2: Momentum integral equation
        self.setup_momentum_equation(s1, s2, cfm, xlog, ulog, tlog, ddlog);

        // Equation 3: Shape parameter equation
        self.setup_shape_equation(s1, s2, &upw, xlog, ulog, hlog, ddlog);
    }

    /// BLDIF (xblsys.f) "turbulent part --> set shear lag equation" (row 1), line for line.
    fn setup_shear_lag_equation(
        &mut self,
        s1: &BLStationState,
        s2: &BLStationState,
        upw: &UpwindParams,
        flow_type: BLFlowType,
    ) {
        let u = upw.upw;
        let sa = (1.0 - u) * s1.ctau + u * s2.ctau;
        let cqa = (1.0 - u) * s1.cq + u * s2.cq;
        let cfa = (1.0 - u) * s1.cf + u * s2.cf;
        let hka = (1.0 - u) * s1.hk + u * s2.hk;
        let usa = 0.5 * (s1.us + s2.us);
        let rta = 0.5 * (s1.rt + s2.rt);
        let dea = 0.5 * (s1.de + s2.de);
        let da = 0.5 * (s1.dstar + s2.dstar);
        // increased dissipation length in wake (decrease its reciprocal)
        let ald = if flow_type == BLFlowType::Wake { DLCON } else { 1.0 };

        // set and linearize  equilibrium 1/Ue dUe/dx   ...  NEW  12 Oct 94
        let (hkc, hkc_hka, hkc_rta) = if flow_type == BLFlowType::Turbulent {
            let gcc = GCCON;
            let mut hkc = hka - 1.0 - gcc / rta;
            let mut hkc_hka = 1.0;
            let mut hkc_rta = gcc / (rta * rta);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hka = 0.0;
                hkc_rta = 0.0;
            }
            (hkc, hkc_hka, hkc_rta)
        } else {
            (hka - 1.0, 1.0, 0.0)
        };
        let hr = hkc / (GACON * ald * hka);
        let hr_hka = hkc_hka / (GACON * ald * hka) - hr / hka;
        let _hr_rta = hkc_rta / (GACON * ald * hka);
        let uq = (0.5 * cfa - hr * hr) / (GBCON * da);
        let uq_hka = -2.0 * hr * hr_hka / (GBCON * da);
        let uq_cfa = 0.5 / (GBCON * da);
        let uq_da = -uq / da;
        // (XFOIL also forms UQ_RTA and UQ_T1..UQ_RE here; none of them enter the Jacobian below)

        let scc = SCCON * 1.333 / (1.0 + usa);
        let scc_usa = -scc / (1.0 + usa);
        let scc_us1 = scc_usa * 0.5;
        let scc_us2 = scc_usa * 0.5;
        let _ = (scc_us1, scc_us2);

        let slog = (s2.ctau / s1.ctau).ln();
        let dxi = s2.x - s1.x;
        let ulog = (s2.u / s1.u).ln();

        let rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog + dea * 2.0 * (uq * dxi - ulog) * DUXCON;

        let z_cfa = dea * 2.0 * uq_cfa * dxi * DUXCON;
        let z_hka = dea * 2.0 * uq_hka * dxi * DUXCON;
        let z_da = dea * 2.0 * uq_da * dxi * DUXCON;
        let z_sl = -dea * 2.0;
        let z_ul = -dea * 2.0 * DUXCON;
        let z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * DUXCON;
        let z_usa = scc_usa * (cqa - sa * ald) * dxi;
        let z_cqa = scc * dxi;
        let z_sa = -scc * dxi * ald;
        let z_dea = 2.0 * ((uq * dxi - ulog) * DUXCON - slog);
        let z_upw =
            z_cqa * (s2.cq - s1.cq) + z_sa * (s2.ctau - s1.ctau) + z_cfa * (s2.cf - s1.cf) + z_hka * (s2.hk - s1.hk);

        let z_de1 = 0.5 * z_dea;
        let z_de2 = 0.5 * z_dea;
        let z_us1 = 0.5 * z_usa;
        let z_us2 = 0.5 * z_usa;
        let z_d1 = 0.5 * z_da;
        let z_d2 = 0.5 * z_da;
        let z_u1 = -z_ul / s1.u;
        let z_u2 = z_ul / s2.u;
        let z_x1 = -z_dxi;
        let z_x2 = z_dxi;
        let z_s1 = (1.0 - u) * z_sa - z_sl / s1.ctau;
        let z_s2 = u * z_sa + z_sl / s2.ctau;
        let z_cq1 = (1.0 - u) * z_cqa;
        let z_cq2 = u * z_cqa;
        let z_cf1 = (1.0 - u) * z_cfa;
        let z_cf2 = u * z_cfa;
        let z_hk1 = (1.0 - u) * z_hka;
        let z_hk2 = u * z_hka;

        self.vs1[0][0] = z_s1;
        self.vs1[0][1] = z_upw * upw.upw_t1 + z_de1 * s1.de_t + z_us1 * s1.us_t;
        self.vs1[0][2] = z_d1 + z_upw * upw.upw_d1 + z_de1 * s1.de_d + z_us1 * s1.us_d;
        self.vs1[0][3] = z_u1 + z_upw * upw.upw_u1 + z_de1 * s1.de_u + z_us1 * s1.us_u;
        self.vs1[0][4] = z_x1;
        self.vs2[0][0] = z_s2;
        self.vs2[0][1] = z_upw * upw.upw_t2 + z_de2 * s2.de_t + z_us2 * s2.us_t;
        self.vs2[0][2] = z_d2 + z_upw * upw.upw_d2 + z_de2 * s2.de_d + z_us2 * s2.us_d;
        self.vs2[0][3] = z_u2 + z_upw * upw.upw_u2 + z_de2 * s2.de_u + z_us2 * s2.us_u;
        self.vs2[0][4] = z_x2;
        self.vsm[0] = z_upw * upw.upw_ms + z_de1 * s1.de_ms + z_us1 * s1.us_ms + z_de2 * s2.de_ms + z_us2 * s2.us_ms;

        self.vs1[0][1] = self.vs1[0][1] + z_cq1 * s1.cq_t + z_cf1 * s1.cf_t + z_hk1 * s1.hk_t;
        self.vs1[0][2] = self.vs1[0][2] + z_cq1 * s1.cq_d + z_cf1 * s1.cf_d + z_hk1 * s1.hk_d;
        self.vs1[0][3] = self.vs1[0][3] + z_cq1 * s1.cq_u + z_cf1 * s1.cf_u + z_hk1 * s1.hk_u;
        self.vs2[0][1] = self.vs2[0][1] + z_cq2 * s2.cq_t + z_cf2 * s2.cf_t + z_hk2 * s2.hk_t;
        self.vs2[0][2] = self.vs2[0][2] + z_cq2 * s2.cq_d + z_cf2 * s2.cf_d + z_hk2 * s2.hk_d;
        self.vs2[0][3] = self.vs2[0][3] + z_cq2 * s2.cq_u + z_cf2 * s2.cf_u + z_hk2 * s2.hk_u;
        self.vsm[0] = self.vsm[0]
            + z_cq1 * s1.cq_ms
            + z_cf1 * s1.cf_ms
            + z_hk1 * s1.hk_ms
            + z_cq2 * s2.cq_ms
            + z_cf2 * s2.cf_ms
            + z_hk2 * s2.hk_ms;
        self.vsr[0] = z_cq1 * s1.cq_re + z_cf1 * s1.cf_re + z_cq2 * s2.cq_re + z_cf2 * s2.cf_re;
        self.vsx[0] = 0.0;
        self.vsrez[0] = -rezc;
    }

    /// Set up the momentum integral equation (row 2)
    fn setup_momentum_equation(
        &mut self,
        s1: &BLStationState,
        s2: &BLStationState,
        cfm: &MidpointCf,
        xlog: f64,
        ulog: f64,
        tlog: f64,
        ddlog: f64,
    ) {
        // Averaged values
        let ha = 0.5 * (s1.h + s2.h);
        let ma = 0.5 * (s1.msq + s2.msq);
        let xa = 0.5 * (s1.x + s2.x);
        let ta = 0.5 * (s1.theta + s2.theta);
        let hwa = 0.5 * (s1.dw / s1.theta + s2.dw / s2.theta);

        // Cf term using central CFM for accuracy
        let cfx = 0.5 * cfm.cfm * xa / ta + 0.25 * (s1.cf * s1.x / s1.theta + s2.cf * s2.x / s2.theta);
        let cfx_xa = 0.5 * cfm.cfm / ta;
        let cfx_ta = -0.5 * cfm.cfm * xa / (ta * ta);
        let cfx_x1 = 0.25 * s1.cf / s1.theta + cfx_xa * 0.5;
        let cfx_x2 = 0.25 * s2.cf / s2.theta + cfx_xa * 0.5;
        let cfx_t1 = -0.25 * s1.cf * s1.x / (s1.theta * s1.theta) + cfx_ta * 0.5;
        let cfx_t2 = -0.25 * s2.cf * s2.x / (s2.theta * s2.theta) + cfx_ta * 0.5;
        let cfx_cf1 = 0.25 * s1.x / s1.theta;
        let cfx_cf2 = 0.25 * s2.x / s2.theta;
        let cfx_cfm = 0.5 * xa / ta;

        let btmp = ha + 2.0 - ma + hwa;

        // Momentum equation residual
        let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;

        // Z coefficients
        let z_cfx = -xlog * 0.5;
        let z_ha = ulog;
        let z_hwa = ulog;
        let z_ma = -ulog;
        let z_xl = -ddlog * 0.5 * cfx;
        let z_ul = ddlog * btmp;
        let z_tl = ddlog;

        let z_cfm = z_cfx * cfx_cfm;
        let z_cf1 = z_cfx * cfx_cf1;
        let z_cf2 = z_cfx * cfx_cf2;

        let z_t1 = -z_tl / s1.theta + z_cfx * cfx_t1 + z_hwa * 0.5 * (-s1.dw / (s1.theta * s1.theta));
        let z_t2 = z_tl / s2.theta + z_cfx * cfx_t2 + z_hwa * 0.5 * (-s2.dw / (s2.theta * s2.theta));
        let z_x1 = -z_xl / s1.x + z_cfx * cfx_x1;
        let z_x2 = z_xl / s2.x + z_cfx * cfx_x2;
        let z_u1 = -z_ul / s1.u;
        let z_u2 = z_ul / s2.u;

        // Jacobian entries for row 2
        self.vs1[1][1] = 0.5 * z_ha * s1.h_t + z_cfm * cfm.cfm_t1 + z_cf1 * s1.cf_t + z_t1;
        self.vs1[1][2] = 0.5 * z_ha * s1.h_d + z_cfm * cfm.cfm_d1 + z_cf1 * s1.cf_d;
        self.vs1[1][3] = 0.5 * z_ma * s1.msq_u + z_cfm * cfm.cfm_u1 + z_cf1 * s1.cf_u + z_u1;
        self.vs1[1][4] = z_x1;

        self.vs2[1][1] = 0.5 * z_ha * s2.h_t + z_cfm * cfm.cfm_t2 + z_cf2 * s2.cf_t + z_t2;
        self.vs2[1][2] = 0.5 * z_ha * s2.h_d + z_cfm * cfm.cfm_d2 + z_cf2 * s2.cf_d;
        self.vs2[1][3] = 0.5 * z_ma * s2.msq_u + z_cfm * cfm.cfm_u2 + z_cf2 * s2.cf_u + z_u2;
        self.vs2[1][4] = z_x2;

        self.vsm[1] =
            0.5 * z_ma * s1.msq_ms + z_cfm * cfm.cfm_ms + z_cf1 * s1.cf_ms + 0.5 * z_ma * s2.msq_ms + z_cf2 * s2.cf_ms;
        self.vsr[1] = z_cfm * cfm.cfm_re + z_cf1 * s1.cf_re + z_cf2 * s2.cf_re;

        self.vsrez[1] = -rezt;
    }

    /// Set up the shape parameter equation (row 3)
    fn setup_shape_equation(
        &mut self,
        s1: &BLStationState,
        s2: &BLStationState,
        upw: &UpwindParams,
        xlog: f64,
        ulog: f64,
        hlog: f64,
        ddlog: f64,
    ) {
        let u = upw.upw;
        let xot1 = s1.x / s1.theta;
        let xot2 = s2.x / s2.theta;

        // Averaged values
        let ha = 0.5 * (s1.h + s2.h);
        let hsa = 0.5 * (s1.hs + s2.hs);
        let hca = 0.5 * (s1.hc + s2.hc);
        let hwa = 0.5 * (s1.dw / s1.theta + s2.dw / s2.theta);

        // Upwind-weighted DI and CF
        let dix = (1.0 - u) * s1.di * xot1 + u * s2.di * xot2;
        let cfx = (1.0 - u) * s1.cf * xot1 + u * s2.cf * xot2;
        let dix_upw = s2.di * xot2 - s1.di * xot1;
        let cfx_upw = s2.cf * xot2 - s1.cf * xot1;

        let btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

        // Shape equation residual
        let rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);

        // Z coefficients
        let z_cfx = xlog * 0.5;
        let z_dix = -xlog;
        let z_hca = 2.0 * ulog / hsa;
        let z_ha = -ulog;
        let z_hwa = -ulog;
        let z_xl = ddlog * (0.5 * cfx - dix);
        let z_ul = ddlog * btmp;
        let z_hl = ddlog;

        let z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

        let z_hs1 = -hca * ulog / (hsa * hsa) - z_hl / s1.hs;
        let z_hs2 = -hca * ulog / (hsa * hsa) + z_hl / s2.hs;

        let z_cf1 = (1.0 - u) * z_cfx * xot1;
        let z_cf2 = u * z_cfx * xot2;
        let z_di1 = (1.0 - u) * z_dix * xot1;
        let z_di2 = u * z_dix * xot2;

        let z_t1 = (1.0 - u) * (z_cfx * s1.cf + z_dix * s1.di) * (-xot1 / s1.theta)
            + z_hwa * 0.5 * (-s1.dw / (s1.theta * s1.theta));
        let z_t2 =
            u * (z_cfx * s2.cf + z_dix * s2.di) * (-xot2 / s2.theta) + z_hwa * 0.5 * (-s2.dw / (s2.theta * s2.theta));
        let z_x1 = (1.0 - u) * (z_cfx * s1.cf + z_dix * s1.di) / s1.theta - z_xl / s1.x;
        let z_x2 = u * (z_cfx * s2.cf + z_dix * s2.di) / s2.theta + z_xl / s2.x;
        let z_u1 = -z_ul / s1.u;
        let z_u2 = z_ul / s2.u;

        // Jacobian entries for row 3
        self.vs1[2][0] = z_di1 * s1.di_s;
        self.vs1[2][1] = z_hs1 * s1.hs_t
            + z_cf1 * s1.cf_t
            + z_di1 * s1.di_t
            + z_t1
            + 0.5 * (z_hca * s1.hc_t + z_ha * s1.h_t)
            + z_upw * upw.upw_t1;
        self.vs1[2][2] = z_hs1 * s1.hs_d
            + z_cf1 * s1.cf_d
            + z_di1 * s1.di_d
            + 0.5 * (z_hca * s1.hc_d + z_ha * s1.h_d)
            + z_upw * upw.upw_d1;
        self.vs1[2][3] =
            z_hs1 * s1.hs_u + z_cf1 * s1.cf_u + z_di1 * s1.di_u + z_u1 + 0.5 * z_hca * s1.hc_u + z_upw * upw.upw_u1;
        self.vs1[2][4] = z_x1;

        self.vs2[2][0] = z_di2 * s2.di_s;
        self.vs2[2][1] = z_hs2 * s2.hs_t
            + z_cf2 * s2.cf_t
            + z_di2 * s2.di_t
            + z_t2
            + 0.5 * (z_hca * s2.hc_t + z_ha * s2.h_t)
            + z_upw * upw.upw_t2;
        self.vs2[2][2] = z_hs2 * s2.hs_d
            + z_cf2 * s2.cf_d
            + z_di2 * s2.di_d
            + 0.5 * (z_hca * s2.hc_d + z_ha * s2.h_d)
            + z_upw * upw.upw_d2;
        self.vs2[2][3] =
            z_hs2 * s2.hs_u + z_cf2 * s2.cf_u + z_di2 * s2.di_u + z_u2 + 0.5 * z_hca * s2.hc_u + z_upw * upw.upw_u2;
        self.vs2[2][4] = z_x2;

        self.vsm[2] = z_hs1 * s1.hs_ms
            + z_cf1 * s1.cf_ms
            + z_di1 * s1.di_ms
            + z_hs2 * s2.hs_ms
            + z_cf2 * s2.cf_ms
            + z_di2 * s2.di_ms
            + 0.5 * (z_hca * s1.hc_ms + z_hca * s2.hc_ms)
            + z_upw * upw.upw_ms;
        self.vsr[2] = z_hs1 * s1.hs_re
            + z_cf1 * s1.cf_re
            + z_di1 * s1.di_re
            + z_hs2 * s2.hs_re
            + z_cf2 * s2.cf_re
            + z_di2 * s2.di_re;

        self.vsrez[2] = -rezh;
    }

    /// Set up Newton system for transition interval (TRDIF equivalent)
    ///
    /// Handles intervals that span laminar-turbulent transition by:
    /// 1. Setting up laminar equations from X1 to XT (transition)
    /// 2. Setting up turbulent equations from XT to X2
    /// 3. Summing the contributions
    ///
    /// # Arguments
    /// * `s1` - Station 1 state (upstream, laminar)
    /// * `s2` - Station 2 state (downstream, turbulent)
    /// * `trans` - Transition location and derivatives
    /// * `acrit` - Critical amplification factor
    /// * `params` - Global BL parameters
    #[allow(clippy::too_many_lines)]
    pub fn trdif(
        &mut self,
        s1: &BLStationState,
        s2: &BLStationState,
        trans: &TransitionLocation,
        acrit: f64,
        params: &BLGlobalParams,
    ) {
        // Weighting factors for linear interpolation to transition point
        let wf2 = (trans.xt - s1.x) / (s2.x - s1.x);
        let wf2_xt = 1.0 / (s2.x - s1.x);
        let wf1 = 1.0 - wf2;

        // Derivatives of weighting factors w.r.t. station variables
        let wf2_a1 = wf2_xt * trans.xt_a1;
        let wf2_x1 = wf2_xt * trans.xt_x1 + (wf2 - 1.0) / (s2.x - s1.x);
        let wf2_x2 = wf2_xt * trans.xt_x2 - wf2 / (s2.x - s1.x);
        let wf2_t1 = wf2_xt * trans.xt_t1;
        let wf2_t2 = wf2_xt * trans.xt_t2;
        let wf2_d1 = wf2_xt * trans.xt_d1;
        let wf2_d2 = wf2_xt * trans.xt_d2;
        let wf2_u1 = wf2_xt * trans.xt_u1;
        let wf2_u2 = wf2_xt * trans.xt_u2;
        let wf2_ms = wf2_xt * trans.xt_ms;
        let wf2_re = wf2_xt * trans.xt_re;
        let wf2_xf = wf2_xt * trans.xt_xf;

        let wf1_a1 = -wf2_a1;
        let wf1_x1 = -wf2_x1;
        let wf1_x2 = -wf2_x2;
        let wf1_t1 = -wf2_t1;
        let wf1_t2 = -wf2_t2;
        let wf1_d1 = -wf2_d1;
        let wf1_d2 = -wf2_d2;
        let wf1_u1 = -wf2_u1;
        let wf1_u2 = -wf2_u2;
        let wf1_ms = -wf2_ms;
        let wf1_re = -wf2_re;
        let wf1_xf = -wf2_xf;

        // *** PART 1: Laminar from X1 to XT ***

        // Interpolate primary variables to transition point
        let tt = s1.theta * wf1 + s2.theta * wf2;
        let tt_a1 = s1.theta * wf1_a1 + s2.theta * wf2_a1;
        let tt_x1 = s1.theta * wf1_x1 + s2.theta * wf2_x1;
        let tt_x2 = s1.theta * wf1_x2 + s2.theta * wf2_x2;
        let tt_t1 = s1.theta * wf1_t1 + s2.theta * wf2_t1 + wf1;
        let tt_t2 = s1.theta * wf1_t2 + s2.theta * wf2_t2 + wf2;
        let tt_d1 = s1.theta * wf1_d1 + s2.theta * wf2_d1;
        let tt_d2 = s1.theta * wf1_d2 + s2.theta * wf2_d2;
        let tt_u1 = s1.theta * wf1_u1 + s2.theta * wf2_u1;
        let tt_u2 = s1.theta * wf1_u2 + s2.theta * wf2_u2;
        let tt_ms = s1.theta * wf1_ms + s2.theta * wf2_ms;
        let tt_re = s1.theta * wf1_re + s2.theta * wf2_re;
        let tt_xf = s1.theta * wf1_xf + s2.theta * wf2_xf;

        let dt = s1.dstar * wf1 + s2.dstar * wf2;
        let dt_a1 = s1.dstar * wf1_a1 + s2.dstar * wf2_a1;
        let dt_x1 = s1.dstar * wf1_x1 + s2.dstar * wf2_x1;
        let dt_x2 = s1.dstar * wf1_x2 + s2.dstar * wf2_x2;
        let dt_t1 = s1.dstar * wf1_t1 + s2.dstar * wf2_t1;
        let dt_t2 = s1.dstar * wf1_t2 + s2.dstar * wf2_t2;
        let dt_d1 = s1.dstar * wf1_d1 + s2.dstar * wf2_d1 + wf1;
        let dt_d2 = s1.dstar * wf1_d2 + s2.dstar * wf2_d2 + wf2;
        let dt_u1 = s1.dstar * wf1_u1 + s2.dstar * wf2_u1;
        let dt_u2 = s1.dstar * wf1_u2 + s2.dstar * wf2_u2;
        let dt_ms = s1.dstar * wf1_ms + s2.dstar * wf2_ms;
        let dt_re = s1.dstar * wf1_re + s2.dstar * wf2_re;
        let dt_xf = s1.dstar * wf1_xf + s2.dstar * wf2_xf;

        let ut = s1.u * wf1 + s2.u * wf2;
        let ut_a1 = s1.u * wf1_a1 + s2.u * wf2_a1;
        let ut_x1 = s1.u * wf1_x1 + s2.u * wf2_x1;
        let ut_x2 = s1.u * wf1_x2 + s2.u * wf2_x2;
        let ut_t1 = s1.u * wf1_t1 + s2.u * wf2_t1;
        let ut_t2 = s1.u * wf1_t2 + s2.u * wf2_t2;
        let ut_d1 = s1.u * wf1_d1 + s2.u * wf2_d1;
        let ut_d2 = s1.u * wf1_d2 + s2.u * wf2_d2;
        let ut_u1 = s1.u * wf1_u1 + s2.u * wf2_u1 + wf1;
        let ut_u2 = s1.u * wf1_u2 + s2.u * wf2_u2 + wf2;
        let ut_ms = s1.u * wf1_ms + s2.u * wf2_ms;
        let ut_re = s1.u * wf1_re + s2.u * wf2_re;
        let ut_xf = s1.u * wf1_xf + s2.u * wf2_xf;

        // Create transition-point state for laminar part
        // set primary "T" variables at XT (really placed into "2" variables): XFOIL overwrites
        // X2/T2/D2/U2/AMPL2/S2 on the saved station-2 COMMON, so U2_UEI, U2_MS and DW2 are
        // those of station 2 — no BLPRV here.
        let mut st = s2.clone();
        st.x = trans.xt;
        st.theta = tt;
        st.dstar = dt;
        st.u = ut;
        st.ampl = acrit;
        st.ctau = 0.0;
        st.blkin(params);
        st.blvar(BLFlowType::Laminar, params);

        // Calculate midpoint Cf for X1-XT
        let cfm_lam = MidpointCf::compute(s1, &st, BLFlowType::Laminar, false);

        // Call BLDIF for laminar part (X1 to XT)
        let mut lam_sys = BLLocalSystem::default();
        lam_sys.bldif(s1, &st, &cfm_lam, BLFlowType::Laminar, false, acrit, params.idampv);

        // Convert laminar system sensitivities from "T" variables to "1" and "2" variables
        // Using chain rule for derivatives
        let mut bl1: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut bl2: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut blrez: [f64; 4] = [0.0; 4];
        let mut blm: [f64; 4] = [0.0; 4];
        let mut blr: [f64; 4] = [0.0; 4];
        let mut blx: [f64; 4] = [0.0; 4];

        for k in 1..3 {
            // Row 2 and 3 (momentum and shape)
            blrez[k] = lam_sys.vsrez[k];
            blm[k] = lam_sys.vsm[k]
                + lam_sys.vs2[k][1] * tt_ms
                + lam_sys.vs2[k][2] * dt_ms
                + lam_sys.vs2[k][3] * ut_ms
                + lam_sys.vs2[k][4] * trans.xt_ms;
            blr[k] = lam_sys.vsr[k]
                + lam_sys.vs2[k][1] * tt_re
                + lam_sys.vs2[k][2] * dt_re
                + lam_sys.vs2[k][3] * ut_re
                + lam_sys.vs2[k][4] * trans.xt_re;
            blx[k] = lam_sys.vsx[k]
                + lam_sys.vs2[k][1] * tt_xf
                + lam_sys.vs2[k][2] * dt_xf
                + lam_sys.vs2[k][3] * ut_xf
                + lam_sys.vs2[k][4] * trans.xt_xf;

            bl1[k][0] = lam_sys.vs1[k][0]
                + lam_sys.vs2[k][1] * tt_a1
                + lam_sys.vs2[k][2] * dt_a1
                + lam_sys.vs2[k][3] * ut_a1
                + lam_sys.vs2[k][4] * trans.xt_a1;
            bl1[k][1] = lam_sys.vs1[k][1]
                + lam_sys.vs2[k][1] * tt_t1
                + lam_sys.vs2[k][2] * dt_t1
                + lam_sys.vs2[k][3] * ut_t1
                + lam_sys.vs2[k][4] * trans.xt_t1;
            bl1[k][2] = lam_sys.vs1[k][2]
                + lam_sys.vs2[k][1] * tt_d1
                + lam_sys.vs2[k][2] * dt_d1
                + lam_sys.vs2[k][3] * ut_d1
                + lam_sys.vs2[k][4] * trans.xt_d1;
            bl1[k][3] = lam_sys.vs1[k][3]
                + lam_sys.vs2[k][1] * tt_u1
                + lam_sys.vs2[k][2] * dt_u1
                + lam_sys.vs2[k][3] * ut_u1
                + lam_sys.vs2[k][4] * trans.xt_u1;
            bl1[k][4] = lam_sys.vs1[k][4]
                + lam_sys.vs2[k][1] * tt_x1
                + lam_sys.vs2[k][2] * dt_x1
                + lam_sys.vs2[k][3] * ut_x1
                + lam_sys.vs2[k][4] * trans.xt_x1;

            bl2[k][0] = 0.0; // No dA2 dependence (A2 is turbulent Ctau)
            bl2[k][1] = lam_sys.vs2[k][1] * tt_t2
                + lam_sys.vs2[k][2] * dt_t2
                + lam_sys.vs2[k][3] * ut_t2
                + lam_sys.vs2[k][4] * trans.xt_t2;
            bl2[k][2] = lam_sys.vs2[k][1] * tt_d2
                + lam_sys.vs2[k][2] * dt_d2
                + lam_sys.vs2[k][3] * ut_d2
                + lam_sys.vs2[k][4] * trans.xt_d2;
            bl2[k][3] = lam_sys.vs2[k][1] * tt_u2
                + lam_sys.vs2[k][2] * dt_u2
                + lam_sys.vs2[k][3] * ut_u2
                + lam_sys.vs2[k][4] * trans.xt_u2;
            bl2[k][4] = lam_sys.vs2[k][1] * tt_x2
                + lam_sys.vs2[k][2] * dt_x2
                + lam_sys.vs2[k][3] * ut_x2
                + lam_sys.vs2[k][4] * trans.xt_x2;
        }

        // *** PART 2: Turbulent from XT to X2 ***

        // Calculate equilibrium shear coefficient CQT at transition
        st.blvar(BLFlowType::Turbulent, params);

        // Set initial shear stress: ST = CTR * CQ
        // where CTR = CTRCON * exp(-CTRCEX/(HK-1))
        let hk_minus_one = st.hk - 1.0;
        let ctr = CTRCON * (-CTRCEX / hk_minus_one).exp();
        let ctr_hk = ctr * CTRCEX / (hk_minus_one * hk_minus_one);

        let s_t = ctr * st.cq;
        let st_tt = ctr * st.cq_t + st.cq * ctr_hk * st.hk_t;
        let st_dt = ctr * st.cq_d + st.cq * ctr_hk * st.hk_d;
        let st_ut = ctr * st.cq_u + st.cq * ctr_hk * st.hk_u;
        let st_ms = ctr * st.cq_ms + st.cq * ctr_hk * st.hk_ms;
        let st_re = ctr * st.cq_re;

        // ST sensitivities w.r.t. actual "1" and "2" variables
        let st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
        let st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
        let st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
        let st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
        let st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
        let st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
        let st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
        let st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
        let st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
        let st_ms_total = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
        let st_re_total = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
        let st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

        // Update transition station with turbulent initial condition
        st.ctau = s_t;

        // Recalculate turbulent secondary variables with proper CTI
        st.blvar(BLFlowType::Turbulent, params);

        // Calculate midpoint Cf for XT-X2
        let cfm_turb = MidpointCf::compute(&st, s2, BLFlowType::Turbulent, false);

        // Call BLDIF for turbulent part (XT to X2)
        let mut turb_sys = BLLocalSystem::default();
        turb_sys.bldif(&st, s2, &cfm_turb, BLFlowType::Turbulent, false, acrit, params.idampv);

        // Convert turbulent system sensitivities from "T" variables to "1" and "2" variables
        let mut bt1: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut bt2: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut btrez: [f64; 4] = [0.0; 4];
        let mut btm: [f64; 4] = [0.0; 4];
        let mut btr: [f64; 4] = [0.0; 4];
        let mut btx: [f64; 4] = [0.0; 4];

        for k in 0..3 {
            btrez[k] = turb_sys.vsrez[k];
            btm[k] = turb_sys.vsm[k]
                + turb_sys.vs1[k][0] * st_ms_total
                + turb_sys.vs1[k][1] * tt_ms
                + turb_sys.vs1[k][2] * dt_ms
                + turb_sys.vs1[k][3] * ut_ms
                + turb_sys.vs1[k][4] * trans.xt_ms;
            btr[k] = turb_sys.vsr[k]
                + turb_sys.vs1[k][0] * st_re_total
                + turb_sys.vs1[k][1] * tt_re
                + turb_sys.vs1[k][2] * dt_re
                + turb_sys.vs1[k][3] * ut_re
                + turb_sys.vs1[k][4] * trans.xt_re;
            btx[k] = turb_sys.vsx[k]
                + turb_sys.vs1[k][0] * st_xf
                + turb_sys.vs1[k][1] * tt_xf
                + turb_sys.vs1[k][2] * dt_xf
                + turb_sys.vs1[k][3] * ut_xf
                + turb_sys.vs1[k][4] * trans.xt_xf;

            bt1[k][0] = turb_sys.vs1[k][0] * st_a1
                + turb_sys.vs1[k][1] * tt_a1
                + turb_sys.vs1[k][2] * dt_a1
                + turb_sys.vs1[k][3] * ut_a1
                + turb_sys.vs1[k][4] * trans.xt_a1;
            bt1[k][1] = turb_sys.vs1[k][0] * st_t1
                + turb_sys.vs1[k][1] * tt_t1
                + turb_sys.vs1[k][2] * dt_t1
                + turb_sys.vs1[k][3] * ut_t1
                + turb_sys.vs1[k][4] * trans.xt_t1;
            bt1[k][2] = turb_sys.vs1[k][0] * st_d1
                + turb_sys.vs1[k][1] * tt_d1
                + turb_sys.vs1[k][2] * dt_d1
                + turb_sys.vs1[k][3] * ut_d1
                + turb_sys.vs1[k][4] * trans.xt_d1;
            bt1[k][3] = turb_sys.vs1[k][0] * st_u1
                + turb_sys.vs1[k][1] * tt_u1
                + turb_sys.vs1[k][2] * dt_u1
                + turb_sys.vs1[k][3] * ut_u1
                + turb_sys.vs1[k][4] * trans.xt_u1;
            bt1[k][4] = turb_sys.vs1[k][0] * st_x1
                + turb_sys.vs1[k][1] * tt_x1
                + turb_sys.vs1[k][2] * dt_x1
                + turb_sys.vs1[k][3] * ut_x1
                + turb_sys.vs1[k][4] * trans.xt_x1;

            bt2[k][0] = turb_sys.vs2[k][0];
            bt2[k][1] = turb_sys.vs2[k][1]
                + turb_sys.vs1[k][0] * st_t2
                + turb_sys.vs1[k][1] * tt_t2
                + turb_sys.vs1[k][2] * dt_t2
                + turb_sys.vs1[k][3] * ut_t2
                + turb_sys.vs1[k][4] * trans.xt_t2;
            bt2[k][2] = turb_sys.vs2[k][2]
                + turb_sys.vs1[k][0] * st_d2
                + turb_sys.vs1[k][1] * tt_d2
                + turb_sys.vs1[k][2] * dt_d2
                + turb_sys.vs1[k][3] * ut_d2
                + turb_sys.vs1[k][4] * trans.xt_d2;
            bt2[k][3] = turb_sys.vs2[k][3]
                + turb_sys.vs1[k][0] * st_u2
                + turb_sys.vs1[k][1] * tt_u2
                + turb_sys.vs1[k][2] * dt_u2
                + turb_sys.vs1[k][3] * ut_u2
                + turb_sys.vs1[k][4] * trans.xt_u2;
            bt2[k][4] = turb_sys.vs2[k][4]
                + turb_sys.vs1[k][0] * st_x2
                + turb_sys.vs1[k][1] * tt_x2
                + turb_sys.vs1[k][2] * dt_x2
                + turb_sys.vs1[k][3] * ut_x2
                + turb_sys.vs1[k][4] * trans.xt_x2;
        }

        // *** COMBINE: Add laminar and turbulent parts ***

        // Row 1: Shear stress (from turbulent part only - laminar row is amplification)
        self.vsrez[0] = btrez[0];
        self.vsm[0] = btm[0];
        self.vsr[0] = btr[0];
        self.vsx[0] = btx[0];
        for l in 0..5 {
            self.vs1[0][l] = bt1[0][l];
            self.vs2[0][l] = bt2[0][l];
        }

        // Rows 2 and 3: Sum laminar and turbulent contributions
        for k in 1..3 {
            self.vsrez[k] = blrez[k] + btrez[k];
            self.vsm[k] = blm[k] + btm[k];
            self.vsr[k] = blr[k] + btr[k];
            self.vsx[k] = blx[k] + btx[k];
            for l in 0..5 {
                self.vs1[k][l] = bl1[k][l] + bt1[k][l];
                self.vs2[k][l] = bl2[k][l] + bt2[k][l];
            }
        }
    }
}

/// Global Newton system for coupled BL equations
///
/// This is the full system that BLSOLV solves.
pub struct BLNewtonSystem {
    /// Number of stations (NSYS in XFOIL)
    pub n_sys: usize,

    /// Diagonal block: VA(3,2,n_sys) - 3 rows, 2 cols per station
    /// In XFOIL: VA(3,2,IZX)
    pub va: Vec<[[f64; 2]; 3]>,

    /// Sub-diagonal block: VB(3,2,n_sys)
    pub vb: Vec<[[f64; 2]; 3]>,

    /// Right-hand side / residual: VDEL(3,2,n_sys)
    /// Column 1: residual, Column 2: Reynolds sensitivity
    pub vdel: Vec<[[f64; 2]; 3]>,

    /// Dense mass coupling: VM(3,n_sys,n_sys)
    /// This contains the DIJ influence terms
    pub vm: DMatrix<f64>,

    /// TE coupling block: VZ(3,2)
    pub vz: [[f64; 2]; 3],

    /// RMS residual
    pub rms_bl: f64,

    /// Max residual
    pub rmx_bl: f64,

    /// Relaxation factor
    pub rlx: f64,
}

impl BLNewtonSystem {
    /// Create new Newton system for given number of stations
    pub fn new(n_sys: usize) -> Self {
        Self {
            n_sys,
            va: vec![[[0.0; 2]; 3]; n_sys],
            vb: vec![[[0.0; 2]; 3]; n_sys],
            vdel: vec![[[0.0; 2]; 3]; n_sys],
            vm: DMatrix::zeros(3 * n_sys, n_sys),
            vz: [[0.0; 2]; 3],
            rms_bl: 0.0,
            rmx_bl: 0.0,
            rlx: 1.0,
        }
    }

    /// Solve the Newton system (BLSOLV equivalent)
    ///
    /// This implements custom block elimination for the coupled BL system:
    /// - Block-tridiagonal from BL equations (VA, VB)
    /// - Dense coupling from mass defect influence (VM)
    ///
    /// For now, this is a simplified solver that ignores the dense coupling
    /// and solves just the block-tridiagonal part. Full XFOIL-style solver
    /// will be implemented later.
    ///
    /// Returns the solution vector [d1, d2, d3, ...] where each di = [dCtau, dTheta, dMass]
    pub fn solve(&mut self) -> DVector<f64> {
        if self.n_sys == 0 {
            return DVector::zeros(0);
        }

        // Build full 3n x 3n matrix and solve directly
        // This is less efficient than XFOIL's custom solver but correct
        let n = self.n_sys;
        let size = 3 * n;
        let mut matrix = DMatrix::zeros(size, size);
        let mut rhs = DVector::zeros(size);

        // Fill matrix from block structure
        for iv in 0..n {
            let row_base = 3 * iv;

            // Diagonal block VA (3x2) and VM diagonal (row 3)
            for k in 0..2 {
                matrix[(row_base + k, row_base)] = self.va[iv][k][0];
                matrix[(row_base + k, row_base + 1)] = self.va[iv][k][1];
            }
            // Row 3 (mass equation) - use VM for coupling
            for l in 0..n {
                matrix[(row_base + 2, 3 * l + 2)] = self.vm[(2, l)];
            }
            // Rows 1,2 also have VM coupling
            for l in 0..n {
                matrix[(row_base, 3 * l + 2)] = self.vm[(0, l)];
                matrix[(row_base + 1, 3 * l + 2)] = self.vm[(1, l)];
            }

            // Sub-diagonal block VB
            if iv > 0 {
                let col_base = 3 * (iv - 1);
                for k in 0..3 {
                    matrix[(row_base + k, col_base)] = self.vb[iv][k][0];
                    matrix[(row_base + k, col_base + 1)] = self.vb[iv][k][1];
                }
            }

            // RHS
            for k in 0..3 {
                rhs[row_base + k] = self.vdel[iv][k][0];
            }
        }

        // Solve using LU decomposition
        match matrix.clone().lu().solve(&rhs) {
            Some(solution) => solution,
            None => {
                // Fallback: return zeros if matrix is singular
                DVector::zeros(size)
            }
        }
    }

    /// Assemble the BL Newton system for one surface (SETBL equivalent)
    ///
    /// This marches through all stations on the surface, setting up the local
    /// BL equation coefficients and assembling them into the global system.
    ///
    /// # Arguments
    /// * `surface` - BL surface state containing all station data
    /// * `params` - Global BL parameters
    /// * `acrit` - Critical amplification factor for transition
    /// * `xiforc` - Forced transition location (set > max_x to disable)
    ///
    /// # Returns
    /// Transition location (station index)
    pub fn assemble_surface(
        &mut self,
        surface: &mut BLSurfaceState,
        params: &BLGlobalParams,
        acrit: f64,
        xiforc: f64,
    ) -> usize {
        let mut itran = 0; // Transition station index
        let mut turb = false; // Currently in turbulent region
        let mut ampl1 = 0.0; // Amplification at previous station

        // Initialize first station residuals
        if surface.nbl < 2 {
            return 0;
        }

        // March through all stations
        for ibl in 1..surface.nbl {
            let iv = surface.isys[ibl]; // Global system index
            let is_simi = ibl == 1; // First station is similarity
            let is_wake = ibl > surface.iblte;
            let is_tran = itran == 0 && !turb && ibl > 1;

            // Get station states
            let s1 = &surface.stations[ibl - 1];
            let s2 = &surface.stations[ibl];

            // Determine flow type
            let flow_type = if is_wake {
                BLFlowType::Wake
            } else if turb {
                BLFlowType::Turbulent
            } else {
                BLFlowType::Laminar
            };

            // Check for transition if laminar
            if is_tran {
                let tr_result = trchek(s1, s2, ampl1, acrit, xiforc, params);
                match tr_result {
                    TransitionResult::FreeTransition { location, ampl2 } => {
                        itran = ibl;
                        turb = true;
                        ampl1 = ampl2;

                        // Handle transition interval with TRDIF
                        let mut local_sys = BLLocalSystem::default();
                        local_sys.trdif(s1, s2, &location, acrit, params);

                        // Copy to global system
                        self.copy_local_to_global(&local_sys, iv);
                        continue;
                    }
                    TransitionResult::ForcedTransition { location } => {
                        itran = ibl;
                        turb = true;

                        // Handle forced transition
                        let mut local_sys = BLLocalSystem::default();
                        local_sys.trdif(s1, s2, &location, acrit, params);

                        // Copy to global system
                        self.copy_local_to_global(&local_sys, iv);
                        continue;
                    }
                    TransitionResult::NoTransition { ampl2 } => {
                        ampl1 = ampl2;
                    }
                }
            }

            // Calculate midpoint skin friction
            let cfm = MidpointCf::compute(s1, s2, flow_type, is_wake);

            // Set up local BL system
            let mut local_sys = BLLocalSystem::default();
            local_sys.bldif(s1, s2, &cfm, flow_type, is_simi, 9.0, params.idampv);

            // Handle similarity station: "1" and "2" variables are the same
            if is_simi {
                for k in 0..4 {
                    for l in 0..5 {
                        local_sys.vs2[k][l] += local_sys.vs1[k][l];
                        local_sys.vs1[k][l] = 0.0;
                    }
                }
            }

            // Copy to global system
            self.copy_local_to_global(&local_sys, iv);
        }

        // Store transition index
        surface.itran = itran;
        itran
    }

    /// Copy local system coefficients to global system
    fn copy_local_to_global(&mut self, local: &BLLocalSystem, iv: usize) {
        // Copy diagonal block (vs2[0:3, 0:2])
        for k in 0..3 {
            self.va[iv][k][0] = local.vs2[k][0]; // dCtau or dAmpl
            self.va[iv][k][1] = local.vs2[k][1]; // dTheta
        }

        // Copy sub-diagonal block (vs1[0:3, 0:2])
        for k in 0..3 {
            self.vb[iv][k][0] = local.vs1[k][0];
            self.vb[iv][k][1] = local.vs1[k][1];
        }

        // Copy residual
        for k in 0..3 {
            self.vdel[iv][k][0] = local.vsrez[k];
            self.vdel[iv][k][1] = 0.0; // Will be filled with Re/Mach sensitivity
        }

        // Track RMS residual
        for k in 0..3 {
            let res = local.vsrez[k].abs();
            self.rms_bl += res * res;
            if res > self.rmx_bl {
                self.rmx_bl = res;
            }
        }
    }

    /// Finalize RMS residual after assembling both surfaces
    pub fn finalize_residual(&mut self) {
        if self.n_sys > 0 {
            self.rms_bl = (self.rms_bl / (3.0 * self.n_sys as f64)).sqrt();
        }
    }
}

/// BL solution state for one surface (upper or lower)
pub struct BLSurfaceState {
    /// Station states
    pub stations: Vec<BLStationState>,

    /// Inviscid edge velocity at each station
    pub uinv: Vec<f64>,

    /// Current edge velocity at each station
    pub uedg: Vec<f64>,

    /// Mass defect at each station (Ue * δ*)
    pub mass: Vec<f64>,

    /// Momentum thickness at each station
    pub thet: Vec<f64>,

    /// Displacement thickness at each station
    pub dstr: Vec<f64>,

    /// Shear coefficient at each station
    pub ctau: Vec<f64>,

    /// Arc length at each station
    pub xssi: Vec<f64>,

    /// Panel index for each BL station
    pub ipan: Vec<usize>,

    /// System index for each BL station
    pub isys: Vec<usize>,

    /// Number of BL stations
    pub nbl: usize,

    /// Transition station index
    pub itran: usize,

    /// TE station index
    pub iblte: usize,
}

impl BLSurfaceState {
    /// Create new surface state with given capacity
    pub fn new(capacity: usize) -> Self {
        Self {
            stations: Vec::with_capacity(capacity),
            uinv: Vec::with_capacity(capacity),
            uedg: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            thet: Vec::with_capacity(capacity),
            dstr: Vec::with_capacity(capacity),
            ctau: Vec::with_capacity(capacity),
            xssi: Vec::with_capacity(capacity),
            ipan: Vec::with_capacity(capacity),
            isys: Vec::with_capacity(capacity),
            nbl: 0,
            itran: 0,
            iblte: 0,
        }
    }
}

// ============================================================================
// UPDATE - Apply Newton Step with Underrelaxation (XFOIL UPDATE equivalent)
// ============================================================================

/// Result from UPDATE
#[derive(Debug, Clone, Default)]
pub struct UpdateResult {
    /// RMS of normalized changes
    pub rms_bl: f64,
    /// Maximum normalized change
    pub rmx_bl: f64,
    /// Variable with maximum change ('C', 'T', 'D', 'U', or 'n')
    pub vmx_bl: char,
    /// Station with maximum change
    pub imx_bl: usize,
    /// Surface with maximum change (0 = upper, 1 = lower)
    pub ismx_bl: usize,
    /// Underrelaxation factor used
    pub rlx: f64,
}

/// Limit displacement thickness to keep Hk above minimum (DSLIM equivalent)
///
/// Adjusts δ* to prevent kinematic shape factor from dropping below HKLIM.
///
/// # Arguments
/// * `dstr` - Displacement thickness (modified in place)
/// * `thet` - Momentum thickness
/// * `uedg` - Edge velocity
/// * `msq` - Edge Mach number squared
/// * `hklim` - Minimum kinematic shape factor
pub fn dslim(dstr: &mut f64, thet: f64, _uedg: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet;
    let (hk, hk_h, _hk_m) = hkin(h, msq);

    // If Hk is below limit, adjust δ*
    let dh = (hklim - hk).max(0.0) / hk_h;
    *dstr += dh * thet;
}

/// Apply Newton deltas to BL variables with underrelaxation (UPDATE equivalent)
///
/// This adds Newton deltas to boundary layer variables, checking for excessive
/// changes and underrelaxing if necessary. Also calculates max and RMS changes.
///
/// # Arguments
/// * `upper` - Upper surface BL state
/// * `lower` - Lower surface BL state
/// * `solution` - Solution vector from BLSOLV [dCtau1, dTheta1, dMass1, ...]
/// * `dij` - DIJ influence matrix (n_panels x n_panels)
/// * `vti` - Velocity sign at each BL station (±1)
/// * `params` - Global BL parameters
///
/// # Returns
/// UpdateResult with convergence information
pub fn update(
    upper: &mut BLSurfaceState,
    lower: &mut BLSurfaceState,
    solution: &DVector<f64>,
    dij: &DMatrix<f64>,
    vti_upper: &[f64],
    vti_lower: &[f64],
    params: &BLGlobalParams,
) -> UpdateResult {
    let mut result = UpdateResult::default();

    // Max allowable changes per iteration (from XFOIL)
    let dhi = 1.5; // Max increase factor
    let dlo = -0.5; // Max decrease factor

    // Hstinv for Mach calculation
    let hstinv = params.hstinv;
    let gm1 = params.gm1;

    // Initialize relaxation factor
    let mut rlx = 1.0;
    let mut rmsbl = 0.0;
    let mut rmxbl = 0.0_f64;
    let mut vmxbl = ' ';
    let mut imxbl = 0;
    let mut ismxbl = 0;

    // =========================================================================
    // Step 1: Calculate new Ue distribution from mass defect changes (DIJ)
    // =========================================================================
    // UNEW(IBL,IS) = UINV(IBL,IS) + sum_j(UE_M * (MASS + dMASS))
    // where UE_M = -VTI(IBL)*VTI(JBL)*DIJ(I,J)

    let mut unew_upper = vec![0.0; upper.nbl];
    let mut unew_lower = vec![0.0; lower.nbl];

    // Calculate new Ue for upper surface
    for ibl in 1..upper.nbl {
        let i = upper.ipan[ibl];
        let mut dui = 0.0;

        // Contribution from upper surface
        for jbl in 1..upper.nbl {
            let j = upper.ipan[jbl];
            let jv = upper.isys[jbl];
            if jv < solution.len() / 3 {
                let ue_m = -vti_upper[ibl] * vti_upper[jbl] * dij[(i, j)];
                let dmass = solution[3 * jv + 2]; // dMass is 3rd variable
                dui += ue_m * (upper.mass[jbl] + dmass);
            }
        }

        // Contribution from lower surface
        for jbl in 1..lower.nbl {
            let j = lower.ipan[jbl];
            let jv = lower.isys[jbl];
            if jv < solution.len() / 3 {
                let ue_m = -vti_upper[ibl] * vti_lower[jbl] * dij[(i, j)];
                let dmass = solution[3 * jv + 2];
                dui += ue_m * (lower.mass[jbl] + dmass);
            }
        }

        unew_upper[ibl] = upper.uinv[ibl] + dui;
    }

    // Calculate new Ue for lower surface
    for ibl in 1..lower.nbl {
        let i = lower.ipan[ibl];
        let mut dui = 0.0;

        // Contribution from upper surface
        for jbl in 1..upper.nbl {
            let j = upper.ipan[jbl];
            let jv = upper.isys[jbl];
            if jv < solution.len() / 3 {
                let ue_m = -vti_lower[ibl] * vti_upper[jbl] * dij[(i, j)];
                let dmass = solution[3 * jv + 2];
                dui += ue_m * (upper.mass[jbl] + dmass);
            }
        }

        // Contribution from lower surface
        for jbl in 1..lower.nbl {
            let j = lower.ipan[jbl];
            let jv = lower.isys[jbl];
            if jv < solution.len() / 3 {
                let ue_m = -vti_lower[ibl] * vti_lower[jbl] * dij[(i, j)];
                let dmass = solution[3 * jv + 2];
                dui += ue_m * (lower.mass[jbl] + dmass);
            }
        }

        unew_lower[ibl] = lower.uinv[ibl] + dui;
    }

    // =========================================================================
    // Step 2: Calculate changes and check for underrelaxation
    // =========================================================================

    let surfaces: [(&mut BLSurfaceState, &[f64]); 2] = [(upper, &unew_upper[..]), (lower, &unew_lower[..])];

    // First pass: determine relaxation factor
    for (is, (surface, unew)) in surfaces.iter().enumerate() {
        for ibl in 1..surface.nbl {
            let iv = surface.isys[ibl];
            if iv * 3 + 2 >= solution.len() {
                continue;
            }

            // Get Newton deltas
            let dctau = solution[3 * iv];
            let dthet = solution[3 * iv + 1];
            let dmass = solution[3 * iv + 2];
            let duedg = unew[ibl] - surface.uedg[ibl];
            let ddstr = (dmass - surface.dstr[ibl] * duedg) / surface.uedg[ibl];

            // Normalize changes
            let dn1 = if ibl < surface.itran {
                dctau / 10.0 // Laminar: normalize by 10
            } else {
                dctau / surface.ctau[ibl].max(1e-10) // Turbulent: normalize by Ctau
            };
            let dn2 = dthet / surface.thet[ibl].max(1e-10);
            let dn3 = ddstr / surface.dstr[ibl].max(1e-10);
            let dn4 = duedg.abs() / 0.25;

            // Accumulate RMS
            rmsbl += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

            // Check Ctau underrelaxation
            let rdn1 = rlx * dn1;
            if dn1.abs() > rmxbl.abs() {
                rmxbl = dn1;
                vmxbl = if ibl < surface.itran { 'n' } else { 'C' };
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn1 > dhi {
                rlx = dhi / dn1;
            }
            if rdn1 < dlo {
                rlx = dlo / dn1;
            }

            // Check Theta underrelaxation
            let rdn2 = rlx * dn2;
            if dn2.abs() > rmxbl.abs() {
                rmxbl = dn2;
                vmxbl = 'T';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn2 > dhi {
                rlx = dhi / dn2;
            }
            if rdn2 < dlo {
                rlx = dlo / dn2;
            }

            // Check Dstar underrelaxation
            let rdn3 = rlx * dn3;
            if dn3.abs() > rmxbl.abs() {
                rmxbl = dn3;
                vmxbl = 'D';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn3 > dhi {
                rlx = dhi / dn3;
            }
            if rdn3 < dlo {
                rlx = dlo / dn3;
            }

            // Check Ue underrelaxation
            let rdn4 = rlx * dn4;
            if dn4.abs() > rmxbl.abs() {
                rmxbl = duedg;
                vmxbl = 'U';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn4 > dhi {
                rlx = dhi / dn4;
            }
            if rdn4 < dlo {
                rlx = dlo / dn4;
            }
        }
    }

    // Finalize RMS
    let total_stations = (upper.nbl + lower.nbl) as f64;
    if total_stations > 0.0 {
        rmsbl = (rmsbl / (4.0 * total_stations)).sqrt();
    }

    // =========================================================================
    // Step 3: Apply underrelaxed updates
    // =========================================================================

    // Update upper surface
    for ibl in 1..upper.nbl {
        let iv = upper.isys[ibl];
        if iv * 3 + 2 >= solution.len() {
            continue;
        }

        let dctau = solution[3 * iv];
        let dthet = solution[3 * iv + 1];
        let dmass = solution[3 * iv + 2];
        let duedg = unew_upper[ibl] - upper.uedg[ibl];
        let ddstr = (dmass - upper.dstr[ibl] * duedg) / upper.uedg[ibl];

        upper.ctau[ibl] += rlx * dctau;
        upper.thet[ibl] += rlx * dthet;
        upper.dstr[ibl] += rlx * ddstr;
        upper.uedg[ibl] += rlx * duedg;

        // Limit Ctau for turbulent region
        if ibl >= upper.itran {
            upper.ctau[ibl] = upper.ctau[ibl].min(0.25);
        }

        // Limit shape factor using DSLIM
        let hklim = if ibl <= upper.iblte { 1.02 } else { 1.00005 };
        let msq = upper.uedg[ibl].powi(2) * hstinv / (gm1 * (1.0 - 0.5 * upper.uedg[ibl].powi(2) * hstinv));
        dslim(&mut upper.dstr[ibl], upper.thet[ibl], upper.uedg[ibl], msq, hklim);

        // Update mass defect (nonlinear)
        upper.mass[ibl] = upper.dstr[ibl] * upper.uedg[ibl];
    }

    // Fix negative Ue islands on upper surface
    for ibl in 2..upper.iblte {
        if upper.uedg[ibl - 1] > 0.0 && upper.uedg[ibl] <= 0.0 {
            upper.uedg[ibl] = upper.uedg[ibl - 1];
            upper.mass[ibl] = upper.dstr[ibl] * upper.uedg[ibl];
        }
    }

    // Update lower surface
    for ibl in 1..lower.nbl {
        let iv = lower.isys[ibl];
        if iv * 3 + 2 >= solution.len() {
            continue;
        }

        let dctau = solution[3 * iv];
        let dthet = solution[3 * iv + 1];
        let dmass = solution[3 * iv + 2];
        let duedg = unew_lower[ibl] - lower.uedg[ibl];
        let ddstr = (dmass - lower.dstr[ibl] * duedg) / lower.uedg[ibl];

        lower.ctau[ibl] += rlx * dctau;
        lower.thet[ibl] += rlx * dthet;
        lower.dstr[ibl] += rlx * ddstr;
        lower.uedg[ibl] += rlx * duedg;

        // Limit Ctau for turbulent region
        if ibl >= lower.itran {
            lower.ctau[ibl] = lower.ctau[ibl].min(0.25);
        }

        // Limit shape factor using DSLIM
        let hklim = if ibl <= lower.iblte { 1.02 } else { 1.00005 };
        let msq = lower.uedg[ibl].powi(2) * hstinv / (gm1 * (1.0 - 0.5 * lower.uedg[ibl].powi(2) * hstinv));
        dslim(&mut lower.dstr[ibl], lower.thet[ibl], lower.uedg[ibl], msq, hklim);

        // Update mass defect (nonlinear)
        lower.mass[ibl] = lower.dstr[ibl] * lower.uedg[ibl];
    }

    // Fix negative Ue islands on lower surface
    for ibl in 2..lower.iblte {
        if lower.uedg[ibl - 1] > 0.0 && lower.uedg[ibl] <= 0.0 {
            lower.uedg[ibl] = lower.uedg[ibl - 1];
            lower.mass[ibl] = lower.dstr[ibl] * lower.uedg[ibl];
        }
    }

    // =========================================================================
    // Step 4: Equate upper wake arrays to lower wake arrays
    // =========================================================================
    let nwake = lower.nbl.saturating_sub(lower.iblte);
    for kbl in 1..=nwake {
        let ibl_upper = upper.iblte + kbl;
        let ibl_lower = lower.iblte + kbl;
        if ibl_upper < upper.ctau.len() && ibl_lower < lower.ctau.len() {
            upper.ctau[ibl_upper] = lower.ctau[ibl_lower];
            upper.thet[ibl_upper] = lower.thet[ibl_lower];
            upper.dstr[ibl_upper] = lower.dstr[ibl_lower];
            upper.uedg[ibl_upper] = lower.uedg[ibl_lower];
            upper.mass[ibl_upper] = lower.mass[ibl_lower];
        }
    }

    // Store results
    result.rms_bl = rmsbl;
    result.rmx_bl = rmxbl;
    result.vmx_bl = vmxbl;
    result.imx_bl = imxbl;
    result.ismx_bl = ismxbl;
    result.rlx = rlx;

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_newton_system_creation() {
        let sys = BLNewtonSystem::new(100);
        assert_eq!(sys.n_sys, 100);
        assert_eq!(sys.va.len(), 100);
        assert_eq!(sys.vb.len(), 100);
        assert_eq!(sys.vdel.len(), 100);
    }

    #[test]
    fn test_station_state_default() {
        let state = BLStationState::default();
        assert_eq!(state.theta, 0.0);
        assert_eq!(state.cf, 0.0);
    }

    // ========================================================================
    // BLGlobalParams Tests
    // ========================================================================

    #[test]
    fn test_global_params_incompressible() {
        let params = BLGlobalParams::incompressible(1e6);

        // At M=0, TKBL should be 0
        assert_eq!(params.tk, 0.0);

        // RST = 1.0 at M=0
        assert_eq!(params.rst, 1.0);

        // HSTINV = 0 at M=0
        assert_eq!(params.hstinv, 0.0);

        // REYBL should equal REINF at M=0 (with some correction factor)
        assert_relative_eq!(params.reybl, 1e6, epsilon = 1.0);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (assumed HVRAT=0.35; XFOIL's analysis path leaves HVRAT=0) — regenerate from the M=0.3 coverage case"]
    fn test_global_params_compressible() {
        let params = BLGlobalParams::new(0.5, 1e6, 1.4);

        // Reference values from Fortran test
        assert_relative_eq!(params.tk, 0.1547005177, epsilon = 1e-6);
        assert_relative_eq!(params.rst, 1.129726171, epsilon = 1e-6);
        assert_relative_eq!(params.hstinv, 0.09523809701, epsilon = 1e-6);
        assert_relative_eq!(params.reybl, 963411.5, epsilon = 10.0);
    }

    // ========================================================================
    // BLPRV Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blprv_incompressible() {
        // Test case: M=0, Re=1e6
        // XFOIL reference values from Fortran test
        let params = BLGlobalParams::incompressible(1e6);
        let mut state = BLStationState::default();

        // Inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // Check primary variables
        assert_eq!(state.x, 0.1);
        assert_eq!(state.ampl, 3.0);
        assert_eq!(state.ctau, 0.015);
        assert_eq!(state.theta, 0.002);
        assert_eq!(state.dstar, 0.005);
        assert_eq!(state.dw, 0.0);

        // At M=0, U2 = Uei (no transformation)
        assert_relative_eq!(state.u, 1.2, epsilon = 1e-10);
        assert_relative_eq!(state.u_uei, 1.0, epsilon = 1e-10);
        // Note: u_ms is the sensitivity d(U2)/d(M²). Even at M=0 this is non-zero because
        // COMSET's TKL_MSQ = 1/(1+BETA)² = 0.25 at M=0 (BETA = 1).
        // U2_MS = (U2*UEI² - UEI) * TKBL_MS = (1.2*1.44 - 1.2) * 0.25 = 0.132
        assert_relative_eq!(state.u_ms, 0.132, epsilon = 1e-6);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (derived with TKBL = 1/beta - 1, not COMSET's TKLAM) — regenerate from the M=0.3 coverage case"]
    fn test_blprv_compressible() {
        // Test case: M=0.5, Re=1e6
        // XFOIL reference values from Fortran test
        let params = BLGlobalParams::new(0.5, 1e6, 1.4);
        let mut state = BLStationState::default();

        // Same inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // XFOIL reference: U2 = 1.305093527
        assert_relative_eq!(state.u, 1.305093527, epsilon = 1e-5);
        // XFOIL reference: U2_UEI = 1.711017609
        assert_relative_eq!(state.u_uei, 1.711017609, epsilon = 1e-5);
        // XFOIL reference: U2_MS = 0.6728398800
        assert_relative_eq!(state.u_ms, 0.6728398800, epsilon = 1e-5);
    }

    // ========================================================================
    // BLKIN Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blkin_incompressible() {
        // Test case: M=0, Re=1e6
        let params = BLGlobalParams::incompressible(1e6);
        let mut state = BLStationState::default();

        // Run BLPRV first
        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.blkin(&params);

        // XFOIL reference values (M=0 case)
        // M2 = 0
        assert_relative_eq!(state.msq, 0.0, epsilon = 1e-10);

        // R2 = 1.0
        assert_relative_eq!(state.r, 1.0, epsilon = 1e-10);

        // H2 = D2/T2 = 0.005/0.002 = 2.5
        assert_relative_eq!(state.h, 2.5, epsilon = 1e-10);
        assert_relative_eq!(state.h_t, -1250.0, epsilon = 1e-6); // -H/T
        assert_relative_eq!(state.h_d, 500.0, epsilon = 1e-6); // 1/T

        // HK2 = 2.5 (same as H at M=0)
        assert_relative_eq!(state.hk, 2.5, epsilon = 1e-6);

        // RT2 = 2400 (R*U*T/V = 1*1.2*0.002/1e-6)
        assert_relative_eq!(state.rt, 2400.0, epsilon = 1.0);

        // RT2_U2 = 2000 (RT/U = 2400/1.2)
        assert_relative_eq!(state.rt_u, 2000.0, epsilon = 1.0);

        // RT2_T2 = 1200000 (RT/T = 2400/0.002)
        assert_relative_eq!(state.rt_t, 1200000.0, epsilon = 100.0);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (assumed HVRAT=0.35; XFOIL's analysis path leaves HVRAT=0) — regenerate from the M=0.3 coverage case"]
    fn test_blkin_compressible() {
        // Test case: M=0.5, Re=1e6
        let params = BLGlobalParams::new(0.5, 1e6, 1.4);
        let mut state = BLStationState::default();

        // Run BLPRV first
        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.blkin(&params);

        // XFOIL reference values (M=0.5 case)
        // M2 = 0.4413362443
        assert_relative_eq!(state.msq, 0.4413362443, epsilon = 1e-5);

        // R2 = 0.9143959880
        assert_relative_eq!(state.r, 0.9143959880, epsilon = 1e-5);

        // HK2 = 2.259336710
        assert_relative_eq!(state.hk, 2.259336710, epsilon = 1e-5);

        // RT2 = 2453.646240
        assert_relative_eq!(state.rt, 2453.646240, epsilon = 1.0);
    }

    #[test]
    fn test_blkin_shape_factor_derivatives() {
        // Verify shape factor derivatives are computed correctly
        let params = BLGlobalParams::incompressible(1e6);
        let mut state = BLStationState::default();

        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);
        state.blkin(&params);

        // H = D/T = 2.5
        // dH/dD = 1/T = 1/0.002 = 500
        // dH/dT = -D/T² = -0.005/0.002² = -1250
        assert_relative_eq!(state.h_d, 1.0 / 0.002, epsilon = 1e-10);
        assert_relative_eq!(state.h_t, -0.005 / 0.002_f64.powi(2), epsilon = 1e-10);
    }

    // ========================================================================
    // BLVAR Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blvar_laminar() {
        // Test case: Laminar BL
        // XFOIL reference values from Fortran test
        let params = BLGlobalParams::incompressible(1e6);
        let mut state = BLStationState::default();

        // Set up: HK=2.5, RT=2400, M=0, H=2.5, T=0.002, D=0.005, S=3.0
        // We need to set the kinematic variables directly since blprv/blkin
        // uses its own computations
        state.hk = 2.5;
        state.rt = 2400.0;
        state.msq = 0.0;
        state.h = 2.5;
        state.theta = 0.002;
        state.dstar = 0.005;
        state.ctau = 3.0; // Amplification for laminar

        // Set the derivatives (simplified - zero for test)
        state.hk_u = 0.0;
        state.hk_t = -1250.0;
        state.hk_d = 500.0;
        state.hk_ms = 0.0;
        state.rt_u = 2000.0;
        state.rt_t = 1200000.0;
        state.rt_ms = 0.0;
        state.rt_re = 0.0024;
        state.h_t = -1250.0;
        state.h_d = 500.0;
        state.msq_u = 0.0;
        state.msq_ms = 0.0;

        state.blvar(BLFlowType::Laminar, &params);

        // XFOIL reference values:
        // HC2 = 0.0 (M=0)
        assert_relative_eq!(state.hc, 0.0, epsilon = 1e-10);

        // HS2 = 1.584867239
        assert_relative_eq!(state.hs, 1.584867239, epsilon = 1e-5);

        // US2 = 0.1584867090
        assert_relative_eq!(state.us, 0.1584867090, epsilon = 1e-5);

        // CQ2 = 0.07772709429
        assert_relative_eq!(state.cq, 0.07772709429, epsilon = 1e-4);

        // CF2 = 0.2045119036e-3
        assert_relative_eq!(state.cf, 0.2045119036e-3, epsilon = 1e-7);

        // DI2 = 0.9419409616e-4
        assert_relative_eq!(state.di, 0.9419409616e-4, epsilon = 1e-7);

        // DE2 = 0.01359333377
        assert_relative_eq!(state.de, 0.01359333377, epsilon = 1e-6);
    }

    #[test]
    fn test_blvar_turbulent() {
        // Test case: Turbulent BL
        // XFOIL reference values from Fortran test
        let params = BLGlobalParams::incompressible(1e6);
        let mut state = BLStationState::default();

        // Set up: HK=1.4, RT=10000, M=0, H=1.4, T=0.005, D=0.007, S=0.015
        state.hk = 1.4;
        state.rt = 10000.0;
        state.msq = 0.0;
        state.h = 1.4;
        state.theta = 0.005;
        state.dstar = 0.007;
        state.ctau = 0.015;

        // Set derivatives (simplified)
        state.hk_u = 0.0;
        state.hk_t = -280.0;
        state.hk_d = 200.0;
        state.hk_ms = 0.0;
        state.rt_u = 0.0;
        state.rt_t = 2000000.0;
        state.rt_ms = 0.0;
        state.rt_re = 0.01;
        state.h_t = -280.0;
        state.h_d = 200.0;
        state.msq_u = 0.0;
        state.msq_ms = 0.0;
        state.hs_u = 0.0;
        state.hs_t = 0.0;
        state.hs_d = 0.0;
        state.hs_ms = 0.0;
        state.hs_re = 0.0;
        state.us_u = 0.0;
        state.us_t = 0.0;
        state.us_d = 0.0;
        state.us_ms = 0.0;
        state.us_re = 0.0;

        state.blvar(BLFlowType::Turbulent, &params);

        // XFOIL reference values:
        // HS2 = 1.755310297
        assert_relative_eq!(state.hs, 1.755310297, epsilon = 1e-4);

        // US2 = 0.5433103442
        assert_relative_eq!(state.us, 0.5433103442, epsilon = 1e-4);

        // CQ2 = 0.03632329032
        assert_relative_eq!(state.cq, 0.03632329032, epsilon = 1e-4);

        // CF2 = 0.2286923816e-2
        assert_relative_eq!(state.cf, 0.2286923816e-2, epsilon = 1e-5);

        // DI2 = 0.8065673755e-3 (turbulent DI is more complex, allow larger tolerance)
        assert_relative_eq!(state.di, 0.8065673755e-3, epsilon = 1e-4);

        // DE2 = 0.04425000027
        assert_relative_eq!(state.de, 0.04425000027, epsilon = 1e-5);
    }

    // ========================================================================
    // BLMID Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blmid_laminar() {
        // Test case: Laminar (ITYP=1)
        // XFOIL reference values from Fortran test
        // HKA = 2.4, RTA = 2200, MA = 0

        // Station 1
        let mut s1 = BLStationState::default();
        s1.hk = 2.3;
        s1.rt = 2000.0;
        s1.msq = 0.0;
        s1.hk_u = 0.0;
        s1.hk_t = -1150.0;
        s1.hk_d = 500.0;
        s1.hk_ms = -0.29;
        s1.rt_u = 1666.67;
        s1.rt_t = 1000000.0;
        s1.rt_ms = 0.0;
        s1.rt_re = 0.002;
        s1.msq_u = 0.0;
        s1.msq_ms = 0.0;

        // Station 2
        let mut s2 = BLStationState::default();
        s2.hk = 2.5;
        s2.rt = 2400.0;
        s2.msq = 0.0;
        s2.hk_u = 0.0;
        s2.hk_t = -1250.0;
        s2.hk_d = 500.0;
        s2.hk_ms = -0.29;
        s2.rt_u = 2000.0;
        s2.rt_t = 1200000.0;
        s2.rt_ms = 0.0;
        s2.rt_re = 0.0024;
        s2.msq_u = 0.0;
        s2.msq_ms = 0.0;

        let result = MidpointCf::compute(&s1, &s2, BLFlowType::Laminar, false);

        // XFOIL reference values
        assert_relative_eq!(result.cfm, 0.2577280102e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_u1, -0.9762444097e-4, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_t1, 0.1515112668, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_d1, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cfm_u2, -0.1171490949e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_t2, 0.1580646932, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_d2, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cfm_ms, 0.1059563219e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_re, -0.2577280056e-9, epsilon = 1e-12);
    }

    #[test]
    fn test_blmid_turbulent() {
        // Test case: Turbulent (ITYP=2)
        // XFOIL reference values from Fortran test
        // HKA = 1.375, RTA = 9000, MA = 0

        // Station 1
        let mut s1 = BLStationState::default();
        s1.hk = 1.35;
        s1.rt = 8000.0;
        s1.msq = 0.0;
        s1.hk_u = 0.0;
        s1.hk_t = -270.0;
        s1.hk_d = 200.0;
        s1.hk_ms = -0.29;
        s1.rt_u = 0.0;
        s1.rt_t = 1600000.0;
        s1.rt_ms = 0.0;
        s1.rt_re = 0.008;
        s1.msq_u = 0.0;
        s1.msq_ms = 0.0;

        // Station 2
        let mut s2 = BLStationState::default();
        s2.hk = 1.40;
        s2.rt = 10000.0;
        s2.msq = 0.0;
        s2.hk_u = 0.0;
        s2.hk_t = -280.0;
        s2.hk_d = 200.0;
        s2.hk_ms = -0.29;
        s2.rt_u = 0.0;
        s2.rt_t = 2000000.0;
        s2.rt_ms = 0.0;
        s2.rt_re = 0.01;
        s2.msq_u = 0.0;
        s2.msq_ms = 0.0;

        let result = MidpointCf::compute(&s1, &s2, BLFlowType::Turbulent, false);

        // XFOIL reference values (turbulent Cf used)
        assert_relative_eq!(result.cfm, 0.2450317144e-2, epsilon = 1e-5);
        assert_relative_eq!(result.cfm_t1, 0.5299984217, epsilon = 1e-2);
        assert_relative_eq!(result.cfm_d1, -0.4310033321, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_t2, 0.5385844707, epsilon = 1e-2);
        assert_relative_eq!(result.cfm_d2, -0.4310033321, epsilon = 1e-3);
    }

    #[test]
    fn test_blmid_wake() {
        // Test case: Wake (ITYP=3)
        // Wake should have zero skin friction
        let s1 = BLStationState::default();
        let s2 = BLStationState::default();

        let result = MidpointCf::compute(&s1, &s2, BLFlowType::Wake, false);

        assert_eq!(result.cfm, 0.0);
        assert_eq!(result.cfm_u1, 0.0);
        assert_eq!(result.cfm_t1, 0.0);
        assert_eq!(result.cfm_d1, 0.0);
    }

    // ========================================================================
    // BLDIF Tests - validate momentum equation against XFOIL Fortran
    // ========================================================================

    #[test]
    fn test_bldif_momentum_equation() {
        // Test case from Fortran test_bldif.f90
        // Typical turbulent BL stations

        // Station 1
        let mut s1 = BLStationState::default();
        s1.x = 0.10;
        s1.u = 1.15;
        s1.theta = 0.0018;
        s1.dstar = 0.0045;
        s1.dw = 0.0;
        s1.h = s1.dstar / s1.theta;
        s1.msq = 0.0;
        s1.hs = 1.75;
        s1.cf = 0.0025;
        s1.hk = 1.40;
        s1.rt = 2070.0;
        s1.de = 0.012;
        s1.ctau = 0.015;
        s1.cq = 0.04;
        s1.us = 0.5;
        s1.di = 0.0008;

        // Station 1 derivatives
        s1.h_t = -s1.h / s1.theta;
        s1.h_d = 1.0 / s1.theta;
        s1.msq_u = 0.0;
        s1.msq_ms = 0.0;
        s1.cf_t = 0.0;
        s1.cf_d = 0.0;
        s1.cf_u = 0.0;
        s1.cf_ms = 0.0;
        s1.cf_re = -s1.cf / 1e6;
        s1.hk_t = 0.0;
        s1.hk_d = 0.0;
        s1.hk_u = 0.0;
        s1.hk_ms = 0.0;
        s1.hs_t = 0.0;
        s1.hs_d = 0.0;
        s1.hs_u = 0.0;
        s1.hs_ms = 0.0;
        s1.hs_re = 0.0;
        s1.de_t = 0.0;
        s1.de_d = 0.0;
        s1.de_u = 0.0;
        s1.de_ms = 0.0;
        s1.us_t = 0.0;
        s1.us_d = 0.0;
        s1.us_u = 0.0;
        s1.us_ms = 0.0;
        s1.us_re = 0.0;
        s1.cq_t = 0.0;
        s1.cq_d = 0.0;
        s1.cq_u = 0.0;
        s1.cq_ms = 0.0;
        s1.cq_re = 0.0;
        s1.di_t = 0.0;
        s1.di_d = 0.0;
        s1.di_u = 0.0;
        s1.di_s = 0.0;
        s1.di_ms = 0.0;
        s1.di_re = 0.0;
        s1.hc = 0.0;
        s1.hc_t = 0.0;
        s1.hc_d = 0.0;
        s1.hc_u = 0.0;
        s1.hc_ms = 0.0;
        s1.rt_t = 0.0;
        s1.rt_u = 0.0;
        s1.rt_ms = 0.0;
        s1.rt_re = 0.0;

        // Station 2
        let mut s2 = BLStationState::default();
        s2.x = 0.12;
        s2.u = 1.12;
        s2.theta = 0.0022;
        s2.dstar = 0.0052;
        s2.dw = 0.0;
        s2.h = s2.dstar / s2.theta;
        s2.msq = 0.0;
        s2.hs = 1.76;
        s2.cf = 0.0024;
        s2.hk = 1.38;
        s2.rt = 2460.0;
        s2.de = 0.015;
        s2.ctau = 0.014;
        s2.cq = 0.038;
        s2.us = 0.52;
        s2.di = 0.00075;

        // Station 2 derivatives
        s2.h_t = -s2.h / s2.theta;
        s2.h_d = 1.0 / s2.theta;
        s2.msq_u = 0.0;
        s2.msq_ms = 0.0;
        s2.cf_t = 0.0;
        s2.cf_d = 0.0;
        s2.cf_u = 0.0;
        s2.cf_ms = 0.0;
        s2.cf_re = -s2.cf / 1e6;
        s2.hk_t = 0.0;
        s2.hk_d = 0.0;
        s2.hk_u = 0.0;
        s2.hk_ms = 0.0;
        s2.hs_t = 0.0;
        s2.hs_d = 0.0;
        s2.hs_u = 0.0;
        s2.hs_ms = 0.0;
        s2.hs_re = 0.0;
        s2.de_t = 0.0;
        s2.de_d = 0.0;
        s2.de_u = 0.0;
        s2.de_ms = 0.0;
        s2.us_t = 0.0;
        s2.us_d = 0.0;
        s2.us_u = 0.0;
        s2.us_ms = 0.0;
        s2.us_re = 0.0;
        s2.cq_t = 0.0;
        s2.cq_d = 0.0;
        s2.cq_u = 0.0;
        s2.cq_ms = 0.0;
        s2.cq_re = 0.0;
        s2.di_t = 0.0;
        s2.di_d = 0.0;
        s2.di_u = 0.0;
        s2.di_s = 0.0;
        s2.di_ms = 0.0;
        s2.di_re = 0.0;
        s2.hc = 0.0;
        s2.hc_t = 0.0;
        s2.hc_d = 0.0;
        s2.hc_u = 0.0;
        s2.hc_ms = 0.0;
        s2.rt_t = 0.0;
        s2.rt_u = 0.0;
        s2.rt_ms = 0.0;
        s2.rt_re = 0.0;

        // Create CFM (simplified average)
        let mut cfm = MidpointCf::default();
        cfm.cfm = 0.5 * (s1.cf + s2.cf);
        cfm.cfm_re = 0.5 * (s1.cf_re + s2.cf_re);

        // Run BLDIF
        let mut sys = BLLocalSystem::default();
        sys.bldif(&s1, &s2, &cfm, BLFlowType::Turbulent, false, 9.0, 0);

        // Check momentum equation residual (row 2)
        // VSREZ[2] = -0.7123274356e-1 from Fortran
        assert_relative_eq!(sys.vsrez[1], -0.07123274356, epsilon = 1e-4);

        // Check Jacobian entries for row 2 (momentum)
        // Note: indexing is [row][col] where col is 0=S, 1=T, 2=D, 3=U, 4=X
        // VS1[2,1] (dT1) = -0.5339051514e+3
        assert_relative_eq!(sys.vs1[1][1], -533.9, epsilon = 1.0);

        // VS1[2,2] (dD1) = -0.7342562675e+1
        assert_relative_eq!(sys.vs1[1][2], -7.34, epsilon = 0.1);

        // VS1[2,3] (dU1) = -0.3853754759e+1
        assert_relative_eq!(sys.vs1[1][3], -3.85, epsilon = 0.1);

        // VS2[2,1] (dT2) = 0.4716367493e+3
        assert_relative_eq!(sys.vs2[1][1], 471.6, epsilon = 1.0);

        // VS2[2,2] (dD2) = -0.6007551670e+1
        assert_relative_eq!(sys.vs2[1][2], -6.0, epsilon = 0.1);

        // VS2[2,3] (dU2) = 0.3956980467e+1
        assert_relative_eq!(sys.vs2[1][3], 3.96, epsilon = 0.1);
    }

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

        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create stations with low Rtheta (below critical for HK ~2.5)
        // Critical log10(RT) for HK=2.5 is about 2.86, so RT ~720 is below critical
        let mut s1 = BLStationState::default();
        s1.x = 0.02;
        s1.u = 1.15;
        s1.theta = 0.0002;
        s1.dstar = 0.0005;
        s1.ampl = 0.0;
        s1.blprv(s1.x, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.u, &params);
        s1.blkin(&params);
        // s1.rt should be around 200 (below critical)

        let mut s2 = BLStationState::default();
        s2.x = 0.04;
        s2.u = 1.12;
        s2.theta = 0.0004;
        s2.dstar = 0.0010;
        s2.ampl = 0.0;
        s2.blprv(s2.x, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.u, &params);
        s2.blkin(&params);

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

        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create stations with moderate Rtheta (above critical, amplification active)
        let mut s1 = BLStationState::default();
        s1.x = 0.08;
        s1.u = 1.10;
        s1.theta = 0.001;
        s1.dstar = 0.0025;
        s1.ampl = 2.0; // Starting amplification
        s1.blprv(s1.x, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.u, &params);
        s1.blkin(&params);

        let mut s2 = BLStationState::default();
        s2.x = 0.10;
        s2.u = 1.08;
        s2.theta = 0.0012;
        s2.dstar = 0.003;
        s2.ampl = 2.5;
        s2.blprv(s2.x, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.u, &params);
        s2.blkin(&params);

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

        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create stations with high Rtheta and amplification close to Ncrit
        let mut s1 = BLStationState::default();
        s1.x = 0.15;
        s1.u = 1.05;
        s1.theta = 0.002;
        s1.dstar = 0.005;
        s1.ampl = 8.0; // Close to Ncrit=9
        s1.blprv(s1.x, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.u, &params);
        s1.blkin(&params);

        let mut s2 = BLStationState::default();
        s2.x = 0.25;
        s2.u = 1.02;
        s2.theta = 0.003;
        s2.dstar = 0.0075;
        s2.ampl = 12.0; // Well above Ncrit
        s2.blprv(s2.x, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.u, &params);
        s2.blkin(&params);

        let result = trchek(&s1, &s2, 8.0, 9.0, 1e6, &params);

        // Should return FreeTransition
        match result {
            TransitionResult::FreeTransition { location, ampl2 } => {
                // Transition should occur between X1 and X2
                assert!(
                    location.xt >= s1.x && location.xt <= s2.x,
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

        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create stations where natural transition wouldn't occur
        let mut s1 = BLStationState::default();
        s1.x = 0.10;
        s1.u = 1.10;
        s1.theta = 0.001;
        s1.dstar = 0.0025;
        s1.ampl = 2.0;
        s1.blprv(s1.x, s1.ampl, 0.0, s1.theta, s1.dstar, 0.0, s1.u, &params);
        s1.blkin(&params);

        let mut s2 = BLStationState::default();
        s2.x = 0.20;
        s2.u = 1.05;
        s2.theta = 0.0015;
        s2.dstar = 0.0038;
        s2.ampl = 3.0;
        s2.blprv(s2.x, s2.ampl, 0.0, s2.theta, s2.dstar, 0.0, s2.u, &params);
        s2.blkin(&params);

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

    // ========================================================================
    // SETBL Tests - validate global Newton system assembly
    // ========================================================================

    #[test]
    fn test_setbl_assemble_laminar() {
        // Test assembling a simple laminar BL surface
        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create a surface with 5 stations (below transition threshold)
        let mut surface = BLSurfaceState::new(5);

        // Set up stations with typical laminar BL development
        let x_vals = [0.01, 0.05, 0.10, 0.15, 0.20];
        let u_vals = [1.20, 1.15, 1.10, 1.08, 1.05];
        let theta_vals = [0.0001, 0.0003, 0.0005, 0.0007, 0.0009];
        let dstar_vals = [0.00025, 0.00075, 0.00125, 0.00175, 0.00225];

        for i in 0..5 {
            let mut s = BLStationState::default();
            s.x = x_vals[i];
            s.u = u_vals[i];
            s.theta = theta_vals[i];
            s.dstar = dstar_vals[i];
            s.ampl = (i as f64) * 0.5; // Low amplification
            s.blprv(s.x, s.ampl, 0.0, s.theta, s.dstar, 0.0, s.u, &params);
            s.blkin(&params);
            s.blvar(BLFlowType::Laminar, &params);

            surface.stations.push(s);
            surface.isys.push(i);
        }
        surface.nbl = 5;
        surface.iblte = 5; // No TE yet (all on surface)

        // Create Newton system
        let mut newton_sys = BLNewtonSystem::new(5);

        // Assemble the surface (no transition expected at low N)
        let itran = newton_sys.assemble_surface(&mut surface, &params, 9.0, 1e6);

        // Check that system was assembled
        // itran should be 0 (no transition) since amplification stays low
        assert_eq!(itran, 0, "Expected no transition for low amplification");

        // Check that VA and VB have been populated
        // Station 0 is similarity, so VB should be zero
        // For other stations, VB should be non-zero
        let sum_vb: f64 = newton_sys.vb[1].iter().flat_map(|r| r.iter()).sum();
        assert!(
            sum_vb.abs() > 0.0 || newton_sys.va[1][1][1].abs() > 0.0,
            "Newton system should have non-zero entries"
        );

        newton_sys.finalize_residual();
        assert!(newton_sys.rms_bl >= 0.0, "RMS residual should be non-negative");
    }

    #[test]
    fn test_blsolv_simple_system() {
        // Test BLSOLV with a simple diagonal-dominant system
        let n = 3;
        let mut sys = BLNewtonSystem::new(n);

        // Set up a simple diagonal-dominant system
        // Each station has 3 equations
        for iv in 0..n {
            // Diagonal block (identity-like)
            sys.va[iv][0][0] = 2.0;
            sys.va[iv][0][1] = 0.1;
            sys.va[iv][1][0] = 0.1;
            sys.va[iv][1][1] = 2.0;

            // Dense coupling (small off-diagonal)
            sys.vm[(0, iv)] = 0.05;
            sys.vm[(1, iv)] = 0.05;
            sys.vm[(2, iv)] = 2.0; // Diagonal for row 3

            // Sub-diagonal (if not first station)
            if iv > 0 {
                sys.vb[iv][0][0] = -0.2;
                sys.vb[iv][0][1] = 0.0;
                sys.vb[iv][1][0] = 0.0;
                sys.vb[iv][1][1] = -0.2;
            }

            // RHS
            sys.vdel[iv][0][0] = 1.0;
            sys.vdel[iv][1][0] = 1.0;
            sys.vdel[iv][2][0] = 1.0;
        }

        // Solve
        let solution = sys.solve();

        // Check that we got a solution of the right size
        assert_eq!(solution.len(), 3 * n, "Solution should have 3*n elements");

        // For a diagonal-dominant system, solution should be finite
        for i in 0..solution.len() {
            assert!(solution[i].is_finite(), "Solution element {} should be finite", i);
        }
    }

    // ========================================================================
    // UPDATE Tests - validate Newton step application
    // ========================================================================

    #[test]
    fn test_dslim_no_change_needed() {
        // Test DSLIM when Hk is already above limit
        let mut dstr = 0.005;
        let thet = 0.002;
        let uedg = 1.0;
        let msq = 0.0;
        let hklim = 1.02;

        // H = 2.5, Hk = 2.5 at M=0 (well above 1.02)
        dslim(&mut dstr, thet, uedg, msq, hklim);

        // Should not change
        assert_relative_eq!(dstr, 0.005, epsilon = 1e-10);
    }

    #[test]
    fn test_dslim_adjustment_needed() {
        // Test DSLIM when Hk is below limit
        let mut dstr = 0.00102; // H = 1.02/2 = 0.51 initially
        let thet = 0.001;
        let uedg = 1.0;
        let msq = 0.0;
        let hklim = 1.5; // Hk = H at M=0, so we need H >= 1.5

        // H = 1.02, Hk = 1.02 at M=0 (below 1.5)
        dslim(&mut dstr, thet, uedg, msq, hklim);

        // Should increase dstr to raise Hk to at least hklim
        let new_h = dstr / thet;
        assert!(new_h >= 1.5, "H should be raised to at least hklim");
    }

    #[test]
    fn test_update_basic() {
        // Basic test for UPDATE function
        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        // Create minimal upper and lower surfaces
        let mut upper = BLSurfaceState::new(3);
        let mut lower = BLSurfaceState::new(3);

        // Initialize with 3 stations each
        for i in 0..3 {
            upper.uinv.push(1.0 - 0.1 * i as f64);
            upper.uedg.push(1.0 - 0.1 * i as f64);
            upper.mass.push(0.001 * (i + 1) as f64);
            upper.thet.push(0.0005 * (i + 1) as f64);
            upper.dstr.push(0.001 * (i + 1) as f64);
            upper.ctau.push(0.01);
            upper.ipan.push(i);
            upper.isys.push(i);

            lower.uinv.push(1.0 - 0.1 * i as f64);
            lower.uedg.push(1.0 - 0.1 * i as f64);
            lower.mass.push(0.001 * (i + 1) as f64);
            lower.thet.push(0.0005 * (i + 1) as f64);
            lower.dstr.push(0.001 * (i + 1) as f64);
            lower.ctau.push(0.01);
            lower.ipan.push(i);
            lower.isys.push(3 + i);
        }
        upper.nbl = 3;
        upper.itran = 2;
        upper.iblte = 3;
        lower.nbl = 3;
        lower.itran = 2;
        lower.iblte = 3;

        // Create a small solution vector (6 stations * 3 vars = 18 elements)
        let solution = DVector::from_element(18, 0.0001);

        // Create dummy DIJ matrix (3x3 for simplicity)
        let dij = DMatrix::zeros(3, 3);

        // VTI signs
        let vti_upper = vec![1.0, 1.0, 1.0];
        let vti_lower = vec![-1.0, -1.0, -1.0];

        // Run update
        let result = update(&mut upper, &mut lower, &solution, &dij, &vti_upper, &vti_lower, &params);

        // Check that update produced reasonable results
        assert!(result.rms_bl.is_finite(), "RMS should be finite");
        assert!(
            result.rlx > 0.0 && result.rlx <= 1.0,
            "Relaxation factor should be in (0,1]"
        );

        // Check that variables were updated
        // With small solution vector and zero DIJ, changes should be small
        for i in 1..upper.nbl {
            assert!(upper.thet[i].is_finite(), "Theta should be finite");
            assert!(upper.dstr[i].is_finite(), "Dstar should be finite");
            assert!(upper.mass[i].is_finite(), "Mass should be finite");
        }
    }

    #[test]
    fn test_update_underrelaxation() {
        // Test that UPDATE applies underrelaxation for large changes
        let params = BLGlobalParams::new(0.0, 1e6, 1.4);

        let mut upper = BLSurfaceState::new(2);
        let mut lower = BLSurfaceState::new(2);

        // Initialize with 2 stations
        for i in 0..2 {
            upper.uinv.push(1.0);
            upper.uedg.push(1.0);
            upper.mass.push(0.001);
            upper.thet.push(0.001);
            upper.dstr.push(0.001);
            upper.ctau.push(0.01);
            upper.ipan.push(i);
            upper.isys.push(i);

            lower.uinv.push(1.0);
            lower.uedg.push(1.0);
            lower.mass.push(0.001);
            lower.thet.push(0.001);
            lower.dstr.push(0.001);
            lower.ctau.push(0.01);
            lower.ipan.push(i);
            lower.isys.push(2 + i);
        }
        upper.nbl = 2;
        upper.itran = 1;
        upper.iblte = 2;
        lower.nbl = 2;
        lower.itran = 1;
        lower.iblte = 2;

        // Create a LARGE solution vector that would require underrelaxation
        // dTheta = 0.1 (100x the original theta) should trigger underrelaxation
        let mut solution_data = vec![0.0; 12];
        solution_data[1] = 0.1; // Large dTheta for station 0
        solution_data[4] = 0.1; // Large dTheta for station 1
        let solution = DVector::from_vec(solution_data);

        let dij = DMatrix::zeros(2, 2);
        let vti_upper = vec![1.0, 1.0];
        let vti_lower = vec![-1.0, -1.0];

        let result = update(&mut upper, &mut lower, &solution, &dij, &vti_upper, &vti_lower, &params);

        // Relaxation factor should be less than 1.0 due to large changes
        assert!(
            result.rlx < 1.0,
            "Expected underrelaxation for large changes, got rlx={}",
            result.rlx
        );

        // The max change variable should be 'T' for theta
        assert_eq!(result.vmx_bl, 'T', "Expected max change in Theta");
    }
}
