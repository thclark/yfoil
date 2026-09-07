//! Boundary layer closure relations
//!
//! This module implements the empirical closure relations used in XFOIL's
//! integral boundary layer formulation. These correlations relate the
//! integral boundary layer parameters (θ, δ*, H) to the skin friction,
//! dissipation, and energy shape factor.
//!
//! References:
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - Drela, M., Giles, M. "Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils"
//!   AIAA Journal, Oct. 1987

/// Result type for closure relations that includes sensitivities
#[derive(Debug, Clone, Copy)]
pub struct ClosureResult {
    /// Primary value
    pub val: f64,
    /// Sensitivity to kinematic shape factor Hk
    pub val_hk: f64,
    /// Sensitivity to momentum thickness Reynolds number Rt
    pub val_rt: f64,
    /// Sensitivity to Mach number squared
    pub val_msq: f64,
}

impl ClosureResult {
    /// Create a new closure result with only the value (zero sensitivities)
    pub fn value(val: f64) -> Self {
        Self {
            val,
            val_hk: 0.0,
            val_rt: 0.0,
            val_msq: 0.0,
        }
    }
}

// ============================================================================
// Shape Factor Conversions
// ============================================================================

/// Calculate kinematic shape factor Hk from H and Mach number
///
/// Hk = (H - 0.29*M²) / (1 + 0.113*M²)
///
/// This accounts for compressibility effects (from Whitfield)
pub fn hkin(h: f64, msq: f64) -> (f64, f64, f64) {
    let denom = 1.0 + 0.113 * msq;
    let hk = (h - 0.29 * msq) / denom;
    let hk_h = 1.0 / denom;
    let hk_msq = (-0.29 - 0.113 * hk) / denom;
    (hk, hk_h, hk_msq)
}

// ============================================================================
// Laminar Closure Relations
// ============================================================================

/// Laminar skin friction coefficient Cf (from Falkner-Skan)
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number Rθ
/// * `msq` - Mach number squared
pub fn cf_lam(hk: f64, rt: f64, _msq: f64) -> ClosureResult {
    let (cf, cf_hk) = if hk < 5.5 {
        let tmp = (5.5 - hk).powi(3) / (hk + 1.0);
        let cf = (0.0727 * tmp - 0.07) / rt;
        let cf_hk = (-0.0727 * tmp * 3.0 / (5.5 - hk) - 0.0727 * tmp / (hk + 1.0)) / rt;
        (cf, cf_hk)
    } else {
        let tmp = 1.0 - 1.0 / (hk - 4.5);
        let cf = (0.015 * tmp.powi(2) - 0.07) / rt;
        let cf_hk = (0.015 * tmp * 2.0 / (hk - 4.5).powi(2)) / rt;
        (cf, cf_hk)
    };

    let cf_rt = -cf / rt;

    ClosureResult {
        val: cf,
        val_hk: cf_hk,
        val_rt: cf_rt,
        val_msq: 0.0,
    }
}

/// Laminar energy shape factor H* correlation
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number (not used in laminar)
/// * `msq` - Mach number squared (not used in laminar)
pub fn hs_lam(hk: f64, _rt: f64, _msq: f64) -> ClosureResult {
    let (hs, hs_hk) = if hk < 4.35 {
        let tmp = hk - 4.35;
        let hs =
            0.0111 * tmp.powi(2) / (hk + 1.0) - 0.0278 * tmp.powi(3) / (hk + 1.0) + 1.528 - 0.0002 * (tmp * hk).powi(2);
        let hs_hk = 0.0111 * (2.0 * tmp - tmp.powi(2) / (hk + 1.0)) / (hk + 1.0)
            - 0.0278 * (3.0 * tmp.powi(2) - tmp.powi(3) / (hk + 1.0)) / (hk + 1.0)
            - 0.0002 * 2.0 * tmp * hk * (tmp + hk);
        (hs, hs_hk)
    } else {
        let hs2 = 0.015;
        let hs = hs2 * (hk - 4.35).powi(2) / hk + 1.528;
        let hs_hk = hs2 * 2.0 * (hk - 4.35) / hk - hs2 * (hk - 4.35).powi(2) / hk.powi(2);
        (hs, hs_hk)
    };

    ClosureResult {
        val: hs,
        val_hk: hs_hk,
        val_rt: 0.0,
        val_msq: 0.0,
    }
}

/// Laminar dissipation coefficient 2*CD/H* (from Falkner-Skan)
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number Rθ
pub fn di_lam(hk: f64, rt: f64) -> ClosureResult {
    let (di, di_hk) = if hk < 4.0 {
        let di = (0.00205 * (4.0 - hk).powf(5.5) + 0.207) / rt;
        let di_hk = (-0.00205 * 5.5 * (4.0 - hk).powf(4.5)) / rt;
        (di, di_hk)
    } else {
        let hkb = hk - 4.0;
        let den = 1.0 + 0.02 * hkb.powi(2);
        let di = (-0.0016 * hkb.powi(2) / den + 0.207) / rt;
        let di_hk = (-0.0016 * 2.0 * hkb * (1.0 / den - 0.02 * hkb.powi(2) / den.powi(2))) / rt;
        (di, di_hk)
    };

    let di_rt = -di / rt;

    ClosureResult {
        val: di,
        val_hk: di_hk,
        val_rt: di_rt,
        val_msq: 0.0,
    }
}

// ============================================================================
// Turbulent Closure Relations
// ============================================================================

/// Turbulent skin friction coefficient Cf (Coles correlation)
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number Rθ
/// * `msq` - Mach number squared
/// * `cffac` - Skin friction factor (typically 1.0)
pub fn cf_turb(hk: f64, rt: f64, msq: f64, cffac: f64) -> ClosureResult {
    const GAM: f64 = 1.4;
    let gm1 = GAM - 1.0;

    let fc = (1.0 + 0.5 * gm1 * msq).sqrt();
    let grt = (rt / fc).ln().max(3.0);

    let gex = -1.74 - 0.31 * hk;
    let arg = (-1.33 * hk).max(-20.0);
    let thk = (4.0 - hk / 0.875).tanh();

    let cfo = cffac * 0.3 * arg.exp() * (grt / 2.3026).powf(gex);
    let cf = (cfo + 1.1e-4 * (thk - 1.0)) / fc;

    let cf_hk = (-1.33 * cfo - 0.31 * (grt / 2.3026).ln() * cfo - 1.1e-4 * (1.0 - thk.powi(2)) / 0.875) / fc;
    let cf_rt = gex * cfo / (fc * grt) / rt;
    let cf_msq = gex * cfo / (fc * grt) * (-0.25 * gm1 / fc.powi(2)) - 0.25 * gm1 * cf / fc.powi(2);

    ClosureResult {
        val: cf,
        val_hk: cf_hk,
        val_rt: cf_rt,
        val_msq: cf_msq,
    }
}

/// Turbulent energy shape factor H* correlation
///
/// # Arguments
/// * `hk` - Kinematic shape factor
/// * `rt` - Momentum thickness Reynolds number Rθ
/// * `msq` - Mach number squared
pub fn hs_turb(hk: f64, rt: f64, msq: f64) -> ClosureResult {
    const HSMIN: f64 = 1.5;
    const DHSINF: f64 = 0.015;

    // Limited Rtheta dependence
    let (ho, ho_rt) = if rt > 400.0 {
        (3.0 + 400.0 / rt, -400.0 / rt.powi(2))
    } else {
        (4.0, 0.0)
    };

    let (rtz, rtz_rt) = if rt > 200.0 { (rt, 1.0) } else { (200.0, 0.0) };

    let (hs, hs_hk, hs_rt) = if hk < ho {
        // Attached branch
        let hr = (ho - hk) / (ho - 1.0);
        let hr_hk = -1.0 / (ho - 1.0);
        let hr_rt = (1.0 - hr) / (ho - 1.0) * ho_rt;

        let hs = (2.0 - HSMIN - 4.0 / rtz) * hr.powi(2) * 1.5 / (hk + 0.5) + HSMIN + 4.0 / rtz;
        let hs_hk = -(2.0 - HSMIN - 4.0 / rtz) * hr.powi(2) * 1.5 / (hk + 0.5).powi(2)
            + (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_hk;
        let hs_rt = (2.0 - HSMIN - 4.0 / rtz) * hr * 2.0 * 1.5 / (hk + 0.5) * hr_rt
            + (hr.powi(2) * 1.5 / (hk + 0.5) - 1.0) * 4.0 / rtz.powi(2) * rtz_rt;
        (hs, hs_hk, hs_rt)
    } else {
        // Separated branch
        let grt = rtz.ln();
        let hdif = hk - ho;
        let rtmp = hk - ho + 4.0 / grt;
        let htmp = 0.007 * grt / rtmp.powi(2) + DHSINF / hk;
        let htmp_hk = -0.014 * grt / rtmp.powi(3) - DHSINF / hk.powi(2);
        let htmp_rt = -0.014 * grt / rtmp.powi(3) * (-ho_rt - 4.0 / grt.powi(2) / rtz * rtz_rt)
            + 0.007 / rtmp.powi(2) / rtz * rtz_rt;

        let hs = hdif.powi(2) * htmp + HSMIN + 4.0 / rtz;
        let hs_hk = hdif * 2.0 * htmp + hdif.powi(2) * htmp_hk;
        let hs_rt = hdif.powi(2) * htmp_rt - 4.0 / rtz.powi(2) * rtz_rt + hdif * 2.0 * htmp * (-ho_rt);
        (hs, hs_hk, hs_rt)
    };

    // Whitfield's compressibility correction
    let fm = 1.0 + 0.014 * msq;
    let hs_final = (hs + 0.028 * msq) / fm;
    let hs_hk_final = hs_hk / fm;
    let hs_rt_final = hs_rt / fm;
    let hs_msq = (0.028 - 0.014 * hs_final) / fm;

    ClosureResult {
        val: hs_final,
        val_hk: hs_hk_final,
        val_rt: hs_rt_final,
        val_msq: hs_msq,
    }
}

/// Density shape parameter (from Whitfield)
pub fn hc_turb(hk: f64, msq: f64) -> (f64, f64, f64) {
    let hc = msq * (0.064 / (hk - 0.8) + 0.251);
    let hc_hk = msq * (-0.064 / (hk - 0.8).powi(2));
    let hc_msq = 0.064 / (hk - 0.8) + 0.251;
    (hc, hc_hk, hc_msq)
}

/// DILW (xblsys.f): laminar wake dissipation function 2*CD/H* and its Hk, Rt sensitivities.
pub fn dilw(hk: f64, rt: f64) -> ClosureResult {
    let msq = 0.0;
    let hs = hs_lam(hk, rt, msq);
    // Laminar wake dissipation function  ( 2 CD/H* )
    let rcd = 1.10 * ((1.0 - 1.0 / hk) * (1.0 - 1.0 / hk)) / hk;
    let rcd_hk = -1.10 * (1.0 - 1.0 / hk) * 2.0 / ((hk * hk) * hk) - rcd / hk;
    let di = 2.0 * rcd / (hs.val * rt);
    let di_hk = 2.0 * rcd_hk / (hs.val * rt) - (di / hs.val) * hs.val_hk;
    let di_rt = -di / rt - (di / hs.val) * hs.val_rt;
    ClosureResult {
        val: di,
        val_hk: di_hk,
        val_rt: di_rt,
        val_msq: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // ========================================================================
    // Shape Factor Tests
    // ========================================================================

    #[test]
    fn test_hkin_incompressible() {
        // At M=0, Hk = H
        let (hk, hk_h, hk_msq) = hkin(2.5, 0.0);
        assert_relative_eq!(hk, 2.5, epsilon = 1e-10);
        assert_relative_eq!(hk_h, 1.0, epsilon = 1e-10);
        // hk_msq = (-0.29 - 0.113*hk) / denom = -0.29 - 0.2825 = -0.5725
        assert_relative_eq!(hk_msq, -0.5725, epsilon = 1e-10);
    }

    #[test]
    fn test_hkin_compressible() {
        // At M=0.5 (M²=0.25), Hk should be reduced
        let (hk, _, _) = hkin(2.5, 0.25);
        assert!(hk < 2.5);
        assert!(hk > 2.0);
    }

    // ========================================================================
    // Laminar Closure Tests
    // ========================================================================

    #[test]
    fn test_cf_lam_blasius() {
        // For Blasius flow, Hk ≈ 2.59, Cf ≈ 0.664/√Re_x
        // At Rθ = 1000, Cf ≈ 0.664/√Rex ≈ 0.0021 for appropriate Rex
        let result = cf_lam(2.59, 1000.0, 0.0);
        assert!(result.val > 0.0);
        assert!(result.val < 0.01);
        // Cf should decrease with increasing Rt
        assert!(result.val_rt < 0.0);
    }

    #[test]
    fn test_cf_lam_separated() {
        // For separated flow (high Hk), Cf should be small or negative
        let result = cf_lam(6.0, 1000.0, 0.0);
        assert!(result.val < 0.0); // Negative Cf indicates separation
    }

    #[test]
    fn test_hs_lam_attached() {
        // For attached laminar flow, H* should be around 1.5-1.6
        let result = hs_lam(2.59, 1000.0, 0.0);
        assert!(result.val > 1.4);
        assert!(result.val < 1.8);
    }

    #[test]
    fn test_di_lam_positive() {
        // Dissipation should always be positive for physical flows
        let result = di_lam(2.59, 1000.0);
        assert!(result.val > 0.0);
    }

    // ========================================================================
    // Turbulent Closure Tests
    // ========================================================================

    #[test]
    fn test_cf_turb_attached() {
        // Turbulent Cf for attached flow (Hk ≈ 1.3-1.5)
        let result = cf_turb(1.4, 10000.0, 0.0, 1.0);
        assert!(result.val > 0.002);
        assert!(result.val < 0.01);
    }

    #[test]
    fn test_cf_turb_reynolds_effect() {
        // Cf should decrease with increasing Reynolds number
        let cf1 = cf_turb(1.4, 10000.0, 0.0, 1.0);
        let cf2 = cf_turb(1.4, 100000.0, 0.0, 1.0);
        assert!(cf2.val < cf1.val);
    }

    #[test]
    fn test_hs_turb_range() {
        // Turbulent H* should be > 1.5 (minimum) and typically < 2.5
        let result = hs_turb(1.4, 10000.0, 0.0);
        assert!(result.val >= 1.5);
        assert!(result.val < 3.0);
    }
}
