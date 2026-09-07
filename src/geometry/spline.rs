//! Cubic spline interpolation
//!
//! This module provides natural cubic spline interpolation for airfoil
//! coordinate representation. All functions are pure with no side effects.
//!
//! Based on the spline routines from XFOIL (spline.f).

use crate::geometry::panel::solve_tridiagonal;

/// Compute cubic spline coefficients for data points
///
/// Given data points (s[i], x[i]), computes the spline derivative coefficients
/// xp[i] = dx/ds at each point using natural boundary conditions (zero second
/// derivative at endpoints).
///
/// # Arguments
/// * `x` - Data values
/// * `s` - Parameter values (typically arc length)
///
/// # Returns
/// Vector of spline derivative coefficients
///
/// # Panics
/// Panics if `x` and `s` have different lengths or fewer than 2 points
pub fn spline_derivatives(x: &[f64], s: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert_eq!(n, s.len(), "x and s must have same length");
    assert!(n >= 2, "Need at least 2 points for spline");

    if n == 2 {
        // Linear interpolation for 2 points
        let slope = (x[1] - x[0]) / (s[1] - s[0]);
        return vec![slope, slope];
    }

    // Set up tridiagonal system for natural spline
    // a[i] * xp[i-1] + b[i] * xp[i] + c[i] * xp[i+1] = d[i]
    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];

    // Natural boundary conditions: x''(s) = 0 at endpoints
    // First point
    let ds0 = s[1] - s[0];
    b[0] = 1.0;
    c[0] = 1.0;
    d[0] = 2.0 * (x[1] - x[0]) / ds0;

    // Interior points
    for i in 1..n - 1 {
        let dsm = s[i] - s[i - 1];
        let dsp = s[i + 1] - s[i];

        a[i] = dsp;
        b[i] = 2.0 * (dsm + dsp);
        c[i] = dsm;
        d[i] = 3.0 * ((x[i + 1] - x[i]) / dsp * dsm + (x[i] - x[i - 1]) / dsm * dsp);
    }

    // Last point
    let dsn1 = s[n - 1] - s[n - 2];
    a[n - 1] = 1.0;
    b[n - 1] = 1.0;
    d[n - 1] = 2.0 * (x[n - 1] - x[n - 2]) / dsn1;

    // Solve tridiagonal system (XFOIL's TRISOL: main diagonal first, then lower, upper, rhs)
    solve_tridiagonal(&mut b, &a, &mut c, &mut d);
    d
}

/// Evaluate spline at parameter value
///
/// # Arguments
/// * `ss` - Parameter value to evaluate at
/// * `x` - Data values
/// * `xp` - Spline derivative coefficients (from `spline()`)
/// * `s` - Parameter values
///
/// # Returns
/// Interpolated value at `ss`
pub fn spline_value(ss: f64, x: &[f64], xp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();

    // Find interval containing ss using binary search
    let i = find_interval(ss, s);
    let i = i.min(n - 2);

    // Cubic Hermite interpolation
    let ds = s[i + 1] - s[i];
    let t = (ss - s[i]) / ds;
    let t2 = t * t;
    let t3 = t2 * t;

    // Hermite basis functions
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;

    h00 * x[i] + h10 * ds * xp[i] + h01 * x[i + 1] + h11 * ds * xp[i + 1]
}

/// Evaluate first derivative of spline at parameter value
///
/// # Returns
/// dx/ds at `ss`
pub fn spline_slope(ss: f64, x: &[f64], xp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();

    let i = find_interval(ss, s);
    let i = i.min(n - 2);

    let ds = s[i + 1] - s[i];
    let t = (ss - s[i]) / ds;
    let t2 = t * t;

    // Derivatives of Hermite basis functions
    let dh00 = (6.0 * t2 - 6.0 * t) / ds;
    let dh10 = 3.0 * t2 - 4.0 * t + 1.0;
    let dh01 = (-6.0 * t2 + 6.0 * t) / ds;
    let dh11 = 3.0 * t2 - 2.0 * t;

    dh00 * x[i] + dh10 * xp[i] + dh01 * x[i + 1] + dh11 * xp[i + 1]
}

/// Evaluate second derivative of spline at parameter value
///
/// # Returns
/// d²x/ds² at `ss`
pub fn spline_second_derivative(ss: f64, x: &[f64], xp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();

    let i = find_interval(ss, s);
    let i = i.min(n - 2);

    let ds = s[i + 1] - s[i];
    let t = (ss - s[i]) / ds;

    // Second derivatives of Hermite basis functions
    let d2h00 = (12.0 * t - 6.0) / (ds * ds);
    let d2h10 = (6.0 * t - 4.0) / ds;
    let d2h01 = (-12.0 * t + 6.0) / (ds * ds);
    let d2h11 = (6.0 * t - 2.0) / ds;

    d2h00 * x[i] + d2h10 * xp[i] + d2h01 * x[i + 1] + d2h11 * xp[i + 1]
}

/// Find interval index containing parameter value using binary search
fn find_interval(ss: f64, s: &[f64]) -> usize {
    if ss <= s[0] {
        return 0;
    }
    if ss >= s[s.len() - 1] {
        return s.len() - 2;
    }

    let mut lo = 0;
    let mut hi = s.len() - 1;

    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if s[mid] > ss {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    lo
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_spline_linear_data() {
        // Spline through linear data should reproduce the line
        let s = vec![0.0, 1.0, 2.0, 3.0];
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let xp = spline_derivatives(&x, &s);

        // Evaluate at midpoints
        assert_relative_eq!(spline_value(0.5, &x, &xp, &s), 0.5, epsilon = 1e-10);
        assert_relative_eq!(spline_value(1.5, &x, &xp, &s), 1.5, epsilon = 1e-10);
        assert_relative_eq!(spline_value(2.5, &x, &xp, &s), 2.5, epsilon = 1e-10);
    }

    #[test]
    fn test_spline_quadratic_data() {
        // Test with quadratic data x = s^2
        let s = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let x: Vec<f64> = s.iter().map(|&si| si * si).collect();
        let xp = spline_derivatives(&x, &s);

        // Derivative should be approximately 2*s
        assert_relative_eq!(spline_slope(2.0, &x, &xp, &s), 4.0, epsilon = 0.1);
    }

    #[test]
    fn test_spline_two_points() {
        let s = vec![0.0, 1.0];
        let x = vec![0.0, 2.0];
        let xp = spline_derivatives(&x, &s);

        assert_relative_eq!(spline_value(0.5, &x, &xp, &s), 1.0, epsilon = 1e-10);
        assert_relative_eq!(spline_slope(0.5, &x, &xp, &s), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_find_interval() {
        let s = vec![0.0, 1.0, 2.0, 3.0, 4.0];

        assert_eq!(find_interval(-1.0, &s), 0);
        assert_eq!(find_interval(0.5, &s), 0);
        assert_eq!(find_interval(1.5, &s), 1);
        assert_eq!(find_interval(2.5, &s), 2);
        assert_eq!(find_interval(5.0, &s), 3);
    }

    #[test]
    fn test_spline_passes_through_data_points() {
        // Spline must pass exactly through all data points
        let s = vec![0.0, 0.5, 1.2, 2.0, 3.5, 4.0];
        let x = vec![1.0, 2.5, 1.8, 3.2, 2.1, 4.0];
        let xp = spline_derivatives(&x, &s);

        for i in 0..s.len() {
            assert_relative_eq!(spline_value(s[i], &x, &xp, &s), x[i], epsilon = 1e-12);
        }
    }

    #[test]
    fn test_spline_derivative_at_data_points() {
        // Derivative at data points should equal xp
        let s = vec![0.0, 1.0, 2.0, 3.0];
        let x = vec![0.0, 1.0, 0.0, 1.0];
        let xp = spline_derivatives(&x, &s);

        for i in 0..s.len() {
            assert_relative_eq!(spline_slope(s[i], &x, &xp, &s), xp[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_spline_sine_wave() {
        // Test with sine wave - common in airfoil upper/lower surface
        let n = 20;
        let s: Vec<f64> = (0..=n).map(|i| i as f64 * std::f64::consts::PI / n as f64).collect();
        let x: Vec<f64> = s.iter().map(|&si| si.sin()).collect();
        let xp = spline_derivatives(&x, &s);

        // Check interpolation at quarter points
        let s_test = std::f64::consts::PI / 4.0;
        let expected = s_test.sin();
        assert_relative_eq!(spline_value(s_test, &x, &xp, &s), expected, epsilon = 0.01);

        // Check derivative (should be cos(s))
        let expected_deriv = s_test.cos();
        assert_relative_eq!(spline_slope(s_test, &x, &xp, &s), expected_deriv, epsilon = 0.05);
    }

    #[test]
    fn test_spline_continuity_at_knots() {
        // First derivative should be continuous at knots (C1 continuity)
        let s = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let x = vec![0.0, 1.5, 1.0, 2.0, 0.5];
        let xp = spline_derivatives(&x, &s);

        // Check continuity at interior knots
        for i in 1..s.len() - 1 {
            let eps = 1e-8;
            let deriv_left = spline_slope(s[i] - eps, &x, &xp, &s);
            let deriv_right = spline_slope(s[i] + eps, &x, &xp, &s);
            assert_relative_eq!(deriv_left, deriv_right, epsilon = 1e-5);
        }
    }

    #[test]
    fn test_second_derivative() {
        // Test d2val for cubic data (should give constant second derivative)
        let s = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let x: Vec<f64> = s.iter().map(|&si| si * si * si).collect(); // x = s^3
        let xp = spline_derivatives(&x, &s);

        // For x = s^3, d²x/ds² = 6s
        // Note: cubic spline won't be exact for cubic data, but should be close
        let d2_mid = spline_second_derivative(2.0, &x, &xp, &s);
        assert_relative_eq!(d2_mid, 12.0, epsilon = 1.0); // 6 * 2 = 12
    }

    #[test]
    fn test_extrapolation_clamped() {
        // Extrapolation should use endpoint intervals
        let s = vec![0.0, 1.0, 2.0];
        let x = vec![0.0, 1.0, 0.0];
        let xp = spline_derivatives(&x, &s);

        // Values outside range should extrapolate using end intervals
        let val_before = spline_value(-0.5, &x, &xp, &s);
        let val_after = spline_value(2.5, &x, &xp, &s);

        // Should not be NaN or infinite
        assert!(val_before.is_finite());
        assert!(val_after.is_finite());
    }

    #[test]
    fn test_non_uniform_spacing() {
        // Test with non-uniformly spaced data (like cosine-spaced airfoil points)
        let n = 10;
        let s: Vec<f64> = (0..=n)
            .map(|i| 0.5 * (1.0 - (std::f64::consts::PI * i as f64 / n as f64).cos()))
            .collect();
        let x: Vec<f64> = s.iter().map(|&si| si * (1.0 - si)).collect(); // Parabola

        let xp = spline_derivatives(&x, &s);

        // Check passes through points
        for i in 0..s.len() {
            assert_relative_eq!(spline_value(s[i], &x, &xp, &s), x[i], epsilon = 1e-12);
        }
    }
}
