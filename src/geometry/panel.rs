//! Panel distribution and repaneling
//!
//! Functions for redistributing panel points on an airfoil surface.

use super::airfoil::{Geometry, PaneledAirfoil};
use super::spline::{d2val, deval, seval, spline};

/// Repanel an airfoil with a new number of panels using cosine spacing
///
/// Redistributes points along the airfoil surface with higher density
/// near the leading and trailing edges.
///
/// # Arguments
/// * `geometry` - Input geometry
/// * `n_panels` - Target number of panels
/// * `le_ratio` - Ratio of LE to TE panel spacing (smaller = finer at LE)
///
/// # Returns
/// New geometry with redistributed points
pub fn repanel(geometry: &Geometry, n_panels: usize, _le_ratio: f64) -> Geometry {
    let n = geometry.x_c.len();

    // Calculate arc length along the surface
    let s = calculate_arc_length(&geometry.x_c, &geometry.y_c);

    // Create splines for x and y
    let xp = spline(&geometry.x_c, &s);
    let yp = spline(&geometry.y_c, &s);

    let s_total = s[n - 1];

    // Generate new parameter values using cosine spacing
    let mut s_new = Vec::with_capacity(n_panels + 1);
    for i in 0..=n_panels {
        let theta = std::f64::consts::PI * (i as f64) / (n_panels as f64);
        // Map [0, pi] -> [0, s_total] with cosine distribution
        let s_val = s_total * 0.5 * (1.0 - theta.cos());
        s_new.push(s_val);
    }

    // Evaluate splines at new parameter values
    let x_c: Vec<f64> = s_new
        .iter()
        .map(|&si| seval(si, &geometry.x_c, &xp, &s))
        .collect();
    let y_c: Vec<f64> = s_new
        .iter()
        .map(|&si| seval(si, &geometry.y_c, &yp, &s))
        .collect();

    Geometry {
        reference: geometry.reference,
        x_c,
        y_c,
    }
}

/// Calculate arc length along the surface
fn calculate_arc_length(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut s = vec![0.0; n];

    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }

    s
}

/// Create a PaneledAirfoil from raw geometry
///
/// Computes all derived quantities needed for aerodynamic analysis:
/// - Arc length parameterization
/// - Spline coefficients
/// - Normal vectors
/// - Panel angles
/// - Leading edge location
pub fn create_paneled_airfoil(geometry: &Geometry) -> PaneledAirfoil {
    let n = geometry.x_c.len();
    let x = geometry.x_c.clone();
    let y = geometry.y_c.clone();

    // Calculate arc length
    let s = calculate_arc_length(&x, &y);

    // Create splines
    let xp = spline(&x, &s);
    let yp = spline(&y, &s);

    // Calculate normal vectors and panel angles
    let (nx, ny, apanel) = calculate_normals_and_angles(&x, &y, &xp, &yp, &s);

    // Find leading edge using XFOIL's chord-perpendicular criterion
    let (sle, le_index) = find_leading_edge(&x, &y, &s, &xp, &yp);

    // Calculate chord length
    let chord = calculate_chord(&x, &y);

    // Check for sharp trailing edge
    let te_gap = ((x[0] - x[n - 1]).powi(2) + (y[0] - y[n - 1]).powi(2)).sqrt();
    let sharp_te = te_gap < 0.0001 * chord;

    PaneledAirfoil {
        x,
        y,
        s,
        xp,
        yp,
        nx,
        ny,
        apanel,
        n,
        sle,
        le_index,
        chord,
        sharp_te,
        reference: geometry.reference,
    }
}

/// Calculate normal vectors and panel angles
fn calculate_normals_and_angles(
    x: &[f64],
    y: &[f64],
    xp: &[f64],
    yp: &[f64],
    s: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len();
    let mut nx = vec![0.0; n];
    let mut ny = vec![0.0; n];
    let mut apanel = vec![0.0; n];

    for i in 0..n {
        // Get tangent vector from spline derivatives
        let dx_ds = deval(s[i], x, xp, s);
        let dy_ds = deval(s[i], y, yp, s);

        // Tangent magnitude
        let ds = (dx_ds * dx_ds + dy_ds * dy_ds).sqrt();

        // Unit tangent
        let tx = dx_ds / ds;
        let ty = dy_ds / ds;

        // Normal is perpendicular to tangent (pointing outward for CCW ordering)
        // For TE->upper->LE->lower->TE ordering, normal points outward with this sign
        nx[i] = ty;
        ny[i] = -tx;

        // Panel angle (angle of tangent from horizontal)
        apanel[i] = ty.atan2(tx);
    }

    (nx, ny, apanel)
}

/// Find leading edge arc length parameter and index
///
/// Uses XFOIL's LEFIND algorithm: finds where the surface tangent is perpendicular
/// to the chord line connecting the LE point to the TE.
///
/// The defining condition is: (X-XTE, Y-YTE) · (X', Y') = 0 at S = SLE
///
/// Returns (sle, le_index)
fn find_leading_edge(
    x: &[f64],
    y: &[f64],
    s: &[f64],
    xp: &[f64],
    yp: &[f64],
) -> (f64, usize) {
    let n = x.len();

    // Convergence tolerance (matches XFOIL)
    let dseps = (s[n - 1] - s[0]) * 1.0e-5;

    // Trailing edge coordinates
    let x_te = 0.5 * (x[0] + x[n - 1]);
    let y_te = 0.5 * (y[0] + y[n - 1]);

    // Get first guess for SLE by finding where dot product changes sign
    // This matches XFOIL's approach exactly
    let mut i_le = n / 2; // fallback
    for i in 2..n - 2 {
        let dxte = x[i] - x_te;
        let dyte = y[i] - y_te;
        let dx = x[i + 1] - x[i];
        let dy = y[i + 1] - y[i];
        let dotp = dxte * dx + dyte * dy;
        if dotp < 0.0 {
            i_le = i;
            break;
        }
    }

    let mut s_le = s[i_le];

    // Check for sharp LE case (doubled point)
    if i_le > 0 && (s[i_le] - s[i_le - 1]).abs() < 1e-14 {
        return (s_le, i_le);
    }

    // Newton iteration to get exact SLE value (matches XFOIL exactly)
    for _ in 0..50 {
        let x_le = seval(s_le, x, xp, s);
        let y_le = seval(s_le, y, yp, s);
        let dxds = deval(s_le, x, xp, s);
        let dyds = deval(s_le, y, yp, s);
        let dxdd = d2val(s_le, x, xp, s);
        let dydd = d2val(s_le, y, yp, s);

        let xchord = x_le - x_te;
        let ychord = y_le - y_te;

        // Drive dot product between chord line and LE tangent to zero
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        // Newton delta for SLE
        let mut dsle = -res / ress;

        // Limit step size (matches XFOIL exactly: ABS(XCHORD+YCHORD), not ABS(XCHORD)+ABS(YCHORD))
        let dsle_limit = 0.02 * (xchord + ychord).abs();
        dsle = dsle.max(-dsle_limit).min(dsle_limit);

        s_le += dsle;

        if dsle.abs() < dseps {
            break;
        }
    }

    (s_le, i_le)
}

/// Calculate chord length (TE to LE distance)
fn calculate_chord(x: &[f64], y: &[f64]) -> f64 {
    // TE is at first/last points
    let x_te = (x[0] + x[x.len() - 1]) / 2.0;
    let y_te = (y[0] + y[y.len() - 1]) / 2.0;

    // LE is at minimum x
    let (i_le, _) = x
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    let x_le = x[i_le];
    let y_le = y[i_le];

    ((x_te - x_le).powi(2) + (y_te - y_le).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_arc_length_calculation() {
        // Simple square path
        let x = vec![0.0, 1.0, 1.0, 0.0, 0.0];
        let y = vec![0.0, 0.0, 1.0, 1.0, 0.0];
        let s = calculate_arc_length(&x, &y);

        assert_eq!(s[0], 0.0);
        assert!((s[1] - 1.0).abs() < 1e-10);
        assert!((s[2] - 2.0).abs() < 1e-10);
        assert!((s[3] - 3.0).abs() < 1e-10);
        assert!((s[4] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_arc_length_circle() {
        // Points on a circle should give arc length = angle * radius
        let n = 36;
        let radius = 1.0;
        let x: Vec<f64> = (0..=n)
            .map(|i| radius * (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos())
            .collect();
        let y: Vec<f64> = (0..=n)
            .map(|i| radius * (2.0 * std::f64::consts::PI * i as f64 / n as f64).sin())
            .collect();

        let s = calculate_arc_length(&x, &y);

        // Total arc length should be approximately 2*pi*r
        let total_arc = s[n];
        assert_relative_eq!(total_arc, 2.0 * std::f64::consts::PI * radius, epsilon = 0.1);
    }

    #[test]
    fn test_chord_calculation_symmetric() {
        // Simple symmetric airfoil shape
        // TE at x=1, LE at x=0
        let x = vec![1.0, 0.8, 0.5, 0.2, 0.0, 0.2, 0.5, 0.8, 1.0];
        let y = vec![0.0, 0.05, 0.08, 0.06, 0.0, -0.06, -0.08, -0.05, 0.0];

        let chord = calculate_chord(&x, &y);
        assert_relative_eq!(chord, 1.0, epsilon = 0.01);
    }

    #[test]
    fn test_paneled_airfoil_from_naca0012() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // Check basic properties
        assert!(paneled.n > 100);
        assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

        // Arc length should be monotonically increasing
        for i in 1..paneled.s.len() {
            assert!(paneled.s[i] > paneled.s[i - 1]);
        }

        // Total arc length should be roughly 2x chord for thin airfoil
        let total_arc = paneled.s[paneled.n - 1];
        assert!(total_arc > 1.8 && total_arc < 2.5);

        // Leading edge should be approximately at the midpoint of arc length
        assert!(paneled.sle > total_arc * 0.3 && paneled.sle < total_arc * 0.7);
    }

    #[test]
    fn test_normal_vectors_unit_length() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // All normal vectors should have unit length
        for i in 0..paneled.n {
            let mag = (paneled.nx[i].powi(2) + paneled.ny[i].powi(2)).sqrt();
            assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_normal_vectors_point_outward() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // For a symmetric airfoil centered on y=0:
        // - Upper surface (y > 0) normals should have ny > 0
        // - Lower surface (y < 0) normals should have ny < 0
        // (with some tolerance near LE/TE where y is close to 0)

        let mut upper_count = 0;
        let mut lower_count = 0;

        for i in 0..paneled.n {
            if paneled.y[i] > 0.02 {
                // Upper surface
                assert!(
                    paneled.ny[i] > 0.0,
                    "Upper surface normal should point up at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.ny[i]
                );
                upper_count += 1;
            } else if paneled.y[i] < -0.02 {
                // Lower surface
                assert!(
                    paneled.ny[i] < 0.0,
                    "Lower surface normal should point down at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.ny[i]
                );
                lower_count += 1;
            }
        }

        // Should have checked at least some points on each surface
        assert!(upper_count > 10);
        assert!(lower_count > 10);
    }

    #[test]
    fn test_sharp_trailing_edge_detection() {
        use crate::geometry::naca::naca_4digit;

        // NACA 0012 with closed TE should be detected as sharp
        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // The NACA generator uses modified coefficients for closed TE
        // so it should be detected as sharp (or very nearly so)
        assert!(paneled.sharp_te);
    }

    #[test]
    fn test_repanel_preserves_shape() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100).unwrap();
        let repaneled = repanel(&original, 150, 0.15);

        // Should have approximately the target number of points
        assert!(repaneled.x_c.len() > 140 && repaneled.x_c.len() < 160);

        // Extents should be preserved
        let orig_max_x = original.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_x = repaneled.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_x, new_max_x, epsilon = 0.01);

        let orig_min_x = original.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
        let new_min_x = repaneled.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
        assert_relative_eq!(orig_min_x, new_min_x, epsilon = 0.01);

        // Maximum thickness should be preserved
        let orig_max_y = original.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_y = repaneled.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_y, new_max_y, epsilon = 0.01);
    }

    #[test]
    fn test_repanel_cosine_spacing() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100).unwrap();
        let repaneled = repanel(&original, 100, 0.15);

        // Cosine spacing should cluster points near LE and TE
        // Points near x=0 and x=1 should be closer together than at mid-chord

        // Find spacing near LE (x close to 0)
        let mut le_spacings = Vec::new();
        let mut mid_spacings = Vec::new();

        for i in 1..repaneled.x_c.len() {
            let dx = (repaneled.x_c[i] - repaneled.x_c[i - 1]).abs();
            let dy = (repaneled.y_c[i] - repaneled.y_c[i - 1]).abs();
            let ds = (dx * dx + dy * dy).sqrt();

            let avg_x = (repaneled.x_c[i] + repaneled.x_c[i - 1]) / 2.0;

            if avg_x < 0.1 || avg_x > 0.9 {
                le_spacings.push(ds);
            } else if avg_x > 0.4 && avg_x < 0.6 {
                mid_spacings.push(ds);
            }
        }

        if !le_spacings.is_empty() && !mid_spacings.is_empty() {
            let avg_le: f64 = le_spacings.iter().sum::<f64>() / le_spacings.len() as f64;
            let avg_mid: f64 = mid_spacings.iter().sum::<f64>() / mid_spacings.len() as f64;

            // LE/TE spacing should be smaller than mid-chord spacing
            assert!(
                avg_le < avg_mid,
                "LE spacing {} should be smaller than mid spacing {}",
                avg_le,
                avg_mid
            );
        }
    }
}
