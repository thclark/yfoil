//! Integration tests for yfoil geometry module

use approx::assert_relative_eq;
use tempfile::NamedTempFile;

use yfoil::geometry::{
    create_paneled_airfoil, naca_4digit, naca_5digit, read_dat_file, read_geometry_from_file,
    repanel, write_dat_file, write_geometry_to_json,
};

/// Test NACA 0012 against the standard NACA formula from external reference.
///
/// Reference data: Standard NACA 0012 formula (Abbott & Von Doenhoff):
///   y = ±(t/0.2) * [0.2969√x - 0.126x - 0.3516x² + 0.2843x³ - 0.1015x⁴]
///
/// Our generator uses a closed-TE modification (0.1036 instead of 0.1015),
/// which produces a small known difference near the trailing edge:
///   Δy = 0.6 * (0.1036 - 0.1015) * x⁴ = 0.00126 * x⁴
///
/// Source: https://turbmodels.larc.nasa.gov/naca0012_val.html
#[test]
fn test_naca_0012_against_standard_formula() {
    // Reference data from standard NACA formula (0.1015 coefficient)
    let (_, reference) = read_dat_file("tests/fixtures/naca0012_reference.dat").unwrap();

    // Helper: compute NACA 0012 thickness using our closed-TE formula
    fn naca_0012_closed_te(x: f64) -> f64 {
        let t = 0.12;
        (t / 0.2)
            * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4))
    }

    // For each reference point, compute what our formula gives and compare
    for i in 0..reference.x_c.len() {
        let ref_x = reference.x_c[i];
        let ref_y = reference.y_c[i];

        // Our formula's expected y at this x (symmetric airfoil)
        let our_y = if ref_y >= 0.0 {
            naca_0012_closed_te(ref_x)
        } else {
            -naca_0012_closed_te(ref_x)
        };

        // Expected difference due to coefficient change: 0.6 * 0.0021 * x^4
        let expected_te_diff = 0.00126 * ref_x.powi(4);

        // Actual difference
        let actual_diff = (our_y - ref_y).abs();

        // Should match the expected TE correction within small tolerance
        // The tolerance accounts for floating point and the formula difference
        assert!(
            (actual_diff - expected_te_diff).abs() < 0.0001,
            "At x={:.4}: ref_y={:.6}, our_y={:.6}, diff={:.6}, expected_te_diff={:.6}",
            ref_x,
            ref_y,
            our_y,
            actual_diff,
            expected_te_diff
        );
    }
}

/// Test that our generator produces coordinates matching our closed-TE formula
#[test]
fn test_naca_generator_matches_formula() {
    let generated = naca_4digit("0012", 200).unwrap();

    // Helper: compute NACA 0012 thickness using our closed-TE formula
    fn naca_0012_closed_te(x: f64) -> f64 {
        let t = 0.12;
        (t / 0.2)
            * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
                - 0.1036 * x.powi(4))
    }

    // For each generated point, verify it matches our formula
    for i in 0..generated.x_c.len() {
        let x = generated.x_c[i];
        let y = generated.y_c[i];

        // For symmetric NACA 0012, y should be ±thickness
        let expected_thickness = naca_0012_closed_te(x);
        let expected_y = if y >= 0.0 {
            expected_thickness
        } else {
            -expected_thickness
        };

        assert!(
            (y - expected_y).abs() < 1e-10,
            "Generator mismatch at x={:.6}: got y={:.10}, expected={:.10}",
            x,
            y,
            expected_y
        );
    }
}

/// Test roundtrip: generate -> write -> read -> compare
#[test]
fn test_json_roundtrip() {
    let original = naca_4digit("4412", 120).unwrap();

    // Write to temp file
    let temp_file = NamedTempFile::with_suffix(".json").unwrap();
    write_geometry_to_json(&original, temp_file.path()).unwrap();

    // Read back
    let loaded = read_geometry_from_file(temp_file.path()).unwrap();

    // Compare
    assert_eq!(original.x_c.len(), loaded.x_c.len());
    for i in 0..original.x_c.len() {
        assert_relative_eq!(original.x_c[i], loaded.x_c[i], epsilon = 1e-10);
        assert_relative_eq!(original.y_c[i], loaded.y_c[i], epsilon = 1e-10);
    }
}

/// Test roundtrip: generate -> write DAT -> read DAT -> compare
#[test]
fn test_dat_roundtrip() {
    let original = naca_5digit("23015", 100).unwrap();

    // Write to temp file
    let temp_file = NamedTempFile::with_suffix(".dat").unwrap();
    write_dat_file(&original, "NACA 23015", temp_file.path()).unwrap();

    // Read back
    let (name, loaded) = read_dat_file(temp_file.path()).unwrap();

    // Compare
    assert_eq!(name, "NACA 23015");
    assert_eq!(original.x_c.len(), loaded.x_c.len());
    for i in 0..original.x_c.len() {
        assert_relative_eq!(original.x_c[i], loaded.x_c[i], epsilon = 1e-5);
        assert_relative_eq!(original.y_c[i], loaded.y_c[i], epsilon = 1e-5);
    }
}

/// Test format conversion: JSON -> DAT -> JSON
#[test]
fn test_format_conversion() {
    let original = naca_4digit("2412", 120).unwrap();

    // Write JSON
    let json_file = NamedTempFile::with_suffix(".json").unwrap();
    write_geometry_to_json(&original, json_file.path()).unwrap();

    // Read JSON
    let from_json = read_geometry_from_file(json_file.path()).unwrap();

    // Write DAT
    let dat_file = NamedTempFile::with_suffix(".dat").unwrap();
    write_dat_file(&from_json, "Test", dat_file.path()).unwrap();

    // Read DAT
    let (_, from_dat) = read_dat_file(dat_file.path()).unwrap();

    // Compare original to final
    assert_eq!(original.x_c.len(), from_dat.x_c.len());
    for i in 0..original.x_c.len() {
        assert_relative_eq!(original.x_c[i], from_dat.x_c[i], epsilon = 1e-5);
        assert_relative_eq!(original.y_c[i], from_dat.y_c[i], epsilon = 1e-5);
    }
}

/// Test repaneling preserves airfoil shape
#[test]
fn test_repanel_shape_preservation() {
    let original = naca_4digit("0012", 100).unwrap();
    let repaneled = repanel(&original, 200, 0.15);

    // Maximum thickness should be preserved
    let orig_max_y = original.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let new_max_y = repaneled.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert_relative_eq!(orig_max_y, new_max_y, epsilon = 0.002);

    // Minimum thickness should be preserved
    let orig_min_y = original.y_c.iter().cloned().fold(f64::INFINITY, f64::min);
    let new_min_y = repaneled.y_c.iter().cloned().fold(f64::INFINITY, f64::min);
    assert_relative_eq!(orig_min_y, new_min_y, epsilon = 0.002);
}

/// Test full pipeline: generate -> repanel -> create paneled airfoil
#[test]
fn test_full_geometry_pipeline() {
    // Generate NACA airfoil
    let raw = naca_5digit("23012", 100).unwrap();

    // Repanel with finer distribution
    let repaneled = repanel(&raw, 160, 0.15);

    // Create paneled airfoil with all derived quantities
    let paneled = create_paneled_airfoil(&repaneled);

    // Verify paneled airfoil has correct properties
    assert!(paneled.n > 150);
    assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

    // Arc length should be monotonically increasing
    for i in 1..paneled.s.len() {
        assert!(
            paneled.s[i] > paneled.s[i - 1],
            "Arc length not monotonic at {}",
            i
        );
    }

    // All normal vectors should be unit length
    for i in 0..paneled.n {
        let mag = (paneled.nx[i].powi(2) + paneled.ny[i].powi(2)).sqrt();
        assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
    }

    // Leading edge should be found
    assert!(paneled.sle > 0.0);
    assert!(paneled.sle < paneled.s[paneled.n - 1]);
}

/// Test that symmetric airfoils are indeed symmetric
#[test]
fn test_airfoil_symmetry() {
    let geom = naca_4digit("0015", 160).unwrap();

    let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

    // Should be symmetric about y=0
    assert_relative_eq!(max_y, -min_y, epsilon = 0.001);
}

/// Test cambered airfoils have positive camber
#[test]
fn test_airfoil_camber() {
    // NACA 4-digit cambered
    let geom_4 = naca_4digit("6412", 160).unwrap();
    let max_y_4 = geom_4.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y_4 = geom_4.y_c.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(max_y_4 > -min_y_4, "4-digit should have positive camber");

    // NACA 5-digit cambered
    let geom_5 = naca_5digit("23018", 160).unwrap();
    let max_y_5 = geom_5.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y_5 = geom_5.y_c.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(max_y_5 > -min_y_5, "5-digit should have positive camber");
}

/// Test different thickness values
#[test]
fn test_thickness_values() {
    let thin = naca_4digit("0006", 100).unwrap();
    let medium = naca_4digit("0012", 100).unwrap();
    let thick = naca_4digit("0024", 100).unwrap();

    let thick_thin = thin.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;
    let thick_med = medium.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;
    let thick_thick = thick.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;

    // Thickness should scale appropriately
    assert!(thick_thin < thick_med);
    assert!(thick_med < thick_thick);

    // Check approximate thickness values
    assert!(thick_thin > 0.05 && thick_thin < 0.08);
    assert!(thick_med > 0.10 && thick_med < 0.14);
    assert!(thick_thick > 0.20 && thick_thick < 0.28);
}

// ============================================================================
// Panel Method Validation Tests
// ============================================================================

use yfoil::forces::integrate_forces;
use yfoil::panel::solve_inviscid;

/// Test inviscid panel method lift curve slope against thin airfoil theory.
///
/// Thin airfoil theory predicts: dCL/dα = 2π ≈ 6.28 per radian
/// For a 12% thick airfoil, the actual value is slightly lower due to
/// thickness effects (typically 5.8-6.2).
#[test]
fn test_panel_method_lift_slope() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);

    // Calculate CL at two angles
    let alpha1 = 0.0_f64;
    let alpha2 = 5.0_f64.to_radians();

    let vel1 = solution.velocity_at_alpha(alpha1);
    let vel2 = solution.velocity_at_alpha(alpha2);
    let cp1: Vec<f64> = vel1.iter().map(|&v| 1.0 - v * v).collect();
    let cp2: Vec<f64> = vel2.iter().map(|&v| 1.0 - v * v).collect();

    let coeffs1 = integrate_forces(&airfoil, &cp1, alpha1);
    let coeffs2 = integrate_forces(&airfoil, &cp2, alpha2);

    let dcl_dalpha = (coeffs2.cl - coeffs1.cl) / (alpha2 - alpha1);

    // Should be in range 5.0 to 7.0 (thin airfoil theory ≈ 6.28)
    assert!(
        dcl_dalpha > 5.0 && dcl_dalpha < 7.0,
        "Lift slope {} per radian should be near 2π ≈ 6.28",
        dcl_dalpha
    );
}

/// Test that symmetric airfoil has zero lift at alpha=0
#[test]
fn test_panel_method_symmetric_zero_lift() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);

    let vel = solution.velocity_at_alpha(0.0);
    let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs = integrate_forces(&airfoil, &cp, 0.0);

    assert!(
        coeffs.cl.abs() < 0.01,
        "CL={} should be ~0 for symmetric airfoil at α=0",
        coeffs.cl
    );
}

/// Test that cambered airfoil has positive lift at alpha=0
#[test]
fn test_panel_method_cambered_lift() {
    let geom = naca_4digit("4412", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);

    let vel = solution.velocity_at_alpha(0.0);
    let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs = integrate_forces(&airfoil, &cp, 0.0);

    // NACA 4412 has 4% camber, should produce positive lift at α=0
    // Thin airfoil theory predicts CL ≈ 2π * 2 * (0.04) ≈ 0.5 for 4% camber
    assert!(
        coeffs.cl > 0.2,
        "CL={} should be positive for cambered airfoil at α=0",
        coeffs.cl
    );
}
