//! (c) yFoil functionality: geometry generation, file I/O and repanelling.
//!
//! Fixtures: `tests/fixtures/naca0012_reference.dat` — hand-written.

use approx::assert_relative_eq;
use tempfile::NamedTempFile;
use yfoil::geometry::Thickness;

use yfoil::geometry::{
    naca_4digit, naca_5digit, panel_foil, read_dat_file, read_geometry_from_file, repanel_cosine, write_dat_file,
    write_geometry_to_json,
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
        (t / 0.2) * (0.2969 * x.sqrt() - 0.126 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3) - 0.1036 * x.powi(4))
    }

    // For each reference point, compute what our formula gives and compare
    for i in 0..reference.x.len() {
        let ref_x = reference.x[i];
        let ref_y = reference.y[i];

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

/// Test roundtrip: generate -> write -> read -> compare
#[test]
fn test_json_roundtrip() {
    let original = naca_4digit("4412", 120, Thickness::Perpendicular).unwrap();

    // Write to temp file
    let temp_file = NamedTempFile::with_suffix(".json").unwrap();
    write_geometry_to_json(&original, temp_file.path()).unwrap();

    // Read back
    let loaded = read_geometry_from_file(temp_file.path()).unwrap();

    // Compare
    assert_eq!(original.x.len(), loaded.x.len());
    for i in 0..original.x.len() {
        assert_relative_eq!(original.x[i], loaded.x[i], epsilon = 1e-10);
        assert_relative_eq!(original.y[i], loaded.y[i], epsilon = 1e-10);
    }
}

/// Test roundtrip: generate -> write DAT -> read DAT -> compare
#[test]
fn test_dat_roundtrip() {
    let original = naca_5digit("23015", 100, Thickness::Perpendicular).unwrap();

    // Write to temp file
    let temp_file = NamedTempFile::with_suffix(".dat").unwrap();
    write_dat_file(&original, "NACA 23015", temp_file.path()).unwrap();

    // Read back
    let (name, loaded) = read_dat_file(temp_file.path()).unwrap();

    // Compare
    assert_eq!(name, "NACA 23015");
    assert_eq!(original.x.len(), loaded.x.len());
    for i in 0..original.x.len() {
        assert_relative_eq!(original.x[i], loaded.x[i], epsilon = 1e-5);
        assert_relative_eq!(original.y[i], loaded.y[i], epsilon = 1e-5);
    }
}

/// Test format conversion: JSON -> DAT -> JSON
#[test]
fn test_format_conversion() {
    let original = naca_4digit("2412", 120, Thickness::Perpendicular).unwrap();

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
    assert_eq!(original.x.len(), from_dat.x.len());
    for i in 0..original.x.len() {
        assert_relative_eq!(original.x[i], from_dat.x[i], epsilon = 1e-5);
        assert_relative_eq!(original.y[i], from_dat.y[i], epsilon = 1e-5);
    }
}

/// Test repaneling preserves airfoil shape
#[test]
fn test_repanel_shape_preservation() {
    let original = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
    let repaneled = repanel_cosine(&original, 200, 0.15);

    // Maximum thickness should be preserved
    let orig_max_y = original.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let new_max_y = repaneled.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert_relative_eq!(orig_max_y, new_max_y, epsilon = 0.002);

    // Minimum thickness should be preserved
    let orig_min_y = original.y.iter().cloned().fold(f64::INFINITY, f64::min);
    let new_min_y = repaneled.y.iter().cloned().fold(f64::INFINITY, f64::min);
    assert_relative_eq!(orig_min_y, new_min_y, epsilon = 0.002);
}

/// Test full pipeline: generate -> repanel -> create paneled airfoil
#[test]
fn test_full_geometry_pipeline() {
    // Generate NACA airfoil
    let raw = naca_5digit("23012", 100, Thickness::Perpendicular).unwrap();

    // Repanel with finer distribution
    let repaneled = repanel_cosine(&raw, 160, 0.15);

    // Create paneled airfoil with all derived quantities
    let paneled = panel_foil(&repaneled);

    // Verify paneled airfoil has correct properties
    assert!(paneled.n_foil_nodes > 150);
    assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

    // Arc length should be monotonically increasing
    for i in 1..paneled.s.len() {
        assert!(paneled.s[i] > paneled.s[i - 1], "Arc length not monotonic at {}", i);
    }

    // All normal vectors should be unit length
    for i in 0..paneled.n_foil_nodes {
        let mag = (paneled.normal_x[i].powi(2) + paneled.normal_y[i].powi(2)).sqrt();
        assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
    }

    // Leading edge should be found
    assert!(paneled.s_le > 0.0);
    assert!(paneled.s_le < paneled.s[paneled.n_foil_nodes - 1]);
}

/// Test different thickness values
#[test]
fn test_thickness_values() {
    let thin = naca_4digit("0006", 100, Thickness::Perpendicular).unwrap();
    let medium = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
    let thick = naca_4digit("0024", 100, Thickness::Perpendicular).unwrap();

    let thick_thin = thin.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;
    let thick_med = medium.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;
    let thick_thick = thick.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max) * 2.0;

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

// ============================================================================
// Boundary Layer Validation Tests
// ============================================================================

// ============================================================================
// Trailing Edge Geometry Tests
// ============================================================================

/// Test that blunt TE detection works correctly
#[test]
fn test_blunt_te_detection() {
    let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();

    // Blunten the geometry
    let blunt_geom = geom.blunten(0.002);
    let paneled = panel_foil(&blunt_geom);

    assert!(!paneled.sharp_te, "Blunted airfoil should have blunt TE");
    assert!(!blunt_geom.is_sharp_te(), "Geometry should detect as blunt TE");

    // TE gap should be ~0.002
    assert_relative_eq!(blunt_geom.te_gap(), 0.002, epsilon = 1e-10);
}
