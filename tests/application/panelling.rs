//! (c) yFoil functionality: yFoil's panelling methods (PANGEN-by-curvature and cosine) produce valid geometry for analysis.
//!
//! Fixtures: none.

use yfoil::geometry::Thickness;

use approx::assert_relative_eq;

use yfoil::geometry::{naca_4digit, panel_foil, repanel_by_curvature, repanel_cosine, PangenConfig};

// ============================================================================
// Repanel method tests (XFOIL PANE vs cosine spacing)
// ============================================================================

/// Test that PANE method produces valid geometry for aerodynamic analysis
#[test]
fn test_pane_method_produces_valid_paneled_airfoil() {
    let original = naca_4digit("0012", 100, Thickness::Perpendicular).expect("Failed to create airfoil");
    let config = PangenConfig::default();
    let paned = repanel_by_curvature(&original, 160, &config);
    let paneled = panel_foil(&paned);

    // Verify panel count
    assert_eq!(paneled.n_foil_nodes, 160, "Panel count should match requested");

    // Verify chord is approximately 1.0
    assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

    // Verify arc lengths are monotonically increasing
    for i in 1..paneled.s.len() {
        assert!(
            paneled.s[i] > paneled.s[i - 1],
            "Arc length should be monotonic at index {}",
            i
        );
    }

    // Verify normal vectors have unit length
    for i in 0..paneled.n_foil_nodes {
        let mag = (paneled.normal_x[i].powi(2) + paneled.normal_y[i].powi(2)).sqrt();
        assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
    }

    // Verify LE index is sensible (approximately in the middle)
    assert!(
        paneled.i_le_node > paneled.n_foil_nodes / 4 && paneled.i_le_node < 3 * paneled.n_foil_nodes / 4,
        "LE index {} should be near middle of {}",
        paneled.i_le_node,
        paneled.n_foil_nodes
    );
}

/// Test that cosine method produces valid geometry for aerodynamic analysis
#[test]
fn test_cosine_method_produces_valid_paneled_airfoil() {
    let original = naca_4digit("0012", 100, Thickness::Perpendicular).expect("Failed to create airfoil");
    let cosined = repanel_cosine(&original, 160, 0.15);
    let paneled = panel_foil(&cosined);

    // Verify panel count is approximately correct (may vary slightly)
    assert!(
        paneled.n_foil_nodes > 155 && paneled.n_foil_nodes < 165,
        "Panel count {} should be near 160",
        paneled.n_foil_nodes
    );

    // Verify chord is approximately 1.0
    assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

    // Verify arc lengths are monotonically increasing
    for i in 1..paneled.s.len() {
        assert!(
            paneled.s[i] > paneled.s[i - 1],
            "Arc length should be monotonic at index {}",
            i
        );
    }

    // Verify normal vectors have unit length
    for i in 0..paneled.n_foil_nodes {
        let mag = (paneled.normal_x[i].powi(2) + paneled.normal_y[i].powi(2)).sqrt();
        assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
    }
}

/// Test that PANE and cosine methods produce different panel distributions
#[test]
fn test_pane_and_cosine_methods_differ() {
    let original = naca_4digit("0012", 200, Thickness::Perpendicular).expect("Failed to create airfoil");
    let config = PangenConfig::default();

    let paned = repanel_by_curvature(&original, 160, &config);
    let cosined = repanel_cosine(&original, 160, 0.15);

    // Both should produce geometry with similar extent
    let paned_max_x = paned.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let cosined_max_x = cosined.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert_relative_eq!(paned_max_x, cosined_max_x, epsilon = 0.01);

    // But interior point distributions should differ
    // Count points where x coordinates differ by more than 0.001
    let min_len = paned.x.len().min(cosined.x.len());
    let mut differences = 0;
    for i in 10..min_len.saturating_sub(10) {
        let dx = (paned.x[i] - cosined.x[i]).abs();
        if dx > 0.001 {
            differences += 1;
        }
    }

    // There should be meaningful differences in the distributions
    assert!(
        differences > 10,
        "PANE and cosine should produce different distributions, but only {} points differ",
        differences
    );
}

/// Test that both repaneling methods preserve airfoil shape (extents)
#[test]
fn test_both_methods_preserve_shape() {
    let original = naca_4digit("4412", 120, Thickness::Perpendicular).expect("Failed to create cambered airfoil");

    let config = PangenConfig::default();
    let paned = repanel_by_curvature(&original, 160, &config);
    let cosined = repanel_cosine(&original, 160, 0.15);

    // Original extents
    let orig_max_x = original.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let orig_min_x = original.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let orig_max_y = original.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let orig_min_y = original.y.iter().cloned().fold(f64::INFINITY, f64::min);

    // PANE extents
    let paned_max_x = paned.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let paned_min_x = paned.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let paned_max_y = paned.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let paned_min_y = paned.y.iter().cloned().fold(f64::INFINITY, f64::min);

    // Cosine extents
    let cosined_max_x = cosined.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let cosined_min_x = cosined.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let cosined_max_y = cosined.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let cosined_min_y = cosined.y.iter().cloned().fold(f64::INFINITY, f64::min);

    // Both should preserve x extents
    assert_relative_eq!(paned_max_x, orig_max_x, epsilon = 0.01);
    assert_relative_eq!(paned_min_x, orig_min_x, epsilon = 0.01);
    assert_relative_eq!(cosined_max_x, orig_max_x, epsilon = 0.01);
    assert_relative_eq!(cosined_min_x, orig_min_x, epsilon = 0.01);

    // Both should preserve y extents (within interpolation tolerance)
    assert_relative_eq!(paned_max_y, orig_max_y, epsilon = 0.01);
    assert_relative_eq!(paned_min_y, orig_min_y, epsilon = 0.01);
    assert_relative_eq!(cosined_max_y, orig_max_y, epsilon = 0.01);
    assert_relative_eq!(cosined_min_y, orig_min_y, epsilon = 0.01);
}

/// Test PANE method clusters panels at leading edge (high curvature)
#[test]
fn test_pane_clusters_at_leading_edge() {
    let original = naca_4digit("0012", 200, Thickness::Perpendicular).expect("Failed to create airfoil");
    let config = PangenConfig::default();
    let paned = repanel_by_curvature(&original, 160, &config);

    // Calculate average panel spacing in different regions
    let mut le_spacings = Vec::new();
    let mut mid_spacings = Vec::new();

    for i in 1..paned.x.len() {
        let dx = (paned.x[i] - paned.x[i - 1]).abs();
        let dy = (paned.y[i] - paned.y[i - 1]).abs();
        let ds = (dx * dx + dy * dy).sqrt();

        let avg_x = (paned.x[i] + paned.x[i - 1]) / 2.0;

        if avg_x < 0.1 {
            le_spacings.push(ds);
        } else if avg_x > 0.4 && avg_x < 0.6 {
            mid_spacings.push(ds);
        }
    }

    if !le_spacings.is_empty() && !mid_spacings.is_empty() {
        let avg_le: f64 = le_spacings.iter().sum::<f64>() / le_spacings.len() as f64;
        let avg_mid: f64 = mid_spacings.iter().sum::<f64>() / mid_spacings.len() as f64;

        // LE spacing should be smaller than mid-chord spacing
        assert!(
            avg_le < avg_mid,
            "LE spacing {:.6} should be smaller than mid spacing {:.6}",
            avg_le,
            avg_mid
        );
    }
}

/// Test that PANE config parameters affect the output
#[test]
fn test_pane_config_affects_output() {
    let original = naca_4digit("0012", 200, Thickness::Perpendicular).expect("Failed to create airfoil");

    let config_low_te = PangenConfig {
        te_curvature_ratio: 0.10,
        ..Default::default()
    };
    let config_high_te = PangenConfig {
        te_curvature_ratio: 0.50,
        ..Default::default()
    };

    let paned_low = repanel_by_curvature(&original, 160, &config_low_te);
    let paned_high = repanel_by_curvature(&original, 160, &config_high_te);

    // Count panels in TE region (x > 0.9)
    let te_count_low = paned_low.x.iter().filter(|&&x| x > 0.9).count();
    let te_count_high = paned_high.x.iter().filter(|&&x| x > 0.9).count();

    // Higher CTERAT should produce more panels at TE (more TE bunching)
    println!(
        "TE panel count: CTERAT=0.10 -> {}, CTERAT=0.50 -> {}",
        te_count_low, te_count_high
    );

    // The distributions should be different
    assert!(
        te_count_low != te_count_high || {
            // If counts are same, check that positions differ
            let mut diff_count = 0;
            for i in 0..paned_low.x.len().min(paned_high.x.len()) {
                if (paned_low.x[i] - paned_high.x[i]).abs() > 0.001 {
                    diff_count += 1;
                }
            }
            diff_count > 10
        },
        "Different CTERAT values should produce different distributions"
    );
}
