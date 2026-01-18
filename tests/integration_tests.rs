//! Integration tests for yfoil geometry module

use approx::assert_relative_eq;
use tempfile::NamedTempFile;

use yfoil::geometry::{
    create_paneled_airfoil, naca_4digit, naca_5digit, read_dat_file, read_geometry_from_file,
    repanel_cosine, write_dat_file, write_geometry_to_json,
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
    let repaneled = repanel_cosine(&original, 200, 0.15);

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
    let repaneled = repanel_cosine(&raw, 160, 0.15);

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

// ============================================================================
// Boundary Layer Validation Tests
// ============================================================================

use yfoil::bl::{cf_lam, dampl, hkin, hs_lam, FlowConditions, FlowRegime};

/// Test Blasius flat plate solution for laminar boundary layer
///
/// The Blasius solution provides the exact laminar BL solution for a flat plate:
/// - θ = 0.664 * sqrt(ν*x/U)   (momentum thickness)
/// - δ* = 1.7208 * sqrt(ν*x/U) (displacement thickness)
/// - H = δ*/θ = 2.591          (shape factor)
/// - Cf = 0.664 / sqrt(Re_x)   (skin friction)
///
/// Reference: Schlichting, "Boundary Layer Theory"
#[test]
fn test_blasius_shape_factor() {
    // Blasius shape factor H = 2.591
    let h_blasius = 2.591;
    let (hk, _, _) = hkin(h_blasius, 0.0);

    // At incompressible conditions, Hk = H
    assert_relative_eq!(hk, h_blasius, epsilon = 1e-10);
}

#[test]
fn test_blasius_skin_friction() {
    // Test Cf correlation at Blasius conditions
    // Blasius: Cf = 0.664 / sqrt(Re_x)
    // For Re_θ, we have Re_x = (Re_θ / 0.664)² approximately
    // So at Re_θ = 1000, Re_x ≈ 2.27e6, and Cf ≈ 0.664/1507 ≈ 0.00044
    // But we test Cf * sqrt(Re_θ) which should be roughly constant

    let hk_blasius = 2.591;

    // Test at different Re_θ values
    for rt in [500.0, 1000.0, 2000.0, 5000.0] {
        let result = cf_lam(hk_blasius, rt, 0.0);
        // Cf should be positive for attached laminar flow
        assert!(result.val > 0.0, "Cf should be positive at Re_θ={}", rt);

        // Cf * sqrt(Re_θ) should be approximately constant for Blasius
        // The Falkner-Skan correlation gives Cf * Re_θ ≈ constant for given Hk
        let cf_rt = result.val * rt;
        assert!(
            cf_rt > 0.1 && cf_rt < 1.0,
            "Cf*Re_θ = {} should be O(1) at Re_θ={}",
            cf_rt,
            rt
        );
    }
}

#[test]
fn test_blasius_energy_shape_factor() {
    // Test H* correlation at Blasius conditions
    // For Blasius flow, H* ≈ 1.573
    let hk_blasius = 2.591;
    let result = hs_lam(hk_blasius, 1000.0, 0.0);

    // H* should be between 1.5 and 1.7 for Blasius-like conditions
    assert!(
        result.val > 1.5 && result.val < 1.7,
        "H* = {} should be ~1.57 for Blasius",
        result.val
    );
}

#[test]
fn test_laminar_amplification_off_below_critical() {
    // Below critical Re_θ, amplification should be zero
    let hk_blasius = 2.591;
    let theta = 0.001; // Small momentum thickness

    // At low Re_θ (below ~200 for Blasius)
    let (ax, _, _, _) = dampl(hk_blasius, theta, 100.0);
    assert_eq!(
        ax, 0.0,
        "Amplification should be 0 below critical Re_θ"
    );
}

#[test]
fn test_laminar_amplification_on_above_critical() {
    // Above critical Re_θ, amplification should be positive
    let hk_blasius = 2.591;
    let theta = 0.001;

    // At high Re_θ (well above critical)
    let (ax, _, _, _) = dampl(hk_blasius, theta, 5000.0);
    assert!(
        ax > 0.0,
        "Amplification should be positive above critical Re_θ"
    );
}

#[test]
fn test_flow_conditions_construction() {
    let cond = FlowConditions::new(1_000_000.0, 0.3, 9.0, 1.0);

    assert_eq!(cond.reynolds, 1_000_000.0);
    assert_eq!(cond.mach, 0.3);
    assert_relative_eq!(cond.msq, 0.09, epsilon = 1e-10);
    assert_eq!(cond.ncrit, 9.0);
    // ν = chord / Re = 1.0 / 1e6 = 1e-6
    assert_relative_eq!(cond.nu, 1e-6, epsilon = 1e-12);
}

#[test]
fn test_flow_regime_enum() {
    let laminar = FlowRegime::Laminar;
    let turbulent = FlowRegime::Turbulent;
    let wake = FlowRegime::Wake;

    assert_eq!(laminar, FlowRegime::Laminar);
    assert_ne!(laminar, turbulent);
    assert_ne!(turbulent, wake);
}

// ============================================================================
// Trailing Edge Geometry Tests
// ============================================================================

/// Test that sharp TE detection works correctly
#[test]
fn test_sharp_te_detection() {
    // NACA generator creates closed TE (sharp)
    let geom = naca_4digit("0012", 160).unwrap();
    let paneled = create_paneled_airfoil(&geom);

    assert!(
        paneled.sharp_te,
        "NACA 0012 from generator should have sharp TE"
    );
    assert!(
        geom.is_sharp_te(),
        "Geometry should also detect as sharp TE"
    );
}

/// Test that blunt TE detection works correctly
#[test]
fn test_blunt_te_detection() {
    let geom = naca_4digit("0012", 160).unwrap();

    // Blunten the geometry
    let blunt_geom = geom.blunten(0.002);
    let paneled = create_paneled_airfoil(&blunt_geom);

    assert!(
        !paneled.sharp_te,
        "Blunted airfoil should have blunt TE"
    );
    assert!(
        !blunt_geom.is_sharp_te(),
        "Geometry should detect as blunt TE"
    );

    // TE gap should be ~0.002
    assert_relative_eq!(blunt_geom.te_gap(), 0.002, epsilon = 1e-10);
}

/// Test panel method with sharp TE gives near-zero CDp for symmetric airfoil at α=0
///
/// This validates that the curvature extrapolation condition (XFOIL equation 9)
/// properly eliminates the TE singularity.
#[test]
fn test_sharp_te_zero_cdp_symmetric() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    assert!(airfoil.sharp_te, "Should use sharp TE handling");

    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs = integrate_forces(&airfoil, &cp, 0.0);

    // For symmetric airfoil at α=0, CDp should be essentially zero
    assert!(
        coeffs.cdp.abs() < 0.001,
        "CDp = {} should be ~0 for symmetric airfoil at α=0 with sharp TE",
        coeffs.cdp
    );

    // CL should also be essentially zero
    assert!(
        coeffs.cl.abs() < 0.01,
        "CL = {} should be ~0 for symmetric airfoil at α=0",
        coeffs.cl
    );
}

/// Test panel method with blunt TE still produces reasonable results
///
/// Blunt TE uses standard flow tangency equations (no curvature extrapolation).
/// The Kutta condition still ensures finite TE velocity, and results should
/// be physically reasonable even if not as clean as sharp TE.
#[test]
fn test_blunt_te_reasonable_results() {
    let geom = naca_4digit("0012", 160).unwrap();
    let blunt_geom = geom.blunten(0.002); // 0.2% gap
    let airfoil = create_paneled_airfoil(&blunt_geom);

    assert!(!airfoil.sharp_te, "Should use blunt TE handling");

    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs = integrate_forces(&airfoil, &cp, 0.0);

    // For symmetric airfoil at α=0, CL should still be near zero
    assert!(
        coeffs.cl.abs() < 0.05,
        "CL = {} should be near 0 for symmetric airfoil at α=0",
        coeffs.cl
    );

    // CDp may not be exactly zero but should be small
    // The blunt TE may cause some small asymmetry in the solution
    assert!(
        coeffs.cdp.abs() < 0.01,
        "CDp = {} should be small for blunt TE",
        coeffs.cdp
    );

    // Velocities should be smooth (no huge spikes)
    let max_vel = vel.iter().map(|v| v.abs()).fold(f64::NEG_INFINITY, f64::max);
    assert!(
        max_vel < 3.0,
        "Max velocity {} should be reasonable (no singularity)",
        max_vel
    );
}

/// Test that lift slope is correct for both sharp and blunt TE
#[test]
fn test_te_type_lift_slope_comparison() {
    // Sharp TE
    let geom_sharp = naca_4digit("0012", 160).unwrap();
    let airfoil_sharp = create_paneled_airfoil(&geom_sharp);
    let solution_sharp = solve_inviscid(&airfoil_sharp);

    // Blunt TE
    let geom_blunt = geom_sharp.blunten(0.002);
    let airfoil_blunt = create_paneled_airfoil(&geom_blunt);
    let solution_blunt = solve_inviscid(&airfoil_blunt);

    // Calculate lift slope for both
    let alpha1 = 0.0_f64;
    let alpha2 = 5.0_f64.to_radians();

    // Sharp TE lift slope
    let vel1_sharp = solution_sharp.velocity_at_alpha(alpha1);
    let vel2_sharp = solution_sharp.velocity_at_alpha(alpha2);
    let cp1_sharp: Vec<f64> = vel1_sharp.iter().map(|&v| 1.0 - v * v).collect();
    let cp2_sharp: Vec<f64> = vel2_sharp.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs1_sharp = integrate_forces(&airfoil_sharp, &cp1_sharp, alpha1);
    let coeffs2_sharp = integrate_forces(&airfoil_sharp, &cp2_sharp, alpha2);
    let dcl_dalpha_sharp = (coeffs2_sharp.cl - coeffs1_sharp.cl) / (alpha2 - alpha1);

    // Blunt TE lift slope
    let vel1_blunt = solution_blunt.velocity_at_alpha(alpha1);
    let vel2_blunt = solution_blunt.velocity_at_alpha(alpha2);
    let cp1_blunt: Vec<f64> = vel1_blunt.iter().map(|&v| 1.0 - v * v).collect();
    let cp2_blunt: Vec<f64> = vel2_blunt.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs1_blunt = integrate_forces(&airfoil_blunt, &cp1_blunt, alpha1);
    let coeffs2_blunt = integrate_forces(&airfoil_blunt, &cp2_blunt, alpha2);
    let dcl_dalpha_blunt = (coeffs2_blunt.cl - coeffs1_blunt.cl) / (alpha2 - alpha1);

    // Both should be close to thin airfoil theory (2π ≈ 6.28)
    assert!(
        dcl_dalpha_sharp > 5.0 && dcl_dalpha_sharp < 7.5,
        "Sharp TE lift slope {} should be near 2π",
        dcl_dalpha_sharp
    );
    assert!(
        dcl_dalpha_blunt > 5.0 && dcl_dalpha_blunt < 7.5,
        "Blunt TE lift slope {} should be near 2π",
        dcl_dalpha_blunt
    );

    // They should be similar (within 10%)
    let diff_pct = ((dcl_dalpha_sharp - dcl_dalpha_blunt) / dcl_dalpha_sharp).abs() * 100.0;
    assert!(
        diff_pct < 15.0,
        "Lift slopes should be similar: sharp={}, blunt={}, diff={}%",
        dcl_dalpha_sharp,
        dcl_dalpha_blunt,
        diff_pct
    );
}

/// Test Kutta condition is satisfied for both sharp and blunt TE
#[test]
fn test_kutta_condition_both_te_types() {
    // Sharp TE
    let geom_sharp = naca_4digit("0012", 160).unwrap();
    let airfoil_sharp = create_paneled_airfoil(&geom_sharp);
    let solution_sharp = solve_inviscid(&airfoil_sharp);

    let kutta_sharp_0 = solution_sharp.gam_0[0] + solution_sharp.gam_0[airfoil_sharp.n - 1];
    let kutta_sharp_90 = solution_sharp.gam_90[0] + solution_sharp.gam_90[airfoil_sharp.n - 1];
    assert!(
        kutta_sharp_0.abs() < 1e-9,
        "Sharp TE Kutta condition violated for α=0°: {}",
        kutta_sharp_0
    );
    assert!(
        kutta_sharp_90.abs() < 1e-9,
        "Sharp TE Kutta condition violated for α=90°: {}",
        kutta_sharp_90
    );

    // Blunt TE
    let geom_blunt = geom_sharp.blunten(0.002);
    let airfoil_blunt = create_paneled_airfoil(&geom_blunt);
    let solution_blunt = solve_inviscid(&airfoil_blunt);

    let kutta_blunt_0 = solution_blunt.gam_0[0] + solution_blunt.gam_0[airfoil_blunt.n - 1];
    let kutta_blunt_90 = solution_blunt.gam_90[0] + solution_blunt.gam_90[airfoil_blunt.n - 1];
    assert!(
        kutta_blunt_0.abs() < 1e-9,
        "Blunt TE Kutta condition violated for α=0°: {}",
        kutta_blunt_0
    );
    assert!(
        kutta_blunt_90.abs() < 1e-9,
        "Blunt TE Kutta condition violated for α=90°: {}",
        kutta_blunt_90
    );
}

/// Test that the gamma distribution is smooth near TE for sharp case
///
/// With curvature extrapolation, gamma should smoothly approach zero at TE
/// from both upper and lower surfaces.
#[test]
fn test_sharp_te_smooth_gamma() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);

    // For α=0° symmetric case, check gamma at TE nodes
    let n = airfoil.n;
    let gam = solution.gamma_at_alpha(0.0);

    // TE gammas should be small (Kutta condition + smooth approach)
    assert!(
        gam[0].abs() < 2.0,
        "Upper TE gamma {} should be small with curvature extrapolation",
        gam[0]
    );
    assert!(
        gam[n - 1].abs() < 2.0,
        "Lower TE gamma {} should be small with curvature extrapolation",
        gam[n - 1]
    );

    // Second derivatives at TE should be similar (curvature condition)
    let curv_upper = gam[2] - 2.0 * gam[1] + gam[0];
    let curv_lower = gam[n - 3] - 2.0 * gam[n - 2] + gam[n - 1];
    assert!(
        (curv_upper - curv_lower).abs() < 0.1,
        "TE curvatures should match: upper={}, lower={}",
        curv_upper,
        curv_lower
    );
}
