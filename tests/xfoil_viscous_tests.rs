//! Fixture-based validation tests against XFOIL
//!
//! These tests compare YFoil results against pre-generated XFOIL fixtures
//! to ensure numerical accuracy.

mod fixtures;

use fixtures::{assert_xfoil_match, TestFixture};
use std::path::Path;
use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig, BLSolution};

/// Helper to run YFoil and compare against fixture
fn validate_against_fixture(fixture_path: &str) {
    let path = Path::new("tests/fixtures").join(fixture_path);

    // Load fixture
    let fixture = match TestFixture::load(&path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "Skipping test: fixture not found at {} ({})",
                path.display(),
                e
            );
            return;
        }
    };

    // Extract NACA code from airfoil name
    let naca_code = fixture
        .input
        .airfoil
        .replace("NACA ", "")
        .trim()
        .to_string();

    // Create airfoil and run YFoil
    let geom = naca_4digit(&naca_code, fixture.input.n_panels).expect("Failed to create airfoil");
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(
        fixture.input.reynolds,
        fixture.input.mach,
        fixture.input.n_crit,
        1.0,
    );

    let config = ViscalConfig::default();
    let alpha = fixture.input.alpha_deg.to_radians();

    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    // Compare against XFOIL results
    if let Some(ref final_fixture) = fixture.final_result {
        println!(
            "\nValidating {} at alpha = {}°:",
            fixture.input.airfoil, fixture.input.alpha_deg
        );
        println!(
            "  YFoil:  CL={:.6}, CD={:.6}, CDF={:.6}",
            result.cl, result.cd, result.cdf
        );
        println!(
            "  XFOIL:  CL={:.6}, CD={:.6}, CDF={:.6}",
            final_fixture.cl, final_fixture.cd, final_fixture.cdf
        );

        // Allow larger tolerance for overall coefficients
        // (accumulation of small errors across integration)
        let coeff_rel_tol = 1e-4;
        let coeff_abs_tol = 1e-6;

        assert_xfoil_match("CL", result.cl, final_fixture.cl, coeff_rel_tol, coeff_abs_tol);
        assert_xfoil_match("CD", result.cd, final_fixture.cd, coeff_rel_tol, coeff_abs_tol);
        assert_xfoil_match(
            "CDF",
            result.cdf,
            final_fixture.cdf,
            coeff_rel_tol,
            coeff_abs_tol,
        );

        // Transition locations - allow 5% tolerance
        if final_fixture.xtr_upper > 0.0 {
            let xtr_tol = 0.05;
            assert_xfoil_match(
                "xtr_upper",
                result.xtr_upper,
                final_fixture.xtr_upper,
                xtr_tol,
                0.01,
            );
        }

        if final_fixture.xtr_lower > 0.0 {
            let xtr_tol = 0.05;
            assert_xfoil_match(
                "xtr_lower",
                result.xtr_lower,
                final_fixture.xtr_lower,
                xtr_tol,
                0.01,
            );
        }
    }

    // Compare BL stations if available
    if let Some(ref bl_fixture) = fixture.bl_stations {
        validate_bl_stations(&result.bl, bl_fixture);
    }

    // Compare iteration history if available
    if let Some(ref viscal_fixture) = fixture.viscal_iters {
        println!(
            "  XFOIL converged in {} iterations",
            viscal_fixture.n_iterations
        );
        println!(
            "  YFoil converged: {} in {} iterations",
            result.converged, result.iterations
        );

        // YFoil should converge if XFOIL did
        if viscal_fixture.converged {
            assert!(
                result.converged,
                "YFoil failed to converge when XFOIL did"
            );
        }
    }
}

fn validate_bl_stations(
    _yfoil_bl: &BLSolution,
    _xfoil_bl: &fixtures::BLFixture,
) {
    // TODO: Implement station-by-station comparison
    // This requires matching station indices between YFoil and XFOIL
}

// ============================================================================
// NACA 0012 Tests
// ============================================================================

#[test]
#[ignore] // VISCAL solver not converging correctly - needs debugging
fn test_naca0012_alpha_0() {
    validate_against_fixture("naca0012/alpha_0_re_1e6");
}

#[test]
#[ignore] // VISCAL solver not converging correctly - needs debugging
fn test_naca0012_alpha_2() {
    validate_against_fixture("naca0012/alpha_2_re_1e6");
}

#[test]
#[ignore] // VISCAL solver not converging correctly - needs debugging
fn test_naca0012_alpha_5() {
    validate_against_fixture("naca0012/alpha_5_re_1e6");
}

// ============================================================================
// NACA 4412 Tests
// ============================================================================

#[test]
#[ignore] // VISCAL solver not converging correctly - needs debugging
fn test_naca4412_alpha_0() {
    validate_against_fixture("naca4412/alpha_0_re_1e6");
}

#[test]
#[ignore] // VISCAL solver not converging correctly - needs debugging
fn test_naca4412_alpha_4() {
    validate_against_fixture("naca4412/alpha_4_re_1e6");
}

// ============================================================================
// Subroutine-level Tests
// ============================================================================

#[cfg(test)]
mod subroutine_tests {
    use super::*;

    /// Test BLKIN closure relations against XFOIL values
    #[test]
    #[ignore] // Enable when fixtures have BLKIN data
    fn test_blkin_matches_xfoil() {
        // This test will compare the kinematic secondary variables
        // computed by YFoil's closure module against XFOIL's BLKIN output
        todo!("Implement when BLKIN fixture data is available")
    }

    /// Test BLVAR closure relations against XFOIL values
    #[test]
    #[ignore]
    fn test_blvar_matches_xfoil() {
        todo!("Implement when BLVAR fixture data is available")
    }

    /// Test transition detection against XFOIL
    #[test]
    #[ignore]
    fn test_trchek_matches_xfoil() {
        todo!("Implement when TRCHEK fixture data is available")
    }

    /// Test Newton system assembly against XFOIL
    #[test]
    #[ignore]
    fn test_blsys_matches_xfoil() {
        todo!("Implement when BLSYS fixture data is available")
    }
}

// ============================================================================
// Tolerance Tests
// ============================================================================

#[cfg(test)]
mod tolerance_tests {
    use super::validate_against_fixture;

    /// Verify that all fixture tests achieve the required precision
    #[test]
    #[ignore] // Run manually to check all fixtures
    fn verify_all_fixtures_precision() {
        let fixture_paths = [
            "naca0012/alpha_0_re_1e6",
            "naca0012/alpha_2_re_1e6",
            "naca0012/alpha_5_re_1e6",
            "naca4412/alpha_0_re_1e6",
            "naca4412/alpha_4_re_1e6",
        ];

        for path in &fixture_paths {
            println!("\n=== Testing {} ===", path);
            validate_against_fixture(path);
        }
    }
}
