//! Tests comparing YFoil's source influence (DZDM) computation against XFOIL
//!
//! XFOIL computes DZDM[j] = dPsi/dSig[j] in the PSILIN subroutine.
//! This is the influence of a source at panel node j on the streamfunction
//! at control point i.
//!
//! The reference data was generated using instrumented XFOIL with a 160-panel
//! NACA 0012 airfoil at control point I=1 (first node).

use std::fs;
use yfoil::geometry::{create_paneled_airfoil, Geometry};
use yfoil::panel::solve_inviscid;

/// Load XFOIL DZDM reference values from fixture file
fn load_xfoil_dzdm_fixture() -> Vec<(usize, f64)> {
    let fixture_path = "tests/fixtures/subroutines/psilin/naca0012_160_dzdm_i1.txt";
    let content = fs::read_to_string(fixture_path).expect("Failed to read DZDM fixture file");

    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let parts: Vec<&str> = line.split_whitespace().collect();
            let idx: usize = parts[0].parse().expect("Failed to parse index");
            let val: f64 = parts[1].parse().expect("Failed to parse value");
            (idx, val)
        })
        .collect()
}

/// Test that YFoil's DZDM computation matches XFOIL for NACA 0012 at control point 1
#[test]
fn test_dzdm_matches_xfoil_naca0012_160() {
    // Load the same geometry used for XFOIL fixture
    let json_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/viscous/naca0012_160.json");
    let json_str = fs::read_to_string(json_path).expect("Failed to read JSON file");
    let geom: Geometry = serde_json::from_str(&json_str).expect("Failed to parse JSON");
    let airfoil = create_paneled_airfoil(&geom);

    assert_eq!(airfoil.n, 160, "Airfoil should have 160 panels");

    // Solve inviscid to compute BIJ = -DZDM
    let inviscid = solve_inviscid(&airfoil);
    let bij = inviscid.get_bij().expect("BIJ matrix not computed");

    // Load XFOIL reference values
    let xfoil_dzdm = load_xfoil_dzdm_fixture();

    // Compare DZDM values for control point I=1 (index 0 in 0-indexed)
    // In YFoil: BIJ[i,j] = -DZDM[j] at control point i
    // So DZDM[j] = -BIJ[0,j]
    let tol = 1e-4; // Allow 0.01% relative error for numerical differences

    let mut max_rel_error = 0.0_f64;
    let mut max_error_idx = 0;

    for (xfoil_idx, xfoil_val) in &xfoil_dzdm {
        // XFOIL uses 1-indexed, YFoil uses 0-indexed
        let yfoil_idx = xfoil_idx - 1;
        if yfoil_idx >= airfoil.n {
            continue;
        }

        let yfoil_dzdm = -bij[(0, yfoil_idx)];
        let xfoil_val = *xfoil_val;

        // Compute relative error (use absolute for small values)
        let abs_error = (yfoil_dzdm - xfoil_val).abs();
        let rel_error = if xfoil_val.abs() > 1e-10 {
            abs_error / xfoil_val.abs()
        } else {
            abs_error
        };

        if rel_error > max_rel_error {
            max_rel_error = rel_error;
            max_error_idx = *xfoil_idx;
        }

        // Print mismatch details for debugging
        if rel_error > tol {
            eprintln!(
                "DZDM[{}] mismatch: XFOIL={:+.10e}, YFoil={:+.10e}, rel_err={:.2e}",
                xfoil_idx, xfoil_val, yfoil_dzdm, rel_error
            );
        }
    }

    eprintln!("Max relative error: {:.2e} at index {}", max_rel_error, max_error_idx);

    // Check first few values explicitly
    let (_, xfoil_dzdm_1) = xfoil_dzdm[0];
    let yfoil_dzdm_1 = -bij[(0, 0)];
    eprintln!("DZDM[1]: XFOIL={:+.10e}, YFoil={:+.10e}", xfoil_dzdm_1, yfoil_dzdm_1);

    // Assert that the max error is within tolerance
    assert!(
        max_rel_error < tol,
        "Max relative error {:.2e} at index {} exceeds tolerance {}",
        max_rel_error,
        max_error_idx,
        tol
    );
}

/// Test that DZDM[1] (self-influence of first panel) is computed correctly
#[test]
fn test_dzdm_self_influence() {
    let json_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/viscous/naca0012_160.json");
    let json_str = fs::read_to_string(json_path).expect("Failed to read JSON file");
    let geom: Geometry = serde_json::from_str(&json_str).expect("Failed to parse JSON");
    let airfoil = create_paneled_airfoil(&geom);

    let inviscid = solve_inviscid(&airfoil);
    let bij = inviscid.get_bij().expect("BIJ matrix not computed");

    // Load XFOIL reference from fixture file (first value)
    let xfoil_dzdm = load_xfoil_dzdm_fixture();
    let (_, xfoil_dzdm_1) = xfoil_dzdm[0];
    let yfoil_dzdm_1 = -bij[(0, 0)];

    let rel_error = (yfoil_dzdm_1 - xfoil_dzdm_1).abs() / xfoil_dzdm_1.abs();

    eprintln!(
        "DZDM[1] self-influence: XFOIL={:+.10e}, YFoil={:+.10e}, rel_err={:.2e}",
        xfoil_dzdm_1, yfoil_dzdm_1, rel_error
    );

    // With identical geometry, values should match to machine precision
    assert!(
        rel_error < 1e-6,
        "DZDM[1] relative error {:.2e} exceeds 1e-6",
        rel_error
    );
}
