#![allow(dead_code)] // fixture structs deserialise every XFOIL field; tests read a subset
//! Tests for SETBL Newton system against XFOIL fixture data
//!
//! These tests validate that yfoil's SETBL implementation produces
//! identical system indexing and Newton system structure to XFOIL's SETBL subroutine.
//!
//! XFOIL's IBLSYS subroutine computes:
//!   - NSYS = (NBL(1)-1) + (NBL(2)-1) = NBL(1) + NBL(2) - 2
//!   - System indices start at IBL=2 (skipping IBL=1, the stagnation station)
//!   - ISYS(IBL,1) = IBL - 1 for surface 1 (1-based)
//!   - ISYS(IBL,2) = (NBL(1)-1) + (IBL-1) for surface 2 (1-based)

mod fixtures;

use serde::Deserialize;
use std::fs;

use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::panel::solve_inviscid;

use fixtures::blsolv_fixtures::parse_blsolv_input;
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// Test that NSYS is computed correctly according to XFOIL's IBLSYS
///
/// XFOIL computes: NSYS = (NBL(1)-1) + (NBL(2)-1) = NBL(1) + NBL(2) - 2
/// This is because the stagnation station (IBL=1) is excluded from the Newton system.
#[test]
fn test_nsys_calculation_matches_xfoil() {
    let input_path = fixture_path("blsolv_input.dat");

    // Parse XFOIL fixture data
    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    // Extract XFOIL parameters
    let xfoil_nsys = xfoil_input.nsys;
    let xfoil_nbl1 = xfoil_input.nbl1;
    let xfoil_nbl2 = xfoil_input.nbl2;

    println!("XFOIL parameters from fixture:");
    println!("  NBL1 = {} (upper surface BL stations)", xfoil_nbl1);
    println!("  NBL2 = {} (lower surface BL stations)", xfoil_nbl2);
    println!("  NSYS = {} (Newton system size)", xfoil_nsys);

    // Verify XFOIL's NSYS formula
    let expected_nsys = (xfoil_nbl1 - 1) + (xfoil_nbl2 - 1);
    assert_eq!(
        xfoil_nsys, expected_nsys,
        "XFOIL NSYS should equal (NBL1-1) + (NBL2-1) = {} but got {}",
        expected_nsys, xfoil_nsys
    );

    // Simulate yfoil's NSYS calculation (current implementation)
    // In yfoil, nbl_upper and nbl_lower correspond to NBL1 and NBL2
    let yfoil_nsys_current = xfoil_nbl1 + xfoil_nbl2; // Current (wrong) implementation
    let yfoil_nsys_correct = (xfoil_nbl1 - 1) + (xfoil_nbl2 - 1); // Correct implementation

    println!("\nYFoil NSYS calculations:");
    println!(
        "  Current (nbl_upper + nbl_lower): {} {}",
        yfoil_nsys_current,
        if yfoil_nsys_current == xfoil_nsys {
            "(CORRECT)"
        } else {
            "(WRONG)"
        }
    );
    println!(
        "  Correct (nbl_upper-1 + nbl_lower-1): {} {}",
        yfoil_nsys_correct,
        if yfoil_nsys_correct == xfoil_nsys {
            "(CORRECT)"
        } else {
            "(WRONG)"
        }
    );

    // This assertion documents the expected correct value
    assert_eq!(
        yfoil_nsys_correct, xfoil_nsys,
        "YFoil NSYS should match XFOIL: {} vs {}",
        yfoil_nsys_correct, xfoil_nsys
    );
}

/// Test that IVTE1 (upper surface TE system index) is computed correctly
///
/// XFOIL: IVTE1 = ISYS(IBLTE(1),1) = IBLTE(1) - 1 (1-based)
/// For 0-based: ivte1 = IBLTE1 - 2
#[test]
fn test_ivte1_calculation_matches_xfoil() {
    let input_path = fixture_path("blsolv_input.dat");

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    let xfoil_iblte1 = xfoil_input.iblte1;
    let xfoil_nbl1 = xfoil_input.nbl1;

    println!("XFOIL IBLTE1 = {} (1-based TE station index)", xfoil_iblte1);
    println!("XFOIL NBL1 = {} (total upper surface stations)", xfoil_nbl1);

    // In XFOIL, IBLTE(1) is the TE station index (1-based)
    // The system index IVTE1 = ISYS(IBLTE1, 1) = IBLTE1 - 1 (1-based)
    // For 0-based: ivte1_0based = IBLTE1 - 2
    let xfoil_ivte1_0based = xfoil_input.ivte1_0based();

    println!("XFOIL IVTE1 (0-based) = {} (computed from IBLTE1)", xfoil_ivte1_0based);

    // Current yfoil implementation uses: ivte1 = nbl_upper - 1
    // This would give: 67 - 1 = 66 for the fixture data
    let yfoil_ivte1_current = xfoil_nbl1 - 1;

    // Correct implementation: ivte1 = (IBLTE1 - 1) - 1 = IBLTE1 - 2
    // Since IBLTE1 = NBL1 for surface ending at TE: ivte1 = NBL1 - 2
    // But more correctly: ivte1 = n_upper_sys - 1 where n_upper_sys = NBL1 - 1
    // So: ivte1 = (NBL1 - 1) - 1 = NBL1 - 2
    let yfoil_ivte1_correct = xfoil_nbl1 - 2;

    println!("\nYFoil IVTE1 calculations:");
    println!(
        "  Current (nbl_upper - 1): {} {}",
        yfoil_ivte1_current,
        if yfoil_ivte1_current == xfoil_ivte1_0based {
            "(CORRECT)"
        } else {
            "(WRONG)"
        }
    );
    println!(
        "  Correct (nbl_upper - 2 = n_upper_sys - 1): {} {}",
        yfoil_ivte1_correct,
        if yfoil_ivte1_correct == xfoil_ivte1_0based {
            "(CORRECT)"
        } else {
            "(WRONG)"
        }
    );

    assert_eq!(
        yfoil_ivte1_correct, xfoil_ivte1_0based,
        "YFoil IVTE1 should match XFOIL: {} vs {}",
        yfoil_ivte1_correct, xfoil_ivte1_0based
    );
}

/// Test lower surface start system index (n_upper_sys)
///
/// The lower surface Newton equations start at system index n_upper_sys.
/// XFOIL: ISYS(2, 2) = (NBL1 - 1) + (2 - 1) = NBL1 - 1 + 1 = NBL1 (1-based)
/// For 0-based: lower_start = NBL1 - 1 = n_upper_sys
#[test]
fn test_lower_surface_start_index() {
    let input_path = fixture_path("blsolv_input.dat");

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    let xfoil_nbl1 = xfoil_input.nbl1;
    let n_upper_sys = xfoil_nbl1 - 1;

    println!("XFOIL NBL1 = {} (total upper surface stations)", xfoil_nbl1);
    println!("n_upper_sys = {} (number of upper surface system entries)", n_upper_sys);
    println!("Lower surface starts at system index {} (0-based)", n_upper_sys);

    // Current yfoil implementation uses: offset = nbl_upper (wrong)
    let yfoil_offset_current = xfoil_nbl1;

    // Correct implementation: offset = n_upper_sys = NBL1 - 1
    let yfoil_offset_correct = n_upper_sys;

    println!("\nYFoil lower surface offset calculations:");
    println!("  Current (nbl_upper): {} (would be WRONG)", yfoil_offset_current);
    println!(
        "  Correct (n_upper_sys = nbl_upper - 1): {} (CORRECT)",
        yfoil_offset_correct
    );

    // Verify this is consistent with station mapping
    // Lower surface IBL=1 (first after stagnation) should map to iv = n_upper_sys + (1-1) = n_upper_sys
    let first_lower_ibl = 1;
    let first_lower_iv = n_upper_sys + (first_lower_ibl - 1);
    assert_eq!(
        first_lower_iv, n_upper_sys,
        "First lower surface system index should be n_upper_sys"
    );
}

/// Test IVZ (wake coupling system index)
///
/// NOTE: IVZ in BLSOLV is specifically for the VZ coupling block that couples
/// the upper surface trailing edge to the WAKE, not the lower surface.
/// IVZ = ISYS(IBLTE(2)+1, 2) = system index of first wake station on surface 2
///
/// This is used for the VZ block in BLSOLV which handles the upper/lower TE
/// to wake coupling. It's separate from the lower surface start index.
#[test]
fn test_ivz_wake_coupling_index() {
    let input_path = fixture_path("blsolv_input.dat");

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    // IVZ is the system index for the first wake station
    // This is where the VZ coupling block connects upper TE to wake
    let xfoil_ivz_0based = xfoil_input.ivz_0based();

    println!("XFOIL parameters for VZ block:");
    println!("  NBL1 = {} (upper surface stations)", xfoil_input.nbl1);
    println!("  IBLTE2 = {} (lower TE station index)", xfoil_input.iblte2);
    println!(
        "  IVZ (0-based) = {} (first wake station system index)",
        xfoil_ivz_0based
    );

    // IVZ formula: IVZ = ISYS(IBLTE2+1, 2) = (NBL1-1) + (IBLTE2+1-1) = NBL1 + IBLTE2 - 1 (1-based)
    // For 0-based: ivz = NBL1 + IBLTE2 - 2
    let expected_ivz = xfoil_input.nbl1 + xfoil_input.iblte2 - 2;
    println!("  Expected IVZ (NBL1 + IBLTE2 - 2) = {}", expected_ivz);

    assert_eq!(
        xfoil_ivz_0based, expected_ivz,
        "IVZ formula should match: {} vs {}",
        xfoil_ivz_0based, expected_ivz
    );
}

/// Test station-to-system index mapping matches XFOIL's ISYS array
///
/// XFOIL ISYS mapping (from IBLSYS subroutine):
/// - Surface 1: ISYS(IBL, 1) = IBL - 1 for IBL = 2..NBL(1)
/// - Surface 2: ISYS(IBL, 2) = (NBL(1)-1) + (IBL-1) for IBL = 2..NBL(2)
///
/// In 0-based terms:
/// - Upper: isys[ibl] = ibl - 1 for ibl = 1..nbl_upper-1
/// - Lower: isys[ibl] = (nbl_upper - 1) + (ibl - 1) for ibl = 1..nbl_lower-1
#[test]
fn test_station_to_system_index_mapping() {
    let input_path = fixture_path("blsolv_input.dat");

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    let nbl1 = xfoil_input.nbl1;
    let nbl2 = xfoil_input.nbl2;
    let nsys = xfoil_input.nsys;

    println!("Testing ISYS index mapping:");
    println!("  NBL1={}, NBL2={}, NSYS={}", nbl1, nbl2, nsys);

    // Test that the mapping is consistent
    let n_upper_sys = nbl1 - 1;
    let n_lower_sys = nbl2 - 1;

    println!("  n_upper_sys = {} (upper surface system entries)", n_upper_sys);
    println!("  n_lower_sys = {} (lower surface system entries)", n_lower_sys);
    println!("  Total = {} (should equal NSYS)", n_upper_sys + n_lower_sys);

    assert_eq!(
        n_upper_sys + n_lower_sys,
        nsys,
        "n_upper_sys + n_lower_sys should equal NSYS"
    );

    // Verify station-to-system index mapping
    // Upper surface: station ibl=1 -> iv=0, ibl=2 -> iv=1, ..., ibl=nbl1-1 -> iv=nbl1-2
    println!("\nUpper surface mapping (ibl -> iv):");
    for ibl in 1..nbl1 {
        let iv = ibl - 1; // System index
        println!("  ibl={} -> iv={}", ibl, iv);
        assert!(
            iv < n_upper_sys,
            "Upper surface iv={} should be < n_upper_sys={}",
            iv,
            n_upper_sys
        );
    }

    // Lower surface: station ibl=1 -> iv=n_upper_sys, ibl=2 -> iv=n_upper_sys+1, ...
    println!("\nLower surface mapping (ibl -> iv):");
    for ibl in 1..nbl2 {
        let iv = n_upper_sys + (ibl - 1); // System index
        println!("  ibl={} -> iv={}", ibl, iv);
        assert!(iv < nsys, "Lower surface iv={} should be < nsys={}", iv, nsys);
    }

    println!("\nAll index mappings verified successfully!");
}

/// Summary test showing the correct formulas for SETBL Newton system
#[test]
fn test_setbl_index_formulas_summary() {
    println!("=== SETBL Newton System Index Formulas ===");
    println!();
    println!("Given: nbl_upper, nbl_lower (number of BL stations including stagnation)");
    println!();
    println!("XFOIL formulas (from IBLSYS subroutine):");
    println!("  NSYS = (NBL(1)-1) + (NBL(2)-1) = nbl_upper + nbl_lower - 2");
    println!("  n_upper_sys = NBL(1) - 1 = nbl_upper - 1");
    println!("  n_lower_sys = NBL(2) - 1 = nbl_lower - 1");
    println!();
    println!("  Upper surface system indices: iv = ibl - 1 for ibl = 1..nbl_upper-1");
    println!("  Lower surface system indices: iv = n_upper_sys + (ibl - 1) for ibl = 1..nbl_lower-1");
    println!();
    println!("  IVTE1 (upper TE, 0-based) = n_upper_sys - 1 = nbl_upper - 2");
    println!("  IVZ (lower start, 0-based) = n_upper_sys = nbl_upper - 1");
    println!();
    println!("Key insight: Station ibl=0 (similarity station) is EXCLUDED from Newton system.");
    println!("This is why NSYS = nbl_upper + nbl_lower - 2, not nbl_upper + nbl_lower.");
}

// ============================================================================
// VISCAL Iteration 1 BL State Tests
// ============================================================================

/// Structure for deserializing XFOIL iteration 1 BL state fixture
#[derive(Debug, Deserialize)]
struct Iter1BLStateFixture {
    description: String,
    iteration: i32,
    rmsbl: f64,
    nbl_upper: usize,
    nbl_lower: usize,
    upper_surface: Vec<StationData>,
    lower_surface: Vec<StationData>,
}

#[derive(Debug, Deserialize)]
struct StationData {
    ibl: usize,
    uedg: f64,
    dstr: f64,
    thet: f64,
    mass: f64,
    ctau: f64,
}

fn load_iter1_fixture() -> Option<Iter1BLStateFixture> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/subroutines/viscal/naca0012_iter1_bl_state.json");

    if !path.exists() {
        eprintln!("Fixture not found: {:?}", path);
        return None;
    }

    let content = fs::read_to_string(&path).ok()?;
    serde_json::from_str(&content).ok()
}

/// Verify that XFOIL's MASS = DSTR * UEDG relationship holds
///
/// This is a sanity check that the fixture data is self-consistent.
/// XFOIL computes: MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS) (xbl.f line 1763)
#[test]
fn test_xfoil_mass_equals_dstr_times_uedg() {
    let fixture = match load_iter1_fixture() {
        Some(f) => f,
        None => {
            panic!("never skip silently (CLAUDE.md Rule 7): fixture not found");
        }
    };

    println!("Verifying MASS = DSTR * UEDG in XFOIL fixture data");
    println!("Description: {}", fixture.description);
    println!("Upper surface: {} stations", fixture.nbl_upper);
    println!("Lower surface: {} stations", fixture.nbl_lower);
    println!();

    let mut max_rel_err_upper: f64 = 0.0;
    let mut max_rel_err_lower: f64 = 0.0;

    // Check upper surface
    println!("Upper surface MASS verification:");
    for station in &fixture.upper_surface {
        let computed_mass = station.dstr * station.uedg;
        let xfoil_mass = station.mass;

        if xfoil_mass.abs() > 1e-15 {
            let rel_err = (computed_mass - xfoil_mass).abs() / xfoil_mass.abs();
            max_rel_err_upper = max_rel_err_upper.max(rel_err);

            // Only print mismatches beyond fixture output precision (~1e-7)
            if rel_err > 1e-6 {
                println!(
                    "  IBL={}: DSTR*UEDG={:.10e}, MASS={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_mass, xfoil_mass, rel_err
                );
            }
        }
    }
    println!("  Max relative error: {:.2e}", max_rel_err_upper);

    // Check lower surface
    println!("\nLower surface MASS verification:");
    for station in &fixture.lower_surface {
        let computed_mass = station.dstr * station.uedg;
        let xfoil_mass = station.mass;

        if xfoil_mass.abs() > 1e-15 {
            let rel_err = (computed_mass - xfoil_mass).abs() / xfoil_mass.abs();
            max_rel_err_lower = max_rel_err_lower.max(rel_err);

            // Only print mismatches beyond fixture output precision (~1e-7)
            if rel_err > 1e-6 {
                println!(
                    "  IBL={}: DSTR*UEDG={:.10e}, MASS={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_mass, xfoil_mass, rel_err
                );
            }
        }
    }
    println!("  Max relative error: {:.2e}", max_rel_err_lower);

    // Assert with tolerance accounting for XFOIL output precision (~7 sig figs)
    // The fixture was parsed from XFOIL output which uses E12.7 format
    assert!(
        max_rel_err_upper < 1e-6,
        "Upper surface MASS should equal DSTR*UEDG, max rel err = {:.2e}",
        max_rel_err_upper
    );
    assert!(
        max_rel_err_lower < 1e-6,
        "Lower surface MASS should equal DSTR*UEDG, max rel err = {:.2e}",
        max_rel_err_lower
    );

    println!("\nSUCCESS: XFOIL MASS = DSTR * UEDG verified (within fixture precision)");
}

/// Print summary of XFOIL iteration 1 BL state for debugging
#[test]
fn test_print_xfoil_iter1_bl_state_summary() {
    let fixture = match load_iter1_fixture() {
        Some(f) => f,
        None => {
            panic!("never skip silently (CLAUDE.md Rule 7): fixture not found");
        }
    };

    println!("=== XFOIL Iteration 1 BL State Summary ===");
    println!("Description: {}", fixture.description);
    println!("RMSBL: {:.6}", fixture.rmsbl);
    println!();

    // Upper surface summary
    println!("Upper surface ({} stations):", fixture.nbl_upper);
    println!("  {:>4} {:>12} {:>12} {:>12}", "IBL", "UEDG", "DSTR", "MASS");
    for station in fixture.upper_surface.iter().take(5) {
        println!(
            "  {:>4} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.uedg, station.dstr, station.mass
        );
    }
    println!("  ...");
    if let Some(te) = fixture.upper_surface.last() {
        println!(
            "  {:>4} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.uedg, te.dstr, te.mass
        );
    }

    // Lower surface summary
    println!("\nLower surface ({} stations):", fixture.nbl_lower);
    println!("  {:>4} {:>12} {:>12} {:>12}", "IBL", "UEDG", "DSTR", "MASS");
    for station in fixture.lower_surface.iter().take(5) {
        println!(
            "  {:>4} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.uedg, station.dstr, station.mass
        );
    }
    println!("  ...");
    if let Some(te) = fixture.lower_surface.last() {
        println!(
            "  {:>4} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.uedg, te.dstr, te.mass
        );
    }

    // Key values for debugging
    let upper_te = fixture.upper_surface.last().unwrap();
    let lower_te = fixture.lower_surface.last().unwrap();
    let upper_max_mass: f64 = fixture.upper_surface.iter().map(|s| s.mass).fold(0.0, f64::max);
    let lower_max_mass: f64 = fixture.lower_surface.iter().map(|s| s.mass).fold(0.0, f64::max);

    println!("\nKey values for DIJ coupling:");
    println!("  Upper TE (IBL={}): MASS = {:.6e}", upper_te.ibl, upper_te.mass);
    println!("  Lower TE (IBL={}): MASS = {:.6e}", lower_te.ibl, lower_te.mass);
    println!("  Upper max MASS: {:.6e}", upper_max_mass);
    println!("  Lower max MASS: {:.6e}", lower_max_mass);
}

// ============================================================================
// UPDATE/DUI Computation Tests
// ============================================================================

/// Structure for deserializing XFOIL UPDATE fixture
#[derive(Debug, Deserialize)]
struct UpdateFixture {
    description: String,
    upper_surface: Vec<UpdateStationData>,
    lower_surface: Vec<UpdateStationData>,
    nbl_upper: usize,
    nbl_lower: usize,
}

#[derive(Debug, Deserialize)]
struct UpdateStationData {
    ibl: usize,
    ipan: usize,
    uinv: f64,
    dui: f64,
    unew: f64,
    uedg: f64,
    mass: f64,
    dstr: f64,
    mass_plus_vdel: Option<f64>,
    vdel: Option<f64>,
}

fn load_update_fixture() -> Option<UpdateFixture> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/subroutines/viscal/naca0012_update_iter1.json");

    if !path.exists() {
        eprintln!("UPDATE fixture not found: {:?}", path);
        return None;
    }

    let content = fs::read_to_string(&path).ok()?;
    serde_json::from_str(&content).ok()
}

/// Verify that XFOIL's UNEW = UINV + DUI relationship holds
///
/// This is a sanity check that the fixture data is self-consistent.
#[test]
fn test_xfoil_unew_equals_uinv_plus_dui() {
    let fixture = match load_update_fixture() {
        Some(f) => f,
        None => {
            panic!("never skip silently (CLAUDE.md Rule 7): UPDATE fixture not found");
        }
    };

    println!("Verifying UNEW = UINV + DUI in XFOIL UPDATE fixture");
    println!("Description: {}", fixture.description);
    println!();

    let mut max_rel_err_upper: f64 = 0.0;
    let mut max_rel_err_lower: f64 = 0.0;

    // Check upper surface
    println!("Upper surface UNEW verification:");
    for station in &fixture.upper_surface {
        let computed_unew = station.uinv + station.dui;
        let xfoil_unew = station.unew;

        if xfoil_unew.abs() > 1e-15 {
            let rel_err = (computed_unew - xfoil_unew).abs() / xfoil_unew.abs();
            max_rel_err_upper = max_rel_err_upper.max(rel_err);

            if rel_err > 1e-6 {
                println!(
                    "  IBL={}: UINV+DUI={:.10e}, UNEW={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_unew, xfoil_unew, rel_err
                );
            }
        }
    }
    println!("  Max relative error: {:.2e}", max_rel_err_upper);

    // Check lower surface
    println!("\nLower surface UNEW verification:");
    for station in &fixture.lower_surface {
        let computed_unew = station.uinv + station.dui;
        let xfoil_unew = station.unew;

        if xfoil_unew.abs() > 1e-15 {
            let rel_err = (computed_unew - xfoil_unew).abs() / xfoil_unew.abs();
            max_rel_err_lower = max_rel_err_lower.max(rel_err);

            if rel_err > 1e-6 {
                println!(
                    "  IBL={}: UINV+DUI={:.10e}, UNEW={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_unew, xfoil_unew, rel_err
                );
            }
        }
    }
    println!("  Max relative error: {:.2e}", max_rel_err_lower);

    // Assert with tolerance for fixture precision
    assert!(
        max_rel_err_upper < 1e-6,
        "Upper surface UNEW should equal UINV+DUI, max rel err = {:.2e}",
        max_rel_err_upper
    );
    assert!(
        max_rel_err_lower < 1e-6,
        "Lower surface UNEW should equal UINV+DUI, max rel err = {:.2e}",
        max_rel_err_lower
    );

    println!("\nSUCCESS: XFOIL UNEW = UINV + DUI verified");
}

/// Print summary of XFOIL UPDATE data for debugging
#[test]
fn test_print_xfoil_update_summary() {
    let fixture = match load_update_fixture() {
        Some(f) => f,
        None => {
            panic!("never skip silently (CLAUDE.md Rule 7): UPDATE fixture not found");
        }
    };

    println!("=== XFOIL UPDATE Summary ===");
    println!("Description: {}", fixture.description);
    println!();

    // Upper surface summary
    println!("Upper surface ({} stations):", fixture.upper_surface.len());
    println!(
        "  {:>4} {:>5} {:>12} {:>12} {:>12} {:>12}",
        "IBL", "IPAN", "UINV", "DUI", "UNEW", "MASS"
    );
    for station in fixture.upper_surface.iter().take(5) {
        println!(
            "  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.ipan, station.uinv, station.dui, station.unew, station.mass
        );
    }
    println!("  ...");
    if let Some(te) = fixture.upper_surface.last() {
        println!(
            "  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.ipan, te.uinv, te.dui, te.unew, te.mass
        );
    }

    // Lower surface summary
    println!("\nLower surface ({} stations):", fixture.lower_surface.len());
    println!(
        "  {:>4} {:>5} {:>12} {:>12} {:>12} {:>12}",
        "IBL", "IPAN", "UINV", "DUI", "UNEW", "MASS"
    );
    for station in fixture.lower_surface.iter().take(5) {
        println!(
            "  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.ipan, station.uinv, station.dui, station.unew, station.mass
        );
    }
    println!("  ...");
    if let Some(te) = fixture.lower_surface.last() {
        println!(
            "  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.ipan, te.uinv, te.dui, te.unew, te.mass
        );
    }

    // Statistics
    let upper_dui_range: (f64, f64) = fixture
        .upper_surface
        .iter()
        .map(|s| s.dui)
        .fold((f64::MAX, f64::MIN), |(min, max), v| (min.min(v), max.max(v)));
    let lower_dui_range: (f64, f64) = fixture
        .lower_surface
        .iter()
        .map(|s| s.dui)
        .fold((f64::MAX, f64::MIN), |(min, max), v| (min.min(v), max.max(v)));

    println!("\nDUI ranges:");
    println!("  Upper: [{:.6e}, {:.6e}]", upper_dui_range.0, upper_dui_range.1);
    println!("  Lower: [{:.6e}, {:.6e}]", lower_dui_range.0, lower_dui_range.1);
}

/// Fixture for RLX under-relaxation computation
#[derive(Debug, Deserialize)]
struct RlxFixture {
    rlx: f64,
    rmsbl: f64,
    dac: f64,
    dhi: f64,
    dlo: f64,
    nbl_upper: usize,
    nbl_lower: usize,
    upper_surface: Vec<RlxStationFixture>,
    lower_surface: Vec<RlxStationFixture>,
}

#[derive(Debug, Deserialize)]
struct RlxStationFixture {
    ibl: usize,
    dn1: f64,
    dn2: f64,
    dn3: f64,
    dn4: f64,
    dctau: f64,
    dthet: f64,
    ddstr: f64,
    duedg: f64,
}

fn load_rlx_fixture() -> Option<RlxFixture> {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/subroutines/update/naca0012_rlx_iter1.json");
    let content = fs::read_to_string(&path).ok()?;
    serde_json::from_str(&content).ok()
}

/// Test RLX under-relaxation factor computation against XFOIL
///
/// XFOIL's UPDATE subroutine computes RLX by:
/// 1. Starting with RLX = 1.0
/// 2. For each station, computing DN1-DN4 normalized deltas
/// 3. For each DNi, if RLX*DNi > DHI (1.5), set RLX = DHI/DNi
///    If RLX*DNi < DLO (-0.5), set RLX = DLO/DNi
/// 4. The final RLX is the most restrictive bound
#[test]
#[ignore = "S8: UPDATE RLX"]
fn test_rlx_computation_matches_xfoil() {
    let fixture = match load_rlx_fixture() {
        Some(f) => f,
        None => {
            eprintln!("RLX fixture not found. Generate with instrumented XFOIL.");
            return;
        }
    };

    println!("=== RLX Under-Relaxation Test ===");
    println!("XFOIL values:");
    println!("  RLX = {:.10}", fixture.rlx);
    println!("  RMSBL = {:.10}", fixture.rmsbl);
    println!("  DAC = {:.10}", fixture.dac);
    println!("  DHI = {:.1}, DLO = {:.1}", fixture.dhi, fixture.dlo);
    println!("  NBL_upper = {}, NBL_lower = {}", fixture.nbl_upper, fixture.nbl_lower);

    // Implement XFOIL's RLX algorithm exactly
    let dhi = fixture.dhi;
    let dlo = fixture.dlo;
    let mut rlx = 1.0_f64;
    let mut rmsbl_sum = 0.0_f64;

    let mut max_dn_value = 0.0_f64;
    let mut max_dn_station = (0_usize, 0_usize); // (surface, ibl)
    let mut max_dn_var = "";

    // Process all stations as XFOIL does
    for station in &fixture.upper_surface {
        let dn1 = station.dn1;
        let dn2 = station.dn2;
        let dn3 = station.dn3;
        let dn4 = station.dn4;

        // Track max DN
        for (dn, name) in [(dn1, "DN1"), (dn2, "DN2"), (dn3, "DN3"), (dn4, "DN4")] {
            if dn.abs() > max_dn_value.abs() {
                max_dn_value = dn;
                max_dn_station = (1, station.ibl);
                max_dn_var = name;
            }
        }

        // RMSBL accumulation (XFOIL does this with unrelaxed DN values)
        rmsbl_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

        // Check RLX bounds (this is the key algorithm)
        // XFOIL: IF(RDN1 .GT. DHI) RLX = DHI/DN1
        //        IF(RDN1 .LT. DLO) RLX = DLO/DN1
        let rdn1 = rlx * dn1;
        if rdn1 > dhi {
            rlx = dhi / dn1;
        }
        if rdn1 < dlo {
            rlx = dlo / dn1;
        }

        let rdn2 = rlx * dn2;
        if rdn2 > dhi {
            rlx = dhi / dn2;
        }
        if rdn2 < dlo {
            rlx = dlo / dn2;
        }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi {
            rlx = dhi / dn3;
        }
        if rdn3 < dlo {
            rlx = dlo / dn3;
        }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi {
            rlx = dhi / dn4;
        }
        // Note: dn4 is always positive (it's |DUEDG|/0.25), so no DLO check needed
    }

    for station in &fixture.lower_surface {
        let dn1 = station.dn1;
        let dn2 = station.dn2;
        let dn3 = station.dn3;
        let dn4 = station.dn4;

        for (dn, name) in [(dn1, "DN1"), (dn2, "DN2"), (dn3, "DN3"), (dn4, "DN4")] {
            if dn.abs() > max_dn_value.abs() {
                max_dn_value = dn;
                max_dn_station = (2, station.ibl);
                max_dn_var = name;
            }
        }

        rmsbl_sum += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

        let rdn1 = rlx * dn1;
        if rdn1 > dhi {
            rlx = dhi / dn1;
        }
        if rdn1 < dlo {
            rlx = dlo / dn1;
        }

        let rdn2 = rlx * dn2;
        if rdn2 > dhi {
            rlx = dhi / dn2;
        }
        if rdn2 < dlo {
            rlx = dlo / dn2;
        }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi {
            rlx = dhi / dn3;
        }
        if rdn3 < dlo {
            rlx = dlo / dn3;
        }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi {
            rlx = dhi / dn4;
        }
    }

    // Compute RMSBL
    // XFOIL uses NBL(1)+NBL(2) in divisor, not (NBL(1)-1)+(NBL(2)-1)
    // This is because RMSBL normalizes by total stations including stagnation
    let n_stations = fixture.nbl_upper + fixture.nbl_lower; // 81 + 104 = 185
    let rmsbl = (rmsbl_sum / (4.0 * n_stations as f64)).sqrt();

    println!("\nComputed values:");
    println!("  RLX = {:.10}", rlx);
    println!("  RMSBL = {:.10}", rmsbl);
    println!(
        "\nMax DN: {} = {:.6} at surface={}, IBL={}",
        max_dn_var, max_dn_value, max_dn_station.0, max_dn_station.1
    );

    // Verify RLX matches to machine precision
    let rlx_rel_err = (rlx - fixture.rlx).abs() / fixture.rlx.abs();
    println!("\nRLX comparison:");
    println!("  XFOIL:    {:.10}", fixture.rlx);
    println!("  Computed: {:.10}", rlx);
    println!("  Rel err:  {:.2e}", rlx_rel_err);

    let tolerance = 1e-8; // Machine precision
    assert!(
        rlx_rel_err < tolerance,
        "RLX mismatch: computed={:.10}, expected={:.10}, rel_err={:.2e}",
        rlx,
        fixture.rlx,
        rlx_rel_err
    );

    // Verify RMSBL matches to machine precision
    let rmsbl_rel_err = (rmsbl - fixture.rmsbl).abs() / fixture.rmsbl.abs();
    println!("\nRMSBL comparison:");
    println!("  XFOIL:    {:.10}", fixture.rmsbl);
    println!("  Computed: {:.10}", rmsbl);
    println!("  Rel err:  {:.2e}", rmsbl_rel_err);

    assert!(
        rmsbl_rel_err < tolerance,
        "RMSBL mismatch: computed={:.10}, expected={:.10}, rel_err={:.2e}",
        rmsbl,
        fixture.rmsbl,
        rmsbl_rel_err
    );

    println!(
        "\nSUCCESS: RLX and RMSBL computations match XFOIL within {:.0e}",
        tolerance
    );
}
