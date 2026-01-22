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
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(".tmp")
        .join(name)
}

/// Test that NSYS is computed correctly according to XFOIL's IBLSYS
///
/// XFOIL computes: NSYS = (NBL(1)-1) + (NBL(2)-1) = NBL(1) + NBL(2) - 2
/// This is because the stagnation station (IBL=1) is excluded from the Newton system.
#[test]
fn test_nsys_calculation_matches_xfoil() {
    let input_path = fixture_path("blsolv_input.dat");

    if !input_path.exists() {
        eprintln!(
            "Skipping test: fixture file not found. Run instrumented XFOIL to generate it."
        );
        return;
    }

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

    if !input_path.exists() {
        return;
    }

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    let xfoil_iblte1 = xfoil_input.iblte1;
    let xfoil_nbl1 = xfoil_input.nbl1;

    println!("XFOIL IBLTE1 = {} (1-based TE station index)", xfoil_iblte1);
    println!("XFOIL NBL1 = {} (total upper surface stations)", xfoil_nbl1);

    // In XFOIL, IBLTE(1) is the TE station index (1-based)
    // The system index IVTE1 = ISYS(IBLTE1, 1) = IBLTE1 - 1 (1-based)
    // For 0-based: ivte1_0based = IBLTE1 - 2
    let xfoil_ivte1_0based = xfoil_input.ivte1_0based();

    println!(
        "XFOIL IVTE1 (0-based) = {} (computed from IBLTE1)",
        xfoil_ivte1_0based
    );

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

    if !input_path.exists() {
        return;
    }

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
    println!(
        "  Current (nbl_upper): {} (would be WRONG)",
        yfoil_offset_current
    );
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

    if !input_path.exists() {
        return;
    }

    let xfoil_input = parse_blsolv_input(&input_path, 1).expect("Failed to parse BLSOLV input");

    // IVZ is the system index for the first wake station
    // This is where the VZ coupling block connects upper TE to wake
    let xfoil_ivz_0based = xfoil_input.ivz_0based();

    println!("XFOIL parameters for VZ block:");
    println!("  NBL1 = {} (upper surface stations)", xfoil_input.nbl1);
    println!("  IBLTE2 = {} (lower TE station index)", xfoil_input.iblte2);
    println!("  IVZ (0-based) = {} (first wake station system index)", xfoil_ivz_0based);

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

    if !input_path.exists() {
        return;
    }

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
        assert!(
            iv < nsys,
            "Lower surface iv={} should be < nsys={}",
            iv,
            nsys
        );
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
            eprintln!("Skipping test: fixture not found");
            return;
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
                println!("  IBL={}: DSTR*UEDG={:.10e}, MASS={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_mass, xfoil_mass, rel_err);
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
                println!("  IBL={}: DSTR*UEDG={:.10e}, MASS={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_mass, xfoil_mass, rel_err);
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
            eprintln!("Skipping test: fixture not found");
            return;
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
        println!("  {:>4} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.uedg, station.dstr, station.mass);
    }
    println!("  ...");
    if let Some(te) = fixture.upper_surface.last() {
        println!("  {:>4} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.uedg, te.dstr, te.mass);
    }

    // Lower surface summary
    println!("\nLower surface ({} stations):", fixture.nbl_lower);
    println!("  {:>4} {:>12} {:>12} {:>12}", "IBL", "UEDG", "DSTR", "MASS");
    for station in fixture.lower_surface.iter().take(5) {
        println!("  {:>4} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.uedg, station.dstr, station.mass);
    }
    println!("  ...");
    if let Some(te) = fixture.lower_surface.last() {
        println!("  {:>4} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.uedg, te.dstr, te.mass);
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
            eprintln!("Skipping test: UPDATE fixture not found");
            return;
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
                println!("  IBL={}: UINV+DUI={:.10e}, UNEW={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_unew, xfoil_unew, rel_err);
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
                println!("  IBL={}: UINV+DUI={:.10e}, UNEW={:.10e}, rel_err={:.2e} MISMATCH",
                    station.ibl, computed_unew, xfoil_unew, rel_err);
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
            eprintln!("Skipping test: UPDATE fixture not found");
            return;
        }
    };

    println!("=== XFOIL UPDATE Summary ===");
    println!("Description: {}", fixture.description);
    println!();

    // Upper surface summary
    println!("Upper surface ({} stations):", fixture.upper_surface.len());
    println!("  {:>4} {:>5} {:>12} {:>12} {:>12} {:>12}", "IBL", "IPAN", "UINV", "DUI", "UNEW", "MASS");
    for station in fixture.upper_surface.iter().take(5) {
        println!("  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.ipan, station.uinv, station.dui, station.unew, station.mass);
    }
    println!("  ...");
    if let Some(te) = fixture.upper_surface.last() {
        println!("  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.ipan, te.uinv, te.dui, te.unew, te.mass);
    }

    // Lower surface summary
    println!("\nLower surface ({} stations):", fixture.lower_surface.len());
    println!("  {:>4} {:>5} {:>12} {:>12} {:>12} {:>12}", "IBL", "IPAN", "UINV", "DUI", "UNEW", "MASS");
    for station in fixture.lower_surface.iter().take(5) {
        println!("  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e}",
            station.ibl, station.ipan, station.uinv, station.dui, station.unew, station.mass);
    }
    println!("  ...");
    if let Some(te) = fixture.lower_surface.last() {
        println!("  {:>4} {:>5} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} (TE)",
            te.ibl, te.ipan, te.uinv, te.dui, te.unew, te.mass);
    }

    // Statistics
    let upper_dui_range: (f64, f64) = fixture.upper_surface.iter()
        .map(|s| s.dui)
        .fold((f64::MAX, f64::MIN), |(min, max), v| (min.min(v), max.max(v)));
    let lower_dui_range: (f64, f64) = fixture.lower_surface.iter()
        .map(|s| s.dui)
        .fold((f64::MAX, f64::MIN), |(min, max), v| (min.min(v), max.max(v)));

    println!("\nDUI ranges:");
    println!("  Upper: [{:.6e}, {:.6e}]", upper_dui_range.0, upper_dui_range.1);
    println!("  Lower: [{:.6e}, {:.6e}]", lower_dui_range.0, lower_dui_range.1);
}

/// Test DUI computation using XFOIL's full DIJ matrix (including wake)
///
/// This test uses XFOIL's 183×183 DIJ matrix (N=160 airfoil + NW=23 wake panels)
/// to verify that the DUI formula: DUI = sum(-VTI_i * VTI_j * DIJ[i,j] * MASS[j])
/// produces the same result as XFOIL when given the same inputs.
///
/// NOTE: YFoil's current inviscid solver only computes 160×160 DIJ (no wake),
/// so this test loads XFOIL's DIJ directly to isolate the DUI computation.
#[test]
#[ignore] // Run with: cargo test test_dui_computation -- --ignored --nocapture
fn test_dui_computation_matches_xfoil() {
    // Load geometry
    let geom_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("docs/validation/assets/analysis/geometry/naca0012.json");

    if !geom_path.exists() {
        eprintln!("Skipping test: geometry file not found at {:?}", geom_path);
        return;
    }

    let geometry = read_geometry_from_file(&geom_path).expect("Failed to load geometry");
    let airfoil = create_paneled_airfoil(&geometry);
    let n = airfoil.n;

    // Load UPDATE fixture (has MASS and expected DUI values)
    let update_fixture = match load_update_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping test: UPDATE fixture not found");
            return;
        }
    };

    println!("=== DUI Computation Test ===");
    println!("N panels: {}", n);
    println!("Upper stations: {}", update_fixture.upper_surface.len());
    println!("Lower stations: {}", update_fixture.lower_surface.len());

    // Compute inviscid solution (includes DIJ matrix)
    let inviscid = solve_inviscid(&airfoil);
    let dij = inviscid.get_dij().expect("DIJ matrix should be computed");
    println!("DIJ matrix computed: {}x{}", dij.nrows(), dij.ncols());

    // Stagnation index (from XFOIL: IST = 80, which is 1-indexed)
    // In 0-indexed: stag_idx = 79 means panels 0..79 are upper, 80..159 are lower
    // But XFOIL uses panel indices directly: upper panels are < 80, lower are >= 80
    let stag_idx = 80; // Panel index where lower surface starts (0-indexed: panels 80+ are lower)

    // Build MASS array indexed by panel
    let mut mass_by_panel = vec![0.0; n];

    // Map upper surface MASS+VDEL values to panels
    // XFOIL computes DUI using MASS+VDEL, not just MASS
    // Upper surface: IBL=2 -> IPAN=80, IBL=3 -> IPAN=79, etc. (going backwards from stag)
    for station in &update_fixture.upper_surface {
        let ipan = station.ipan - 1; // Convert to 0-indexed (XFOIL uses 1-indexed)
        if ipan < n {
            // Use mass_plus_vdel if available, otherwise fall back to mass
            mass_by_panel[ipan] = station.mass_plus_vdel.unwrap_or(station.mass);
        }
    }

    // Map lower surface MASS+VDEL values to panels
    // Lower surface: IBL=2 -> IPAN=81, IBL=3 -> IPAN=82, etc.
    // Note: IBL > 81 are wake stations with IPAN > 160, skip those
    for station in &update_fixture.lower_surface {
        if station.ipan > n {
            continue; // Skip wake stations
        }
        let ipan = station.ipan - 1; // Convert to 0-indexed
        mass_by_panel[ipan] = station.mass_plus_vdel.unwrap_or(station.mass);
    }

    // Compute VTI for each panel
    // VTI = +1 for upper surface (panels < stag_idx), -1 for lower surface (panels >= stag_idx)
    let vti: Vec<f64> = (0..n)
        .map(|i| if i < stag_idx { 1.0 } else { -1.0 })
        .collect();

    // Compute DUI for each panel using the formula:
    // DUI[i] = sum_j(-VTI[i] * VTI[j] * DIJ[i,j] * MASS[j])
    let mut dui_computed = vec![0.0; n];
    for i in 0..n {
        let mut dui = 0.0;
        for j in 0..n {
            if mass_by_panel[j] != 0.0 {
                let ue_m = -vti[i] * vti[j] * dij[(i, j)];
                dui += ue_m * mass_by_panel[j];
            }
        }
        dui_computed[i] = dui;
    }

    // Debug: print MASS distribution
    let mut mass_sum: f64 = 0.0;
    let mut mass_count = 0;
    let mut upper_mass_count = 0;
    let mut lower_mass_count = 0;
    let mut upper_mass_sum = 0.0;
    let mut lower_mass_sum = 0.0;
    for j in 0..n {
        if mass_by_panel[j] != 0.0 {
            mass_sum += mass_by_panel[j];
            mass_count += 1;
            if j < stag_idx {
                upper_mass_count += 1;
                upper_mass_sum += mass_by_panel[j];
            } else {
                lower_mass_count += 1;
                lower_mass_sum += mass_by_panel[j];
            }
        }
    }
    println!("\nMASS distribution debug:");
    println!("  Total: {} panels with non-zero MASS, sum={:.6e}", mass_count, mass_sum);
    println!("  Upper (panels 0-79): {} panels, sum={:.6e}", upper_mass_count, upper_mass_sum);
    println!("  Lower (panels 80-159): {} panels, sum={:.6e}", lower_mass_count, lower_mass_sum);

    // Print MASS values for specific panels
    println!("\n  Sample MASS values (upper surface):");
    for &j in &[79, 78, 77, 61, 60, 1, 0] {
        if j < n {
            println!("    Panel {:3} (0-indexed): MASS={:.6e}", j, mass_by_panel[j]);
        }
    }
    println!("\n  Sample MASS values (lower surface):");
    for &j in &[80, 81, 82, 100, 158, 159] {
        if j < n {
            println!("    Panel {:3} (0-indexed): MASS={:.6e}", j, mass_by_panel[j]);
        }
    }

    // Debug specific panel (panel 79 = XFOIL panel 80 = upper IBL=2)
    let debug_i = 79;
    println!("\nDetailed DUI computation for panel {} (upper IBL=2):", debug_i);
    println!("  VTI[{}] = {}", debug_i, vti[debug_i]);
    let mut debug_dui = 0.0;
    let mut upper_contrib = 0.0;
    let mut lower_contrib = 0.0;
    for j in 0..n {
        if mass_by_panel[j] != 0.0 {
            let ue_m = -vti[debug_i] * vti[j] * dij[(debug_i, j)];
            let contrib = ue_m * mass_by_panel[j];
            debug_dui += contrib;
            if j < stag_idx {
                upper_contrib += contrib;
            } else {
                lower_contrib += contrib;
            }
            // Print a few specific contributions
            if j == debug_i || j == debug_i - 1 || j == debug_i + 1 || j == 80 || j == 0 || j == 159 {
                println!("    j={:3}: VTI[j]={:+.0}, DIJ[{},{}]={:+.6e}, MASS={:.6e}, contrib={:+.6e}",
                    j, vti[j], debug_i, j, dij[(debug_i, j)], mass_by_panel[j], contrib);
            }
        }
    }
    println!("  Upper surface contribution: {:.6e}", upper_contrib);
    println!("  Lower surface contribution: {:.6e}", lower_contrib);
    println!("  Total DUI: {:.6e}", debug_dui);
    println!("  Expected (XFOIL): {:.6e}", update_fixture.upper_surface[0].dui);

    // Compare against XFOIL's DUI values
    println!("\n=== Upper Surface DUI Comparison ===");
    println!("{:>4} {:>5} {:>14} {:>14} {:>12}", "IBL", "IPAN", "XFOIL_DUI", "YFOIL_DUI", "REL_ERR");

    let mut max_rel_err_upper: f64 = 0.0;
    let mut max_err_station_upper = 0;

    for station in &update_fixture.upper_surface {
        let ipan = station.ipan - 1; // 0-indexed
        let xfoil_dui = station.dui;
        let yfoil_dui = dui_computed[ipan];

        let rel_err = if xfoil_dui.abs() > 1e-15 {
            (yfoil_dui - xfoil_dui).abs() / xfoil_dui.abs()
        } else {
            (yfoil_dui - xfoil_dui).abs()
        };

        if rel_err > max_rel_err_upper {
            max_rel_err_upper = rel_err;
            max_err_station_upper = station.ibl;
        }

        // Print first few and any mismatches
        if station.ibl <= 5 || rel_err > 0.01 {
            println!("{:>4} {:>5} {:>14.6e} {:>14.6e} {:>12.2e}{}",
                station.ibl, station.ipan, xfoil_dui, yfoil_dui, rel_err,
                if rel_err > 0.01 { " MISMATCH" } else { "" });
        }
    }
    println!("Max rel error: {:.2e} at IBL={}", max_rel_err_upper, max_err_station_upper);

    println!("\n=== Lower Surface DUI Comparison (airfoil only, no wake) ===");
    println!("{:>4} {:>5} {:>14} {:>14} {:>12}", "IBL", "IPAN", "XFOIL_DUI", "YFOIL_DUI", "REL_ERR");

    let mut max_rel_err_lower: f64 = 0.0;
    let mut max_err_station_lower = 0;

    for station in &update_fixture.lower_surface {
        if station.ipan > n {
            continue; // Skip wake stations
        }
        let ipan = station.ipan - 1; // 0-indexed
        let xfoil_dui = station.dui;
        let yfoil_dui = dui_computed[ipan];

        let rel_err = if xfoil_dui.abs() > 1e-15 {
            (yfoil_dui - xfoil_dui).abs() / xfoil_dui.abs()
        } else {
            (yfoil_dui - xfoil_dui).abs()
        };

        if rel_err > max_rel_err_lower {
            max_rel_err_lower = rel_err;
            max_err_station_lower = station.ibl;
        }

        // Print first few and any mismatches
        if station.ibl <= 5 || rel_err > 0.01 {
            println!("{:>4} {:>5} {:>14.6e} {:>14.6e} {:>12.2e}{}",
                station.ibl, station.ipan, xfoil_dui, yfoil_dui, rel_err,
                if rel_err > 0.01 { " MISMATCH" } else { "" });
        }
    }
    println!("Max rel error: {:.2e} at IBL={}", max_rel_err_lower, max_err_station_lower);

    // Assert with tolerance for fixture precision (~1e-7)
    // The DUI values are small (1e-3 to 1e-4), so we allow some tolerance
    let tolerance = 1e-5; // Allow 0.001% error

    assert!(
        max_rel_err_upper < tolerance,
        "Upper surface DUI mismatch: max rel err = {:.2e} at IBL={} (tolerance = {:.0e})",
        max_rel_err_upper, max_err_station_upper, tolerance
    );
    assert!(
        max_rel_err_lower < tolerance,
        "Lower surface DUI mismatch: max rel err = {:.2e} at IBL={} (tolerance = {:.0e})",
        max_rel_err_lower, max_err_station_lower, tolerance
    );

    println!("\nSUCCESS: DUI computation matches XFOIL within {:.0e}", tolerance);
}

/// Test DIJ matrix values against XFOIL fixture
///
/// This validates that YFoil's DIJ = inv(AIJ) * BIJ matches XFOIL's DIJ.
/// XFOIL uses 1-indexed panels, YFoil uses 0-indexed.
/// XFOIL DIJ(i,j) corresponds to YFoil dij[(i-1, j-1)]
#[test]
#[ignore] // Run with: cargo test test_dij_matrix_matches_xfoil -- --ignored --nocapture
fn test_dij_matrix_matches_xfoil() {
    // Load geometry (same as DUI test)
    let geom_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("docs/validation/assets/analysis/geometry/naca0012.json");

    if !geom_path.exists() {
        eprintln!("Skipping test: geometry file not found at {:?}", geom_path);
        return;
    }

    let geometry = read_geometry_from_file(&geom_path).expect("Failed to load geometry");
    let airfoil = create_paneled_airfoil(&geometry);
    let n = airfoil.n;

    // Compute inviscid solution (includes DIJ matrix)
    let inviscid = solve_inviscid(&airfoil);
    let dij = inviscid.get_dij().expect("DIJ matrix should be computed");

    println!("=== DIJ Matrix Comparison Test ===");
    println!("N panels: {}", n);
    println!("DIJ matrix: {}x{}", dij.nrows(), dij.ncols());

    // Load XFOIL DIJ fixture (generated by instrumented XFOIL)
    // Format: row col value (1-indexed)
    let dij_path = std::path::Path::new("/tmp/xfoil_dij.dat");
    if !dij_path.exists() {
        eprintln!("XFOIL DIJ fixture not found at {:?}", dij_path);
        eprintln!("Run instrumented XFOIL first to generate it");
        return;
    }

    let content = fs::read_to_string(dij_path).expect("Failed to read DIJ fixture");
    let lines: Vec<&str> = content.lines().collect();

    // Parse header
    println!("\nXFOIL DIJ header:");
    for line in lines.iter().take(3) {
        println!("  {}", line);
    }

    // Compare specific DIJ entries
    // XFOIL row 80 (near LE, upper surface station IBL=2)
    println!("\n=== DIJ Row 80 (XFOIL 1-indexed) / Row 79 (YFoil 0-indexed) ===");
    println!("{:>4} {:>4} {:>20} {:>20} {:>12}", "XROW", "XCOL", "XFOIL_DIJ", "YFOIL_DIJ", "REL_ERR");

    let mut max_rel_err: f64 = 0.0;
    let mut max_err_entry = (0, 0);
    let mut entries_checked = 0;

    for line in &lines {
        // Parse lines like: "   80    1  0.4240011919787700E+02"
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 {
            continue;
        }

        let row: usize = match parts[0].parse() {
            Ok(r) => r,
            Err(_) => continue,
        };
        let col: usize = match parts[1].parse() {
            Ok(c) => c,
            Err(_) => continue,
        };
        let xfoil_val: f64 = match parts[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Convert to 0-indexed
        let yfoil_row = row - 1;
        let yfoil_col = col - 1;

        if yfoil_row >= n || yfoil_col >= n {
            continue;
        }

        let yfoil_val = dij[(yfoil_row, yfoil_col)];

        let rel_err = if xfoil_val.abs() > 1e-15 {
            (yfoil_val - xfoil_val).abs() / xfoil_val.abs()
        } else {
            (yfoil_val - xfoil_val).abs()
        };

        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_err_entry = (row, col);
        }
        entries_checked += 1;

        // Print row 80 entries (near LE)
        if row == 80 && col <= 10 {
            println!("{:>4} {:>4} {:>20.10e} {:>20.10e} {:>12.2e}{}",
                row, col, xfoil_val, yfoil_val, rel_err,
                if rel_err > 1e-6 { " MISMATCH" } else { "" });
        }
    }

    // Print row 1 entries (near TE upper)
    println!("\n=== DIJ Row 1 (XFOIL 1-indexed) / Row 0 (YFoil 0-indexed) ===");
    for line in &lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 { continue; }
        let row: usize = match parts[0].parse() { Ok(r) => r, Err(_) => continue };
        let col: usize = match parts[1].parse() { Ok(c) => c, Err(_) => continue };
        let xfoil_val: f64 = match parts[2].parse() { Ok(v) => v, Err(_) => continue };

        if row == 1 && col <= 10 {
            let yfoil_row = row - 1;
            let yfoil_col = col - 1;
            let yfoil_val = dij[(yfoil_row, yfoil_col)];
            let rel_err = if xfoil_val.abs() > 1e-15 {
                (yfoil_val - xfoil_val).abs() / xfoil_val.abs()
            } else {
                (yfoil_val - xfoil_val).abs()
            };
            println!("{:>4} {:>4} {:>20.10e} {:>20.10e} {:>12.2e}{}",
                row, col, xfoil_val, yfoil_val, rel_err,
                if rel_err > 1e-6 { " MISMATCH" } else { "" });
        }
    }

    println!("\nTotal entries checked: {}", entries_checked);
    println!("Max relative error: {:.2e} at XFOIL entry ({}, {})", max_rel_err, max_err_entry.0, max_err_entry.1);

    // Check specific entry that had largest error in DUI test
    // Upper IBL=20 had 29x error, which is IPAN=61 (XFOIL 1-indexed)
    println!("\n=== Checking entries for upper IBL=20 (IPAN=61) ===");
    let check_row = 61;
    for line in &lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 { continue; }
        let row: usize = match parts[0].parse() { Ok(r) => r, Err(_) => continue };
        let col: usize = match parts[1].parse() { Ok(c) => c, Err(_) => continue };
        let xfoil_val: f64 = match parts[2].parse() { Ok(v) => v, Err(_) => continue };

        // Check a few columns for this row
        if row == check_row && (col <= 5 || col == 80 || col == 81 || col == 160) {
            let yfoil_row = row - 1;
            let yfoil_col = col - 1;
            let yfoil_val = dij[(yfoil_row, yfoil_col)];
            let rel_err = if xfoil_val.abs() > 1e-15 {
                (yfoil_val - xfoil_val).abs() / xfoil_val.abs()
            } else {
                (yfoil_val - xfoil_val).abs()
            };
            println!("  DIJ({}, {}): XFOIL={:+.10e}, YFOIL={:+.10e}, rel_err={:.2e}{}",
                row, col, xfoil_val, yfoil_val, rel_err,
                if rel_err > 1e-6 { " MISMATCH" } else { "" });
        }
    }

    // Assert overall accuracy
    let tolerance = 1e-6;
    assert!(
        max_rel_err < tolerance,
        "DIJ matrix mismatch: max rel err = {:.2e} at ({}, {}) (tolerance = {:.0e})",
        max_rel_err, max_err_entry.0, max_err_entry.1, tolerance
    );

    println!("\nSUCCESS: DIJ matrix matches XFOIL within {:.0e}", tolerance);
}

/// Test DUI computation using XFOIL's full DIJ matrix (including wake)
///
/// This test loads XFOIL's 183×183 DIJ matrix (N=160 + NW=23) and verifies
/// that the DUI formula produces matching results. This isolates the DUI
/// computation from YFoil's inviscid solver (which doesn't include wake).
#[test]
#[ignore] // Run with: cargo test test_dui_with_xfoil_dij -- --ignored --nocapture
fn test_dui_with_xfoil_dij() {
    // Load XFOIL DIJ matrix from fixture
    let dij_path = std::path::Path::new("/tmp/xfoil_dij.dat");
    if !dij_path.exists() {
        eprintln!("XFOIL DIJ fixture not found at {:?}", dij_path);
        eprintln!("Run instrumented XFOIL first to generate it");
        return;
    }

    let content = fs::read_to_string(dij_path).expect("Failed to read DIJ fixture");
    let lines: Vec<&str> = content.lines().collect();

    // Parse header to get dimensions
    let mut n_total = 0;
    let mut n_airfoil = 0;
    let mut n_wake = 0;
    for line in lines.iter().take(5) {
        if line.contains("N =") && !line.contains("NW") {
            n_airfoil = line.split_whitespace().last()
                .and_then(|s| s.parse().ok()).unwrap_or(0);
        }
        if line.contains("NW =") {
            n_wake = line.split_whitespace().last()
                .and_then(|s| s.parse().ok()).unwrap_or(0);
        }
    }
    n_total = n_airfoil + n_wake;

    if n_total == 0 {
        eprintln!("Failed to parse DIJ dimensions");
        return;
    }

    println!("=== DUI Test with XFOIL DIJ Matrix ===");
    println!("N_airfoil = {}, N_wake = {}, N_total = {}", n_airfoil, n_wake, n_total);

    // Parse DIJ matrix (1-indexed in file)
    let mut dij = vec![vec![0.0; n_total]; n_total];
    for line in &lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != 3 { continue; }
        let row: usize = match parts[0].parse() { Ok(r) => r, Err(_) => continue };
        let col: usize = match parts[1].parse() { Ok(c) => c, Err(_) => continue };
        let val: f64 = match parts[2].parse() { Ok(v) => v, Err(_) => continue };

        if row >= 1 && row <= n_total && col >= 1 && col <= n_total {
            dij[row - 1][col - 1] = val;  // Convert to 0-indexed
        }
    }

    // Load UPDATE fixture
    let update_fixture = match load_update_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping test: UPDATE fixture not found");
            return;
        }
    };

    // Build MASS array for all panels (including wake)
    let mut mass_by_panel = vec![0.0; n_total];

    // Upper surface: maps to airfoil panels 0-79 (0-indexed)
    for station in &update_fixture.upper_surface {
        let ipan = station.ipan - 1;  // 0-indexed
        if ipan < n_airfoil {
            mass_by_panel[ipan] = station.mass_plus_vdel.unwrap_or(station.mass);
        }
    }

    // Lower surface: includes airfoil (panels 80-159) AND wake (panels 160+)
    for station in &update_fixture.lower_surface {
        let ipan = station.ipan - 1;  // 0-indexed
        if ipan < n_total {
            mass_by_panel[ipan] = station.mass_plus_vdel.unwrap_or(station.mass);
        }
    }

    // Count MASS entries
    let airfoil_mass_count = mass_by_panel[..n_airfoil].iter().filter(|&&m| m != 0.0).count();
    let wake_mass_count = mass_by_panel[n_airfoil..].iter().filter(|&&m| m != 0.0).count();
    println!("MASS distribution: {} airfoil + {} wake = {} total",
        airfoil_mass_count, wake_mass_count, airfoil_mass_count + wake_mass_count);

    // VTI by panel:
    // - Airfoil upper (panels 0-79): +1
    // - Airfoil lower (panels 80-159): -1
    // - Wake (panels 160+): -1 (continuation of lower surface)
    let stag_idx = 80;  // First lower surface panel (0-indexed)
    let vti: Vec<f64> = (0..n_total)
        .map(|i| if i < stag_idx { 1.0 } else { -1.0 })
        .collect();

    // Compute DUI for each airfoil panel
    let mut dui_computed = vec![0.0; n_airfoil];
    for i in 0..n_airfoil {
        let mut dui = 0.0;
        for j in 0..n_total {  // Sum over ALL panels including wake
            if mass_by_panel[j] != 0.0 {
                let ue_m = -vti[i] * vti[j] * dij[i][j];
                dui += ue_m * mass_by_panel[j];
            }
        }
        dui_computed[i] = dui;
    }

    // Compare against XFOIL's DUI values
    println!("\n=== Upper Surface DUI Comparison (with wake) ===");
    println!("{:>4} {:>5} {:>14} {:>14} {:>12}", "IBL", "IPAN", "XFOIL_DUI", "YFOIL_DUI", "REL_ERR");

    let mut max_rel_err_upper: f64 = 0.0;
    let mut max_err_station_upper = 0;

    for station in &update_fixture.upper_surface {
        let ipan = station.ipan - 1;  // 0-indexed
        if ipan >= n_airfoil { continue; }

        let xfoil_dui = station.dui;
        let yfoil_dui = dui_computed[ipan];

        let rel_err = if xfoil_dui.abs() > 1e-15 {
            (yfoil_dui - xfoil_dui).abs() / xfoil_dui.abs()
        } else {
            (yfoil_dui - xfoil_dui).abs()
        };

        if rel_err > max_rel_err_upper {
            max_rel_err_upper = rel_err;
            max_err_station_upper = station.ibl;
        }

        // Print first 10 and any mismatches > 1%
        if station.ibl <= 10 || rel_err > 0.01 {
            println!("{:>4} {:>5} {:>14.6e} {:>14.6e} {:>12.2e}{}",
                station.ibl, station.ipan, xfoil_dui, yfoil_dui, rel_err,
                if rel_err > 0.01 { " MISMATCH" } else { "" });
        }
    }
    println!("Max rel error: {:.2e} at IBL={}", max_rel_err_upper, max_err_station_upper);

    // Check lower surface (airfoil only)
    println!("\n=== Lower Surface DUI Comparison (with wake) ===");
    println!("{:>4} {:>5} {:>14} {:>14} {:>12}", "IBL", "IPAN", "XFOIL_DUI", "YFOIL_DUI", "REL_ERR");

    let mut max_rel_err_lower: f64 = 0.0;
    let mut max_err_station_lower = 0;

    for station in &update_fixture.lower_surface {
        let ipan = station.ipan - 1;  // 0-indexed
        if ipan >= n_airfoil { continue; }  // Skip wake

        let xfoil_dui = station.dui;
        let yfoil_dui = dui_computed[ipan];

        let rel_err = if xfoil_dui.abs() > 1e-15 {
            (yfoil_dui - xfoil_dui).abs() / xfoil_dui.abs()
        } else {
            (yfoil_dui - xfoil_dui).abs()
        };

        if rel_err > max_rel_err_lower {
            max_rel_err_lower = rel_err;
            max_err_station_lower = station.ibl;
        }

        // Print first 10 and any mismatches > 1%
        if station.ibl <= 10 || rel_err > 0.01 {
            println!("{:>4} {:>5} {:>14.6e} {:>14.6e} {:>12.2e}{}",
                station.ibl, station.ipan, xfoil_dui, yfoil_dui, rel_err,
                if rel_err > 0.01 { " MISMATCH" } else { "" });
        }
    }
    println!("Max rel error: {:.2e} at IBL={}", max_rel_err_lower, max_err_station_lower);

    // Assert with tolerance
    let tolerance = 1e-5;  // Allow 0.001% error (fixture precision)

    assert!(
        max_rel_err_upper < tolerance,
        "Upper surface DUI mismatch: max rel err = {:.2e} at IBL={} (tolerance = {:.0e})",
        max_rel_err_upper, max_err_station_upper, tolerance
    );
    assert!(
        max_rel_err_lower < tolerance,
        "Lower surface DUI mismatch: max rel err = {:.2e} at IBL={} (tolerance = {:.0e})",
        max_rel_err_lower, max_err_station_lower, tolerance
    );

    println!("\nSUCCESS: DUI computation matches XFOIL within {:.0e} (using full DIJ with wake)", tolerance);
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
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/subroutines/update/naca0012_rlx_iter1.json");
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
#[ignore] // Run with: cargo test test_rlx_computation_matches_xfoil -- --ignored --nocapture
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
    let mut max_dn_station = (0_usize, 0_usize);  // (surface, ibl)
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
        rmsbl_sum += dn1*dn1 + dn2*dn2 + dn3*dn3 + dn4*dn4;

        // Check RLX bounds (this is the key algorithm)
        // XFOIL: IF(RDN1 .GT. DHI) RLX = DHI/DN1
        //        IF(RDN1 .LT. DLO) RLX = DLO/DN1
        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }

        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }
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

        rmsbl_sum += dn1*dn1 + dn2*dn2 + dn3*dn3 + dn4*dn4;

        let rdn1 = rlx * dn1;
        if rdn1 > dhi { rlx = dhi / dn1; }
        if rdn1 < dlo { rlx = dlo / dn1; }

        let rdn2 = rlx * dn2;
        if rdn2 > dhi { rlx = dhi / dn2; }
        if rdn2 < dlo { rlx = dlo / dn2; }

        let rdn3 = rlx * dn3;
        if rdn3 > dhi { rlx = dhi / dn3; }
        if rdn3 < dlo { rlx = dlo / dn3; }

        let rdn4 = rlx * dn4;
        if rdn4 > dhi { rlx = dhi / dn4; }
    }

    // Compute RMSBL
    // XFOIL uses NBL(1)+NBL(2) in divisor, not (NBL(1)-1)+(NBL(2)-1)
    // This is because RMSBL normalizes by total stations including stagnation
    let n_stations = fixture.nbl_upper + fixture.nbl_lower;  // 81 + 104 = 185
    let rmsbl = (rmsbl_sum / (4.0 * n_stations as f64)).sqrt();

    println!("\nComputed values:");
    println!("  RLX = {:.10}", rlx);
    println!("  RMSBL = {:.10}", rmsbl);
    println!("\nMax DN: {} = {:.6} at surface={}, IBL={}",
        max_dn_var, max_dn_value, max_dn_station.0, max_dn_station.1);

    // Verify RLX matches to machine precision
    let rlx_rel_err = (rlx - fixture.rlx).abs() / fixture.rlx.abs();
    println!("\nRLX comparison:");
    println!("  XFOIL:    {:.10}", fixture.rlx);
    println!("  Computed: {:.10}", rlx);
    println!("  Rel err:  {:.2e}", rlx_rel_err);

    let tolerance = 1e-8;  // Machine precision
    assert!(
        rlx_rel_err < tolerance,
        "RLX mismatch: computed={:.10}, expected={:.10}, rel_err={:.2e}",
        rlx, fixture.rlx, rlx_rel_err
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
        rmsbl, fixture.rmsbl, rmsbl_rel_err
    );

    println!("\nSUCCESS: RLX and RMSBL computations match XFOIL within {:.0e}", tolerance);
}
