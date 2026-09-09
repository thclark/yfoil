//! Check DIJ matrix against XFOIL

use std::fs;
use yfoil::geometry::{panel_foil, Geometry};
use yfoil::panel::solve_inviscid;

fn main() {
    // Load the same geometry used for XFOIL
    let json_path = "/tmp/naca0012.json";
    let json_str = fs::read_to_string(json_path).expect("Failed to read JSON file");
    let geom: Geometry = serde_json::from_str(&json_str).expect("Failed to parse JSON");
    let airfoil = panel_foil(&geom);

    println!("=== yFoil DIJ Matrix Check ===\n");
    println!("N = {}", airfoil.n);

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);

    // Get DIJ matrix
    let dij = inviscid.get_dij().expect("DIJ matrix not computed");

    println!("\nDIJ matrix shape: {} x {}", dij.nrows(), dij.ncols());

    // Print some DIJ values for comparison with XFOIL
    println!("\n=== Sample DIJ values (first 5x5) ===");
    for i in 0..5.min(dij.nrows()) {
        for j in 0..5.min(dij.ncols()) {
            println!("DIJ[{},{}] = {:24.16e}", i+1, j+1, dij[(i, j)]);
        }
    }

    // Print diagonal values
    println!("\n=== DIJ diagonal values ===");
    for i in 0..dij.nrows().min(10) {
        println!("DIJ[{},{}] = {:24.16e}", i+1, i+1, dij[(i, i)]);
    }

    // Check row sums to see overall influence pattern
    println!("\n=== Row sums (total influence on each panel) ===");
    for i in 0..dij.nrows().min(10) {
        let row_sum: f64 = (0..dij.ncols()).map(|j| dij[(i, j)]).sum();
        println!("Row {} sum = {:+.6e}", i+1, row_sum);
    }

    // Print a specific value near stagnation
    let stag = airfoil.n / 2;  // Approximate stagnation
    println!("\n=== DIJ near stagnation (index {}) ===", stag);
    for i in (stag-2).max(0)..(stag+3).min(dij.nrows()) {
        for j in (stag-2).max(0)..(stag+3).min(dij.ncols()) {
            println!("DIJ[{},{}] = {:+.6e}", i+1, j+1, dij[(i, j)]);
        }
    }

    // Get BIJ matrix (= -DZDM) for comparison with XFOIL
    if let Some(bij) = inviscid.get_bij() {
        println!("\n=== BIJ (=-DZDM) for I=1 (first control point) ===");
        println!("Compare with XFOIL's /tmp/xfoil_dzdm_i1.dat");
        for j in 0..10.min(bij.ncols()) {
            // BIJ = -DZDM, so DZDM = -BIJ
            let dzdm = -bij[(0, j)];
            println!("  DZDM[{:3}] = {:+24.16e}", j+1, dzdm);
        }
        println!("  ...");

        // Print BIJ diagonal for comparison
        println!("\n=== BIJ diagonal (first 5) ===");
        for i in 0..5.min(bij.nrows()).min(bij.ncols()) {
            println!("  BIJ[{},{}] = {:+24.16e}", i+1, i+1, bij[(i, i)]);
        }
    } else {
        println!("\n=== BIJ matrix not available ===");
    }
}
