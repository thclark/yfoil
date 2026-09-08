//! Compare yfoil's BIJ matrix against XFOIL's output
//!
//! Run XFOIL debug version first to generate /tmp/xfoil_bij.dat

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::geometry::{panel_foil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    println!("=== BIJ Matrix Comparison: yfoil vs XFOIL ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;

    println!("Airfoil: {} panels", n);

    // Solve inviscid to get BIJ
    let inviscid = solve_inviscid(&airfoil);
    let bij = inviscid.get_bij().expect("BIJ should be computed");

    // Print sample yfoil BIJ values (matching XFOIL debug output indices)
    // Note: XFOIL uses 1-based indexing, Rust uses 0-based
    println!("\nyfoil BIJ values:");
    println!("BIJ[0,0] (XFOIL 1,1) = {:e}", bij[(0, 0)]);
    println!("BIJ[0,1] (XFOIL 1,2) = {:e}", bij[(0, 1)]);
    println!("BIJ[1,0] (XFOIL 2,1) = {:e}", bij[(1, 0)]);
    println!("BIJ[39,39] (XFOIL 40,40) = {:e}", bij[(39, 39)]);
    println!("BIJ[39,40] (XFOIL 40,41) = {:e}", bij[(39, 40)]);
    println!("BIJ[79,79] (XFOIL 80,80) = {:e}", bij[(79, 79)]);
    println!("BIJ[79,78] (XFOIL 80,79) = {:e}", bij[(79, 78)]);

    // Read XFOIL's BIJ values from file
    println!("\nXFOIL BIJ values (from debug output):");
    println!("BIJ(1,1) = -2.207e-02");
    println!("BIJ(1,2) = 9.218e-05");
    println!("BIJ(2,1) = -4.500e-01");
    println!("BIJ(40,40) = 3.757e-01");
    println!("BIJ(40,41) = 6.297e-02");
    println!("BIJ(80,80) = 3.898e-01");
    println!("BIJ(80,79) = 6.964e-02");

    // Load XFOIL's BIJ file and compare specific values
    if let Ok(file) = File::open("/tmp/xfoil_bij.dat") {
        println!("\n=== Loading full comparison from /tmp/xfoil_bij.dat ===");
        let reader = BufReader::new(file);

        let mut max_abs_diff: f64 = 0.0;
        let mut max_diff_i = 0;
        let mut max_diff_j = 0;
        let mut max_rel_diff: f64 = 0.0;
        let mut count = 0;

        for line in reader.lines() {
            if let Ok(line) = line {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 3 {
                    if let (Ok(i), Ok(j), Ok(xfoil_val)) = (
                        parts[0].parse::<usize>(),
                        parts[1].parse::<usize>(),
                        parts[2].parse::<f64>(),
                    ) {
                        // Convert to 0-based indexing
                        let ri = i - 1;
                        let rj = j - 1;
                        if ri < n && rj < n {
                            let yfoil_val = bij[(ri, rj)];
                            let abs_diff = (yfoil_val - xfoil_val).abs();

                            if abs_diff > max_abs_diff {
                                max_abs_diff = abs_diff;
                                max_diff_i = i;
                                max_diff_j = j;
                            }

                            if xfoil_val.abs() > 1e-10 {
                                let rel_diff = abs_diff / xfoil_val.abs();
                                if rel_diff > max_rel_diff {
                                    max_rel_diff = rel_diff;
                                }
                            }
                            count += 1;
                        }
                    }
                }
            }
        }

        println!("Compared {} values", count);
        println!("Max absolute difference: {:e} at ({}, {})", max_abs_diff, max_diff_i, max_diff_j);
        println!("Max relative difference: {:.2}%", max_rel_diff * 100.0);

        // Print the values at max diff location
        let ri = max_diff_i - 1;
        let rj = max_diff_j - 1;
        println!("\nAt max diff location ({}, {}):", max_diff_i, max_diff_j);
        println!("  yfoil BIJ: {:e}", bij[(ri, rj)]);

    } else {
        println!("\nCould not open /tmp/xfoil_bij.dat for comparison");
    }

    // Also compare DIJ
    let dij = inviscid.get_dij().expect("DIJ should be computed");
    println!("\n=== DIJ comparison ===");
    println!("\nyfoil DIJ values:");
    println!("DIJ[0,0] (XFOIL 1,1) = {:e}", dij[(0, 0)]);
    println!("DIJ[0,1] (XFOIL 1,2) = {:e}", dij[(0, 1)]);
    println!("DIJ[1,0] (XFOIL 2,1) = {:e}", dij[(1, 0)]);
    println!("DIJ[39,39] (XFOIL 40,40) = {:e}", dij[(39, 39)]);
    println!("DIJ[79,79] (XFOIL 80,80) = {:e}", dij[(79, 79)]);

    println!("\nXFOIL DIJ values (from debug output):");
    println!("DIJ(1,1) = -4.951e+01");
    println!("DIJ(1,2) = 3.902e+01");
    println!("DIJ(2,1) = 9.657e+01");
    println!("DIJ(40,40) = -7.455e+01");
    println!("DIJ(80,80) = -6.741e+02");
}
