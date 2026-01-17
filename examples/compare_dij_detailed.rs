//! Detailed DIJ comparison between yfoil and XFOIL
//!
//! Run XFOIL debug version first to generate /tmp/xfoil_dij.dat

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    println!("=== Detailed DIJ Comparison ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;
    let le = airfoil.le_index;

    println!("Airfoil: {} panels, LE at index {}", n, le);

    let inviscid = solve_inviscid(&airfoil);
    let dij = inviscid.get_dij().expect("DIJ should be computed");

    // Load XFOIL's DIJ values
    let file = File::open("/tmp/xfoil_dij.dat").expect("Run XFOIL debug first");
    let reader = BufReader::new(file);

    let mut xfoil_dij = vec![vec![0.0; n]; n];
    for line in reader.lines() {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Ok(i), Ok(j), Ok(val)) = (
                    parts[0].parse::<usize>(),
                    parts[1].parse::<usize>(),
                    parts[2].parse::<f64>(),
                ) {
                    if i <= n && j <= n {
                        xfoil_dij[i - 1][j - 1] = val;
                    }
                }
            }
        }
    }

    // Compare diagonal elements
    println!("\nDiagonal comparison (DIJ[i,i]):");
    println!("{:>6} {:>12} {:>12} {:>10} {:>8}", "i", "yfoil", "XFOIL", "ratio", "region");
    for i in [0, 1, 2, 3, 10, 20, 40, 60, le - 1, le, le + 1, n - 3, n - 2, n - 1] {
        if i < n {
            let yf = dij[(i, i)];
            let xf = xfoil_dij[i][i];
            let ratio = if xf.abs() > 1e-10 { yf / xf } else { f64::NAN };
            let region = if i <= le { "upper" } else { "lower" };
            println!("{:6} {:12.4} {:12.4} {:10.4} {:>8}", i, yf, xf, ratio, region);
        }
    }

    // Compare first few rows/columns
    println!("\nRow 0 (TE region, upper surface):");
    println!("{:>6} {:>12} {:>12} {:>10}", "j", "yfoil", "XFOIL", "ratio");
    for j in 0..10.min(n) {
        let yf = dij[(0, j)];
        let xf = xfoil_dij[0][j];
        let ratio = if xf.abs() > 1e-10 { yf / xf } else { f64::NAN };
        println!("{:6} {:12.4} {:12.4} {:10.4}", j, yf, xf, ratio);
    }

    println!("\nRow 40 (middle upper surface):");
    println!("{:>6} {:>12} {:>12} {:>10}", "j", "yfoil", "XFOIL", "ratio");
    for j in 35..45.min(n) {
        let yf = dij[(40, j)];
        let xf = xfoil_dij[40][j];
        let ratio = if xf.abs() > 1e-10 { yf / xf } else { f64::NAN };
        println!("{:6} {:12.4} {:12.4} {:10.4}", j, yf, xf, ratio);
    }

    println!("\nRow 79/80 (LE region):");
    println!("{:>6} {:>12} {:>12} {:>10}", "j", "yfoil", "XFOIL", "ratio");
    for j in 75..85.min(n) {
        let yf = dij[(80, j)];
        let xf = xfoil_dij[80][j];
        let ratio = if xf.abs() > 1e-10 { yf / xf } else { f64::NAN };
        println!("{:6} {:12.4} {:12.4} {:10.4}", j, yf, xf, ratio);
    }

    // Overall statistics
    println!("\n=== Overall Statistics ===");
    let mut sum_abs_diff = 0.0;
    let mut max_abs_diff = 0.0;
    let mut max_i = 0;
    let mut max_j = 0;
    let mut count_good = 0;
    let mut count_bad = 0;

    for i in 0..n {
        for j in 0..n {
            let yf = dij[(i, j)];
            let xf = xfoil_dij[i][j];
            let abs_diff = (yf - xf).abs();
            sum_abs_diff += abs_diff;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
                max_i = i;
                max_j = j;
            }

            // Count how many are within 10% (relative to XFOIL)
            if xf.abs() > 1e-6 {
                let rel_diff = abs_diff / xf.abs();
                if rel_diff < 0.1 {
                    count_good += 1;
                } else {
                    count_bad += 1;
                }
            }
        }
    }

    println!("Mean absolute difference: {:e}", sum_abs_diff / (n * n) as f64);
    println!("Max absolute difference: {:e} at ({}, {})", max_abs_diff, max_i, max_j);
    println!("Values within 10% of XFOIL: {} ({:.1}%)", count_good, 100.0 * count_good as f64 / (n * n) as f64);
    println!("Values outside 10%: {} ({:.1}%)", count_bad, 100.0 * count_bad as f64 / (n * n) as f64);

    // Check which rows have significant differences
    println!("\n=== Row-by-row max difference ===");
    println!("{:>6} {:>12} {:>6} {:>8}", "row", "max_diff", "at_col", "region");
    for i in 0..n {
        let mut max_diff_row = 0.0;
        let mut max_col = 0;
        for j in 0..n {
            let diff = (dij[(i, j)] - xfoil_dij[i][j]).abs();
            if diff > max_diff_row {
                max_diff_row = diff;
                max_col = j;
            }
        }
        if max_diff_row > 10.0 {  // Only show rows with significant differences
            let region = if i <= le { "upper" } else { "lower" };
            println!("{:6} {:12.4} {:6} {:>8}", i, max_diff_row, max_col, region);
        }
    }
}
