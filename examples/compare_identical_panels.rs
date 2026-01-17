//! Compare yfoil and XFOIL using identical panel coordinates
//!
//! This example loads the exact paneled coordinates that XFOIL used,
//! runs yfoil's BL solver on them, and compares results.
//!
//! First run XFOIL with:
//!   xfoil < /tmp/xfoil_repanel.txt
//!
//! This generates:
//!   /tmp/xfoil_naca0012_paneled.dat - The paneled coordinates
//!   /tmp/xfoil_bl_paneled.dat - XFOIL's BL solution

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::bl::{march_newton, FlowConditions, NewtonConfig};
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, find_stagnation_point};

fn main() {
    println!("=== Loading XFOIL-paneled Coordinates ===");

    // Load the exact coordinates that XFOIL used
    let (name, geom) = read_dat_file("/tmp/xfoil_naca0012_paneled.dat")
        .expect("Failed to load XFOIL paneled coordinates. Run XFOIL first.");

    println!("Loaded: {}", name);
    println!("Number of points: {}", geom.x_c.len());

    // Show first few coordinates
    println!("\nFirst 5 coordinates:");
    for i in 0..5.min(geom.x_c.len()) {
        println!("  [{:3}] x = {:10.6}, y = {:10.6}", i, geom.x_c[i], geom.y_c[i]);
    }

    // Create paneled airfoil from these exact coordinates
    let airfoil = create_paneled_airfoil(&geom);
    println!("\nAirfoil properties:");
    println!("  Chord: {:.6}", airfoil.chord);
    println!("  LE index: {}", airfoil.le_index);
    println!("  LE arc length: {:.6}", airfoil.sle);
    println!("  Sharp TE: {}", airfoil.sharp_te);

    // Solve inviscid flow
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    println!("\nStagnation point: index = {}, x = {:.4}", stag, airfoil.x[stag]);

    // Extract upper surface
    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);

    println!("\nUpper surface:");
    println!("  Points: {}", x_upper.len());
    println!("  x range: {:.4} to {:.4}", x_upper.first().unwrap(), x_upper.last().unwrap());
    println!("  s range: {:.4} to {:.4}", s_upper.first().unwrap(), s_upper.last().unwrap());

    // Run BL solver
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    // Find transition
    let trans_idx = bl_upper.iter().position(|r| r.n_amp >= cond.ncrit);
    let trans_xc = trans_idx.map(|i| x_upper[i]).unwrap_or(f64::INFINITY);
    println!("\nyfoil transition: x/c = {:.4}", trans_xc);

    // Load XFOIL BL dump for comparison
    println!("\n=== Loading XFOIL BL Data ===");
    let file = File::open("/tmp/xfoil_bl_paneled.dat").expect("Failed to open XFOIL BL dump");
    let reader = BufReader::new(file);

    let mut xfoil_data: Vec<(f64, f64, f64, f64, f64, f64, f64)> = Vec::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 8 {
            let s: f64 = parts[0].parse().unwrap_or(0.0);
            let x: f64 = parts[1].parse().unwrap_or(0.0);
            let y: f64 = parts[2].parse().unwrap_or(0.0);
            let ue: f64 = parts[3].parse().unwrap_or(0.0);
            let _dstar: f64 = parts[4].parse().unwrap_or(0.0);
            let theta: f64 = parts[5].parse().unwrap_or(0.0);
            let cf: f64 = parts[6].parse().unwrap_or(0.0);
            let h: f64 = parts[7].parse().unwrap_or(0.0);

            if y >= 0.0 {
                xfoil_data.push((s, x, ue, theta, cf, h, y));
            }
        }
    }

    // Sort by x ascending (LE to TE)
    xfoil_data.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Find stagnation point in XFOIL data (min |Ue|)
    let stag_idx = xfoil_data
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.2.abs().partial_cmp(&b.2.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    let xfoil_upper: Vec<_> = xfoil_data[stag_idx..].to_vec();

    println!("XFOIL upper surface points: {}", xfoil_upper.len());

    // Check station by station - find XFOIL transition
    let mut xfoil_trans_x = 1.0;
    for i in 1..xfoil_upper.len() {
        // Transition occurs where H jumps suddenly or Cf drops
        if xfoil_upper[i].5 > 2.8 && xfoil_upper[i - 1].5 < 2.5 {
            xfoil_trans_x = xfoil_upper[i].1;
            break;
        }
    }
    println!("XFOIL transition (from H jump): x/c ≈ {:.4}", xfoil_trans_x);

    // Compare panel coordinates
    println!("\n=== Panel Coordinate Verification ===");
    println!(
        "{:>6} {:>12} {:>12} {:>12} {:>12}",
        "idx", "x_xfoil", "x_yfoil", "y_xfoil", "y_yfoil"
    );
    println!("{:-<60}", "");

    let mut max_x_diff: f64 = 0.0;
    let mut max_y_diff: f64 = 0.0;

    for (i, xf) in xfoil_upper.iter().take(10).enumerate() {
        // Find closest yfoil point
        if let Some((yi, _)) = x_upper
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (*a - xf.1).abs().partial_cmp(&(*b - xf.1).abs()).unwrap())
        {
            let x_diff = (x_upper[yi] - xf.1).abs();
            let y_diff = (airfoil.y[stag + yi] - xf.6).abs();

            max_x_diff = max_x_diff.max(x_diff);
            max_y_diff = max_y_diff.max(y_diff);

            println!(
                "{:>6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}",
                i, xf.1, x_upper[yi], xf.6, airfoil.y[stag + yi]
            );
        }
    }
    println!("...");
    println!("Max x difference: {:.2e}", max_x_diff);
    println!("Max y difference: {:.2e}", max_y_diff);

    // Compare edge velocity
    println!("\n=== Edge Velocity Comparison ===");
    println!(
        "{:>8} {:>12} {:>12} {:>10}",
        "x/c", "Ue_xfoil", "Ue_yfoil", "diff %"
    );
    println!("{:-<50}", "");

    let check_xc = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90];

    for &target_xc in &check_xc {
        if let Some(xf) = xfoil_upper
            .iter()
            .min_by(|a, b| (a.1 - target_xc).abs().partial_cmp(&(b.1 - target_xc).abs()).unwrap())
        {
            if let Some((yi, _)) = x_upper
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| (*a - target_xc).abs().partial_cmp(&(*b - target_xc).abs()).unwrap())
            {
                if yi < ue_upper.len() {
                    let diff_pct = 100.0 * (ue_upper[yi] - xf.2) / xf.2;
                    println!(
                        "{:>8.2} {:>12.6} {:>12.6} {:>10.3}",
                        target_xc, xf.2, ue_upper[yi], diff_pct
                    );
                }
            }
        }
    }

    // Compare BL solution
    println!("\n=== Boundary Layer Comparison ===");
    println!(
        "{:>8} {:>12} {:>12} {:>10} {:>10} {:>10} {:>10}",
        "x/c", "θ_xfoil", "θ_yfoil", "θ ratio", "H_xfoil", "H_yfoil", "H diff"
    );
    println!("{:-<90}", "");

    for &target_xc in &check_xc {
        if let Some(xf) = xfoil_upper
            .iter()
            .min_by(|a, b| (a.1 - target_xc).abs().partial_cmp(&(b.1 - target_xc).abs()).unwrap())
        {
            if let Some((yi, _)) = x_upper
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| (*a - target_xc).abs().partial_cmp(&(*b - target_xc).abs()).unwrap())
            {
                if yi < bl_upper.len() {
                    let theta_ratio = bl_upper[yi].theta / xf.3;
                    let h_diff = bl_upper[yi].hk - xf.5;
                    println!(
                        "{:>8.2} {:>12.2e} {:>12.2e} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
                        target_xc, xf.3, bl_upper[yi].theta, theta_ratio, xf.5, bl_upper[yi].hk, h_diff
                    );
                }
            }
        }
    }

    // Focus on the divergence region
    println!("\n=== Detailed H Comparison (x/c = 0.50 to 0.75) ===");
    println!(
        "{:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "x/c", "H_xfoil", "H_yfoil", "H_diff", "Cf_xfoil", "Cf_yfoil"
    );
    println!("{:-<70}", "");

    for (yi, &xc) in x_upper.iter().enumerate() {
        if xc >= 0.50 && xc <= 0.75 && yi < bl_upper.len() {
            if let Some(xf) = xfoil_upper
                .iter()
                .min_by(|a, b| (a.1 - xc).abs().partial_cmp(&(b.1 - xc).abs()).unwrap())
            {
                if (xf.1 - xc).abs() < 0.02 {
                    let h_diff = bl_upper[yi].hk - xf.5;
                    println!(
                        "{:>8.4} {:>10.4} {:>10.4} {:>10.4} {:>10.6} {:>10.6}",
                        xc, xf.5, bl_upper[yi].hk, h_diff, xf.4, bl_upper[yi].cf
                    );
                }
            }
        }
    }

    // Summary
    println!("\n=== Summary ===");
    println!("yfoil transition: x/c = {:.4}", trans_xc);
    println!("XFOIL transition: x/c ≈ 0.687 (from XFOIL output)");

    let converged_count = bl_upper.iter().filter(|r| r.converged).count();
    let total_count = bl_upper.len();
    println!(
        "Converged stations: {}/{} ({:.1}%)",
        converged_count,
        total_count,
        100.0 * converged_count as f64 / total_count as f64
    );
}
