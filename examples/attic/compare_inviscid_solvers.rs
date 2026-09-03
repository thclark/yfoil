//! Compare inviscid solvers: yfoil panel method vs XFOIL
//!
//! This investigates why yfoil produces ~0.5% different edge velocities

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, find_stagnation_point};

fn main() {
    println!("=== Inviscid Solver Comparison ===\n");

    // Load XFOIL-paneled coordinates (both use same geometry)
    let (_, geom) = read_dat_file("/tmp/xfoil_naca0012_paneled.dat")
        .expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_yfoil) = extract_upper_surface(&airfoil, &vel, stag);

    println!("yfoil panel method:");
    println!("  Panels: {}", airfoil.n);
    println!("  Stagnation point: index {}, x = {:.6}", stag, airfoil.x[stag]);
    println!("  Upper surface points: {}", x_upper.len());

    // Load XFOIL's Ue distribution
    let file = File::open("/tmp/xfoil_bl_paneled.dat").expect("XFOIL BL dump not found");
    let reader = BufReader::new(file);

    let mut xfoil_data: Vec<(f64, f64, f64)> = Vec::new(); // x, ue, y

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 4 {
            let x: f64 = parts[1].parse().unwrap_or(0.0);
            let y: f64 = parts[2].parse().unwrap_or(0.0);
            let ue: f64 = parts[3].parse().unwrap_or(0.0);
            if y >= 0.0 {
                xfoil_data.push((x, ue, y));
            }
        }
    }

    xfoil_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Find stagnation point in XFOIL data
    let stag_idx_xfoil = xfoil_data
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    let xfoil_upper: Vec<_> = xfoil_data[stag_idx_xfoil..].to_vec();

    println!("\nXFOIL inviscid:");
    println!("  Upper surface points: {}", xfoil_upper.len());

    // Compare Ue distributions
    println!("\n=== Edge Velocity Comparison ===");
    println!("{:>8} {:>12} {:>12} {:>10} {:>10}",
             "x/c", "Ue_xfoil", "Ue_yfoil", "diff", "diff %");
    println!("{:-<60}", "");

    let mut max_diff_pct: f64 = 0.0;
    let mut max_diff_x: f64 = 0.0;

    for (yi, &xc) in x_upper.iter().enumerate() {
        if xc >= 0.02 && xc <= 0.98 && yi < ue_yfoil.len() {
            // Find closest XFOIL point
            if let Some(xf) = xfoil_upper
                .iter()
                .min_by(|a, b| (a.0 - xc).abs().partial_cmp(&(b.0 - xc).abs()).unwrap())
            {
                if (xf.0 - xc).abs() < 0.02 {
                    let diff = ue_yfoil[yi] - xf.1;
                    let diff_pct = 100.0 * diff / xf.1;

                    if diff_pct.abs() > max_diff_pct.abs() {
                        max_diff_pct = diff_pct;
                        max_diff_x = xc;
                    }

                    // Print every 10th point
                    if (xc * 10.0) as i32 % 1 == 0 && xc > 0.05 {
                        println!("{:>8.3} {:>12.6} {:>12.6} {:>10.6} {:>10.3}",
                                 xc, xf.1, ue_yfoil[yi], diff, diff_pct);
                    }
                }
            }
        }
    }

    println!("\nMaximum difference: {:.3}% at x/c = {:.3}", max_diff_pct, max_diff_x);

    // Analyze the pattern
    println!("\n=== Difference Pattern Analysis ===");
    println!("Checking systematic bias in Ue:");

    let mut sum_diff_pct = 0.0;
    let mut count = 0;

    for (yi, &xc) in x_upper.iter().enumerate() {
        if xc >= 0.1 && xc <= 0.9 && yi < ue_yfoil.len() {
            if let Some(xf) = xfoil_upper
                .iter()
                .min_by(|a, b| (a.0 - xc).abs().partial_cmp(&(b.0 - xc).abs()).unwrap())
            {
                if (xf.0 - xc).abs() < 0.02 {
                    let diff_pct = 100.0 * (ue_yfoil[yi] - xf.1) / xf.1;
                    sum_diff_pct += diff_pct;
                    count += 1;
                }
            }
        }
    }

    let avg_diff_pct = sum_diff_pct / count as f64;
    println!("Average Ue difference (x/c = 0.1 to 0.9): {:.3}%", avg_diff_pct);

    // Check if it's a systematic offset or varying
    println!("\n=== Regional Analysis ===");

    for region in [(0.0, 0.2, "LE region"), (0.2, 0.5, "Forward"), (0.5, 0.8, "Aft"), (0.8, 1.0, "TE region")] {
        let (x_min, x_max, name) = region;
        let mut sum = 0.0;
        let mut n = 0;

        for (yi, &xc) in x_upper.iter().enumerate() {
            if xc >= x_min && xc <= x_max && yi < ue_yfoil.len() {
                if let Some(xf) = xfoil_upper
                    .iter()
                    .min_by(|a, b| (a.0 - xc).abs().partial_cmp(&(b.0 - xc).abs()).unwrap())
                {
                    if (xf.0 - xc).abs() < 0.02 {
                        sum += 100.0 * (ue_yfoil[yi] - xf.1) / xf.1;
                        n += 1;
                    }
                }
            }
        }

        if n > 0 {
            println!("  {:12}: avg diff = {:+.3}%", name, sum / n as f64);
        }
    }

    // Check Cp at stagnation point and TE
    println!("\n=== Key Points ===");
    println!("Stagnation point:");
    println!("  yfoil: x = {:.6}, Ue = {:.6}", x_upper[0], ue_yfoil[0]);
    println!("  XFOIL: x = {:.6}, Ue = {:.6}", xfoil_upper[0].0, xfoil_upper[0].1);

    println!("\nTrailing edge:");
    let last_yfoil = x_upper.len() - 1;
    println!("  yfoil: x = {:.6}, Ue = {:.6}", x_upper[last_yfoil], ue_yfoil[last_yfoil]);
    let last_xfoil = xfoil_upper.len() - 1;
    println!("  XFOIL: x = {:.6}, Ue = {:.6}", xfoil_upper[last_xfoil].0, xfoil_upper[last_xfoil].1);

    // Check maximum Ue (should be near LE)
    let max_ue_yfoil = ue_yfoil.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let max_ue_idx_yfoil = ue_yfoil.iter().position(|&u| u == max_ue_yfoil).unwrap();
    let max_ue_xfoil = xfoil_upper.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);

    println!("\nMaximum Ue:");
    println!("  yfoil: {:.6} at x/c = {:.4}", max_ue_yfoil, x_upper[max_ue_idx_yfoil]);
    println!("  XFOIL: {:.6}", max_ue_xfoil);
    println!("  Ratio: {:.4}", max_ue_yfoil / max_ue_xfoil);
}
