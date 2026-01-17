//! Compare BL solvers using XFOIL's edge velocity distribution
//!
//! This eliminates the inviscid solver difference by using XFOIL's
//! Ue(s) directly as input to yfoil's BL solver.

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::bl::{march_newton, FlowConditions, NewtonConfig};

fn main() {
    println!("=== BL Solver Comparison with Identical Ue ===\n");

    // Load XFOIL's BL dump which contains Ue(s) distribution
    let file = File::open("/tmp/xfoil_bl_paneled.dat")
        .expect("Run XFOIL first to generate BL dump");
    let reader = BufReader::new(file);

    // Parse XFOIL BL data
    // Format: s, x, y, Ue/Vinf, Dstar, Theta, Cf, H, ...
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

            // Upper surface only (y >= 0)
            if y >= 0.0 {
                xfoil_data.push((s, x, ue, theta, cf, h, y));
            }
        }
    }

    // Sort by x to get from LE to TE
    xfoil_data.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Find stagnation point (min |Ue|)
    let stag_idx = xfoil_data
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.2.abs().partial_cmp(&b.2.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Extract upper surface from stagnation to TE
    let xfoil_upper: Vec<_> = xfoil_data[stag_idx..].to_vec();

    println!("XFOIL upper surface: {} stations", xfoil_upper.len());
    println!("  x range: {:.4} to {:.4}", xfoil_upper.first().unwrap().1, xfoil_upper.last().unwrap().1);

    // Build s and Ue arrays from XFOIL data
    // Note: XFOIL's s is from the TE, we need to adjust
    let s_te = xfoil_upper.last().unwrap().0;
    let s_stag = xfoil_upper.first().unwrap().0;

    // Create monotonically increasing s from stagnation point
    let mut s_dist: Vec<f64> = Vec::with_capacity(xfoil_upper.len());
    let mut ue_dist: Vec<f64> = Vec::with_capacity(xfoil_upper.len());
    let mut x_dist: Vec<f64> = Vec::with_capacity(xfoil_upper.len());

    for (i, pt) in xfoil_upper.iter().enumerate() {
        // s from stagnation point
        let s_from_stag = if i == 0 { 0.0 } else { pt.0 - s_stag };
        s_dist.push(s_from_stag.abs());
        ue_dist.push(pt.2);
        x_dist.push(pt.1);
    }

    // Ensure s is monotonically increasing
    for i in 1..s_dist.len() {
        if s_dist[i] <= s_dist[i-1] {
            s_dist[i] = s_dist[i-1] + 0.001;
        }
    }

    println!("  s range: {:.4} to {:.4}", s_dist.first().unwrap(), s_dist.last().unwrap());
    println!("  Ue range: {:.4} to {:.4}\n", ue_dist.iter().cloned().fold(f64::INFINITY, f64::min),
             ue_dist.iter().cloned().fold(f64::NEG_INFINITY, f64::max));

    // Run yfoil BL solver with XFOIL's Ue distribution
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_yfoil = march_newton(&ue_dist, &s_dist, &cond, &config);

    // Find transition in yfoil
    let trans_idx_yfoil = bl_yfoil.iter().position(|r| r.n_amp >= cond.ncrit);
    let trans_x_yfoil = trans_idx_yfoil.map(|i| x_dist[i]).unwrap_or(f64::INFINITY);

    println!("yfoil transition (with XFOIL's Ue): x/c = {:.4}", trans_x_yfoil);
    println!("XFOIL transition: x/c ≈ 0.687\n");

    // Compare BL parameters station by station
    println!("=== Station-by-Station Comparison ===");
    println!("(Using identical Ue distribution from XFOIL)");
    println!();
    println!("{:>8} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10}",
             "x/c", "H_xfoil", "H_yfoil", "H_diff", "θ_xfoil", "θ_yfoil", "θ_ratio");
    println!("{:-<90}", "");

    for (i, pt) in xfoil_upper.iter().enumerate() {
        let x = pt.1;
        if i < bl_yfoil.len() && x >= 0.05 && x <= 0.75 && (x * 20.0) as i32 % 2 == 0 {
            let h_xfoil = pt.5;
            let h_yfoil = bl_yfoil[i].hk;
            let h_diff = h_yfoil - h_xfoil;

            let theta_xfoil = pt.3;
            let theta_yfoil = bl_yfoil[i].theta;
            let theta_ratio = theta_yfoil / theta_xfoil;

            println!("{:>8.4} {:>10.4} {:>10.4} {:>10.4} {:>12.2e} {:>12.2e} {:>10.4}",
                     x, h_xfoil, h_yfoil, h_diff, theta_xfoil, theta_yfoil, theta_ratio);
        }
    }

    // Detailed comparison near transition region
    println!("\n=== Detailed Comparison (x/c = 0.50 to 0.70) ===");
    println!("{:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
             "x/c", "H_xfoil", "H_yfoil", "H_diff", "n_yfoil", "Ue", "conv");
    println!("{:-<80}", "");

    for (i, pt) in xfoil_upper.iter().enumerate() {
        let x = pt.1;
        if i < bl_yfoil.len() && x >= 0.50 && x <= 0.70 {
            let h_xfoil = pt.5;
            let h_yfoil = bl_yfoil[i].hk;
            let h_diff = h_yfoil - h_xfoil;
            let conv = if bl_yfoil[i].converged { "yes" } else { "NO" };

            println!("{:>8.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10}",
                     x, h_xfoil, h_yfoil, h_diff, bl_yfoil[i].n_amp, ue_dist[i], conv);
        }
    }

    // Summary
    println!("\n=== Summary ===");
    println!("With XFOIL's exact Ue distribution:");
    println!("  yfoil transition: x/c = {:.4}", trans_x_yfoil);
    println!("  XFOIL transition: x/c ≈ 0.687");

    if trans_x_yfoil < 0.687 {
        let diff_pct = (0.687 - trans_x_yfoil) / 0.687 * 100.0;
        println!("  yfoil transitions {:.1}% earlier", diff_pct);
        println!("\n  This isolates the issue to the BL marching algorithm,");
        println!("  not the inviscid solver.");
    } else {
        println!("  Transition matches! The difference was in the inviscid solver.");
    }

    // Check convergence
    let converged_count = bl_yfoil.iter().filter(|r| r.converged).count();
    println!("\nConverged stations: {}/{}", converged_count, bl_yfoil.len());
}
