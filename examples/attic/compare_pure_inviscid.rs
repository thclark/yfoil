//! Compare pure inviscid solutions: yfoil vs XFOIL
//!
//! Uses XFOIL's CPWR output (pure inviscid) to compare against yfoil's panel method.

use std::fs::File;
use std::io::{BufRead, BufReader};
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, find_stagnation_point};

fn main() {
    println!("=== Pure Inviscid Solution Comparison ===\n");

    // Load XFOIL-paneled coordinates
    let (_, geom) = read_dat_file("/tmp/xfoil_naca0012_paneled.dat")
        .expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, _, ue_yfoil) = extract_upper_surface(&airfoil, &vel, stag);

    println!("yfoil panel method:");
    println!("  Panels: {}", airfoil.n);
    println!("  Upper surface points: {}", x_upper.len());

    // Load XFOIL's pure inviscid Cp distribution
    let file = File::open("/tmp/xfoil_cp_inviscid.dat")
        .expect("Run XFOIL in inviscid mode first");
    let reader = BufReader::new(file);

    // Parse XFOIL Cp data - format: x, Cp
    let mut xfoil_data: Vec<(f64, f64, f64)> = Vec::new(); // x, Cp, Ue

    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let x: f64 = parts[0].parse().unwrap_or(0.0);
            let cp: f64 = parts[1].parse().unwrap_or(0.0);
            // Ue/U∞ = sqrt(1 - Cp)
            let ue = if cp < 1.0 { (1.0 - cp).sqrt() } else { 0.0 };
            xfoil_data.push((x, cp, ue));
        }
    }

    println!("\nXFOIL pure inviscid:");
    println!("  Points: {}", xfoil_data.len());

    // XFOIL Cp file is in airfoil order (TE -> upper -> LE -> lower -> TE)
    // Find LE (minimum x) to split upper/lower
    let le_idx = xfoil_data.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.0.partial_cmp(&b.0).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Upper surface is from index 0 to le_idx (TE to LE, decreasing x)
    let xfoil_upper: Vec<_> = xfoil_data[0..=le_idx].to_vec();

    println!("  Upper surface points: {} (to LE at index {})", xfoil_upper.len(), le_idx);

    // Compare Cp and Ue distributions
    println!("\n=== Cp Comparison (Upper Surface) ===");
    println!("{:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
             "x/c", "Cp_xfoil", "Cp_yfoil", "Cp_diff", "Ue_xfoil", "Ue_yfoil");
    println!("{:-<70}", "");

    let mut sum_cp_diff_sq = 0.0;
    let mut sum_ue_diff_sq = 0.0;
    let mut count = 0;

    // yfoil goes from LE to TE (increasing x), XFOIL goes TE to LE (decreasing x)
    for (yi, &xc) in x_upper.iter().enumerate() {
        if xc >= 0.02 && xc <= 0.98 && yi < ue_yfoil.len() {
            // Find closest XFOIL point
            if let Some(xf) = xfoil_upper
                .iter()
                .min_by(|a, b| (a.0 - xc).abs().partial_cmp(&(b.0 - xc).abs()).unwrap())
            {
                if (xf.0 - xc).abs() < 0.02 {
                    // Compute yfoil Cp from Ue: Cp = 1 - Ue²
                    let cp_yfoil = 1.0 - ue_yfoil[yi].powi(2);
                    let cp_diff = cp_yfoil - xf.1;
                    let ue_diff = ue_yfoil[yi] - xf.2;

                    sum_cp_diff_sq += cp_diff.powi(2);
                    sum_ue_diff_sq += ue_diff.powi(2);
                    count += 1;

                    // Print every few points
                    if (xc * 10.0) as i32 % 1 == 0 {
                        println!("{:>8.3} {:>10.5} {:>10.5} {:>10.5} {:>10.5} {:>10.5}",
                                 xc, xf.1, cp_yfoil, cp_diff, xf.2, ue_yfoil[yi]);
                    }
                }
            }
        }
    }

    let rms_cp = (sum_cp_diff_sq / count as f64).sqrt();
    let rms_ue = (sum_ue_diff_sq / count as f64).sqrt();

    println!("\n=== Statistics ===");
    println!("RMS Cp difference: {:.6}", rms_cp);
    println!("RMS Ue difference: {:.6}", rms_ue);

    // Regional analysis
    println!("\n=== Regional Ue Difference (% of XFOIL) ===");
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
                    if (xf.0 - xc).abs() < 0.02 && xf.2.abs() > 0.01 {
                        sum += 100.0 * (ue_yfoil[yi] - xf.2) / xf.2;
                        n += 1;
                    }
                }
            }
        }

        if n > 0 {
            println!("  {:12}: {:+.4}%", name, sum / n as f64);
        }
    }

    // Peak velocity comparison
    println!("\n=== Peak Velocity ===");
    let max_ue_yfoil = ue_yfoil.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let max_ue_xfoil = xfoil_upper.iter().map(|p| p.2).fold(f64::NEG_INFINITY, f64::max);
    println!("yfoil max Ue: {:.6}", max_ue_yfoil);
    println!("XFOIL max Ue: {:.6}", max_ue_xfoil);
    println!("Difference: {:.4}%", 100.0 * (max_ue_yfoil - max_ue_xfoil) / max_ue_xfoil);

    // TE velocity comparison
    println!("\n=== Trailing Edge ===");
    let te_ue_yfoil = ue_yfoil.last().unwrap();
    let te_ue_xfoil = xfoil_upper.first().unwrap().2; // XFOIL starts at TE
    println!("yfoil TE Ue: {:.6}", te_ue_yfoil);
    println!("XFOIL TE Ue: {:.6}", te_ue_xfoil);
    println!("Difference: {:.4}%", 100.0 * (te_ue_yfoil - te_ue_xfoil) / te_ue_xfoil);
}
