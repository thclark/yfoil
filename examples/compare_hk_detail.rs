use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    // Run yfoil with same conditions as XFOIL
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    // Load XFOIL dump data
    let file = File::open("/tmp/xfoil_bl.dat").expect("Failed to open XFOIL dump");
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

    xfoil_data.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    xfoil_data.reverse();

    let stag_idx = xfoil_data.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.2.abs().partial_cmp(&b.2.abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);

    let xfoil_upper: Vec<_> = xfoil_data[stag_idx..].to_vec();

    println!("=== Detailed Hk Comparison (x/c = 0.40 to 0.70) ===");
    println!();
    println!("{:>8} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10} {:>10} {:>10}",
        "x/c", "H_xfoil", "H_yfoil", "H_diff", "θ_xfoil", "θ_yfoil", "Cf_xfoil", "Cf_yfoil", "converged");
    println!("{:-<120}", "");

    // Find all yfoil stations in range 0.40 to 0.70
    for (yi, &xc) in x_upper.iter().enumerate() {
        if xc >= 0.40 && xc <= 0.70 && yi < bl_upper.len() {
            // Find closest XFOIL point
            if let Some(xf) = xfoil_upper.iter()
                .min_by(|a, b| (a.1 - xc).abs().partial_cmp(&(b.1 - xc).abs()).unwrap())
            {
                if (xf.1 - xc).abs() < 0.02 {
                    let h_diff = bl_upper[yi].hk - xf.5;
                    let conv = if bl_upper[yi].converged { "yes" } else { "NO" };

                    println!("{:>8.4} {:>10.4} {:>10.4} {:>10.4} {:>12.2e} {:>12.2e} {:>10.6} {:>10.6} {:>10}",
                        xc, xf.5, bl_upper[yi].hk, h_diff,
                        xf.3, bl_upper[yi].theta,
                        xf.4, bl_upper[yi].cf, conv);
                }
            }
        }
    }

    // Check where inverse mode kicks in for yfoil
    println!();
    println!("=== yfoil Inverse Mode Check ===");
    println!("Hmax (laminar) = {}", config.hlmax);
    println!();

    let mut first_inverse = None;
    for (i, r) in bl_upper.iter().enumerate() {
        if !r.converged && first_inverse.is_none() {
            first_inverse = Some(i);
        }
        if i > 0 && x_upper[i] >= 0.50 && x_upper[i] <= 0.70 {
            let conv = if r.converged { " " } else { "*" };
            println!("  Stn {:3}: x/c = {:.4}, Hk = {:.4}, n_amp = {:.4} {}",
                i, x_upper[i], r.hk, r.n_amp, conv);
        }
    }

    if let Some(idx) = first_inverse {
        println!();
        println!("First non-converged station: {} at x/c = {:.4}, Hk = {:.4}",
            idx, x_upper[idx], bl_upper[idx].hk);
    }

    // Look at pressure gradient in this region
    println!();
    println!("=== Pressure Gradient (dUe/ds) in Region ===");
    println!("{:>8} {:>12} {:>12} {:>12}", "x/c", "Ue", "dUe/ds", "type");
    println!("{:-<50}", "");

    for i in 1..x_upper.len() {
        if x_upper[i] >= 0.40 && x_upper[i] <= 0.70 {
            let ds = s_upper[i] - s_upper[i-1];
            let due_ds = (ue_upper[i] - ue_upper[i-1]) / ds;
            let grad_type = if due_ds > 0.01 { "favorable" }
                else if due_ds < -0.01 { "ADVERSE" }
                else { "~zero" };
            println!("{:>8.4} {:>12.6} {:>12.6} {:>12}", x_upper[i], ue_upper[i], due_ds, grad_type);
        }
    }
}
