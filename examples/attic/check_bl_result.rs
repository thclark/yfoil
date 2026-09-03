//! Debug BL result cf values

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    let result = solve_viscous(&airfoil, 0.0, &cond, &config);

    println!("CL = {}", result.cl);
    println!("CD = {}", result.cd);
    println!("CDf = {}", result.cdf);
    println!("CDp = {}", result.cdp);
    println!("Converged: {}", result.converged);
    println!("Iterations: {}", result.iterations);
    println!("Residual: {}", result.residual);
    println!();

    println!("Upper surface BL stations: {}", result.bl.upper.len());
    println!("Lower surface BL stations: {}", result.bl.lower.len());
    println!();

    println!("Upper surface cf values (first 10):");
    for (i, r) in result.bl.upper.iter().enumerate().take(10) {
        println!("  [{}] theta={:.6e}, dstar={:.6e}, h={:.3}, hk={:.3}, cf={:.6e}, ue={:.4}",
                 i, r.theta, r.dstar, r.h, r.hk, r.cf, r.ue);
    }
    println!();

    println!("Lower surface cf values (first 10):");
    for (i, r) in result.bl.lower.iter().enumerate().take(10) {
        println!("  [{}] theta={:.6e}, dstar={:.6e}, h={:.3}, hk={:.3}, cf={:.6e}, ue={:.4}",
                 i, r.theta, r.dstar, r.h, r.hk, r.cf, r.ue);
    }
    println!();

    // Check for negative or NaN cf
    let upper_cf_issues: Vec<_> = result.bl.upper.iter().enumerate()
        .filter(|(_, r)| r.cf.is_nan() || r.cf < 0.0)
        .collect();
    let lower_cf_issues: Vec<_> = result.bl.lower.iter().enumerate()
        .filter(|(_, r)| r.cf.is_nan() || r.cf < 0.0)
        .collect();

    if !upper_cf_issues.is_empty() {
        println!("ISSUES in upper surface cf:");
        for (i, r) in upper_cf_issues {
            println!("  [{}] cf={:.6e}, hk={:.3}, theta={:.6e}", i, r.cf, r.hk, r.theta);
        }
    }
    if !lower_cf_issues.is_empty() {
        println!("ISSUES in lower surface cf:");
        for (i, r) in lower_cf_issues {
            println!("  [{}] cf={:.6e}, hk={:.3}, theta={:.6e}", i, r.cf, r.hk, r.theta);
        }
    }

    // Calculate CDf manually
    let cdf_upper: f64 = (1..result.bl.upper.len())
        .map(|i| {
            let dx = result.bl.x_upper[i] - result.bl.x_upper[i - 1];
            let cf_avg = 0.5 * (result.bl.upper[i - 1].cf + result.bl.upper[i].cf);
            cf_avg * dx
        })
        .sum();
    let cdf_lower: f64 = (1..result.bl.lower.len())
        .map(|i| {
            let dx = result.bl.x_lower[i] - result.bl.x_lower[i - 1];
            let cf_avg = 0.5 * (result.bl.lower[i - 1].cf + result.bl.lower[i].cf);
            cf_avg * dx
        })
        .sum();

    println!("\nManual CDf calculation:");
    println!("  Upper: {:.6e}", cdf_upper);
    println!("  Lower: {:.6e}", cdf_lower);
    println!("  Total: {:.6e}", cdf_upper + cdf_lower);
}
