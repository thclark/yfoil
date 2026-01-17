//! Debug theta growth rate in the BL
//!
//! Compare theta growth vs theoretical expectations

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Debug Theta Growth ===\n");

    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let mut config = ViscalConfig::default();
    config.tol_cl = 5e-4;

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("Converged: {}\n", result.converged);

    // Analyze theta growth rate
    println!("=== Upper Surface Theta Growth ===");
    println!("{:>4} {:>8} {:>10} {:>10} {:>10} {:>8} {:>8}",
             "i", "x/c", "theta", "d_theta", "growth%", "H", "regime");
    println!("{}", "-".repeat(78));

    let mut laminar = true;
    for i in 1..result.bl.upper.len() {
        let prev = &result.bl.upper[i - 1];
        let curr = &result.bl.upper[i];
        let x = result.bl.x_upper.get(i).copied().unwrap_or(0.0);

        let d_theta = curr.theta - prev.theta;
        let growth_pct = if prev.theta > 0.0 {
            100.0 * d_theta / prev.theta
        } else {
            0.0
        };

        // Check for transition (n_amp crossing ncrit=9)
        if laminar && curr.n_amp >= 9.0 {
            laminar = false;
            println!("  ** TRANSITION **");
        }

        let regime = if laminar { "LAM" } else { "TURB" };

        // Flag anomalies
        let flag = if d_theta < 0.0 {
            " << DECREASE!"
        } else if growth_pct > 100.0 {
            " << LARGE JUMP!"
        } else {
            ""
        };

        println!("{:4} {:8.4} {:10.6} {:10.6} {:10.1} {:8.3} {:>8}{}",
                 i, x, curr.theta, d_theta, growth_pct, curr.h, regime, flag);
    }

    // Check trailing edge values vs expected
    println!("\n=== Trailing Edge Analysis ===");
    let upper_te = result.bl.upper.last().unwrap();
    let lower_te = result.bl.lower.last().unwrap();

    println!("Upper TE: θ = {:.6}, H = {:.3}", upper_te.theta, upper_te.h);
    println!("Lower TE: θ = {:.6}, H = {:.3}", lower_te.theta, lower_te.h);
    println!("\nExpected (XFOIL): θ_upper ≈ 0.0015, θ_lower ≈ 0.0015");
    println!("Ratio: θ_actual/θ_expected = {:.2}", upper_te.theta / 0.0015);

    // Sum up where the excess theta comes from
    println!("\n=== Cumulative Theta Growth ===");
    let mut theta_lam_end = 0.0;
    let mut theta_trans = 0.0;
    let mut idx_trans = 0;
    for i in 0..result.bl.upper.len() {
        if result.bl.upper[i].n_amp >= 9.0 && theta_trans == 0.0 {
            theta_trans = result.bl.upper[i].theta;
            idx_trans = i;
        }
        if result.bl.upper[i].n_amp < 9.0 {
            theta_lam_end = result.bl.upper[i].theta;
        }
    }

    let theta_turb_growth = upper_te.theta - theta_trans;

    println!("θ at end of laminar (before trans): {:.6}", theta_lam_end);
    println!("θ at transition point (station {}): {:.6}", idx_trans, theta_trans);
    println!("θ at trailing edge: {:.6}", upper_te.theta);
    println!("θ growth in turbulent region: {:.6} ({:.1}% of total)",
             theta_turb_growth, 100.0 * theta_turb_growth / upper_te.theta);

    // Estimate expected values based on Blasius
    // At x=0.5 (roughly transition), theta_blasius ≈ 0.00033
    // From x=0.5 to x=1.0 in turbulent: theta grows roughly 5x
    // So expected theta_TE ≈ 0.0015
    println!("\n=== Theoretical Estimates ===");
    let x_trans = result.bl.x_upper.get(idx_trans).copied().unwrap_or(0.5);
    let theta_blasius_trans = 0.664 * x_trans / (1e6 * x_trans).sqrt();
    println!("Blasius θ at x={:.2}: {:.6}", x_trans, theta_blasius_trans);
    println!("Actual θ at transition: {:.6}", theta_trans);
    println!("Ratio actual/Blasius at transition: {:.2}", theta_trans / theta_blasius_trans);
}
