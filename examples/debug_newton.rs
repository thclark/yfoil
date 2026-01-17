//! Debug Newton iteration at specific stations
//!
//! Trace what happens during turbulent BL marching

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Debug Newton Iteration ===\n");

    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let mut config = ViscalConfig::default();
    config.tol_cl = 5e-4;

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("Converged: {}\n", result.converged);

    // Look at the momentum equation balance at key stations
    println!("=== Momentum Balance Check ===");
    println!("At each station, momentum eq should be satisfied:");
    println!("dθ/ds + (H+2)*θ/Ue * dUe/ds = Cf/2\n");

    let upper = &result.bl.upper;
    let s = &result.bl.s_upper;
    let x = &result.bl.x_upper;
    let nu = 1.0 / 1_000_000.0;

    println!("{:>4} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
             "i", "x/c", "dtheta/ds", "shape_term", "sum", "Cf/2", "residual");
    println!("{}", "-".repeat(78));

    for i in 26..upper.len().min(40) {
        let prev = &upper[i - 1];
        let curr = &upper[i];
        let ds = s[i] - s[i - 1];

        let dtheta_ds = (curr.theta - prev.theta) / ds;
        let ue_avg = 0.5 * (prev.ue + curr.ue);
        let h_avg = 0.5 * (prev.h + curr.h);
        let theta_avg = 0.5 * (prev.theta + curr.theta);
        let due_ds = (curr.ue - prev.ue) / ds;

        let shape_term = (h_avg + 2.0) * theta_avg / ue_avg * due_ds;
        let lhs = dtheta_ds + shape_term;
        let cf_avg = 0.5 * (prev.cf + curr.cf);
        let rhs = cf_avg / 2.0;
        let residual = lhs - rhs;

        let flag = if residual.abs() > 0.001 { " << LARGE" } else { "" };

        println!("{:4} {:8.4} {:10.4e} {:10.4e} {:10.4e} {:10.4e} {:10.4e}{}",
                 i, x[i], dtheta_ds, shape_term, lhs, rhs, residual, flag);
    }

    // Check derivative consistency
    println!("\n=== Check: Is θ decrease possible? ===");
    println!("For θ to decrease: dθ/ds < 0");
    println!("From mom eq: dθ/ds = Cf/2 - (H+2)*θ/Ue * dUe/ds");
    println!("For adverse pressure gradient, dUe/ds < 0, so shape_term > 0");
    println!("dθ/ds < 0 requires Cf/2 < shape_term\n");

    for i in 26..upper.len().min(40) {
        let prev = &upper[i - 1];
        let curr = &upper[i];
        let ds = s[i] - s[i - 1];

        let ue_avg = 0.5 * (prev.ue + curr.ue);
        let h_avg = 0.5 * (prev.h + curr.h);
        let theta_avg = 0.5 * (prev.theta + curr.theta);
        let due_ds = (curr.ue - prev.ue) / ds;

        let shape_term = (h_avg + 2.0) * theta_avg / ue_avg * due_ds;
        let cf_avg = 0.5 * (prev.cf + curr.cf);

        let could_decrease = cf_avg / 2.0 < shape_term;
        let did_decrease = curr.theta < prev.theta;

        if could_decrease || did_decrease {
            println!("Station {}: Cf/2={:.4e}, shape={:.4e}, could_dec={}, did_dec={}",
                     i, cf_avg / 2.0, shape_term, could_decrease, did_decrease);
        }
    }

    // Summary: where is the extra θ coming from?
    println!("\n=== Θ Growth Summary ===");
    let mut total_growth = 0.0;
    let mut pos_growth = 0.0;
    let mut neg_growth = 0.0;
    for i in 23..upper.len() {
        let prev = &upper[i - 1];
        let curr = &upper[i];
        let delta = curr.theta - prev.theta;
        total_growth += delta;
        if delta > 0.0 {
            pos_growth += delta;
        } else {
            neg_growth += delta;
        }
    }
    println!("Turbulent region (stations 23 to TE):");
    println!("  Total θ growth: {:.6}", total_growth);
    println!("  Positive growth: {:.6}", pos_growth);
    println!("  Negative growth: {:.6}", neg_growth);
    println!("  Net = total: {:.6}", pos_growth + neg_growth);
}
