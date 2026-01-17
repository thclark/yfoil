use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::{FlowConditions, NewtonConfig, march_newton, cf_lam, hkin};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    println!("=== Cf Analysis ===");
    println!();
    println!("{:>4} {:>8} {:>10} {:>10} {:>8} {:>10} {:>10} {:>10}",
        "Stn", "x/c", "θ", "Hk", "Rθ", "Cf_stored", "Cf_calc", "Diff(%)");
    println!("{:-<90}", "");
    
    for i in 0..20.min(bl_upper.len()) {
        let theta = bl_upper[i].theta;
        let ue = ue_upper[i];
        let hk = bl_upper[i].hk;
        let rt = ue * theta / cond.nu;
        let cf_stored = bl_upper[i].cf;
        
        // Recalculate Cf from the closure relation
        let cf_calc = cf_lam(hk, rt, cond.msq);
        
        let diff_pct = (cf_stored - cf_calc.val) / cf_calc.val.abs().max(1e-10) * 100.0;
        
        println!("{:>4} {:>8.4} {:>10.2e} {:>8.4} {:>10.1} {:>10.4} {:>10.4} {:>10.2}",
            i, x_upper[i], theta, hk, rt, cf_stored, cf_calc.val, diff_pct);
    }
    
    // Print momentum balance for a few stations
    println!();
    println!("=== Momentum Balance Check ===");
    println!("dθ/ds + (H+2-M²) * θ/Ue * dUe/ds = Cf/2");
    println!();
    println!("{:>4} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Stn", "x/c", "dθ/ds", "(H+2)θ/Ue*dUe/ds", "Cf/2", "LHS", "RHS-LHS");
    println!("{:-<90}", "");
    
    for i in 1..15.min(bl_upper.len()) {
        let ds = s_upper[i] - s_upper[i-1];
        let dtheta_ds = (bl_upper[i].theta - bl_upper[i-1].theta) / ds;
        
        let theta_avg = 0.5 * (bl_upper[i].theta + bl_upper[i-1].theta);
        let h_avg = 0.5 * (bl_upper[i].h + bl_upper[i-1].h);
        let ue_avg = 0.5 * (ue_upper[i] + ue_upper[i-1]);
        let due_ds = (ue_upper[i] - ue_upper[i-1]) / ds;
        let cf_avg = 0.5 * (bl_upper[i].cf + bl_upper[i-1].cf);
        
        let pressure_term = (h_avg + 2.0) * theta_avg / ue_avg * due_ds;
        let lhs = dtheta_ds + pressure_term;
        let rhs = cf_avg / 2.0;
        let residual = rhs - lhs;
        
        println!("{:>4} {:>8.4} {:>10.2e} {:>10.2e} {:>10.4} {:>10.4} {:>10.2e}",
            i, x_upper[i], dtheta_ds, pressure_term, rhs, lhs, residual);
    }
}
