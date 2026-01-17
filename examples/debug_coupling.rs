//! Debug detailed coupling iterations

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{
    compute_mass_defect, compute_source_velocity_dij, extract_lower_surface,
    extract_upper_surface, find_stagnation_point, solve_boundary_layer, ViscalConfig,
};

fn main() {
    println!("=== Debug coupling iterations ===\n");

    let geom = naca_4digit("0012", 80).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    // Inviscid solution
    let inviscid = solve_inviscid(&airfoil);
    let alpha = 0.0_f64.to_radians();
    let qinv = inviscid.velocity_at_alpha(alpha);
    let gamma = inviscid.gamma_at_alpha(alpha);
    let gamma_total: f64 = gamma.iter().sum::<f64>() / airfoil.n as f64 * airfoil.chord;

    println!("Inviscid: qinv[LE]={:.4}, gamma_total={:.4}", qinv[airfoil.le_index], gamma_total);

    // Start with inviscid velocity
    let mut velocity = qinv.clone();

    // Run a few coupling iterations manually
    for iter in 0..5 {
        println!("\n=== Iteration {} ===", iter);

        // Find stagnation point
        let stag_idx = find_stagnation_point(&airfoil, &velocity);
        println!("Stagnation index: {} (LE={})", stag_idx, airfoil.le_index);

        // Extract surfaces
        let (x_upper, _, _, ue_upper) = extract_upper_surface(&airfoil, &velocity, stag_idx);
        let (x_lower, _, _, ue_lower) = extract_lower_surface(&airfoil, &velocity, stag_idx);

        println!("Upper: {} stations, ue[0]={:.4}, ue[-1]={:.4}",
                 x_upper.len(), ue_upper[0], ue_upper.last().unwrap());
        println!("Lower: {} stations, ue[0]={:.4}, ue[-1]={:.4}",
                 x_lower.len(), ue_lower[0], ue_lower.last().unwrap());

        // Solve BL
        let bl = solve_boundary_layer(
            &airfoil,
            &velocity,
            stag_idx,
            alpha,
            gamma_total,
            &cond,
            &config.newton,
            &config.wake,
        );

        // Print BL TE values
        if !bl.upper.is_empty() {
            let u = bl.upper.last().unwrap();
            println!("Upper TE: θ={:.6}, δ*={:.6}, H={:.3}, Hk={:.3}, Cf={:.6}",
                     u.theta, u.dstar, u.h, u.hk, u.cf);
        }
        if !bl.lower.is_empty() {
            let l = bl.lower.last().unwrap();
            println!("Lower TE: θ={:.6}, δ*={:.6}, H={:.3}, Hk={:.3}, Cf={:.6}",
                     l.theta, l.dstar, l.h, l.hk, l.cf);
        }

        // Compute mass defect
        let mass = compute_mass_defect(&airfoil, &bl.upper, &bl.lower, &ue_upper, &ue_lower, stag_idx);

        // Check mass defect magnitude
        let max_mass = mass.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_mass = mass.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("Mass defect: min={:.6}, max={:.6}", min_mass, max_mass);

        // Compute source velocity correction
        let dq = compute_source_velocity_dij(&inviscid, &mass, airfoil.le_index);
        let max_dq = dq.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_dq = dq.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("dq correction: min={:.4}, max={:.4}", min_dq, max_dq);

        // Update velocity
        for i in 0..airfoil.n {
            velocity[i] = qinv[i] + dq[i];
        }

        // Check velocity after update
        let v_max = velocity.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let v_min = velocity.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("Velocity after: min={:.4}, max={:.4}", v_min, v_max);

        // Check TE velocities specifically
        println!("Velocity[0] (TE upper)={:.4}", velocity[0]);
        println!("Velocity[{}] (TE lower)={:.4}", airfoil.n - 1, velocity[airfoil.n - 1]);
        println!("Velocity[{}] (LE)={:.4}", airfoil.le_index, velocity[airfoil.le_index]);
    }
}
