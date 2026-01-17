use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};

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

    println!("=== θ Comparison with Blasius ===");
    println!("Blasius: θ = 0.664 * sqrt(ν*s/Ue)");
    println!();
    println!("{:>8} {:>10} {:>12} {:>12} {:>10}", "x/c", "s", "θ_ours", "θ_blasius", "Ratio");
    
    let check_points = [0.10, 0.20, 0.30, 0.40, 0.50];
    for xc_target in check_points {
        // Find closest station
        let mut closest_idx = 0;
        let mut min_diff = f64::INFINITY;
        for (i, &x) in x_upper.iter().enumerate() {
            let diff = (x - xc_target).abs();
            if diff < min_diff {
                min_diff = diff;
                closest_idx = i;
            }
        }
        
        if closest_idx < bl_upper.len() {
            let s = s_upper[closest_idx];
            let ue = ue_upper[closest_idx];
            let theta_ours = bl_upper[closest_idx].theta;
            // Blasius: θ = 0.664 * sqrt(ν*x/Ue) where x is streamwise coordinate
            // For arc length s and varying Ue, use local Blasius-like scaling
            let theta_blasius = 0.664 * (cond.nu * s / ue).sqrt();
            let ratio = theta_ours / theta_blasius;
            
            println!("{:>8.2} {:>10.4} {:>12.2e} {:>12.2e} {:>10.3}", 
                x_upper[closest_idx], s, theta_ours, theta_blasius, ratio);
        }
    }
    
    println!();
    println!("Note: Blasius formula is for flat plate (dUe/ds = 0).");
    println!("On an airfoil with favorable gradient, θ is typically smaller.");
    
    // Show the actual dUe/ds at each point
    println!();
    println!("=== Pressure Gradient Effect ===");
    println!("{:>8} {:>12} {:>12}", "x/c", "dUe/ds", "Type");
    for xc_target in check_points {
        let mut closest_idx = 0;
        let mut min_diff = f64::INFINITY;
        for (i, &x) in x_upper.iter().enumerate() {
            let diff = (x - xc_target).abs();
            if diff < min_diff && i > 0 && i < ue_upper.len() - 1 {
                min_diff = diff;
                closest_idx = i;
            }
        }
        
        if closest_idx > 0 && closest_idx < ue_upper.len() - 1 {
            let due_ds = (ue_upper[closest_idx + 1] - ue_upper[closest_idx - 1]) 
                / (s_upper[closest_idx + 1] - s_upper[closest_idx - 1]);
            let grad_type = if due_ds > 0.01 { "favorable" } 
                else if due_ds < -0.01 { "adverse" } 
                else { "flat" };
            println!("{:>8.2} {:>12.4} {:>12}", x_upper[closest_idx], due_ds, grad_type);
        }
    }
}
