use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface, extract_lower_surface, compute_mass_defect, solve_boundary_layer, ViscalConfig};
use yfoil::bl::{NewtonConfig, march_newton};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    
    println!("Stagnation index: {}", stag_idx);
    
    // Extract surfaces
    let (_, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (_, _, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);
    
    println!("Upper surface: {} stations", s_upper.len());
    println!("Lower surface: {} stations", s_lower.len());
    
    // March BL
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &config);
    
    // Print some BL values
    if let Some(te) = bl_upper.last() {
        println!("\nUpper TE: theta={:.6e}, dstar={:.6e}, ue={:.6}", te.theta, te.dstar, te.ue);
        println!("Mass defect at upper TE: ue*dstar = {:.6e}", te.ue * te.dstar);
    }
    if let Some(te) = bl_lower.last() {
        println!("Lower TE: theta={:.6e}, dstar={:.6e}, ue={:.6}", te.theta, te.dstar, te.ue);
        println!("Mass defect at lower TE: ue*dstar = {:.6e}", te.ue * te.dstar);
    }
    
    // Compute mass defect
    let mass = compute_mass_defect(&airfoil, &bl_upper, &bl_lower, &ue_upper, &ue_lower, stag_idx);
    
    println!("\n=== Mass defect array ===");
    println!("mass[0] (upper TE): {:.6e}", mass[0]);
    println!("mass[1]: {:.6e}", mass[1]);
    println!("mass[{}] (near LE): {:.6e}", stag_idx, mass[stag_idx]);
    println!("mass[{}] (lower TE): {:.6e}", mass.len()-1, mass[mass.len()-1]);
    
    let max_mass = mass.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_mass = mass.iter().cloned().fold(f64::INFINITY, f64::min);
    let sum_mass: f64 = mass.iter().sum();
    println!("\nMax mass: {:.6e}", max_mass);
    println!("Min mass: {:.6e}", min_mass);
    println!("Sum mass: {:.6e}", sum_mass);
    
    // Check how many are nonzero
    let nonzero = mass.iter().filter(|&&m| m.abs() > 1e-10).count();
    println!("Nonzero entries: {}/{}", nonzero, mass.len());
}
