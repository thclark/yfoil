use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface, extract_lower_surface, compute_mass_defect};
use yfoil::bl::{NewtonConfig, march_newton};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    
    // Extract surfaces and march BL
    let (_, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (_, _, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);
    
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &config);
    
    // Compute mass defect
    let mass = compute_mass_defect(&airfoil, &bl_upper, &bl_lower, &ue_upper, &ue_lower, stag_idx);
    
    // Use inviscid method to compute dq
    if let Some(dq) = inviscid.velocity_from_mass_defect_unsigned(&mass) {
        println!("=== dq from velocity_from_mass_defect_unsigned ===");
        println!("dq[0] (upper TE): {:.6e}", dq[0]);
        println!("dq[1]: {:.6e}", dq[1]);
        println!("dq[10]: {:.6e}", dq[10]);
        println!("dq[{}] (near LE): {:.6e}", stag_idx, dq[stag_idx]);
        println!("dq[{}] (lower TE): {:.6e}", dq.len()-1, dq[dq.len()-1]);
        
        let max_dq = dq.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_dq = dq.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("\nMax dq: {:.6e}", max_dq);
        println!("Min dq: {:.6e}", min_dq);
        
        // What should the TE velocity become?
        let qinv_te = qinv[0].abs();
        let ue_new_te = qinv_te + dq[0];
        println!("\nInviscid TE velocity: {:.6}", qinv_te);
        println!("dq at TE: {:.6}", dq[0]);
        println!("New Ue at TE would be: {:.6}", ue_new_te);
        println!("XFOIL viscous TE: 0.886");
        
    } else {
        println!("Failed to compute dq!");
    }
    
    // Also manually compute dq[0] to verify
    if let Some(dij) = inviscid.get_dij() {
        let mut dq0_manual = 0.0;
        for j in 0..mass.len() {
            dq0_manual -= dij[(0, j)] * mass[j];
        }
        println!("\nManual dq[0]: {:.6e}", dq0_manual);
        
        // Show breakdown by contributions
        println!("\nTop 5 contributions to dq[0]:");
        let mut contribs: Vec<(usize, f64)> = (0..mass.len())
            .map(|j| (j, -dij[(0, j)] * mass[j]))
            .collect();
        contribs.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());
        for i in 0..5 {
            let (j, c) = contribs[i];
            println!("  j={}: DIJ[0,{}]={:.3}, mass[{}]={:.6e}, contrib={:.6e}",
                     j, j, dij[(0, j)], j, mass[j], c);
        }
    }
}
