//! Debug viscous-inviscid coupling

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{
    extract_lower_surface, extract_upper_surface, find_stagnation_point,
    solve_boundary_layer, solve_viscous, ViscalConfig,
};

fn main() {
    println!("=== Debug viscous-inviscid coupling ===\n");

    let geom = naca_4digit("0012", 80).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    // First check inviscid solution
    let inviscid = solve_inviscid(&airfoil);

    let alpha = 0.0_f64.to_radians();
    let vel_inv = inviscid.velocity_at_alpha(alpha);

    println!("Inviscid solution at α=0:");
    println!("  Panels: {}", airfoil.n);
    println!("  LE index: {}", airfoil.le_index);
    println!("  Velocity range: {:.4} to {:.4}",
             vel_inv.iter().cloned().fold(f64::INFINITY, f64::min),
             vel_inv.iter().cloned().fold(f64::NEG_INFINITY, f64::max));

    let stag_idx = find_stagnation_point(&airfoil, &vel_inv);
    println!("  Stagnation index: {} (expected ~{})", stag_idx, airfoil.le_index);

    // Check stagnation point velocity
    println!("  Velocity at stag: {:.6}", vel_inv[stag_idx]);

    // Extract surfaces
    let (x_upper, _y_upper, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel_inv, stag_idx);
    let (x_lower, _y_lower, s_lower, ue_lower) = extract_lower_surface(&airfoil, &vel_inv, stag_idx);

    println!("\nUpper surface: {} stations", x_upper.len());
    println!("  x: {:.4} to {:.4}", x_upper.first().unwrap(), x_upper.last().unwrap());
    println!("  s: {:.4} to {:.4}", s_upper.first().unwrap(), s_upper.last().unwrap());
    println!("  ue[0]={:.4}, ue[-1]={:.4}", ue_upper[0], ue_upper.last().unwrap());

    println!("\nLower surface: {} stations", x_lower.len());
    println!("  x: {:.4} to {:.4}", x_lower.first().unwrap(), x_lower.last().unwrap());
    println!("  s: {:.4} to {:.4}", s_lower.first().unwrap(), s_lower.last().unwrap());
    println!("  ue[0]={:.4}, ue[-1]={:.4}", ue_lower[0], ue_lower.last().unwrap());

    // Check first few Ue values
    println!("\nFirst 5 upper surface ue:");
    for i in 0..5.min(ue_upper.len()) {
        println!("  s={:.4}, x={:.4}, ue={:.4}", s_upper[i], x_upper[i], ue_upper[i]);
    }

    println!("\nLast 5 upper surface ue:");
    let n = ue_upper.len();
    for i in (n-5).max(0)..n {
        println!("  s={:.4}, x={:.4}, ue={:.4}", s_upper[i], x_upper[i], ue_upper[i]);
    }

    // Run full solve
    let config = ViscalConfig {
        max_iter: 50,
        ..Default::default()
    };

    println!("\n=== Running solve_viscous ===");
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("\nResult:");
    println!("  Converged: {}", result.converged);
    println!("  Iterations: {}", result.iterations);
    println!("  CL: {:.4}", result.cl);
    println!("  CD: {:.5}", result.cd);
    println!("  CDf: {:.5}", result.cdf);
    println!("  CDp: {:.5}", result.cdp);
    println!("  xtr_upper: {:.3}", result.xtr_upper);
    println!("  xtr_lower: {:.3}", result.xtr_lower);

    // Check BL results
    println!("\nBL upper surface:");
    println!("  Stations: {}", result.bl.upper.len());
    if !result.bl.upper.is_empty() {
        let last = result.bl.upper.last().unwrap();
        println!("  θ_TE: {:.6}", last.theta);
        println!("  δ*_TE: {:.6}", last.dstar);
        println!("  H_TE: {:.3}", last.h);
    }

    println!("\nBL lower surface:");
    println!("  Stations: {}", result.bl.lower.len());
    if !result.bl.lower.is_empty() {
        let last = result.bl.lower.last().unwrap();
        println!("  θ_TE: {:.6}", last.theta);
        println!("  δ*_TE: {:.6}", last.dstar);
        println!("  H_TE: {:.3}", last.h);
    }

    println!("\nWake:");
    println!("  Stations: {}", result.bl.wake.len());
    if !result.bl.wake.is_empty() {
        let last = result.bl.wake.last().unwrap();
        println!("  θ_end: {:.6}", last.theta);
        println!("  H_end: {:.3}", last.h);
        println!("  ue_end: {:.4}", last.ue);
    }
}
