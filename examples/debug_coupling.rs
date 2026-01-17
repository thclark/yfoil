//! Debug viscous-inviscid coupling to understand source correction magnitude

use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{analyze_airfoil, FlowConditions, SolverConfig};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    println!("=== Analyzing coupling behavior ===\n");

    // Just run inviscid first
    let inviscid = solve_inviscid(&airfoil);
    let qinv_nodes = inviscid.velocity_at_nodes(0.0);

    println!("Inviscid velocities near LE (first 5 upper surface nodes):");
    let le_idx = airfoil.le_index;
    for j in 0..5 {
        let i = le_idx - j;
        println!("  Node {}: x={:.6}, qinv={:+.6}", i, airfoil.x[i], qinv_nodes[i]);
    }

    // Get DIJ matrix and check its magnitude
    if let Some(dij) = inviscid.get_dij() {
        println!("\nDIJ matrix sample (row 79, columns 78-82):");
        for j in 78..=82.min(dij.ncols() - 1) {
            println!("  DIJ[79,{}] = {:+.6e}", j, dij[(79, j)]);
        }

        // Check DIJ row sum for panel 79
        let mut row_sum = 0.0;
        for j in 0..dij.ncols() {
            row_sum += dij[(79, j)].abs();
        }
        println!("\n  |DIJ[79,:]| sum = {:.6}", row_sum);
    }

    println!("\n=== Running full coupled analysis ===\n");

    // Run with full convergence
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = SolverConfig::default();

    let result = analyze_airfoil(&airfoil, 0.0f64.to_radians(), &cond, &config, None);

    println!("Final result:");
    println!("  CL = {:.6}", result.cl);
    println!("  CD = {:.6}", result.cd);
    println!("  xtr_upper = {:.4}", result.xtr_upper);
    println!("  xtr_lower = {:.4}", result.xtr_lower);
    println!("  converged = {}", result.converged);
    println!("  iterations = {}", result.iterations);

    println!("\nXFOIL reference (same geometry):");
    println!("  CD = 0.00532");
    println!("  xtr = 0.687");
}
