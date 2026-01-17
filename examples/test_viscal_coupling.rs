//! Test viscous-inviscid coupling against XFOIL
//!
//! This compares yfoil's coupled solution against XFOIL's viscous results
//! to verify that the DIJ-based coupling improves agreement.

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Viscous-Inviscid Coupling Test ===\n");

    // Load XFOIL-paneled coordinates (use same geometry as XFOIL)
    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first to generate paneled coordinates");

    let airfoil = create_paneled_airfoil(&geom);

    println!("Airfoil: {} panels", airfoil.n);

    // Flow conditions matching XFOIL
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    // Run yfoil viscous solver with DIJ coupling
    let config = ViscalConfig {
        max_iter: 50,
        tol_cl: 1e-5,
        relax: 0.5,
        ..Default::default()
    };

    println!("\nRunning yfoil viscous solver...");
    let result = solve_viscous(&airfoil, 0.0, &cond, &config);

    println!("\n=== Results ===");
    println!("Converged: {}", result.converged);
    println!("Iterations: {}", result.iterations);
    println!("CL: {:.6}", result.cl);
    println!("CD: {:.6}", result.cd);
    println!("CDf: {:.6}", result.cdf);
    println!("CDp: {:.6}", result.cdp);
    println!("Transition upper (x/c): {:.4}", result.xtr_upper);
    println!("Transition lower (x/c): {:.4}", result.xtr_lower);

    // Compare against XFOIL values (from previous runs)
    println!("\n=== Comparison with XFOIL ===");
    println!("XFOIL transition: x/c = 0.687");
    println!(
        "yfoil transition: x/c = {:.4}",
        result.xtr_upper.min(result.xtr_lower)
    );

    let xfoil_trans = 0.687;
    let yfoil_trans = result.xtr_upper.min(result.xtr_lower);
    let trans_diff = 100.0 * (yfoil_trans - xfoil_trans) / xfoil_trans;
    println!("Difference: {:.2}%", trans_diff);

    // Check if coupling improved the solution
    if trans_diff.abs() < 5.0 {
        println!("\n✓ Good agreement with XFOIL transition location");
    } else if trans_diff.abs() < 10.0 {
        println!("\n○ Moderate agreement - coupling is working but may need tuning");
    } else {
        println!("\n✗ Significant difference - investigate coupling implementation");
    }

    // Print velocity correction magnitude
    let max_dq: f64 = result
        .dq_source
        .iter()
        .map(|x| x.abs())
        .fold(0.0, f64::max);
    println!("\nMax velocity correction: {:.6}", max_dq);
}
