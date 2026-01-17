//! Debug what happens at transition

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Debug Transition Region ===\n");

    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let mut config = ViscalConfig::default();
    config.tol_cl = 5e-4;

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("=== Upper Surface Around Transition ===");
    println!("{:>4} {:>10} {:>10} {:>8} {:>8} {:>10} {:>8} {:>10}",
             "i", "theta", "dstar", "H", "Hk", "Cf", "Ue", "n_amp");

    for i in 18..result.bl.upper.len().min(35) {
        let st = &result.bl.upper[i];
        println!("{:4} {:10.6} {:10.6} {:8.3} {:8.3} {:10.6} {:8.4} {:10.4}",
                 i, st.theta, st.dstar, st.h, st.hk, st.cf, st.ue, st.n_amp);
    }

    // Find transition index
    println!("\n=== Transition Analysis ===");
    for i in 1..result.bl.upper.len() {
        if result.bl.upper[i].n_amp >= 9.0 && result.bl.upper[i-1].n_amp < 9.0 {
            println!("Transition detected between stations {} and {}", i-1, i);
            println!("  Station {}: n={:.2}, H={:.3}, theta={:.6}",
                     i-1, result.bl.upper[i-1].n_amp, result.bl.upper[i-1].h, result.bl.upper[i-1].theta);
            println!("  Station {}: n={:.2}, H={:.3}, theta={:.6}",
                     i, result.bl.upper[i].n_amp, result.bl.upper[i].h, result.bl.upper[i].theta);
            break;
        }
    }

    // Check where H drops below 1.5 (indicating turbulent BL)
    println!("\n=== H Jump Analysis ===");
    for i in 1..result.bl.upper.len() {
        if result.bl.upper[i].h < 1.5 && result.bl.upper[i-1].h > 2.0 {
            println!("H jump at station {} -> {}", i-1, i);
            println!("  H dropped from {:.3} to {:.3}", result.bl.upper[i-1].h, result.bl.upper[i].h);
            println!("  theta: {:.6} -> {:.6}", result.bl.upper[i-1].theta, result.bl.upper[i].theta);
            break;
        }
    }

    // Summary of problem
    println!("\n=== Problem Summary ===");
    println!("When H jumps from laminar (~3.8 in inverse mode) to turbulent,");
    println!("it's dropping to H=1.0 instead of H~1.4-1.6 (equilibrium turbulent).");
    println!("This causes theta to be incorrectly calculated in the turbulent region.");

    // What XFOIL does at transition
    println!("\n=== Expected Behavior ===");
    println!("XFOIL resets H to equilibrium turbulent value at transition.");
    println!("For turbulent BL: H_eq ≈ 1.4-1.6 depending on pressure gradient.");
    println!("The Newton solver should start with a good initial guess for H.");
}
