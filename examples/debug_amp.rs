use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::bl::{FlowConditions, NewtonConfig, hkin, dampl};
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::march_newton;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (_, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    println!("=== Arc Length Analysis ===\n");
    println!("Total arc length (upper surface): {:.6}", s_upper.last().unwrap());
    println!("Chord: 1.0");
    println!();

    println!("=== Step Size (ds) and Amplification Integration ===");
    println!("{:>5} {:>8} {:>8} {:>8} {:>8} {:>8} {:>10} {:>10}",
        "Stn", "s", "ds", "Hk", "ax", "ax*ds", "n_amp", "dn_actual");

    let mut prev_n = 0.0;
    for i in 25..55 {
        if i >= bl_upper.len() { break; }
        let ds = if i > 0 { s_upper[i] - s_upper[i-1] } else { 0.0 };
        let rt = ue_upper[i] * bl_upper[i].theta / cond.nu;
        let (hk, _, _) = hkin(bl_upper[i].h, cond.msq);
        let (ax, _, _, _) = dampl(hk, bl_upper[i].theta, rt);
        let ax_ds = ax * ds;
        let dn_actual = bl_upper[i].n_amp - prev_n;
        prev_n = bl_upper[i].n_amp;

        println!("{:>5} {:>8.4} {:>8.5} {:>8.4} {:>8.2} {:>8.4} {:>10.4} {:>10.4}",
            i, s_upper[i], ds, hk, ax, ax_ds, bl_upper[i].n_amp, dn_actual);
    }

    // Check if our dampl matches XFOIL's dampl at specific points
    println!("\n=== DAMPL Verification ===");
    println!("Testing dampl at known points:");

    // Test point 1: Hk = 2.5, θ = 0.001, Rθ = 500
    let hk_test = 2.5;
    let theta_test = 0.001;
    let rt_test = 500.0;
    let (ax, _, _, _) = dampl(hk_test, theta_test, rt_test);
    println!("  Hk=2.5, θ=0.001, Rθ=500 -> ax = {:.4}", ax);

    // Test point 2: Hk = 3.0, θ = 0.001, Rθ = 700
    let hk_test = 3.0;
    let rt_test = 700.0;
    let (ax, _, _, _) = dampl(hk_test, theta_test, rt_test);
    println!("  Hk=3.0, θ=0.001, Rθ=700 -> ax = {:.4}", ax);

    // Test point 3: Hk = 3.5, θ = 0.001, Rθ = 800
    let hk_test = 3.5;
    let rt_test = 800.0;
    let (ax, _, _, _) = dampl(hk_test, theta_test, rt_test);
    println!("  Hk=3.5, θ=0.001, Rθ=800 -> ax = {:.4}", ax);

    // Expected XFOIL values at these points (from XFOIL code/papers)
    println!("\nExpected XFOIL values (approximate):");
    println!("  Hk=2.5, Rθ=500 -> ax ≈ 1-3");
    println!("  Hk=3.0, Rθ=700 -> ax ≈ 5-10");
    println!("  Hk=3.5, Rθ=800 -> ax ≈ 15-25");
}
