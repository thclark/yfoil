use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::bl::{FlowConditions, NewtonConfig};
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::march_newton;

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

    // Expected XFOIL values for NACA 0012 at Re=1e6, α=0°
    // These are approximate values from running XFOIL
    println!("=== BL Evolution Comparison ===\n");
    println!("NACA 0012, Re=1e6, α=0°\n");

    println!("{:>6} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10}",
        "x/c", "ue/U∞", "Hk_our", "Hk_exp", "θ_our", "θ_exp", "δ*_our");
    println!("{:-<85}", "");

    // Expected XFOIL values at select x/c locations
    // (approximate, from typical XFOIL runs)
    let expected = [
        (0.05, 1.25, 2.3, 1.5e-4),
        (0.10, 1.18, 2.35, 2.5e-4),
        (0.20, 1.12, 2.45, 4.0e-4),
        (0.30, 1.08, 2.55, 5.5e-4),
        (0.40, 1.05, 2.70, 7.0e-4),
        (0.50, 1.02, 2.90, 8.5e-4),
        (0.60, 0.98, 3.20, 1.0e-3),
    ];

    for &(xc_exp, ue_exp, hk_exp, th_exp) in &expected {
        // Find closest station to this x/c
        let mut closest_idx = 0;
        let mut min_diff = f64::INFINITY;
        for (i, &x) in x_upper.iter().enumerate() {
            let diff = (x - xc_exp).abs();
            if diff < min_diff {
                min_diff = diff;
                closest_idx = i;
            }
        }

        if closest_idx < bl_upper.len() {
            let ue_our = ue_upper[closest_idx];
            let hk_our = bl_upper[closest_idx].hk;
            let th_our = bl_upper[closest_idx].theta;
            let ds_our = bl_upper[closest_idx].dstar;

            println!("{:>6.2} {:>8.4} {:>8.4} {:>8.2} {:>10.2e} {:>10.2e} {:>10.2e}",
                x_upper[closest_idx], ue_our, hk_our, hk_exp, th_our, th_exp, ds_our);
        }
    }

    // Show pressure gradient
    println!("\n=== Pressure Gradient (dCp/dx proxy: due/dx) ===");
    println!("{:>6} {:>10} {:>10}", "x/c", "ue", "due/ds");
    for i in [10, 20, 30, 40, 50, 60, 70] {
        if i + 1 < bl_upper.len() {
            let due_ds = (ue_upper[i+1] - ue_upper[i]) / (s_upper[i+1] - s_upper[i]);
            println!("{:>6.3} {:>10.4} {:>10.4}",
                x_upper[i], ue_upper[i], due_ds);
        }
    }

    // Show XFOIL expected transition info
    println!("\n=== Expected XFOIL Transition ===");
    println!("For NACA 0012 at Re=1e6, α=0°, XFOIL typically gives:");
    println!("  Transition x/c ≈ 0.45-0.50");
    println!("  CD ≈ 0.0065-0.0070");
    println!();
    println!("Our current results:");

    // Find our transition
    for (i, r) in bl_upper.iter().enumerate() {
        if r.n_amp >= 9.0 {
            println!("  Transition at x/c = {:.4} (station {})", x_upper[i], i);
            break;
        }
    }

    // Key observation: compare Hk at the same x/c
    println!("\n=== Key Observation ===");
    println!("At x/c ≈ 0.50:");
    let idx_50 = x_upper.iter().position(|&x| x > 0.49 && x < 0.51).unwrap_or(40);
    println!("  Our Hk = {:.4}", bl_upper[idx_50].hk);
    println!("  Expected XFOIL Hk ≈ 2.8-3.0 (approaching separation)");
    println!("  Our n_amp = {:.4}", bl_upper[idx_50].n_amp);
    println!("  Expected XFOIL n_amp ≈ 8-9 (near transition)");
}
