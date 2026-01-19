use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface, extract_lower_surface};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    println!("Stagnation index: {} (LE index: {})", stag_idx, airfoil.le_index);
    println!("Stagnation x: {:.6}", airfoil.x[stag_idx]);

    // Extract surfaces
    let (x_upper, _y, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (x_lower, _, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);

    println!("\nUpper surface: {} stations", x_upper.len());
    println!("  Start: x={:.6}, s={:.6}, ue={:.6}", x_upper[0], s_upper[0], ue_upper[0]);
    println!("  End (TE): x={:.6}, s={:.6}, ue={:.6}",
             x_upper.last().unwrap(), s_upper.last().unwrap(), ue_upper.last().unwrap());

    println!("\nLower surface: {} stations", x_lower.len());
    println!("  Start: x={:.6}, s={:.6}, ue={:.6}", x_lower[0], s_lower[0], ue_lower[0]);
    println!("  End (TE): x={:.6}, s={:.6}, ue={:.6}",
             x_lower.last().unwrap(), s_lower.last().unwrap(), ue_lower.last().unwrap());

    // March BL
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();

    let upper_bl = march_newton(&ue_upper, &s_upper, &cond, &config);
    let lower_bl = march_newton(&ue_lower, &s_lower, &cond, &config);

    println!("\n=== Upper BL March ===");
    println!("Station 0 (near stag): x={:.6}, theta={:.6e}, H={:.4}, Ue={:.4}",
             x_upper[0], upper_bl[0].theta, upper_bl[0].h, upper_bl[0].ue);
    println!("Station 5: x={:.6}, theta={:.6e}, H={:.4}",
             x_upper[5], upper_bl[5].theta, upper_bl[5].h);
    println!("Station 10: x={:.6}, theta={:.6e}, H={:.4}",
             x_upper[10], upper_bl[10].theta, upper_bl[10].h);
    if let Some(last) = upper_bl.last() {
        let last_idx = upper_bl.len()-1;
        println!("Station {} (TE): x={:.6}, theta={:.6e}, H={:.4}, Cf={:.6}",
                 last_idx, x_upper[last_idx], last.theta, last.h, last.cf);
    }

    println!("\n=== Lower BL March ===");
    println!("Station 0 (stag): x={:.6}, theta={:.6e}, H={:.4}, Ue={:.4}",
             x_lower[0], lower_bl[0].theta, lower_bl[0].h, lower_bl[0].ue);
    if let Some(last) = lower_bl.last() {
        let last_idx = lower_bl.len()-1;
        println!("Station {} (TE): x={:.6}, theta={:.6e}, H={:.4}, Cf={:.6}",
                 last_idx, x_lower[last_idx], last.theta, last.h, last.cf);
    }

    // Compare with XFOIL (from dump file)
    println!("\n=== XFOIL Reference (from dump) ===");
    println!("Upper TE: theta=0.001998, H=1.59, Cf=0.002044, Ue=0.886");
    println!("Near LE (s~1): theta~0.000033");
    println!("Combined wake theta: 0.003996");

    println!("\n=== Ratio ===");
    if let Some(last) = upper_bl.last() {
        println!("YFoil/XFOIL theta ratio at TE: {:.2}x", last.theta / 0.001998);
    }

    // Print some stations for detailed comparison
    println!("\n=== Detailed Upper Surface ===");
    println!("{:>8} {:>10} {:>12} {:>8} {:>8} {:>8}", "idx", "x", "theta", "H", "Cf", "Ue");
    for (i, bl) in upper_bl.iter().enumerate().step_by(10) {
        if i < x_upper.len() {
            println!("{:>8} {:>10.5} {:>12.6e} {:>8.4} {:>8.5} {:>8.4}",
                     i, x_upper[i], bl.theta, bl.h, bl.cf, bl.ue);
        }
    }
    if let Some(last) = upper_bl.last() {
        let i = upper_bl.len() - 1;
        println!("{:>8} {:>10.5} {:>12.6e} {:>8.4} {:>8.5} {:>8.4}",
                 i, x_upper[i], last.theta, last.h, last.cf, last.ue);
    }
}
