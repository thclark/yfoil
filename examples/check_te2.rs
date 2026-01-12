use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    println!("Airfoil structure:");
    println!("  airfoil.x.len() = {} (nodes)", airfoil.x.len());
    println!("  airfoil.n = {} (panels)", airfoil.n);
    println!("  sharp_te = {}", airfoil.sharp_te);

    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let gamma = solution.gamma_at_alpha(0.0);

    println!("\nSolution arrays:");
    println!("  vel.len() = {} (panel midpoints)", vel.len());
    println!("  gamma.len() = {} (nodes)", gamma.len());

    // Check TE panel velocity vs adjacent panels
    println!("\n=== Panel velocities near TE ===");
    println!("(Panel i goes from node i to node i+1)");
    println!();

    // Upper surface near TE (panels 0, 1, 2 start from TE going to LE)
    for i in 0..5 {
        let ip1 = if i == airfoil.n - 1 { 0 } else { i + 1 };
        println!("Panel {:>3}: vel={:+.4}, node {} to {}, x: {:.4} to {:.4}",
            i, vel[i], i, ip1, airfoil.x[i], airfoil.x[ip1]);
    }
    println!("...");

    // Lower surface near TE (panels n-5 to n-1, ending at TE)
    for i in (airfoil.n - 5)..airfoil.n {
        let ip1 = if i == airfoil.n - 1 { 0 } else { i + 1 };
        println!("Panel {:>3}: vel={:+.4}, node {} to {}, x: {:.4} to {:.4}",
            i, vel[i], i, ip1, airfoil.x[i], airfoil.x[ip1]);
    }

    // Check gamma (node values) near TE
    println!("\n=== Node gamma values near TE ===");
    let n = gamma.len();
    for i in 0..3 {
        println!("Node {:>3}: gamma={:+.4}, x={:.4}, y={:.6}",
            i, gamma[i], airfoil.x[i], airfoil.y[i]);
    }
    println!("...");
    for i in (n - 3)..n {
        println!("Node {:>3}: gamma={:+.4}, x={:.4}, y={:.6}",
            i, gamma[i], airfoil.x[i], airfoil.y[i]);
    }

    // Force integration check
    use yfoil::forces::{calculate_cp, integrate_forces};
    let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();
    let coeffs = integrate_forces(&airfoil, &cp, 0.0);

    println!("\n=== Force integration ===");
    println!("CL  = {:.6}", coeffs.cl);
    println!("CM  = {:.6}", coeffs.cm);
    println!("CDp = {:.6} (should be 0 for symmetric airfoil at α=0)", coeffs.cdp);

    // Check which panels contribute most to CDp
    println!("\n=== Panel CDp contributions ===");
    let mut panel_cdp = vec![];
    for i in 0..airfoil.n {
        let ip1 = if i == airfoil.n - 1 { 0 } else { i + 1 };
        let dx = airfoil.x[ip1] - airfoil.x[i];
        let dy = airfoil.y[ip1] - airfoil.y[i];
        let cp_i = cp[i];
        let cn_i = cp_i * dx;
        let ca_i = -cp_i * dy;
        let cdp_i = cn_i * 0.0 + ca_i * 1.0; // at α=0, cdp = ca
        panel_cdp.push((i, cdp_i, dx, dy, cp_i));
    }

    // Sort by |CDp| contribution
    panel_cdp.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());

    println!("Top 5 panels by |CDp| contribution:");
    for &(i, cdp_i, dx, dy, cp_i) in panel_cdp.iter().take(5) {
        let ip1 = if i == airfoil.n - 1 { 0 } else { i + 1 };
        println!("  Panel {}: CDp={:+.6}, Cp={:+.4}, dx={:+.4}, dy={:+.4}, x:{:.3}->{:.3}",
            i, cdp_i, cp_i, dx, dy, airfoil.x[i], airfoil.x[ip1]);
    }
}
