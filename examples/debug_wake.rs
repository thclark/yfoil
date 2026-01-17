//! Debug wake values

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Debug wake values ===\n");

    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("Converged: {}", result.converged);
    println!("Iterations: {}", result.iterations);
    println!("CL: {:.4}", result.cl);
    println!("CD: {:.4e}", result.cd);
    println!("CDf: {:.4e}", result.cdf);
    println!("CDp: {:.4e}", result.cdp);

    println!("\n=== Upper TE ===");
    if let Some(u) = result.bl.upper.last() {
        println!("θ={:.6}, δ*={:.6}, H={:.3}, Hk={:.3}", u.theta, u.dstar, u.h, u.hk);
    }

    println!("\n=== Lower TE ===");
    if let Some(l) = result.bl.lower.last() {
        println!("θ={:.6}, δ*={:.6}, H={:.3}, Hk={:.3}", l.theta, l.dstar, l.h, l.hk);
    }

    println!("\n=== Wake ===");
    println!("Wake stations: {}", result.bl.wake.len());
    for (i, w) in result.bl.wake.iter().enumerate() {
        println!("  {} θ={:.6}, δ*={:.6}, H={:.3}, Hk={:.3}, ue={:.4}",
                 i, w.theta, w.dstar, w.h, w.hk, w.ue);
    }

    // Compute what Squire-Young gives
    if let Some(wake_end) = result.bl.wake.last() {
        let theta_wake = wake_end.theta;
        let h_wake = wake_end.h;
        let ue_wake = wake_end.ue;

        let exponent = 0.5 * (5.0 + h_wake);
        let theta_inf = theta_wake * ue_wake.powf(exponent);
        let cd_sq = 2.0 * theta_inf / airfoil.chord;

        println!("\n=== Squire-Young ===");
        println!("θ_wake = {:.6}", theta_wake);
        println!("H_wake = {:.3}", h_wake);
        println!("Ue_wake = {:.4}", ue_wake);
        println!("Exponent = {:.3}", exponent);
        println!("Ue^exp = {:.4e}", ue_wake.powf(exponent));
        println!("θ_inf = {:.4e}", theta_inf);
        println!("CD_SY = {:.4e}", cd_sq);
    }
}
