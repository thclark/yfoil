//! Debug wake values for drag calculation

use yfoil::bl::FlowConditions;
use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = panel_foil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("=== Upper surface TE ===");
    if let Some(te) = result.bl.upper.last() {
        println!("  theta = {:.6e}", te.theta);
        println!("  H = {:.4}", te.h);
        println!("  Ue = {:.6}", te.ue);
    }

    println!("\n=== Lower surface TE ===");
    if let Some(te) = result.bl.lower.last() {
        println!("  theta = {:.6e}", te.theta);
        println!("  H = {:.4}", te.h);
        println!("  Ue = {:.6}", te.ue);
    }

    println!("\n=== Wake ({} stations) ===", result.bl.wake.len());
    if !result.bl.wake.is_empty() {
        println!("First wake station:");
        let first = &result.bl.wake[0];
        println!("  theta = {:.6e}", first.theta);
        println!("  H = {:.4}", first.h);
        println!("  Ue = {:.6}", first.ue);

        println!("\nLast wake station:");
        let last = result.bl.wake.last().unwrap();
        println!("  theta = {:.6e}", last.theta);
        println!("  H = {:.4}", last.h);
        println!("  Ue = {:.6}", last.ue);
    }

    println!("\n=== Drag calculation ===");
    println!("CD (Squire-Young) = {:.6}", result.cd);
    println!("CDf = {:.6}", result.cdf);
    println!("CDp = {:.6}", result.cdp);

    // Manual verification
    if let Some(wake_end) = result.bl.wake.last() {
        let theta = wake_end.theta;
        let h = wake_end.h;
        let ue = wake_end.ue;
        let exponent = 0.5 * (5.0 + h);
        let theta_inf = theta * ue.powf(exponent);
        let cd_manual = 2.0 * theta_inf / airfoil.chord;
        println!("\nManual Squire-Young:");
        println!("  exponent = 0.5 * (5 + {:.4}) = {:.4}", h, exponent);
        println!("  theta_inf = theta * Ue^exp = {:.6e} * {:.4}^{:.4} = {:.6e}",
                 theta, ue, exponent, theta_inf);
        println!("  CD = 2 * theta_inf / c = 2 * {:.6e} = {:.6e}",
                 theta_inf, cd_manual);
    }
}
