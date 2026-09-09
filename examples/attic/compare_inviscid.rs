//! Compare raw inviscid velocities between XFOIL and yfoil

use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    // Generate NACA 0012 with 160 panels using yfoil
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = panel_foil(&geom);
    let solution = solve_inviscid(&airfoil);

    // Get inviscid velocity at alpha = 0
    let qinv = solution.velocity_at_alpha(0.0);
    let gam = solution.gamma_at_alpha(0.0);

    println!("yFoil INVISCID SOLUTION");
    println!("N = {}", airfoil.n);
    println!("LE index = {}", airfoil.le_index);
    println!();
    println!("--- GAM (node-based) and QINV (midpoint) ---");
    println!("I, X, Y, GAM, QINV");

    for i in 0..airfoil.n {
        println!("{:4}  {:12.7}  {:12.7}  {:14.10}  {:14.10}",
                 i+1, airfoil.x[i], airfoil.y[i], gam[i], qinv[i]);
    }

    // Focus on LE region (around index 80)
    println!();
    println!("=== LE REGION (panels 75-86) ===");
    for i in 74..86.min(airfoil.n) {
        println!("{:4}  x={:12.7}  GAM={:14.10}  QINV={:14.10}",
                 i+1, airfoil.x[i], gam[i], qinv[i]);
    }
}
