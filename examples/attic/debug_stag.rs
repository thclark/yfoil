//! Debug stagnation point detection

use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = panel_foil(&geom);
    let solution = solve_inviscid(&airfoil);

    let qinv = solution.velocity_at_nodes(0.0);

    println!("LE index = {}", airfoil.le_index);
    println!("\nVelocities near LE:");
    for i in 77..85.min(airfoil.n) {
        println!("  Node {}: x={:.6}, vel={:+.6}", i, airfoil.x[i], qinv[i]);
    }

    // Find sign change
    println!("\nLooking for sign change:");
    for i in 77..84 {
        if qinv[i] > 0.0 && qinv[i+1] <= 0.0 {
            println!("  Sign change at i={}: vel[{}]={:+.6}, vel[{}]={:+.6}",
                     i, i, qinv[i], i+1, qinv[i+1]);
        }
    }
}
