//! Compare initial inviscid state between yfoil approaches

use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);

    // Get both types of velocities
    let qinv_nodes = solution.velocity_at_nodes(0.0);
    let qinv_midpt = solution.velocity_at_alpha(0.0);

    let le_idx = airfoil.le_index;

    println!("LE index = {}", le_idx);
    println!("=== Raw inviscid velocities near LE ===\n");
    println!("Node   X          Node_vel   Midpt_vel   Ratio");
    for i in (le_idx - 5).max(0)..=(le_idx + 5).min(airfoil.n - 1) {
        let ratio = if qinv_nodes[i].abs() > 0.001 {
            qinv_midpt[i] / qinv_nodes[i]
        } else {
            0.0
        };
        println!("{:4}  {:10.6}  {:+10.6}  {:+10.6}  {:8.3}",
                 i, airfoil.x[i], qinv_nodes[i], qinv_midpt[i], ratio);
    }

    println!();
    println!("=== Upper surface stations (from LE toward TE) ===");
    println!("Station  Node  X          Node_vel   Midpt_vel");
    // Upper surface goes from LE (index le_idx) toward TE (index 0)
    for j in 0..10 {
        let i = le_idx - j;
        if i > 0 {
            println!("{:4}    {:4}  {:10.6}  {:+10.6}  {:+10.6}",
                     j+2, i, airfoil.x[i], qinv_nodes[i], qinv_midpt[i]);
        }
    }

    println!();
    println!("XFOIL reference (inviscid GAM at same nodes):");
    println!("  Panel 80 (x=0.000385): GAM = 0.285");
    println!("  Panel 79 (x=0.00154):  GAM = 0.534");
    println!("  Panel 78 (x=0.00347):  GAM = 0.727");
}
