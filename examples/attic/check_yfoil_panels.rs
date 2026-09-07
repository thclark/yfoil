//! Check yfoil panel coordinates near LE

use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = panel_foil(&geom);
    let inviscid = solve_inviscid(&airfoil);

    let qinv = inviscid.velocity_at_nodes(0.0);

    println!("yfoil panel coordinates near LE (le_index = {})", airfoil.le_index);
    println!("Panel    X              Y              gamma");

    let le = airfoil.le_index;
    for i in (le.saturating_sub(5))..=(le + 5).min(airfoil.n - 1) {
        let x = airfoil.x[i];
        let y = airfoil.y[i];
        let g = qinv[i];
        let marker = if i == le { " <-- LE" } else { "" };
        println!("{:4}  {:14.10e}  {:14.10e}  {:+12.8}{}", i, x, y, g, marker);
    }

    // Also show midpoint velocities for comparison
    println!("\nyfoil MIDPOINT velocities (velocity_at_alpha):");
    let qmid = inviscid.velocity_at_alpha(0.0);
    for i in (le.saturating_sub(5))..=(le + 5).min(airfoil.n - 1) {
        let marker = if i == le { " <-- LE" } else { "" };
        println!("{:4}  qmid = {:+12.8}{}", i, qmid[i], marker);
    }
}
