//! Debug script to check airfoil symmetry and inviscid solution

use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    // Create NACA 0012 airfoil with 160 panels
    let geom = naca_4digit("0012", 160).expect("Failed to create airfoil");
    let airfoil = create_paneled_airfoil(&geom);

    println!("=== Geometry Symmetry Check ===");
    println!("N panels: {}", airfoil.n);
    println!("LE index: {}", airfoil.le_index);

    // Check symmetry of geometry
    println!("\n=== First and Last 5 Nodes ===");
    println!("{:>6} {:>14} {:>14} (first 5)", "Idx", "x", "y");
    for i in 0..5 {
        println!("{:>6} {:>14.10} {:>14.10}", i, airfoil.x[i], airfoil.y[i]);
    }
    println!("...");
    for i in (airfoil.n - 5)..airfoil.n {
        println!("{:>6} {:>14.10} {:>14.10}", i, airfoil.x[i], airfoil.y[i]);
    }

    // Check symmetry around LE
    println!("\n=== Symmetry around LE (index {}) ===", airfoil.le_index);
    let le = airfoil.le_index;
    println!("{:>6} {:>14} {:>14} {:>14} {:>14}", "Offset", "x_upper", "y_upper", "x_lower", "y_lower");
    for offset in 0..5 {
        let upper = le - offset - 1;
        let lower = le + offset;
        if upper > 0 && lower < airfoil.n {
            println!("{:>6} {:>14.10} {:>14.10} {:>14.10} {:>14.10}",
                     offset, airfoil.x[upper], airfoil.y[upper],
                     airfoil.x[lower], airfoil.y[lower]);
            println!("       y diff: {:.2e}", (airfoil.y[upper] + airfoil.y[lower]).abs());
        }
    }

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let alpha = 0.0_f64;
    let qinv = inviscid.velocity_at_nodes(alpha);

    println!("\n=== Inviscid Velocity at TE ===");
    println!("Node 0 (upper TE):   x={:.10}, y={:.10}, qinv={:.10}",
             airfoil.x[0], airfoil.y[0], qinv[0]);
    println!("Node {} (lower TE): x={:.10}, y={:.10}, qinv={:.10}",
             airfoil.n - 1, airfoil.x[airfoil.n - 1], airfoil.y[airfoil.n - 1], qinv[airfoil.n - 1]);

    println!("\nTE velocity asymmetry: {:.4e}", (qinv[0].abs() - qinv[airfoil.n - 1].abs()));

    // Check velocity symmetry around LE
    println!("\n=== Velocity Symmetry around LE ===");
    for offset in 0..5 {
        let upper = le - offset - 1;
        let lower = le + offset;
        if upper > 0 && lower < airfoil.n {
            println!("Offset {}: v_upper={:.10}, v_lower={:.10}, diff={:.4e}",
                     offset, qinv[upper], qinv[lower],
                     (qinv[upper].abs() - qinv[lower].abs()));
        }
    }

    // Check velocity at nodes 1-5 and n-5 to n-1
    println!("\n=== Velocity at TE region ===");
    println!("{:>6} {:>14} {:>14}", "Idx", "qinv", "x");
    for i in 0..5 {
        println!("{:>6} {:>14.10} {:>14.10}", i, qinv[i], airfoil.x[i]);
    }
    println!("...");
    for i in (airfoil.n - 5)..airfoil.n {
        println!("{:>6} {:>14.10} {:>14.10}", i, qinv[i], airfoil.x[i]);
    }

    // XFOIL comparison
    println!("\n=== XFOIL Comparison ===");
    println!("XFOIL QINV at station 1 (upper TE): 0.7670558229");
    println!("XFOIL QINV at station 160 (lower TE): -0.7670558229");
    println!("YFoil qinv[0] (upper TE): {:.10}", qinv[0]);
    println!("YFoil qinv[{}] (lower TE): {:.10}", airfoil.n - 1, qinv[airfoil.n - 1]);
}
