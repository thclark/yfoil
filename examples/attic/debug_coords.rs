//! Debug coordinate indexing

use std::fs::File;
use std::io::BufReader;
use yfoil::geometry::{panel_foil, Geometry};
use yfoil::panel::solve_inviscid;
use yfoil::solver::find_stagnation_point;

fn load_geometry(path: &str) -> Geometry {
    let file = File::open(path).expect("Failed to open geometry file");
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).expect("Failed to parse geometry JSON")
}

fn main() {
    let geom = load_geometry(".tmp/naca0012_xfoil_paneled.json");
    let airfoil = panel_foil(&geom);

    println!("Total nodes: {}", airfoil.n);
    println!();

    // Print coordinates near LE
    println!("Coordinates near LE:");
    println!("{:>4} {:>14} {:>14}", "i", "x", "y");
    for i in 75..85 {
        println!("{:>4} {:>14.8e} {:>14.8e}", i, airfoil.x[i], airfoil.y[i]);
    }

    println!();

    // Find stagnation point
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    println!("Stagnation index: {}", stag_idx);
    println!("Stagnation x: {:.10e}", airfoil.x[stag_idx]);
    println!("Stagnation y: {:.10e}", airfoil.y[stag_idx]);
    println!();

    // Check what extract_upper_surface would do
    let start_idx = stag_idx.saturating_sub(1);
    println!("start_idx (stag_idx - 1): {}", start_idx);
    println!("x[start_idx]: {:.10e}", airfoil.x[start_idx]);
    println!("y[start_idx]: {:.10e}", airfoil.y[start_idx]);

    // Arc length from stag to start
    let dx = airfoil.x[start_idx] - airfoil.x[stag_idx];
    let dy = airfoil.y[start_idx] - airfoil.y[stag_idx];
    let arc_len = (dx * dx + dy * dy).sqrt();
    println!("Initial arc length: {:.10e}", arc_len);

    println!();
    println!("Velocities near LE:");
    for i in 75..85 {
        println!("{:>4} velocity: {:>14.8e}", i, qinv[i]);
    }
}
