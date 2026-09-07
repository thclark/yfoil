//! Check station coordinates and velocities against XFOIL

use std::fs::File;
use std::io::BufReader;
use yfoil::geometry::{panel_foil, Geometry};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, find_stagnation_point};

fn load_geometry(path: &str) -> Geometry {
    let file = File::open(path).expect("Failed to open geometry file");
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).expect("Failed to parse geometry JSON")
}

fn main() {
    // Load same geometry as XFOIL
    let geom = load_geometry(".tmp/naca0012_xfoil_paneled.json");
    let airfoil = panel_foil(&geom);

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    println!("Stagnation point index: {}", stag_idx);
    println!("Total nodes: {}", airfoil.n);
    println!();

    let (x_upper, y_upper, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);

    println!("Upper surface: {} points", x_upper.len());
    println!();

    // Print first 15 stations matching XFOIL format
    // XFOIL IBL starts at 2 (IBL=1 is stagnation point itself)
    println!("{:>4} {:>12} {:>12} {:>12} {:>12}", "IBL", "X", "Y", "S", "UE");
    println!("{:->4} {:->12} {:->12} {:->12} {:->12}", "", "", "", "", "");

    for i in 0..x_upper.len().min(15) {
        println!(
            "{:>4} {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e}",
            i + 2,  // IBL index (1-based, starting from 2)
            x_upper[i],
            y_upper[i],
            s_upper[i],
            ue_upper[i]
        );
    }

    println!();
    println!("XFOIL station 2 reference:");
    println!("X = 0.000203, U = 0.01677, S = 0.03");
}
