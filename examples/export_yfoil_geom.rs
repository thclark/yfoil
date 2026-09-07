//! Export yfoil geometry to XFOIL format

use yfoil::geometry::{naca_4digit, panel_foil};

fn main() {
    // Generate NACA 0012 with 160 panels using yfoil
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = panel_foil(&geom);

    // Output in XFOIL .dat format (Selig format)
    // First line: name
    // Then: x y pairs
    println!("NACA 0012 from yfoil");
    for i in 0..airfoil.n_foil_nodes {
        println!("{:12.7} {:12.7}", airfoil.x[i], airfoil.y[i]);
    }

    // Also output to stderr for debugging
    eprintln!("Exported {} nodes", airfoil.n_foil_nodes);
    eprintln!("First node: ({:.6}, {:.6})", airfoil.x[0], airfoil.y[0]);
    eprintln!(
        "Last node: ({:.6}, {:.6})",
        airfoil.x[airfoil.n_foil_nodes - 1],
        airfoil.y[airfoil.n_foil_nodes - 1]
    );
    eprintln!(
        "LE node {}: ({:.6}, {:.6})",
        airfoil.i_le_node, airfoil.x[airfoil.i_le_node], airfoil.y[airfoil.i_le_node]
    );
}
