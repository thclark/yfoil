//! Export yfoil geometry to XFOIL format

use yfoil::geometry::{create_paneled_airfoil, naca_4digit};

fn main() {
    // Generate NACA 0012 with 160 panels using yfoil
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    // Output in XFOIL .dat format (Selig format)
    // First line: name
    // Then: x y pairs
    println!("NACA 0012 from yfoil");
    for i in 0..airfoil.n {
        println!("{:12.7} {:12.7}", airfoil.x[i], airfoil.y[i]);
    }

    // Also output to stderr for debugging
    eprintln!("Exported {} nodes", airfoil.n);
    eprintln!("First node: ({:.6}, {:.6})", airfoil.x[0], airfoil.y[0]);
    eprintln!("Last node: ({:.6}, {:.6})", airfoil.x[airfoil.n-1], airfoil.y[airfoil.n-1]);
    eprintln!("LE node {}: ({:.6}, {:.6})", airfoil.le_index,
              airfoil.x[airfoil.le_index], airfoil.y[airfoil.le_index]);
}
