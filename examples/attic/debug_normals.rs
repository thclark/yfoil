//! Debug normal vectors at TE panels

use yfoil::geometry::{panel_foil, read_dat_file};

fn main() {
    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;

    println!("=== Normal vectors and panel angles at TE ===\n");
    println!("{:>4} {:>10} {:>10} {:>10} {:>10} {:>10}", "j", "nx", "ny", "apanel", "x", "y");

    for j in [0, 1, 2, 3, n - 3, n - 2, n - 1] {
        let apanel = (-airfoil.ny[j]).atan2(-airfoil.nx[j]);
        println!("{:4} {:10.6} {:10.6} {:10.6} {:10.6} {:10.6}",
                 j, airfoil.nx[j], airfoil.ny[j], apanel, airfoil.x[j], airfoil.y[j]);
    }

    // Check panel directions
    println!("\n=== Panel directions (from node to next node) ===");
    for jo in [0, 1, 2, n - 2, n - 1] {
        let jp = if jo == n - 1 { 0 } else { jo + 1 };
        let dx = airfoil.x[jp] - airfoil.x[jo];
        let dy = airfoil.y[jp] - airfoil.y[jo];
        let ds = (dx * dx + dy * dy).sqrt();
        let tx = dx / ds;  // tangent x
        let ty = dy / ds;  // tangent y
        println!("Panel {} → {}: tangent=({:.6}, {:.6}), ds={:.6}",
                 jo, jp, tx, ty, ds);
    }

    // Verify normals are perpendicular to tangent
    println!("\n=== Normal · tangent (should be ~0) ===");
    for jo in [0, 1, 2, n - 2, n - 1] {
        let jp = if jo == n - 1 { 0 } else { jo + 1 };
        let dx = airfoil.x[jp] - airfoil.x[jo];
        let dy = airfoil.y[jp] - airfoil.y[jo];
        let ds = (dx * dx + dy * dy).sqrt();
        let tx = dx / ds;
        let ty = dy / ds;
        let dot = airfoil.nx[jo] * tx + airfoil.ny[jo] * ty;
        println!("Panel {}: nx*tx + ny*ty = {:.6}", jo, dot);
    }
}
