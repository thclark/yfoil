//! Debug the trailing edge panel handling

use yfoil::geometry::{create_paneled_airfoil, read_dat_file};

fn main() {
    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;

    println!("Airfoil: {} panels, sharp_te: {}", n, airfoil.sharp_te);
    println!();

    // TE geometry
    println!("=== Trailing Edge Geometry ===");
    println!("Node 0 (upper TE): ({:.6}, {:.6})", airfoil.x[0], airfoil.y[0]);
    println!("Node {} (lower TE): ({:.6}, {:.6})", n - 1, airfoil.x[n - 1], airfoil.y[n - 1]);

    let te_gap = ((airfoil.x[n - 1] - airfoil.x[0]).powi(2)
        + (airfoil.y[n - 1] - airfoil.y[0]).powi(2))
    .sqrt();
    println!("TE gap: {:.10}", te_gap);

    // Total arc length
    let s_total = airfoil.s[n - 1] - airfoil.s[0];
    let seps = s_total * 1.0e-5;
    println!("s_total: {:.6}", s_total);
    println!("seps (tolerance): {:.10}", seps);
    println!("te_gap < seps: {}", te_gap < seps);

    // Panel 0 and panel n-1
    println!();
    println!("=== Panel lengths near TE ===");
    let ds0 = ((airfoil.x[1] - airfoil.x[0]).powi(2) + (airfoil.y[1] - airfoil.y[0]).powi(2)).sqrt();
    let ds_n1 = ((airfoil.x[0] - airfoil.x[n - 1]).powi(2) + (airfoil.y[0] - airfoil.y[n - 1]).powi(2)).sqrt();
    let ds_n2 = ((airfoil.x[n - 1] - airfoil.x[n - 2]).powi(2) + (airfoil.y[n - 1] - airfoil.y[n - 2]).powi(2)).sqrt();
    println!("Panel 0 length (node 0 to 1): {:.6}", ds0);
    println!("Panel n-1 length (node n-1 to 0): {:.6}", ds_n1);
    println!("Panel n-2 length (node n-2 to n-1): {:.6}", ds_n2);

    // Check what XFOIL uses for JM, JQ at boundaries
    println!();
    println!("=== Boundary node indices (yfoil 0-based) ===");
    println!("At jo=0: jm=0, jp=1, jq=2");
    println!("At jo=1: jm=0, jp=2, jq=3");
    println!("At jo=n-2={}: jm={}, jp={}, jq={}", n - 2, n - 3, n - 1, n - 1);
    println!("At jo=n-1={}: jm={}, jp=0, jq=1", n - 1, n - 2);

    // LE geometry for comparison
    let le = airfoil.le_index;
    println!();
    println!("=== Leading Edge (index {}) ===", le);
    println!("LE node: ({:.6}, {:.6})", airfoil.x[le], airfoil.y[le]);
    let ds_le = ((airfoil.x[le + 1] - airfoil.x[le]).powi(2)
        + (airfoil.y[le + 1] - airfoil.y[le]).powi(2))
    .sqrt();
    println!("LE panel length: {:.6}", ds_le);
}
