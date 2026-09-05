use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);

    println!("Velocity at key stations (alpha=0):");
    println!("{:>6} {:>10} {:>10} {:>12}", "idx", "x", "y", "qinv");

    // TE (index 0)
    println!("{:>6} {:>10.5} {:>10.5} {:>12.6}", 0, airfoil.x[0], airfoil.y[0], qinv[0]);

    // Near upper TE
    for i in [1, 2, 3, 5, 10].iter() {
        if *i < qinv.len() {
            println!("{:>6} {:>10.5} {:>10.5} {:>12.6}", i, airfoil.x[*i], airfoil.y[*i], qinv[*i]);
        }
    }

    // LE region
    let le = airfoil.le_index;
    println!("LE index: {}", le);
    for i in [le-2, le-1, le, le+1, le+2].iter() {
        if *i < qinv.len() {
            println!("{:>6} {:>10.5} {:>10.5} {:>12.6}", i, airfoil.x[*i], airfoil.y[*i], qinv[*i]);
        }
    }

    // Near lower TE
    let n = qinv.len();
    for i in [n-10, n-5, n-3, n-2, n-1].iter() {
        if *i < n {
            println!("{:>6} {:>10.5} {:>10.5} {:>12.6}", i, airfoil.x[*i], airfoil.y[*i], qinv[*i]);
        }
    }

    // Find max and min velocities
    let max_q = qinv.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_q = qinv.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("\nMax velocity: {:.6}", max_q);
    println!("Min velocity: {:.6}", min_q);

    // Check TE velocity (should be ~0.88 for NACA 0012)
    println!("\nExpected TE velocity from XFOIL: 0.8856");
    println!("YFoil TE velocity (index 0): {:.6}", qinv[0].abs());
    println!("YFoil TE velocity (index {}): {:.6}", n-1, qinv[n-1].abs());
}
