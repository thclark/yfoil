//! Find where the large DIJ values occur

use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;
    let le = airfoil.le_index;

    println!("Airfoil: {} panels, LE at index {}", n, le);

    let inviscid = solve_inviscid(&airfoil);
    let dij = inviscid.get_dij().expect("DIJ should be computed");

    // Find the 10 largest |DIJ| values
    let mut entries: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            entries.push((i, j, dij[(i, j)].abs()));
        }
    }
    entries.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());

    println!("\nTop 20 largest |DIJ| values:");
    println!("{:>6} {:>6} {:>12} {:>8} {:>8} {:>8} {:>8}",
             "i", "j", "|DIJ|", "x_i", "x_j", "surf_i", "surf_j");
    for (i, j, val) in entries.iter().take(20) {
        let surf_i = if *i <= le { "upper" } else { "lower" };
        let surf_j = if *j <= le { "upper" } else { "lower" };
        println!("{:6} {:6} {:12.4} {:8.4} {:8.4} {:>8} {:>8}",
                 i, j, val, airfoil.x[*i], airfoil.x[*j], surf_i, surf_j);
    }

    // Check diagonal elements specifically
    println!("\nDiagonal elements (self-influence):");
    println!("{:>6} {:>12} {:>8} {:>8}", "i", "DIJ[i,i]", "x", "surface");
    for i in [0, 1, 2, le-1, le, le+1, n-3, n-2, n-1] {
        let surf = if i <= le { "upper" } else { "lower" };
        println!("{:6} {:12.4} {:8.4} {:>8}", i, dij[(i, i)], airfoil.x[i], surf);
    }

    // Check row sums (total influence at each panel)
    println!("\nRow sums (total influence at each panel):");
    let mut row_sums: Vec<(usize, f64)> = Vec::new();
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            sum += dij[(i, j)].abs();
        }
        row_sums.push((i, sum));
    }
    row_sums.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    println!("{:>6} {:>12} {:>8} {:>8}", "i", "row_sum", "x", "surface");
    for (i, sum) in row_sums.iter().take(10) {
        let surf = if *i <= le { "upper" } else { "lower" };
        println!("{:6} {:12.4} {:8.4} {:>8}", i, sum, airfoil.x[*i], surf);
    }
}
