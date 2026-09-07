//! Check DIJ matrix symmetry properties for symmetric airfoil

use yfoil::geometry::{panel_foil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    println!("=== DIJ Matrix Symmetry Check ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;
    let le = airfoil.le_index;

    println!("Airfoil: {} panels, LE at index {}", n, le);

    let inviscid = solve_inviscid(&airfoil);
    let dij = inviscid.get_dij().expect("DIJ should be computed");

    // For a symmetric airfoil, DIJ should have certain symmetry properties
    // Check if DIJ[i,j] relates to DIJ[n-1-i, n-1-j] by symmetry

    println!("\n=== Checking upper/lower symmetry ===");
    println!("For symmetric airfoil: DIJ[i,j] should relate to DIJ[n-1-i, n-1-j]\n");

    println!("{:>6} {:>6} {:>12} {:>12} {:>12} {:>10}",
             "i", "j", "DIJ[i,j]", "DIJ[i',j']", "diff", "ratio");

    for (i, j) in [(0, 0), (0, 1), (0, 5), (1, 0), (1, 1), (5, 5), (10, 10), (20, 20), (40, 40)] {
        // Symmetric counterpart
        let i_sym = n - 1 - i;
        let j_sym = n - 1 - j;

        if i < n && j < n && i_sym < n && j_sym < n {
            let d1 = dij[(i, j)];
            let d2 = dij[(i_sym, j_sym)];
            let diff = d1 - d2;
            let ratio = if d2.abs() > 1e-10 { d1 / d2 } else { f64::NAN };
            println!("{:6} {:6} {:12.4} {:12.4} {:12.4} {:10.4}",
                     i, j, d1, d2, diff, ratio);
        }
    }

    // Check row sums - for symmetric airfoil, row sums should be symmetric
    println!("\n=== Row sum symmetry ===");
    println!("{:>6} {:>12} {:>12} {:>12}",
             "i", "row_sum", "sym_row_sum", "diff");

    for i in [0, 1, 2, 5, 10, 20, le, le + 1, n - 2, n - 1] {
        let i_sym = n - 1 - i;
        if i < n && i_sym < n {
            let sum_i: f64 = (0..n).map(|j| dij[(i, j)]).sum();
            let sum_sym: f64 = (0..n).map(|j| dij[(i_sym, j)]).sum();
            println!("{:6} {:12.4} {:12.4} {:12.4}", i, sum_i, sum_sym, sum_i - sum_sym);
        }
    }

    // Check influence of upper panels on lower panels vs vice versa
    println!("\n=== Cross-surface influence ===");
    println!("Panel 0 (upper TE) influenced by:");
    let upper_influence: f64 = (0..=le).map(|j| dij[(0, j)]).sum();
    let lower_influence: f64 = (le + 1..n).map(|j| dij[(0, j)]).sum();
    println!("  Upper surface panels (0..{}): {:.4}", le, upper_influence);
    println!("  Lower surface panels ({}..{}): {:.4}", le + 1, n - 1, lower_influence);
    println!("  Ratio (lower/upper): {:.4}", lower_influence / upper_influence);

    println!("\nPanel {} (lower TE) influenced by:", n - 1);
    let upper_influence: f64 = (0..=le).map(|j| dij[(n - 1, j)]).sum();
    let lower_influence: f64 = (le + 1..n).map(|j| dij[(n - 1, j)]).sum();
    println!("  Upper surface panels (0..{}): {:.4}", le, upper_influence);
    println!("  Lower surface panels ({}..{}): {:.4}", le + 1, n - 1, lower_influence);
    println!("  Ratio (lower/upper): {:.4}", lower_influence / upper_influence);
}
