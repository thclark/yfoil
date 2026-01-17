//! Debug DIJ matrix values

use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let inviscid = solve_inviscid(&airfoil);

    if let Some(dij) = inviscid.get_dij() {
        println!("=== yfoil DIJ[79,78-82] ===");
        for j in 78..=82.min(dij.ncols() - 1) {
            println!("   79   {}  {:+.16E}", j, dij[(79, j)]);
        }

        println!("\n=== yfoil DIJ row 79 stats ===");
        let mut row_sum = 0.0;
        let mut row_abs_sum = 0.0;
        for j in 0..dij.ncols() {
            row_sum += dij[(79, j)];
            row_abs_sum += dij[(79, j)].abs();
        }
        println!("Row sum: {:+.6e}", row_sum);
        println!("Row abs sum: {:.6e}", row_abs_sum);
    }
}
