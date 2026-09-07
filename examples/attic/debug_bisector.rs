//! Debug bisector condition computation

use std::f64::consts::PI;
use yfoil::geometry::{panel_foil, read_dat_file};

fn main() {
    println!("=== Debug Bisector Condition ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;

    println!("N = {}, sharp_te = {}", n, airfoil.sharp_te);

    // TE point
    let x_te = 0.5 * (airfoil.x[0] + airfoil.x[n - 1]);
    let y_te = 0.5 * (airfoil.y[0] + airfoil.y[n - 1]);
    println!("\nTE point: ({:.6}, {:.6})", x_te, y_te);
    println!("Node 0: ({:.6}, {:.6})", airfoil.x[0], airfoil.y[0]);
    println!("Node {}: ({:.6}, {:.6})", n - 1, airfoil.x[n - 1], airfoil.y[n - 1]);

    // Panel tangent directions at TE
    // Upper surface: from node 0 to node 1
    let dx1 = airfoil.x[1] - airfoil.x[0];
    let dy1 = airfoil.y[1] - airfoil.y[0];
    let ds1 = (dx1 * dx1 + dy1 * dy1).sqrt();

    // Lower surface: from node n-2 to node n-1
    let dxn = airfoil.x[n - 1] - airfoil.x[n - 2];
    let dyn_te = airfoil.y[n - 1] - airfoil.y[n - 2];
    let dsn = (dxn * dxn + dyn_te * dyn_te).sqrt();

    println!("\nUpper TE tangent: ({:.6}, {:.6}), ds={:.6}", dx1 / ds1, dy1 / ds1, ds1);
    println!("Lower TE tangent: ({:.6}, {:.6}), ds={:.6}", dxn / dsn, dyn_te / dsn, dsn);

    // XFOIL-style angles
    // AG1 = ATAN2(-YP(1), -XP(1)) - reversed tangent at upper TE
    let ag1 = (-dy1 / ds1).atan2(-dx1 / ds1);

    // AG2 = ATANC(YP(N), XP(N), AG1) - tangent at lower TE
    let mut ag2 = (dyn_te / dsn).atan2(dxn / dsn);

    // Adjust for branch cut
    while ag2 - ag1 > PI {
        ag2 -= 2.0 * PI;
    }
    while ag2 - ag1 < -PI {
        ag2 += 2.0 * PI;
    }

    println!("\nAG1 = {:.6} ({:.2} deg)", ag1, ag1.to_degrees());
    println!("AG2 = {:.6} ({:.2} deg)", ag2, ag2.to_degrees());

    // Bisector angle
    let a_bis = 0.5 * (ag1 + ag2);
    let c_bis = a_bis.cos();
    let s_bis = a_bis.sin();

    println!("\nA_BIS = {:.6} ({:.2} deg)", a_bis, a_bis.to_degrees());
    println!("C_BIS = {:.6}", c_bis);
    println!("S_BIS = {:.6}", s_bis);

    // Bisector control point
    let ds_min = ds1.min(dsn);
    let bwt = 0.1;
    let x_bis = x_te - bwt * ds_min * c_bis;
    let y_bis = y_te - bwt * ds_min * s_bis;

    println!("\nDS_MIN = {:.6}", ds_min);
    println!("X_BIS = {:.6}", x_bis);
    println!("Y_BIS = {:.6}", y_bis);

    // The normal direction at bisector (perpendicular to bisector)
    let nx_bis = -s_bis;
    let ny_bis = c_bis;
    println!("\nBisector normal: ({:.6}, {:.6})", nx_bis, ny_bis);

    // Compare with XFOIL's expected values
    // For NACA 0012, the TE panels should be nearly horizontal on upper surface
    // and nearly horizontal on lower surface, so bisector should point along x-axis
    println!("\n=== Expected values (approximate for symmetric NACA 0012) ===");
    println!("For symmetric airfoil at zero y-offset TE:");
    println!("  AG1 ~ -pi/2 (upper TE points backward-up) or ~pi (backward)");
    println!("  AG2 ~ pi/2 (lower TE points forward-up) or ~-pi (backward)");
    println!("  A_BIS ~ 0 or ~pi (horizontal bisector)");
}
