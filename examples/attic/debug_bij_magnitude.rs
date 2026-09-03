//! Debug BIJ vs AIJ magnitude comparison
//!
//! Compare the magnitudes of the vortex influence matrix (AIJ) and
//! source influence matrix (BIJ) to understand the 100x scaling issue.

use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    println!("=== BIJ vs AIJ Magnitude Analysis ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;

    println!("Airfoil: {} panels", n);

    // Get typical panel length
    let mut total_length = 0.0;
    for i in 0..n {
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        let ds = ((airfoil.x[ip1] - airfoil.x[i]).powi(2)
            + (airfoil.y[ip1] - airfoil.y[i]).powi(2))
        .sqrt();
        total_length += ds;
    }
    let avg_panel_length = total_length / n as f64;
    println!("Average panel length: {:.6}", avg_panel_length);
    println!("1/avg_panel_length: {:.1}", 1.0 / avg_panel_length);

    // Solve inviscid to get DIJ
    let inviscid = solve_inviscid(&airfoil);

    // Check DIJ statistics
    let dij = inviscid.get_dij().expect("DIJ should be computed");

    let mut dij_sum = 0.0;
    let mut dij_max: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let val = dij[(i, j)].abs();
            dij_sum += val;
            if val > dij_max {
                dij_max = val;
            }
        }
    }
    let mean_dij = dij_sum / (n * n) as f64;

    println!("\nDIJ matrix (dQtan/dSig):");
    println!("  Mean |DIJ|: {:.4}", mean_dij);
    println!("  Max |DIJ|: {:.4}", dij_max);

    // Check inviscid velocity (gamma) magnitudes for reference
    let qinv = inviscid.velocity_at_alpha(0.0);
    let mut qinv_max: f64 = 0.0;
    for &q in &qinv {
        if q.abs() > qinv_max {
            qinv_max = q.abs();
        }
    }
    println!("\nInviscid velocity (gamma) at alpha=0:");
    println!("  Max |Qinv|: {:.4}", qinv_max);

    // Analysis
    println!("\n=== Analysis ===");
    println!("For typical mass defect m* = 0.002:");
    println!("  Expected dQ = DIJ * m* ~ 0.01 (1% of Ue)");
    println!("  This implies DIJ ~ 5");
    println!("  Actual mean DIJ: {:.1}", mean_dij);
    println!("  Ratio: {:.0}x too large", mean_dij / 5.0);

    // Key insight: The extra factor might be related to 1/panel_length
    println!("\n=== Hypothesis ===");
    println!("1/avg_panel_length = {:.1}", 1.0 / avg_panel_length);
    println!("If DIJ is scaled by panel length:");
    println!("  Scaled DIJ = {:.2}", mean_dij * avg_panel_length);
}
