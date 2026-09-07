//! Debug DZDM computation at trailing edge (control point 0)
//!
//! This traces through the exact computation to compare with XFOIL

use std::f64::consts::PI;
use yfoil::geometry::{panel_foil, read_dat_file};

const QOPI: f64 = 0.25 / PI;

fn main() {
    println!("=== Debug DZDM at TE (i=0) ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;

    println!("N = {}, sharp_te = {}", n, airfoil.sharp_te);

    // Compute panel angles
    let mut apanel = vec![0.0; n];
    for j in 0..n {
        apanel[j] = (-airfoil.ny[j]).atan2(-airfoil.nx[j]);
    }

    // Control point i=0 (first node, upper TE)
    let i = 0;
    let xi = airfoil.x[i];
    let yi = airfoil.y[i];
    println!("\nControl point i={}: ({:.6}, {:.6})", i, xi, yi);

    let mut dzdm = vec![0.0; n];

    // Print detailed info for first few panels and last few panels
    let debug_panels = vec![0, 1, 2, 3, n-3, n-2, n-1];

    println!("\n{:>4} {:>4} {:>4} {:>4} {:>10} {:>10} {:>10} {:>10}",
             "jo", "jp", "jm", "jq", "DSO", "DSM", "DSP", "contrib");

    for jo in 0..n {
        let jp = if jo == n - 1 { 0 } else { jo + 1 };
        let jm = if jo == 0 { jo } else { jo - 1 };
        let jq = if jo == n - 2 {
            n - 1
        } else {
            (jp + 1) % n
        };

        // Panel length
        let dso = ((airfoil.x[jp] - airfoil.x[jo]).powi(2)
            + (airfoil.y[jp] - airfoil.y[jo]).powi(2))
        .sqrt();
        if dso < 1e-14 {
            continue;
        }
        let dsio = 1.0 / dso;

        let apan = apanel[jo];

        // Vectors from panel nodes to control point
        let rx1 = xi - airfoil.x[jo];
        let ry1 = yi - airfoil.y[jo];
        let rx2 = xi - airfoil.x[jp];
        let ry2 = yi - airfoil.y[jp];

        // Unit tangent along panel
        let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
        let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

        // Transform to panel-local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        // Squared distances
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn = 1.0_f64;

        // Log and arctan terms at panel endpoints
        let (g1, t1) = if i == jo || rs1 < 1e-24 {
            (0.0, 0.0)
        } else {
            (rs1.ln(), (sgn * x1).atan2(sgn * yy))
        };

        let (g2, t2) = if i == jp || rs2 < 1e-24 {
            (0.0, 0.0)
        } else {
            (rs2.ln(), (sgn * x2).atan2(sgn * yy))
        };

        // Midpoint quantities
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = if rs0 > 1e-24 { rs0.ln() } else { 0.0 };
        let t0 = (sgn * x0).atan2(sgn * yy);

        // First half-panel (1 to 0)
        let dxinv_10 = if (x1 - x0).abs() > 1e-14 {
            1.0 / (x1 - x0)
        } else {
            0.0
        };

        let psum_10 = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
        let pdif_10 = if dxinv_10.abs() > 1e-14 {
            ((x1 + x0) * psum_10 + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy)
                * dxinv_10
        } else {
            0.0
        };

        let dsm = ((airfoil.x[jp] - airfoil.x[jm]).powi(2)
            + (airfoil.y[jp] - airfoil.y[jm]).powi(2))
        .sqrt();
        let dsim = if dsm > 1e-14 { 1.0 / dsm } else { 0.0 };

        // dPsi/dm contributions for first half-panel
        let contrib_jm_10 = QOPI * (-psum_10 * dsim + pdif_10 * dsim);
        let contrib_jo_10 = QOPI * (-psum_10 * dsio - pdif_10 * dsio);
        let contrib_jp_10 = QOPI * (psum_10 * (dsio + dsim) + pdif_10 * (dsio - dsim));

        dzdm[jm] += contrib_jm_10;
        dzdm[jo] += contrib_jo_10;
        dzdm[jp] += contrib_jp_10;

        // Second half-panel (0 to 2)
        let dxinv_02 = if (x0 - x2).abs() > 1e-14 {
            1.0 / (x0 - x2)
        } else {
            0.0
        };

        let psum_02 = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
        let pdif_02 = if dxinv_02.abs() > 1e-14 {
            ((x0 + x2) * psum_02 + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy)
                * dxinv_02
        } else {
            0.0
        };

        let dsp = ((airfoil.x[jq] - airfoil.x[jo]).powi(2)
            + (airfoil.y[jq] - airfoil.y[jo]).powi(2))
        .sqrt();
        let dsip = if dsp > 1e-14 { 1.0 / dsp } else { 0.0 };

        // dPsi/dm contributions for second half-panel
        let contrib_jo_02 = QOPI * (-psum_02 * (dsip + dsio) - pdif_02 * (dsip - dsio));
        let contrib_jp_02 = QOPI * (psum_02 * dsio - pdif_02 * dsio);
        let contrib_jq_02 = QOPI * (psum_02 * dsip + pdif_02 * dsip);

        dzdm[jo] += contrib_jo_02;
        dzdm[jp] += contrib_jp_02;
        dzdm[jq] += contrib_jq_02;

        // Total contribution from this panel
        let total_contrib = contrib_jm_10.abs() + contrib_jo_10.abs() + contrib_jp_10.abs()
            + contrib_jo_02.abs() + contrib_jp_02.abs() + contrib_jq_02.abs();

        if debug_panels.contains(&jo) {
            println!("{:4} {:4} {:4} {:4} {:10.6} {:10.6} {:10.6} {:10.6}",
                     jo, jp, jm, jq, dso, dsm, dsp, total_contrib);
            if jo == 0 || jo == n - 1 {
                println!("     g1={:.4}, t1={:.4}, g2={:.4}, t2={:.4}", g1, t1, g2, t2);
                println!("     psum_10={:.6}, pdif_10={:.6}", psum_10, pdif_10);
                println!("     psum_02={:.6}, pdif_02={:.6}", psum_02, pdif_02);
            }
        }
    }

    // BIJ[i,j] = -DZDM[j]
    println!("\n=== Final BIJ values (= -DZDM) ===");
    println!("BIJ[0,0] = {:.6} (XFOIL: -0.0221)", -dzdm[0]);
    println!("BIJ[0,1] = {:.6} (XFOIL: 0.0001)", -dzdm[1]);
    println!("BIJ[0,159] = {:.6}", -dzdm[n - 1]);
    println!("BIJ[0,158] = {:.6}", -dzdm[n - 2]);

    // Check self-influence (when i==jo)
    println!("\n=== Self-influence check ===");
    println!("Control point i=0 is on panel jo=0 (from node 0 to node 1)");
    println!("When jo=0, the contribution to DZDM[0] comes from:");
    println!("  - First half-panel contrib_jo_10");
    println!("  - Second half-panel contrib_jo_02");
    println!("These use g1=0, t1=0 when i==jo (coincident point handling)");
}
