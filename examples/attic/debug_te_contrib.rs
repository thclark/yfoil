//! Debug TE panel contribution to AIJ

use std::f64::consts::PI;
use yfoil::geometry::{panel_foil, read_dat_file};

const HOPI: f64 = 0.5 / PI;

fn main() {
    println!("=== Debug TE Panel Contribution ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let n = airfoil.n;

    // TE gap
    let te_gap_x = airfoil.x[0] - airfoil.x[n - 1];
    let te_gap_y = airfoil.y[0] - airfoil.y[n - 1];
    let te_gap = (te_gap_x * te_gap_x + te_gap_y * te_gap_y).sqrt();

    println!("N = {}", n);
    println!("TE gap: ({:.6}, {:.6}), magnitude = {:.6}", te_gap_x, te_gap_y, te_gap);

    // Compute bisector direction
    let dx1 = airfoil.x[1] - airfoil.x[0];
    let dy1 = airfoil.y[1] - airfoil.y[0];
    let ds1 = (dx1 * dx1 + dy1 * dy1).sqrt();

    let dxn = airfoil.x[n - 1] - airfoil.x[n - 2];
    let dyn_te = airfoil.y[n - 1] - airfoil.y[n - 2];
    let dsn = (dxn * dxn + dyn_te * dyn_te).sqrt();

    println!("\nUpper panel (0->1): ({:.6}, {:.6}), ds = {:.6}", dx1/ds1, dy1/ds1, ds1);
    println!("Lower panel (n-2->n-1): ({:.6}, {:.6}), ds = {:.6}", dxn/dsn, dyn_te/dsn, dsn);

    // Bisector direction
    let dxs = 0.5 * (-dx1 / ds1 + dxn / dsn);
    let dys = 0.5 * (-dy1 / ds1 + dyn_te / dsn);
    println!("\nBisector direction (DXS, DYS): ({:.6}, {:.6})", dxs, dys);

    // ANTE, ASTE
    let ante = dxs * te_gap_y - dys * te_gap_x;
    let aste = dxs * te_gap_x + dys * te_gap_y;
    println!("ANTE = {:.6}, ASTE = {:.6}", ante, aste);

    // SCS, SDS
    let scs = ante / te_gap;
    let sds = aste / te_gap;
    println!("SCS = {:.6}, SDS = {:.6}", scs, sds);

    // For control point i=0, compute TE panel contribution
    let i = 0;
    let xi = airfoil.x[i];
    let yi = airfoil.y[i];

    println!("\n=== Control point i={}: ({:.6}, {:.6}) ===", i, xi, yi);

    // TE panel geometry
    let jo = n - 1;
    let jp = 0;
    let dso = te_gap;
    let dsio = 1.0 / dso;

    // Vectors from TE panel nodes to control point
    let rx1 = xi - airfoil.x[jo];
    let ry1 = yi - airfoil.y[jo];
    let rx2 = xi - airfoil.x[jp];
    let ry2 = yi - airfoil.y[jp];

    println!("r1 = ({:.6}, {:.6})", rx1, ry1);
    println!("r2 = ({:.6}, {:.6})", rx2, ry2);

    // Unit tangent along TE panel
    let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
    let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;
    println!("TE panel tangent: ({:.6}, {:.6})", sx, sy);

    // Transform to panel-local coordinates
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;

    println!("Local coords: x1={:.6}, x2={:.6}, yy={:.6}", x1, x2, yy);

    // Squared distances
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;
    println!("rs1={:.6}, rs2={:.6}", rs1, rs2);

    // Log and arctan terms
    let apan = airfoil.apanel[jo];
    let (g1, t1) = if i == jo || rs1 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs1.ln(), x1.atan2(yy))
    };
    let (g2, t2) = if i == jp || rs2 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs2.ln(), x2.atan2(yy))
    };

    println!("apan = {:.6} ({:.2} deg)", apan, apan.to_degrees());
    println!("g1={:.6}, t1={:.6}, g2={:.6}, t2={:.6}", g1, t1, g2, t2);
    println!("i==jo: {}, i==jp: {}", i == jo, i == jp);

    // PSIG and PGAM
    let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
    let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
    println!("\nPSIG = {:.6}, PGAM = {:.6}", psig, pgam);

    // TE panel contribution to DZDG
    let contrib_jo = -HOPI * psig * scs * 0.5 + HOPI * pgam * sds * 0.5;
    let contrib_jp = HOPI * psig * scs * 0.5 - HOPI * pgam * sds * 0.5;
    println!("\nTE contribution to DZDG:");
    println!("  DZDG[{}] += {:.6}", jo, contrib_jo);
    println!("  DZDG[{}] += {:.6}", jp, contrib_jp);
}
