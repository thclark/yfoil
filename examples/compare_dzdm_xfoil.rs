//! Compare yfoil's DZDM array at i=0 against XFOIL's output
//!
//! Run XFOIL debug version first to generate /tmp/xfoil_dzdm_i1.dat

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::f64::consts::PI;
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};

const QOPI: f64 = 0.25 / PI;

fn main() {
    println!("=== DZDM Comparison: yfoil vs XFOIL ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;

    println!("N = {}", n);

    // Compute panel angles - XFOIL's APCALC formula
    let mut apanel = vec![0.0; n];
    for j in 0..n {
        let jp = if j == n - 1 { 0 } else { j + 1 };
        let sx = airfoil.x[jp] - airfoil.x[j];
        let sy = airfoil.y[jp] - airfoil.y[j];

        if sx == 0.0 && sy == 0.0 {
            apanel[j] = (-airfoil.ny[j]).atan2(-airfoil.nx[j]);
        } else if j == n - 1 {
            // TE panel
            if airfoil.sharp_te {
                apanel[j] = std::f64::consts::PI;
            } else {
                apanel[j] = (-sx).atan2(sy) + std::f64::consts::PI;
            }
        } else {
            apanel[j] = sx.atan2(-sy);
        }
    }

    // Compute DZDM at control point i=0 (same as XFOIL I=1)
    let i = 0;
    let xi = airfoil.x[i];
    let yi = airfoil.y[i];

    let mut dzdm = vec![0.0; n];

    // Loop over all panels - same as solver.rs
    // NOTE: XFOIL skips the TE panel (JO=N) for source influence
    for jo in 0..n {
        let jp = if jo == n - 1 { 0 } else { jo + 1 };

        // Skip the TE panel - XFOIL line 245: IF(JO.EQ.N) GO TO 11
        if jo == n - 1 {
            continue;
        }

        let jm = if jo == 0 { jo } else { jo - 1 };
        let jq = if jo == n - 2 {
            n - 1
        } else {
            (jp + 1) % n
        };

        let dso = ((airfoil.x[jp] - airfoil.x[jo]).powi(2)
            + (airfoil.y[jp] - airfoil.y[jo]).powi(2))
        .sqrt();
        if dso < 1e-14 {
            continue;
        }
        let dsio = 1.0 / dso;

        let apan = apanel[jo];

        let rx1 = xi - airfoil.x[jo];
        let ry1 = yi - airfoil.y[jo];
        let rx2 = xi - airfoil.x[jp];
        let ry2 = yi - airfoil.y[jp];

        let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
        let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn = 1.0_f64;

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

        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = if rs0 > 1e-24 { rs0.ln() } else { 0.0 };
        let t0 = (sgn * x0).atan2(sgn * yy);

        // First half-panel
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

        dzdm[jm] += QOPI * (-psum_10 * dsim + pdif_10 * dsim);
        dzdm[jo] += QOPI * (-psum_10 * dsio - pdif_10 * dsio);
        dzdm[jp] += QOPI * (psum_10 * (dsio + dsim) + pdif_10 * (dsio - dsim));

        // Second half-panel
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

        dzdm[jo] += QOPI * (-psum_02 * (dsip + dsio) - pdif_02 * (dsip - dsio));
        dzdm[jp] += QOPI * (psum_02 * dsio - pdif_02 * dsio);
        dzdm[jq] += QOPI * (psum_02 * dsip + pdif_02 * dsip);
    }

    // Save yfoil's DZDM to file
    let mut file = File::create("/tmp/yfoil_dzdm_i0.dat").expect("create file");
    for j in 0..n {
        writeln!(file, "{:5} {:20.10e}", j + 1, dzdm[j]).expect("write");
    }
    println!("Wrote yfoil DZDM to /tmp/yfoil_dzdm_i0.dat");

    // Load XFOIL's DZDM
    let xfoil_file = File::open("/tmp/xfoil_dzdm_i1.dat").expect("Run XFOIL debug first");
    let reader = BufReader::new(xfoil_file);
    let mut xfoil_dzdm = vec![0.0; n];

    for line in reader.lines() {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(j), Ok(val)) = (
                    parts[0].parse::<usize>(),
                    parts[1].parse::<f64>(),
                ) {
                    if j >= 1 && j <= n {
                        xfoil_dzdm[j - 1] = val;
                    }
                }
            }
        }
    }

    // Compare
    println!("\n{:>6} {:>14} {:>14} {:>10}", "j", "yfoil", "XFOIL", "ratio");
    for j in [0, 1, 2, 3, 4, 5, 10, 20, 40, 60, 79, 80, 155, 156, 157, 158, 159] {
        if j < n {
            let yf = dzdm[j];
            let xf = xfoil_dzdm[j];
            let ratio = if xf.abs() > 1e-10 { yf / xf } else { f64::NAN };
            println!("{:6} {:14.6e} {:14.6e} {:10.4}", j, yf, xf, ratio);
        }
    }

    // Statistics
    let mut max_abs_diff = 0.0;
    let mut max_idx = 0;
    let mut sum_diff = 0.0;

    for j in 0..n {
        let diff = (dzdm[j] - xfoil_dzdm[j]).abs();
        sum_diff += diff;
        if diff > max_abs_diff {
            max_abs_diff = diff;
            max_idx = j;
        }
    }

    println!("\nMean absolute diff: {:e}", sum_diff / n as f64);
    println!("Max absolute diff: {:e} at j={}", max_abs_diff, max_idx);
    println!("yfoil DZDM[{}] = {:e}", max_idx, dzdm[max_idx]);
    println!("XFOIL DZDM[{}] = {:e}", max_idx + 1, xfoil_dzdm[max_idx]);
}
