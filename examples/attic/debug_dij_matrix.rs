//! Debug DIJ matrix computation
//!
//! Compare the DIJ matrix influence against expected values

use yfoil::geometry::{panel_foil, read_dat_file};
use yfoil::panel::solve_inviscid;

fn main() {
    println!("=== DIJ Matrix Diagnostic ===\n");

    // Load XFOIL-paneled coordinates
    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = panel_foil(&geom);
    let inviscid = solve_inviscid(&airfoil);

    println!("Airfoil: {} panels", airfoil.n);

    // Check DIJ matrix
    if let Some(dij) = inviscid.get_dij() {
        let n = inviscid.num_panels();

        // Statistics on DIJ
        let mut min_val = f64::MAX;
        let mut max_val = f64::MIN;
        let mut sum_abs = 0.0;

        for i in 0..n {
            for j in 0..n {
                let val = dij[(i, j)];
                if val < min_val {
                    min_val = val;
                }
                if val > max_val {
                    max_val = val;
                }
                sum_abs += val.abs();
            }
        }

        println!("\nDIJ matrix statistics:");
        println!("  Min value: {:.6}", min_val);
        println!("  Max value: {:.6}", max_val);
        println!("  Mean |DIJ|: {:.6}", sum_abs / (n * n) as f64);

        // Check diagonal elements
        println!("\nDiagonal elements (self-influence):");
        for i in [0, n / 4, n / 2, 3 * n / 4, n - 1] {
            println!("  DIJ[{}, {}] = {:.6}", i, i, dij[(i, i)]);
        }

        // Check a row near LE
        let le_idx = airfoil.le_index;
        println!("\nRow {} (near LE) - first few non-diagonal elements:", le_idx);
        let mut vals: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != le_idx)
            .map(|j| (j, dij[(le_idx, j)]))
            .collect();
        vals.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());
        for (j, val) in vals.iter().take(5) {
            println!("  DIJ[{}, {}] = {:.6}", le_idx, j, val);
        }

        // Check what velocity correction we get from a typical mass defect distribution
        println!("\n=== Simulated Mass Defect Test ===");

        // Create a synthetic mass defect distribution (typical δ* ~ 0.002c at mid-chord)
        let mut mass_defect = vec![0.0; n];
        for i in 0..n {
            // Simple model: δ* grows with arc length from LE
            let x = airfoil.x[i];
            if x > 0.1 && x < 0.9 {
                // δ* ~ 0.002 * x^0.5, Ue ~ 1.1
                let dstar = 0.002 * x.sqrt();
                let ue = 1.1;
                mass_defect[i] = ue * dstar;
            }
        }

        // Compute velocity correction
        if let Some(dq) = inviscid.velocity_from_mass_defect(&mass_defect, airfoil.le_index) {
            let mut max_dq: f64 = 0.0;
            let mut max_dq_idx = 0;

            for (i, &dq_i) in dq.iter().enumerate() {
                if dq_i.abs() > max_dq.abs() {
                    max_dq = dq_i;
                    max_dq_idx = i;
                }
            }

            println!("Max velocity correction: {:.6} at panel {} (x={:.3})",
                     max_dq, max_dq_idx, airfoil.x[max_dq_idx]);

            // Print dQ at a few stations
            println!("\nVelocity correction at selected stations:");
            for i in [n / 4, n / 2, 3 * n / 4] {
                println!("  x={:.3}: mass={:.6}, dQ={:.6}",
                         airfoil.x[i], mass_defect[i], dq[i]);
            }

            // Compare to expected: XFOIL shows about +1.4% Ue increase at transition
            // With Ue ~ 1.07, that's dQ ~ 0.015
            // Our mass defect at mid-chord is ~ 0.002 * 0.7 * 1.1 ~ 0.0015
            println!("\nExpected: ~0.015 dQ for typical BL");
            println!("Actual max: {:.6} ({:.2}% of Ue)", max_dq, 100.0 * max_dq / 1.07);
        }
    } else {
        println!("DIJ matrix not computed!");
    }
}
