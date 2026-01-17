//! Debug first iteration to match XFOIL's viscal iteration log format
//!
//! Outputs values in same format as /tmp/xfoil_viscal_iter.dat for comparison

use yfoil::bl::{march_newton, FlowConditions, NewtonConfig};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{compute_mass_defect, extract_lower_surface, extract_upper_surface, find_stagnation_point};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let inviscid = solve_inviscid(&airfoil);

    let alpha = 0.0_f64.to_radians();
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let newton_config = NewtonConfig::default();

    // Get inviscid velocities using MIDPOINT averaging (matching what viscal now uses)
    // This gives smoother velocities near stagnation
    let qinv = inviscid.velocity_at_alpha(alpha);

    // Find stagnation point
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    println!("yfoil IST (stag panel) = {} (LE index = {})", stag_idx, airfoil.le_index);

    // Extract surfaces using inviscid velocity (no source correction yet - iteration 1)
    let (x_upper, _y_upper, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (x_lower, _y_lower, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);

    println!("NBL(1) equivalent = {} (upper surface stations)", ue_upper.len());
    println!("NBL(2) equivalent = {} (lower surface stations)", ue_lower.len());

    // March BL on each surface
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &newton_config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &newton_config);

    // Output in XFOIL format
    println!("\n=== YFOIL ITERATION 1 ===");
    println!("--- Upper surface BL (IS=1) ---");
    println!("IBL, UEDG, DSTR, THET, MASS, HK");

    // Station 0 is stagnation (IBL=1 in XFOIL)
    println!("{:5} {:17.12e} {:17.12e} {:17.12e} {:17.12e} {:8.4}",
             1, ue_upper[0], 0.0, 0.0, 0.0, 0.0);

    // BL results start at station 1 (IBL=2 in XFOIL)
    for (j, res) in bl_upper.iter().enumerate() {
        let ibl = j + 2; // XFOIL IBL numbering
        let mass = ue_upper[j + 1] * res.dstar;
        println!("{:5} {:17.12e} {:17.12e} {:17.12e} {:17.12e} {:8.4}",
                 ibl, ue_upper[j + 1], res.dstar, res.theta, mass, res.hk);
    }

    // Compare specific stations with XFOIL
    println!("\n=== COMPARISON WITH XFOIL (iteration 1) ===");
    println!("IBL   yfoil_UEDG    xfoil_UEDG    ratio     yfoil_DSTR    xfoil_DSTR    ratio");

    // XFOIL iteration 1 values (from the log):
    let xfoil_uedg = [
        0.0746371488911,  // IBL=2
        0.225868323226,   // IBL=3
        0.372890806171,   // IBL=4
        0.508297598941,   // IBL=5
        0.627694871234,   // IBL=6
        0.729399795164,   // IBL=7
        0.814115973980,   // IBL=8
        0.883790878804,   // IBL=9
        0.940809427252,   // IBL=10
    ];
    let xfoil_dstr = [
        0.714537319852e-4,  // IBL=2
        0.725335634810e-4,  // IBL=3
        0.753823908733e-4,  // IBL=4
        0.798659080465e-4,  // IBL=5
        0.860846051510e-4,  // IBL=6
        0.938142958877e-4,  // IBL=7
        1.03001750037e-4,   // IBL=8
        1.13462466152e-4,   // IBL=9
        1.25142725527e-4,   // IBL=10
    ];

    for i in 0..9.min(bl_upper.len()) {
        let yfoil_ue = ue_upper[i + 1];
        let yfoil_ds = bl_upper[i].dstar;
        let ue_ratio = yfoil_ue / xfoil_uedg[i];
        let ds_ratio = yfoil_ds / xfoil_dstr[i];
        println!("{:3}   {:12.6e}  {:12.6e}  {:6.4}    {:12.6e}  {:12.6e}  {:6.4}",
                 i + 2, yfoil_ue, xfoil_uedg[i], ue_ratio, yfoil_ds, xfoil_dstr[i], ds_ratio);
    }

    // Compute mass defect
    let mass = compute_mass_defect(&airfoil, &bl_upper, &bl_lower, &ue_upper, &ue_lower, stag_idx);

    println!("\n=== MASS DEFECT at key panels ===");
    println!("Panel   yfoil_MASS     xfoil_MASS (approx)   ratio");

    // XFOIL MASS values from iteration 1:
    // IBL=2, panel=80: MASS = 5.333e-6
    // IBL=3, panel=79: MASS = 1.638e-5
    let xfoil_mass = [5.33310283300e-6, 1.63830343610e-5, 2.81094005039e-5];

    for (i, &xm) in xfoil_mass.iter().enumerate() {
        let panel = stag_idx - i - 1;
        if panel < airfoil.n {
            let ym = mass[panel];
            let ratio = ym / xm;
            println!("{:4}    {:12.6e}  {:12.6e}           {:6.4}",
                     panel, ym, xm, ratio);
        }
    }

    // Now compute source velocity correction
    println!("\n=== SOURCE VELOCITY CORRECTION ===");

    if let Some(dij) = inviscid.get_dij() {
        // Compute dq = -DIJ * MASS (unsigned formula)
        let mut dq = vec![0.0; airfoil.n];
        for i in 0..airfoil.n {
            for j in 0..airfoil.n {
                dq[i] -= dij[(i, j)] * mass[j];
            }
        }

        println!("Panel   qinv_abs      dq_source     ue_new (qinv+dq)    xfoil_QVIS");

        // XFOIL QVIS at panels (from iteration 1 log):
        // Panel 80: 0.0746 (near stag)
        // Panel 79: 0.2259
        // Panel 78: 0.3729
        let xfoil_qvis = [
            (80, 0.746371488911e-1),
            (79, 0.225868323226),
            (78, 0.372890806171),
            (77, 0.508297598941),
            (76, 0.627694871234),
            (75, 0.729399795164),
        ];

        for (panel, xfoil_q) in &xfoil_qvis {
            let p = *panel as usize;
            if p < airfoil.n {
                let qinv_abs = qinv[p].abs();
                let ue_new = qinv_abs + dq[p];
                println!("{:4}    {:12.6}  {:+12.6e}  {:12.6}            {:12.6}",
                         p, qinv_abs, dq[p], ue_new, xfoil_q);
            }
        }

        // Also check accumulated effect
        println!("\n=== ACCUMULATED SOURCE CORRECTION MAGNITUDE ===");
        let mut total_dq_upper = 0.0;
        let mut count_upper = 0;
        for i in 0..stag_idx {
            total_dq_upper += dq[i].abs();
            count_upper += 1;
        }
        println!("Sum |dq| on upper surface (panels 0..{}): {:.6e}",
                 stag_idx - 1, total_dq_upper);
        println!("Average |dq| on upper surface: {:.6e}",
                 total_dq_upper / count_upper as f64);

        // Check relative correction at key stations
        println!("\nRelative correction dq/qinv at key stations:");
        for i in [stag_idx - 1, stag_idx - 3, stag_idx - 5, stag_idx - 10, 10, 0] {
            if i < airfoil.n {
                let qinv_abs = qinv[i].abs();
                let rel = if qinv_abs > 1e-6 { dq[i] / qinv_abs } else { 0.0 };
                println!("  Panel {:3}: qinv={:.4}, dq={:+.4e}, dq/qinv={:+.2}%",
                         i, qinv_abs, dq[i], rel * 100.0);
            }
        }
    }
}
