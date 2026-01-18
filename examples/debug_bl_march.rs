//! Debug script to compare BL march station-by-station against XFOIL
//!
//! Outputs detailed BL state at each station for comparison with instrumented XFOIL.

use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, extract_lower_surface, find_stagnation_point};

fn main() {
    // Create NACA 0012 airfoil with 160 panels (XFOIL default)
    let geom = naca_4digit("0012", 160).expect("Failed to create airfoil");
    let airfoil = create_paneled_airfoil(&geom);

    println!("=== BL March Debug for NACA 0012 at alpha=0, Re=1e6 ===");
    println!("N panels: {}", airfoil.n);

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let alpha = 0.0_f64;
    let qinv = inviscid.velocity_at_nodes(alpha);

    // Find stagnation point
    let velocity: Vec<f64> = qinv.iter().map(|&q| if q >= 0.0 { q } else { -q }).collect();
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    println!("\nStagnation point: index {}, x={:.6}, y={:.6}",
             stag_idx, airfoil.x[stag_idx], airfoil.y[stag_idx]);

    // Extract surfaces
    let (x_upper, y_upper, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (x_lower, y_lower, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);

    println!("\nUpper surface: {} stations", ue_upper.len());
    println!("Lower surface: {} stations", ue_lower.len());

    // Print edge velocity comparison
    println!("\n=== Edge Velocity at Key Stations ===");
    println!("{:>6} {:>10} {:>14} {:>14}", "Idx", "x/c", "Ue_upper", "Ue_lower");

    let n_print = 5.min(ue_upper.len()).min(ue_lower.len());
    for i in 0..n_print {
        println!("{:>6} {:>10.6} {:>14.8e} {:>14.8e}",
                 i, x_upper[i], ue_upper[i], ue_lower[i]);
    }

    // Print first few and last few
    println!("...");
    let start = ue_upper.len().saturating_sub(5);
    for i in start..ue_upper.len() {
        if i < ue_lower.len() {
            println!("{:>6} {:>10.6} {:>14.8e} {:>14.8e}",
                     i, x_upper[i], ue_upper[i], ue_lower[i]);
        }
    }

    // Flow conditions
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();

    // March BL on upper surface
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &config);

    println!("\n=== BL State at First Few Stations (Upper) ===");
    println!("{:>6} {:>10} {:>14} {:>14} {:>14} {:>10} {:>10}",
             "Idx", "s", "theta", "dstar", "Ue", "Hk", "N");
    for (i, r) in bl_upper.iter().take(10).enumerate() {
        println!("{:>6} {:>10.6} {:>14.8e} {:>14.8e} {:>14.8e} {:>10.4} {:>10.4}",
                 i, s_upper[i], r.theta, r.dstar, r.ue, r.hk, r.n_amp);
    }

    println!("\n=== BL State at TE (Last 5 Stations) ===");
    println!("\nUpper surface:");
    let start = bl_upper.len().saturating_sub(5);
    for i in start..bl_upper.len() {
        let r = &bl_upper[i];
        println!("{:>6} {:>10.6} {:>14.8e} {:>14.8e} {:>14.8e} {:>10.4} {:>10.4}",
                 i, s_upper[i], r.theta, r.dstar, r.ue, r.hk, r.n_amp);
    }

    println!("\nLower surface:");
    let start = bl_lower.len().saturating_sub(5);
    for i in start..bl_lower.len() {
        let r = &bl_lower[i];
        println!("{:>6} {:>10.6} {:>14.8e} {:>14.8e} {:>14.8e} {:>10.4} {:>10.4}",
                 i, s_lower[i], r.theta, r.dstar, r.ue, r.hk, r.n_amp);
    }

    // Final TE values
    let upper_te = bl_upper.last().unwrap();
    let lower_te = bl_lower.last().unwrap();

    println!("\n=== Summary ===");
    println!("Upper TE: theta={:.8e}, dstar={:.8e}, Ue={:.8e}",
             upper_te.theta, upper_te.dstar, upper_te.ue);
    println!("Lower TE: theta={:.8e}, dstar={:.8e}, Ue={:.8e}",
             lower_te.theta, lower_te.dstar, lower_te.ue);
    println!("Combined theta (TTE): {:.8e}", upper_te.theta + lower_te.theta);
    println!("Combined dstar (DTE): {:.8e}", upper_te.dstar + lower_te.dstar);

    // XFOIL expected values
    println!("\n=== XFOIL Expected (from fixture) ===");
    println!("XFOIL TTE: 3.9961520077e-3");
    println!("XFOIL DTE: 8.8492753897e-3 (includes ANTE)");
    println!("XFOIL ANTE: 2.4955581871e-3");
    println!("XFOIL theta per surface: ~2.0e-3");

    // Calculate error
    let xfoil_tte = 3.9961520077e-3;
    let yfoil_tte = upper_te.theta + lower_te.theta;
    let error = (yfoil_tte - xfoil_tte).abs() / xfoil_tte * 100.0;
    println!("\nYFoil TTE error: {:.1}%", error);
}
