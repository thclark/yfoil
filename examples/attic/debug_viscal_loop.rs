//! Debug VISCAL loop - trace multiple iterations

use std::fs;
use yfoil::bl::{blsolv, FlowConditions};
use yfoil::geometry::{panel_foil, Geometry};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{
    find_stagnation_point, SetblConfig, SetblState, build_newton_system,
    apply_newton_update_with_ue,
};

fn main() {
    // Load geometry from JSON file (same as used by XFOIL for comparison)
    let json_path = "/tmp/naca0012.json";
    let json_str = fs::read_to_string(json_path)
        .expect("Failed to read geometry JSON - run: yfoil geometry naca 0012 -n 160 -o /tmp/naca0012.json");
    let geom: Geometry = serde_json::from_str(&json_str).expect("Failed to parse JSON");
    let airfoil = panel_foil(&geom);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    println!("=== Debug VISCAL Loop (Combined Update) ===\n");

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let qinv_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();

    // Find stagnation
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    let sst = airfoil.s[stag_idx];
    println!("Stagnation: idx={}, sst={:.6}", stag_idx, sst);

    // Initialize SETBL state
    let mut setbl_state = SetblState::new(&airfoil, stag_idx, sst, &cond);
    setbl_state.init_from_velocity(&airfoil, &qinv_mag);

    // Initialize edge velocities
    let mut ue_mag = qinv_mag.clone();
    let config = SetblConfig::default();

    // Run 15 iterations with detailed logging
    for iter in 0..15 {
        println!("\n{}", "=".repeat(60));
        println!("ITERATION {}", iter + 1);
        println!("{}", "=".repeat(60));

        // Update edge velocities in SETBL state
        for (ibl, &ipan) in setbl_state.ipan_upper.iter().enumerate() {
            if ibl < setbl_state.nbl_upper {
                setbl_state.upper.uedg[ibl] = ue_mag[ipan];
            }
        }
        for (ibl, &ipan) in setbl_state.ipan_lower.iter().enumerate() {
            if ibl < setbl_state.nbl_lower {
                setbl_state.lower.uedg[ibl] = ue_mag[ipan];
            }
        }

        println!("\nEdge velocities (station 0-2):");
        for i in 0..3.min(setbl_state.nbl_upper) {
            println!("  Upper[{}]: uedg={:.6}", i, setbl_state.upper.uedg[i]);
        }

        // Build Newton system (pass qinv_mag for USAV computation)
        let mut blsolv_input = build_newton_system(
            &mut setbl_state,
            &airfoil,
            &inviscid,
            &qinv_mag,
            &config,
        );

        println!("\nAFTER build_newton_system:");
        println!("  NSYS = {}", blsolv_input.nsys);

        // VA diagonal values (critical for BLSOLV pivots)
        println!("\n  VA diagonal (first 3 stations):");
        for iv in 0..3.min(blsolv_input.nsys) {
            println!("    [{}] VA[0][0]={:.6e}, VA[1][1]={:.6e}",
                     iv, blsolv_input.va[iv][0][0], blsolv_input.va[iv][1][1]);
        }

        // VM diagonal values
        println!("\n  VM diagonal (first 3 stations):");
        for iv in 0..3.min(blsolv_input.nsys) {
            println!("    [{}] VM[iv][iv]=[{:.6e}, {:.6e}, {:.6e}]",
                     iv, blsolv_input.vm[iv][iv][0], blsolv_input.vm[iv][iv][1], blsolv_input.vm[iv][iv][2]);
        }

        // Newton residuals
        println!("\n  Residuals (first 3):");
        for iv in 0..3.min(blsolv_input.nsys) {
            println!("    [{}] res0={:.6e}, res1={:.6e}, res2={:.6e}",
                     iv, blsolv_input.vdel[iv][0][0], blsolv_input.vdel[iv][1][0], blsolv_input.vdel[iv][2][0]);
        }

        // Check station 76 (lower surface near TE)
        let iv76 = 76.min(blsolv_input.nsys - 1);
        println!("\n  Station {} (lower surface near TE):", iv76);
        println!("    VA[0][0]={:.6e}, VA[1][1]={:.6e}",
                 blsolv_input.va[iv76][0][0], blsolv_input.va[iv76][1][1]);
        println!("    VM[iv][iv]=[{:.6e}, {:.6e}, {:.6e}]",
                 blsolv_input.vm[iv76][iv76][0], blsolv_input.vm[iv76][iv76][1], blsolv_input.vm[iv76][iv76][2]);
        println!("    res=[{:.6e}, {:.6e}, {:.6e}]",
                 blsolv_input.vdel[iv76][0][0], blsolv_input.vdel[iv76][1][0], blsolv_input.vdel[iv76][2][0]);

        // Solve
        blsolv(&mut blsolv_input);

        // Newton deltas
        println!("\n  Deltas (first 3):");
        for iv in 0..3.min(blsolv_input.nsys) {
            println!("    [{}] dCtau={:.6e}, dTheta={:.6e}, dMass={:.6e}",
                     iv, blsolv_input.vdel[iv][0][0], blsolv_input.vdel[iv][1][0], blsolv_input.vdel[iv][2][0]);
        }

        // Apply combined Newton update (BL variables + edge velocity)
        let (result, new_ue_mag) = apply_newton_update_with_ue(
            &mut setbl_state,
            &blsolv_input.vdel,
            &inviscid,
            &qinv_mag,
            &config,
        );

        println!("\nAFTER apply_newton_update_with_ue:");
        println!("  RMSBL = {:.6e}", result.rmsbl);

        println!("\n  Upper stations (0-2):");
        for i in 0..3.min(setbl_state.nbl_upper) {
            let s = &setbl_state.stations_upper[i];
            println!("    [{}] theta={:.6e}, dstar={:.6e}, h={:.3}, cf={:.6e}, u={:.4}",
                     i, s.theta, s.dstar, s.h, s.cf, s.u);
        }

        // Show edge velocity changes
        let old_ue1 = ue_mag[setbl_state.ipan_upper[1]];
        let new_ue1 = new_ue_mag[setbl_state.ipan_upper[1]];
        println!("\n  Edge velocity change:");
        println!("    Upper[1] panel {}: ue {:.6} -> {:.6}",
                 setbl_state.ipan_upper[1], old_ue1, new_ue1);

        // Update for next iteration
        ue_mag = new_ue_mag;

        // Check mass defect values being used
        println!("\n  Mass defect (upper stations 0-2):");
        for i in 0..3.min(setbl_state.nbl_upper) {
            let station = &setbl_state.stations_upper[i];
            println!("    [{}] mass={:.6e}", i, station.dstar * station.u);
        }

        // Check for convergence
        if result.rmsbl < 1e-4 {
            println!("\n*** CONVERGED at iteration {} ***", iter + 1);
            break;
        }

        if result.rmsbl > 100.0 {
            println!("\n*** DIVERGED at iteration {} ***", iter + 1);
            break;
        }
    }

    // Final result
    println!("\n{}", "=".repeat(60));
    println!("FINAL STATE");
    println!("{}", "=".repeat(60));

    println!("\nUpper surface (first 10 stations):");
    for (i, s) in setbl_state.stations_upper.iter().enumerate().take(10) {
        let cf_status = if s.cf.is_nan() { "NaN" } else if s.cf < 0.0 { "NEG" } else { "OK" };
        println!("  [{}] theta={:.6e}, cf={:.6e} ({}), h={:.3}", i, s.theta, s.cf, cf_status, s.h);
    }

    // Calculate CDf
    let cdf_upper: f64 = (1..setbl_state.nbl_upper)
        .map(|i| {
            let ipan = setbl_state.ipan_upper[i];
            let ipan_prev = setbl_state.ipan_upper[i - 1];
            let dx = (airfoil.x[ipan] - airfoil.x[ipan_prev]).abs();
            let cf_avg = 0.5 * (setbl_state.stations_upper[i - 1].cf + setbl_state.stations_upper[i].cf);
            cf_avg * dx
        })
        .sum();

    println!("\nCalculated CDf (upper) = {:.6e}", cdf_upper);
}
