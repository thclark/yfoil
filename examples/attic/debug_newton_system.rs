//! Debug VISCAL first iteration

use yfoil::bl::{blsolv, FlowConditions};
use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, SetblConfig, SetblState, build_newton_system};

fn main() {
    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = panel_foil(&geom);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    println!("=== Debug VISCAL First Iteration ===\n");

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let qinv_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();

    // Find stagnation
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    let sst = airfoil.s[stag_idx];
    println!("Stagnation: idx={}, sst={:.6}", stag_idx, sst);
    println!("nbl_upper={}, nbl_lower={}", stag_idx, airfoil.n - stag_idx);

    // Initialize SETBL state
    let mut setbl_state = SetblState::new(&airfoil, stag_idx, sst, &cond);
    setbl_state.init_from_velocity(&airfoil, &qinv_mag);

    println!("\n=== Initial Edge Velocities ===");
    println!("Upper surface (first 5):");
    for i in 0..5.min(setbl_state.nbl_upper) {
        println!("  [{}] xssi={:.6e}, uedg={:.6}", i,
                 setbl_state.upper.xssi[i], setbl_state.upper.uedg[i]);
    }
    println!("Lower surface (first 5):");
    for i in 0..5.min(setbl_state.nbl_lower) {
        println!("  [{}] xssi={:.6e}, uedg={:.6}", i,
                 setbl_state.lower.xssi[i], setbl_state.lower.uedg[i]);
    }

    // Build Newton system
    let config = SetblConfig::default();
    let mut blsolv_input = build_newton_system(&mut setbl_state, &airfoil, &inviscid, &qinv_mag, &config);

    println!("\n=== After build_newton_system ===");
    println!("NSYS = {}", blsolv_input.nsys);
    println!("\nUpper surface stations (first 5):");
    for (i, s) in setbl_state.stations_upper.iter().enumerate().take(5) {
        println!("  [{}] theta={:.6e}, dstar={:.6e}, h={:.3}, hk={:.3}, u={:.4}, cf={:.6e}",
                 i, s.theta, s.dstar, s.h, s.hk, s.u, s.cf);
    }
    println!("\nLower surface stations (first 5):");
    for (i, s) in setbl_state.stations_lower.iter().enumerate().take(5) {
        println!("  [{}] theta={:.6e}, dstar={:.6e}, h={:.3}, hk={:.3}, u={:.4}, cf={:.6e}",
                 i, s.theta, s.dstar, s.h, s.hk, s.u, s.cf);
    }

    // Newton residuals (VDEL before solve)
    println!("\n=== Newton Residuals (first 5 system entries) ===");
    for iv in 0..5.min(blsolv_input.nsys) {
        println!("  [{}] res0={:.6e}, res1={:.6e}, res2={:.6e}",
                 iv, blsolv_input.vdel[iv][0][0], blsolv_input.vdel[iv][1][0], blsolv_input.vdel[iv][2][0]);
    }

    // Solve
    blsolv(&mut blsolv_input);

    println!("\n=== Newton Deltas (after blsolv, first 5 entries) ===");
    for iv in 0..5.min(blsolv_input.nsys) {
        println!("  [{}] dCtau={:.6e}, dTheta={:.6e}, dMass={:.6e}",
                 iv, blsolv_input.vdel[iv][0][0], blsolv_input.vdel[iv][1][0], blsolv_input.vdel[iv][2][0]);
    }
}
