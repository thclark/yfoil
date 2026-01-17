use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::forces::{calculate_cp, integrate_forces};
use yfoil::solver::{
    find_stagnation_point, solve_boundary_layer, compute_mass_defect,
    compute_source_velocity_dij, extract_upper_surface, extract_lower_surface,
    ViscalConfig,
};

fn main() {
    let geom = naca_4digit("0012", 160).expect("Failed");
    let airfoil = create_paneled_airfoil(&geom);

    let alpha_deg = 2.0f64;
    let alpha_rad = alpha_deg.to_radians();
    let cond = FlowConditions::new(1e6, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_alpha(alpha_rad);

    // Calculate inviscid CL
    let cp_inv = calculate_cp(&qinv, cond.mach);
    let coeffs_inv = integrate_forces(&airfoil, &cp_inv, alpha_rad);
    println!("Inviscid CL: {:.6}", coeffs_inv.cl);

    // Find stagnation
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    // Get gamma total for wake
    let gamma_total = inviscid.gam_0.iter().sum::<f64>() * alpha_rad.cos()
        + inviscid.gam_90.iter().sum::<f64>() * alpha_rad.sin();

    // First BL solution (with inviscid velocity)
    let bl = solve_boundary_layer(
        &airfoil, &qinv, stag_idx, alpha_rad, gamma_total,
        &cond, &config.newton, &config.wake,
    );

    // Get edge velocities for mass defect
    let (_, _, _, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (_, _, _, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);

    // Compute mass defect
    let mass_defect = compute_mass_defect(
        &airfoil, &bl.upper, &bl.lower, &ue_upper, &ue_lower, stag_idx,
    );

    // Compute source velocity
    let dq_new = compute_source_velocity_dij(&inviscid, &mass_defect, airfoil.le_index);

    println!("\n--- Mass defect (first 10, around LE) ---");
    let le = airfoil.le_index;
    for i in (le.saturating_sub(5))..=(le + 5).min(airfoil.n - 1) {
        println!("i={:3} mass_defect={:+12.6e}", i, mass_defect[i]);
    }

    println!("\n--- Source velocity dq_new (first 10, around LE) ---");
    for i in (le.saturating_sub(5))..=(le + 5).min(airfoil.n - 1) {
        println!("i={:3} dq_new={:+12.6e}", i, dq_new[i]);
    }

    // Compute corrected velocity
    let mut velocity_corrected = qinv.clone();
    for i in 0..airfoil.n {
        // Note: viscal subtracts, so we do: v = qinv - dq_new
        // But let's check both ways
        velocity_corrected[i] = qinv[i] - dq_new[i];
    }

    // Calculate CL with corrected velocity
    let cp_corr = calculate_cp(&velocity_corrected, cond.mach);
    let coeffs_corr = integrate_forces(&airfoil, &cp_corr, alpha_rad);
    println!("\nCL with v = qinv - dq_new: {:.6}", coeffs_corr.cl);

    // Try the other sign
    let mut velocity_corrected2 = qinv.clone();
    for i in 0..airfoil.n {
        velocity_corrected2[i] = qinv[i] + dq_new[i];
    }
    let cp_corr2 = calculate_cp(&velocity_corrected2, cond.mach);
    let coeffs_corr2 = integrate_forces(&airfoil, &cp_corr2, alpha_rad);
    println!("CL with v = qinv + dq_new: {:.6}", coeffs_corr2.cl);

    println!("\nXFOIL viscous CL at alpha=2°: 0.2009");

    // Summary of velocity statistics
    let dq_sum: f64 = dq_new.iter().sum();
    let dq_abs_sum: f64 = dq_new.iter().map(|x| x.abs()).sum();
    let dq_min: f64 = dq_new.iter().copied().fold(f64::INFINITY, f64::min);
    let dq_max: f64 = dq_new.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    println!("\n--- dq_new statistics ---");
    println!("Sum: {:+.6e}", dq_sum);
    println!("Abs Sum: {:.6e}", dq_abs_sum);
    println!("Min: {:+.6e}", dq_min);
    println!("Max: {:+.6e}", dq_max);
}
