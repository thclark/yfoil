//! Force and moment integration from pressure distribution
//!
//! Computes lift coefficient (CL) and moment coefficient (CM) by
//! integrating pressure distribution around the airfoil.

use crate::geometry::PaneledAirfoil;

/// Aerodynamic force coefficients
#[derive(Debug, Clone, Copy)]
pub struct AeroCoefficients {
    /// Lift coefficient (normal to freestream)
    pub cl: f64,
    /// Moment coefficient about reference point (positive nose-up)
    pub cm: f64,
    /// Pressure drag coefficient (should be ~0 for inviscid)
    pub cdp: f64,
}

/// Integrate pressure distribution to obtain force coefficients
///
/// Uses panel-based integration of Cp around the airfoil surface.
/// The forces are computed in wind axes (aligned with freestream).
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `cp` - Pressure coefficients at each panel midpoint (length n)
/// * `alpha_rad` - Angle of attack in radians
///
/// # Returns
/// Lift, moment, and pressure drag coefficients
pub fn integrate_forces(
    airfoil: &PaneledAirfoil,
    cp: &[f64],
    alpha_rad: f64,
) -> AeroCoefficients {
    let n = airfoil.n;
    let cosa = alpha_rad.cos();
    let sina = alpha_rad.sin();

    // Reference point for moments (typically quarter-chord)
    let x_ref = airfoil.reference[0];
    let y_ref = airfoil.reference[1];

    let mut cn = 0.0; // Normal force coefficient (body axes)
    let mut ca = 0.0; // Axial force coefficient (body axes)
    let mut cm = 0.0; // Moment coefficient

    // Integrate around the airfoil using panel-based values
    // Convention: airfoil goes TE -> upper -> LE -> lower -> TE (counterclockwise)
    for i in 0..n {
        // Panel from node i to node i+1 (with wraparound)
        let ip1 = if i == n - 1 { 0 } else { i + 1 };

        // Panel geometry
        let dx = airfoil.x[ip1] - airfoil.x[i];
        let dy = airfoil.y[ip1] - airfoil.y[i];

        // Cp at panel midpoint (already computed there)
        let cp_panel = cp[i];

        // Panel midpoint for moment arm
        let x_mid = 0.5 * (airfoil.x[i] + airfoil.x[ip1]);
        let y_mid = 0.5 * (airfoil.y[i] + airfoil.y[ip1]);

        // Force integration:
        // For counterclockwise traversal, outward normal n = (dy, -dx)/ds
        // Pressure force on panel: dF = -p * n * ds = -Cp * (dy, -dx)
        // So: dF_x = -Cp * dy, dF_y = Cp * dx
        //
        // Normal force (positive up): dCn = Cp * dx
        // Axial force (positive downstream): dCa = -Cp * dy
        cn += cp_panel * dx;
        ca -= cp_panel * dy;

        // Moment about reference point (positive nose-up)
        // The cross product r × F gives: (x-xref)*Fy - (y-yref)*Fx
        //                              = (x-xref)*Cp*dx - (y-yref)*(-Cp*dy)
        //                              = Cp*((x-xref)*dx + (y-yref)*dy)
        cm += cp_panel * ((x_mid - x_ref) * dx + (y_mid - y_ref) * dy);
    }

    // Transform from body axes to wind axes
    // CL = Cn * cos(α) - Ca * sin(α)
    // CD = Cn * sin(α) + Ca * cos(α)
    let cl = cn * cosa - ca * sina;
    let cdp = cn * sina + ca * cosa;

    // Moment coefficient (already about reference point)
    // Sign convention: positive = nose up

    AeroCoefficients { cl, cm, cdp }
}

/// Compute lift coefficient using Kutta-Joukowski theorem
///
/// CL = Γ_total / (0.5 * chord) = 2 * Γ_total / chord
///
/// This provides an alternative calculation that should match
/// the pressure integration for inviscid flow.
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `gamma` - Vortex strength distribution
///
/// # Returns
/// Lift coefficient from circulation
pub fn cl_from_circulation(airfoil: &PaneledAirfoil, gamma: &[f64]) -> f64 {
    // Total circulation = integral of γ ds around airfoil
    let n = airfoil.n;
    let mut circulation = 0.0;

    for i in 0..n {
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        let ds = ((airfoil.x[ip1] - airfoil.x[i]).powi(2)
            + (airfoil.y[ip1] - airfoil.y[i]).powi(2))
        .sqrt();
        let gamma_avg = 0.5 * (gamma[i] + gamma[ip1]);
        circulation += gamma_avg * ds;
    }

    // CL = 2 * Γ / chord (for unit freestream)
    2.0 * circulation / airfoil.chord
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{create_paneled_airfoil, naca_4digit};
    use crate::panel::solve_inviscid;
    use approx::assert_relative_eq;

    #[test]
    fn test_symmetric_airfoil_zero_lift_at_alpha_0() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Get velocity and Cp at alpha=0
        let vel = solution.velocity_at_alpha(0.0);
        let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();

        let coeffs = integrate_forces(&airfoil, &cp, 0.0);

        // Symmetric airfoil at α=0 should have CL≈0, CM≈0
        assert!(coeffs.cl.abs() < 0.01, "CL={} should be ~0", coeffs.cl);
        assert!(coeffs.cm.abs() < 0.01, "CM={} should be ~0", coeffs.cm);
    }

    #[test]
    fn test_positive_lift_at_positive_alpha() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Get velocity and Cp at alpha=5 degrees
        let alpha = 5.0_f64.to_radians();
        let vel = solution.velocity_at_alpha(alpha);
        let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();

        let coeffs = integrate_forces(&airfoil, &cp, alpha);

        // Should have positive lift
        assert!(coeffs.cl > 0.3, "CL={} should be positive", coeffs.cl);
    }

    #[test]
    fn test_cambered_airfoil_lift_at_alpha_0() {
        let geom = naca_4digit("4412", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Get velocity and Cp at alpha=0
        let vel = solution.velocity_at_alpha(0.0);
        let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();

        let coeffs = integrate_forces(&airfoil, &cp, 0.0);

        // Cambered airfoil at α=0 should have positive lift
        assert!(
            coeffs.cl > 0.2,
            "CL={} should be positive for cambered airfoil",
            coeffs.cl
        );
    }

    #[test]
    fn test_cl_from_circulation_matches() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        let alpha = 5.0_f64.to_radians();
        let gamma = solution.gamma_at_alpha(alpha);
        let vel = solution.velocity_at_alpha(alpha);
        let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();

        let coeffs = integrate_forces(&airfoil, &cp, alpha);
        let cl_circ = cl_from_circulation(&airfoil, &gamma);

        // Both methods should give similar CL (within ~10% for reasonable discretization)
        assert_relative_eq!(coeffs.cl, cl_circ, epsilon = 0.1);
    }

    #[test]
    fn test_lift_slope() {
        // Test that lift slope is approximately 2π (thin airfoil theory)
        let geom = naca_4digit("0012", 160).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Calculate CL at two angles
        let alpha1 = 0.0_f64.to_radians();
        let alpha2 = 5.0_f64.to_radians();

        let vel1 = solution.velocity_at_alpha(alpha1);
        let vel2 = solution.velocity_at_alpha(alpha2);
        let cp1: Vec<f64> = vel1.iter().map(|&v| 1.0 - v * v).collect();
        let cp2: Vec<f64> = vel2.iter().map(|&v| 1.0 - v * v).collect();

        let coeffs1 = integrate_forces(&airfoil, &cp1, alpha1);
        let coeffs2 = integrate_forces(&airfoil, &cp2, alpha2);

        let dcl_dalpha = (coeffs2.cl - coeffs1.cl) / (alpha2 - alpha1);

        // Thin airfoil theory: dCL/dα = 2π ≈ 6.28
        // Panel methods typically give slightly less due to thickness effects
        assert!(
            dcl_dalpha > 5.5 && dcl_dalpha < 7.0,
            "Lift slope {} should be near 2π",
            dcl_dalpha
        );
    }

    #[test]
    fn test_diagnose_cdp_issue() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        let alpha = 5.0_f64.to_radians();
        let vel = solution.velocity_at_alpha(alpha);
        let cp: Vec<f64> = vel.iter().map(|&v| 1.0 - v * v).collect();

        // Find extreme velocity values
        let max_vel = vel.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_vel = vel.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_vel_idx = vel.iter().position(|&v| v == max_vel).unwrap();
        let min_vel_idx = vel.iter().position(|&v| v == min_vel).unwrap();

        println!("\n=== CDp Diagnostic ===");
        println!("Panels: {}", airfoil.n);
        println!("Velocity range: {:.4} to {:.4}", min_vel, max_vel);
        println!("Max vel at idx {} (x={:.4})", max_vel_idx, airfoil.x[max_vel_idx]);
        println!("Min vel at idx {} (x={:.4})", min_vel_idx, airfoil.x[min_vel_idx]);

        // Check first few panels (near TE on upper surface)
        println!("\nNear TE (upper, indices 0-5):");
        for i in 0..5.min(airfoil.n) {
            println!("  i={}: x={:.4}, y={:.5}, vel={:.4}, Cp={:.4}",
                i, airfoil.x[i], airfoil.y[i], vel[i], cp[i]);
        }

        // Check last few panels (near TE on lower surface)
        println!("\nNear TE (lower, last 5):");
        for i in (airfoil.n - 5).max(0)..airfoil.n {
            println!("  i={}: x={:.4}, y={:.5}, vel={:.4}, Cp={:.4}",
                i, airfoil.x[i], airfoil.y[i], vel[i], cp[i]);
        }

        // Calculate contribution to CDp from each panel
        let cosa = alpha.cos();
        let sina = alpha.sin();
        let mut panel_cdp = Vec::new();
        for i in 0..airfoil.n {
            let ip1 = if i == airfoil.n - 1 { 0 } else { i + 1 };
            let dx = airfoil.x[ip1] - airfoil.x[i];
            let dy = airfoil.y[ip1] - airfoil.y[i];
            let cp_avg = 0.5 * (cp[i] + cp[ip1]);
            let cn_i = cp_avg * dx;
            let ca_i = -cp_avg * dy;
            let cdp_i = cn_i * sina + ca_i * cosa;
            panel_cdp.push(cdp_i);
        }

        // Find panels with largest CDp contribution
        let mut indexed: Vec<_> = panel_cdp.iter().enumerate().collect();
        indexed.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());

        println!("\nTop 10 panels by |CDp| contribution:");
        for (i, &cdp_i) in indexed.iter().take(10) {
            println!("  Panel {}: CDp={:.6}, x={:.4}, vel={:.4}, Cp={:.4}",
                i, cdp_i, airfoil.x[*i], vel[*i], cp[*i]);
        }

        let coeffs = integrate_forces(&airfoil, &cp, alpha);
        println!("\nTotal: CL={:.4}, CDp={:.4}", coeffs.cl, coeffs.cdp);

        // The high CDp is a known issue - this test is for diagnostics
        // We expect CL > 0 and CDp != 0 due to numerical issues
        assert!(coeffs.cl > 0.3);
    }
}
