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

/// Integrate forces using surface velocity (XFOIL-style CLCALC)
///
/// This method uses surface velocity at nodes to compute Cp with
/// Karman-Tsien compressibility correction, then integrates in wind
/// axes using trapezoidal averaging. This matches XFOIL's CLCALC
/// subroutine exactly.
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `velocity` - Surface velocity at nodes (normalized by freestream, length n)
/// * `alpha_rad` - Angle of attack in radians
/// * `mach` - Mach number for compressibility correction
///
/// # Returns
/// Lift, moment, and pressure drag coefficients
pub fn integrate_forces(
    airfoil: &PaneledAirfoil,
    velocity: &[f64],
    alpha_rad: f64,
    mach: f64,
) -> AeroCoefficients {
    let n = airfoil.n;
    let cosa = alpha_rad.cos();
    let sina = alpha_rad.sin();

    // Compressibility parameters (Karman-Tsien)
    let m2 = mach * mach;
    let beta = (1.0 - m2).sqrt().max(0.001);
    let bfac = 0.5 * m2 / (1.0 + beta);

    // Reference point for moments
    let x_ref = airfoil.reference[0];
    let y_ref = airfoil.reference[1];

    let mut cl = 0.0;
    let mut cm = 0.0;
    let mut cdp = 0.0;

    // Unit freestream velocity
    let qinf = 1.0;

    // Compute Cp at first node
    let cginc1 = 1.0 - (velocity[0] / qinf).powi(2);
    let mut cpg1 = if mach < 0.001 {
        cginc1
    } else {
        cginc1 / (beta + bfac * cginc1)
    };

    // Integrate around airfoil using trapezoidal averaging (like XFOIL's CLCALC)
    for i in 0..n {
        let ip1 = (i + 1) % (n + 1);
        let ip1_idx = if ip1 == n { 0 } else { ip1 }; // Handle wraparound for velocity

        // Cp at node ip1
        let cginc2 = 1.0 - (velocity[ip1_idx] / qinf).powi(2);
        let cpg2 = if mach < 0.001 {
            cginc2
        } else {
            cginc2 / (beta + bfac * cginc2)
        };

        // Panel geometry in wind axes
        let i_node = i;
        let ip1_node = if i == n - 1 { 0 } else { i + 1 };

        let dx_body = airfoil.x[ip1_node] - airfoil.x[i_node];
        let dy_body = airfoil.y[ip1_node] - airfoil.y[i_node];

        // Transform to wind axes
        let dx = dx_body * cosa + dy_body * sina;
        let dy = dy_body * cosa - dx_body * sina;

        // Average Cp over panel (trapezoidal)
        let ag = 0.5 * (cpg2 + cpg1);

        // Midpoint in body axes for moment calculation
        let ax_body = 0.5 * (airfoil.x[i_node] + airfoil.x[ip1_node]) - x_ref;
        let ay_body = 0.5 * (airfoil.y[i_node] + airfoil.y[ip1_node]) - y_ref;

        // Transform midpoint to wind axes
        let ax = ax_body * cosa + ay_body * sina;
        let ay = ay_body * cosa - ax_body * sina;

        // Linear Cp variation for moment (XFOIL includes this refinement)
        let dg = cpg2 - cpg1;

        // Force integration
        cl += dx * ag;
        cdp -= dy * ag;
        cm -= dx * (ag * ax + dg * dx / 12.0) + dy * (ag * ay + dg * dy / 12.0);

        // Move to next panel
        cpg1 = cpg2;
    }

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
/// * `velocity` - Vortex strength distribution
///
/// # Returns
/// Lift coefficient from circulation
pub fn cl_from_circulation(airfoil: &PaneledAirfoil, velocity: &[f64]) -> f64 {
    // Total circulation = integral of γ ds around airfoil
    let n = airfoil.n;
    let mut circulation = 0.0;

    for i in 0..n {
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        let ds = ((airfoil.x[ip1] - airfoil.x[i]).powi(2)
            + (airfoil.y[ip1] - airfoil.y[i]).powi(2))
        .sqrt();
        let velocity_avg = 0.5 * (velocity[i] + velocity[ip1]);
        circulation += velocity_avg * ds;
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

        // Get velocity at nodes for alpha=0
        let velocity = solution.velocity_at_nodes(0.0);
        let coeffs = integrate_forces(&airfoil, &velocity, 0.0, 0.0);

        // Symmetric airfoil at α=0 should have CL≈0, CM≈0 (to machine precision)
        assert!(coeffs.cl.abs() < 1e-6, "CL={} should be ~0", coeffs.cl);
        assert!(coeffs.cm.abs() < 1e-6, "CM={} should be ~0", coeffs.cm);
    }

    #[test]
    fn test_positive_lift_at_positive_alpha() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Get velocity at nodes for alpha=5 degrees
        let alpha = 5.0_f64.to_radians();
        let velocity = solution.velocity_at_nodes(alpha);
        let coeffs = integrate_forces(&airfoil, &velocity, alpha, 0.0);

        // Should have positive lift
        assert!(coeffs.cl > 0.3, "CL={} should be positive", coeffs.cl);
    }

    #[test]
    fn test_cambered_airfoil_lift_at_alpha_0() {
        let geom = naca_4digit("4412", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Get velocity at nodes for alpha=0
        let velocity = solution.velocity_at_nodes(0.0);
        let coeffs = integrate_forces(&airfoil, &velocity, 0.0, 0.0);

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
        let velocity = solution.velocity_at_nodes(alpha);

        let coeffs = integrate_forces(&airfoil, &velocity, alpha, 0.0);
        let cl_circ = cl_from_circulation(&airfoil, &velocity);

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

        let vel1 = solution.velocity_at_nodes(alpha1);
        let vel2 = solution.velocity_at_nodes(alpha2);

        let coeffs1 = integrate_forces(&airfoil, &vel1, alpha1, 0.0);
        let coeffs2 = integrate_forces(&airfoil, &vel2, alpha2, 0.0);

        let dcl_dalpha = (coeffs2.cl - coeffs1.cl) / (alpha2 - alpha1);

        // Thin airfoil theory: dCL/dα = 2π ≈ 6.28
        // Panel methods typically give slightly less due to thickness effects
        assert!(
            dcl_dalpha > 5.5 && dcl_dalpha < 7.0,
            "Lift slope {} should be near 2π",
            dcl_dalpha
        );
    }
}
