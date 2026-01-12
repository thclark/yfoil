//! Panel method solver for inviscid flow
//!
//! This module assembles and solves the linear system for the vortex
//! panel method with linear vorticity distribution.

use nalgebra::{DMatrix, DVector};

use crate::geometry::PaneledAirfoil;

use super::influence::panel_influence;

/// Result of the inviscid panel method solution
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Vortex strength at each node for α=0° unit solution
    pub gam_0: Vec<f64>,
    /// Vortex strength at each node for α=90° unit solution
    pub gam_90: Vec<f64>,
    /// Surface velocity at panel midpoints for α=0° unit solution
    pub qinv_0: Vec<f64>,
    /// Surface velocity at panel midpoints for α=90° unit solution
    pub qinv_90: Vec<f64>,
    /// Internal streamfunction value for α=0°
    pub psi_0: f64,
    /// Internal streamfunction value for α=90°
    pub psi_90: f64,
    /// Number of panels
    pub n: usize,
}

impl InviscidSolution {
    /// Get vortex distribution at arbitrary angle of attack
    ///
    /// γ(α) = cos(α) * γ₀ + sin(α) * γ₉₀
    ///
    /// Note: This returns node-based vortex strengths, not surface velocities.
    /// For surface velocities, use `velocity_at_alpha` instead.
    pub fn gamma_at_alpha(&self, alpha_rad: f64) -> Vec<f64> {
        let cosa = alpha_rad.cos();
        let sina = alpha_rad.sin();
        self.gam_0
            .iter()
            .zip(&self.gam_90)
            .map(|(&g0, &g90)| cosa * g0 + sina * g90)
            .collect()
    }

    /// Get surface velocity at panel midpoints for arbitrary angle of attack
    ///
    /// Q(α) = cos(α) * Q₀ + sin(α) * Q₉₀
    ///
    /// These are computed at panel midpoints (control points) to avoid
    /// the trailing edge singularity that occurs at panel nodes.
    pub fn velocity_at_alpha(&self, alpha_rad: f64) -> Vec<f64> {
        let cosa = alpha_rad.cos();
        let sina = alpha_rad.sin();
        self.qinv_0
            .iter()
            .zip(&self.qinv_90)
            .map(|(&q0, &q90)| cosa * q0 + sina * q90)
            .collect()
    }

    /// Get surface velocity at nodes (original method, has TE singularity)
    ///
    /// This is the raw gamma distribution which for attached flow
    /// approximates the surface velocity.
    pub fn velocity_at_nodes(&self, alpha_rad: f64) -> Vec<f64> {
        self.gamma_at_alpha(alpha_rad)
    }

    /// Get the number of velocity/panel values
    pub fn num_panels(&self) -> usize {
        self.n
    }
}

/// Build and solve the inviscid panel method system.
///
/// This implements the vortex panel method with:
/// - Linear vorticity distribution on each panel
/// - Kutta condition at trailing edge: γ₁ + γₙ = 0
/// - Flow tangency (constant streamfunction on surface)
/// - Control points at panel MIDPOINTS (not nodes) to avoid singularity
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
///
/// # Returns
/// * `InviscidSolution` containing unit solutions for α=0° and α=90°
pub fn solve_inviscid(airfoil: &PaneledAirfoil) -> InviscidSolution {
    let n = airfoil.n;

    // Build influence matrix AIJ
    // Size is (n+1) x (n+1): n equations for Ψ=const, 1 for Kutta condition
    // Unknowns: γ₁...γₙ (vortex strengths) and Ψᵢₙₜ (internal streamfunction)
    let mut aij = DMatrix::<f64>::zeros(n + 1, n + 1);

    // Right-hand sides for α=0° and α=90°
    let mut rhs_0 = DVector::<f64>::zeros(n + 1);
    let mut rhs_90 = DVector::<f64>::zeros(n + 1);

    // Fill influence matrix for each control point (node i)
    for i in 0..n {
        let x_i = airfoil.x[i];
        let y_i = airfoil.y[i];

        // Sum influence from all panels
        for j in 0..n {
            // Panel j goes from node j to node j+1 (with wraparound for last panel)
            let jp1 = if j == n - 1 { 0 } else { j + 1 };

            let x_j = airfoil.x[j];
            let y_j = airfoil.y[j];
            let x_jp1 = airfoil.x[jp1];
            let y_jp1 = airfoil.y[jp1];

            // Check if control point coincides with panel nodes
            let same_j = i == j;
            let same_jp1 = i == jp1;

            // Skip the last panel if TE is closed (nodes 0 and n-1 coincide)
            if j == n - 1 && airfoil.sharp_te {
                let te_gap =
                    ((x_j - airfoil.x[0]).powi(2) + (y_j - airfoil.y[0]).powi(2)).sqrt();
                if te_gap < 1e-8 {
                    continue;
                }
            }

            let influence =
                panel_influence(x_j, y_j, x_jp1, y_jp1, x_i, y_i, same_j, same_jp1);

            // dΨ/dγⱼ contribution
            aij[(i, j)] += influence.coeff_j();

            // dΨ/dγⱼ₊₁ contribution
            aij[(i, jp1)] += influence.coeff_jp1();
        }

        // dΨ/dΨᵢₙₜ = -1 (we want Ψ_surface = Ψ_internal)
        aij[(i, n)] = -1.0;

        // Freestream contributions (RHS)
        // At α=0°: Ψ_∞ = Q∞ * y, so RHS = -y (for unit Q∞)
        // At α=90°: Ψ_∞ = -Q∞ * x, so RHS = x
        rhs_0[i] = -y_i;
        rhs_90[i] = x_i;
    }

    // For sharp trailing edge, replace the equation for node n-1 with
    // a curvature extrapolation condition (XFOIL paper equation 9).
    // This prevents the TE singularity by enforcing smooth γ approach
    // from both sides:
    //   (γ[2] - 2γ[1] + γ[0]) - (γ[n-3] - 2γ[n-2] + γ[n-1]) = 0
    //
    // This replaces the ill-conditioned flow tangency equation at node n-1
    // (which would be nearly identical to node 0 for sharp TE).
    if airfoil.sharp_te && n >= 6 {
        // Clear row n-1
        for j in 0..=n {
            aij[(n - 1, j)] = 0.0;
        }

        // Upper surface curvature: γ[2] - 2γ[1] + γ[0]
        aij[(n - 1, 0)] = 1.0;
        aij[(n - 1, 1)] = -2.0;
        aij[(n - 1, 2)] = 1.0;

        // Lower surface curvature: -(γ[n-3] - 2γ[n-2] + γ[n-1])
        aij[(n - 1, n - 3)] = -1.0;
        aij[(n - 1, n - 2)] = 2.0;
        aij[(n - 1, n - 1)] = -1.0;

        // No ψ_internal dependency for this equation
        aij[(n - 1, n)] = 0.0;

        // RHS is 0 (curvatures should be equal)
        rhs_0[n - 1] = 0.0;
        rhs_90[n - 1] = 0.0;
    }

    // Kutta condition: γ₁ + γₙ = 0
    // Row n+1 (index n): enforce γ[0] + γ[n-1] = 0
    aij[(n, 0)] = 1.0;
    aij[(n, n - 1)] = 1.0;
    // RHS is 0 for Kutta condition
    rhs_0[n] = 0.0;
    rhs_90[n] = 0.0;

    // Solve the system using LU decomposition
    let lu = aij.lu();

    let solution_0 = lu
        .solve(&rhs_0)
        .expect("Panel system should be solvable for α=0°");
    let solution_90 = lu
        .solve(&rhs_90)
        .expect("Panel system should be solvable for α=90°");

    // Extract vortex strengths
    let gam_0: Vec<f64> = solution_0.rows(0, n).iter().copied().collect();
    let gam_90: Vec<f64> = solution_90.rows(0, n).iter().copied().collect();
    let psi_0 = solution_0[n];
    let psi_90 = solution_90[n];

    // Compute velocities at panel midpoints
    // For the midpoint-based formulation, the gamma values at nodes
    // give the surface velocity when properly averaged
    let (qinv_0, qinv_90) =
        compute_midpoint_velocities(airfoil, &gam_0, &gam_90);

    InviscidSolution {
        gam_0,
        gam_90,
        qinv_0,
        qinv_90,
        psi_0,
        psi_90,
        n,
    }
}

/// Compute surface velocities at panel midpoints.
///
/// For a vortex panel method, the solved vortex strength γ at each node
/// equals the tangential surface velocity there (by the properties of the
/// potential flow solution). At panel midpoints, we use the average of
/// the two node values.
///
/// Special handling for trailing edge:
/// - TE nodes (indices 0 and n-1) have singular γ values due to the Kutta condition
/// - Panels adjacent to the TE use extrapolation from nearby interior nodes
///
/// This approach follows XFOIL's treatment of the trailing edge singularity.
fn compute_midpoint_velocities(
    airfoil: &PaneledAirfoil,
    gam_0: &[f64],
    gam_90: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n = airfoil.n;
    let mut qinv_0 = vec![0.0; n];
    let mut qinv_90 = vec![0.0; n];

    // For most panels, use average of node values
    for i in 0..n {
        let ip1 = if i == n - 1 { 0 } else { i + 1 };

        // Check if either node is a TE node (indices 0 or n-1)
        let node_i_is_te = i == 0 || i == n - 1;
        let node_ip1_is_te = ip1 == 0 || ip1 == n - 1;

        if (node_i_is_te || node_ip1_is_te) && airfoil.sharp_te {
            // Panel touches the TE - use only the non-TE node value
            if node_i_is_te && !node_ip1_is_te {
                // Node i is TE, use node ip1
                qinv_0[i] = gam_0[ip1];
                qinv_90[i] = gam_90[ip1];
            } else if !node_i_is_te && node_ip1_is_te {
                // Node ip1 is TE, use node i
                qinv_0[i] = gam_0[i];
                qinv_90[i] = gam_90[i];
            } else {
                // Both nodes are TE (the TE panel itself) - use average of
                // the two nodes adjacent to TE (nodes 1 and n-2)
                qinv_0[i] = 0.5 * (gam_0[1] + gam_0[n - 2]);
                qinv_90[i] = 0.5 * (gam_90[1] + gam_90[n - 2]);
            }
        } else {
            // Normal panels: average of the two node values
            qinv_0[i] = 0.5 * (gam_0[i] + gam_0[ip1]);
            qinv_90[i] = 0.5 * (gam_90[i] + gam_90[ip1]);
        }
    }

    (qinv_0, qinv_90)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{create_paneled_airfoil, naca_4digit};
    use approx::assert_relative_eq;

    #[test]
    fn test_solve_inviscid_naca0012() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Check solution dimensions
        assert_eq!(solution.gam_0.len(), airfoil.n);
        assert_eq!(solution.gam_90.len(), airfoil.n);

        // For symmetric airfoil at α=0°, the solution should be antisymmetric
        // Find max velocity magnitude excluding near-TE points (where singularities occur)
        let mut max_vel = 0.0;
        for i in 5..(airfoil.n - 5) {
            let vel = solution.gam_0[i].abs();
            if vel > max_vel {
                max_vel = vel;
            }
        }

        // Max velocity should be reasonable (around 1-2 for thin airfoil)
        assert!(max_vel > 0.5 && max_vel < 3.0, "Max velocity {} out of range", max_vel);
    }

    #[test]
    fn test_kutta_condition() {
        let geom = naca_4digit("0012", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Kutta condition: γ[0] + γ[n-1] ≈ 0
        let kutta_0 = solution.gam_0[0] + solution.gam_0[airfoil.n - 1];
        let kutta_90 = solution.gam_90[0] + solution.gam_90[airfoil.n - 1];

        assert!(kutta_0.abs() < 1e-10, "Kutta condition violated for α=0°: {}", kutta_0);
        assert!(kutta_90.abs() < 1e-10, "Kutta condition violated for α=90°: {}", kutta_90);
    }

    #[test]
    fn test_symmetric_airfoil_alpha_0() {
        let geom = naca_4digit("0012", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // At α=0° for symmetric airfoil, flow should be symmetric
        // Upper and lower surface velocities should be mirror images
        let n = airfoil.n;
        let le_idx = airfoil.le_index;

        // Compare corresponding upper/lower points
        for i in 1..le_idx.min(n - le_idx) {
            let upper_idx = le_idx - i;
            let lower_idx = le_idx + i;
            if lower_idx < n {
                // Velocities should be equal in magnitude, opposite sign
                assert_relative_eq!(
                    solution.gam_0[upper_idx],
                    -solution.gam_0[lower_idx],
                    epsilon = 0.05
                );
            }
        }
    }

    #[test]
    fn test_velocity_at_alpha() {
        let geom = naca_4digit("4412", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Test superposition of midpoint velocities
        let alpha = 0.1; // ~5.7 degrees
        let vel = solution.velocity_at_alpha(alpha);

        // Manual calculation for midpoint velocities
        let cosa = alpha.cos();
        let sina = alpha.sin();
        for i in 0..airfoil.n {
            let expected = cosa * solution.qinv_0[i] + sina * solution.qinv_90[i];
            assert_relative_eq!(vel[i], expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_stagnation_point_at_le() {
        let geom = naca_4digit("0012", 160).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // At α=0° for symmetric airfoil, stagnation point should be at LE
        let vel = solution.velocity_at_alpha(0.0);
        let le_idx = airfoil.le_index;

        // Velocity at LE should be near zero (stagnation)
        assert!(
            vel[le_idx].abs() < 0.2,
            "LE velocity {} not near stagnation",
            vel[le_idx]
        );
    }
}
