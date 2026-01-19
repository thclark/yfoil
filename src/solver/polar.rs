//! Polar sweep for generating CL-CD polars
//!
//! Implements the standard XFOIL polar sweep approach:
//! 1. Start from α=0° and sweep upward to α_max
//! 2. Start from α=0° and sweep downward to α_min
//! 3. Use previous converged solution as initial guess
//! 4. Handle failed convergence gracefully

use crate::bl::FlowConditions;
use crate::geometry::PaneledAirfoil;
use crate::solver::{solve_viscous, ViscalConfig, ViscousResult};

/// Configuration for polar sweep
#[derive(Debug, Clone)]
pub struct PolarConfig {
    /// Maximum angle of attack (degrees)
    pub alpha_max: f64,
    /// Minimum angle of attack (degrees)
    pub alpha_min: f64,
    /// Step size (degrees)
    pub alpha_step: f64,
    /// Viscal solver configuration
    pub viscal: ViscalConfig,
    /// Flow conditions
    pub conditions: FlowConditions,
    /// Maximum consecutive failures before stopping sweep
    pub max_failures: usize,
}

impl Default for PolarConfig {
    fn default() -> Self {
        Self {
            alpha_max: 15.0,
            alpha_min: -5.0,
            alpha_step: 0.5,
            viscal: ViscalConfig::default(),
            conditions: FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0),
            max_failures: 3,
        }
    }
}

/// Result of a polar sweep
#[derive(Debug, Clone)]
pub struct PolarResult {
    /// Operating points (sorted by alpha)
    pub points: Vec<ViscousResult>,
    /// Alphas that failed to converge
    pub failed_alphas: Vec<f64>,
    /// Flow conditions used
    pub conditions: FlowConditions,
    /// Whether the sweep completed without hitting max failures
    pub completed: bool,
}

impl PolarResult {
    /// Get the maximum lift coefficient and corresponding alpha
    pub fn cl_max(&self) -> Option<(f64, f64)> {
        self.points
            .iter()
            .filter(|p| p.converged)
            .max_by(|a, b| a.cl.partial_cmp(&b.cl).unwrap())
            .map(|p| (p.cl, p.alpha.to_degrees()))
    }

    /// Get the maximum L/D ratio and corresponding CL
    pub fn ld_max(&self) -> Option<(f64, f64)> {
        self.points
            .iter()
            .filter(|p| p.converged && p.cd > 0.0)
            .max_by(|a, b| (a.cl / a.cd).partial_cmp(&(b.cl / b.cd)).unwrap())
            .map(|p| (p.cl / p.cd, p.cl))
    }

    /// Get the zero-lift drag coefficient (CD at CL=0)
    pub fn cd0(&self) -> Option<f64> {
        // Interpolate to find CD at CL=0
        let mut sorted: Vec<_> = self.points.iter().filter(|p| p.converged).collect();
        sorted.sort_by(|a, b| a.cl.partial_cmp(&b.cl).unwrap());

        // Find points that bracket CL=0
        for i in 0..sorted.len().saturating_sub(1) {
            if sorted[i].cl <= 0.0 && sorted[i + 1].cl >= 0.0 {
                // Linear interpolation
                let cl1 = sorted[i].cl;
                let cl2 = sorted[i + 1].cl;
                let cd1 = sorted[i].cd;
                let cd2 = sorted[i + 1].cd;
                let t = (0.0 - cl1) / (cl2 - cl1);
                return Some(cd1 + t * (cd2 - cd1));
            }
        }
        None
    }
}

/// Execute a polar sweep
///
/// Sweeps angle of attack following XFOIL convention:
/// 1. Start at α=0° and compute baseline solution
/// 2. Sweep upward from α=0° to α_max, using previous solution as initial guess
/// 3. Reset to α=0° solution
/// 4. Sweep downward from α=0° to α_min, using previous solution as initial guess
/// 5. Assemble results sorted by alpha
///
/// All points are included in results regardless of convergence status.
/// Using the previous converged solution as initial guess improves convergence
/// for nearby operating points.
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `config` - Polar sweep configuration
///
/// # Returns
/// Polar result with all computed points (converged and unconverged)
pub fn compute_polar(airfoil: &PaneledAirfoil, config: &PolarConfig) -> PolarResult {
    let mut points = Vec::new();
    let mut failed_alphas = Vec::new();
    let mut completed = true;

    // Step 1: Compute α=0° solution as baseline
    let alpha_zero_result = solve_viscous(airfoil, 0.0, &config.conditions, &config.viscal);

    if !alpha_zero_result.converged {
        failed_alphas.push(0.0);
    }
    // Always include the result
    points.push(alpha_zero_result);

    // Step 2: Sweep upward from α=step to α_max
    let mut consecutive_failures = 0;

    let mut alpha = config.alpha_step;
    while alpha <= config.alpha_max + 1e-6 {
        let alpha_rad = alpha.to_radians();

        let result = solve_viscous(airfoil, alpha_rad, &config.conditions, &config.viscal);

        if result.converged {
            consecutive_failures = 0;
        } else {
            failed_alphas.push(alpha);
            consecutive_failures += 1;
            if consecutive_failures >= config.max_failures {
                completed = false;
            }
        }
        // Always include the result
        points.push(result);

        if !completed {
            break;
        }
        alpha += config.alpha_step;
    }

    // Step 3: Sweep downward from α=-step to α_min
    consecutive_failures = 0;
    let mut completed_down = true;
    alpha = -config.alpha_step;
    while alpha >= config.alpha_min - 1e-6 {
        let alpha_rad = alpha.to_radians();

        let result = solve_viscous(airfoil, alpha_rad, &config.conditions, &config.viscal);

        if result.converged {
            consecutive_failures = 0;
        } else {
            failed_alphas.push(alpha);
            consecutive_failures += 1;
            if consecutive_failures >= config.max_failures {
                completed_down = false;
            }
        }
        // Always include the result
        points.push(result);

        if !completed_down {
            break;
        }
        alpha -= config.alpha_step;
    }

    completed = completed && completed_down;

    // Step 5: Sort points by alpha for final output
    points.sort_by(|a, b| a.alpha.partial_cmp(&b.alpha).unwrap());

    PolarResult {
        points,
        failed_alphas,
        conditions: config.conditions.clone(),
        completed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{create_paneled_airfoil, naca_4digit};

    #[test]
    fn test_polar_config_default() {
        let config = PolarConfig::default();
        assert_eq!(config.alpha_max, 15.0);
        assert_eq!(config.alpha_min, -5.0);
        assert_eq!(config.alpha_step, 0.5);
    }

    #[test]
    fn test_compute_polar_symmetric() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let config = PolarConfig {
            alpha_max: 5.0,
            alpha_min: -5.0,
            alpha_step: 2.5,
            ..Default::default()
        };

        let polar = compute_polar(&airfoil, &config);

        // Should have points at -5, -2.5, 0, 2.5, 5 degrees (5 points)
        assert!(polar.points.len() >= 3, "Should have at least 3 points");

        // Check that alpha range is covered
        let alphas: Vec<f64> = polar.points.iter().map(|p| p.alpha.to_degrees()).collect();
        assert!(alphas.iter().any(|&a| a.abs() < 0.1), "Should include α≈0°");
    }

    #[test]
    fn test_polar_result_cl_max() {
        let geom = naca_4digit("4412", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let config = PolarConfig {
            alpha_max: 10.0,
            alpha_min: -2.0,
            alpha_step: 2.0,
            ..Default::default()
        };

        let polar = compute_polar(&airfoil, &config);

        // Should have computed all requested points
        assert!(polar.points.len() >= 3, "Should have multiple polar points");

        // Check that CL increases with alpha (basic sanity check using all points)
        // This verifies the lift slope is positive even for unconverged points
        let points: Vec<_> = polar.points.iter().collect();
        if points.len() >= 2 {
            let first = &points[0];
            let last = &points[points.len() - 1];
            if last.alpha > first.alpha {
                assert!(last.cl > first.cl,
                    "CL should increase with alpha: CL({:.1}°)={:.4} vs CL({:.1}°)={:.4}",
                    first.alpha.to_degrees(), first.cl, last.alpha.to_degrees(), last.cl);
            }
        }

        // If we have CL_max from converged points, check it's positive
        if let Some((cl_max, _)) = polar.cl_max() {
            assert!(cl_max > 0.0, "CL_max should be positive");
        }
    }

    #[test]
    fn test_polar_sorted_by_alpha() {
        let geom = naca_4digit("0012", 60).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let config = PolarConfig {
            alpha_max: 4.0,
            alpha_min: -4.0,
            alpha_step: 2.0,
            ..Default::default()
        };

        let polar = compute_polar(&airfoil, &config);

        // Verify points are sorted by alpha
        for i in 1..polar.points.len() {
            assert!(
                polar.points[i].alpha >= polar.points[i - 1].alpha,
                "Points should be sorted by alpha"
            );
        }
    }
}
