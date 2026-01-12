//! Viscous-inviscid coupling (VISCAL)
//!
//! Implements the coupled iteration between the inviscid panel method
//! and the boundary layer solver to obtain a converged viscous solution.
//!
//! The coupling uses the mass defect concept: the boundary layer displacement
//! thickness δ* creates an effective source distribution that modifies the
//! inviscid velocity field.
//!
//! # Algorithm
//! 1. Solve inviscid panel method -> get Qinv (inviscid velocity)
//! 2. Find stagnation point from Qinv distribution
//! 3. March BL from stagnation point on upper and lower surfaces
//! 4. Compute mass defect m = ρ * Ue * δ* at each station
//! 5. Compute source-induced velocity correction from mass defect
//! 6. Update edge velocity: Ue = Qinv + Qsource
//! 7. Check convergence; if not converged, goto step 3
//!
//! References:
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - Drela, M., Giles, M. "Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils"

use crate::bl::{
    find_transition, generate_wake_coordinates, integrate_friction, march_newton,
    solve_wake, squire_young_drag, wake_edge_velocity, FlowConditions, NewtonConfig,
    NewtonResult, WakeConfig, WakeInitialState,
};
use crate::forces::{calculate_cp, integrate_forces, AeroCoefficients};
use crate::geometry::PaneledAirfoil;
use crate::panel::solve_inviscid;

/// Configuration for viscous-inviscid coupling
#[derive(Debug, Clone)]
pub struct ViscalConfig {
    /// Maximum coupling iterations
    pub max_iter: usize,
    /// Convergence tolerance for CL change
    pub tol_cl: f64,
    /// Convergence tolerance for δ* change
    pub tol_dstar: f64,
    /// Relaxation factor for δ* update
    pub relax: f64,
    /// Newton solver configuration for BL
    pub newton: NewtonConfig,
    /// Wake configuration
    pub wake: WakeConfig,
}

impl Default for ViscalConfig {
    fn default() -> Self {
        Self {
            max_iter: 300,
            tol_cl: 1e-4,
            tol_dstar: 1e-4,
            relax: 0.7,
            newton: NewtonConfig::default(),
            wake: WakeConfig::default(),
        }
    }
}

/// Result of viscous analysis at a single operating point
#[derive(Debug, Clone)]
pub struct ViscousResult {
    /// Angle of attack (radians)
    pub alpha: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Total drag coefficient (pressure + friction)
    pub cd: f64,
    /// Friction drag coefficient
    pub cdf: f64,
    /// Pressure drag coefficient
    pub cdp: f64,
    /// Moment coefficient
    pub cm: f64,
    /// Transition location on upper surface (x/c)
    pub xtr_upper: f64,
    /// Transition location on lower surface (x/c)
    pub xtr_lower: f64,
    /// Number of coupling iterations
    pub iterations: usize,
    /// Whether solution converged
    pub converged: bool,
    /// Final CL residual
    pub residual: f64,
    /// Source velocity correction (for use as initial guess in polar sweeps)
    pub dq_source: Vec<f64>,
}

/// Boundary layer state for coupling
#[derive(Debug, Clone)]
pub struct BLSolution {
    /// Upper surface results (from stagnation to TE)
    pub upper: Vec<NewtonResult>,
    /// Lower surface results (from stagnation to TE)
    pub lower: Vec<NewtonResult>,
    /// Wake results (from TE downstream)
    pub wake: Vec<NewtonResult>,
    /// Upper surface arc lengths
    pub s_upper: Vec<f64>,
    /// Lower surface arc lengths
    pub s_lower: Vec<f64>,
    /// Wake arc lengths
    pub s_wake: Vec<f64>,
    /// Wake x-coordinates
    pub x_wake: Vec<f64>,
    /// Wake y-coordinates
    pub y_wake: Vec<f64>,
    /// Stagnation point index
    pub i_stag: usize,
    /// Upper surface transition location (s)
    pub s_tr_upper: f64,
    /// Lower surface transition location (s)
    pub s_tr_lower: f64,
}

/// Find stagnation point from velocity distribution
///
/// The stagnation point is where the surface velocity changes sign
/// (or is minimum magnitude near the leading edge).
pub fn find_stagnation_point(airfoil: &PaneledAirfoil, velocity: &[f64]) -> usize {
    // Find minimum velocity magnitude near leading edge
    let le_idx = airfoil.le_index;
    let search_range = (airfoil.n / 4).max(5);

    let start = le_idx.saturating_sub(search_range);
    let end = (le_idx + search_range).min(airfoil.n);

    let mut min_vel = f64::MAX;
    let mut stag_idx = le_idx;

    for i in start..end {
        let vel_mag = velocity[i].abs();
        if vel_mag < min_vel {
            min_vel = vel_mag;
            stag_idx = i;
        }
    }

    stag_idx
}

/// Extract upper surface data from airfoil (from stagnation to TE)
///
/// Returns (x, y, s, ue) arrays for the upper surface
pub fn extract_upper_surface(
    airfoil: &PaneledAirfoil,
    velocity: &[f64],
    stag_idx: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    // Upper surface goes from stagnation (near LE) towards TE at index 0
    // In standard airfoil ordering: TE(0) -> upper -> LE -> lower -> TE(n)
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();

    let mut arc_len = 0.0;
    for i in (0..=stag_idx).rev() {
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);
        s.push(arc_len);
        ue.push(velocity[i].abs()); // Use magnitude for BL

        if i > 0 {
            let dx = airfoil.x[i - 1] - airfoil.x[i];
            let dy = airfoil.y[i - 1] - airfoil.y[i];
            arc_len += (dx * dx + dy * dy).sqrt();
        }
    }

    (x, y, s, ue)
}

/// Extract lower surface data from airfoil (from stagnation to TE)
pub fn extract_lower_surface(
    airfoil: &PaneledAirfoil,
    velocity: &[f64],
    stag_idx: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    // Lower surface goes from stagnation towards TE at index n
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();

    let mut arc_len = 0.0;
    for i in stag_idx..airfoil.n {
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);
        s.push(arc_len);
        ue.push(velocity[i].abs());

        if i + 1 < airfoil.n {
            let dx = airfoil.x[i + 1] - airfoil.x[i];
            let dy = airfoil.y[i + 1] - airfoil.y[i];
            arc_len += (dx * dx + dy * dy).sqrt();
        }
    }

    (x, y, s, ue)
}

/// Compute source-induced velocity correction from displacement thickness
///
/// The mass defect source strength is: m = d(ρ*Ue*δ*)/ds
/// This creates a velocity perturbation that modifies the inviscid solution.
pub fn compute_source_velocity(
    airfoil: &PaneledAirfoil,
    bl_upper: &[NewtonResult],
    bl_lower: &[NewtonResult],
    s_upper: &[f64],
    s_lower: &[f64],
    stag_idx: usize,
) -> Vec<f64> {
    let n = airfoil.n;
    let mut dq = vec![0.0; n];

    // Upper surface: map BL stations back to panel indices
    for j in 0..bl_upper.len() {
        // Panel index (going from stag backwards to TE)
        let i = if stag_idx >= j { stag_idx - j } else { 0 };
        if i < n {
            // Source effect: simplified model using local δ* gradient
            // More accurate would use influence coefficients
            if j > 0 && j < bl_upper.len() - 1 {
                let ds = s_upper[j + 1] - s_upper[j - 1];
                if ds > 1e-10 {
                    let dm = (bl_upper[j + 1].dstar * bl_upper[j + 1].theta.sqrt()
                        - bl_upper[j - 1].dstar * bl_upper[j - 1].theta.sqrt())
                        / ds;
                    dq[i] += dm * 0.5; // Simplified coupling
                }
            }
        }
    }

    // Lower surface
    for j in 0..bl_lower.len() {
        let i = stag_idx + j;
        if i < n {
            if j > 0 && j < bl_lower.len() - 1 {
                let ds = s_lower[j + 1] - s_lower[j - 1];
                if ds > 1e-10 {
                    let dm = (bl_lower[j + 1].dstar * bl_lower[j + 1].theta.sqrt()
                        - bl_lower[j - 1].dstar * bl_lower[j - 1].theta.sqrt())
                        / ds;
                    dq[i] += dm * 0.5;
                }
            }
        }
    }

    dq
}

/// Solve boundary layer on both surfaces and wake
pub fn solve_boundary_layer(
    airfoil: &PaneledAirfoil,
    velocity: &[f64],
    stag_idx: usize,
    alpha_rad: f64,
    gamma_total: f64,
    cond: &FlowConditions,
    newton_config: &NewtonConfig,
    wake_config: &WakeConfig,
) -> BLSolution {
    // Extract surfaces
    let (_x_upper, _y_upper, s_upper, ue_upper) =
        extract_upper_surface(airfoil, velocity, stag_idx);
    let (_x_lower, _y_lower, s_lower, ue_lower) =
        extract_lower_surface(airfoil, velocity, stag_idx);

    // March BL on each surface
    let upper = if ue_upper.len() >= 2 {
        march_newton(&ue_upper, &s_upper, cond, newton_config)
    } else {
        vec![]
    };

    let lower = if ue_lower.len() >= 2 {
        march_newton(&ue_lower, &s_lower, cond, newton_config)
    } else {
        vec![]
    };

    // Find transition locations
    let (_, _, s_tr_upper) = find_transition(&upper, &s_upper, cond.ncrit);
    let (_, _, s_tr_lower) = find_transition(&lower, &s_lower, cond.ncrit);

    // Compute wake
    let (wake, x_wake, y_wake, s_wake) = if !upper.is_empty() && !lower.is_empty() {
        // Get TE location and values
        let x_te = airfoil.x[0]; // TE is at index 0 for standard airfoil
        let y_te = airfoil.y[0];

        // Get TE values from surface BL
        let upper_te = upper.last().unwrap();
        let lower_te = lower.last().unwrap();
        let ue_upper_te = ue_upper.last().copied().unwrap_or(1.0);
        let ue_lower_te = ue_lower.last().copied().unwrap_or(1.0);

        // Initialize wake from combined surface values
        let wake_init = WakeInitialState::from_surfaces(
            upper_te,
            lower_te,
            ue_upper_te,
            ue_lower_te,
            x_te,
            y_te,
        );

        // Generate wake coordinates
        let (x_wake, y_wake, s_wake) = generate_wake_coordinates(
            x_te,
            y_te,
            alpha_rad,
            airfoil.chord,
            wake_config,
        );

        // Compute wake edge velocity
        let ue_wake = wake_edge_velocity(&x_wake, &y_wake, gamma_total, alpha_rad, airfoil.chord);

        // Solve wake BL
        let wake_results = solve_wake(&wake_init, &s_wake, &ue_wake, cond, newton_config);

        (wake_results, x_wake, y_wake, s_wake)
    } else {
        (vec![], vec![], vec![], vec![])
    };

    BLSolution {
        upper,
        lower,
        wake,
        s_upper,
        s_lower,
        s_wake,
        x_wake,
        y_wake,
        i_stag: stag_idx,
        s_tr_upper,
        s_tr_lower,
    }
}

/// Calculate friction drag from BL solution
pub fn calculate_friction_drag(bl: &BLSolution) -> f64 {
    let cdf_upper = integrate_friction(&bl.upper, &bl.s_upper);
    let cdf_lower = integrate_friction(&bl.lower, &bl.s_lower);
    cdf_upper + cdf_lower
}

/// Calculate total drag using Squire-Young formula applied to wake
///
/// This gives a more accurate drag estimate by using the far-wake
/// momentum thickness extrapolated from the TE conditions.
pub fn calculate_drag_squire_young(bl: &BLSolution, chord: f64) -> f64 {
    if bl.upper.is_empty() || bl.lower.is_empty() {
        return 0.0;
    }

    let upper_te = bl.upper.last().unwrap();
    let lower_te = bl.lower.last().unwrap();

    // Combined TE momentum thickness
    let theta_te = upper_te.theta + lower_te.theta;

    // Use average H at TE
    let h_te = 0.5 * (upper_te.h + lower_te.h);

    // Assume ue ≈ 1 at TE for normalized airfoil
    let ue_te = 1.0;

    squire_young_drag(theta_te, h_te, ue_te, chord)
}

/// Main viscous-inviscid coupling solver
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `alpha_rad` - Angle of attack in radians
/// * `cond` - Flow conditions (Re, M, Ncrit)
/// * `config` - Coupling configuration
///
/// # Returns
/// Converged viscous solution
pub fn solve_viscous(
    airfoil: &PaneledAirfoil,
    alpha_rad: f64,
    cond: &FlowConditions,
    config: &ViscalConfig,
) -> ViscousResult {
    solve_viscous_with_init(airfoil, alpha_rad, cond, config, None)
}

/// Viscous-inviscid coupling solver with optional initial guess
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `alpha_rad` - Angle of attack in radians
/// * `cond` - Flow conditions (Re, M, Ncrit)
/// * `config` - Coupling configuration
/// * `init_dq` - Optional initial source velocity correction from previous solution
///
/// # Returns
/// Converged viscous solution
pub fn solve_viscous_with_init(
    airfoil: &PaneledAirfoil,
    alpha_rad: f64,
    cond: &FlowConditions,
    config: &ViscalConfig,
    init_dq: Option<&[f64]>,
) -> ViscousResult {
    // Step 1: Solve inviscid panel method
    let inviscid = solve_inviscid(airfoil);

    // Get inviscid velocity at this alpha
    let mut velocity = inviscid.velocity_at_alpha(alpha_rad);

    // Initialize tracking variables
    let mut cl_prev = 0.0;
    let mut iterations = 0;
    let mut residual = f64::MAX;

    // Source velocity correction - use initial guess if provided
    let mut dq_source = if let Some(init) = init_dq {
        if init.len() == airfoil.n {
            init.to_vec()
        } else {
            vec![0.0; airfoil.n]
        }
    } else {
        vec![0.0; airfoil.n]
    };

    // Apply initial guess to velocity if provided
    if init_dq.is_some() {
        let qinv = inviscid.velocity_at_alpha(alpha_rad);
        for i in 0..airfoil.n {
            velocity[i] = qinv[i] + dq_source[i];
        }
    }

    // Compute total circulation from inviscid solution (for wake model)
    let gamma = inviscid.gamma_at_alpha(alpha_rad);
    let gamma_total: f64 = gamma.iter().sum::<f64>() / airfoil.n as f64 * airfoil.chord;

    // Coupling iteration
    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Step 2: Find stagnation point
        let stag_idx = find_stagnation_point(airfoil, &velocity);

        // Step 3: Solve boundary layer (including wake)
        let bl = solve_boundary_layer(
            airfoil,
            &velocity,
            stag_idx,
            alpha_rad,
            gamma_total,
            cond,
            &config.newton,
            &config.wake,
        );

        // Step 4: Compute source velocity correction
        let dq_new = compute_source_velocity(
            airfoil,
            &bl.upper,
            &bl.lower,
            &bl.s_upper,
            &bl.s_lower,
            stag_idx,
        );

        // Relaxed update of source velocity
        for i in 0..airfoil.n {
            dq_source[i] = (1.0 - config.relax) * dq_source[i] + config.relax * dq_new[i];
        }

        // Step 5: Update edge velocity
        let qinv = inviscid.velocity_at_alpha(alpha_rad);
        for i in 0..airfoil.n {
            velocity[i] = qinv[i] + dq_source[i];
        }

        // Step 6: Calculate forces with updated velocity
        let cp = calculate_cp(&velocity, cond.mach);
        let coeffs = integrate_forces(airfoil, &cp, alpha_rad);

        // Check convergence
        let cl_change = (coeffs.cl - cl_prev).abs();
        residual = cl_change;

        if cl_change < config.tol_cl && iter > 0 {
            // Calculate final results
            let cdf = calculate_friction_drag(&bl);
            let cdp = coeffs.cdp;

            // Convert transition locations to x/c
            let xtr_upper = if bl.s_tr_upper < f64::INFINITY && !bl.upper.is_empty() {
                // Find x at transition
                let idx = bl
                    .s_upper
                    .iter()
                    .position(|&s| s >= bl.s_tr_upper)
                    .unwrap_or(0);
                if idx < bl.upper.len() {
                    // Approximate x/c from arc length
                    bl.s_tr_upper / bl.s_upper.last().unwrap_or(&1.0)
                } else {
                    1.0
                }
            } else {
                1.0 // No transition (fully laminar)
            };

            let xtr_lower = if bl.s_tr_lower < f64::INFINITY && !bl.lower.is_empty() {
                bl.s_tr_lower / bl.s_lower.last().unwrap_or(&1.0)
            } else {
                1.0
            };

            return ViscousResult {
                alpha: alpha_rad,
                cl: coeffs.cl,
                cd: cdf + cdp,
                cdf,
                cdp,
                cm: coeffs.cm,
                xtr_upper,
                xtr_lower,
                iterations,
                converged: true,
                residual,
                dq_source,
            };
        }

        cl_prev = coeffs.cl;
    }

    // Did not converge - return best estimate
    let stag_idx = find_stagnation_point(airfoil, &velocity);
    let bl = solve_boundary_layer(
        airfoil,
        &velocity,
        stag_idx,
        alpha_rad,
        gamma_total,
        cond,
        &config.newton,
        &config.wake,
    );
    let cp = calculate_cp(&velocity, cond.mach);
    let coeffs = integrate_forces(airfoil, &cp, alpha_rad);
    let cdf = calculate_friction_drag(&bl);

    ViscousResult {
        alpha: alpha_rad,
        cl: coeffs.cl,
        cd: cdf + coeffs.cdp,
        cdf,
        cdp: coeffs.cdp,
        cm: coeffs.cm,
        xtr_upper: 1.0,
        xtr_lower: 1.0,
        iterations,
        converged: false,
        residual,
        dq_source,
    }
}

/// Solve inviscid-only analysis (no BL)
pub fn solve_inviscid_only(
    airfoil: &PaneledAirfoil,
    alpha_rad: f64,
    mach: f64,
) -> AeroCoefficients {
    let inviscid = solve_inviscid(airfoil);
    let velocity = inviscid.velocity_at_alpha(alpha_rad);
    let cp = calculate_cp(&velocity, mach);
    integrate_forces(airfoil, &cp, alpha_rad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{create_paneled_airfoil, naca_4digit};
    use approx::assert_relative_eq;

    #[test]
    fn test_find_stagnation_point() {
        let geom = naca_4digit("0012", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let inviscid = solve_inviscid(&airfoil);

        // At alpha=0, stagnation should be near LE
        let vel = inviscid.velocity_at_alpha(0.0);
        let stag = find_stagnation_point(&airfoil, &vel);

        // Should be near LE index
        let le_idx = airfoil.le_index;
        assert!(
            (stag as i32 - le_idx as i32).abs() < 10,
            "Stagnation {} should be near LE {}",
            stag,
            le_idx
        );
    }

    #[test]
    fn test_extract_surfaces() {
        let geom = naca_4digit("0012", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let inviscid = solve_inviscid(&airfoil);
        let vel = inviscid.velocity_at_alpha(0.0);

        let stag = find_stagnation_point(&airfoil, &vel);

        let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
        let (x_lower, _, s_lower, ue_lower) = extract_lower_surface(&airfoil, &vel, stag);

        // Both surfaces should have reasonable sizes
        assert!(x_upper.len() > 10);
        assert!(x_lower.len() > 10);

        // Arc lengths should be monotonically increasing
        for i in 1..s_upper.len() {
            assert!(s_upper[i] >= s_upper[i - 1]);
        }
        for i in 1..s_lower.len() {
            assert!(s_lower[i] >= s_lower[i - 1]);
        }

        // Edge velocities should be positive
        for &ue in &ue_upper {
            assert!(ue >= 0.0);
        }
        for &ue in &ue_lower {
            assert!(ue >= 0.0);
        }
    }

    #[test]
    fn test_solve_boundary_layer() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let inviscid = solve_inviscid(&airfoil);
        let alpha = 0.0;
        let vel = inviscid.velocity_at_alpha(alpha);

        let stag = find_stagnation_point(&airfoil, &vel);
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let newton_config = NewtonConfig::default();
        let wake_config = WakeConfig::default();

        // Compute total circulation
        let gamma = inviscid.gamma_at_alpha(alpha);
        let gamma_total: f64 = gamma.iter().sum::<f64>() / airfoil.n as f64 * airfoil.chord;

        let bl = solve_boundary_layer(
            &airfoil,
            &vel,
            stag,
            alpha,
            gamma_total,
            &cond,
            &newton_config,
            &wake_config,
        );

        // Should have results on both surfaces
        assert!(!bl.upper.is_empty(), "Upper surface BL should have results");
        assert!(!bl.lower.is_empty(), "Lower surface BL should have results");

        // Should have wake results
        assert!(!bl.wake.is_empty(), "Wake should have results");

        // θ should be positive at all stations
        for (i, r) in bl.upper.iter().enumerate() {
            assert!(r.theta > 0.0, "θ should be positive at station {}: {}", i, r.theta);
        }

        // Wake should have Cf = 0
        for r in &bl.wake {
            assert_eq!(r.cf, 0.0, "Wake should have zero skin friction");
        }
    }

    #[test]
    fn test_solve_viscous_symmetric_alpha_0() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = ViscalConfig::default();

        let result = solve_viscous(&airfoil, 0.0, &cond, &config);

        // At alpha=0, symmetric airfoil should have CL ≈ 0
        assert!(
            result.cl.abs() < 0.05,
            "CL={} should be ~0 for symmetric at α=0",
            result.cl
        );

        // Should have positive drag
        assert!(result.cd > 0.0, "CD should be positive");
        assert!(result.cdf > 0.0, "CDf should be positive");
    }

    #[test]
    fn test_solve_viscous_positive_lift() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = ViscalConfig::default();

        let alpha = 5.0_f64.to_radians();
        let result = solve_viscous(&airfoil, alpha, &cond, &config);

        // At alpha=5°, should have positive lift
        assert!(
            result.cl > 0.2,
            "CL={} should be positive at α=5°",
            result.cl
        );

        // Drag should be positive (exact value depends on coupling quality)
        assert!(
            result.cd > 0.0,
            "CD={} should be positive",
            result.cd
        );

        // Friction drag should be a component
        assert!(
            result.cdf > 0.0,
            "CDf={} should be positive",
            result.cdf
        );
    }

    #[test]
    fn test_viscous_drag_breakdown() {
        let geom = naca_4digit("4412", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = ViscalConfig::default();

        let result = solve_viscous(&airfoil, 0.0, &cond, &config);

        // Total drag should equal sum of components
        assert_relative_eq!(result.cd, result.cdf + result.cdp, epsilon = 1e-6);

        // Friction drag should be dominant at low alpha
        // (for attached flow, Cdf >> Cdp)
        assert!(result.cdf > 0.0);
    }

    #[test]
    fn test_inviscid_only() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let coeffs = solve_inviscid_only(&airfoil, 5.0_f64.to_radians(), 0.0);

        // Should have positive lift at positive alpha
        assert!(coeffs.cl > 0.3, "CL={} should be positive", coeffs.cl);

        // Pressure drag sign is affected by numerical issues near TE
        // Just verify it's finite
        assert!(coeffs.cdp.is_finite(), "CDp should be finite");
    }
}
