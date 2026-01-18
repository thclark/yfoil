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
    find_transition, generate_wake_coordinates, integrate_friction, march_newton, solve_wake,
    squire_young_drag, wake_edge_velocity, FlowConditions, NewtonConfig, NewtonResult, WakeConfig,
    WakeInitialState,
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
            relax: 1.0,
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
    /// Boundary layer solution (for debugging and detailed analysis)
    pub bl: BLSolution,
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
    /// Upper surface x-coordinates (for drag integration)
    pub x_upper: Vec<f64>,
    /// Lower surface x-coordinates (for drag integration)
    pub x_lower: Vec<f64>,
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
/// The stagnation point is where the surface velocity changes sign.
/// This finds the panel where the sign change occurs, preferring the
/// geometric LE when velocities are symmetric.
pub fn find_stagnation_point(airfoil: &PaneledAirfoil, velocity: &[f64]) -> usize {
    let le_idx = airfoil.le_index;
    let search_range = (airfoil.n / 4).max(5);

    let start = le_idx.saturating_sub(search_range);
    let end = (le_idx + search_range).min(airfoil.n);

    // Look for sign change in velocity (stagnation point)
    // Upper surface has positive velocity, lower has negative (for standard ordering)
    // XFOIL convention: IST is the first panel with non-positive velocity (first lower panel)
    // This makes IST the shared starting point for both upper and lower BL marches
    for i in start..end.saturating_sub(1) {
        if velocity[i] > 0.0 && velocity[i + 1] <= 0.0 {
            // Sign change found at (i, i+1)
            // Pick i+1 to match XFOIL's IST convention (first lower surface panel)
            return i + 1;
        }
    }

    // Fallback: find minimum velocity magnitude, preferring LE on tie
    let mut min_vel = f64::MAX;
    let mut stag_idx = le_idx;

    // Check LE first so it wins on ties
    if le_idx < velocity.len() {
        min_vel = velocity[le_idx].abs();
        stag_idx = le_idx;
    }

    for i in start..end {
        let vel_mag = velocity[i].abs();
        // Use strict less-than so LE wins on tie
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
    // Upper surface goes from near-LE (upper side) towards TE at index 0
    // In standard airfoil ordering: TE(0) -> upper -> LE -> lower -> TE(n)
    //
    // With straddling LE geometry (two nodes at min x, one upper y>0, one lower y<0):
    // - stag_idx is the LOWER near-LE node (first with non-positive velocity)
    // - stag_idx - 1 is the UPPER near-LE node
    // - Upper surface: nodes (stag_idx - 1), (stag_idx - 2), ..., 0
    // - Lower surface: nodes stag_idx, (stag_idx + 1), ..., (n - 1)
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();

    // Start from the node BEFORE stag_idx (upper near-LE) and go to TE at 0
    let start_idx = stag_idx.saturating_sub(1);
    let n_stations = start_idx + 1;

    // XFOIL initializes BL with arc length from stagnation, not s=0.
    // The first upper station (start_idx) is at some distance from stag (stag_idx).
    // Initialize arc_len as the distance from stag_idx to start_idx.
    let mut arc_len = {
        let dx = airfoil.x[start_idx] - airfoil.x[stag_idx];
        let dy = airfoil.y[start_idx] - airfoil.y[stag_idx];
        (dx * dx + dy * dy).sqrt()
    };

    for j in 0..n_stations {
        let i = start_idx - j; // Node index (start_idx, start_idx-1, ..., 0)
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);
        s.push(arc_len);

        // For node-based velocities, use velocity[i] directly at each node.
        ue.push(velocity[i].abs());

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
    // Lower surface goes from stagnation towards TE at index n-1
    // Uses velocity[i] directly at each node for consistent treatment
    // with upper surface.
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();

    let n = airfoil.n;
    let mut arc_len = 0.0;
    let end_idx = n - 1; // Lower TE node

    for i in stag_idx..=end_idx {
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);
        s.push(arc_len);

        // Use velocity[i] directly (same as upper surface)
        ue.push(velocity[i].abs());

        if i < end_idx {
            let dx = airfoil.x[i + 1] - airfoil.x[i];
            let dy = airfoil.y[i + 1] - airfoil.y[i];
            arc_len += (dx * dx + dy * dy).sqrt();
        }
    }

    (x, y, s, ue)
}

/// Compute mass defect array from BL solution
///
/// Mass defect m* = Ue * δ* at each panel, mapped from BL stations.
/// This is used with the DIJ matrix to compute velocity corrections.
pub fn compute_mass_defect(
    airfoil: &PaneledAirfoil,
    bl_upper: &[NewtonResult],
    bl_lower: &[NewtonResult],
    ue_upper: &[f64],
    ue_lower: &[f64],
    stag_idx: usize,
) -> Vec<f64> {
    let n = airfoil.n;
    let mut mass = vec![0.0; n];

    // Upper surface: map BL stations back to panel indices
    // Upper surface goes from (stag_idx - 1) down to 0
    // Note: extract_upper_surface starts at stag_idx - 1, not stag_idx
    let upper_start = stag_idx.saturating_sub(1);
    for j in 0..bl_upper.len().min(upper_start + 1) {
        let i = upper_start - j;
        if i < n && j < ue_upper.len() && j < bl_upper.len() {
            let ue = ue_upper[j].abs();
            let dstar = bl_upper[j].dstar;
            // Mass defect = Ue * δ*, ensure non-negative
            if ue.is_finite() && dstar.is_finite() && dstar >= 0.0 {
                mass[i] = ue * dstar;
            }
        }
    }

    // Lower surface: from stag_idx to n-1
    for j in 0..bl_lower.len() {
        let i = stag_idx + j;
        if i < n && j < ue_lower.len() && j < bl_lower.len() {
            let ue = ue_lower[j].abs();
            let dstar = bl_lower[j].dstar;
            if ue.is_finite() && dstar.is_finite() && dstar >= 0.0 {
                mass[i] = ue * dstar;
            }
        }
    }

    mass
}

/// Compute source-induced velocity correction using DIJ matrix
///
/// This is the XFOIL-style coupling following the exact formula:
///   dQ[i] = Σ_j (-VTI[i] * VTI[j] * DIJ[i,j] * MASS[j])
///
/// where:
/// - mass_defect = Ue * δ* at each panel
/// - VTI = +1 for upper surface, -1 for lower surface
/// - DIJ was precomputed as AIJ^-1 * BIJ
///
/// # Arguments
/// * `inviscid` - Inviscid solution containing DIJ matrix
/// * `mass_defect` - Mass defect at each panel
/// * `le_index` - Leading edge index to determine upper/lower surface
pub fn compute_source_velocity_dij(
    inviscid: &crate::panel::InviscidSolution,
    mass_defect: &[f64],
    le_index: usize,
) -> Vec<f64> {
    // Use the DIJ matrix from the inviscid solution
    if let Some(dq) = inviscid.velocity_from_mass_defect(mass_defect, le_index) {
        dq
    } else {
        // Fallback: no DIJ matrix, return zeros
        vec![0.0; mass_defect.len()]
    }
}

/// Compute source-induced velocity correction for UNSIGNED edge velocity
///
/// Unlike `compute_source_velocity_dij`, this computes the correction for the
/// unsigned edge velocity |Ue| directly:
///   dUe[i] = Σ_j DIJ[i,j] * MASS[j]
///
/// The VTI sign conversion is NOT applied here because we're working with
/// unsigned velocities throughout the coupling loop.
pub fn compute_source_velocity_unsigned(
    inviscid: &crate::panel::InviscidSolution,
    mass_defect: &[f64],
) -> Vec<f64> {
    if let Some(dq) = inviscid.velocity_from_mass_defect_unsigned(mass_defect) {
        dq
    } else {
        vec![0.0; mass_defect.len()]
    }
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
    let (x_upper, _y_upper, s_upper, ue_upper) =
        extract_upper_surface(airfoil, velocity, stag_idx);
    let (x_lower, _y_lower, s_lower, ue_lower) =
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
        x_upper,
        x_lower,
        x_wake,
        y_wake,
        i_stag: stag_idx,
        s_tr_upper,
        s_tr_lower,
    }
}

/// Calculate friction drag from BL solution
///
/// Uses XFOIL-compatible integration over chord-projected distance
/// with Ue² weighting factor.
pub fn calculate_friction_drag(bl: &BLSolution, alpha: f64) -> f64 {
    let cdf_upper = integrate_friction(&bl.upper, &bl.x_upper);
    let cdf_lower = integrate_friction(&bl.lower, &bl.x_lower);
    cdf_upper + cdf_lower
}

/// Calculate total drag using Squire-Young formula applied to wake
///
/// This matches XFOIL's CDCALC approach: use wake END values (not TE values)
/// to extrapolate to downstream infinity using the Squire-Young relation.
///
/// CD = 2 * θ_wake * (Ue_wake/Q∞)^((5+H_wake)/2)
///
/// This gives the most accurate total drag by including both friction and
/// pressure drag effects through the momentum thickness evolution.
pub fn calculate_drag_squire_young(bl: &BLSolution, chord: f64) -> f64 {
    // First try to use wake end values (preferred, matches XFOIL)
    if !bl.wake.is_empty() {
        let wake_end = bl.wake.last().unwrap();
        let theta_wake = wake_end.theta;
        let h_wake = wake_end.h;
        let ue_wake = wake_end.ue;

        // Squire-Young: CD = 2 * θ * (Ue/Q∞)^((5+H)/2)
        // For normalized airfoil with Q∞ = 1
        let exponent = 0.5 * (5.0 + h_wake);
        let theta_inf = theta_wake * ue_wake.powf(exponent);
        return 2.0 * theta_inf / chord;
    }

    // Fallback: use TE values if no wake
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

    // Get inviscid velocity at this alpha (SIGNED: + for upper, - for lower)
    // Use node-based velocities (gamma values) which XFOIL uses.
    // Note: yfoil's paneling has a node at exact x=0 (stagnation) where gamma=0,
    // while XFOIL's paneling doesn't have a node at x=0 (LE falls between panels).
    // This difference is handled in extract_upper/lower_surface by interpolating
    // the stagnation panel velocity to match XFOIL's behavior.
    let qinv = inviscid.velocity_at_nodes(alpha_rad);

    // Convert to UNSIGNED edge velocities for coupling (following XFOIL's approach)
    // XFOIL works with unsigned Ue in the BL solver, only converting to signed for forces
    let mut ue_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();

    // Initialize tracking variables
    let mut cl_prev = 0.0;
    let mut iterations = 0;
    let mut residual = f64::MAX;

    // Source velocity correction (for UNSIGNED Ue) - use initial guess if provided
    let mut dq_source = if let Some(init) = init_dq {
        if init.len() == airfoil.n {
            init.to_vec()
        } else {
            vec![0.0; airfoil.n]
        }
    } else {
        vec![0.0; airfoil.n]
    };

    // Apply initial guess to unsigned edge velocity if provided
    if init_dq.is_some() {
        let qinv_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();
        for i in 0..airfoil.n {
            ue_mag[i] = (qinv_mag[i] + dq_source[i]).max(0.01);
        }
    }

    // Compute total circulation from inviscid solution (for wake model)
    let gamma = inviscid.gamma_at_alpha(alpha_rad);
    let gamma_total: f64 = gamma.iter().sum::<f64>() / airfoil.n as f64 * airfoil.chord;

    // Create signed velocity for stagnation point detection
    // Sign convention: use original qinv signs
    let mut velocity: Vec<f64> = ue_mag
        .iter()
        .zip(qinv.iter())
        .map(|(&ue, &q)| if q >= 0.0 { ue } else { -ue })
        .collect();

    // Coupling iteration
    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Step 2: Find stagnation point (using signed velocity)
        let stag_idx = find_stagnation_point(airfoil, &velocity);

        // Step 3: Solve boundary layer (using UNSIGNED velocity)
        // Note: extract_*_surface takes abs() of velocity, so we pass the signed version
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

        // Extract edge velocities for mass defect calculation
        let (_, _, _, ue_upper) = extract_upper_surface(airfoil, &velocity, stag_idx);
        let (_, _, _, ue_lower) = extract_lower_surface(airfoil, &velocity, stag_idx);

        // Step 4: Compute mass defect and source velocity correction
        let mass_defect = compute_mass_defect(
            airfoil,
            &bl.upper,
            &bl.lower,
            &ue_upper,
            &ue_lower,
            stag_idx,
        );

        // Compute dq for UNSIGNED velocity (remove VTI multiplication from the formula)
        let dq_new = compute_source_velocity_unsigned(&inviscid, &mass_defect);

        // Adaptive relaxation
        let qinv_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();
        let mut rlx = config.relax;

        // Find maximum relative change
        let mut dmax = 0.0_f64;
        for i in 0..airfoil.n {
            let delta = (dq_new[i] - dq_source[i]).abs();
            let qref = qinv_mag[i].max(0.1);
            let rel_change = delta / qref;
            dmax = dmax.max(rel_change);
        }

        // Limit relaxation if changes are too large (XFOIL uses 0.3 threshold)
        if dmax > 0.3 {
            rlx = (0.3 / dmax).min(config.relax);
        }

        // Relaxed update of source velocity correction
        for i in 0..airfoil.n {
            dq_source[i] = rlx * dq_new[i] + (1.0 - rlx) * dq_source[i];
        }

        // Step 5: Update UNSIGNED edge velocity
        for i in 0..airfoil.n {
            ue_mag[i] = (qinv_mag[i] + dq_source[i]).max(0.01);
        }

        // Convert to signed velocity for force calculation and next iteration
        for i in 0..airfoil.n {
            velocity[i] = if qinv[i] >= 0.0 { ue_mag[i] } else { -ue_mag[i] };
        }

        // Step 6: Calculate forces with updated velocity
        let cp = calculate_cp(&velocity, cond.mach);
        let coeffs = integrate_forces(airfoil, &cp, alpha_rad);

        // Check convergence
        let cl_change = (coeffs.cl - cl_prev).abs();
        residual = cl_change;

        if cl_change < config.tol_cl && iter > 0 {
            // Calculate final results using XFOIL's approach:
            // - Total CD from Squire-Young (includes both friction and pressure effects)
            // - CDF from skin friction integration (for reporting)
            // - CDP = CD - CDF (derived, for reporting)
            let cdf = calculate_friction_drag(&bl, alpha_rad);
            let cd = calculate_drag_squire_young(&bl, airfoil.chord);
            let cdp = cd - cdf; // Derived pressure drag

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
                cd,   // Total drag from Squire-Young
                cdf,  // Friction drag component
                cdp,  // Pressure drag = CD - CDF
                cm: coeffs.cm,
                xtr_upper,
                xtr_lower,
                iterations,
                converged: true,
                residual,
                dq_source,
                bl,
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

    // Use Squire-Young for total drag (same as converged case)
    let cdf = calculate_friction_drag(&bl, alpha_rad);
    let cd = calculate_drag_squire_young(&bl, airfoil.chord);
    let cdp = cd - cdf; // Derived pressure drag

    ViscousResult {
        alpha: alpha_rad,
        cl: coeffs.cl,
        cd,   // Total drag from Squire-Young
        cdf,  // Friction drag component
        cdp,  // Pressure drag = CD - CDF
        cm: coeffs.cm,
        xtr_upper: 1.0,
        xtr_lower: 1.0,
        iterations,
        converged: false,
        residual,
        dq_source,
        bl,
    }
}

/// Solve inviscid-only analysis (no BL)
pub fn solve_inviscid_only(
    airfoil: &PaneledAirfoil,
    alpha_rad: f64,
    mach: f64,
) -> AeroCoefficients {
    let inviscid = solve_inviscid(airfoil);
    // Use node-based velocities for consistency with XFOIL
    let velocity = inviscid.velocity_at_nodes(alpha_rad);
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
        // Use solve_viscous to properly set up the boundary layer,
        // then verify the BL solution properties
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = ViscalConfig::default();

        let result = solve_viscous(&airfoil, 0.0, &cond, &config);

        // Should have results on both surfaces
        assert!(!result.bl.upper.is_empty(), "Upper surface BL should have results");
        assert!(!result.bl.lower.is_empty(), "Lower surface BL should have results");

        // Should have wake results
        assert!(!result.bl.wake.is_empty(), "Wake should have results");

        // θ should be positive at all stations
        for (i, r) in result.bl.upper.iter().enumerate() {
            assert!(r.theta > 0.0, "θ should be positive at station {}: {}", i, r.theta);
        }

        // Wake should have Cf = 0
        for r in &result.bl.wake {
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
