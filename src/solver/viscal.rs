//! Viscous-inviscid coupling (VISCAL)
//!
//! Implements XFOIL's coupled viscous-inviscid iteration using the
//! SETBL + BLSOLV + UPDATE algorithm.
//!
//! # Algorithm (XFOIL-compatible)
//!
//! Each iteration:
//! 1. **SETBL** - Build linearized Newton system from current BL state
//!    - March BL on both surfaces using current edge velocities
//!    - Compute local Jacobians (VA, VB matrices)
//!    - Build global coupling matrix (VM) from DIJ influence
//! 2. **BLSOLV** - Solve the 3N×3N block-tridiagonal system
//! 3. **UPDATE** - Apply Newton deltas with adaptive per-variable relaxation
//!    - Compute RMSBL convergence metric
//!
//! References:
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - Drela, M., Giles, M. "Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils"

use crate::bl::{blsolv, integrate_friction, squire_young_drag, system::BLStationState, FlowConditions, NewtonResult};
use crate::forces::{integrate_forces, AeroCoefficients};
use crate::geometry::PaneledAirfoil;
use crate::panel::solve_inviscid;
use crate::solver::setbl::{apply_newton_update_with_ue, build_newton_system, SetblConfig, SetblState};

/// Configuration for viscous-inviscid coupling
#[derive(Debug, Clone)]
pub struct ViscalConfig {
    /// Maximum coupling iterations
    pub max_iter: usize,
    /// Convergence tolerance for RMSBL (RMS of normalized variable changes)
    pub tol_rmsbl: f64,
    /// BLSOLV acceleration parameter (XFOIL's VACCEL)
    pub vaccel: f64,
    /// SETBL configuration
    pub setbl: SetblConfig,
}

impl Default for ViscalConfig {
    fn default() -> Self {
        Self {
            max_iter: 300,
            tol_rmsbl: 1e-4,
            vaccel: 0.01,
            setbl: SetblConfig::default(),
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

/// Compute interpolated stagnation arc length (SST) matching XFOIL's STFIND
///
/// XFOIL computes an interpolated arc length at the exact stagnation point
/// by linearly interpolating between the two panels where the velocity changes sign.
/// This is used as the origin for BL arc length computation on both surfaces.
///
/// # Arguments
/// * `airfoil` - Paneled airfoil with arc lengths
/// * `velocity` - Velocity distribution at panel nodes
/// * `stag_idx` - Stagnation point index (first panel with non-positive velocity)
///
/// # Returns
/// Interpolated arc length at stagnation point
pub fn compute_stagnation_arc_length(airfoil: &PaneledAirfoil, velocity: &[f64], stag_idx: usize) -> f64 {
    // XFOIL convention:
    // - stag_idx (IST in XFOIL, 1-based) is where GAM(I) >= 0 and GAM(I+1) < 0
    // - In YFoil, stag_idx is the first panel with non-positive velocity (IST+1 in XFOIL terms)
    // - So stag_idx-1 has positive velocity (upper side), stag_idx has non-positive (lower side)

    if stag_idx == 0 || stag_idx >= airfoil.n {
        return airfoil.s[stag_idx.min(airfoil.n - 1)];
    }

    // XFOIL uses panels IST and IST+1, which in our 0-based indexing is stag_idx-1 and stag_idx
    let i = stag_idx - 1; // Panel with positive velocity (XFOIL's IST)
    let gam_i = velocity[i];
    let gam_i1 = velocity[stag_idx];
    let dgam = gam_i1 - gam_i;
    let ds = airfoil.s[stag_idx] - airfoil.s[i];

    // XFOIL formula: evaluate to minimize roundoff for very small GAM values
    // See STFIND in xpanel.f
    let mut sst = if gam_i < -gam_i1 {
        airfoil.s[i] - ds * (gam_i / dgam)
    } else {
        airfoil.s[stag_idx] - ds * (gam_i1 / dgam)
    };

    // Tweak stagnation point if it falls right on a node (very unlikely)
    if sst <= airfoil.s[i] {
        sst = airfoil.s[i] + 1.0e-7;
    }
    if sst >= airfoil.s[stag_idx] {
        sst = airfoil.s[stag_idx] - 1.0e-7;
    }

    sst
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
/// Returns (x, y, s, ue) arrays for the upper surface.
/// Arc lengths are computed as XSSI = SST - S(i), matching XFOIL's XICALC.
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

    // Compute XEPS: minimum arc length near stagnation (XFOIL uses 1e-7 * total)
    let total_arc = airfoil.s[airfoil.n - 1] - airfoil.s[0];
    let xeps = 1e-7 * total_arc;

    // Compute SST: interpolated stagnation arc length
    let sst = compute_stagnation_arc_length(airfoil, velocity, stag_idx);

    // XFOIL XICALC: XSSI(IBL,IS) = MAX( SST - S(I) , XEPS )
    // Upper surface goes backward from stagnation, so S(I) < SST and XSSI > 0
    for j in 0..n_stations {
        let i = start_idx - j; // Node index (start_idx, start_idx-1, ..., 0)
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);

        // Arc length from stagnation point (SST - S[i])
        let arc_len = (sst - airfoil.s[i]).max(xeps);
        s.push(arc_len);

        // For node-based velocities, use velocity[i] directly at each node.
        ue.push(velocity[i].abs());
    }

    (x, y, s, ue)
}

/// Extract lower surface data from airfoil (from stagnation to TE)
///
/// Returns (x, y, s, ue) arrays for the lower surface.
/// Arc lengths are computed as XSSI = S(i) - SST, matching XFOIL's XICALC.
pub fn extract_lower_surface(
    airfoil: &PaneledAirfoil,
    velocity: &[f64],
    stag_idx: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    // Lower surface goes from stagnation towards TE at index n-1
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();

    // Compute XEPS: minimum arc length near stagnation (XFOIL uses 1e-7 * total)
    let total_arc = airfoil.s[airfoil.n - 1] - airfoil.s[0];
    let xeps = 1e-7 * total_arc;

    // Compute SST: interpolated stagnation arc length
    let sst = compute_stagnation_arc_length(airfoil, velocity, stag_idx);

    let n = airfoil.n;
    let end_idx = n - 1; // Lower TE node

    // XFOIL XICALC: XSSI(IBL,IS) = MAX( S(I) - SST , XEPS )
    // Lower surface goes forward from stagnation, so S(I) > SST and XSSI > 0
    for i in stag_idx..=end_idx {
        x.push(airfoil.x[i]);
        y.push(airfoil.y[i]);

        // Arc length from stagnation point (S[i] - SST)
        let arc_len = (airfoil.s[i] - sst).max(xeps);
        s.push(arc_len);

        // Use velocity[i] directly (same as upper surface)
        ue.push(velocity[i].abs());
    }

    (x, y, s, ue)
}

/// Calculate friction drag from BL solution
///
/// Uses XFOIL-compatible integration over chord-projected distance
/// with Ue² weighting factor.
pub fn calculate_friction_drag(bl: &BLSolution, _alpha: f64) -> f64 {
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
/// Uses XFOIL-compatible SETBL + BLSOLV + UPDATE algorithm.
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
    solve_viscous_impl(airfoil, alpha_rad, cond, config)
}

/// Internal implementation using SETBL+BLSOLV+UPDATE loop
fn solve_viscous_impl(
    airfoil: &PaneledAirfoil,
    alpha_rad: f64,
    cond: &FlowConditions,
    config: &ViscalConfig,
) -> ViscousResult {
    // Step 1: Solve inviscid panel method
    // TODO: Wake DIJ (solve_coupled_system) causes divergence - needs debugging
    // For now, use basic inviscid solver without wake
    let inviscid = solve_inviscid(airfoil);

    // Get inviscid velocity (unsigned magnitude)
    let qinv = inviscid.velocity_at_nodes(alpha_rad);
    let qinv_mag: Vec<f64> = qinv.iter().map(|&q| q.abs()).collect();

    // Find stagnation point and compute interpolated arc length
    let stag_idx = find_stagnation_point(airfoil, &qinv);
    let sst = compute_stagnation_arc_length(airfoil, &qinv, stag_idx);

    // Initialize SETBL state
    let mut setbl_state = SetblState::new(airfoil, stag_idx, sst, cond);
    setbl_state.init_from_velocity(airfoil, &qinv_mag);

    // Initialize edge velocities
    let mut ue_mag = qinv_mag.clone();

    // Tracking variables
    let mut iterations = 0;
    let mut rmsbl = f64::MAX;
    let mut converged = false;

    // Main VISCAL iteration loop (SETBL + BLSOLV + UPDATE)
    for iter in 0..config.max_iter {
        iterations = iter + 1;

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

        // Step 1: SETBL - Build Newton system
        // Pass qinv_mag (inviscid velocities) so USAV can be computed from QINV + DIJ*MASS
        // DUE2 = UEDG - USAV is the mismatch that drives Newton convergence
        let mut blsolv_input = build_newton_system(
            &mut setbl_state,
            airfoil,
            &inviscid,
            &qinv_mag, // Pass inviscid velocities for USAV computation
            &config.setbl,
        );

        // Step 2: BLSOLV - Solve the block system
        blsolv(&mut blsolv_input);

        // Step 3: UPDATE - Apply Newton deltas with relaxation and update edge velocities
        // This function properly computes UNEW using MASS + VDEL (Newton delta) before
        // applying under-relaxation, matching XFOIL's UPDATE subroutine exactly.
        let (result, new_ue) = apply_newton_update_with_ue(
            &mut setbl_state,
            &blsolv_input.vdel,
            &inviscid,
            &qinv_mag,
            &config.setbl,
        );

        rmsbl = result.rmsbl;
        ue_mag = new_ue;

        log::debug!(
            "VISCAL iter {}: RMSBL={:.16e} RMXBL={:.16e}",
            iter + 1,
            rmsbl,
            result.dmax
        );

        // XFOIL xoper.f: IF(RMSBL .LT. EPS1) -> converged. There is no divergence
        // check in VISCAL; it runs NITER iterations and reports "Convergence failed".
        if rmsbl < config.tol_rmsbl {
            converged = true;
            break;
        }
    }

    // Convert SETBL state to BLSolution for result
    let bl = convert_setbl_to_bl_solution(&setbl_state, airfoil, stag_idx);

    // Calculate forces
    let velocity = inviscid.velocity_at_nodes(alpha_rad);
    let coeffs = integrate_forces(airfoil, &velocity, alpha_rad, cond.mach);

    // Calculate drag
    let cdf = calculate_friction_drag(&bl, alpha_rad);
    let cd = calculate_drag_squire_young(&bl, airfoil.chord);
    let cdp = cd - cdf;

    // Transition locations
    let xtr_upper = find_transition_x(&setbl_state.stations_upper, &setbl_state.upper.xssi, cond.ncrit);
    let xtr_lower = find_transition_x(&setbl_state.stations_lower, &setbl_state.lower.xssi, cond.ncrit);

    ViscousResult {
        alpha: alpha_rad,
        cl: coeffs.cl,
        cd,
        cdf,
        cdp,
        cm: coeffs.cm,
        xtr_upper,
        xtr_lower,
        iterations,
        converged,
        residual: rmsbl,
        dq_source: vec![0.0; airfoil.n], // Not used in SETBL approach
        bl,
    }
}

/// Find transition x-coordinate from station states
fn find_transition_x(stations: &[BLStationState], xssi: &[f64], ncrit: f64) -> f64 {
    for (i, station) in stations.iter().enumerate() {
        if station.ampl >= ncrit && i > 0 {
            // Interpolate transition location
            let ampl_prev = stations[i - 1].ampl;
            let ampl_curr = station.ampl;
            let frac = (ncrit - ampl_prev) / (ampl_curr - ampl_prev);
            return xssi[i - 1] + frac * (xssi[i] - xssi[i - 1]);
        }
    }
    // No transition found - return large value
    xssi.last().copied().unwrap_or(1.0)
}

/// Convert SETBL state to BLSolution format for result output
fn convert_setbl_to_bl_solution(state: &SetblState, airfoil: &PaneledAirfoil, stag_idx: usize) -> BLSolution {
    // Extract upper surface coordinates from airfoil
    let (x_upper, _, s_upper, _) = extract_upper_surface(
        airfoil,
        &vec![0.0; airfoil.n], // Dummy velocity
        stag_idx,
    );
    let (x_lower, _, s_lower, _) = extract_lower_surface(airfoil, &vec![0.0; airfoil.n], stag_idx);

    // Convert station states to NewtonResult format
    let upper: Vec<NewtonResult> = state
        .stations_upper
        .iter()
        .enumerate()
        .map(|(i, s)| NewtonResult {
            theta: s.theta,
            dstar: s.dstar,
            h: s.h,
            hk: s.hk,
            cf: s.cf,
            cd: s.di * s.hs / 2.0,
            hs: s.hs,
            n_amp: s.ampl,
            ue: s.u,
            x: if i < x_upper.len() { x_upper[i] } else { 0.0 },
            iterations: 0,
            residual: 0.0,
            converged: true,
        })
        .collect();

    let lower: Vec<NewtonResult> = state
        .stations_lower
        .iter()
        .enumerate()
        .map(|(i, s)| NewtonResult {
            theta: s.theta,
            dstar: s.dstar,
            h: s.h,
            hk: s.hk,
            cf: s.cf,
            cd: s.di * s.hs / 2.0,
            hs: s.hs,
            n_amp: s.ampl,
            ue: s.u,
            x: if i < x_lower.len() { x_lower[i] } else { 0.0 },
            iterations: 0,
            residual: 0.0,
            converged: true,
        })
        .collect();

    BLSolution {
        upper,
        lower,
        wake: Vec::new(),
        s_upper,
        s_lower,
        s_wake: Vec::new(),
        x_upper,
        x_lower,
        x_wake: Vec::new(),
        y_wake: Vec::new(),
        i_stag: stag_idx,
        s_tr_upper: state.march_upper.itran as f64 * 0.01, // Approximate
        s_tr_lower: state.march_lower.itran as f64 * 0.01,
    }
}

/// Solve inviscid-only analysis (no BL)
pub fn solve_inviscid_only(airfoil: &PaneledAirfoil, alpha_rad: f64, mach: f64) -> AeroCoefficients {
    let inviscid = solve_inviscid(airfoil);
    // Use node-based velocities for consistency with XFOIL
    let velocity = inviscid.velocity_at_nodes(alpha_rad);
    integrate_forces(airfoil, &velocity, alpha_rad, mach)
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

        // θ should be positive at all stations
        for (i, r) in result.bl.upper.iter().enumerate() {
            assert!(r.theta > 0.0, "θ should be positive at station {}: {}", i, r.theta);
        }

        // Note: Wake handling not yet implemented in SETBL-based solver
        // Wake tests will be added when wake marching is implemented
    }

    #[test]
    #[ignore = "S9: VISCAL loop closure"]
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
    #[ignore = "S9: VISCAL loop closure"]
    fn test_solve_viscous_positive_lift() {
        let geom = naca_4digit("0012", 80).unwrap();
        let airfoil = create_paneled_airfoil(&geom);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = ViscalConfig::default();

        let alpha = 5.0_f64.to_radians();
        let result = solve_viscous(&airfoil, alpha, &cond, &config);

        // At alpha=5°, should have positive lift
        assert!(result.cl > 0.2, "CL={} should be positive at α=5°", result.cl);

        // Drag should be positive (exact value depends on coupling quality)
        assert!(result.cd > 0.0, "CD={} should be positive", result.cd);

        // Friction drag should be a component
        assert!(result.cdf > 0.0, "CDf={} should be positive", result.cdf);
    }

    #[test]
    #[ignore = "S9: VISCAL loop closure"]
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

    /// Test that stagnation index matches XFOIL's IST
    ///
    /// XFOIL reference (from instrumented run with 160-panel NACA 0012):
    /// - IST = 80 (1-based), meaning panel 80 has positive velocity
    /// - stag_idx in YFoil should be the first lower surface panel (IST+1 equivalent)
    #[test]
    fn test_stagnation_index_matches_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        // Generate geometry matching the XFOIL test case
        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let inviscid = solve_inviscid(&airfoil);
        let velocity = inviscid.velocity_at_nodes(0.0);

        let stag_idx = find_stagnation_point(&airfoil, &velocity);

        // XFOIL IST = 80 (1-based), meaning index 79 (0-based) has positive velocity
        // YFoil stag_idx should be 80 (0-based) = first panel with non-positive velocity
        // This is IST+1 in XFOIL's 1-based indexing = 81
        // So stag_idx should be 80 (0-based)
        assert_eq!(
            stag_idx, 80,
            "stag_idx should be 80 (0-based), matching XFOIL IST+1=81 (1-based)"
        );

        // Verify the sign change: velocity[stag_idx-1] > 0 and velocity[stag_idx] <= 0
        assert!(
            velocity[stag_idx - 1] > 0.0,
            "Panel before stag should have positive velocity: {}",
            velocity[stag_idx - 1]
        );
        assert!(
            velocity[stag_idx] <= 0.0,
            "Panel at stag should have non-positive velocity: {}",
            velocity[stag_idx]
        );
    }

    /// Test SST computation matches XFOIL's STFIND
    ///
    /// XFOIL reference:
    /// - SST = 1.019621842270175
    /// - S(IST=80) = 1.019418602014701
    /// - S(IST+1=81) = 1.021246928567578
    /// - GAM(IST) = 0.01676756915... (positive)
    /// - GAM(IST+1) = -0.13407160... (negative)
    #[test]
    fn test_stagnation_arc_length_matches_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let inviscid = solve_inviscid(&airfoil);
        let velocity = inviscid.velocity_at_nodes(0.0);

        let stag_idx = find_stagnation_point(&airfoil, &velocity);
        let sst = compute_stagnation_arc_length(&airfoil, &velocity, stag_idx);

        // XFOIL reference value (first iteration, inviscid)
        let xfoil_sst = 1.019621842270175;

        // Compare - should match within floating point tolerance
        // Note: The exact match depends on having identical geometry
        let rel_err = (sst - xfoil_sst).abs() / xfoil_sst;
        assert!(
            rel_err < 1e-6,
            "SST should match XFOIL: yfoil={:.15e} xfoil={:.15e} rel_err={:.2e}",
            sst,
            xfoil_sst,
            rel_err
        );
    }

    /// Test upper surface arc lengths match XFOIL's XICALC
    ///
    /// XFOIL reference (first few stations):
    /// IBL  IPAN        S(IPAN)              XSSI
    ///  2    80   1.019418602014701e+00   2.032402554739132e-04
    ///  3    79   1.017517698063360e+00   2.104144206815040e-03
    ///  4    78   1.015530615921399e+00   4.091226348776011e-03
    #[test]
    fn test_upper_surface_arc_lengths_match_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let inviscid = solve_inviscid(&airfoil);
        let velocity = inviscid.velocity_at_nodes(0.0);

        let stag_idx = find_stagnation_point(&airfoil, &velocity);
        let (_, _, s_upper, _) = extract_upper_surface(&airfoil, &velocity, stag_idx);

        // XFOIL reference values for first few stations
        let xfoil_xssi_upper = [
            2.032402554739132e-04, // IBL=2, IPAN=80
            2.104144206815040e-03, // IBL=3, IPAN=79
            4.091226348776011e-03, // IBL=4, IPAN=78
            6.157629882446392e-03, // IBL=5, IPAN=77
            8.319496927608583e-03, // IBL=6, IPAN=76
        ];

        // Note: s_upper[0] corresponds to IBL=2 (first actual station after stagnation)
        for (i, &xfoil_val) in xfoil_xssi_upper.iter().enumerate() {
            let rel_err = (s_upper[i] - xfoil_val).abs() / xfoil_val;
            assert!(
                rel_err < 1e-6,
                "Upper XSSI[{}] mismatch: yfoil={:.15e} xfoil={:.15e} rel_err={:.2e}",
                i,
                s_upper[i],
                xfoil_val,
                rel_err
            );
        }
    }

    /// Test lower surface arc lengths match XFOIL's XICALC
    ///
    /// XFOIL reference (first few stations):
    /// IBL  IPAN        S(IPAN)              XSSI
    ///  2    81   1.021246928567578e+00   1.625086297402545e-03
    ///  3    82   1.023074108325968e+00   3.452266055793185e-03
    ///  4    83   1.024980294773205e+00   5.358452503029687e-03
    #[test]
    fn test_lower_surface_arc_lengths_match_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let inviscid = solve_inviscid(&airfoil);
        let velocity = inviscid.velocity_at_nodes(0.0);

        let stag_idx = find_stagnation_point(&airfoil, &velocity);
        let (_, _, s_lower, _) = extract_lower_surface(&airfoil, &velocity, stag_idx);

        // XFOIL reference values for first few stations
        let xfoil_xssi_lower = [
            1.625086297402545e-03, // IBL=2, IPAN=81
            3.452266055793185e-03, // IBL=3, IPAN=82
            5.358452503029687e-03, // IBL=4, IPAN=83
            7.369101817904955e-03, // IBL=5, IPAN=84
            9.478965571405595e-03, // IBL=6, IPAN=85
        ];

        // Note: s_lower[0] corresponds to IBL=2 (first station at stag_idx)
        for (i, &xfoil_val) in xfoil_xssi_lower.iter().enumerate() {
            let rel_err = (s_lower[i] - xfoil_val).abs() / xfoil_val;
            assert!(
                rel_err < 1e-6,
                "Lower XSSI[{}] mismatch: yfoil={:.15e} xfoil={:.15e} rel_err={:.2e}",
                i,
                s_lower[i],
                xfoil_val,
                rel_err
            );
        }
    }

    /// Test XEPS calculation matches XFOIL
    ///
    /// XFOIL reference: XEPS = 2.039242600293527e-07
    #[test]
    fn test_xeps_matches_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        // Compute XEPS the same way as in extract_upper_surface
        let total_arc = airfoil.s[airfoil.n - 1] - airfoil.s[0];
        let xeps = 1e-7 * total_arc;

        let xfoil_xeps = 2.039242600293527e-07;

        let rel_err = (xeps - xfoil_xeps).abs() / xfoil_xeps;
        assert!(
            rel_err < 1e-6,
            "XEPS mismatch: yfoil={:.15e} xfoil={:.15e} rel_err={:.2e}",
            xeps,
            xfoil_xeps,
            rel_err
        );
    }

    /// Test that IBLTE counts match XFOIL
    ///
    /// XFOIL reference:
    /// - IBLTE(1) = 81 (upper surface stations including stagnation)
    /// - IBLTE(2) = 81 (lower surface stations including stagnation)
    #[test]
    fn test_surface_station_counts_match_xfoil() {
        use crate::geometry::{repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 160).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 160, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let inviscid = solve_inviscid(&airfoil);
        let velocity = inviscid.velocity_at_nodes(0.0);

        let stag_idx = find_stagnation_point(&airfoil, &velocity);
        let (_, _, s_upper, _) = extract_upper_surface(&airfoil, &velocity, stag_idx);
        let (_, _, s_lower, _) = extract_lower_surface(&airfoil, &velocity, stag_idx);

        // XFOIL has IBLTE(1) = 81 and IBLTE(2) = 81
        // This includes IBL=1 (stagnation) plus the surface stations
        // Our arrays don't include the stagnation point explicitly,
        // so we should have IBLTE - 1 = 80 stations
        // But actually XFOIL's IPAN starts at IBL=2 with IPAN=IST, so there are 80 actual panels

        // Upper surface: from stag_idx-1 down to 0, that's stag_idx stations
        // XFOIL: IBLTE(1) = IST + 1 = 81 (includes IBL=1 which has XSSI=0)
        // YFoil: s_upper has stag_idx entries
        assert_eq!(
            s_upper.len(),
            80,
            "Upper surface should have 80 stations (IST=80, from panel 79 to 0)"
        );

        // Lower surface: from stag_idx to n-1
        // XFOIL: IBLTE(2) = N - IST = 160 - 80 = 80, plus 1 for IBL=1 = 81
        // But wait, let's check: in IBLPAN, lower goes from IST+1 to N, that's N-IST panels
        // N=160, IST=80, so 160-80=80 panels, plus IBL=1 = 81 total
        // YFoil: from stag_idx to n-1, that's n - stag_idx = 160 - 80 = 80 stations
        assert_eq!(
            s_lower.len(),
            80,
            "Lower surface should have 80 stations (from panel 80 to 159)"
        );
    }
}
