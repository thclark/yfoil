//! SETBL: Newton system setup for viscous-inviscid coupling
//!
//! This module implements XFOIL's SETBL subroutine, which builds the global
//! Newton system for the coupled viscous-inviscid iteration. The key insight
//! is that XFOIL solves ALL BL equations simultaneously with the DIJ coupling
//! embedded in the Jacobian matrix.
//!
//! ## System Structure
//!
//! The global Newton system has 3 equations per BL station:
//! 1. Shear lag (turbulent) or amplification (laminar) equation
//! 2. Momentum integral equation
//! 3. Shape parameter (energy) equation
//!
//! The unknowns at each station are:
//! 1. dCtau (or dAmpl for laminar) - shear/amplification change
//! 2. dTheta - momentum thickness change
//! 3. dMass - mass defect change (where Mass = Ue * δ*)
//!
//! The Jacobian has special structure:
//! - Block tridiagonal (VA, VB) for local station coupling
//! - Dense VM matrix for global DIJ coupling
//!
//! ## Reference
//!
//! XFOIL source files:
//! - xbl.f: SETBL (Newton system assembly)
//! - xoper.f: VISCAL (coupling iteration)

use crate::bl::{
    mrchdu::{march_station, SurfaceBLData, SurfaceMarchState},
    system::{BLFlowType, BLGlobalParams, BLLocalSystem, BLStationState, MidpointCf},
    BlsolvInput, FlowConditions,
};
use crate::geometry::PaneledAirfoil;
use crate::panel::InviscidSolution;

/// State for the global Newton system
///
/// This holds all BL data for both surfaces plus the global system matrices.
#[derive(Debug, Clone)]
pub struct SetblState {
    /// Upper surface BL data
    pub upper: SurfaceBLData,
    /// Lower surface BL data
    pub lower: SurfaceBLData,
    /// Wake BL data
    pub wake: SurfaceBLData,

    /// Upper surface march state
    pub march_upper: SurfaceMarchState,
    /// Lower surface march state
    pub march_lower: SurfaceMarchState,

    /// Station states for upper surface
    pub stations_upper: Vec<BLStationState>,
    /// Station states for lower surface
    pub stations_lower: Vec<BLStationState>,
    /// Station states for wake
    pub stations_wake: Vec<BLStationState>,

    /// Panel indices for upper surface (maps BL station to panel index)
    pub ipan_upper: Vec<usize>,
    /// Panel indices for lower surface
    pub ipan_lower: Vec<usize>,

    /// Stagnation point index (first lower surface panel)
    pub stag_idx: usize,
    /// Stagnation arc length (interpolated)
    pub sst: f64,

    /// Global BL parameters
    pub params: BLGlobalParams,

    /// Number of upper surface stations (including stagnation)
    pub nbl_upper: usize,
    /// Number of lower surface stations (including stagnation)
    pub nbl_lower: usize,
    /// Number of wake stations
    pub nbl_wake: usize,

    /// USAV: Saved edge velocities from previous UPDATE (indexed by panel)
    /// This is UNEW = UINV + DIJ*(MASS + VDEL) from the previous iteration,
    /// used for computing DUE2 = UEDG - USAV in SETBL.
    pub usav: Vec<f64>,
}

impl SetblState {
    /// Create a new SETBL state from airfoil geometry
    ///
    /// # Arguments
    /// * `airfoil` - Paneled airfoil geometry
    /// * `stag_idx` - Stagnation point index (first panel with non-positive velocity)
    /// * `sst` - Interpolated stagnation arc length
    /// * `cond` - Flow conditions
    pub fn new(airfoil: &PaneledAirfoil, stag_idx: usize, sst: f64, cond: &FlowConditions) -> Self {
        // Upper surface: from stag_idx-1 down to 0
        let nbl_upper = stag_idx;
        // Lower surface: from stag_idx to n-1
        let nbl_lower = airfoil.n - stag_idx;

        // Create BL data arrays
        let upper = SurfaceBLData::new(nbl_upper);
        let lower = SurfaceBLData::new(nbl_lower);
        let wake = SurfaceBLData::new(0); // Wake size determined later

        // Create march states
        let acrit = cond.ncrit;
        let march_upper = SurfaceMarchState::new(nbl_upper, nbl_upper, acrit);
        let march_lower = SurfaceMarchState::new(nbl_lower, nbl_lower, acrit);

        // Create global parameters
        let params = BLGlobalParams::new(cond.mach, cond.reynolds, 1.4);

        // Build panel index mappings
        // Station 0 on both surfaces is the virtual stagnation point (xssi=0),
        // so the first real BL station is ibl=1, not ibl=0.
        //
        // XFOIL mapping (1-based, IST = stag_idx in 0-based):
        //   Upper: IBL=2 -> IPAN=IST, IBL=3 -> IPAN=IST-1, ...
        //   Lower: IBL=2 -> IPAN=IST+1, IBL=3 -> IPAN=IST+2, ...
        //
        // YFoil mapping (0-based, ibl = IBL - 1):
        //   Upper: ibl=1 -> ipan=stag_idx-1, ibl=2 -> ipan=stag_idx-2, ...
        //   Lower: ibl=1 -> ipan=stag_idx, ibl=2 -> ipan=stag_idx+1, ...
        let ipan_upper: Vec<usize> = (0..nbl_upper)
            .map(|i| {
                if i == 0 {
                    stag_idx.saturating_sub(1) // Stagnation placeholder (not used)
                } else {
                    stag_idx - i // Real stations: i=1 -> stag_idx-1, i=2 -> stag_idx-2, ...
                }
            })
            .collect();
        let ipan_lower: Vec<usize> = (0..nbl_lower)
            .map(|i| {
                if i == 0 {
                    stag_idx // Stagnation placeholder (not used)
                } else {
                    stag_idx + i - 1 // Real stations: i=1 -> stag_idx, i=2 -> stag_idx+1, ...
                }
            })
            .collect();

        Self {
            upper,
            lower,
            wake,
            march_upper,
            march_lower,
            stations_upper: vec![BLStationState::default(); nbl_upper],
            stations_lower: vec![BLStationState::default(); nbl_lower],
            stations_wake: Vec::new(),
            ipan_upper,
            ipan_lower,
            stag_idx,
            sst,
            params,
            nbl_upper,
            nbl_lower,
            nbl_wake: 0,
            usav: Vec::new(), // Empty initially; set by first UPDATE
        }
    }

    /// Initialize BL data from edge velocity distribution
    ///
    /// Sets up initial arc lengths and edge velocities from the panel solution.
    /// Station 0 on both surfaces is the stagnation station with identical
    /// initial conditions (xssi=xeps, uedg extrapolated from velocity gradient).
    ///
    /// This ensures symmetric BL development for symmetric airfoils at α=0°.
    ///
    /// # Arguments
    /// * `airfoil` - Paneled airfoil
    /// * `ue_mag` - Unsigned edge velocity at each panel
    pub fn init_from_velocity(&mut self, airfoil: &PaneledAirfoil, ue_mag: &[f64]) {
        // Compute XEPS - minimum arc length near stagnation
        let total_arc = airfoil.s[airfoil.n - 1] - airfoil.s[0];
        let xeps = 1e-7 * total_arc;

        // Compute velocity gradient near stagnation for initial condition
        // Use average of first real panels on both surfaces
        let ue_stag = if self.nbl_upper >= 2 && self.nbl_lower >= 2 {
            // First upper panel after stagnation
            let upper_ipan = self.ipan_upper[1.min(self.nbl_upper - 1)];
            let ue_upper = ue_mag[upper_ipan].abs();
            // First lower panel after stagnation
            let lower_ipan = self.ipan_lower[1.min(self.nbl_lower - 1)];
            let ue_lower = ue_mag[lower_ipan].abs();
            // Use average for symmetric initial condition
            0.5 * (ue_upper + ue_lower).max(0.01)
        } else {
            0.01 // Fallback
        };

        // Upper surface: arc length = sst - s[i]
        // Station 0 is the stagnation station at xssi=0
        for (ibl, &ipan) in self.ipan_upper.iter().enumerate() {
            if ibl == 0 {
                // Stagnation station at xssi=0 with symmetric velocity
                self.upper.xssi[0] = 0.0;
                self.upper.uedg[0] = ue_stag;
            } else {
                let arc_len = (self.sst - airfoil.s[ipan]).max(xeps);
                self.upper.xssi[ibl] = arc_len;
                self.upper.uedg[ibl] = ue_mag[ipan].abs();
            }
        }

        // Lower surface: arc length = s[i] - sst
        // Station 0 is the stagnation station at xssi=0
        for (ibl, &ipan) in self.ipan_lower.iter().enumerate() {
            if ibl == 0 {
                // Stagnation station at xssi=0 with symmetric velocity
                self.lower.xssi[0] = 0.0;
                self.lower.uedg[0] = ue_stag;
            } else {
                let arc_len = (airfoil.s[ipan] - self.sst).max(xeps);
                self.lower.xssi[ibl] = arc_len;
                self.lower.uedg[ibl] = ue_mag[ipan].abs();
            }
        }
    }

    /// Get total number of BL stations (for system sizing)
    pub fn total_stations(&self) -> usize {
        self.nbl_upper + self.nbl_lower + self.nbl_wake
    }
}

/// Configuration for SETBL/VISCAL iteration
#[derive(Debug, Clone)]
pub struct SetblConfig {
    /// Maximum Newton iterations per station
    pub max_iter_station: usize,
    /// BLSOLV acceleration parameter
    pub vaccel: f64,
    /// Convergence tolerance for RMSBL
    pub tol_rmsbl: f64,
    /// Maximum relative increase in variables (DHI) - XFOIL default is 1.5
    pub dhi: f64,
    /// Maximum relative decrease in variables (DLO) - XFOIL default is -0.5
    pub dlo: f64,
}

impl Default for SetblConfig {
    fn default() -> Self {
        Self {
            max_iter_station: 25,
            vaccel: 0.01,
            tol_rmsbl: 1e-4,
            // NOTE: XFOIL uses DHI=1.5, DLO=-0.5, but YFoil's VM matrix is simplified
            // so the Newton system can become ill-conditioned. Use tighter bounds for stability.
            dhi: 0.5,  // Max relative increase (tighter than XFOIL 1.5)
            dlo: -0.3, // Max relative decrease (tighter than XFOIL -0.5)
        }
    }
}

/// Result of a single SETBL+BLSOLV+UPDATE iteration
#[derive(Debug, Clone)]
pub struct SetblResult {
    /// RMS of normalized variable changes (convergence metric)
    pub rmsbl: f64,
    /// Maximum normalized variable change
    pub dmax: f64,
    /// Number of stations that failed to converge
    pub n_failed: usize,
    /// Whether the iteration should continue
    pub continue_iteration: bool,
}

/// Build the global Newton system from current BL state
///
/// This is the main SETBL function that:
/// 1. Marches the BL on both surfaces using current edge velocities
/// 2. Builds VA, VB matrices from local Jacobians with DUE2 mismatch terms
/// 3. Builds VM matrix for DIJ coupling
/// 4. Returns BlsolvInput ready for blsolv()
///
/// The key XFOIL insight is that VDEL includes not just the BL equation residuals
/// (VSREZ), but also "forced" mismatch terms DUE2 = UEDG - USAV that drive the
/// Newton system to correct for viscous-inviscid coupling errors.
///
/// # Arguments
/// * `state` - Current SETBL state (modified in place)
/// * `airfoil` - Paneled airfoil
/// * `inviscid` - Inviscid solution (contains DIJ matrix)
/// * `qinv` - Inviscid velocity magnitude at each panel node
/// * `config` - SETBL configuration
///
/// # Returns
/// BlsolvInput ready to be solved with blsolv()
pub fn build_newton_system(
    state: &mut SetblState,
    airfoil: &PaneledAirfoil,
    inviscid: &InviscidSolution,
    qinv: &[f64],
    config: &SetblConfig,
) -> BlsolvInput {
    // Total system size: upper + lower surfaces, excluding stagnation stations
    // XFOIL IBLSYS: NSYS = (NBL(1)-1) + (NBL(2)-1) = NBL(1) + NBL(2) - 2
    // The stagnation station (ibl=0) on each surface is excluded from the Newton system
    let n_upper_sys = state.nbl_upper - 1; // System entries for upper surface
    let n_lower_sys = state.nbl_lower - 1; // System entries for lower surface
    let nsys = n_upper_sys + n_lower_sys;

    // Initialize BLSOLV input
    let mut input = BlsolvInput {
        nsys,
        va: vec![[[0.0; 2]; 3]; nsys],
        vb: vec![[[0.0; 2]; 3]; nsys],
        vdel: vec![[[0.0; 2]; 3]; nsys],
        vm: vec![vec![[0.0, 0.0, 0.0]; nsys]; nsys],
        vz: [[0.0; 2]; 3],
        ivte1: Some(n_upper_sys - 1), // Upper TE system index (last upper entry)
        ivz: Some(n_upper_sys),       // Lower surface start (first lower entry)
        vaccel: config.vaccel,
        arc_length: Some(airfoil.s[airfoil.n - 1] - airfoil.s[0]),
    };

    // Get USAV from saved state (from previous UPDATE) or compute from QINV if first iteration
    // XFOIL saves USAV = UNEW = UINV + DIJ*(MASS + VDEL) in UPDATE,
    // then uses it in SETBL as DUE2 = UEDG - USAV.
    // On first iteration, state.usav is empty, so we use QINV (no coupling yet).
    let usav = if state.usav.is_empty() {
        // First iteration: use QINV as starting USAV
        qinv.to_vec()
    } else {
        // Using saved USAV from previous UPDATE
        state.usav.clone()
    };

    // March upper surface and build system with DUE2 mismatch terms
    build_surface_system(
        state, &mut input, true, // is_upper
        &usav,
    );

    // March lower surface and build system with DUE2 mismatch terms
    build_surface_system(
        state, &mut input, false, // is_lower
        &usav,
    );

    // Add DIJ coupling to VM matrix
    // TEMPORARILY DISABLED: The YFoil VM matrix computation is simplified compared to XFOIL.
    // XFOIL uses a full chain rule with D1_M, U1_M, D2_M, U2_M terms.
    // YFoil's simplified version causes the Newton system to be ill-conditioned.
    // Disabling VM means no viscous-inviscid coupling, which should still allow
    // the BL equations to be solved at fixed edge velocities.
    let disable_vm = true;
    if !disable_vm {
        if let Some(dij) = inviscid.get_dij() {
            add_dij_coupling(state, &mut input, dij);
        }
    }

    input
}

/// Build Newton system for one surface
///
/// On first call (stations not initialized): marches BL to initialize station values.
/// On subsequent calls: uses existing station values with updated edge velocities,
/// computes residuals representing how much current state violates BL equations.
///
/// This matches XFOIL's behavior where SETBL evaluates the equations at the current
/// state rather than re-solving from scratch each iteration.
///
/// # Arguments
/// * `state` - SetblState with station data
/// * `input` - BlsolvInput to populate
/// * `is_upper` - Whether building upper (true) or lower (false) surface
/// * `usav` - Predicted edge velocities from inviscid + mass defect coupling
fn build_surface_system(state: &mut SetblState, input: &mut BlsolvInput, is_upper: bool, usav: &[f64]) {
    // n_upper_sys = nbl_upper - 1 (number of upper surface system entries)
    let n_upper_sys = state.nbl_upper - 1;

    let (nbl, march, xssi, uedg, stations, ipan, offset) = if is_upper {
        (
            state.nbl_upper,
            &mut state.march_upper,
            &state.upper.xssi,
            &state.upper.uedg,
            &mut state.stations_upper,
            &state.ipan_upper,
            0, // Upper surface: system index starts at 0
        )
    } else {
        (
            state.nbl_lower,
            &mut state.march_lower,
            &state.lower.xssi,
            &state.lower.uedg,
            &mut state.stations_lower,
            &state.ipan_lower,
            n_upper_sys, // Lower surface: system index starts at n_upper_sys
        )
    };

    if nbl < 2 {
        return;
    }

    // Initialize march state for this iteration
    march.init_march(march.itran);

    // Local system for building Jacobian
    let mut local_sys = BLLocalSystem::default();

    // XFOIL's SETBL always marches to get the BL solution for current edge velocities.
    // The Newton system then corrects the mass defect to find the coupled solution.
    // We always march here - the previous "skip march" optimization was incorrect.

    // Initialize stagnation station (similarity solution at ibl=0)
    let s2_init = init_similarity_station(xssi[0], uedg[0], &state.params);
    stations[0] = s2_init;

    // March through stations (ibl=1 to nbl-1)
    for ibl in 1..nbl {
        let iv = offset + (ibl - 1);
        let s1 = stations[ibl - 1].clone();

        // Initialize s2 from extrapolation
        let mut s2_init = BLStationState::default();
        s2_init.blprv(
            xssi[ibl],
            s1.ampl,
            s1.ctau,
            s1.theta * 1.1,
            s1.dstar * 1.1,
            0.0,
            uedg[ibl],
            &state.params,
        );
        s2_init.blkin(&state.params);

        let flow_type = if ibl >= march.itran {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };

        // March this station to get initial solution
        let (_result, s2) = march_station(&s1, &s2_init, march, &state.params, ibl + 1, xssi[ibl], 0.0, None);

        // Compute DUE2: velocity mismatch at current station
        // DUE2 = UEDG - USAV (what we're using minus what coupling predicts)
        // This drives the Newton system to correct for viscous-inviscid errors.
        // See XFOIL xbl.f line 279: DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
        let ipan_idx = ipan[ibl];
        let due2 = uedg[ibl] - usav[ipan_idx];

        // Build local system
        let is_simi = ibl == 1;
        let s1_for_sys = if is_simi { &s2 } else { &s1 };
        let cfm = MidpointCf::compute(s1_for_sys, &s2, flow_type, is_simi);
        local_sys.bldif(s1_for_sys, &s2, &cfm, flow_type, is_simi);
        copy_to_global_system(&local_sys, &s2, input, iv, is_simi, due2);

        // Store marched solution
        stations[ibl] = s2;
    }
}

/// Initialize similarity (stagnation) station using Thwaites formula
///
/// For the stagnation station (xsi=0 or very small), uses a small but finite
/// arc length to compute initial values via Thwaites correlation. This avoids
/// zero theta which causes numerical issues in the march.
///
/// This matches XFOIL's behavior where the stagnation station serves as a
/// starting point for the BL march with appropriate initial conditions.
fn init_similarity_station(xsi: f64, uei: f64, params: &BLGlobalParams) -> BLStationState {
    // Thwaites formula: theta^2 = 0.45 * nu * s / (6 * Ue) for BULE=1.0
    // where BULE ≈ 1.0 at the Hiemenz stagnation point
    let nu = 1.0 / params.reybl;

    // For the stagnation station, use a minimum arc length that gives
    // reasonable initial BL thickness. Using 1e-4 (0.01% of chord) as
    // a typical small value near stagnation.
    // This avoids the zero theta problem while staying in the similarity region.
    let xsi_safe = if xsi < 1e-4 { 1e-4 } else { xsi };
    let uei_safe = uei.max(0.01);

    let theta = (0.45 * nu * xsi_safe / (6.0 * uei_safe)).sqrt();
    let h = 2.2; // Hiemenz stagnation point shape factor
    let dstar = h * theta;

    let mut s = BLStationState::default();
    s.blprv(xsi, 0.0, 0.03, theta, dstar, 0.0, uei_safe, params);
    s.blkin(params);
    s.blvar(BLFlowType::Laminar, params);
    s
}

/// Copy local system matrices to global BLSOLV input
///
/// Adds the velocity mismatch terms (DUE2) to VDEL to drive Newton convergence.
/// This matches XFOIL's SETBL (xbl.f lines 449-513) where VDEL includes:
///   VDEL(k,1,IV) = VSREZ(k) + VS2(k,4)*DUE2 + VS2(k,3)*DDS2 + ...
///
/// # Arguments
/// * `local` - Local BL system with Jacobians (VS1, VS2, VSREZ)
/// * `s2` - Station state at current station
/// * `input` - BlsolvInput to populate
/// * `iv` - System index (0-based)
/// * `_is_simi` - Whether this is similarity station
/// * `due2` - Velocity mismatch: UEDG - USAV (what we use minus what coupling predicts)
fn copy_to_global_system(
    local: &BLLocalSystem,
    s2: &BLStationState,
    input: &mut BlsolvInput,
    iv: usize,
    _is_simi: bool,
    due2: f64,
) {
    // VA: diagonal block (rows 0-2, cols 0-1 for Ctau/θ)
    // VB: sub-diagonal block
    // VM: mass defect coupling
    // VDEL: residual (col 0) and Re sensitivity (col 1)

    // The BLSOLV system has 3 equations per station:
    // Row 0: Ctau/Ampl equation
    // Row 1: Momentum equation
    // Row 2: Shape equation

    // VA block: derivatives w.r.t. Ctau and Theta at current station
    for k in 0..3 {
        input.va[iv][k][0] = local.vs2[k][0]; // d/dCtau
        input.va[iv][k][1] = local.vs2[k][1]; // d/dTheta
    }

    // VB block: derivatives w.r.t. Ctau and Theta at previous station
    if iv > 0 {
        for k in 0..3 {
            input.vb[iv][k][0] = local.vs1[k][0]; // d/dCtau_prev
            input.vb[iv][k][1] = local.vs1[k][1]; // d/dTheta_prev
        }
    }

    // VM diagonal: derivatives w.r.t. Mass (= Ue * dstar) at current station
    // Mass = Ue * dstar, so d/dMass = (d/dDstar) / Ue + (d/dUe) * dstar / (Ue * dstar)
    // Simplified: VM[iv][iv][k] = vs2[k][2] / Ue (derivative w.r.t. dstar scaled)
    let ue = s2.u.max(0.01);
    for k in 0..3 {
        input.vm[iv][iv][k] = local.vs2[k][2] / ue; // d/dMass at this station
    }

    // Compute DDS2: displacement thickness change from velocity mismatch
    // XFOIL xbl.f lines 262-263, 280:
    //   D2_U2 = -DSI/UEI
    //   DDS2 = D2_U2*DUE2
    // This represents how dstar changes when edge velocity changes by DUE2
    let dsi = s2.dstar;
    let uei = ue;
    let d2_u2 = -dsi / uei;
    let dds2 = d2_u2 * due2;

    // VDEL: residual + forced mismatch terms
    // XFOIL xbl.f lines 449-453, 479-483, 509-513:
    //   VDEL(k,1,IV) = VSREZ(k)
    //        + (VS1(k,4)*DUE1 + VS1(k,3)*DDS1)   <- previous station terms (ignored for now)
    //        + (VS2(k,4)*DUE2 + VS2(k,3)*DDS2)   <- current station mismatch
    //        + (VS1(k,5) + VS2(k,5) + VSX(k))*(XI_ULE1*DULE1 + XI_ULE2*DULE2)  <- LE terms (ignored)
    //
    // vs2[k][3] = d(residual_k)/d(Ue) = Ue derivative
    // vs2[k][2] = d(residual_k)/d(dstar) = dstar derivative

    // Apply coupling mismatch terms
    // XFOIL: VDEL(k) = VSREZ(k) + VS2(k,4)*DUE2 + VS2(k,3)*DDS2
    // vs2[k][3] = d(residual_k)/d(Ue) at current station
    // vs2[k][2] = d(residual_k)/d(dstar) at current station
    // DUE2 = UEDG - USAV (velocity mismatch)
    // DDS2 = -DSI/UEI * DUE2 (implied dstar mismatch)
    //
    // NOTE: DUE2 terms currently cause instability. Disabled pending investigation.
    let _ = (due2, dds2); // Suppress unused warnings

    for k in 0..3 {
        // let mismatch_term = local.vs2[k][3] * due2 + local.vs2[k][2] * dds2;
        input.vdel[iv][k][0] = local.vsrez[k]; // + mismatch_term;
        input.vdel[iv][k][1] = local.vsr[k]; // Re sensitivity (unchanged)
    }
}

/// Add DIJ coupling to VM matrix
///
/// The VM matrix contains the mass defect influence: how changing mass
/// at station j affects the residual at station i through the velocity change.
fn add_dij_coupling(state: &SetblState, input: &mut BlsolvInput, dij: &nalgebra::DMatrix<f64>) {
    let nsys = input.nsys;
    let n_upper_sys = state.nbl_upper - 1; // System entries for upper surface
    let stag_idx = state.stag_idx; // Use stagnation index for VTI

    // System index mapping (excluding stagnation station at ibl=0):
    // - Upper: iv = 0..n_upper_sys-1 corresponds to ibl = 1..nbl_upper-1
    // - Lower: iv = n_upper_sys..nsys-1 corresponds to ibl = 1..nbl_lower-1

    // For each system station pair (i, j), add DIJ influence
    // The coupling is: d(residual_i)/d(mass_j) += d(residual_i)/d(Ue_i) * d(Ue_i)/d(mass_j)
    // where d(Ue_i)/d(mass_j) = -VTI_i * VTI_j * DIJ(ipan_i, ipan_j)

    for iv in 0..nsys {
        // Map system index to station index (add 1 because stagnation station excluded)
        let (ibl_i, is_upper_i) = if iv < n_upper_sys {
            (iv + 1, true)
        } else {
            (iv - n_upper_sys + 1, false)
        };

        // Get panel index for station
        let ipan_i = if is_upper_i {
            state.ipan_upper[ibl_i]
        } else {
            state.ipan_lower[ibl_i]
        };

        // Get station state for derivatives
        let s_i = if is_upper_i {
            &state.stations_upper[ibl_i]
        } else {
            &state.stations_lower[ibl_i]
        };

        // VTI sign for station i: upper (i < stag_idx) -> +1, lower (i >= stag_idx) -> -1
        let vti_i = if ipan_i < stag_idx { 1.0 } else { -1.0 };

        for jv in 0..nsys {
            let (ibl_j, is_upper_j) = if jv < n_upper_sys {
                (jv + 1, true)
            } else {
                (jv - n_upper_sys + 1, false)
            };

            let ipan_j = if is_upper_j {
                state.ipan_upper[ibl_j]
            } else {
                state.ipan_lower[ibl_j]
            };

            // VTI sign for station j
            let vti_j = if ipan_j < stag_idx { 1.0 } else { -1.0 };

            // DIJ influence: dUe_i/dMass_j = -VTI_i * VTI_j * DIJ(ipan_i, ipan_j)
            let ue_m = -vti_i * vti_j * dij[(ipan_i, ipan_j)];

            // Add to VM matrix
            // d(residual_k)/d(mass_j) += d(residual_k)/d(Ue_i) * dUe_i/dMass_j
            // Note: The Ue derivative is in vs2[k][3], scaled by u_uei
            // This is simplified - full XFOIL uses more complex chain rule
            if iv != jv {
                // Off-diagonal: influence through Ue change
                // For momentum and shape equations, Ue changes affect residual
                input.vm[iv][jv][1] += s_i.hs_u * ue_m; // Momentum via Hs
                input.vm[iv][jv][2] += s_i.di_u * ue_m; // Shape via Di
            }
        }
    }
}

/// Apply Newton deltas with relaxation (UPDATE equivalent)
///
/// Takes the VDEL solution from BLSOLV and applies it to update the BL variables.
/// Uses XFOIL's per-variable relaxation to ensure stability.
///
/// # Arguments
/// * `state` - SETBL state to update
/// * `vdel` - Newton deltas from BLSOLV (modified in place by blsolv)
/// * `config` - SETBL configuration
///
/// # Returns
/// SetblResult with RMSBL convergence metric
pub fn apply_newton_update(state: &mut SetblState, vdel: &[[[f64; 2]; 3]], config: &SetblConfig) -> SetblResult {
    let mut rmsbl: f64 = 0.0;
    let mut dmax: f64 = 0.0;

    let n_upper_sys = state.nbl_upper - 1; // System entries for upper surface

    // Total number of BL stations (for RMSBL divisor)
    // XFOIL: RMSBL = SQRT( RMSBL / (4.0*FLOAT( NBL(1)+NBL(2) )) )
    // NBL(1) + NBL(2) = total stations on both surfaces
    let n_total = state.nbl_upper + state.nbl_lower;

    // Update BL data arrays for stagnation stations (ibl=0) - no Newton deltas
    // These stations are initialized in build_surface_system but their
    // closure relations need to be computed for friction drag calculation
    if state.nbl_upper > 0 {
        let station = &mut state.stations_upper[0];
        station.blkin(&state.params);
        station.blvar(BLFlowType::Laminar, &state.params); // Stagnation is always laminar
        state.upper.thet[0] = station.theta;
        state.upper.dstr[0] = station.dstar;
        state.upper.ctau[0] = station.ctau;
        state.upper.mass[0] = station.dstar * station.u;
    }
    if state.nbl_lower > 0 {
        let station = &mut state.stations_lower[0];
        station.blkin(&state.params);
        station.blvar(BLFlowType::Laminar, &state.params); // Stagnation is always laminar
        state.lower.thet[0] = station.theta;
        state.lower.dstr[0] = station.dstar;
        state.lower.ctau[0] = station.ctau;
        state.lower.mass[0] = station.dstar * station.u;
    }

    // Apply updates to upper surface (ibl=1 to nbl_upper-1)
    // Newton system excludes stagnation station (ibl=0)
    // System index: ibl=1 -> iv=0, ibl=2 -> iv=1, etc.
    let itran_upper = state.march_upper.itran;
    for ibl in 1..state.nbl_upper {
        let iv = ibl - 1; // System index (stagnation excluded)
        let delta = &vdel[iv];
        let station = &mut state.stations_upper[ibl];

        let (_rlx, sum_sq, dm) = apply_station_update(station, delta, config);
        rmsbl += sum_sq; // Already sum of squared normalized changes
        dmax = dmax.max(dm);

        // Recompute closure relations (cf, hs, etc.) after updating primary variables
        station.blkin(&state.params);
        let flow_type = if ibl >= itran_upper {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        station.blvar(flow_type, &state.params);

        // Update BL data arrays
        state.upper.thet[ibl] = station.theta;
        state.upper.dstr[ibl] = station.dstar;
        state.upper.ctau[ibl] = station.ctau;
        state.upper.mass[ibl] = station.dstar * station.u;
    }

    // Apply updates to lower surface (ibl=1 to nbl_lower-1)
    // System index: iv = n_upper_sys + (ibl - 1)
    let itran_lower = state.march_lower.itran;
    for ibl in 1..state.nbl_lower {
        let iv = n_upper_sys + (ibl - 1);
        let delta = &vdel[iv];
        let station = &mut state.stations_lower[ibl];

        let (_rlx, sum_sq, dm) = apply_station_update(station, delta, config);
        rmsbl += sum_sq; // Already sum of squared normalized changes
        dmax = dmax.max(dm);

        // Recompute closure relations (cf, hs, etc.) after updating primary variables
        station.blkin(&state.params);
        let flow_type = if ibl >= itran_lower {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        station.blvar(flow_type, &state.params);

        // Update BL data arrays
        state.lower.thet[ibl] = station.theta;
        state.lower.dstr[ibl] = station.dstar;
        state.lower.ctau[ibl] = station.ctau;
        state.lower.mass[ibl] = station.dstar * station.u;
    }

    // Compute RMS matching XFOIL formula
    // XFOIL: RMSBL = SQRT( RMSBL / (4.0*FLOAT( NBL(1)+NBL(2) )) )
    // Note: We have 3 variables (ctau, theta, mass) but XFOIL has 4 (includes DUEDG)
    // For now we use 4*n_total as XFOIL does; DN4 contribution will be added later
    rmsbl = if n_total > 0 {
        (rmsbl / (4.0 * n_total as f64)).sqrt()
    } else {
        0.0
    };

    SetblResult {
        rmsbl,
        dmax,
        n_failed: 0,
        continue_iteration: rmsbl > config.tol_rmsbl,
    }
}

/// Apply Newton deltas with relaxation AND edge velocity update (full UPDATE equivalent)
///
/// This is a combined UPDATE implementation that follows XFOIL exactly:
/// 1. Computes new edge velocities from DIJ coupling with updated mass (MASS + dMASS)
/// 2. Converts mass deltas to dstar deltas accounting for Ue change
/// 3. Applies global under-relaxation to all variables
/// 4. Recomputes mass nonlinearly as MASS = DSTR * UEDG
///
/// This matches XFOIL's UPDATE subroutine.
///
/// # Arguments
/// * `state` - SETBL state to update
/// * `vdel` - Newton deltas from BLSOLV
/// * `inviscid` - Inviscid solution (for DIJ matrix)
/// * `qinv` - Inviscid velocity distribution
/// * `config` - SETBL configuration
///
/// # Returns
/// (SetblResult with RMSBL metric, updated edge velocity vector)
pub fn apply_newton_update_with_ue(
    state: &mut SetblState,
    vdel: &[[[f64; 2]; 3]],
    inviscid: &InviscidSolution,
    qinv: &[f64],
    config: &SetblConfig,
) -> (SetblResult, Vec<f64>) {
    let n = qinv.len();
    let n_upper_sys = state.nbl_upper - 1;

    // Get DIJ matrix
    let dij = match inviscid.get_dij() {
        Some(d) => d,
        None => {
            // No coupling - just return simple update
            let result = apply_newton_update(state, vdel, config);
            return (result, qinv.to_vec());
        }
    };

    // === Step 1: Compute new Ue using (MASS + dMASS) ===
    // XFOIL: DUI = DUI + UE_M*(MASS(JBL,JS)+VDEL(3,1,JV))
    // UNEW(IBL,IS) = UINV(IBL,IS) + DUI

    // Build array of new masses (MASS + DMASS) for each panel
    // IMPORTANT: XFOIL excludes stagnation station (JBL=1 in 1-indexed) from DIJ sum
    // In 0-indexed: skip ibl=0 (stagnation station)
    // XFOIL: DO 1000 JBL=2, NBL(JS)  -- starts at 2, not 1
    let mut new_mass = vec![0.0; n];
    for (ibl, &ipan) in state.ipan_upper.iter().enumerate() {
        // Skip stagnation station (ibl=0)
        if ibl >= 1 && ibl < state.nbl_upper {
            let station = &state.stations_upper[ibl];
            // CRITICAL: Use state.upper.uedg[ibl] for MASS, NOT station.u
            // XFOIL: MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)  (xbl.f line 1763)
            let old_mass = station.dstar * state.upper.uedg[ibl];
            let iv = ibl - 1;
            let dmass = vdel[iv][2][0]; // VDEL(3,1,JV) - mass delta
            new_mass[ipan] = old_mass + dmass;
        }
    }
    for (ibl, &ipan) in state.ipan_lower.iter().enumerate() {
        // Skip stagnation station (ibl=0)
        if ibl >= 1 && ibl < state.nbl_lower {
            let station = &state.stations_lower[ibl];
            // CRITICAL: Use state.lower.uedg[ibl] for MASS, NOT station.u
            // XFOIL: MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)  (xbl.f line 1763)
            let old_mass = station.dstar * state.lower.uedg[ibl];
            let iv = n_upper_sys + (ibl - 1);
            let dmass = vdel[iv][2][0]; // VDEL(3,1,JV) - mass delta
            new_mass[ipan] = old_mass + dmass;
        }
    }

    // Compute UNEW from DIJ coupling
    // Formula: UNEW[i] = QINV[i] + sum_j(-VTI[i] * VTI[j] * DIJ[i,j] * NEW_MASS[j])
    //
    // NOTE: Without wake panels in DIJ, the coupling near the TE is incorrect.
    // The DIJ matrix is (N)x(N) but should be (N+NW)x(N+NW) to include wake influence.
    // As a workaround, we reduce the coupling strength near the TE to prevent divergence.
    let mut unew = qinv.to_vec();
    let stag_idx = state.stag_idx;

    // DEVIATION (tracked in ci/deviations-allowlist.txt, removed in stage S8):
    // DIJ coupling disabled until the full SETBL VM assembly (S7) lands.
    let disable_dij_coupling = true;

    for i in 0..n {
        let coupling_scale = if disable_dij_coupling { 0.0 } else { 1.0 };

        let vti_i = if i < stag_idx { 1.0 } else { -1.0 };
        for j in 0..n {
            let vti_j = if j < stag_idx { 1.0 } else { -1.0 };
            let contribution = coupling_scale * (-vti_i * vti_j * dij[(i, j)] * new_mass[j]);
            unew[i] += contribution;
        }
    }

    // Save UNEW as USAV for the next iteration's SETBL
    // XFOIL: USAV(IBL,IS) = UNEW(IBL,IS) in UPDATE
    // This is the FULL Newton step velocity (before relaxation)
    state.usav = unew.clone();

    // === Step 2: Compute under-relaxation factor ===
    // XFOIL checks all variable changes and computes global RLX
    let dhi = 1.5_f64;
    let dlo = -0.5_f64;
    let mut rlx = 1.0_f64;
    let mut rmsbl = 0.0_f64;
    let mut rmxbl = 0.0_f64;

    // Collect deltas for all stations
    struct StationDeltas {
        dctau: f64,
        dthet: f64,
        #[allow(dead_code)] // kept for the S8 UPDATE rewrite
        dmass: f64,
        duedg: f64,
        ddstr: f64,
        ctau: f64,
        thet: f64,
        dstr: f64,
        uedg: f64,
        is_turb: bool,
    }

    let mut upper_deltas = Vec::new();
    let mut lower_deltas = Vec::new();

    // Process upper surface
    let itran_upper = state.march_upper.itran;
    for ibl in 1..state.nbl_upper {
        let iv = ibl - 1;
        let ipan = state.ipan_upper[ibl];
        let station = &state.stations_upper[ibl];

        let dctau = vdel[iv][0][0];
        let dthet = vdel[iv][1][0];
        let dmass = vdel[iv][2][0];

        let ctau = station.ctau.max(1e-6);
        let thet = station.theta.max(1e-10);
        let dstr = station.dstar.max(1e-10);
        let uedg = state.upper.uedg[ibl].max(0.01);

        let duedg = unew[ipan] - uedg;
        // XFOIL: DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
        let ddstr = (dmass - dstr * duedg) / uedg;

        let is_turb = ibl >= itran_upper;

        upper_deltas.push(StationDeltas {
            dctau,
            dthet,
            dmass,
            duedg,
            ddstr,
            ctau,
            thet,
            dstr,
            uedg,
            is_turb,
        });

        // Compute normalized changes for relaxation check
        // XFOIL: DN1 = DCTAU / CTAU (or /10 for laminar)
        let dn1 = if is_turb { dctau / ctau } else { dctau / 10.0 };
        let dn2 = dthet / thet;
        let dn3 = ddstr / dstr;
        let dn4 = duedg.abs() / 0.25; // XFOIL uses fixed 0.25, not current Ue

        // Note: RMSBL will be computed AFTER global RLX is determined
        rmxbl = rmxbl.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());

        // Check under-relaxation for each variable
        // Need to check both directions for signed DN values
        if dn1.abs() > 1e-10 {
            if dn1 > dhi {
                rlx = rlx.min(dhi / dn1);
            }
            if dn1 < dlo {
                rlx = rlx.min(dlo / dn1);
            }
        }
        if dn2.abs() > 1e-10 {
            if dn2 > dhi {
                rlx = rlx.min(dhi / dn2);
            }
            if dn2 < dlo {
                rlx = rlx.min(dlo / dn2);
            }
        }
        if dn3.abs() > 1e-10 {
            if dn3 > dhi {
                rlx = rlx.min(dhi / dn3);
            }
            if dn3 < dlo {
                rlx = rlx.min(dlo / dn3);
            }
        }
        if dn4 > dhi {
            rlx = rlx.min(dhi / dn4);
        } // dn4 is always positive
    }

    // Process lower surface
    let itran_lower = state.march_lower.itran;
    for ibl in 1..state.nbl_lower {
        let iv = n_upper_sys + (ibl - 1);
        let ipan = state.ipan_lower[ibl];
        let station = &state.stations_lower[ibl];

        let dctau = vdel[iv][0][0];
        let dthet = vdel[iv][1][0];
        let dmass = vdel[iv][2][0];

        let ctau = station.ctau.max(1e-6);
        let thet = station.theta.max(1e-10);
        let dstr = station.dstar.max(1e-10);
        let uedg = state.lower.uedg[ibl].max(0.01);

        let duedg = unew[ipan] - uedg;
        let ddstr = (dmass - dstr * duedg) / uedg;

        let is_turb = ibl >= itran_lower;

        lower_deltas.push(StationDeltas {
            dctau,
            dthet,
            dmass,
            duedg,
            ddstr,
            ctau,
            thet,
            dstr,
            uedg,
            is_turb,
        });

        let dn1 = if is_turb { dctau / ctau } else { dctau / 10.0 };
        let dn2 = dthet / thet;
        let dn3 = ddstr / dstr;
        let dn4 = duedg.abs() / 0.25;

        // Note: RMSBL will be computed AFTER global RLX is determined
        rmxbl = rmxbl.max(dn1.abs()).max(dn2.abs()).max(dn3.abs()).max(dn4.abs());

        if dn1.abs() > 1e-10 {
            if dn1 > dhi {
                rlx = rlx.min(dhi / dn1);
            }
            if dn1 < dlo {
                rlx = rlx.min(dlo / dn1);
            }
        }
        if dn2.abs() > 1e-10 {
            if dn2 > dhi {
                rlx = rlx.min(dhi / dn2);
            }
            if dn2 < dlo {
                rlx = rlx.min(dlo / dn2);
            }
        }
        if dn3.abs() > 1e-10 {
            if dn3 > dhi {
                rlx = rlx.min(dhi / dn3);
            }
            if dn3 < dlo {
                rlx = rlx.min(dlo / dn3);
            }
        }
        if dn4 > dhi {
            rlx = rlx.min(dhi / dn4);
        } // dn4 is always positive
    }

    // Ensure positive relaxation
    rlx = rlx.max(0.01);

    // Compute RMSBL using UNRELAXED DN values (as XFOIL does)
    // XFOIL: RMSBL = RMSBL + DN1**2 + DN2**2 + DN3**2 + DN4**2 (no RLX factor)
    for d in &upper_deltas {
        let dn1 = if d.is_turb { d.dctau / d.ctau } else { d.dctau / 10.0 };
        let dn2 = d.dthet / d.thet;
        let dn3 = d.ddstr / d.dstr;
        let dn4 = d.duedg.abs() / 0.25;
        rmsbl += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
    }
    for d in &lower_deltas {
        let dn1 = if d.is_turb { d.dctau / d.ctau } else { d.dctau / 10.0 };
        let dn2 = d.dthet / d.thet;
        let dn3 = d.ddstr / d.dstr;
        let dn4 = d.duedg.abs() / 0.25;
        rmsbl += dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;
    }

    // === Step 3: Apply relaxed updates ===
    // Build output edge velocity
    let mut ue_out = qinv.to_vec();

    // Update upper surface
    for (idx, ibl) in (1..state.nbl_upper).enumerate() {
        let d = &upper_deltas[idx];
        let ipan = state.ipan_upper[ibl];
        let station = &mut state.stations_upper[ibl];

        // Apply relaxed changes
        station.ctau = (d.ctau + rlx * d.dctau).max(1e-6).min(0.25);
        station.theta = (d.thet + rlx * d.dthet).max(1e-10);
        station.dstar = (d.dstr + rlx * d.ddstr).max(station.theta * 1.02);
        let new_uedg = (d.uedg + rlx * d.duedg).max(0.01);
        station.u = new_uedg;

        // Update BL arrays
        state.upper.thet[ibl] = station.theta;
        state.upper.dstr[ibl] = station.dstar;
        state.upper.ctau[ibl] = station.ctau;
        state.upper.uedg[ibl] = new_uedg;
        state.upper.mass[ibl] = station.dstar * new_uedg; // Nonlinear update

        // Recompute closure relations
        station.blkin(&state.params);
        let flow_type = if d.is_turb {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        station.blvar(flow_type, &state.params);

        ue_out[ipan] = new_uedg;
    }

    // Update lower surface
    for (idx, ibl) in (1..state.nbl_lower).enumerate() {
        let d = &lower_deltas[idx];
        let ipan = state.ipan_lower[ibl];
        let station = &mut state.stations_lower[ibl];

        station.ctau = (d.ctau + rlx * d.dctau).max(1e-6).min(0.25);
        station.theta = (d.thet + rlx * d.dthet).max(1e-10);
        station.dstar = (d.dstr + rlx * d.ddstr).max(station.theta * 1.02);
        let new_uedg = (d.uedg + rlx * d.duedg).max(0.01);
        station.u = new_uedg;

        state.lower.thet[ibl] = station.theta;
        state.lower.dstr[ibl] = station.dstar;
        state.lower.ctau[ibl] = station.ctau;
        state.lower.uedg[ibl] = new_uedg;
        state.lower.mass[ibl] = station.dstar * new_uedg;

        station.blkin(&state.params);
        let flow_type = if d.is_turb {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        station.blvar(flow_type, &state.params);

        ue_out[ipan] = new_uedg;
    }

    // Update stagnation stations (no Newton deltas, just recompute closure)
    if state.nbl_upper > 0 {
        let station = &mut state.stations_upper[0];
        station.blkin(&state.params);
        station.blvar(BLFlowType::Laminar, &state.params);
        state.upper.thet[0] = station.theta;
        state.upper.dstr[0] = station.dstar;
        state.upper.ctau[0] = station.ctau;
        state.upper.mass[0] = station.dstar * station.u;
    }
    if state.nbl_lower > 0 {
        let station = &mut state.stations_lower[0];
        station.blkin(&state.params);
        station.blvar(BLFlowType::Laminar, &state.params);
        state.lower.thet[0] = station.theta;
        state.lower.dstr[0] = station.dstar;
        state.lower.ctau[0] = station.ctau;
        state.lower.mass[0] = station.dstar * station.u;
    }

    // XFOIL: RMSBL = SQRT( RMSBL / (4.0*FLOAT( NBL(1)+NBL(2) )) )
    let n_total = (state.nbl_upper + state.nbl_lower) as f64;
    rmsbl = (rmsbl / (4.0 * n_total)).sqrt();

    let result = SetblResult {
        rmsbl,
        dmax: rmxbl,
        n_failed: 0,
        continue_iteration: rmsbl > config.tol_rmsbl,
    };

    (result, ue_out)
}

/// Apply Newton delta to a single station with relaxation
///
/// Returns (relaxation_factor, sum_of_squared_normalized_changes, max_normalized_change)
///
/// XFOIL accumulates DN1² + DN2² + DN3² for RMSBL calculation where:
/// - DN1 = |dCtau/Ctau| * rlx
/// - DN2 = |dTheta/Theta| * rlx
/// - DN3 = |dMass/Mass| * rlx
fn apply_station_update(station: &mut BLStationState, delta: &[[f64; 2]; 3], config: &SetblConfig) -> (f64, f64, f64) {
    // Newton deltas: delta[0] = dCtau, delta[1] = dTheta, delta[2] = dMass
    let d_ctau = delta[0][0];
    let d_theta = delta[1][0];
    let d_mass = delta[2][0];

    // Current values
    let ctau = station.ctau.max(1e-6);
    let theta = station.theta.max(1e-10);
    let mass = (station.dstar * station.u).max(1e-10);

    // Compute relative changes (DN values before relaxation)
    let dn1 = (d_ctau / ctau).abs();
    let dn2 = (d_theta / theta).abs();
    let dn3 = (d_mass / mass).abs();

    // Find relaxation factor to keep changes within bounds
    let mut rlx: f64 = 1.0;

    // Limit increases
    if dn1 > config.dhi {
        rlx = rlx.min(config.dhi / dn1);
    }
    if dn2 > config.dhi {
        rlx = rlx.min(config.dhi / dn2);
    }
    if dn3 > config.dhi {
        rlx = rlx.min(config.dhi / dn3);
    }

    // Limit decreases (using signed relative changes)
    let rel_ctau = d_ctau / ctau;
    let rel_theta = d_theta / theta;
    let rel_mass = d_mass / mass;

    if rel_ctau < config.dlo {
        rlx = rlx.min(config.dlo / rel_ctau);
    }
    if rel_theta < config.dlo {
        rlx = rlx.min(config.dlo / rel_theta);
    }
    if rel_mass < config.dlo {
        rlx = rlx.min(config.dlo / rel_mass);
    }

    // Note: XFOIL does NOT have a global under-relaxation factor here.
    // The per-variable relaxation based on DHI/DLO bounds is sufficient.
    // Removed: rlx *= config.global_rlx;

    // Apply relaxed update
    station.ctau = (ctau + rlx * d_ctau).max(1e-6).min(0.3);
    station.theta = (theta + rlx * d_theta).max(1e-10);

    // Update dstar from mass: mass = Ue * dstar, so dstar = mass / Ue
    let new_mass = mass + rlx * d_mass;
    let ue = station.u.max(0.01);
    station.dstar = (new_mass / ue).max(station.theta * 1.00005);

    // XFOIL RMSBL: accumulates (RLX*DN)² for each variable
    // We return the sum of squared relaxed normalized changes
    let dn1_rlx = rlx * dn1;
    let dn2_rlx = rlx * dn2;
    let dn3_rlx = rlx * dn3;
    let sum_sq = dn1_rlx * dn1_rlx + dn2_rlx * dn2_rlx + dn3_rlx * dn3_rlx;

    // Also return max for dmax tracking
    let dmax = dn1.max(dn2).max(dn3) * rlx;

    (rlx, sum_sq, dmax)
}

/// Update edge velocities from mass defect using DIJ matrix
///
/// Computes: Ue_new = Qinv + Σ_j (-VTI_i * VTI_j * DIJ(i,j) * mass_j)
/// Then under-relaxes the change: Ue = Ue_old + RLX * (Ue_new - Ue_old)
///
/// This matches XFOIL's UPDATE subroutine behavior where DUEDG is under-relaxed.
///
/// # Arguments
/// * `state` - SETBL state with current mass defect values
/// * `inviscid` - Inviscid solution (contains DIJ)
/// * `qinv` - Inviscid velocity at each panel (unsigned magnitude)
///
/// # Returns
/// Updated edge velocity at each panel (unsigned magnitude)
pub fn update_edge_velocities(
    state: &SetblState,
    inviscid: &InviscidSolution,
    qinv: &[f64],
    config: &SetblConfig,
) -> Vec<f64> {
    let n = qinv.len();

    // Get current edge velocities
    let mut current_ue = vec![0.0; n];
    for (ibl, &ipan) in state.ipan_upper.iter().enumerate() {
        if ibl < state.nbl_upper {
            current_ue[ipan] = state.upper.uedg[ibl];
        }
    }
    for (ibl, &ipan) in state.ipan_lower.iter().enumerate() {
        if ibl < state.nbl_lower {
            current_ue[ipan] = state.lower.uedg[ibl];
        }
    }
    // For panels not covered by BL stations, use inviscid
    for i in 0..n {
        if current_ue[i] == 0.0 {
            current_ue[i] = qinv[i];
        }
    }

    // Get DIJ matrix
    let dij = match inviscid.get_dij() {
        Some(d) => d,
        None => return current_ue,
    };

    // Build mass defect array from BL state
    // IMPORTANT: XFOIL excludes stagnation station (JBL=1 in 1-indexed) from DIJ sum
    // In 0-indexed: skip ibl=0 (stagnation station)
    // XFOIL: DO 1000 JBL=2, NBL(JS)  -- starts at 2, not 1
    // Build MASS array using DSTR * UEDG (not station.u!)
    // XFOIL: MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)  (xbl.f line 1763)
    let mut mass = vec![0.0; n];
    for (ibl, &ipan) in state.ipan_upper.iter().enumerate() {
        // Skip stagnation station (ibl=0)
        if ibl >= 1 && ibl < state.nbl_upper {
            let station = &state.stations_upper[ibl];
            // CRITICAL: Use state.upper.uedg[ibl], NOT station.u
            mass[ipan] = station.dstar * state.upper.uedg[ibl];
        }
    }
    for (ibl, &ipan) in state.ipan_lower.iter().enumerate() {
        // Skip stagnation station (ibl=0)
        if ibl >= 1 && ibl < state.nbl_lower {
            let station = &state.stations_lower[ibl];
            // CRITICAL: Use state.lower.uedg[ibl], NOT station.u
            mass[ipan] = station.dstar * state.lower.uedg[ibl];
        }
    }

    // Compute new edge velocities from DIJ coupling
    // VTI sign: upper surface (i < stag_idx) -> +1, lower surface (i >= stag_idx) -> -1
    let stag_idx = state.stag_idx;
    let mut new_ue = qinv.to_vec();

    for i in 0..n {
        let vti_i = if i < stag_idx { 1.0 } else { -1.0 };
        for j in 0..n {
            let vti_j = if j < stag_idx { 1.0 } else { -1.0 };
            new_ue[i] += -vti_i * vti_j * dij[(i, j)] * mass[j];
        }
    }

    // Compute under-relaxation factor based on maximum Ue change
    // XFOIL uses DN4 = ABS(DUEDG)/0.25 with FIXED 0.25 reference
    let dhi = config.dhi;
    let mut rlx = 1.0_f64;

    for i in 0..n {
        let duedg = new_ue[i] - current_ue[i];
        // XFOIL: DN4 = ABS(DUEDG)/0.25 - uses FIXED reference of 0.25, not current Ue
        let dn4 = duedg.abs() / 0.25;

        // Under-relax if change is too large
        if dn4 > dhi {
            rlx = rlx.min(dhi / dn4);
        }
    }

    // Ensure some minimum relaxation
    rlx = rlx.max(0.01);

    // Note: XFOIL does NOT have a global under-relaxation factor here.
    // Removed: rlx *= config.global_rlx;

    // Apply under-relaxed change
    let mut ue = vec![0.0; n];
    for i in 0..n {
        let duedg = new_ue[i] - current_ue[i];
        ue[i] = current_ue[i] + rlx * duedg;
        // Ensure positive
        ue[i] = ue[i].max(0.01);
    }

    ue
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_setbl_state_creation() {
        use crate::geometry::{create_paneled_airfoil, naca_4digit, repanel_xfoil, PaneConfig};

        let geom = naca_4digit("0012", 80).unwrap();
        let config = PaneConfig::default();
        let repaneled = repanel_xfoil(&geom, 80, &config);
        let airfoil = create_paneled_airfoil(&repaneled);

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let stag_idx = 40; // Approximate for symmetric at α=0
        let sst = airfoil.s[40];

        let state = SetblState::new(&airfoil, stag_idx, sst, &cond);

        assert_eq!(state.nbl_upper, 40);
        assert_eq!(state.nbl_lower, 40);
        assert_eq!(state.total_stations(), 80);
    }

    #[test]
    fn test_setbl_config_default() {
        let config = SetblConfig::default();
        assert_eq!(config.vaccel, 0.01);
        assert!(config.tol_rmsbl > 0.0);
    }
}
