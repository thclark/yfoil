//! MRCHDU - Mixed-mode BL marching with Ue-Hk characteristic line
//!
//! This module implements XFOIL's MRCHDU subroutine, which marches the
//! boundary layer equations downstream while avoiding the Goldstein
//! singularity by using a quasi-normal approach in Ue-Hk space.
//!
//! ## Algorithm Overview
//!
//! MRCHDU solves the BL equations station-by-station, starting from
//! the similarity solution at the stagnation point. At each station:
//!
//! 1. Set primary variables from current state
//! 2. Newton iteration loop:
//!    a. Call BLPRV, BLKIN to compute secondary variables
//!    b. Check for transition (TRCHEK)
//!    c. Assemble local Newton system (BLSYS or TESYS)
//!    d. Calculate Ue-Hk characteristic slope
//!    e. Replace 4th equation with Ue-Hk constraint
//!    f. Solve 4x4 system with GAUSS
//!    g. Underrelax if needed
//!    h. Update variables
//! 3. Store converged values
//! 4. Copy "2" → "1" for next station
//!
//! ## Reference
//!
//! XFOIL source: xbl.f, SUBROUTINE MRCHDU (lines ~1003-1321)

use super::closure::hkin;
use super::gauss::gauss_solve_4x4;
use super::system::{
    BLFlowType, BLGlobalParams, BLLocalSystem, BLStationState, MidpointCf, TransitionLocation,
    TransitionResult,
};

/// Convergence tolerance for Newton iteration
const DEPS: f64 = 5.0e-6;

/// Weight for Ue-Hk sensitivity (controls how far Hk can deviate)
const SENSWT: f64 = 1000.0;

/// Maximum Newton iterations per station
const MAX_ITER: usize = 25;

/// Result of a single station march
#[derive(Debug, Clone)]
pub enum MarchResult {
    /// Station converged successfully
    Converged,
    /// Station converged with slightly high residual (< 0.1)
    PartiallyConverged { dmax: f64 },
    /// Station failed to converge
    Failed { dmax: f64 },
}

/// State tracking for a single surface during marching
#[derive(Debug, Clone)]
pub struct SurfaceMarchState {
    /// Number of BL stations on this surface
    pub nbl: usize,
    /// Trailing edge station index (1-based, as in XFOIL)
    pub iblte: usize,
    /// Critical amplification factor for this surface
    pub acrit: f64,
    /// Current transition station index (1-based)
    pub itran: usize,
    /// Old transition station index (from previous VISCAL)
    pub itrold: usize,
    /// Whether currently in transition interval
    pub tran: bool,
    /// Whether currently turbulent
    pub turb: bool,
    /// Forced transition arc-length location
    pub xiforc: f64,
    /// Transition arc-length location after march
    pub xssitr: f64,
    /// Whether transition was forced
    pub tforce: bool,
}

impl SurfaceMarchState {
    /// Create new surface march state
    pub fn new(nbl: usize, iblte: usize, acrit: f64) -> Self {
        Self {
            nbl,
            iblte,
            acrit,
            itran: iblte, // Initially set to TE
            itrold: iblte,
            tran: false,
            turb: false,
            xiforc: f64::MAX, // Disabled by default
            xssitr: 0.0,
            tforce: false,
        }
    }

    /// Initialize for a new march (called at start of MRCHDU surface loop)
    ///
    /// # Arguments
    /// * `itrold` - Transition station from previous VISCAL iteration.
    ///              For first iteration, set to iblte (no previous transition).
    ///              For subsequent iterations, set to the detected transition station.
    pub fn init_march(&mut self, itrold: usize) {
        self.itrold = itrold;
        self.tran = false;
        self.turb = false;
        // ITRAN is initialized to ITROLD (from previous iteration)
        // It may be updated during the march if transition moves
        self.itran = itrold;
    }
}

/// BL data arrays for one surface
#[derive(Debug, Clone)]
pub struct SurfaceBLData {
    /// Arc length XSSI(IBL)
    pub xssi: Vec<f64>,
    /// Edge velocity UEDG(IBL)
    pub uedg: Vec<f64>,
    /// Momentum thickness THET(IBL)
    pub thet: Vec<f64>,
    /// Displacement thickness DSTR(IBL)
    pub dstr: Vec<f64>,
    /// Shear stress coeff or amplification CTAU(IBL)
    pub ctau: Vec<f64>,
    /// Mass defect MASS(IBL)
    pub mass: Vec<f64>,
    /// Wall shear stress TAU(IBL)
    pub tau: Vec<f64>,
    /// Dissipation integral DIS(IBL)
    pub dis: Vec<f64>,
    /// Equilibrium Ctau CTQ(IBL)
    pub ctq: Vec<f64>,
    /// Energy thickness DELT(IBL)
    pub delt: Vec<f64>,
    /// Entrainment thickness TSTR(IBL)
    pub tstr: Vec<f64>,
    /// Wake gap (for wake stations) WGAP(IW)
    pub wgap: Vec<f64>,
}

impl SurfaceBLData {
    /// Create BL data arrays with given size
    pub fn new(nbl: usize) -> Self {
        Self {
            xssi: vec![0.0; nbl],
            uedg: vec![0.0; nbl],
            thet: vec![0.0; nbl],
            dstr: vec![0.0; nbl],
            ctau: vec![0.0; nbl],
            mass: vec![0.0; nbl],
            tau: vec![0.0; nbl],
            dis: vec![0.0; nbl],
            ctq: vec![0.0; nbl],
            delt: vec![0.0; nbl],
            tstr: vec![0.0; nbl],
            wgap: Vec::new(),
        }
    }

    /// Get station at index (0-based)
    pub fn get_station(&self, ibl: usize) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.xssi[ibl],
            self.uedg[ibl],
            self.thet[ibl],
            self.dstr[ibl],
            self.ctau[ibl],
            self.mass[ibl],
        )
    }

    /// Set station at index (0-based)
    pub fn set_station(&mut self, ibl: usize, thi: f64, dsi: f64, uei: f64, cti_or_ami: f64) {
        self.thet[ibl] = thi;
        self.dstr[ibl] = dsi;
        self.uedg[ibl] = uei;
        self.ctau[ibl] = cti_or_ami;
        self.mass[ibl] = dsi * uei;
    }
}

/// Set up trailing edge system (TESYS equivalent)
///
/// This creates a "dummy" BL system at the trailing edge that simply
/// enforces continuity of Ctau, theta, and delta* from the TE to the
/// first wake point.
///
/// # Arguments
/// * `sys` - Local system to populate (VS1, VS2, VSREZ)
/// * `s2` - Wake station state (already set up with BLPRV, BLKIN, BLVAR)
/// * `cte` - Trailing edge Ctau (combined from upper and lower)
/// * `tte` - Trailing edge theta (sum of upper and lower)
/// * `dte` - Trailing edge delta* (sum of upper and lower + ANTE)
pub fn tesys(sys: &mut BLLocalSystem, s2: &BLStationState, cte: f64, tte: f64, dte: f64) {
    // Initialize to zero
    for k in 0..4 {
        sys.vsrez[k] = 0.0;
        sys.vsm[k] = 0.0;
        sys.vsr[k] = 0.0;
        sys.vsx[k] = 0.0;
        for l in 0..5 {
            sys.vs1[k][l] = 0.0;
            sys.vs2[k][l] = 0.0;
        }
    }

    // Ctau continuity: CTE = S2
    sys.vs1[0][0] = -1.0;
    sys.vs2[0][0] = 1.0;
    sys.vsrez[0] = cte - s2.ctau;

    // Theta continuity: TTE = T2
    sys.vs1[1][1] = -1.0;
    sys.vs2[1][1] = 1.0;
    sys.vsrez[1] = tte - s2.theta;

    // Delta* continuity: DTE = D2 + DW2
    sys.vs1[2][2] = -1.0;
    sys.vs2[2][2] = 1.0;
    sys.vsrez[2] = dte - s2.dstar - s2.dw;
}

/// Calculate Ue-Hk characteristic slope for the 4th equation
///
/// This computes the sensitivity dUe/dHk along the characteristic line,
/// which is used to prescribe a Ue-Hk combination that avoids the
/// Goldstein singularity.
///
/// # Arguments
/// * `sys` - Local system (VS2, VSREZ)
/// * `s2` - Current station state
/// * `hkref` - Reference Hk value
/// * `ueref` - Reference Ue value
/// * `sens_old` - Previous sensitivity (for averaging)
/// * `itbl` - Newton iteration number (1-based)
///
/// # Returns
/// (new_sens, vs2_row4, vsrez4) - Updated sensitivity and row 4 of system
pub fn calc_ue_hk_characteristic(
    sys: &BLLocalSystem,
    s2: &BLStationState,
    hkref: f64,
    ueref: f64,
    sens_old: f64,
    itbl: usize,
) -> (f64, [f64; 5], f64) {
    // Make copies for GAUSS (it destroys the matrix)
    let mut vtmp = [[0.0f64; 4]; 4];
    let mut vztmp = [0.0f64; 4];

    // Copy the 4x4 portion of VS2
    for k in 0..4 {
        vztmp[k] = sys.vsrez[k];
        for l in 0..4 {
            vtmp[k][l] = sys.vs2[k][l];
        }
    }

    // Set unit dHk in 4th equation
    // HK2 = HK2(T2, D2, U2) via HK2_T2, HK2_D2, HK2_U2
    vtmp[3][0] = 0.0;
    vtmp[3][1] = s2.hk_t;
    vtmp[3][2] = s2.hk_d;
    vtmp[3][3] = s2.hk_u * s2.u_uei;
    vztmp[3] = 1.0;

    // Solve to get dUe response to unit dHk
    gauss_solve_4x4(&mut vtmp, &mut vztmp);

    // VZTMP[3] is now dUe for unit dHk
    // Set SENSWT * (normalized dUe/dHk)
    let sennew = SENSWT * vztmp[3] * hkref / ueref;

    // Average sensitivity for stability
    let sens = if itbl <= 5 {
        sennew
    } else if itbl <= 15 {
        0.5 * (sens_old + sennew)
    } else {
        sens_old
    };

    // Set prescribed Ue-Hk combination (4th row of system)
    let mut vs2_row4 = [0.0f64; 5];
    vs2_row4[0] = 0.0;
    vs2_row4[1] = s2.hk_t * hkref;
    vs2_row4[2] = s2.hk_d * hkref;
    vs2_row4[3] = (s2.hk_u * hkref + sens / ueref) * s2.u_uei;
    vs2_row4[4] = 0.0;

    // Residual: drive Hk and Ue toward reference values
    let vsrez4 = -(hkref * hkref) * (s2.hk / hkref - 1.0) - sens * (s2.u / ueref - 1.0);

    (sens, vs2_row4, vsrez4)
}

/// Set up prescribed Ue equation for similarity station or first wake point
///
/// For these special stations, we simply prescribe Ue to its reference value.
///
/// # Arguments
/// * `s2` - Current station state
/// * `ueref` - Reference Ue value
///
/// # Returns
/// (vs2_row4, vsrez4) - Row 4 of system
pub fn prescribed_ue_equation(s2: &BLStationState, ueref: f64) -> ([f64; 5], f64) {
    let mut vs2_row4 = [0.0f64; 5];
    vs2_row4[0] = 0.0;
    vs2_row4[1] = 0.0;
    vs2_row4[2] = 0.0;
    vs2_row4[3] = s2.u_uei;
    vs2_row4[4] = 0.0;

    let vsrez4 = ueref - s2.u;

    (vs2_row4, vsrez4)
}

/// Limit delta* to maintain minimum Hk (DSLIM equivalent)
///
/// Ensures the kinematic shape parameter Hk doesn't go below the
/// specified limit, which would indicate non-physical separation.
///
/// # Arguments
/// * `dsw` - Delta* without wake gap
/// * `thi` - Momentum thickness
/// * `uei` - Edge velocity
/// * `msq` - Edge Mach number squared
/// * `hklim` - Minimum Hk limit (1.02 on surface, 1.00005 in wake)
///
/// # Returns
/// Limited delta* value
pub fn dslim(dsw: f64, thi: f64, _uei: f64, msq: f64, hklim: f64) -> f64 {
    // Get current H and Hk
    let h = dsw / thi;
    let (hk, _hk_h, _hk_msq) = hkin(h, msq);

    if hk < hklim {
        // Calculate H corresponding to Hklim
        // H = (Hk + 0.028*(1+0.5*M²))*Hk / (1 - 0.014*M²/(1-0.4*M²))
        // Solve for H given Hk=Hklim

        // Simplified: just enforce minimum delta*
        let h_min = hklim + 0.028 * (1.0 + 0.5 * msq) * hklim;
        let dsw_min = h_min * thi;
        dsw_min.max(dsw)
    } else {
        dsw
    }
}

/// Perform Newton iteration at a single BL station
///
/// This is the core of MRCHDU - it iterates to convergence at one station
/// using the local 4x4 Newton system.
///
/// # Arguments
/// * `s1` - Previous station state (already converged)
/// * `s2` - Current station state (to be updated)
/// * `march` - Surface march state
/// * `params` - Global BL parameters
/// * `ibl` - Station index (1-based, as in XFOIL)
/// * `xsi` - Arc length at current station
/// * `dswaki` - Wake gap at current station
/// * `te_values` - Optional (cte, tte, dte) for first wake point
///
/// # Returns
/// (MarchResult, updated s2)
pub fn march_station(
    s1: &BLStationState,
    s2_init: &BLStationState,
    march: &mut SurfaceMarchState,
    params: &BLGlobalParams,
    ibl: usize,
    xsi: f64,
    dswaki: f64,
    te_values: Option<(f64, f64, f64)>,
) -> (MarchResult, BLStationState) {
    let is_simi = ibl == 2;
    let is_wake = ibl > march.iblte;
    let is_first_wake = ibl == march.iblte + 1;

    // At start of each station: reset TRAN and set TURB based on ITRAN
    // (Matches XFOIL MRCHDU: TRAN = .FALSE., TURB = IBL .GE. ITRAN(IS))
    // ITRAN is updated during the march when transition is detected.
    march.tran = false;
    march.turb = ibl >= march.itran;

    // Initialize working variables from input state
    let xsi_i = xsi;
    let mut uei = s2_init.u / s2_init.u_uei; // Convert back to incompressible
    let mut thi = s2_init.theta;
    let mut dsi = s2_init.dstar + s2_init.dw;

    // For laminar stations, amplification comes from the PREVIOUS station (s1.ampl),
    // not from the input fixture. For turbulent stations, use Ctau from fixture.
    let mut ami = if ibl < march.itrold {
        s1.ampl // Accumulated amplification from previous station
    } else {
        0.0 // Not used for turbulent (CTAU is used instead)
    };

    // Handle Ctau initialization based on transition state
    let mut cti = if ibl < march.itrold {
        0.03 // CTI not used for laminar, but initialize to reasonable value
    } else {
        let c = s2_init.ctau;
        if c <= 0.0 { 0.03 } else { c }
    };

    // Enforce minimum delta* for H > 1.02
    let hklim = if !is_wake { 1.02000 } else { 1.00005 };
    dsi = (dsi - dswaki).max(hklim * thi) + dswaki;

    // Create working station state
    let mut s2 = BLStationState::default();
    let mut local_sys = BLLocalSystem::default();

    // Reference values for Ue-Hk characteristic (set on first iteration)
    let mut ueref = 0.0;
    let mut hkref = 0.0;
    let mut sens = 0.0;

    // Transition location (stored when transition detected)
    let mut trans_loc: Option<TransitionLocation> = None;

    // Newton iteration loop
    for itbl in 1..=MAX_ITER {
        // Set up station state
        s2.blprv(xsi_i, ami, cti, thi, dsi, dswaki, uei, params);
        s2.blkin(params);

        // Check for transition (if not similarity and not already turbulent)
        if !is_simi && !march.turb {
            let result = super::system::trchek(
                s1,
                &s2,
                s1.ampl,
                march.acrit,
                march.xiforc,
                params,
            );

            match result {
                TransitionResult::NoTransition { ampl2 } => {
                    ami = ampl2;
                    // ITRAN is NOT updated on NoTransition - it stays at its current value
                    // (either IBLTE for first iteration, or the known transition location
                    // from a previous iteration)
                    trans_loc = None;
                }
                TransitionResult::FreeTransition { ampl2, location } => {
                    ami = ampl2;
                    march.tran = true;
                    march.itran = ibl;
                    trans_loc = Some(location);
                }
                TransitionResult::ForcedTransition { location } => {
                    march.tran = true;
                    march.itran = ibl;
                    trans_loc = Some(location);
                }
            }

            // CRITICAL: Update s2.ampl with the value computed by TRCHEK
            // This ensures BLDIF uses the correct amplification when computing
            // the amplification equation residual. Without this, s2.ampl would
            // still contain the old value (s1.ampl) from the blprv call, causing
            // BLDIF to compute a non-zero residual that effectively doubles
            // the amplification accumulation.
            s2.ampl = ami;
        }

        // Determine flow type for this station
        let flow_type = if is_wake {
            BLFlowType::Wake
        } else if ibl >= march.itran {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };

        // Calculate secondary variables
        s2.blvar(flow_type, params);

        // For similarity station, "station 1" is really "station 2" (XFOIL BLSYS lines 627-631)
        let s1_for_bldif = if is_simi { &s2 } else { s1 };

        // Assemble local Newton system
        if let Some((cte, tte, dte)) = te_values {
            if is_first_wake {
                tesys(&mut local_sys, &s2, cte, tte, dte);
            } else if march.tran {
                // Transition interval: use TRDIF
                if let Some(ref trans) = trans_loc {
                    local_sys.trdif(s1, &s2, trans, march.acrit, params);
                } else {
                    // Fallback to turbulent if no transition location
                    let cfm = MidpointCf::compute(s1_for_bldif, &s2, flow_type, is_simi);
                    local_sys.bldif(s1_for_bldif, &s2, &cfm, flow_type, is_simi);
                }
            } else {
                // Calculate midpoint Cf
                let cfm = MidpointCf::compute(s1_for_bldif, &s2, flow_type, is_simi);
                local_sys.bldif(s1_for_bldif, &s2, &cfm, flow_type, is_simi);
            }
        } else if march.tran {
            // Transition interval: use TRDIF
            if let Some(ref trans) = trans_loc {
                local_sys.trdif(s1, &s2, trans, march.acrit, params);
            } else {
                // Fallback to turbulent if no transition location
                let cfm = MidpointCf::compute(s1_for_bldif, &s2, flow_type, is_simi);
                local_sys.bldif(s1_for_bldif, &s2, &cfm, flow_type, is_simi);
            }
        } else {
            let cfm = MidpointCf::compute(s1_for_bldif, &s2, flow_type, is_simi);
            local_sys.bldif(s1_for_bldif, &s2, &cfm, flow_type, is_simi);
        }

        // For similarity station, combine Jacobians: VS2 = VS1 + VS2, VS1 = 0 (XFOIL BLSYS lines 646-654)
        if is_simi {
            for k in 0..4 {
                for l in 0..5 {
                    local_sys.vs2[k][l] += local_sys.vs1[k][l];
                    local_sys.vs1[k][l] = 0.0;
                }
            }
        }

        // Set reference values on first iteration
        if itbl == 1 {
            ueref = s2.u;
            hkref = s2.hk;

            // If current point was turbulent but is now laminar, extrapolate Hk
            if ibl < march.itran && ibl >= march.itrold {
                // Use s1 to extrapolate baseline Hk
                hkref = s1.hk;
            }

            // If point was laminar but is now turbulent, reinit Ctau
            if ibl < march.itrold {
                if march.tran {
                    cti = 0.03;
                } else if march.turb {
                    cti = s1.ctau;
                }
                if march.tran || march.turb {
                    s2.ctau = cti;
                }
            }
        }

        // Set up 4th equation
        let (vs2_row4, vsrez4) = if is_simi || is_first_wake {
            // Prescribe Ue
            prescribed_ue_equation(&s2, ueref)
        } else {
            // Use Ue-Hk characteristic
            let (sens_new, row4, rez4) =
                calc_ue_hk_characteristic(&local_sys, &s2, hkref, ueref, sens, itbl);
            sens = sens_new;
            (row4, rez4)
        };

        // Copy row 4 into system
        for l in 0..5 {
            local_sys.vs2[3][l] = vs2_row4[l];
        }
        local_sys.vsrez[3] = vsrez4;

        // Extract 4x4 system for GAUSS
        let mut z = [[0.0f64; 4]; 4];
        let mut r = [0.0f64; 4];
        for k in 0..4 {
            r[k] = local_sys.vsrez[k];
            for l in 0..4 {
                z[k][l] = local_sys.vs2[k][l];
            }
        }

        // Solve Newton system
        gauss_solve_4x4(&mut z, &mut r);

        // Determine max change and underrelax if needed
        let mut dmax = (r[1] / thi).abs().max((r[2] / dsi).abs()).max((r[3] / uei).abs());
        if ibl >= march.itran {
            dmax = dmax.max((r[0] / (10.0 * cti)).abs());
        }

        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

        // Update variables
        if ibl < march.itran {
            ami += rlx * r[0];
        } else {
            cti += rlx * r[0];
        }
        thi += rlx * r[1];
        dsi += rlx * r[2];
        uei += rlx * r[3];

        // Clamp Ctau to reasonable range
        if ibl >= march.itran {
            cti = cti.min(0.30).max(0.0000001);
        }

        // Apply delta* limit
        let msq = s2.msq;
        let dsw = dsi - dswaki;
        let dsw_lim = dslim(dsw, thi, uei, msq, hklim);
        dsi = dsw_lim + dswaki;

        // Check convergence
        if dmax <= DEPS {
            // Final update of station state
            s2.blprv(xsi_i, ami, cti, thi, dsi, dswaki, uei, params);
            s2.blkin(params);
            s2.blvar(flow_type, params);

            // Store final values
            if ibl < march.itran {
                s2.ctau = ami; // Store amplification in ctau slot
            } else {
                s2.ctau = cti;
            }

            return (MarchResult::Converged, s2);
        }
    }

    // Failed to converge
    let dmax_final = (s2.theta / thi - 1.0)
        .abs()
        .max(((s2.dstar + s2.dw) / dsi - 1.0).abs())
        .max((s2.u / uei - 1.0).abs());

    if dmax_final <= 0.1 {
        // Partially converged - use current values
        s2.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, params);
        s2.blkin(params);
        let flow_type = if is_wake {
            BLFlowType::Wake
        } else if ibl >= march.itran {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        s2.blvar(flow_type, params);
        if ibl < march.itran {
            s2.ctau = ami;
        } else {
            s2.ctau = cti;
        }

        (MarchResult::PartiallyConverged { dmax: dmax_final }, s2)
    } else {
        // Station failed to converge (dmax > 0.1 after 25 iterations).
        // XFOIL implements fallback extrapolation here using:
        //   - Surface: sqrt(arc-length) scaling of theta and dstar
        //   - First wake point: trailing edge values
        //   - Wake continuation: blending formula with ratlen
        // See xfoil/xfoil6.99/src/xbl.f lines 1318-1336 (MRCHDU).
        //
        // Current status: All test cases converge adequately without this.
        // Implement only if high-alpha or extreme separation cases require it.
        (MarchResult::Failed { dmax: dmax_final }, s2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tesys_structure() {
        let mut sys = BLLocalSystem::default();
        let mut s2 = BLStationState::default();
        s2.ctau = 0.03;
        s2.theta = 0.001;
        s2.dstar = 0.002;
        s2.dw = 0.0001;

        let cte = 0.025;
        let tte = 0.0012;
        let dte = 0.0025;

        tesys(&mut sys, &s2, cte, tte, dte);

        // Check VS1 diagonal
        assert_eq!(sys.vs1[0][0], -1.0);
        assert_eq!(sys.vs1[1][1], -1.0);
        assert_eq!(sys.vs1[2][2], -1.0);

        // Check VS2 diagonal
        assert_eq!(sys.vs2[0][0], 1.0);
        assert_eq!(sys.vs2[1][1], 1.0);
        assert_eq!(sys.vs2[2][2], 1.0);

        // Check residuals
        assert!((sys.vsrez[0] - (cte - s2.ctau)).abs() < 1e-15);
        assert!((sys.vsrez[1] - (tte - s2.theta)).abs() < 1e-15);
        assert!((sys.vsrez[2] - (dte - s2.dstar - s2.dw)).abs() < 1e-15);
    }

    #[test]
    fn test_dslim() {
        let thi = 0.001;
        let uei = 1.0;
        let msq = 0.0;

        // Test that delta* below limit gets clamped
        let dsw_low = 1.01 * thi; // H = 1.01, Hk < 1.02
        let dsw_result = dslim(dsw_low, thi, uei, msq, 1.02);
        assert!(dsw_result >= dsw_low);

        // Test that delta* above limit is unchanged
        let dsw_high = 2.5 * thi; // H = 2.5, Hk > 1.02
        let dsw_result2 = dslim(dsw_high, thi, uei, msq, 1.02);
        assert_eq!(dsw_result2, dsw_high);
    }
}
