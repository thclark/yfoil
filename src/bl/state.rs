//! Boundary layer state definitions
//!
//! This module defines the data structures for boundary layer variables
//! at each station along the airfoil surface.

/// Flow regime indicator
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// Laminar flow
    Laminar,
    /// Turbulent flow
    Turbulent,
    /// Wake region (behind trailing edge)
    Wake,
}

/// Boundary layer state at a single station
///
/// Contains all BL variables needed for the integral formulation:
/// - Primary variables: θ, δ*, Cf, and amplification factor
/// - Derived quantities: H, H*, Ue, etc.
#[derive(Debug, Clone)]
pub struct BLStation {
    /// Arc length position along surface (from stagnation point)
    pub s: f64,

    /// x-coordinate
    pub x: f64,

    /// y-coordinate
    pub y: f64,

    /// Momentum thickness θ
    pub theta: f64,

    /// Displacement thickness δ*
    pub dstar: f64,

    /// Shape factor H = δ*/θ
    pub h: f64,

    /// Kinematic shape factor Hk (compressibility-corrected)
    pub hk: f64,

    /// Energy shape factor H*
    pub hs: f64,

    /// Skin friction coefficient Cf
    pub cf: f64,

    /// Edge velocity (normalized by freestream)
    pub ue: f64,

    /// Edge velocity squared for Mach correction
    pub ue_sq: f64,

    /// Momentum thickness Reynolds number Rθ = ue*θ/ν
    pub rtheta: f64,

    /// Amplification factor N (for eN transition)
    pub n_amp: f64,

    /// Dissipation coefficient CD
    pub cd: f64,

    /// Flow regime at this station
    pub regime: FlowRegime,

    /// Mass defect ρ*ue*δ* (for coupling)
    pub mass_defect: f64,
}

impl BLStation {
    /// Create a new station with laminar initialization
    pub fn new_laminar(s: f64, x: f64, y: f64, ue: f64) -> Self {
        // Typical laminar initialization near stagnation
        Self {
            s,
            x,
            y,
            theta: 0.0,
            dstar: 0.0,
            h: 2.5, // Typical laminar value
            hk: 2.5,
            hs: 1.55,
            cf: 0.0,
            ue,
            ue_sq: ue * ue,
            rtheta: 0.0,
            n_amp: 0.0,
            cd: 0.0,
            regime: FlowRegime::Laminar,
            mass_defect: 0.0,
        }
    }

    /// Create a wake station
    pub fn new_wake(s: f64, x: f64, y: f64, ue: f64) -> Self {
        Self {
            s,
            x,
            y,
            theta: 0.0,
            dstar: 0.0,
            h: 1.1, // Wake starts near H=1
            hk: 1.1,
            hs: 1.5,
            cf: 0.0, // No skin friction in wake
            ue,
            ue_sq: ue * ue,
            rtheta: 0.0,
            n_amp: 0.0,
            cd: 0.0,
            regime: FlowRegime::Wake,
            mass_defect: 0.0,
        }
    }

    /// Update derived quantities from primary variables
    pub fn update_derived(&mut self, msq: f64) {
        self.ue_sq = self.ue * self.ue;
        if self.theta > 0.0 {
            self.h = self.dstar / self.theta;
        }
        // Update kinematic shape factor
        let (hk, _, _) = crate::bl::hkin(self.h, msq);
        self.hk = hk;
    }
}

/// Complete boundary layer solution along one side of airfoil
#[derive(Debug, Clone)]
pub struct BLSide {
    /// Stations from stagnation point to trailing edge
    pub stations: Vec<BLStation>,

    /// Index of transition location (None if fully laminar)
    pub transition_index: Option<usize>,

    /// Arc length at transition
    pub s_transition: f64,

    /// Total number of stations
    pub n_stations: usize,
}

impl BLSide {
    /// Create a new BL side with given number of stations
    pub fn new(n_stations: usize) -> Self {
        Self {
            stations: Vec::with_capacity(n_stations),
            transition_index: None,
            s_transition: f64::INFINITY,
            n_stations,
        }
    }

    /// Add a station
    pub fn push(&mut self, station: BLStation) {
        self.stations.push(station);
    }

    /// Get station at index
    pub fn get(&self, i: usize) -> Option<&BLStation> {
        self.stations.get(i)
    }

    /// Get mutable station at index
    pub fn get_mut(&mut self, i: usize) -> Option<&mut BLStation> {
        self.stations.get_mut(i)
    }

    /// Number of stations
    pub fn len(&self) -> usize {
        self.stations.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.stations.is_empty()
    }

    /// Mark transition at given index
    pub fn set_transition(&mut self, index: usize) {
        self.transition_index = Some(index);
        if let Some(station) = self.stations.get(index) {
            self.s_transition = station.s;
        }
        // Mark all subsequent stations as turbulent
        for i in index..self.stations.len() {
            if let Some(station) = self.stations.get_mut(i) {
                station.regime = FlowRegime::Turbulent;
            }
        }
    }
}

/// Complete boundary layer state for the airfoil
#[derive(Debug, Clone)]
pub struct BLState {
    /// Upper surface BL (from stagnation to TE)
    pub upper: BLSide,

    /// Lower surface BL (from stagnation to TE)
    pub lower: BLSide,

    /// Wake BL (from TE downstream)
    pub wake: BLSide,

    /// Stagnation point arc length
    pub s_stag: f64,

    /// Stagnation point index on panel array
    pub i_stag: usize,

    /// Total drag coefficient (pressure + friction)
    pub cd: f64,

    /// Friction drag coefficient
    pub cdf: f64,

    /// Pressure drag coefficient
    pub cdp: f64,
}

impl BLState {
    /// Create a new empty BL state
    pub fn new() -> Self {
        Self {
            upper: BLSide::new(0),
            lower: BLSide::new(0),
            wake: BLSide::new(0),
            s_stag: 0.0,
            i_stag: 0,
            cd: 0.0,
            cdf: 0.0,
            cdp: 0.0,
        }
    }

    /// Initialize BL state for a given airfoil geometry
    pub fn initialize(n_upper: usize, n_lower: usize, n_wake: usize) -> Self {
        Self {
            upper: BLSide::new(n_upper),
            lower: BLSide::new(n_lower),
            wake: BLSide::new(n_wake),
            s_stag: 0.0,
            i_stag: 0,
            cd: 0.0,
            cdf: 0.0,
            cdp: 0.0,
        }
    }
}

impl Default for BLState {
    fn default() -> Self {
        Self::new()
    }
}

/// Flow conditions for BL calculation
#[derive(Debug, Clone, Copy)]
pub struct FlowConditions {
    /// Freestream Mach number
    pub mach: f64,

    /// Mach number squared (cached for efficiency)
    pub msq: f64,

    /// Reynolds number based on chord
    pub reynolds: f64,

    /// Critical amplification factor for transition
    pub ncrit: f64,

    /// Kinematic viscosity (computed from Re)
    pub nu: f64,
}

impl FlowConditions {
    /// Create flow conditions from basic parameters
    ///
    /// # Arguments
    /// * `reynolds` - Reynolds number based on chord
    /// * `mach` - Mach number
    /// * `ncrit` - Critical amplification factor (typically 9)
    /// * `chord` - Airfoil chord length
    pub fn new(reynolds: f64, mach: f64, ncrit: f64, chord: f64) -> Self {
        Self {
            mach,
            msq: mach * mach,
            reynolds,
            ncrit,
            nu: chord / reynolds, // ν = c/Re (assuming unit freestream velocity)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bl_station_new_laminar() {
        let station = BLStation::new_laminar(0.1, 0.5, 0.05, 1.2);
        assert_eq!(station.s, 0.1);
        assert_eq!(station.x, 0.5);
        assert_eq!(station.ue, 1.2);
        assert_eq!(station.regime, FlowRegime::Laminar);
    }

    #[test]
    fn test_bl_side_transition() {
        let mut side = BLSide::new(5);
        for i in 0..5 {
            side.push(BLStation::new_laminar(i as f64 * 0.1, 0.0, 0.0, 1.0));
        }
        assert_eq!(side.len(), 5);

        // Set transition at index 2
        side.set_transition(2);
        assert_eq!(side.transition_index, Some(2));
        assert_eq!(side.stations[0].regime, FlowRegime::Laminar);
        assert_eq!(side.stations[1].regime, FlowRegime::Laminar);
        assert_eq!(side.stations[2].regime, FlowRegime::Turbulent);
        assert_eq!(side.stations[3].regime, FlowRegime::Turbulent);
        assert_eq!(side.stations[4].regime, FlowRegime::Turbulent);
    }

    #[test]
    fn test_flow_conditions() {
        let cond = FlowConditions::new(1_000_000.0, 0.3, 9.0, 1.0);
        assert_eq!(cond.reynolds, 1_000_000.0);
        assert_eq!(cond.mach, 0.3);
        assert_eq!(cond.msq, 0.09);
        assert_eq!(cond.ncrit, 9.0);
        assert_eq!(cond.nu, 1e-6);
    }
}
