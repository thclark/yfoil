//! Foil output: the geometry, the wake and every per-station boundary-layer quantity of an
//! operating point, plus the significant locations (stagnation, transition, separation), in a
//! serialisable form that `yfoil plot foil` consumes.
//!
//! # Live versus stored closure quantities
//!
//! XFOIL's SETBL calls MRCHDU at the top of *every* Newton iteration (`xbl.f:93`); MRCHDU and
//! SETBL store `TAU, DIS, CTQ, DELT, USLP, TSTR` from the state *entering* that iteration
//! (`xbl.f:277-282, 1157-1167`). BLSOLV/UPDATE then correct the primaries `THET, DSTR, UEDG, CTAU,
//! MASS`, and nothing refreshes the closure arrays afterwards. So XFOIL's `DUMP` prints `Ue, Dstar,
//! Theta, H, HK` from the converged state but `Cf, H*, K, tau, Di` one Newton correction behind
//! (`H* = TSTR/THET` and `K = TSTR*Ue^3` mix a lagged numerator with a current denominator), and the
//! VPLO `CF`, `CD` and `DELT` plots are lagged the same way. The mismatch is the size of the final
//! Newton correction (bounded by the RMSBL < 1e-4 convergence test on a converged point; unbounded
//! on an unconverged one).
//!
//! YFoil reproduces those stored arrays exactly (they are fixture-gated). The `closures` of
//! [`SideStations`] are therefore computed *live* here by running XFOIL's own BLPRV → BLKIN → BLVAR
//! on the converged `primaries`, and the lagged arrays are emitted verbatim under
//! [`SideStations::lagged_closures`] only when asked for (`--include-lagged-closures`), so a
//! column-for-column comparison with an XFOIL `DUMP` stays possible. A side-by-side of the live
//! `hstar` against DUMP's `H*` differs at that ~1e-4 level by construction.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

use crate::bl::system::{FlowParameters, FlowRegime, StationState};
use crate::geometry::{spline_value, PanelledFoil};
use crate::solver::blstate::SolverState;

// ============================================================================
// Geometry
// ============================================================================

/// Wake node geometry: XYWAKE's nodes `N+1..=N+NW` and their normals. The wake normal is
/// `-grad(psi)/|grad(psi)|` (`xpanel.f:1311-1330`), which points towards the *lower* side.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WakeNodes {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Arc length continued from the airfoil's S(N)
    pub s: Vec<f64>,
    pub normal_x: Vec<f64>,
    pub normal_y: Vec<f64>,
}

/// Airfoil node geometry as the solver holds it (XFOIL's `X, Y, S, NX, NY` for `I = 1..N`), with
/// the TECALC scalars and, after a viscous solve, the wake.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FoilNodes {
    /// Node coordinates, TE → upper → LE → lower → TE
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Spline arc length at each node
    pub s: Vec<f64>,
    /// Outward unit normals at the nodes (NCALC: spline normals, corner-averaged at doubled nodes)
    pub normal_x: Vec<f64>,
    pub normal_y: Vec<f64>,
    pub chord: f64,
    pub x_le: f64,
    pub y_le: f64,
    pub x_te: f64,
    pub y_te: f64,
    /// Arc length of the leading edge
    pub s_le: f64,
    /// The node nearest the spline leading edge (0-based; splits the upper and lower surfaces)
    pub i_le_node: usize,
    /// TE gap area projected normal to the TE bisector (TECALC's ANTE)
    pub te_thickness_normal: f64,
    pub sharp_te: bool,
    /// Present once XYWAKE has run (any viscous point)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wake: Option<WakeNodes>,
}

impl FoilNodes {
    /// Airfoil nodes `1..=n` and, if `st.lwake`, wake nodes `n+1..=n+nw`.
    pub fn from_state(state: &SolverState) -> Self {
        let n = state.n_foil_nodes;
        let nw = state.n_wake_nodes;
        let wake = if state.wake_built && nw > 0 {
            Some(WakeNodes {
                x: state.x[n + 1..=n + nw].to_vec(),
                y: state.y[n + 1..=n + nw].to_vec(),
                s: state.s[n + 1..=n + nw].to_vec(),
                normal_x: state.normal_x[n + 1..=n + nw].to_vec(),
                normal_y: state.normal_y[n + 1..=n + nw].to_vec(),
            })
        } else {
            None
        };
        Self {
            x: state.x[1..=n].to_vec(),
            y: state.y[1..=n].to_vec(),
            s: state.s[1..=n].to_vec(),
            normal_x: state.normal_x[1..=n].to_vec(),
            normal_y: state.normal_y[1..=n].to_vec(),
            chord: state.chord,
            x_le: state.x_le,
            y_le: state.y_le,
            x_te: state.x_te,
            y_te: state.y_te,
            s_le: state.s_le,
            i_le_node: (1..=n)
                .min_by(|&i, &j| {
                    (state.s[i] - state.s_le)
                        .abs()
                        .partial_cmp(&(state.s[j] - state.s_le).abs())
                        .unwrap()
                })
                .map(|i| i - 1)
                .unwrap_or(0),
            te_thickness_normal: state.te_thickness_normal,
            sharp_te: state.sharp_te,
            wake,
        }
    }

    /// Geometry only (no wake), for plotting panelings without a solve.
    pub fn from_panelled(airfoil: &PanelledFoil) -> Self {
        Self::from_state(&SolverState::from_foil(airfoil, 0))
    }

    /// Number of airfoil nodes
    pub fn n(&self) -> usize {
        self.x.len()
    }

    /// Bitwise comparison of the panel nodes (what "same geometry" means for overlaying design
    /// points: identical panels, so the BL stations coincide).
    pub fn same_panels(&self, other: &Self) -> bool {
        self.x.len() == other.x.len()
            && self.x.iter().zip(&other.x).all(|(a, b)| a.to_bits() == b.to_bits())
            && self.y.iter().zip(&other.y).all(|(a, b)| a.to_bits() == b.to_bits())
    }
}

// ============================================================================
// BL quantities
// ============================================================================

/// A per-station boundary-layer quantity that can be plotted normal to the surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BlQuantity {
    /// Displacement thickness δ* (DSTR; in the wake the total including the wake gap)
    Dstar,
    /// Momentum thickness θ (THET)
    Theta,
    /// Boundary-layer thickness δ (Green's correlation, BLVAR's DE)
    Delta,
    /// Shape factor H = δ*/θ
    H,
    /// Kinematic shape factor Hk (HKIN; unclamped, as DUMP prints it)
    Hk,
    /// Kinetic-energy shape factor H* (live BLVAR value; see module docs for DUMP's lagged `H*`)
    Hstar,
    /// Edge velocity Ue/Vinf after the Karman–Tsien transformation (unsigned on both sides)
    Ue,
    /// Skin-friction coefficient Cf = τ/(½ρ∞V∞²) (zero in the wake)
    Cf,
    /// Dissipation coefficient CD = Di/(ρ∞V∞³) (DUMP's `CDIS`)
    Cdiss,
    /// Shear-stress coefficient √Cτ (turbulent) or amplification N (laminar): CTAU
    Sqrtctau,
    /// Equilibrium shear-stress coefficient √Cτ_eq: CTQ
    Sqrtctaueq,
    /// Slip-velocity parameter 1.6/(1+Us): USLP
    Us,
    /// Viscous pressure coefficient at the node: CPV
    Cp,
    /// Mass defect m = Ue·δ*: MASS
    MassDefect,
}

impl BlQuantity {
    /// Every quantity, in the order the CLI documents them
    pub const ALL: [BlQuantity; 14] = [
        BlQuantity::Dstar,
        BlQuantity::Theta,
        BlQuantity::Delta,
        BlQuantity::H,
        BlQuantity::Hk,
        BlQuantity::Hstar,
        BlQuantity::Ue,
        BlQuantity::Cf,
        BlQuantity::Cdiss,
        BlQuantity::Sqrtctau,
        BlQuantity::Sqrtctaueq,
        BlQuantity::Us,
        BlQuantity::Cp,
        BlQuantity::MassDefect,
    ];

    /// CLI / JSON name
    pub fn name(self) -> &'static str {
        match self {
            BlQuantity::Dstar => "dstar",
            BlQuantity::Theta => "theta",
            BlQuantity::Delta => "delta",
            BlQuantity::H => "h",
            BlQuantity::Hk => "hk",
            BlQuantity::Hstar => "hstar",
            BlQuantity::Ue => "ue",
            BlQuantity::Cf => "cf",
            BlQuantity::Cdiss => "cdiss",
            BlQuantity::Sqrtctau => "sqrtctau",
            BlQuantity::Sqrtctaueq => "sqrtctaueq",
            BlQuantity::Us => "us",
            BlQuantity::Cp => "cp",
            BlQuantity::MassDefect => "mass_defect",
        }
    }

    /// Legend label
    pub fn label(self) -> &'static str {
        match self {
            BlQuantity::Dstar => "δ*",
            BlQuantity::Theta => "θ",
            BlQuantity::Delta => "δ",
            BlQuantity::H => "H",
            BlQuantity::Hk => "Hk",
            BlQuantity::Hstar => "H*",
            BlQuantity::Ue => "Ue/V∞",
            BlQuantity::Cf => "Cf",
            BlQuantity::Cdiss => "CD",
            BlQuantity::Sqrtctau => "√Cτ / N",
            BlQuantity::Sqrtctaueq => "√Cτ_eq",
            BlQuantity::Us => "1.6/(1+Us)",
            BlQuantity::Cp => "Cp",
            BlQuantity::MassDefect => "m",
        }
    }

    /// True for the quantities that are lengths in chord units (δ*, θ, δ), so that a scale
    /// factor of 1 draws them at true geometric size.
    pub fn is_length(self) -> bool {
        matches!(self, BlQuantity::Dstar | BlQuantity::Theta | BlQuantity::Delta)
    }
}

impl fmt::Display for BlQuantity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for BlQuantity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let key = s.trim().to_ascii_lowercase();
        BlQuantity::ALL
            .iter()
            .copied()
            .find(|q| q.name() == key)
            .ok_or_else(|| {
                let names: Vec<&str> = BlQuantity::ALL.iter().map(|q| q.name()).collect();
                format!("unknown BL quantity '{}' (expected one of: {})", s, names.join(", "))
            })
    }
}

/// XFOIL's closure arrays exactly as stored after the final SETBL/MRCHDU call, i.e. evaluated on
/// the state entering the last Newton iteration and therefore one Newton correction behind the
/// primaries (see the module docs). These are what XFOIL's DUMP `Cf`, `H*`, `K`, `tau`, `Di` columns
/// and its VPLO `CF`/`CD`/`DELT` plots show.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct LaggedClosures {
    /// TAU = ½·ρ·Ue²·Cf
    pub tau: Vec<f64>,
    /// DIS = ½·ρ·Ue³·Di·H*
    pub dissipation: Vec<f64>,
    /// CTQ
    pub sqrtctaueq: Vec<f64>,
    /// DELT
    pub delta: Vec<f64>,
    /// USLP
    pub us_plot_scale: Vec<f64>,
    /// TSTR = H*·θ from the march
    pub thetastar: Vec<f64>,
    /// DUMP's `H*` = TSTR/THET (1.0 where THET == 0, as XFOIL)
    pub hstar_dump: Vec<f64>,
    /// DUMP's `Cf` = TAU/(½QINF²)
    pub cf_dump: Vec<f64>,
}

/// The converged primary variables of one side, one entry per station (the solver's state
/// after the last UPDATE)
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Primaries {
    /// UEDG: incompressible edge velocity
    pub ue: Vec<f64>,
    /// THET
    pub theta: Vec<f64>,
    /// DSTR (in the wake: total, including the wake gap WGAP)
    pub dstar: Vec<f64>,
    /// CTAU (Cτ^½ turbulent, amplification N laminar)
    pub sqrtctau: Vec<f64>,
    /// MASS
    pub mass_defect: Vec<f64>,
}

/// The closure quantities evaluated live on the converged primaries (BLPRV → BLKIN → BLVAR)
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Closures {
    /// Ue/q∞ after the Karman–Tsien transformation (unsigned; DUMP signs it by GAM)
    pub ue_compressible: Vec<f64>,
    /// H = δ*/θ (wake: without the wake gap, as BLKIN)
    pub h: Vec<f64>,
    /// Hk from HKIN (unclamped)
    pub hk: Vec<f64>,
    /// H* from BLVAR
    pub hstar: Vec<f64>,
    /// Cf = ρ·Ue²·cf/q∞² (BLVAR's cf scaled as SETBL forms TAU, then DUMP's normalisation)
    pub cf: Vec<f64>,
    /// CD = DIS/q∞³ with DIS as SETBL forms it
    pub cdiss: Vec<f64>,
    /// δ from BLVAR (DE)
    pub delta: Vec<f64>,
    /// Cτ_eq^½ from BLVAR (CQ)
    pub sqrtctaueq: Vec<f64>,
    /// Us/Ue from BLVAR (the normalised wall-slip velocity; XFOIL's USLP is 1.6/(1+Us))
    pub us: Vec<f64>,
    /// Rθ from BLKIN
    pub retheta: Vec<f64>,
    /// Edge Mach² from BLKIN
    pub machsqd_edge: Vec<f64>,
}

/// One side's stations (upper, lower or wake) as a struct of arrays; every vector has one entry
/// per station, in marching order.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct SideStations {
    /// IBL (1-based station index on this side, as XFOIL)
    pub i_station: Vec<usize>,
    /// IPAN (1-based panel/wake node index, as XFOIL)
    pub i_node: Vec<usize>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// BL arc coordinate ξ from the stagnation point (XSSI)
    pub xi: Vec<f64>,
    /// CPV at the station's node
    pub cp: Vec<f64>,
    pub primaries: Primaries,
    pub closures: Closures,
    /// XFOIL's lagged closure arrays, only with `--include-lagged-closures`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lagged_closures: Option<LaggedClosures>,
}

impl SideStations {
    /// Number of stations
    pub fn len(&self) -> usize {
        self.i_station.len()
    }

    /// True when the side has no stations
    pub fn is_empty(&self) -> bool {
        self.i_station.is_empty()
    }

    /// The column for a quantity
    pub fn column(&self, q: BlQuantity) -> &[f64] {
        match q {
            BlQuantity::Dstar => &self.primaries.dstar,
            BlQuantity::Theta => &self.primaries.theta,
            BlQuantity::Delta => &self.closures.delta,
            BlQuantity::H => &self.closures.h,
            BlQuantity::Hk => &self.closures.hk,
            BlQuantity::Hstar => &self.closures.hstar,
            BlQuantity::Ue => &self.closures.ue_compressible,
            BlQuantity::Cf => &self.closures.cf,
            BlQuantity::Cdiss => &self.closures.cdiss,
            BlQuantity::Sqrtctau => &self.primaries.sqrtctau,
            BlQuantity::Sqrtctaueq => &self.closures.sqrtctaueq,
            BlQuantity::Us => &self.closures.us,
            BlQuantity::Cp => &self.cp,
            BlQuantity::MassDefect => &self.primaries.mass_defect,
        }
    }

    /// Every column has exactly one entry per station
    pub fn check_lengths(&self) -> Result<(), String> {
        let n = self.len();
        let (p, c) = (&self.primaries, &self.closures);
        let mut cols: Vec<(&str, usize)> = vec![
            ("i_node", self.i_node.len()),
            ("x", self.x.len()),
            ("y", self.y.len()),
            ("xi", self.xi.len()),
            ("cp", self.cp.len()),
            ("primaries.ue", p.ue.len()),
            ("primaries.theta", p.theta.len()),
            ("primaries.dstar", p.dstar.len()),
            ("primaries.sqrtctau", p.sqrtctau.len()),
            ("primaries.mass_defect", p.mass_defect.len()),
            ("closures.ue_compressible", c.ue_compressible.len()),
            ("closures.h", c.h.len()),
            ("closures.hk", c.hk.len()),
            ("closures.hstar", c.hstar.len()),
            ("closures.cf", c.cf.len()),
            ("closures.cdiss", c.cdiss.len()),
            ("closures.delta", c.delta.len()),
            ("closures.sqrtctaueq", c.sqrtctaueq.len()),
            ("closures.us", c.us.len()),
            ("closures.retheta", c.retheta.len()),
            ("closures.machsqd_edge", c.machsqd_edge.len()),
        ];
        if let Some(l) = &self.lagged_closures {
            cols.extend([
                ("lagged_closures.tau", l.tau.len()),
                ("lagged_closures.dissipation", l.dissipation.len()),
                ("lagged_closures.sqrtctaueq", l.sqrtctaueq.len()),
                ("lagged_closures.delta", l.delta.len()),
                ("lagged_closures.us_plot_scale", l.us_plot_scale.len()),
                ("lagged_closures.thetastar", l.thetastar.len()),
                ("lagged_closures.hstar_dump", l.hstar_dump.len()),
                ("lagged_closures.cf_dump", l.cf_dump.len()),
            ]);
        }
        for (name, len) in cols {
            if len != n {
                return Err(format!("column {name} has {len} entries, expected {n}"));
            }
        }
        Ok(())
    }
}

// ============================================================================
// Significant locations
// ============================================================================

/// The stagnation point (STFIND): arc length SST and its position on the spline
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StagnationMarker {
    /// IST: the node ahead of the stagnation point
    pub i_stagnation_node: usize,
    /// SST: spline arc length of the stagnation point
    pub s_stagnation: f64,
    pub x: f64,
    pub y: f64,
}

/// Transition on one side
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TransitionMarker {
    /// ITRAN: first turbulent station
    pub i_station: usize,
    /// TFORCE: transition was forced (X_TRIP)
    pub forced: bool,
    /// XOCTR, YOCTR: the transition point as chord fractions
    pub x_transition: f64,
    pub y_transition: f64,
    /// Foil arc coordinate of the transition point (SST ∓ XSSITR)
    pub s_transition: f64,
}

/// Direction of a wall-shear sign change
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeparationKind {
    /// Cf goes from ≥ 0 to < 0 marching downstream
    Separation,
    /// Cf returns to ≥ 0
    Reattachment,
}

/// A separation or reattachment point. **Derived by YFoil for plotting**: XFOIL reports no such
/// location anywhere. It is the sign change of the live `cf` between two consecutive surface
/// stations, interpolated linearly in arc length and placed on the spline.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SeparationMarker {
    /// 1 = upper, 2 = lower
    pub side: usize,
    pub kind: SeparationKind,
    /// IBL of the last station before the sign change
    pub station_before: usize,
    pub s: f64,
    pub x: f64,
    pub y: f64,
}

/// The boundary layer of one operating point: both sides, the wake, and the markers
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundaryLayerOutput {
    pub upper: SideStations,
    pub lower: SideStations,
    pub wake: SideStations,
    /// IBLTE(1..2): last airfoil station on each side
    pub i_te_station: [usize; 2],
    /// ITRAN(1..2)
    pub i_transition_station: [usize; 2],
    /// NW
    pub n_wake_nodes: usize,
    /// QINF
    pub qinf: f64,
    /// ANTE (TECALC)
    pub te_thickness_normal: f64,
    pub stagnation: StagnationMarker,
    /// Upper, lower
    pub transition: [TransitionMarker; 2],
    /// Derived (see [`SeparationMarker`]); empty when Cf never changes sign
    pub derived_separation: Vec<SeparationMarker>,
    /// CPDISP's DSF1, DSF2: the fraction of the wake δ* drawn on the upper and lower side of the
    /// wake, `(DSTR(IBLTE(is),is) + ½·ANTE) / DSTR(IBLTE(2)+1, 2)`, or 0.5/0.5 when the first wake
    /// δ* is exactly zero (`xplots.f:714-721`)
    pub wake_split: [f64; 2],
}

/// Live closure values at one station
struct Live {
    ue_compressible: f64,
    h: f64,
    hk: f64,
    hstar: f64,
    cf: f64,
    cdiss: f64,
    delta: f64,
    sqrtctaueq: f64,
    us: f64,
    retheta: f64,
    machsqd_edge: f64,
}

/// BLPRV → BLKIN → BLVAR on the converged primaries of station (is, ibl), with the flow type
/// SETBL would use there (laminar ahead of ITRAN, turbulent from it, wake past IBLTE).
fn live_closures(state: &SolverState, params: &FlowParameters, side: usize, i_station: usize) -> Live {
    let wake = i_station > state.i_te_station[side];
    let turb = i_station >= state.i_transition_station[side];
    let ctau = state.sqrtctau[side][i_station];
    let (ami, cti) = if turb { (0.0, ctau) } else { (ctau, 0.0) };
    let uei = state.ue[side][i_station];
    let thi = state.theta[side][i_station];
    if thi == 0.0 || uei == 0.0 {
        // an unsolved station (DUMP prints H = H* = 1 there)
        return Live {
            ue_compressible: 0.0,
            h: 1.0,
            hk: 1.0,
            hstar: 1.0,
            cf: 0.0,
            cdiss: 0.0,
            delta: 0.0,
            sqrtctaueq: 0.0,
            us: 0.0,
            retheta: 0.0,
            machsqd_edge: 0.0,
        };
    }
    // DSI = MDI/UEI, DSWAKI = WGAP(IW), as SETBL sets the "2" station
    let dsi = state.mass_defect[side][i_station] / uei;
    let dswaki = if wake {
        state.wake_gap[i_station - state.i_te_station[side]]
    } else {
        0.0
    };

    let mut s = StationState::default();
    s.set_primary_variables(state.xi[side][i_station], ami, cti, thi, dsi, dswaki, uei, params);
    s.set_kinematic_variables(params);
    let (h, hk, retheta, machsqd_edge) = (s.h, s.hk, s.retheta, s.machsqd_edge);
    let flow = if wake {
        FlowRegime::Wake
    } else if turb {
        FlowRegime::Turbulent
    } else {
        FlowRegime::Laminar
    };
    s.set_closure_variables(flow, params);

    let qinf = state.qinf;
    Live {
        ue_compressible: s.ue / qinf,
        h,
        hk,
        hstar: s.hstar,
        // TAU = ½·R2·U2²·CF2 (SETBL), Cf = TAU/(½·QINF²) (DUMP)
        cf: s.rho * s.ue * s.ue * s.cf / (qinf * qinf),
        // DIS = R2·U2³·DI2·HS2·½ (SETBL), CDIS = DIS/QINF³ (DUMP)
        cdiss: s.rho * s.ue * s.ue * s.ue * s.cdiss * s.hstar * 0.5 / (qinf * qinf * qinf),
        delta: s.delta,
        sqrtctaueq: s.sqrtctaueq,
        us: s.us,
        retheta,
        machsqd_edge,
    }
}

fn side_output(
    state: &SolverState,
    params: &FlowParameters,
    side: usize,
    stations: impl Iterator<Item = usize>,
    include_lagged_closures: bool,
) -> SideStations {
    let mut out = SideStations::default();
    let mut lagged = LaggedClosures::default();
    let qinf = state.qinf;
    for i_station in stations {
        let i = state.i_node[side][i_station];
        let live = live_closures(state, params, side, i_station);
        let theta = state.theta[side][i_station];
        out.i_station.push(i_station);
        out.i_node.push(i);
        out.x.push(state.x[i]);
        out.y.push(state.y[i]);
        out.xi.push(state.xi[side][i_station]);
        out.cp.push(state.cp_viscous[i]);
        let p = &mut out.primaries;
        p.ue.push(state.ue[side][i_station]);
        p.theta.push(theta);
        p.dstar.push(state.dstar[side][i_station]);
        p.sqrtctau.push(state.sqrtctau[side][i_station]);
        p.mass_defect.push(state.mass_defect[side][i_station]);
        let c = &mut out.closures;
        c.ue_compressible.push(live.ue_compressible);
        c.h.push(live.h);
        c.hk.push(live.hk);
        c.hstar.push(live.hstar);
        c.cf.push(live.cf);
        c.cdiss.push(live.cdiss);
        c.delta.push(live.delta);
        c.sqrtctaueq.push(live.sqrtctaueq);
        c.us.push(live.us);
        c.retheta.push(live.retheta);
        c.machsqd_edge.push(live.machsqd_edge);
        if include_lagged_closures {
            lagged.tau.push(state.tau[side][i_station]);
            lagged.dissipation.push(state.dissipation[side][i_station]);
            lagged.sqrtctaueq.push(state.sqrtctaueq[side][i_station]);
            lagged.delta.push(state.delta[side][i_station]);
            lagged.us_plot_scale.push(state.us_plot_scale[side][i_station]);
            lagged.thetastar.push(state.thetastar[side][i_station]);
            // BLDUMP: IF(TH.EQ.0.0) HS = 1.0 ELSE HS = TS/TH
            lagged.hstar_dump.push(if theta == 0.0 {
                1.0
            } else {
                state.thetastar[side][i_station] / theta
            });
            lagged.cf_dump.push(state.tau[side][i_station] / (0.5 * qinf * qinf));
        }
    }
    if include_lagged_closures {
        out.lagged_closures = Some(lagged);
    }
    out
}

/// Position on the airfoil spline at arc length `s`
fn spline_point(state: &SolverState, s: f64) -> (f64, f64) {
    let n = state.n_foil_nodes;
    let x = spline_value(s, &state.x[1..=n], &state.dxds[1..=n], &state.s[1..=n]);
    let y = spline_value(s, &state.y[1..=n], &state.dyds[1..=n], &state.s[1..=n]);
    (x, y)
}

/// Sign changes of the live `cf` along one side's surface stations, from the pair (2, 3) on
fn separation_markers(state: &SolverState, stations: &SideStations, side: usize) -> Vec<SeparationMarker> {
    let mut out = Vec::new();
    for k in 1..stations.len() {
        let (prev, cur) = (stations.closures.cf[k - 1], stations.closures.cf[k]);
        let kind = if prev >= 0.0 && cur < 0.0 {
            SeparationKind::Separation
        } else if prev < 0.0 && cur >= 0.0 {
            SeparationKind::Reattachment
        } else {
            continue;
        };
        let (s0, s1) = (state.s[stations.i_node[k - 1]], state.s[stations.i_node[k]]);
        let frac = if cur == prev { 0.0 } else { prev / (prev - cur) };
        let s = s0 + frac * (s1 - s0);
        let (x, y) = spline_point(state, s);
        out.push(SeparationMarker {
            side,
            kind,
            station_before: stations.i_station[k - 1],
            s,
            x,
            y,
        });
    }
    out
}

impl BoundaryLayerOutput {
    /// Every station of a solved state. Requires the pointer layer and BL arrays (any state after
    /// a viscous VISCAL call, converged or not). `include_lagged_closures` adds XFOIL's lagged
    /// closure arrays to each side.
    pub fn from_state(state: &SolverState, include_lagged_closures: bool) -> Self {
        let params = FlowParameters::new(state.mach, state.re, state.gamma_gas);
        let lag = include_lagged_closures;
        let upper = side_output(state, &params, 1, 2..=state.i_te_station[1], lag);
        let lower = side_output(state, &params, 2, 2..=state.i_te_station[2], lag);
        let wake = side_output(state, &params, 2, state.i_te_station[2] + 1..=state.n_stations[2], lag);

        let (sx, sy) = spline_point(state, state.s_stagnation);
        let stagnation = StagnationMarker {
            i_stagnation_node: state.i_stagnation_node,
            s_stagnation: state.s_stagnation,
            x: sx,
            y: sy,
        };
        let transition_marker = |side: usize| {
            let s = if side == 1 {
                state.s_stagnation - state.xi_transition[side]
            } else {
                state.s_stagnation + state.xi_transition[side]
            };
            TransitionMarker {
                i_station: state.i_transition_station[side],
                forced: state.transition_forced[side],
                x_transition: state.x_transition[side],
                y_transition: state.y_transition[side],
                s_transition: s,
            }
        };
        let transition = [transition_marker(1), transition_marker(2)];

        let mut derived_separation = separation_markers(state, &upper, 1);
        derived_separation.extend(separation_markers(state, &lower, 2));

        // CPDISP: upper and lower wake Dstar fractions from the first wake point
        let dstrte = state.dstar[2][state.i_te_station[2] + 1];
        let wake_split = if dstrte != 0.0 {
            [
                (state.dstar[1][state.i_te_station[1]] + 0.5 * state.te_thickness_normal) / dstrte,
                (state.dstar[2][state.i_te_station[2]] + 0.5 * state.te_thickness_normal) / dstrte,
            ]
        } else {
            [0.5, 0.5]
        };

        Self {
            upper,
            lower,
            wake,
            i_te_station: [state.i_te_station[1], state.i_te_station[2]],
            i_transition_station: [state.i_transition_station[1], state.i_transition_station[2]],
            n_wake_nodes: state.n_wake_nodes,
            qinf: state.qinf,
            te_thickness_normal: state.te_thickness_normal,
            stagnation,
            transition,
            derived_separation,
            wake_split,
        }
    }

    /// The largest |value| of a quantity over both sides and the wake
    pub fn max_abs(&self, q: BlQuantity) -> f64 {
        [&self.upper, &self.lower, &self.wake]
            .iter()
            .flat_map(|s| s.column(q).iter().copied())
            .filter(|v| v.is_finite())
            .fold(0.0_f64, |m, v| m.max(v.abs()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantity_names_round_trip() {
        for q in BlQuantity::ALL {
            assert_eq!(q.name().parse::<BlQuantity>().unwrap(), q);
            assert_eq!(q.name().to_ascii_uppercase().parse::<BlQuantity>().unwrap(), q);
            let json = serde_json::to_string(&q).unwrap();
            assert_eq!(json, format!("\"{}\"", q.name()));
        }
        assert!("dstar,theta".parse::<BlQuantity>().is_err());
        assert!("nope".parse::<BlQuantity>().unwrap_err().contains("dstar"));
    }

    #[test]
    fn geometry_from_paneled_has_no_wake_and_matches_nodes() {
        use crate::geometry::{naca_4digit, panel_foil, Thickness};
        let geom = naca_4digit("2412", 100, Thickness::Perpendicular).unwrap();
        let airfoil = panel_foil(&geom);
        let g = FoilNodes::from_panelled(&airfoil);
        assert_eq!(g.n(), 100);
        assert!(g.wake.is_none());
        assert_eq!(g.x, airfoil.x);
        assert_eq!(g.normal_x, airfoil.normal_x);
        assert_eq!(g.sharp_te, airfoil.sharp_te);
        assert!(g.same_panels(&g.clone()));
        let mut h = g.clone();
        h.y[10] = f64::from_bits(h.y[10].to_bits() + 1);
        assert!(!g.same_panels(&h));
    }
}
