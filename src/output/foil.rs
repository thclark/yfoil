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
//! YFoil reproduces those stored arrays exactly (they are fixture-gated). The canonical columns of
//! [`BlSideOutput`] are therefore computed *live* here by running XFOIL's own BLPRV → BLKIN → BLVAR
//! on the converged primaries, and the lagged arrays are also emitted verbatim under
//! [`BlSideOutput::stored`] so a column-for-column comparison with an XFOIL `DUMP` is possible. A
//! side-by-side of the live `hs` against DUMP's `H*` differs at that ~1e-4 level by construction.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

use crate::bl::system::{BLStationState, FlowParameters, FlowRegime};
use crate::geometry::{seval, PaneledAirfoil};
use crate::solver::blstate::BlState;

// ============================================================================
// Geometry
// ============================================================================

/// Wake node geometry: XYWAKE's nodes `N+1..=N+NW` and their normals. The wake normal is
/// `-grad(psi)/|grad(psi)|` (`xpanel.f:1311-1330`), which points towards the *lower* side.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WakeGeometryOutput {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Arc length continued from the airfoil's S(N)
    pub s: Vec<f64>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
}

/// Airfoil node geometry as the solver holds it (XFOIL's `X, Y, S, NX, NY` for `I = 1..N`), with
/// the TECALC scalars and, after a viscous solve, the wake.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FoilGeometryOutput {
    /// Node coordinates, TE → upper → LE → lower → TE
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Spline arc length at each node
    pub s: Vec<f64>,
    /// Outward unit normals at the nodes (NCALC: spline normals, corner-averaged at doubled nodes)
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub chord: f64,
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    /// Arc length of the leading edge
    pub sle: f64,
    /// TE gap area projected normal to the TE bisector (TECALC's ANTE)
    pub ante: f64,
    pub sharp_te: bool,
    /// Present once XYWAKE has run (any viscous point)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wake: Option<WakeGeometryOutput>,
}

impl FoilGeometryOutput {
    /// Airfoil nodes `1..=n` and, if `st.lwake`, wake nodes `n+1..=n+nw`.
    pub fn from_state(st: &BlState) -> Self {
        let n = st.n;
        let nw = st.nw;
        let wake = if st.lwake && nw > 0 {
            Some(WakeGeometryOutput {
                x: st.x[n + 1..=n + nw].to_vec(),
                y: st.y[n + 1..=n + nw].to_vec(),
                s: st.s[n + 1..=n + nw].to_vec(),
                nx: st.nx[n + 1..=n + nw].to_vec(),
                ny: st.ny[n + 1..=n + nw].to_vec(),
            })
        } else {
            None
        };
        Self {
            x: st.x[1..=n].to_vec(),
            y: st.y[1..=n].to_vec(),
            s: st.s[1..=n].to_vec(),
            nx: st.nx[1..=n].to_vec(),
            ny: st.ny[1..=n].to_vec(),
            chord: st.chord,
            xle: st.xle,
            yle: st.yle,
            xte: st.xte,
            yte: st.yte,
            sle: st.sle,
            ante: st.ante,
            sharp_te: st.sharp,
            wake,
        }
    }

    /// Geometry only (no wake), for plotting panelings without a solve.
    pub fn from_paneled(airfoil: &PaneledAirfoil) -> Self {
        Self::from_state(&BlState::from_airfoil(airfoil, 0))
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
    Hs,
    /// Edge velocity Ue/Vinf after the Karman–Tsien transformation (unsigned on both sides)
    Ue,
    /// Skin-friction coefficient Cf = τ/(½ρ∞V∞²) (zero in the wake)
    Cf,
    /// Dissipation coefficient CD = Di/(ρ∞V∞³) (DUMP's `CDIS`)
    Cdis,
    /// Shear-stress coefficient √Cτ (turbulent) or amplification N (laminar): CTAU
    Ctau,
    /// Equilibrium shear-stress coefficient √Cτ_eq: CTQ
    Ctq,
    /// Slip-velocity parameter 1.6/(1+Us): USLP
    Uslp,
    /// Viscous pressure coefficient at the node: CPV
    Cp,
    /// Mass defect m = Ue·δ*: MASS
    Mass,
}

impl BlQuantity {
    /// Every quantity, in the order the CLI documents them
    pub const ALL: [BlQuantity; 14] = [
        BlQuantity::Dstar,
        BlQuantity::Theta,
        BlQuantity::Delta,
        BlQuantity::H,
        BlQuantity::Hk,
        BlQuantity::Hs,
        BlQuantity::Ue,
        BlQuantity::Cf,
        BlQuantity::Cdis,
        BlQuantity::Ctau,
        BlQuantity::Ctq,
        BlQuantity::Uslp,
        BlQuantity::Cp,
        BlQuantity::Mass,
    ];

    /// CLI / JSON name
    pub fn name(self) -> &'static str {
        match self {
            BlQuantity::Dstar => "dstar",
            BlQuantity::Theta => "theta",
            BlQuantity::Delta => "delta",
            BlQuantity::H => "h",
            BlQuantity::Hk => "hk",
            BlQuantity::Hs => "hs",
            BlQuantity::Ue => "ue",
            BlQuantity::Cf => "cf",
            BlQuantity::Cdis => "cdis",
            BlQuantity::Ctau => "ctau",
            BlQuantity::Ctq => "ctq",
            BlQuantity::Uslp => "uslp",
            BlQuantity::Cp => "cp",
            BlQuantity::Mass => "mass",
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
            BlQuantity::Hs => "H*",
            BlQuantity::Ue => "Ue/V∞",
            BlQuantity::Cf => "Cf",
            BlQuantity::Cdis => "CD",
            BlQuantity::Ctau => "√Cτ / N",
            BlQuantity::Ctq => "√Cτ_eq",
            BlQuantity::Uslp => "1.6/(1+Us)",
            BlQuantity::Cp => "Cp",
            BlQuantity::Mass => "m",
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
pub struct StoredClosures {
    /// TAU = ½·ρ·Ue²·Cf
    pub tau: Vec<f64>,
    /// DIS = ½·ρ·Ue³·Di·H*
    pub dis: Vec<f64>,
    /// CTQ
    pub ctq: Vec<f64>,
    /// DELT
    pub delt: Vec<f64>,
    /// USLP
    pub uslp: Vec<f64>,
    /// TSTR = H*·θ from the march
    pub tstr: Vec<f64>,
    /// DUMP's `H*` = TSTR/THET (1.0 where THET == 0, as XFOIL)
    pub hs_dump: Vec<f64>,
    /// DUMP's `Cf` = TAU/(½QINF²)
    pub cf_dump: Vec<f64>,
}

/// One side's stations (upper, lower or wake) as a struct of arrays; every vector has one entry
/// per station, in marching order.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct BlSideOutput {
    /// IBL (1-based station index on this side, as XFOIL)
    pub ibl: Vec<usize>,
    /// IPAN (1-based panel/wake node index, as XFOIL)
    pub node: Vec<usize>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// BL arc length from the stagnation point (XSSI)
    pub xssi: Vec<f64>,

    // ---- converged primaries (post-UPDATE) ----
    /// UEDG: incompressible edge velocity
    pub uedg: Vec<f64>,
    /// THET
    pub thet: Vec<f64>,
    /// DSTR (in the wake: total, including the wake gap WGAP)
    pub dstr: Vec<f64>,
    /// CTAU (√Cτ turbulent, amplification N laminar)
    pub ctau: Vec<f64>,
    /// MASS
    pub mass: Vec<f64>,

    // ---- live closure values on the primaries above (BLPRV → BLKIN → BLVAR) ----
    /// Ue/Vinf after the Karman–Tsien transformation (unsigned; DUMP signs it by GAM)
    pub ue: Vec<f64>,
    /// H = δ*/θ (wake: without the wake gap, as BLKIN)
    pub h: Vec<f64>,
    /// Hk from HKIN (unclamped)
    pub hk: Vec<f64>,
    /// H* from BLVAR
    pub hs: Vec<f64>,
    /// Cf = ρ·Ue²·cf/QINF² (BLVAR's cf scaled as SETBL forms TAU, then DUMP's normalisation)
    pub cf: Vec<f64>,
    /// CD = DIS/QINF³ with DIS as SETBL forms it
    pub cdis: Vec<f64>,
    /// δ from BLVAR (DE)
    pub delta: Vec<f64>,
    /// √Cτ_eq from BLVAR (CQ)
    pub ctq: Vec<f64>,
    /// 1.6/(1+Us) from BLVAR
    pub uslp: Vec<f64>,
    /// Rθ from BLKIN
    pub rt: Vec<f64>,
    /// Edge Mach² from BLKIN
    pub msq: Vec<f64>,
    /// CPV at the node
    pub cp: Vec<f64>,

    /// XFOIL's stored (lagged) closure arrays
    pub stored: StoredClosures,
}

impl BlSideOutput {
    /// Number of stations
    pub fn len(&self) -> usize {
        self.ibl.len()
    }

    /// True when the side has no stations
    pub fn is_empty(&self) -> bool {
        self.ibl.is_empty()
    }

    /// The column for a quantity
    pub fn column(&self, q: BlQuantity) -> &[f64] {
        match q {
            BlQuantity::Dstar => &self.dstr,
            BlQuantity::Theta => &self.thet,
            BlQuantity::Delta => &self.delta,
            BlQuantity::H => &self.h,
            BlQuantity::Hk => &self.hk,
            BlQuantity::Hs => &self.hs,
            BlQuantity::Ue => &self.ue,
            BlQuantity::Cf => &self.cf,
            BlQuantity::Cdis => &self.cdis,
            BlQuantity::Ctau => &self.ctau,
            BlQuantity::Ctq => &self.ctq,
            BlQuantity::Uslp => &self.uslp,
            BlQuantity::Cp => &self.cp,
            BlQuantity::Mass => &self.mass,
        }
    }

    /// Every column has exactly one entry per station
    pub fn check_lengths(&self) -> Result<(), String> {
        let n = self.len();
        let cols: [(&str, usize); 30] = [
            ("node", self.node.len()),
            ("x", self.x.len()),
            ("y", self.y.len()),
            ("xssi", self.xssi.len()),
            ("uedg", self.uedg.len()),
            ("thet", self.thet.len()),
            ("dstr", self.dstr.len()),
            ("ctau", self.ctau.len()),
            ("mass", self.mass.len()),
            ("ue", self.ue.len()),
            ("h", self.h.len()),
            ("hk", self.hk.len()),
            ("hs", self.hs.len()),
            ("cf", self.cf.len()),
            ("cdis", self.cdis.len()),
            ("delta", self.delta.len()),
            ("ctq", self.ctq.len()),
            ("uslp", self.uslp.len()),
            ("rt", self.rt.len()),
            ("msq", self.msq.len()),
            ("cp", self.cp.len()),
            ("stored.tau", self.stored.tau.len()),
            ("stored.dis", self.stored.dis.len()),
            ("stored.ctq", self.stored.ctq.len()),
            ("stored.delt", self.stored.delt.len()),
            ("stored.uslp", self.stored.uslp.len()),
            ("stored.tstr", self.stored.tstr.len()),
            ("stored.hs_dump", self.stored.hs_dump.len()),
            ("stored.cf_dump", self.stored.cf_dump.len()),
            ("ibl", n),
        ];
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
    pub ist: usize,
    /// SST: spline arc length of the stagnation point
    pub sst: f64,
    pub x: f64,
    pub y: f64,
}

/// Transition on one side
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TransitionMarker {
    /// ITRAN: first turbulent station
    pub station: usize,
    /// TFORCE: transition was forced (XSTRIP)
    pub forced: bool,
    /// XOCTR, YOCTR: chord-projected position
    pub x_c: f64,
    pub y_c: f64,
    /// Spline arc length (SST ∓ XSSITR) and position on the surface
    pub s: f64,
    pub x: f64,
    pub y: f64,
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
    pub upper: BlSideOutput,
    pub lower: BlSideOutput,
    pub wake: BlSideOutput,
    /// IBLTE(1..2): last airfoil station on each side
    pub iblte: [usize; 2],
    /// ITRAN(1..2)
    pub itran: [usize; 2],
    /// NW
    pub nw: usize,
    /// QINF
    pub qinf: f64,
    /// ANTE (TECALC)
    pub ante: f64,
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
    ue: f64,
    h: f64,
    hk: f64,
    hs: f64,
    cf: f64,
    cdis: f64,
    delta: f64,
    ctq: f64,
    uslp: f64,
    rt: f64,
    msq: f64,
}

/// BLPRV → BLKIN → BLVAR on the converged primaries of station (is, ibl), with the flow type
/// SETBL would use there (laminar ahead of ITRAN, turbulent from it, wake past IBLTE).
fn live_closures(st: &BlState, params: &FlowParameters, is: usize, ibl: usize) -> Live {
    let wake = ibl > st.iblte[is];
    let turb = ibl >= st.itran[is];
    let ctau = st.ctau[is][ibl];
    let (ami, cti) = if turb { (0.0, ctau) } else { (ctau, 0.0) };
    let uei = st.uedg[is][ibl];
    let thi = st.thet[is][ibl];
    if thi == 0.0 || uei == 0.0 {
        // an unsolved station (DUMP prints H = H* = 1 there)
        return Live {
            ue: 0.0,
            h: 1.0,
            hk: 1.0,
            hs: 1.0,
            cf: 0.0,
            cdis: 0.0,
            delta: 0.0,
            ctq: 0.0,
            uslp: 0.0,
            rt: 0.0,
            msq: 0.0,
        };
    }
    // DSI = MDI/UEI, DSWAKI = WGAP(IW), as SETBL sets the "2" station
    let dsi = st.mass[is][ibl] / uei;
    let dswaki = if wake { st.wgap[ibl - st.iblte[is]] } else { 0.0 };

    let mut s = BLStationState::default();
    s.blprv(st.xssi[is][ibl], ami, cti, thi, dsi, dswaki, uei, params);
    s.blkin(params);
    let (h, hk, rt, msq) = (s.h, s.hk, s.rt, s.msq);
    let flow = if wake {
        FlowRegime::Wake
    } else if turb {
        FlowRegime::Turbulent
    } else {
        FlowRegime::Laminar
    };
    s.blvar(flow, params);

    let qinf = st.qinf;
    Live {
        ue: s.u / qinf,
        h,
        hk,
        hs: s.hs,
        // TAU = ½·R2·U2²·CF2 (SETBL), Cf = TAU/(½·QINF²) (DUMP)
        cf: s.r * s.u * s.u * s.cf / (qinf * qinf),
        // DIS = R2·U2³·DI2·HS2·½ (SETBL), CDIS = DIS/QINF³ (DUMP)
        cdis: s.r * s.u * s.u * s.u * s.di * s.hs * 0.5 / (qinf * qinf * qinf),
        delta: s.de,
        ctq: s.cq,
        uslp: 1.60 / (1.0 + s.us),
        rt,
        msq,
    }
}

fn side_output(
    st: &BlState,
    params: &FlowParameters,
    is: usize,
    stations: impl Iterator<Item = usize>,
) -> BlSideOutput {
    let mut out = BlSideOutput::default();
    let qinf = st.qinf;
    for ibl in stations {
        let i = st.ipan[is][ibl];
        let live = live_closures(st, params, is, ibl);
        let thet = st.thet[is][ibl];
        out.ibl.push(ibl);
        out.node.push(i);
        out.x.push(st.x[i]);
        out.y.push(st.y[i]);
        out.xssi.push(st.xssi[is][ibl]);
        out.uedg.push(st.uedg[is][ibl]);
        out.thet.push(thet);
        out.dstr.push(st.dstr[is][ibl]);
        out.ctau.push(st.ctau[is][ibl]);
        out.mass.push(st.mass[is][ibl]);
        out.ue.push(live.ue);
        out.h.push(live.h);
        out.hk.push(live.hk);
        out.hs.push(live.hs);
        out.cf.push(live.cf);
        out.cdis.push(live.cdis);
        out.delta.push(live.delta);
        out.ctq.push(live.ctq);
        out.uslp.push(live.uslp);
        out.rt.push(live.rt);
        out.msq.push(live.msq);
        out.cp.push(st.cpv[i]);
        let stored = &mut out.stored;
        stored.tau.push(st.tau[is][ibl]);
        stored.dis.push(st.dis[is][ibl]);
        stored.ctq.push(st.ctq[is][ibl]);
        stored.delt.push(st.delt[is][ibl]);
        stored.uslp.push(st.uslp[is][ibl]);
        stored.tstr.push(st.tstr[is][ibl]);
        // BLDUMP: IF(TH.EQ.0.0) HS = 1.0 ELSE HS = TS/TH
        stored
            .hs_dump
            .push(if thet == 0.0 { 1.0 } else { st.tstr[is][ibl] / thet });
        stored.cf_dump.push(st.tau[is][ibl] / (0.5 * qinf * qinf));
    }
    out
}

/// Position on the airfoil spline at arc length `s`
fn spline_point(st: &BlState, s: f64) -> (f64, f64) {
    let n = st.n;
    let x = seval(s, &st.x[1..=n], &st.xp[1..=n], &st.s[1..=n]);
    let y = seval(s, &st.y[1..=n], &st.yp[1..=n], &st.s[1..=n]);
    (x, y)
}

/// Sign changes of the live `cf` along one side's surface stations, from the pair (2, 3) on
fn separation_markers(st: &BlState, side: &BlSideOutput, is: usize) -> Vec<SeparationMarker> {
    let mut out = Vec::new();
    for k in 1..side.len() {
        let (prev, cur) = (side.cf[k - 1], side.cf[k]);
        let kind = if prev >= 0.0 && cur < 0.0 {
            SeparationKind::Separation
        } else if prev < 0.0 && cur >= 0.0 {
            SeparationKind::Reattachment
        } else {
            continue;
        };
        let (s0, s1) = (st.s[side.node[k - 1]], st.s[side.node[k]]);
        let frac = if cur == prev { 0.0 } else { prev / (prev - cur) };
        let s = s0 + frac * (s1 - s0);
        let (x, y) = spline_point(st, s);
        out.push(SeparationMarker {
            side: is,
            kind,
            station_before: side.ibl[k - 1],
            s,
            x,
            y,
        });
    }
    out
}

impl BoundaryLayerOutput {
    /// Every station of a solved state. Requires the pointer layer and BL arrays (any state after
    /// a viscous VISCAL call, converged or not).
    pub fn from_state(st: &BlState) -> Self {
        let params = FlowParameters::new(st.minf, st.reinf, st.gamma);
        let upper = side_output(st, &params, 1, 2..=st.iblte[1]);
        let lower = side_output(st, &params, 2, 2..=st.iblte[2]);
        let wake = side_output(st, &params, 2, st.iblte[2] + 1..=st.nbl[2]);

        let (sx, sy) = spline_point(st, st.sst);
        let stagnation = StagnationMarker {
            ist: st.ist,
            sst: st.sst,
            x: sx,
            y: sy,
        };
        let transition_marker = |is: usize| {
            let s = if is == 1 {
                st.sst - st.xssitr[is]
            } else {
                st.sst + st.xssitr[is]
            };
            let (x, y) = spline_point(st, s);
            TransitionMarker {
                station: st.itran[is],
                forced: st.tforce[is],
                x_c: st.xoctr[is],
                y_c: st.yoctr[is],
                s,
                x,
                y,
            }
        };
        let transition = [transition_marker(1), transition_marker(2)];

        let mut derived_separation = separation_markers(st, &upper, 1);
        derived_separation.extend(separation_markers(st, &lower, 2));

        // CPDISP: upper and lower wake Dstar fractions from the first wake point
        let dstrte = st.dstr[2][st.iblte[2] + 1];
        let wake_split = if dstrte != 0.0 {
            [
                (st.dstr[1][st.iblte[1]] + 0.5 * st.ante) / dstrte,
                (st.dstr[2][st.iblte[2]] + 0.5 * st.ante) / dstrte,
            ]
        } else {
            [0.5, 0.5]
        };

        Self {
            upper,
            lower,
            wake,
            iblte: [st.iblte[1], st.iblte[2]],
            itran: [st.itran[1], st.itran[2]],
            nw: st.nw,
            qinf: st.qinf,
            ante: st.ante,
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
        use crate::geometry::{create_paneled_airfoil, naca_4digit};
        let geom = naca_4digit("2412", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let g = FoilGeometryOutput::from_paneled(&airfoil);
        assert_eq!(g.n(), 100);
        assert!(g.wake.is_none());
        assert_eq!(g.x, airfoil.x);
        assert_eq!(g.nx, airfoil.nx);
        assert_eq!(g.sharp_te, airfoil.sharp_te);
        assert!(g.same_panels(&g.clone()));
        let mut h = g.clone();
        h.y[10] = f64::from_bits(h.y[10].to_bits() + 1);
        assert!(!g.same_panels(&h));
    }
}
