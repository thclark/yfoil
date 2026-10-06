//! The OPER analysis driver on `SolverState`: one persistent session per airfoil and flow
//! specification, `alfa` = XFOIL's `ALFA` command (SPECAL, then VISCAL when viscous),
//! `init` = XFOIL's `INIT`, and the polar sweep as CLAUDE.md prescribes it (0° → max,
//! reinitialise and re-solve 0°, 0° → min, stitched ascending).

use crate::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
use serde::{Deserialize, Serialize};

use crate::geometry::PanelledFoil;
use crate::solver::blstate::SolverState;
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::point_id::IdHasher;
use crate::solver::specal::{alpha_command, cl_command, sequence_command};
use crate::solver::viscal::{solve_viscous, IterationRecord};

/// The flow conditions shared by every point of a polar: the OPER settings that must be pinned
/// explicitly. Also the `conditions` block of the analysis and polar JSON outputs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowConditions {
    /// Reynolds number REINF1; `None` for an inviscid analysis
    pub re: Option<f64>,
    /// Mach number MINF1
    pub mach: f64,
    /// Critical amplification ACRIT (both sides)
    pub ncrit: f64,
    /// ITMAX: VISCAL iteration limit
    pub max_iterations: usize,
    /// WAKLEN: wake length in chords
    pub wake_length: f64,
    /// VACCEL: BLSOLV sparse-elimination threshold
    pub elimination_threshold: f64,
    /// XSTRIP(1..2): forced-transition x/c per side (1.0 = free transition)
    pub x_trip: [f64; 2],
    /// MATYP / RETYP: Mach and Re dependence on CL (1 = fixed)
    pub mach_cl_dependence: MachClDependence,
    pub re_cl_dependence: ReClDependence,
    /// OPER `DAMP`: modified envelope e^n amplification (IDAMP = 1, DAMPL2)
    pub amplification_model: AmplificationModel,
}

impl Default for FlowConditions {
    fn default() -> Self {
        Self {
            re: Some(1.0e6),
            mach: 0.0,
            ncrit: 9.0,
            max_iterations: 20,
            wake_length: 1.0,
            elimination_threshold: 0.01,
            x_trip: [1.0, 1.0],
            mach_cl_dependence: MachClDependence::Fixed,
            re_cl_dependence: ReClDependence::Fixed,
            amplification_model: AmplificationModel::Envelope,
        }
    }
}

/// One converged (or not) operating point.
#[derive(Debug, Clone, Default)]
pub struct PointResult {
    /// Angle of attack (radians)
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
    pub cd_friction: f64,
    pub cd_pressure: f64,
    pub cm: f64,
    pub cl_d_alpha: f64,
    /// (x, y) of the transition point on the upper and lower side, chord fractions
    /// (XOCTR(IS), YOCTR(IS))
    pub transition_upper: [f64; 2],
    pub transition_lower: [f64; 2],
    pub i_transition_station: [usize; 3],
    /// LVCONV after VISCAL (true for an inviscid point)
    pub converged: bool,
    /// VISCAL iterations performed
    pub iterations: usize,
    /// SPECCL's exit iteration for an OPER `CL` point (XFOIL's ITAL: 21 when the 20-iteration
    /// alpha Newton was exhausted), 0 for an `ALFA` point
    pub inviscid_cl_iterations: usize,
    /// The evidence for whether this point's values are valid: Kármán–Tsien margins, the CL(M)
    /// Newton's exhaustion and MRCL's substitutions (`crate::solver::validity`). Recorded by the
    /// solver, never read by it; the output layer classifies the point from this.
    pub validity: crate::solver::validity::ValidityRecord,
    /// This point's stable reference: a random 8-character id, unique to this solve. Random
    /// rather than positional so that points from different runs can be held together — a polar
    /// restarted from a point of an earlier one, or two sweeps merged — without the references
    /// colliding, which a per-run counter would guarantee.
    pub id: String,
    /// The `id` of the point whose converged boundary layer seeded this solve, or `None` when the
    /// BL was marched fresh (MRCHUE from the inviscid solution). A polar is a state machine, so
    /// following this back gives the chain of states a point depends on — two chains rooted at
    /// the 0° solve, one per leg — which recovers the order of execution and lets a sweep be
    /// restarted part way through rather than from the beginning.
    pub initialised_from: Option<String>,
    /// RMSBL of the last iteration (0 for an inviscid point)
    pub residual: f64,
    /// Per-iteration record
    pub iteration_records: Vec<IterationRecord>,
}

/// A persistent analysis session: the geometry, inviscid system and BL state that XFOIL keeps
/// in COMMON between OPER commands.
#[derive(Debug, Clone)]
pub struct Session {
    state: SolverState,
    inviscid: Option<InviscidSystem>,
    conditions: FlowConditions,
    /// Everything a point id depends on that does not change within a session: the yFoil version,
    /// the panel geometry and the flow conditions (`crate::solver::point_id`). Hashed once here,
    /// then each point mixes in its own operating point and predecessor.
    inputs_hash: IdHasher,
    /// The id of the last point solved, which is the predecessor of the next one whenever the
    /// boundary layer is carried across.
    last_point_id: Option<String>,
}

impl Session {
    /// The solver state (XFOIL's COMMON blocks) as the last command left it
    pub fn state(&self) -> &SolverState {
        &self.state
    }

    /// The flow conditions the session was opened with
    pub fn conditions(&self) -> &FlowConditions {
        &self.conditions
    }

    /// The inviscid system (AIJ factors, BIJ), once the first SPECAL/SPECCL has built it
    pub fn inviscid(&self) -> Option<&InviscidSystem> {
        self.inviscid.as_ref()
    }

    /// Mutable state, for driving the translated subroutines directly (fixture tests). The
    /// OPER commands on `Session` keep XFOIL's flag set consistent; a caller of this does not.
    pub fn state_mut(&mut self) -> &mut SolverState {
        &mut self.state
    }

    /// Mutable inviscid-system slot, for the same purpose as [`Session::state_mut`]
    pub fn inviscid_mut(&mut self) -> &mut Option<InviscidSystem> {
        &mut self.inviscid
    }

    /// Both mutable parts at once, as the translated SPECAL/VISCAL signatures take them
    pub fn parts_mut(&mut self) -> (&mut SolverState, &mut Option<InviscidSystem>) {
        (&mut self.state, &mut self.inviscid)
    }

    /// LOAD + OPER settings: geometry in, VISC/MACH/N/ITER/VACCEL/XTR pinned.
    pub fn new(airfoil: &PanelledFoil, spec: FlowConditions) -> Self {
        // NW = N/12 + 10*INT(WAKLEN)
        let nw = airfoil.n_foil_nodes / 12 + 10 * (spec.wake_length as usize);
        let mut state = SolverState::from_foil(airfoil, nw);
        state.re_cl1 = spec.re.unwrap_or(0.0);
        state.re = spec.re.unwrap_or(0.0);
        state.mach_cl1 = spec.mach;
        state.mach = spec.mach;
        state.mach_cl_dependence = spec.mach_cl_dependence;
        state.re_cl_dependence = spec.re_cl_dependence;
        state.amplification_model = spec.amplification_model;
        state.ncrit = [0.0, spec.ncrit, spec.ncrit];
        state.elimination_threshold = spec.elimination_threshold;
        state.x_trip = [0.0, spec.x_trip[0], spec.x_trip[1]];
        state.viscous = spec.re.is_some();
        state.alpha_specified = true;
        state.qinf = 1.0;
        let mut inputs_hash = IdHasher::new();
        inputs_hash
            .str("yfoil_version", env!("CARGO_PKG_VERSION"))
            .f64_slice("x", &airfoil.x)
            .f64_slice("y", &airfoil.y)
            .opt_f64("re", spec.re)
            .f64("mach", spec.mach)
            .f64("ncrit", spec.ncrit)
            .usize("max_iterations", spec.max_iterations)
            .f64("wake_length", spec.wake_length)
            .f64("elimination_threshold", spec.elimination_threshold)
            .f64("x_trip_1", spec.x_trip[0])
            .f64("x_trip_2", spec.x_trip[1])
            .u8("mach_cl_dependence", spec.mach_cl_dependence as u8)
            .u8("re_cl_dependence", spec.re_cl_dependence as u8)
            .u8("amplification_model", spec.amplification_model as u8);
        Self {
            state,
            inviscid: None,
            inputs_hash,
            last_point_id: None,
            conditions: spec,
        }
    }

    /// OPER `INIT`: BL initialisation flag toggled off so the next VISCAL re-marches with
    /// MRCHUE, and the pointer layer is rebuilt.
    #[doc(alias = "INIT")]
    pub fn init(&mut self) {
        self.state.bl_initialised = !self.state.bl_initialised;
        if !self.state.bl_initialised {
            // 'BLs will be initialized on next point'
            self.state.pointers_built = false;
        }
    }

    /// OPER `ALFA`: SPECAL for the new angle, then VISCAL(ITMAX) when viscous.
    #[doc(alias = "ALFA")]
    pub fn alpha(&mut self, alpha: f64) -> PointResult {
        alpha_command(&mut self.state, &mut self.inviscid, alpha);
        self.solve_point(self.conditions.max_iterations, 0)
    }

    /// OPER `CL`: SPECCL for the specified CL (alpha is the unknown), then VISCAL(ITMAX) when
    /// viscous — UPDATE then drives alpha so that the viscous CL meets CLSPEC.
    pub fn cl(&mut self, clspec: f64) -> PointResult {
        let ital = cl_command(&mut self.state, &mut self.inviscid, clspec);
        self.solve_point(self.conditions.max_iterations, ital)
    }

    /// One point of OPER `ASEQ`: SPECAL for the new angle, then VISCAL(ITMAX + 5) when viscous.
    #[doc(alias = "ASEQ")]
    pub fn sequence_point(&mut self, alpha: f64) -> PointResult {
        sequence_command(&mut self.state, &mut self.inviscid, alpha);
        self.solve_point(self.conditions.max_iterations + 5, 0)
    }

    fn solve_point(&mut self, niter: usize, inviscid_cl_iterations: usize) -> PointResult {
        // The predecessor is the point whose converged BL seeds this solve. There is one only if
        // the BL is being carried across: an inviscid point depends on nothing before it (SPECAL
        // resets CLM and GAMU is alpha-independent), and after INIT the march starts fresh.
        let initialised_from = if self.state.viscous && self.state.bl_initialised {
            self.last_point_id.clone()
        } else {
            None
        };
        let mut id_hash = self.inputs_hash.clone();
        id_hash
            .u8("alpha_specified", u8::from(self.state.alpha_specified))
            .f64(
                "operating_point",
                if self.state.alpha_specified {
                    self.state.alpha
                } else {
                    self.state.cl_specified
                },
            )
            .usize("niter", niter)
            .opt_str("initialised_from", initialised_from.as_deref());
        let id = id_hash.finish();
        self.last_point_id = Some(id.clone());

        let mut trace = Vec::new();
        let converged = if self.state.viscous {
            solve_viscous(
                &mut self.state,
                self.inviscid.as_mut(),
                niter,
                self.conditions.wake_length,
                Some(&mut trace),
            )
        } else {
            true
        };
        let state = &self.state;
        PointResult {
            alpha: state.alpha,
            cl: state.cl,
            cd: state.cd,
            cd_friction: state.cd_friction,
            cd_pressure: state.cd_pressure,
            cm: state.cm,
            cl_d_alpha: state.cl_d_alpha,
            transition_upper: [state.x_transition[1], state.y_transition[1]],
            transition_lower: [state.x_transition[2], state.y_transition[2]],
            i_transition_station: state.i_transition_station,
            converged,
            iterations: trace.len(),
            inviscid_cl_iterations,
            validity: state.validity,
            id,
            initialised_from,
            residual: trace.last().map(|t| t.residual).unwrap_or(0.0),
            iteration_records: trace,
        }
    }
}

/// Translates XFOIL's `ALFA`.
///
/// Single operating point from scratch (fresh session).
#[doc(alias = "ALFA")]
pub fn analyse(airfoil: &PanelledFoil, alpha: f64, spec: &FlowConditions) -> PointResult {
    Session::new(airfoil, spec.clone()).alpha(alpha)
}

/// Polar sweep configuration.
#[derive(Debug, Clone)]
pub struct PolarConfig {
    /// Maximum angle of attack (degrees)
    pub alpha_max: f64,
    /// Minimum angle of attack (degrees)
    pub alpha_min: f64,
    /// Step size (degrees, positive)
    pub alpha_step: f64,
    pub conditions: FlowConditions,
    /// NSEQEX: an ASEQ sequence halts once this many consecutive points fail to converge
    pub max_consecutive_failures: usize,
}

impl Default for PolarConfig {
    fn default() -> Self {
        Self {
            alpha_max: 15.0,
            alpha_min: -5.0,
            alpha_step: 0.5,
            conditions: FlowConditions::default(),
            max_consecutive_failures: 4,
        }
    }
}

/// Result of a polar sweep: **every** alpha the sweep was asked for, sorted ascending. A point
/// that did not converge is kept here with its whole record — the caller decides what to do with
/// it — and an alpha the sweep halted before reaching is listed in `not_attempted`, so nothing
/// the caller asked for goes missing from the result.
#[derive(Debug, Clone)]
pub struct PolarResult {
    /// Every point the sweep solved, converged or not, ascending in alpha
    pub results: Vec<PointResult>,
    /// Alphas (radians) the sweep halted before reaching, ascending. The halt is XFOIL's own
    /// NSEQEX rule and is not changed by recording what it skipped.
    pub not_attempted: Vec<f64>,
    pub conditions: FlowConditions,
    /// Neither sequence was halted by NSEQEX consecutive failures
    pub completed: bool,
}

/// The polar as XFOIL's OPER script runs it (CLAUDE.md):
/// `ALFA 0` / `ASEQ step alpha_max step` / `INIT` / `ALFA 0` / `ASEQ -step alpha_min -step`,
/// with one persistent session so each point starts from the previous point's BL, and each
/// ASEQ halting after NSEQEX consecutive non-converged points. Both legs start from the 0°
/// solution: OPER cannot store a BL state, so the second leg re-solves 0° from a fresh march
/// (`INIT` / `ALFA 0`) and sweeps down from it. That re-solve seeds the leg and is not a polar
/// point — the polar carries the first 0° solve — and the polar fixture test asserts the two
/// solves are identical. Points are stitched ascending.
pub fn compute_polar(airfoil: &PanelledFoil, config: &PolarConfig) -> PolarResult {
    compute_polar_with(airfoil, config, &mut |_, _| {})
}

/// [`compute_polar`] with an observer called after every point the sweep visits, converged or
/// not, with the session in that point's state (before the next SPECAL moves it on). This is how
/// `yfoil polar --distributions` captures the BL of each point: a point reached inside a sweep
/// starts from the previous alpha's BL and is not the same solve as a fresh `analyze`.
pub fn compute_polar_with(
    airfoil: &PanelledFoil,
    config: &PolarConfig,
    observe: &mut dyn FnMut(&Session, &PointResult),
) -> PolarResult {
    let step = config.alpha_step.abs();
    let mut session = Session::new(airfoil, config.conditions.clone());
    let mut points = Vec::new();
    let mut not_attempted = Vec::new();
    let mut completed = true;

    // ASEQ: NPOINT = INT((A2-A1)/DA + 0.5) + 1 points from A1 in steps of DA
    let aseq_alphas = |a1: f64, a2: f64, da: f64| -> Vec<f64> {
        if da == 0.0 || (a2 - a1) * da < 0.0 {
            return Vec::new();
        }
        let npoint = ((a2 - a1) / da + 0.5).floor() as usize + 1;
        (0..npoint).map(|i| a1 + da * i as f64).collect()
    };
    let aseq = |session: &mut Session,
                alphas: Vec<f64>,
                points: &mut Vec<PointResult>,
                not_attempted: &mut Vec<f64>,
                observe: &mut dyn FnMut(&Session, &PointResult)|
     -> bool {
        let mut iseqex = 0;
        for (k, adeg) in alphas.iter().enumerate() {
            let p = session.sequence_point(adeg.to_radians());
            observe(session, &p);
            let conv = p.converged;
            points.push(p);
            if session.state.viscous && !conv {
                iseqex += 1;
                if iseqex >= config.max_consecutive_failures {
                    // 'Sequence halted since previous N points did not converge'. The halt is
                    // XFOIL's; recording the alphas it skipped is not a change to it.
                    not_attempted.extend(alphas[k + 1..].iter().map(|a| a.to_radians()));
                    return false;
                }
            } else {
                iseqex = 0;
            }
        }
        true
    };

    // ALFA 0, ASEQ step alpha_max step
    // ALFA 0: a fresh march, so it is seeded by nothing and roots both legs
    let p = session.alpha(0.0);
    observe(&session, &p);
    points.push(p);
    if config.alpha_max >= step {
        completed &= aseq(
            &mut session,
            aseq_alphas(step, config.alpha_max, step),
            &mut points,
            &mut not_attempted,
            observe,
        );
    }

    // INIT, ALFA 0 (re-solve to seed the leg; not recorded), ASEQ -step alpha_min -step
    if config.alpha_min <= -step {
        session.init();
        session.alpha(0.0);
        // The re-solve seeds the downward leg and is not a polar point (OPER cannot store a BL
        // state, so the leg starts from a fresh MRCHUE march at 0°). Its inputs — geometry,
        // conditions, operating point, iteration limit, and no predecessor — are those of the
        // first 0° solve, so it takes the same content-addressed id, and the leg's first point
        // cites the recorded 0° point rather than a phantom. That the two solves really are
        // identical is what the polar fixture test asserts.
        completed &= aseq(
            &mut session,
            aseq_alphas(-step, config.alpha_min, -step),
            &mut points,
            &mut not_attempted,
            observe,
        );
    }

    points.sort_by(|p, q| p.alpha.partial_cmp(&q.alpha).unwrap());
    not_attempted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    PolarResult {
        results: points,
        not_attempted,
        conditions: config.conditions.clone(),
        completed,
    }
}
