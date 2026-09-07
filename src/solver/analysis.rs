//! The OPER analysis driver on `SolverState`: one persistent session per airfoil and flow
//! specification, `alfa` = XFOIL's `ALFA` command (SPECAL, then VISCAL when viscous),
//! `init` = XFOIL's `INIT`, and the polar sweep as CLAUDE.md prescribes it (0° → max,
//! reinitialise, −step → min, stitched ascending).

use crate::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
use serde::{Deserialize, Serialize};

use crate::geometry::PanelledFoil;
use crate::solver::blstate::SolverState;
use crate::solver::ggcalc::InviscidSystem;
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
        let mut st = SolverState::from_foil(airfoil, nw);
        st.re_cl1 = spec.re.unwrap_or(0.0);
        st.re = spec.re.unwrap_or(0.0);
        st.mach_cl1 = spec.mach;
        st.mach = spec.mach;
        st.mach_cl_dependence = spec.mach_cl_dependence;
        st.re_cl_dependence = spec.re_cl_dependence;
        st.amplification_model = spec.amplification_model;
        st.ncrit = [0.0, spec.ncrit, spec.ncrit];
        st.elimination_threshold = spec.elimination_threshold;
        st.x_trip = [0.0, spec.x_trip[0], spec.x_trip[1]];
        st.viscous = spec.re.is_some();
        st.alpha_specified = true;
        st.qinf = 1.0;
        Self {
            state: st,
            inviscid: None,
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
        self.solve_point(self.conditions.max_iterations)
    }

    /// OPER `CL`: SPECCL for the specified CL (alpha is the unknown), then VISCAL(ITMAX) when
    /// viscous — UPDATE then drives alpha so that the viscous CL meets CLSPEC.
    pub fn cl(&mut self, clspec: f64) -> PointResult {
        cl_command(&mut self.state, &mut self.inviscid, clspec);
        self.solve_point(self.conditions.max_iterations)
    }

    /// One point of OPER `ASEQ`: SPECAL for the new angle, then VISCAL(ITMAX + 5) when viscous.
    #[doc(alias = "ASEQ")]
    pub fn sequence_point(&mut self, alpha: f64) -> PointResult {
        sequence_command(&mut self.state, &mut self.inviscid, alpha);
        self.solve_point(self.conditions.max_iterations + 5)
    }

    fn solve_point(&mut self, niter: usize) -> PointResult {
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
        let st = &self.state;
        PointResult {
            alpha: st.alpha,
            cl: st.cl,
            cd: st.cd,
            cd_friction: st.cd_friction,
            cd_pressure: st.cd_pressure,
            cm: st.cm,
            cl_d_alpha: st.cl_d_alpha,
            transition_upper: [st.x_transition[1], st.y_transition[1]],
            transition_lower: [st.x_transition[2], st.y_transition[2]],
            i_transition_station: st.i_transition_station,
            converged,
            iterations: trace.len(),
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

/// Result of a polar sweep: points sorted by ascending alpha, in the order XFOIL's `PACC`
/// would keep them (unconverged viscous points are recorded in `failed_alphas`, not in `points`).
#[derive(Debug, Clone)]
pub struct PolarResult {
    pub results: Vec<PointResult>,
    /// Alphas (radians) that did not converge
    pub failed_alphas: Vec<f64>,
    pub conditions: FlowConditions,
    /// Neither sequence was halted by NSEQEX consecutive failures
    pub completed: bool,
}

impl PolarResult {
    /// (CL_max, alpha in degrees at CL_max)
    pub fn cl_max(&self) -> Option<(f64, f64)> {
        self.results
            .iter()
            .map(|p| (p.cl, p.alpha.to_degrees()))
            .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
    }
    /// (max L/D, CL at max L/D)
    pub fn ldratio_max(&self) -> Option<(f64, f64)> {
        self.results
            .iter()
            .filter(|p| p.cd > 1e-10)
            .map(|p| (p.cl / p.cd, p.cl))
            .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
    }
    /// CD at the point with the smallest |CL|
    pub fn cd0(&self) -> Option<f64> {
        self.results
            .iter()
            .min_by(|a, b| a.cl.abs().partial_cmp(&b.cl.abs()).unwrap())
            .map(|p| p.cd)
    }
}

/// The polar as XFOIL's OPER script runs it (CLAUDE.md):
/// `ALFA 0` / `ASEQ step alpha_max step` / `INIT` / `ALFA -step` / `ASEQ -2step alpha_min -step`,
/// with one persistent session so each point starts from the previous point's BL, and each
/// ASEQ halting after NSEQEX consecutive non-converged points. Points are stitched ascending.
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
    let mut failed = Vec::new();
    let mut completed = true;

    // ASEQ: NPOINT = INT((A2-A1)/DA + 0.5) + 1 points from A1 in steps of DA
    let aseq_alphas = |a1: f64, a2: f64, da: f64| -> Vec<f64> {
        if da == 0.0 || (a2 - a1) * da < 0.0 {
            return Vec::new();
        }
        let npoint = ((a2 - a1) / da + 0.5).floor() as usize + 1;
        (0..npoint).map(|i| a1 + da * i as f64).collect()
    };
    let record = |p: PointResult, points: &mut Vec<PointResult>, failed: &mut Vec<f64>| {
        if p.converged {
            points.push(p);
        } else {
            failed.push(p.alpha);
        }
    };
    let aseq = |session: &mut Session,
                alphas: Vec<f64>,
                points: &mut Vec<PointResult>,
                failed: &mut Vec<f64>,
                observe: &mut dyn FnMut(&Session, &PointResult)|
     -> bool {
        let mut iseqex = 0;
        for adeg in alphas {
            let p = session.sequence_point(adeg.to_radians());
            observe(session, &p);
            let conv = p.converged;
            record(p, points, failed);
            if session.state.viscous && !conv {
                iseqex += 1;
                if iseqex >= config.max_consecutive_failures {
                    // 'Sequence halted since previous N points did not converge'
                    return false;
                }
            } else {
                iseqex = 0;
            }
        }
        true
    };

    // ALFA 0, ASEQ step alpha_max step
    let p = session.alpha(0.0);
    observe(&session, &p);
    record(p, &mut points, &mut failed);
    if config.alpha_max >= step {
        completed &= aseq(
            &mut session,
            aseq_alphas(step, config.alpha_max, step),
            &mut points,
            &mut failed,
            observe,
        );
    }

    // INIT, ALFA -step, ASEQ -2step alpha_min -step
    if config.alpha_min <= -step {
        session.init();
        let p = session.alpha(-step.to_radians());
        observe(&session, &p);
        record(p, &mut points, &mut failed);
        if config.alpha_min <= -2.0 * step {
            completed &= aseq(
                &mut session,
                aseq_alphas(-2.0 * step, config.alpha_min, -step),
                &mut points,
                &mut failed,
                observe,
            );
        }
    }

    points.sort_by(|p, q| p.alpha.partial_cmp(&q.alpha).unwrap());
    PolarResult {
        results: points,
        failed_alphas: failed,
        conditions: config.conditions.clone(),
        completed,
    }
}
