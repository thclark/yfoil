//! The OPER analysis driver on `BlState`: one persistent session per airfoil and flow
//! specification, `alfa` = XFOIL's `ALFA` command (SPECAL, then VISCAL when viscous),
//! `init` = XFOIL's `INIT`, and the polar sweep as CLAUDE.md prescribes it (0° → max,
//! reinitialise, −step → min, stitched ascending).

use crate::geometry::PaneledAirfoil;
use crate::solver::blstate::BlState;
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::specal::{alfa_command, aseq_point, cl_command};
use crate::solver::viscal::{solve_viscous, IterationRecord};

/// Flow specification (the OPER settings that must be pinned explicitly).
#[derive(Debug, Clone)]
pub struct FlowSpec {
    /// Reynolds number REINF1 (0 → inviscid analysis)
    pub re: f64,
    /// Mach number MINF1
    pub mach: f64,
    /// Critical amplification ACRIT (both sides)
    pub ncrit: f64,
    /// ITMAX: VISCAL iteration limit
    pub itmax: usize,
    /// WAKLEN: wake length in chords
    pub waklen: f64,
    /// VACCEL: BLSOLV sparse-elimination threshold
    pub vaccel: f64,
    /// XSTRIP(1..2): forced-transition x/c per side (1.0 = free transition)
    pub xstrip: [f64; 2],
    /// MATYP / RETYP: Mach and Re dependence on CL (1 = fixed)
    pub matyp: usize,
    pub retyp: usize,
    /// OPER `DAMP`: modified envelope e^n amplification (IDAMP = 1, DAMPL2)
    pub idamp: bool,
}

impl Default for FlowSpec {
    fn default() -> Self {
        Self {
            re: 1.0e6,
            mach: 0.0,
            ncrit: 9.0,
            itmax: 20,
            waklen: 1.0,
            vaccel: 0.01,
            xstrip: [1.0, 1.0],
            matyp: 1,
            retyp: 1,
            idamp: false,
        }
    }
}

/// One converged (or not) operating point.
#[derive(Debug, Clone, Default)]
pub struct OperatingPoint {
    /// Angle of attack (radians)
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cdp: f64,
    pub cm: f64,
    pub cl_alf: f64,
    /// XOCTR(1), XOCTR(2): transition x/c on the upper and lower side
    pub xtr_upper: f64,
    pub xtr_lower: f64,
    pub itran: [usize; 3],
    /// LVCONV after VISCAL (true for an inviscid point)
    pub converged: bool,
    /// VISCAL iterations performed
    pub iterations: usize,
    /// RMSBL of the last iteration (0 for an inviscid point)
    pub rmsbl: f64,
    /// Per-iteration record
    pub trace: Vec<IterationRecord>,
}

/// A persistent analysis session: the geometry, inviscid system and BL state that XFOIL keeps
/// in COMMON between OPER commands.
#[derive(Debug, Clone)]
pub struct Session {
    pub st: BlState,
    pub sys: Option<InviscidSystem>,
    pub spec: FlowSpec,
}

impl Session {
    /// LOAD + OPER settings: geometry in, VISC/MACH/N/ITER/VACCEL/XTR pinned.
    pub fn new(airfoil: &PaneledAirfoil, spec: FlowSpec) -> Self {
        // NW = N/12 + 10*INT(WAKLEN)
        let nw = airfoil.n / 12 + 10 * (spec.waklen as usize);
        let mut st = BlState::from_airfoil(airfoil, nw);
        st.reinf1 = spec.re;
        st.reinf = spec.re;
        st.minf1 = spec.mach;
        st.minf = spec.mach;
        st.matyp = spec.matyp;
        st.retyp = spec.retyp;
        st.idamp = usize::from(spec.idamp);
        st.acrit = [0.0, spec.ncrit, spec.ncrit];
        st.vaccel = spec.vaccel;
        st.xstrip = [0.0, spec.xstrip[0], spec.xstrip[1]];
        st.lvisc = spec.re > 0.0;
        st.lalfa = true;
        st.qinf = 1.0;
        Self { st, sys: None, spec }
    }

    /// OPER `INIT`: BL initialisation flag toggled off so the next VISCAL re-marches with
    /// MRCHUE, and the pointer layer is rebuilt.
    pub fn init(&mut self) {
        self.st.lblini = !self.st.lblini;
        if !self.st.lblini {
            // 'BLs will be initialized on next point'
            self.st.lipan = false;
        }
    }

    /// OPER `ALFA`: SPECAL for the new angle, then VISCAL(ITMAX) when viscous.
    pub fn alfa(&mut self, alpha: f64) -> OperatingPoint {
        alfa_command(&mut self.st, &mut self.sys, alpha);
        self.run_viscal(self.spec.itmax)
    }

    /// OPER `CL`: SPECCL for the specified CL (alpha is the unknown), then VISCAL(ITMAX) when
    /// viscous — UPDATE then drives alpha so that the viscous CL meets CLSPEC.
    pub fn cl(&mut self, clspec: f64) -> OperatingPoint {
        cl_command(&mut self.st, &mut self.sys, clspec);
        self.run_viscal(self.spec.itmax)
    }

    /// One point of OPER `ASEQ`: SPECAL for the new angle, then VISCAL(ITMAX + 5) when viscous.
    pub fn aseq(&mut self, alpha: f64) -> OperatingPoint {
        aseq_point(&mut self.st, &mut self.sys, alpha);
        self.run_viscal(self.spec.itmax + 5)
    }

    fn run_viscal(&mut self, niter: usize) -> OperatingPoint {
        let mut trace = Vec::new();
        let converged = if self.st.lvisc {
            solve_viscous(
                &mut self.st,
                self.sys.as_mut(),
                niter,
                self.spec.waklen,
                Some(&mut trace),
            )
        } else {
            true
        };
        let st = &self.st;
        OperatingPoint {
            alpha: st.alfa,
            cl: st.cl,
            cd: st.cd,
            cdf: st.cdf,
            cdp: st.cdp,
            cm: st.cm,
            cl_alf: st.cl_alf,
            xtr_upper: st.xoctr[1],
            xtr_lower: st.xoctr[2],
            itran: st.itran,
            converged,
            iterations: trace.len(),
            rmsbl: trace.last().map(|t| t.residual).unwrap_or(0.0),
            trace,
        }
    }
}

/// Single operating point from scratch (fresh session).
pub fn analyze(airfoil: &PaneledAirfoil, alpha: f64, spec: &FlowSpec) -> OperatingPoint {
    Session::new(airfoil, spec.clone()).alfa(alpha)
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
    pub spec: FlowSpec,
    /// NSEQEX: an ASEQ sequence halts once this many consecutive points fail to converge
    pub nseqex: usize,
}

impl Default for PolarConfig {
    fn default() -> Self {
        Self {
            alpha_max: 15.0,
            alpha_min: -5.0,
            alpha_step: 0.5,
            spec: FlowSpec::default(),
            nseqex: 4,
        }
    }
}

/// Result of a polar sweep: points sorted by ascending alpha, in the order XFOIL's `PACC`
/// would keep them (unconverged viscous points are recorded in `failed_alphas`, not in `points`).
#[derive(Debug, Clone)]
pub struct PolarResult {
    pub points: Vec<OperatingPoint>,
    /// Alphas (radians) that did not converge
    pub failed_alphas: Vec<f64>,
    pub spec: FlowSpec,
    /// Neither sequence was halted by NSEQEX consecutive failures
    pub completed: bool,
}

impl PolarResult {
    /// (CL_max, alpha in degrees at CL_max)
    pub fn cl_max(&self) -> Option<(f64, f64)> {
        self.points
            .iter()
            .map(|p| (p.cl, p.alpha.to_degrees()))
            .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
    }
    /// (max L/D, CL at max L/D)
    pub fn ld_max(&self) -> Option<(f64, f64)> {
        self.points
            .iter()
            .filter(|p| p.cd > 1e-10)
            .map(|p| (p.cl / p.cd, p.cl))
            .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
    }
    /// CD at the point with the smallest |CL|
    pub fn cd0(&self) -> Option<f64> {
        self.points
            .iter()
            .min_by(|a, b| a.cl.abs().partial_cmp(&b.cl.abs()).unwrap())
            .map(|p| p.cd)
    }
}

/// The polar as XFOIL's OPER script runs it (CLAUDE.md):
/// `ALFA 0` / `ASEQ step alpha_max step` / `INIT` / `ALFA -step` / `ASEQ -2step alpha_min -step`,
/// with one persistent session so each point starts from the previous point's BL, and each
/// ASEQ halting after NSEQEX consecutive non-converged points. Points are stitched ascending.
pub fn compute_polar(airfoil: &PaneledAirfoil, config: &PolarConfig) -> PolarResult {
    compute_polar_with(airfoil, config, &mut |_, _| {})
}

/// [`compute_polar`] with an observer called after every point the sweep visits, converged or
/// not, with the session in that point's state (before the next SPECAL moves it on). This is how
/// `yfoil polar --distributions` captures the BL of each point: a point reached inside a sweep
/// starts from the previous alpha's BL and is not the same solve as a fresh `analyze`.
pub fn compute_polar_with(
    airfoil: &PaneledAirfoil,
    config: &PolarConfig,
    observe: &mut dyn FnMut(&Session, &OperatingPoint),
) -> PolarResult {
    let step = config.alpha_step.abs();
    let mut session = Session::new(airfoil, config.spec.clone());
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
    let record = |p: OperatingPoint, points: &mut Vec<OperatingPoint>, failed: &mut Vec<f64>| {
        if p.converged {
            points.push(p);
        } else {
            failed.push(p.alpha);
        }
    };
    let aseq = |session: &mut Session,
                alphas: Vec<f64>,
                points: &mut Vec<OperatingPoint>,
                failed: &mut Vec<f64>,
                observe: &mut dyn FnMut(&Session, &OperatingPoint)|
     -> bool {
        let mut iseqex = 0;
        for adeg in alphas {
            let p = session.aseq(adeg.to_radians());
            observe(session, &p);
            let conv = p.converged;
            record(p, points, failed);
            if session.st.lvisc && !conv {
                iseqex += 1;
                if iseqex >= config.nseqex {
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
    let p = session.alfa(0.0);
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
        let p = session.alfa(-step.to_radians());
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
        points,
        failed_alphas: failed,
        spec: config.spec.clone(),
        completed,
    }
}
