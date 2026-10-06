//! The serialised results: one analysed point (`AnalysisOutput`), a polar (`PolarOutput`) and the
//! geometry report (`GeometryInfo`). Keys are the variable names of `docs/conventions/naming.md`.

use serde::{Deserialize, Serialize};

/// Why a point is not valid. Each variant is an exact, single-run fact — a formula evaluated
/// outside its analytic domain, or an iteration that ran out before meeting its tolerance.
/// Nothing here is inferred from the size or the smoothness of a returned number
/// (`docs/xfoil-known-issues.md` §7.10).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Reason {
    /// The Kármán–Tsien denominator `β + BFAC·Cp_inc` reached zero or below, so the
    /// compressibility correction has passed its pole and Cp has changed sign
    KarmanTsienOutOfDomain,
    /// SPECAL's CL(M) Newton used all 20 iterations without reaching `|DCLM| ≤ 1e-6`
    MachClNewtonExhausted,
    /// SPECCL's alpha Newton used all 20 iterations without reaching `|DALFA| ≤ 1e-6`
    ClNewtonExhausted,
    /// VISCAL finished without LVCONV: the viscous-inviscid iteration did not converge
    ViscousNotConverged,
    /// A station Newton failed with a non-finite residual and the march was carried over it by
    /// the garbage extrapolation: the run went through a NaN
    NonFiniteStationFailure,
    /// The sweep halted before reaching this alpha, so it was never attempted
    SequenceHalted,
}

/// Whether a point's numeric results can be believed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PointStatus {
    /// Solved, in domain, and every iteration met its tolerance
    Valid,
    /// Solved, but at least one `Reason` applies, so the numbers are not to be believed
    Invalid,
    /// Never solved: the sweep halted before reaching this alpha
    NotAttempted,
}

impl std::fmt::Display for Reason {
    /// The variant's own name, which is also what it serialises as: someone who reads a reason in
    /// a terminal or a JSON file and searches for it lands on the definition above.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Reason::KarmanTsienOutOfDomain => "KarmanTsienOutOfDomain",
            Reason::MachClNewtonExhausted => "MachClNewtonExhausted",
            Reason::ClNewtonExhausted => "ClNewtonExhausted",
            Reason::ViscousNotConverged => "ViscousNotConverged",
            Reason::NonFiniteStationFailure => "NonFiniteStationFailure",
            Reason::SequenceHalted => "SequenceHalted",
        };
        f.write_str(s)
    }
}

impl std::fmt::Display for PointStatus {
    /// As for [`Reason`]: the printed token is the serialised token is the variant name.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            PointStatus::Valid => "Valid",
            PointStatus::Invalid => "Invalid",
            PointStatus::NotAttempted => "NotAttempted",
        })
    }
}

/// The evidence for a point's status. Present whenever the point was attempted, and present
/// whether or not the numbers were withheld — these are facts about the computation, not
/// physical results, and withholding them would make the verdict unfalsifiable.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Diagnostics {
    /// Smallest Kármán–Tsien denominator over the force integration; at or below zero the
    /// correction has passed its pole. `None` if it was never evaluated
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub karman_tsien_margin_forces: Option<f64>,
    /// The same over the stored Cp arrays (CPCALC's `DENNEG` test, kept as a margin)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub karman_tsien_margin_pressure: Option<f64>,
    /// SPECAL's CL(M) Newton was exhausted
    pub mach_cl_newton_exhausted: bool,
    /// MRCL floored the lift coefficient when setting the reported Mach and Reynolds number
    pub cl_floored: bool,
    /// MRCL limited the reported Mach to 0.99
    pub mach_limited: bool,
    /// MRCL limited the reported Reynolds number to 100 × REINF1
    pub re_limited: bool,
    /// LVCONV: the viscous solution converged
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converged: Option<bool>,
    /// VISCAL iterations performed
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub iterations: Option<usize>,
    /// RMSBL of the last iteration
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub residual: Option<f64>,
    /// SPECCL's exit iteration for an OPER `CL` point (21 when its Newton was exhausted)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inviscid_cl_iterations: Option<usize>,
    /// Station Newton failures with a non-finite residual in this point's marches
    pub nonfinite_station_failures: usize,
}

/// One operating point's numbers, as XFOIL reports them after `ALFA`/`CL`. The viscous-only
/// fields are `None` for an inviscid point and are then omitted from the JSON, so an inviscid
/// record is a strict subset of a viscous one.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointValues {
    pub cl: f64,
    /// Moment coefficient about `cm_ref`
    pub cm: f64,
    /// Pressure drag (CDP): the only drag an inviscid solve has
    pub cd_pressure: f64,
    /// Total drag from the wake momentum defect (CD)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cd: Option<f64>,
    /// Friction drag (CDF)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cd_friction: Option<f64>,
    /// CL/CD
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ldratio: Option<f64>,
    /// (x, y) of the transition point on the upper side, chord fractions (XOCTR(1), YOCTR(1))
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transition_upper: Option<[f64; 2]>,
    /// (x, y) of the transition point on the lower side (XOCTR(2), YOCTR(2))
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transition_lower: Option<[f64; 2]>,
}

/// One entry of a result set: an operating point that was asked for. Every requested alpha
/// produces one of these, so a point is never silently missing — `status` says what happened to
/// it, `diagnostics` is the evidence, and `values` carries the numbers only when the point is
/// valid or the caller passed `--allow-invalid`. The same record is one polar entry and the
/// `results` of a single `analyse`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PointRecord {
    /// Angle of attack (degrees)
    pub alpha_deg: f64,
    /// This point's stable reference: a random 8-character id, unique to this solve. Absent when
    /// the point was never attempted. Random rather than positional so that points from different
    /// runs can be held together without their references colliding.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    /// The `id` of the record whose converged boundary layer seeded this solve, or absent when
    /// the BL was marched fresh. Records are presented ascending in alpha, which is not the order
    /// they were solved in; following this chain back recovers that order, and lets a sweep be
    /// restarted part way through instead of from the beginning.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initialised_from: Option<String>,
    pub status: PointStatus,
    /// Empty when the status is valid
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reasons: Vec<Reason>,
    /// Absent only when the point was never attempted
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub diagnostics: Option<Diagnostics>,
    /// The numbers, withheld unless the point is valid or `--allow-invalid` was given
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub values: Option<PointValues>,
}

impl PointRecord {
    /// Classify a solved point. `viscous` selects whether the viscous-only fields are carried;
    /// `allow_invalid` populates `values` even when the status is invalid, and never changes the
    /// status itself, so a consumer cannot read a withheld point as a sound one.
    pub fn from_point(p: &crate::solver::analysis::PointResult, viscous: bool, allow_invalid: bool) -> Self {
        let v = viscous;
        let mut reasons = Vec::new();
        if p.validity.out_of_domain() {
            reasons.push(Reason::KarmanTsienOutOfDomain);
        }
        if p.validity.mach_cl_newton_exhausted {
            reasons.push(Reason::MachClNewtonExhausted);
        }
        if p.inviscid_cl_iterations == 21 {
            reasons.push(Reason::ClNewtonExhausted);
        }
        if v && !p.converged {
            reasons.push(Reason::ViscousNotConverged);
        }
        if p.nonfinite_station_failures > 0 {
            reasons.push(Reason::NonFiniteStationFailure);
        }
        let status = if reasons.is_empty() {
            PointStatus::Valid
        } else {
            PointStatus::Invalid
        };
        let finite = |x: f64| x.is_finite().then_some(x);
        let values = (status == PointStatus::Valid || allow_invalid).then(|| PointValues {
            cl: p.cl,
            cm: p.cm,
            cd_pressure: p.cd_pressure,
            cd: v.then_some(p.cd),
            cd_friction: v.then_some(p.cd_friction),
            ldratio: v.then_some(if p.cd > 1e-10 { p.cl / p.cd } else { 0.0 }),
            transition_upper: v.then_some(p.transition_upper),
            transition_lower: v.then_some(p.transition_lower),
        });
        Self {
            alpha_deg: p.alpha.to_degrees(),
            id: Some(p.id.clone()),
            initialised_from: p.initialised_from.clone(),
            status,
            reasons,
            diagnostics: Some(Diagnostics {
                karman_tsien_margin_forces: finite(p.validity.karman_tsien_margin_forces),
                karman_tsien_margin_pressure: finite(p.validity.karman_tsien_margin_pressure),
                mach_cl_newton_exhausted: p.validity.mach_cl_newton_exhausted,
                cl_floored: p.validity.cl_floored,
                mach_limited: p.validity.mach_limited,
                re_limited: p.validity.re_limited,
                converged: v.then_some(p.converged),
                iterations: v.then_some(p.iterations),
                residual: v.then_some(p.residual),
                inviscid_cl_iterations: (p.inviscid_cl_iterations > 0).then_some(p.inviscid_cl_iterations),
                nonfinite_station_failures: p.nonfinite_station_failures,
            }),
            values,
        }
    }

    /// An alpha the sweep halted before reaching. It was asked for, so it is reported; it was
    /// never solved, so it has neither values nor evidence.
    pub fn not_attempted(alpha_deg: f64) -> Self {
        Self {
            alpha_deg,
            id: None,
            initialised_from: None,
            status: PointStatus::NotAttempted,
            reasons: vec![Reason::SequenceHalted],
            diagnostics: None,
            values: None,
        }
    }

    /// The point was solved and nothing disqualifies it.
    pub fn is_valid(&self) -> bool {
        self.status == PointStatus::Valid
    }
}

/// Polar sweep result (multiple operating points)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarOutput {
    /// Aerofoil name (the geometry file stem)
    pub foil: String,
    /// Display label for this polar (legend entry when plotted); falls back to `foil`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    /// The flow conditions every point shares
    pub conditions: crate::solver::analysis::FlowConditions,
    /// One record per alpha the sweep was asked for, ascending — including the ones that came
    /// back invalid and the ones it halted before reaching. Nothing requested is omitted.
    pub results: Vec<PointRecord>,
    /// Summary statistics
    pub summary: PolarSummary,
    /// Whether the sweep completed without excessive failures
    pub completed: bool,
    /// Full analysis records (geometry, wake and boundary layer) at every point the sweep
    /// visited, converged or not, ascending in alpha; written by `yfoil polar --distributions`.
    /// A point reached inside a sweep starts from the previous alpha's BL and is not the same
    /// solve as `yfoil analyse` at that alpha.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub distributions: Vec<AnalysisOutput>,
}

/// Summary statistics for a polar
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarSummary {
    /// Maximum lift coefficient
    pub cl_max: Option<f64>,
    /// Alpha at CL_max (degrees)
    pub alpha_at_cl_max: Option<f64>,
    /// Maximum L/D
    pub ldratio_max: Option<f64>,
    /// CL at maximum L/D
    pub cl_at_ldratio_max: Option<f64>,
    /// CD at the point of smallest |CL|
    pub cd0: Option<f64>,
    /// Number of points whose status is valid
    pub n_valid: usize,
    /// Number of points that were solved but came back invalid
    pub n_invalid: usize,
    /// Number of alphas the sweep halted before reaching
    pub n_not_attempted: usize,
}

impl PolarOutput {
    /// Create from PolarResult. `allow_invalid` carries the numbers of invalid points through
    /// as well; it never changes a status, and the summary is always computed from the valid
    /// points alone — a Kármán–Tsien artefact must not be able to become the reported CL_max.
    pub fn from_polar(result: &crate::solver::analysis::PolarResult, foil_name: &str, allow_invalid: bool) -> Self {
        let viscous = result.conditions.re.is_some();
        let mut results: Vec<PointRecord> = result
            .results
            .iter()
            .map(|p| PointRecord::from_point(p, viscous, allow_invalid))
            .chain(
                result
                    .not_attempted
                    .iter()
                    .map(|a| PointRecord::not_attempted(a.to_degrees())),
            )
            .collect();
        results.sort_by(|a, b| a.alpha_deg.partial_cmp(&b.alpha_deg).unwrap());

        // every statistic below is over the valid points only
        let valid: Vec<&PointRecord> = results.iter().filter(|r| r.is_valid()).collect();
        fn value(r: &PointRecord) -> &PointValues {
            r.values.as_ref().expect("a valid point carries its values")
        }
        let best = |f: &dyn Fn(&PointRecord) -> Option<(f64, f64)>| -> (Option<f64>, Option<f64>) {
            valid
                .iter()
                .filter_map(|r| f(r))
                .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
                .map_or((None, None), |(x, y)| (Some(x), Some(y)))
        };
        let (cl_max, alpha_at_cl_max) = best(&|r| Some((value(r).cl, r.alpha_deg)));
        let (ldratio_max, cl_at_ldratio_max) = best(&|r| {
            let v = value(r);
            v.cd.filter(|cd| *cd > 1e-10).map(|cd| (v.cl / cd, v.cl))
        });
        let cd0 = valid
            .iter()
            .min_by(|a, b| value(a).cl.abs().partial_cmp(&value(b).cl.abs()).unwrap())
            .map(|r| value(r).cd_pressure);

        let summary = PolarSummary {
            cl_max,
            alpha_at_cl_max,
            ldratio_max,
            cl_at_ldratio_max,
            cd0,
            n_valid: valid.len(),
            n_invalid: results.iter().filter(|r| r.status == PointStatus::Invalid).count(),
            n_not_attempted: result.not_attempted.len(),
        };

        Self {
            foil: foil_name.to_string(),
            label: None,
            conditions: result.conditions.clone(),
            results,
            summary,
            completed: result.completed,
            distributions: Vec::new(),
        }
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Per-node surface distributions of one solve: the surface speed q and Cp at every foil node,
/// inviscid or viscous according to `conditions.re`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SurfaceDistributions {
    /// Surface speed q/q∞ (QINV for an inviscid solve, QVIS for a viscous one)
    pub q: Vec<f64>,
    /// Pressure coefficient (CPI or CPV)
    pub cp: Vec<f64>,
}

/// One analysed point: the conditions, the results, the geometry as solved (with the wake after
/// a viscous solve), the surface distributions and, after a viscous solve, every per-station
/// boundary-layer quantity (see [`crate::output::BoundaryLayerOutput`]). An inviscid point has
/// no `geometry.wake` and no `boundary_layer`; everything else is the same shape.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisOutput {
    /// Aerofoil name (the geometry file stem)
    pub foil: String,
    pub conditions: crate::solver::analysis::FlowConditions,
    /// The point's status, evidence and — when it is valid, or `--allow-invalid` was given —
    /// its numbers. Exactly the record one polar entry carries.
    pub results: PointRecord,
    /// Panel nodes, normals and (after a viscous solve) the wake. Geometry is an *input*, not a
    /// result, so it is present whatever the status.
    pub geometry: crate::output::FoilNodes,
    /// Surface q and Cp. These are results — on an out-of-domain point the Cp field is exactly
    /// what passed through the Kármán–Tsien pole — so they are withheld on the same rule as
    /// `results.values`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub surface: Option<SurfaceDistributions>,
    /// Boundary-layer distributions and markers; absent for an inviscid point
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub boundary_layer: Option<crate::output::BoundaryLayerOutput>,
}

impl AnalysisOutput {
    /// From a session's state after an operating point. `include_lagged_closures` also emits
    /// XFOIL's lagged closure arrays under each side's `lagged_closures`
    /// (`yfoil analyse --include-lagged-closures`).
    pub fn from_session(
        session: &crate::solver::analysis::Session,
        p: &crate::solver::analysis::PointResult,
        foil_name: &str,
        include_lagged_closures: bool,
        allow_invalid: bool,
    ) -> Self {
        let state = session.state();
        let viscous = session.conditions().re.is_some();
        let n = state.n_foil_nodes;
        let surface = if viscous {
            SurfaceDistributions {
                q: state.q_viscous[1..=n].to_vec(),
                cp: state.cp_viscous[1..=n].to_vec(),
            }
        } else {
            SurfaceDistributions {
                q: state.q_inviscid[1..=n].to_vec(),
                cp: state.cp_inviscid[1..=n].to_vec(),
            }
        };
        let boundary_layer = if viscous && state.bl_initialised {
            Some(crate::output::BoundaryLayerOutput::from_state(
                state,
                include_lagged_closures,
            ))
        } else {
            None
        };
        let results = PointRecord::from_point(p, viscous, allow_invalid);
        // the numeric payload is gated as one: values, surface distributions and BL alike
        let disclose = results.is_valid() || allow_invalid;
        Self {
            foil: foil_name.to_string(),
            conditions: session.conditions().clone(),
            results,
            geometry: crate::output::FoilNodes::from_state(state),
            surface: disclose.then_some(surface),
            boundary_layer: boundary_layer.filter(|_| disclose),
        }
    }

    /// Upper-surface (x, cp, q), TE to LE (nodes `0..=i_le_node`). `None` when the surface
    /// distributions were withheld because the point is not valid.
    pub fn upper_surface(&self) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let le = self.geometry.i_le_node;
        let s = self.surface.as_ref()?;
        Some((
            self.geometry.x[0..=le].to_vec(),
            s.cp[0..=le].to_vec(),
            s.q[0..=le].to_vec(),
        ))
    }

    /// Lower-surface (x, cp, q), LE to TE (nodes `i_le_node..`). `None` as above.
    pub fn lower_surface(&self) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let le = self.geometry.i_le_node;
        let s = self.surface.as_ref()?;
        Some((self.geometry.x[le..].to_vec(), s.cp[le..].to_vec(), s.q[le..].to_vec()))
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Geometry information output
///
/// Contains comprehensive geometric properties of a paneled airfoil,
/// including both top-level summary statistics and detailed distributions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryInfo {
    /// Summary statistics (always included)
    pub summary: GeometrySummary,
    /// Detailed distributions (included in JSON output)
    pub distributions: GeometryDistributions,
}

/// Summary statistics for airfoil geometry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometrySummary {
    /// Number of panel nodes
    pub n_foil_nodes: usize,
    /// Chord length
    pub chord: f64,
    /// X-coordinate range [min, max]
    pub x_range: [f64; 2],
    /// Y-coordinate range [min, max]
    pub y_range: [f64; 2],
    /// Maximum thickness (max_y - min_y)
    pub y_extent: f64,
    /// Trailing edge gap (distance between first and last points)
    pub te_gap: f64,
    /// Whether trailing edge is sharp (gap < 0.01% chord)
    pub sharp_te: bool,
    /// Reference point for moment calculation [x/c, y/c]
    pub cm_ref: [f64; 2],
    /// Leading edge node index
    pub i_le_node: usize,
    /// Leading edge arc length parameter
    pub s_le: f64,
    /// Total arc length around the airfoil
    pub s_total: f64,
    /// Maximum curvature (typically at leading edge)
    pub max_curvature: f64,
    /// First point coordinates [x, y] (trailing edge upper)
    pub first_point: [f64; 2],
    /// Last point coordinates [x, y] (trailing edge lower)
    pub last_point: [f64; 2],
}

/// Detailed distributions along the airfoil surface
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryDistributions {
    /// X-coordinates at each node
    pub x: Vec<f64>,
    /// Y-coordinates at each node
    pub y: Vec<f64>,
    /// Arc length parameter at each node
    pub s: Vec<f64>,
    /// Curvature at each node
    pub curvature: Vec<f64>,
    /// Panel angle at each node (radians)
    pub panel_angle: Vec<f64>,
    /// Normal vector x-component at each node
    pub normal_x: Vec<f64>,
    /// Normal vector y-component at each node
    pub normal_y: Vec<f64>,
}

impl GeometryInfo {
    /// Create from a PanelledFoil
    pub fn from_panelled(airfoil: &crate::geometry::PanelledFoil) -> Self {
        let n = airfoil.n_foil_nodes;

        // Calculate ranges
        let min_x = airfoil.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_x = airfoil.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = airfoil.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_y = airfoil.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Calculate TE gap
        let te_gap = ((airfoil.x[0] - airfoil.x[n - 1]).powi(2) + (airfoil.y[0] - airfoil.y[n - 1]).powi(2)).sqrt();

        // Calculate curvature at each point
        let curvature = Self::calculate_curvature(airfoil);
        let max_curvature = curvature.iter().cloned().fold(0.0_f64, |m, v| m.max(v.abs()));

        let total_arc_length = airfoil.s[n - 1];

        let summary = GeometrySummary {
            n_foil_nodes: n,
            chord: airfoil.chord,
            x_range: [min_x, max_x],
            y_range: [min_y, max_y],
            y_extent: max_y - min_y,
            te_gap,
            sharp_te: airfoil.sharp_te,
            cm_ref: airfoil.cm_ref,
            i_le_node: airfoil.i_le_node,
            s_le: airfoil.s_le,
            s_total: total_arc_length,
            max_curvature,
            first_point: [airfoil.x[0], airfoil.y[0]],
            last_point: [airfoil.x[n - 1], airfoil.y[n - 1]],
        };

        let distributions = GeometryDistributions {
            x: airfoil.x.clone(),
            y: airfoil.y.clone(),
            s: airfoil.s.clone(),
            curvature,
            panel_angle: airfoil.panel_angle.clone(),
            normal_x: airfoil.normal_x.clone(),
            normal_y: airfoil.normal_y.clone(),
        };

        Self { summary, distributions }
    }

    /// Calculate curvature at each node
    ///
    /// Uses κ = dθ/ds where θ is the panel angle
    fn calculate_curvature(airfoil: &crate::geometry::PanelledFoil) -> Vec<f64> {
        let n = airfoil.n_foil_nodes;
        let mut curvature = vec![0.0; n];

        // Central differences for interior points
        for i in 1..n - 1 {
            let ds = airfoil.s[i + 1] - airfoil.s[i - 1];
            if ds > 1e-12 {
                // Handle angle wrap-around
                let mut dtheta = airfoil.panel_angle[i + 1] - airfoil.panel_angle[i - 1];
                if dtheta > std::f64::consts::PI {
                    dtheta -= 2.0 * std::f64::consts::PI;
                } else if dtheta < -std::f64::consts::PI {
                    dtheta += 2.0 * std::f64::consts::PI;
                }
                curvature[i] = dtheta / ds;
            }
        }

        // Forward difference for first point
        if n > 1 {
            let ds = airfoil.s[1] - airfoil.s[0];
            if ds > 1e-12 {
                let mut dtheta = airfoil.panel_angle[1] - airfoil.panel_angle[0];
                if dtheta > std::f64::consts::PI {
                    dtheta -= 2.0 * std::f64::consts::PI;
                } else if dtheta < -std::f64::consts::PI {
                    dtheta += 2.0 * std::f64::consts::PI;
                }
                curvature[0] = dtheta / ds;
            }
        }

        // Backward difference for last point
        if n > 1 {
            let ds = airfoil.s[n - 1] - airfoil.s[n - 2];
            if ds > 1e-12 {
                let mut dtheta = airfoil.panel_angle[n - 1] - airfoil.panel_angle[n - 2];
                if dtheta > std::f64::consts::PI {
                    dtheta -= 2.0 * std::f64::consts::PI;
                } else if dtheta < -std::f64::consts::PI {
                    dtheta += 2.0 * std::f64::consts::PI;
                }
                curvature[n - 1] = dtheta / ds;
            }
        }

        curvature
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_record_serialization() {
        let point = PointRecord {
            alpha_deg: 5.0,
            id: Some("k3n8qz1p".to_string()),
            initialised_from: Some("a7f2m0xd".to_string()),
            status: PointStatus::Valid,
            reasons: vec![],
            diagnostics: Some(Diagnostics {
                karman_tsien_margin_forces: Some(0.81),
                karman_tsien_margin_pressure: Some(0.81),
                mach_cl_newton_exhausted: false,
                cl_floored: false,
                mach_limited: false,
                re_limited: false,
                converged: Some(true),
                iterations: Some(12),
                residual: Some(1.2e-5),
                inviscid_cl_iterations: None,
                nonfinite_station_failures: 0,
            }),
            values: Some(PointValues {
                cl: 0.55,
                cm: -0.05,
                cd_pressure: 0.0035,
                cd: Some(0.0085),
                cd_friction: Some(0.005),
                ldratio: Some(64.7),
                transition_upper: Some([0.15, 0.03]),
                transition_lower: Some([0.45, -0.02]),
            }),
        };

        let json = serde_json::to_string(&point).unwrap();
        assert!(json.contains("\"alpha_deg\":5.0"));
        assert!(json.contains("\"status\":\"Valid\""));
        assert!(json.contains("\"cl\":0.55"));
        assert!(json.contains("\"residual\":"));

        let restored: PointRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.status, PointStatus::Valid);
        assert!((restored.values.unwrap().cl - 0.55).abs() < 1e-10);
    }

    #[test]
    fn test_withheld_point_still_reports_status_and_evidence() {
        // the contract: a point that was asked for comes back, carrying its verdict and the
        // evidence for it, with only the numbers withheld
        let point = PointRecord {
            alpha_deg: -1.0,
            id: Some("q9w4ez2r".to_string()),
            initialised_from: Some("t6y1u8io".to_string()),
            status: PointStatus::Invalid,
            reasons: vec![Reason::KarmanTsienOutOfDomain, Reason::MachClNewtonExhausted],
            diagnostics: Some(Diagnostics {
                karman_tsien_margin_forces: Some(-0.0296),
                karman_tsien_margin_pressure: Some(-0.0296),
                mach_cl_newton_exhausted: true,
                cl_floored: true,
                mach_limited: true,
                re_limited: false,
                converged: None,
                iterations: None,
                residual: None,
                inviscid_cl_iterations: None,
                nonfinite_station_failures: 0,
            }),
            values: None,
        };
        let json = serde_json::to_string(&point).unwrap();
        assert!(json.contains("\"status\":\"Invalid\""));
        assert!(json.contains("KarmanTsienOutOfDomain"));
        assert!(json.contains("MachClNewtonExhausted"));
        assert!(json.contains("-0.0296"), "the margin is evidence and is not withheld");
        assert!(!json.contains("\"values\""), "the numbers are withheld");
    }

    #[test]
    fn test_not_attempted_point_is_still_reported() {
        let p = PointRecord::not_attempted(24.0);
        assert_eq!(p.status, PointStatus::NotAttempted);
        assert_eq!(p.reasons, vec![Reason::SequenceHalted]);
        assert!(p.values.is_none() && p.diagnostics.is_none());
        let json = serde_json::to_string(&p).unwrap();
        assert!(json.contains("\"alpha_deg\":24.0"));
        assert!(json.contains("NotAttempted"));
    }

    #[test]
    fn test_polar_summary_serialization() {
        let summary = PolarSummary {
            cl_max: Some(1.2),
            alpha_at_cl_max: Some(12.0),
            ldratio_max: Some(80.0),
            cl_at_ldratio_max: Some(0.6),
            cd0: Some(0.006),
            n_valid: 25,
            n_invalid: 2,
            n_not_attempted: 0,
        };

        let json = serde_json::to_string_pretty(&summary).unwrap();
        assert!(json.contains("\"cl_max\""));
        assert!(json.contains("1.2"));
    }

    #[test]
    fn test_geometry_info_from_paneled() {
        use crate::geometry::{naca_4digit, panel_foil, Thickness};

        let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
        let airfoil = panel_foil(&geom);
        let info = GeometryInfo::from_panelled(&airfoil);

        // Check summary fields
        assert_eq!(info.summary.n_foil_nodes, 160);
        assert!((info.summary.chord - 1.0).abs() < 0.05);
        assert!(info.summary.x_range[0] < 0.01); // LE near x=0
        assert!(info.summary.x_range[1] > 0.99); // TE near x=1
        assert!(info.summary.y_extent > 0.10); // NACA 0012 has 12% thickness
        assert!(info.summary.y_extent < 0.14);
        assert!(info.summary.s_total > 1.8); // Arc length > chord
        assert!(info.summary.max_curvature > 0.0); // Should have positive curvature

        // Check distributions
        assert_eq!(info.distributions.x.len(), 160);
        assert_eq!(info.distributions.y.len(), 160);
        assert_eq!(info.distributions.s.len(), 160);
        assert_eq!(info.distributions.curvature.len(), 160);
        assert_eq!(info.distributions.panel_angle.len(), 160);
        assert_eq!(info.distributions.normal_x.len(), 160);
        assert_eq!(info.distributions.normal_y.len(), 160);

        // Arc length should be monotonically increasing
        for i in 1..info.distributions.s.len() {
            assert!(info.distributions.s[i] > info.distributions.s[i - 1]);
        }
    }

    #[test]
    fn test_geometry_info_serialization() {
        use crate::geometry::{naca_4digit, panel_foil, Thickness};

        let geom = naca_4digit("0012", 120, Thickness::Perpendicular).unwrap();
        let airfoil = panel_foil(&geom);
        let info = GeometryInfo::from_panelled(&airfoil);

        // Serialize to JSON
        let json = info.to_json().unwrap();

        // Check that key fields are present
        assert!(json.contains("\"n_foil_nodes\""));
        assert!(json.contains("\"chord\""));
        assert!(json.contains("\"x_range\""));
        assert!(json.contains("\"y_range\""));
        assert!(json.contains("\"y_extent\""));
        assert!(json.contains("\"te_gap\""));
        assert!(json.contains("\"sharp_te\""));
        assert!(json.contains("\"i_le_node\""));
        assert!(json.contains("\"s_le\""));
        assert!(json.contains("\"s_total\""));
        assert!(json.contains("\"max_curvature\""));
        assert!(json.contains("\"curvature\""));
        assert!(json.contains("\"panel_angle\""));

        // Deserialize back
        let restored: GeometryInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.summary.n_foil_nodes, info.summary.n_foil_nodes);
        assert!((restored.summary.chord - info.summary.chord).abs() < 1e-10);
        assert_eq!(restored.distributions.x.len(), info.distributions.x.len());
    }

    #[test]
    fn test_geometry_summary_serialization() {
        let summary = GeometrySummary {
            n_foil_nodes: 160,
            chord: 1.0,
            x_range: [0.0, 1.0],
            y_range: [-0.06, 0.06],
            y_extent: 0.12,
            te_gap: 0.00252,
            sharp_te: false,
            cm_ref: [0.25, 0.0],
            i_le_node: 80,
            s_le: 1.05,
            s_total: 2.1,
            max_curvature: 50.0,
            first_point: [1.0, 0.00126],
            last_point: [1.0, -0.00126],
        };

        let json = serde_json::to_string_pretty(&summary).unwrap();
        assert!(json.contains("\"n_foil_nodes\": 160"));
        assert!(json.contains("\"chord\": 1.0"));
        assert!(json.contains("\"sharp_te\": false"));

        // Deserialize back
        let restored: GeometrySummary = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.n_foil_nodes, 160);
        assert!((restored.y_extent - 0.12).abs() < 1e-10);
    }
}
