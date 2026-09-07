//! The serialised results: one analysed point (`AnalysisOutput`), a polar (`PolarOutput`) and the
//! geometry report (`GeometryInfo`). Keys are the variable names of `docs/conventions/naming.md`.

use serde::{Deserialize, Serialize};

/// One operating point's results, as XFOIL reports them after `ALFA`/`CL` and keeps them in a
/// polar. The viscous-only fields are `None` for an inviscid point and are then omitted from the
/// JSON, so an inviscid record is a strict subset of a viscous one.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarPoint {
    /// Angle of attack (degrees)
    pub alpha_deg: f64,
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
    /// LVCONV: the viscous solution converged
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converged: Option<bool>,
    /// VISCAL iterations performed
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub iterations: Option<usize>,
    /// RMSBL of the last iteration
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub residual: Option<f64>,
}

impl PolarPoint {
    /// From a solved point; `viscous` selects whether the viscous-only fields are carried
    pub fn from_point(p: &crate::solver::analysis::PointResult, viscous: bool) -> Self {
        let v = viscous;
        Self {
            alpha_deg: p.alpha.to_degrees(),
            cl: p.cl,
            cm: p.cm,
            cd_pressure: p.cd_pressure,
            cd: v.then_some(p.cd),
            cd_friction: v.then_some(p.cd_friction),
            ldratio: v.then_some(if p.cd > 1e-10 { p.cl / p.cd } else { 0.0 }),
            transition_upper: v.then_some(p.transition_upper),
            transition_lower: v.then_some(p.transition_lower),
            converged: v.then_some(p.converged),
            iterations: v.then_some(p.iterations),
            residual: v.then_some(p.residual),
        }
    }

    /// True unless the point is a viscous one that did not converge
    pub fn is_converged(&self) -> bool {
        self.converged.unwrap_or(true)
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
    /// One record per converged point, ascending in alpha
    pub results: Vec<PolarPoint>,
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
    /// Number of converged points
    pub n_converged: usize,
    /// Number of failed points
    pub n_failed: usize,
}

impl PolarOutput {
    /// Create from PolarResult
    pub fn from_polar(result: &crate::solver::analysis::PolarResult, foil_name: &str) -> Self {
        let viscous = result.conditions.re.is_some();
        let results: Vec<PolarPoint> = result
            .results
            .iter()
            .map(|p| PolarPoint::from_point(p, viscous))
            .collect();

        let (cl_max, alpha_at_cl_max) = result.cl_max().map_or((None, None), |(cl, a)| (Some(cl), Some(a)));

        let (ldratio_max, cl_at_ldratio_max) = result
            .ldratio_max()
            .map_or((None, None), |(ld, cl)| (Some(ld), Some(cl)));

        let summary = PolarSummary {
            cl_max,
            alpha_at_cl_max,
            ldratio_max,
            cl_at_ldratio_max,
            cd0: result.cd0(),
            n_converged: result.results.iter().filter(|p| p.converged).count(),
            n_failed: result.failed_alphas.len(),
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
    pub results: PolarPoint,
    /// Panel nodes, normals and (after a viscous solve) the wake
    pub geometry: crate::output::FoilNodes,
    pub surface: SurfaceDistributions,
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
        Self {
            foil: foil_name.to_string(),
            conditions: session.conditions().clone(),
            results: PolarPoint::from_point(p, viscous),
            geometry: crate::output::FoilNodes::from_state(state),
            surface,
            boundary_layer,
        }
    }

    /// Upper-surface (x, cp, q), TE to LE (nodes `0..=i_le_node`)
    pub fn upper_surface(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let le = self.geometry.i_le_node;
        (
            self.geometry.x[0..=le].to_vec(),
            self.surface.cp[0..=le].to_vec(),
            self.surface.q[0..=le].to_vec(),
        )
    }

    /// Lower-surface (x, cp, q), LE to TE (nodes `i_le_node..`)
    pub fn lower_surface(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let le = self.geometry.i_le_node;
        (
            self.geometry.x[le..].to_vec(),
            self.surface.cp[le..].to_vec(),
            self.surface.q[le..].to_vec(),
        )
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
    fn test_operating_point_serialization() {
        let point = PolarPoint {
            alpha_deg: 5.0,
            cl: 0.55,
            cm: -0.05,
            cd_pressure: 0.0035,
            cd: Some(0.0085),
            cd_friction: Some(0.005),
            ldratio: Some(64.7),
            transition_upper: Some([0.15, 0.03]),
            transition_lower: Some([0.45, -0.02]),
            converged: Some(true),
            iterations: Some(12),
            residual: Some(1.2e-5),
        };

        let json = serde_json::to_string(&point).unwrap();
        assert!(json.contains("\"alpha_deg\":5.0"));
        assert!(json.contains("\"cl\":0.55"));
        assert!(json.contains("\"residual\":"));

        // Deserialize back
        let restored: PolarPoint = serde_json::from_str(&json).unwrap();
        assert!((restored.cl - 0.55).abs() < 1e-10);
        assert!((restored.residual.unwrap() - 1.2e-5).abs() < 1e-10);
    }

    #[test]
    fn test_polar_summary_serialization() {
        let summary = PolarSummary {
            cl_max: Some(1.2),
            alpha_at_cl_max: Some(12.0),
            ldratio_max: Some(80.0),
            cl_at_ldratio_max: Some(0.6),
            cd0: Some(0.006),
            n_converged: 25,
            n_failed: 2,
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
