//! Result structures for JSON serialization
//!
//! These structs provide a clean, serializable interface for analysis results
//! that can be easily exported to JSON or other formats.

use serde::{Deserialize, Serialize};

/// Single operating point result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarPoint {
    /// Angle of attack (degrees)
    pub alpha_deg: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Moment coefficient (about quarter chord)
    pub cm: f64,
    /// Friction drag coefficient
    pub cd_friction: f64,
    /// Pressure drag coefficient
    pub cd_pressure: f64,
    /// Lift-to-drag ratio
    pub ldratio: f64,
    /// Transition location on upper surface (x/c)
    pub transition_upper: f64,
    /// Transition location on lower surface (x/c)
    pub transition_lower: f64,
    /// Whether solution converged
    pub converged: bool,
    /// Number of iterations to converge
    pub iterations: usize,
    /// Final RMSBL (rms BL Newton change) of the last VISCAL iteration
    pub residual: f64,
}

impl PolarPoint {
    /// Create from an analysis operating point
    pub fn from_point(p: &crate::solver::analysis::PointResult) -> Self {
        Self {
            alpha_deg: p.alpha.to_degrees(),
            cl: p.cl,
            cd: p.cd,
            cm: p.cm,
            cd_friction: p.cd_friction,
            cd_pressure: p.cd_pressure,
            ldratio: if p.cd > 1e-10 { p.cl / p.cd } else { 0.0 },
            transition_upper: p.transition_upper[0],
            transition_lower: p.transition_lower[0],
            converged: p.converged,
            iterations: p.iterations,
            residual: p.residual,
        }
    }
}

/// Flow conditions for analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowConditionsOutput {
    /// Reynolds number
    pub re: f64,
    /// Mach number
    pub mach: f64,
    /// Critical amplification factor (Ncrit)
    pub ncrit: f64,
}

impl FlowConditionsOutput {
    /// Create from the flow specification
    pub fn from_spec(spec: &crate::solver::analysis::FlowConditions) -> Self {
        Self {
            re: spec.re,
            mach: spec.mach,
            ncrit: spec.ncrit,
        }
    }
}

/// Polar sweep result (multiple operating points)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarOutput {
    /// Airfoil name/description
    pub foil: String,
    /// Display label for this polar (legend entry when plotted); falls back to `airfoil`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    /// Flow conditions
    pub conditions: FlowConditionsOutput,
    /// Operating points
    pub results: Vec<PolarPoint>,
    /// Summary statistics
    pub summary: PolarSummary,
    /// Whether sweep completed without excessive failures
    pub completed: bool,
    /// Full point records (geometry, wake and BL distributions) at every point the sweep visited,
    /// converged or not, ascending in alpha; written by `yfoil polar --distributions`. A point
    /// reached inside a sweep starts from the previous alpha's BL and is not the same solve as
    /// `yfoil analyze` at that alpha.
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
    /// Maximum L/D ratio
    pub ldratio_max: Option<f64>,
    /// CL at maximum L/D
    pub cl_at_ldratio_max: Option<f64>,
    /// Zero-lift drag coefficient
    pub cd0: Option<f64>,
    /// Number of converged points
    pub n_converged: usize,
    /// Number of failed points
    pub n_failed: usize,
}

impl PolarOutput {
    /// Create from PolarResult
    pub fn from_polar(result: &crate::solver::analysis::PolarResult, airfoil_name: &str) -> Self {
        let points: Vec<PolarPoint> = result.results.iter().map(PolarPoint::from_point).collect();

        let (cl_max, alpha_cl_max) = result.cl_max().map_or((None, None), |(cl, a)| (Some(cl), Some(a)));

        let (ld_max, cl_at_ld_max) = result
            .ldratio_max()
            .map_or((None, None), |(ld, cl)| (Some(ld), Some(cl)));

        let summary = PolarSummary {
            cl_max,
            alpha_at_cl_max: alpha_cl_max,
            ldratio_max: ld_max,
            cl_at_ldratio_max: cl_at_ld_max,
            cd0: result.cd0(),
            n_converged: result.results.iter().filter(|p| p.converged).count(),
            n_failed: result.failed_alphas.len(),
        };

        Self {
            foil: airfoil_name.to_string(),
            label: None,
            conditions: FlowConditionsOutput::from_spec(&result.conditions),
            results: points,
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

/// Single-point analysis result: forces, the geometry as solved (with the wake), and every
/// per-station boundary-layer quantity (see [`crate::output::BoundaryLayerOutput`]).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisOutput {
    /// Airfoil name/description
    pub foil: String,
    /// Flow conditions
    pub conditions: FlowConditionsOutput,
    /// Operating point result
    pub results: PolarPoint,
    /// Inviscid-only mode
    pub inviscid_only: bool,
    /// Panel nodes, normals and (after a viscous solve) the wake
    pub geometry: crate::output::FoilNodes,
    /// Boundary-layer distributions and markers; `None` for an inviscid point
    pub boundary_layer: Option<crate::output::BoundaryLayerOutput>,
}

impl AnalysisOutput {
    /// Create from a session's state after an operating point
    pub fn from_session(
        session: &crate::solver::analysis::Session,
        p: &crate::solver::analysis::PointResult,
        airfoil_name: &str,
        spec: &crate::solver::analysis::FlowConditions,
        inviscid_only: bool,
    ) -> Self {
        let st = &session.state;
        let boundary_layer = if st.viscous && st.bl_initialised {
            Some(crate::output::BoundaryLayerOutput::from_state(st))
        } else {
            None
        };
        Self {
            foil: airfoil_name.to_string(),
            conditions: FlowConditionsOutput::from_spec(spec),
            results: PolarPoint::from_point(p),
            inviscid_only,
            geometry: crate::output::FoilNodes::from_state(st),
            boundary_layer,
        }
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Inviscid analysis result with velocity distributions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InviscidAnalysisOutput {
    /// Airfoil name/description
    pub airfoil: String,
    /// Angle of attack in degrees
    pub alpha_deg: f64,
    /// Mach number
    pub mach: f64,
    /// Number of stations (panel nodes)
    pub n_stations: usize,
    /// Leading edge index (for splitting upper/lower surfaces)
    pub le_index: usize,
    /// Force coefficients
    pub coefficients: InviscidCoefficients,
    /// Station distributions
    pub stations: SurfaceDistributions,
}

/// Inviscid force coefficients
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InviscidCoefficients {
    /// Lift coefficient
    pub cl: f64,
    /// Moment coefficient (about quarter chord)
    pub cm: f64,
    /// Pressure drag coefficient
    pub cdp: f64,
}

/// Distributions at each station
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceDistributions {
    /// X-coordinates (x/c)
    pub x: Vec<f64>,
    /// Y-coordinates (y/c)
    pub y: Vec<f64>,
    /// Arc length parameter
    pub s: Vec<f64>,
    /// Surface velocity (normalized by freestream)
    pub velocity: Vec<f64>,
    /// Pressure coefficient
    pub cp: Vec<f64>,
}

impl InviscidAnalysisOutput {
    /// Create from inviscid solution and airfoil geometry
    pub fn new(
        airfoil: &crate::geometry::PanelledFoil,
        velocity: &[f64],
        cp: &[f64],
        (cl, cm, cdp): (f64, f64, f64),
        alpha_deg: f64,
        mach: f64,
        airfoil_name: &str,
    ) -> Self {
        Self {
            airfoil: airfoil_name.to_string(),
            alpha_deg,
            mach,
            n_stations: airfoil.n_foil_nodes,
            le_index: airfoil.i_le_node,
            coefficients: InviscidCoefficients { cl, cm, cdp },
            stations: SurfaceDistributions {
                x: airfoil.x.clone(),
                y: airfoil.y.clone(),
                s: airfoil.s.clone(),
                velocity: velocity.to_vec(),
                cp: cp.to_vec(),
            },
        }
    }

    /// Get upper surface data (TE to LE, indices 0..=le_index)
    ///
    /// Returns (x, cp, velocity) vectors for the upper surface
    pub fn upper_surface(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let le = self.le_index;
        let x: Vec<f64> = self.stations.x[0..=le].to_vec();
        let cp: Vec<f64> = self.stations.cp[0..=le].to_vec();
        let vel: Vec<f64> = self.stations.velocity[0..=le].to_vec();
        (x, cp, vel)
    }

    /// Get lower surface data (LE to TE, indices le_index..n)
    ///
    /// Returns (x, cp, velocity) vectors for the lower surface
    pub fn lower_surface(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let le = self.le_index;
        let x: Vec<f64> = self.stations.x[le..].to_vec();
        let cp: Vec<f64> = self.stations.cp[le..].to_vec();
        let vel: Vec<f64> = self.stations.velocity[le..].to_vec();
        (x, cp, vel)
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
            cd: 0.0085,
            cm: -0.05,
            cd_friction: 0.005,
            cd_pressure: 0.0035,
            ldratio: 64.7,
            transition_upper: 0.15,
            transition_lower: 0.45,
            converged: true,
            iterations: 12,
            residual: 1.2e-5,
        };

        let json = serde_json::to_string(&point).unwrap();
        assert!(json.contains("\"alpha_deg\":5.0"));
        assert!(json.contains("\"cl\":0.55"));
        assert!(json.contains("\"residual\":"));

        // Deserialize back
        let restored: PolarPoint = serde_json::from_str(&json).unwrap();
        assert!((restored.cl - 0.55).abs() < 1e-10);
        assert!((restored.residual - 1.2e-5).abs() < 1e-10);
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
        use crate::geometry::{naca_4digit, panel_foil};

        let geom = naca_4digit("0012", 160).unwrap();
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
        use crate::geometry::{naca_4digit, panel_foil};

        let geom = naca_4digit("0012", 120).unwrap();
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
