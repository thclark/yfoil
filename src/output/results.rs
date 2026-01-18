//! Result structures for JSON serialization
//!
//! These structs provide a clean, serializable interface for analysis results
//! that can be easily exported to JSON or other formats.

use serde::{Deserialize, Serialize};

/// Single operating point result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OperatingPoint {
    /// Angle of attack (degrees)
    pub alpha_deg: f64,
    /// Lift coefficient
    pub cl: f64,
    /// Drag coefficient
    pub cd: f64,
    /// Moment coefficient (about quarter chord)
    pub cm: f64,
    /// Friction drag coefficient
    pub cdf: f64,
    /// Pressure drag coefficient
    pub cdp: f64,
    /// Lift-to-drag ratio
    pub ld: f64,
    /// Transition location on upper surface (x/c)
    pub xtr_upper: f64,
    /// Transition location on lower surface (x/c)
    pub xtr_lower: f64,
    /// Whether solution converged
    pub converged: bool,
    /// Number of iterations to converge
    pub iterations: usize,
    /// Final convergence residual (CL change)
    pub residual: f64,
}

impl OperatingPoint {
    /// Create from ViscousResult
    pub fn from_viscous(result: &crate::solver::ViscousResult) -> Self {
        Self {
            alpha_deg: result.alpha.to_degrees(),
            cl: result.cl,
            cd: result.cd,
            cm: result.cm,
            cdf: result.cdf,
            cdp: result.cdp,
            ld: if result.cd > 1e-10 {
                result.cl / result.cd
            } else {
                0.0
            },
            xtr_upper: result.xtr_upper,
            xtr_lower: result.xtr_lower,
            converged: result.converged,
            iterations: result.iterations,
            residual: result.residual,
        }
    }
}

/// Flow conditions for analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowConditionsOutput {
    /// Reynolds number
    pub reynolds: f64,
    /// Mach number
    pub mach: f64,
    /// Critical amplification factor (Ncrit)
    pub ncrit: f64,
}

impl FlowConditionsOutput {
    /// Create from FlowConditions
    pub fn from_conditions(cond: &crate::bl::FlowConditions) -> Self {
        Self {
            reynolds: cond.reynolds,
            mach: cond.mach,
            ncrit: cond.ncrit,
        }
    }
}

/// Polar sweep result (multiple operating points)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarOutput {
    /// Airfoil name/description
    pub airfoil: String,
    /// Flow conditions
    pub conditions: FlowConditionsOutput,
    /// Operating points
    pub points: Vec<OperatingPoint>,
    /// Summary statistics
    pub summary: PolarSummary,
    /// Whether sweep completed without excessive failures
    pub completed: bool,
}

/// Summary statistics for a polar
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolarSummary {
    /// Maximum lift coefficient
    pub cl_max: Option<f64>,
    /// Alpha at CL_max (degrees)
    pub alpha_cl_max: Option<f64>,
    /// Maximum L/D ratio
    pub ld_max: Option<f64>,
    /// CL at maximum L/D
    pub cl_at_ld_max: Option<f64>,
    /// Zero-lift drag coefficient
    pub cd0: Option<f64>,
    /// Number of converged points
    pub num_converged: usize,
    /// Number of failed points
    pub num_failed: usize,
}

impl PolarOutput {
    /// Create from PolarResult
    pub fn from_polar(
        result: &crate::solver::PolarResult,
        airfoil_name: &str,
    ) -> Self {
        let points: Vec<OperatingPoint> = result
            .points
            .iter()
            .map(OperatingPoint::from_viscous)
            .collect();

        let (cl_max, alpha_cl_max) = result.cl_max().map_or((None, None), |(cl, a)| {
            (Some(cl), Some(a))
        });

        let (ld_max, cl_at_ld_max) = result.ld_max().map_or((None, None), |(ld, cl)| {
            (Some(ld), Some(cl))
        });

        let summary = PolarSummary {
            cl_max,
            alpha_cl_max,
            ld_max,
            cl_at_ld_max,
            cd0: result.cd0(),
            num_converged: result.points.iter().filter(|p| p.converged).count(),
            num_failed: result.failed_alphas.len(),
        };

        Self {
            airfoil: airfoil_name.to_string(),
            conditions: FlowConditionsOutput::from_conditions(&result.conditions),
            points,
            summary,
            completed: result.completed,
        }
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Single-point analysis result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisOutput {
    /// Airfoil name/description
    pub airfoil: String,
    /// Flow conditions
    pub conditions: FlowConditionsOutput,
    /// Operating point result
    pub result: OperatingPoint,
    /// Inviscid-only mode
    pub inviscid_only: bool,
}

impl AnalysisOutput {
    /// Create from ViscousResult
    pub fn from_viscous(
        result: &crate::solver::ViscousResult,
        airfoil_name: &str,
        conditions: &crate::bl::FlowConditions,
        inviscid_only: bool,
    ) -> Self {
        Self {
            airfoil: airfoil_name.to_string(),
            conditions: FlowConditionsOutput::from_conditions(conditions),
            result: OperatingPoint::from_viscous(result),
            inviscid_only,
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
    /// Force coefficients
    pub coefficients: InviscidCoefficients,
    /// Station distributions
    pub stations: StationDistributions,
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
pub struct StationDistributions {
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
        airfoil: &crate::geometry::PaneledAirfoil,
        velocity: &[f64],
        cp: &[f64],
        coeffs: &crate::forces::AeroCoefficients,
        alpha_deg: f64,
        mach: f64,
        airfoil_name: &str,
    ) -> Self {
        Self {
            airfoil: airfoil_name.to_string(),
            alpha_deg,
            mach,
            n_stations: airfoil.n,
            coefficients: InviscidCoefficients {
                cl: coeffs.cl,
                cm: coeffs.cm,
                cdp: coeffs.cdp,
            },
            stations: StationDistributions {
                x: airfoil.x.clone(),
                y: airfoil.y.clone(),
                s: airfoil.s.clone(),
                velocity: velocity.to_vec(),
                cp: cp.to_vec(),
            },
        }
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
    pub n_points: usize,
    /// Chord length
    pub chord: f64,
    /// X-coordinate range [min, max]
    pub x_range: [f64; 2],
    /// Y-coordinate range [min, max]
    pub y_range: [f64; 2],
    /// Maximum thickness (max_y - min_y)
    pub max_thickness: f64,
    /// Trailing edge gap (distance between first and last points)
    pub te_gap: f64,
    /// Whether trailing edge is sharp (gap < 0.01% chord)
    pub sharp_te: bool,
    /// Reference point for moment calculation [x/c, y/c]
    pub reference: [f64; 2],
    /// Leading edge node index
    pub le_index: usize,
    /// Leading edge arc length parameter
    pub sle: f64,
    /// Total arc length around the airfoil
    pub total_arc_length: f64,
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
    pub apanel: Vec<f64>,
    /// Normal vector x-component at each node
    pub nx: Vec<f64>,
    /// Normal vector y-component at each node
    pub ny: Vec<f64>,
}

impl GeometryInfo {
    /// Create from a PaneledAirfoil
    pub fn from_paneled(airfoil: &crate::geometry::PaneledAirfoil) -> Self {
        let n = airfoil.n;

        // Calculate ranges
        let min_x = airfoil.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_x = airfoil.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = airfoil.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_y = airfoil.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Calculate TE gap
        let te_gap = ((airfoil.x[0] - airfoil.x[n - 1]).powi(2)
            + (airfoil.y[0] - airfoil.y[n - 1]).powi(2))
        .sqrt();

        // Calculate curvature at each point
        let curvature = Self::calculate_curvature(airfoil);
        let max_curvature = curvature
            .iter()
            .cloned()
            .fold(0.0_f64, |m, v| m.max(v.abs()));

        let total_arc_length = airfoil.s[n - 1];

        let summary = GeometrySummary {
            n_points: n,
            chord: airfoil.chord,
            x_range: [min_x, max_x],
            y_range: [min_y, max_y],
            max_thickness: max_y - min_y,
            te_gap,
            sharp_te: airfoil.sharp_te,
            reference: airfoil.reference,
            le_index: airfoil.le_index,
            sle: airfoil.sle,
            total_arc_length,
            max_curvature,
            first_point: [airfoil.x[0], airfoil.y[0]],
            last_point: [airfoil.x[n - 1], airfoil.y[n - 1]],
        };

        let distributions = GeometryDistributions {
            x: airfoil.x.clone(),
            y: airfoil.y.clone(),
            s: airfoil.s.clone(),
            curvature,
            apanel: airfoil.apanel.clone(),
            nx: airfoil.nx.clone(),
            ny: airfoil.ny.clone(),
        };

        Self {
            summary,
            distributions,
        }
    }

    /// Calculate curvature at each node
    ///
    /// Uses κ = dθ/ds where θ is the panel angle
    fn calculate_curvature(airfoil: &crate::geometry::PaneledAirfoil) -> Vec<f64> {
        let n = airfoil.n;
        let mut curvature = vec![0.0; n];

        // Central differences for interior points
        for i in 1..n - 1 {
            let ds = airfoil.s[i + 1] - airfoil.s[i - 1];
            if ds > 1e-12 {
                // Handle angle wrap-around
                let mut dtheta = airfoil.apanel[i + 1] - airfoil.apanel[i - 1];
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
                let mut dtheta = airfoil.apanel[1] - airfoil.apanel[0];
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
                let mut dtheta = airfoil.apanel[n - 1] - airfoil.apanel[n - 2];
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
        let point = OperatingPoint {
            alpha_deg: 5.0,
            cl: 0.55,
            cd: 0.0085,
            cm: -0.05,
            cdf: 0.005,
            cdp: 0.0035,
            ld: 64.7,
            xtr_upper: 0.15,
            xtr_lower: 0.45,
            converged: true,
            iterations: 12,
            residual: 1.2e-5,
        };

        let json = serde_json::to_string(&point).unwrap();
        assert!(json.contains("\"alpha_deg\":5.0"));
        assert!(json.contains("\"cl\":0.55"));
        assert!(json.contains("\"residual\":"));

        // Deserialize back
        let restored: OperatingPoint = serde_json::from_str(&json).unwrap();
        assert!((restored.cl - 0.55).abs() < 1e-10);
        assert!((restored.residual - 1.2e-5).abs() < 1e-10);
    }

    #[test]
    fn test_polar_summary_serialization() {
        let summary = PolarSummary {
            cl_max: Some(1.2),
            alpha_cl_max: Some(12.0),
            ld_max: Some(80.0),
            cl_at_ld_max: Some(0.6),
            cd0: Some(0.006),
            num_converged: 25,
            num_failed: 2,
        };

        let json = serde_json::to_string_pretty(&summary).unwrap();
        assert!(json.contains("\"cl_max\""));
        assert!(json.contains("1.2"));
    }

    #[test]
    fn test_geometry_info_from_paneled() {
        use crate::geometry::{create_paneled_airfoil, naca_4digit};

        let geom = naca_4digit("0012", 160).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let info = GeometryInfo::from_paneled(&airfoil);

        // Check summary fields
        assert_eq!(info.summary.n_points, 160);
        assert!((info.summary.chord - 1.0).abs() < 0.05);
        assert!(info.summary.x_range[0] < 0.01); // LE near x=0
        assert!(info.summary.x_range[1] > 0.99); // TE near x=1
        assert!(info.summary.max_thickness > 0.10); // NACA 0012 has 12% thickness
        assert!(info.summary.max_thickness < 0.14);
        assert!(info.summary.total_arc_length > 1.8); // Arc length > chord
        assert!(info.summary.max_curvature > 0.0); // Should have positive curvature

        // Check distributions
        assert_eq!(info.distributions.x.len(), 160);
        assert_eq!(info.distributions.y.len(), 160);
        assert_eq!(info.distributions.s.len(), 160);
        assert_eq!(info.distributions.curvature.len(), 160);
        assert_eq!(info.distributions.apanel.len(), 160);
        assert_eq!(info.distributions.nx.len(), 160);
        assert_eq!(info.distributions.ny.len(), 160);

        // Arc length should be monotonically increasing
        for i in 1..info.distributions.s.len() {
            assert!(info.distributions.s[i] > info.distributions.s[i - 1]);
        }
    }

    #[test]
    fn test_geometry_info_serialization() {
        use crate::geometry::{create_paneled_airfoil, naca_4digit};

        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let info = GeometryInfo::from_paneled(&airfoil);

        // Serialize to JSON
        let json = info.to_json().unwrap();

        // Check that key fields are present
        assert!(json.contains("\"n_points\""));
        assert!(json.contains("\"chord\""));
        assert!(json.contains("\"x_range\""));
        assert!(json.contains("\"y_range\""));
        assert!(json.contains("\"max_thickness\""));
        assert!(json.contains("\"te_gap\""));
        assert!(json.contains("\"sharp_te\""));
        assert!(json.contains("\"le_index\""));
        assert!(json.contains("\"sle\""));
        assert!(json.contains("\"total_arc_length\""));
        assert!(json.contains("\"max_curvature\""));
        assert!(json.contains("\"curvature\""));
        assert!(json.contains("\"apanel\""));

        // Deserialize back
        let restored: GeometryInfo = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.summary.n_points, info.summary.n_points);
        assert!((restored.summary.chord - info.summary.chord).abs() < 1e-10);
        assert_eq!(restored.distributions.x.len(), info.distributions.x.len());
    }

    #[test]
    fn test_geometry_summary_serialization() {
        let summary = GeometrySummary {
            n_points: 160,
            chord: 1.0,
            x_range: [0.0, 1.0],
            y_range: [-0.06, 0.06],
            max_thickness: 0.12,
            te_gap: 0.00252,
            sharp_te: false,
            reference: [0.25, 0.0],
            le_index: 80,
            sle: 1.05,
            total_arc_length: 2.1,
            max_curvature: 50.0,
            first_point: [1.0, 0.00126],
            last_point: [1.0, -0.00126],
        };

        let json = serde_json::to_string_pretty(&summary).unwrap();
        assert!(json.contains("\"n_points\": 160"));
        assert!(json.contains("\"chord\": 1.0"));
        assert!(json.contains("\"sharp_te\": false"));

        // Deserialize back
        let restored: GeometrySummary = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.n_points, 160);
        assert!((restored.max_thickness - 0.12).abs() < 1e-10);
    }
}
