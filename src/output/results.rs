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
}
