//! Influence coefficient calculations for vortex panel method
//!
//! This module implements the computation of influence coefficients
//! that relate vortex strengths on panels to streamfunction values
//! at control points.
//!
//! The approach follows XFOIL's PSILIN methodology but with a cleaner,
//! functional interface.

use std::f64::consts::PI;

/// Coefficient 1/(4*PI) used in influence calculations
const QOPI: f64 = 0.25 / PI;

/// Result of influence coefficient calculation for a single panel
/// on a single control point.
///
/// For a panel with linearly varying vorticity from γ₀ to γ₁,
/// the contribution to streamfunction is:
///   ΔΨ = QOPI * (psi_const * (γ₀ + γ₁) + psi_linear * (γ₁ - γ₀))
///
/// This gives:
///   dΨ/dγ₀ = QOPI * (psi_const - psi_linear)
///   dΨ/dγ₁ = QOPI * (psi_const + psi_linear)
#[derive(Debug, Clone, Copy)]
pub struct PanelInfluence {
    /// Constant vorticity contribution (PSIS in XFOIL)
    pub psi_const: f64,
    /// Linear vorticity variation contribution (PSID in XFOIL)
    pub psi_linear: f64,
}

impl PanelInfluence {
    /// Influence coefficient for the panel's first node (j)
    #[inline]
    pub fn coeff_j(&self) -> f64 {
        QOPI * (self.psi_const - self.psi_linear)
    }

    /// Influence coefficient for the panel's second node (j+1)
    #[inline]
    pub fn coeff_jp1(&self) -> f64 {
        QOPI * (self.psi_const + self.psi_linear)
    }
}

/// Calculate the influence of a vortex panel on a control point.
///
/// The panel extends from (x_j, y_j) to (x_jp1, y_jp1).
/// The control point is at (x_i, y_i).
///
/// # Arguments
/// * `x_j`, `y_j` - First node of panel (node j)
/// * `x_jp1`, `y_jp1` - Second node of panel (node j+1)
/// * `x_i`, `y_i` - Control point location
/// * `same_point_j` - True if control point coincides with node j
/// * `same_point_jp1` - True if control point coincides with node j+1
///
/// # Returns
/// Panel influence coefficients (psi_const, psi_linear)
pub fn panel_influence(
    x_j: f64,
    y_j: f64,
    x_jp1: f64,
    y_jp1: f64,
    x_i: f64,
    y_i: f64,
    same_point_j: bool,
    same_point_jp1: bool,
) -> PanelInfluence {
    // Panel length and unit tangent vector
    let dx = x_jp1 - x_j;
    let dy = y_jp1 - y_j;
    let ds = (dx * dx + dy * dy).sqrt();

    // Skip null panels
    if ds < 1e-14 {
        return PanelInfluence {
            psi_const: 0.0,
            psi_linear: 0.0,
        };
    }

    let ds_inv = 1.0 / ds;

    // Unit tangent vector along panel (s direction)
    let sx = dx * ds_inv;
    let sy = dy * ds_inv;

    // Vector from panel node j to control point
    let rx1 = x_i - x_j;
    let ry1 = y_i - y_j;

    // Vector from panel node j+1 to control point
    let rx2 = x_i - x_jp1;
    let ry2 = y_i - y_jp1;

    // Transform to panel-local coordinates
    // x1, x2: projections along panel
    // yy: perpendicular distance to panel line
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;

    // Squared distances to panel endpoints
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;

    // Log and arctan terms, handling singularities
    // Note: using x.atan2(y) gives atan2(x, y) in Rust
    // XFOIL uses ATAN2(Y, X) for angle of (X, Y) from positive x-axis
    let (g1, t1) = if same_point_j || rs1 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs1.ln(), x1.atan2(yy))
    };

    let (g2, t2) = if same_point_jp1 || rs2 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs2.ln(), x2.atan2(yy))
    };

    // Vortex panel contribution to streamfunction
    // For linearly varying vorticity γ = γ₀ + (γ₁-γ₀)*s/ds along panel
    //
    // PSIS: contribution from constant (average) part
    // PSID: contribution from linear variation part
    let psi_const = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

    // Linear variation term
    let dx_inv = x1 - x2;
    let psi_linear = if dx_inv.abs() < 1e-14 {
        // Panel is perpendicular to line from control point - use limit
        0.0
    } else {
        ((x1 + x2) * psi_const + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) / dx_inv
    };

    PanelInfluence {
        psi_const,
        psi_linear,
    }
}

/// Calculate the freestream contribution to streamfunction.
///
/// Ψ_∞ = Q_∞ * (cos(α)*y - sin(α)*x)
///
/// For unit freestream:
/// - At α=0°: Ψ = y
/// - At α=90°: Ψ = -x
#[inline]
pub fn freestream_psi(x: f64, y: f64, alpha_rad: f64, q_inf: f64) -> f64 {
    q_inf * (alpha_rad.cos() * y - alpha_rad.sin() * x)
}

/// Result of source panel influence calculation for a single panel
/// on a single control point.
///
/// For a panel with linearly varying source strength from σ₀ to σ₁,
/// the contribution to streamfunction is:
///   ΔΨ = QOPI * (psi_const * (σ₀ + σ₁) + psi_linear * (σ₁ - σ₀))
#[derive(Debug, Clone, Copy)]
pub struct SourcePanelInfluence {
    /// Constant source contribution
    pub psi_const: f64,
    /// Linear source variation contribution
    pub psi_linear: f64,
}

impl SourcePanelInfluence {
    /// Influence coefficient for the panel's first node (j)
    #[inline]
    pub fn coeff_j(&self) -> f64 {
        QOPI * (self.psi_const - self.psi_linear)
    }

    /// Influence coefficient for the panel's second node (j+1)
    #[inline]
    pub fn coeff_jp1(&self) -> f64 {
        QOPI * (self.psi_const + self.psi_linear)
    }
}

/// Calculate the influence of a source panel on a control point.
///
/// The panel extends from (x_j, y_j) to (x_jp1, y_jp1).
/// The control point is at (x_i, y_i).
///
/// This is the source analog of panel_influence() for vortex panels.
/// Source panels induce a flow perpendicular to the panel.
///
/// # Arguments
/// * `x_j`, `y_j` - First node of panel (node j)
/// * `x_jp1`, `y_jp1` - Second node of panel (node j+1)
/// * `x_i`, `y_i` - Control point location
/// * `same_point_j` - True if control point coincides with node j
/// * `same_point_jp1` - True if control point coincides with node j+1
///
/// # Returns
/// Source panel influence coefficients
pub fn source_panel_influence(
    x_j: f64,
    y_j: f64,
    x_jp1: f64,
    y_jp1: f64,
    x_i: f64,
    y_i: f64,
    same_point_j: bool,
    same_point_jp1: bool,
) -> SourcePanelInfluence {
    // Panel length and unit tangent vector
    let dx = x_jp1 - x_j;
    let dy = y_jp1 - y_j;
    let ds = (dx * dx + dy * dy).sqrt();

    // Skip null panels
    if ds < 1e-14 {
        return SourcePanelInfluence {
            psi_const: 0.0,
            psi_linear: 0.0,
        };
    }

    let ds_inv = 1.0 / ds;

    // Unit tangent vector along panel (s direction)
    let sx = dx * ds_inv;
    let sy = dy * ds_inv;

    // Vector from panel node j to control point
    let rx1 = x_i - x_j;
    let ry1 = y_i - y_j;

    // Vector from panel node j+1 to control point
    let rx2 = x_i - x_jp1;
    let ry2 = y_i - y_jp1;

    // Transform to panel-local coordinates
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;

    // Squared distances to panel endpoints
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;

    // Log and arctan terms, handling singularities
    let (g1, t1) = if same_point_j || rs1 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs1.ln(), x1.atan2(yy))
    };

    let (g2, t2) = if same_point_jp1 || rs2 < 1e-24 {
        (0.0, 0.0)
    } else {
        (rs2.ln(), x2.atan2(yy))
    };

    // Source panel contribution to streamfunction
    // For sources, the contribution is related to the log term differently
    // than for vortices (it's the "cross" term)
    //
    // XFOIL's PSWLIN: source influence on streamfunction
    // PSI = (yy*(g1-g2)/2 + x2*t2 - x1*t1) / (2*pi)
    let psi_const = 0.5 * yy * (g1 - g2) + x2 * t2 - x1 * t1;

    // Linear variation term
    let dx_inv = x1 - x2;
    let psi_linear = if dx_inv.abs() < 1e-14 {
        0.0
    } else {
        ((x1 + x2) * psi_const + 0.5 * (rs2 * t2 - rs1 * t1 + yy * (g1 - g2) * 0.5)) / dx_inv
    };

    SourcePanelInfluence {
        psi_const,
        psi_linear,
    }
}

/// Calculate the tangential velocity induced by a source panel at a field point.
///
/// This computes dQ/dσ for a source panel, which is needed for the DIJ matrix.
///
/// Returns (coeff_j, coeff_jp1) for the tangential velocity induced at (x_i, y_i)
/// in direction (tx, ty) by source strengths at nodes j and j+1.
pub fn source_panel_velocity(
    x_j: f64,
    y_j: f64,
    x_jp1: f64,
    y_jp1: f64,
    x_i: f64,
    y_i: f64,
    tx: f64,
    ty: f64,
    same_point_j: bool,
    same_point_jp1: bool,
) -> (f64, f64) {
    // Panel length and unit tangent vector
    let dx = x_jp1 - x_j;
    let dy = y_jp1 - y_j;
    let ds = (dx * dx + dy * dy).sqrt();

    // Skip null panels
    if ds < 1e-14 {
        return (0.0, 0.0);
    }

    let ds_inv = 1.0 / ds;

    // Unit tangent vector along panel (s direction)
    let sx = dx * ds_inv;
    let sy = dy * ds_inv;

    // Vector from panel nodes to control point
    let rx1 = x_i - x_j;
    let ry1 = y_i - y_j;
    let rx2 = x_i - x_jp1;
    let ry2 = y_i - y_jp1;

    // Transform to panel-local coordinates
    let x1 = sx * rx1 + sy * ry1;
    let x2 = sx * rx2 + sy * ry2;
    let yy = sx * ry1 - sy * rx1;

    // Squared distances to panel endpoints
    let rs1 = rx1 * rx1 + ry1 * ry1;
    let rs2 = rx2 * rx2 + ry2 * ry2;

    // Log and arctan terms, handling singularities
    let (g1, t1) = if same_point_j || rs1 < 1e-24 {
        (0.0, PI / 2.0) // Limit value when point coincides with node
    } else {
        (rs1.ln(), x1.atan2(yy))
    };

    let (g2, t2) = if same_point_jp1 || rs2 < 1e-24 {
        (0.0, PI / 2.0)
    } else {
        (rs2.ln(), x2.atan2(yy))
    };

    // Source panel induces velocity components:
    // u_s = (1/2π) * (θ₁ - θ₂)  (tangent to panel direction)
    // u_n = (1/4π) * ln(r₁²/r₂²) (normal to panel direction)
    //
    // The velocity in direction (tx, ty) is:
    // Q = u_panel_tangent * (tx*sx + ty*sy) + u_panel_normal * (tx*(-sy) + ty*sx)

    // Transform velocity direction to panel coordinates
    let cos_angle = tx * sx + ty * sy;
    let sin_angle = -tx * sy + ty * sx;

    // For constant source:
    // u_s = (t1 - t2) / (2π)
    // u_n = (g1 - g2) / (4π)
    let u_s_const = t1 - t2;
    let u_n_const = 0.5 * (g1 - g2);

    // Total tangential velocity from constant source part
    let q_const = QOPI * (u_s_const * cos_angle + u_n_const * sin_angle);

    // For linear source variation, compute the additional terms
    let dx_inv = x1 - x2;
    let q_linear = if dx_inv.abs() < 1e-14 {
        0.0
    } else {
        // Linear variation contribution (more complex - simplified here)
        let extra_s = if yy.abs() > 1e-14 {
            (x2 / rs2.max(1e-20) - x1 / rs1.max(1e-20)) * yy
        } else {
            0.0
        };
        let extra_n = if yy.abs() > 1e-14 {
            (yy / rs2.max(1e-20) - yy / rs1.max(1e-20)) * yy * 0.5 + (t2 - t1)
        } else {
            0.0
        };
        QOPI * ((x1 + x2) / dx_inv) * (extra_s * cos_angle + extra_n * sin_angle)
    };

    // Coefficients for node j and j+1
    let coeff_j = q_const - q_linear;
    let coeff_jp1 = q_const + q_linear;

    (coeff_j, coeff_jp1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_panel_influence_symmetric() {
        // Control point directly above panel midpoint
        let influence = panel_influence(0.0, 0.0, 1.0, 0.0, 0.5, 1.0, false, false);

        // psi_const should be non-zero
        assert!(influence.psi_const.abs() > 0.01);
        // For symmetric configuration, psi_linear should be near zero
        assert!(influence.psi_linear.abs() < 0.01);
    }

    #[test]
    fn test_panel_influence_at_node() {
        // Control point at panel node j
        let influence = panel_influence(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, true, false);

        // Should handle singularity gracefully
        assert!(influence.psi_const.is_finite());
        assert!(influence.psi_linear.is_finite());
    }

    #[test]
    fn test_panel_influence_null_panel() {
        // Null panel (both nodes same)
        let influence = panel_influence(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, false, false);

        assert_eq!(influence.psi_const, 0.0);
        assert_eq!(influence.psi_linear, 0.0);
    }

    #[test]
    fn test_freestream_psi_alpha_0() {
        // At alpha=0, Psi = Q*y
        let psi = freestream_psi(0.5, 0.3, 0.0, 1.0);
        assert_relative_eq!(psi, 0.3, epsilon = 1e-10);
    }

    #[test]
    fn test_freestream_psi_alpha_90() {
        // At alpha=90°, Psi = -Q*x
        let psi = freestream_psi(0.5, 0.3, PI / 2.0, 1.0);
        assert_relative_eq!(psi, -0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_coefficient_sum() {
        // For a panel, coeff_j + coeff_jp1 should equal 2*QOPI*psi_const
        let influence = panel_influence(0.0, 0.0, 1.0, 0.0, 0.5, 0.5, false, false);
        let sum = influence.coeff_j() + influence.coeff_jp1();
        assert_relative_eq!(sum, 2.0 * QOPI * influence.psi_const, epsilon = 1e-12);
    }

    #[test]
    fn test_coefficient_difference() {
        // coeff_jp1 - coeff_j should equal 2*QOPI*psi_linear
        let influence = panel_influence(0.0, 0.0, 1.0, 0.0, 0.3, 0.5, false, false);
        let diff = influence.coeff_jp1() - influence.coeff_j();
        assert_relative_eq!(diff, 2.0 * QOPI * influence.psi_linear, epsilon = 1e-12);
    }
}
