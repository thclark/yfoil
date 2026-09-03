//! Wake panel generation and influence calculations
//!
//! This module generates wake panels behind the trailing edge and computes
//! their influence on the airfoil and on each other. This is required for
//! accurate viscous-inviscid coupling.
//!
//! XFOIL Reference:
//! - XYWAKE in xpanel.f: Wake trajectory calculation
//! - QDCALC in xpanel.f: Wake source influence matrix
//! - PSWLIN in xpanel.f: Wake panel influence coefficients

use std::f64::consts::PI;

use crate::geometry::PaneledAirfoil;

/// Coefficient 1/(4*PI) used in influence calculations
const QOPI: f64 = 0.25 / PI;

/// Configuration for wake panel generation
#[derive(Debug, Clone, Copy)]
pub struct WakePanelConfig {
    /// Wake length in chord lengths (WAKLEN in XFOIL, default 1.0)
    pub wake_length: f64,
    /// Maximum number of wake panels (IWX in XFOIL)
    pub max_panels: usize,
}

impl Default for WakePanelConfig {
    fn default() -> Self {
        Self {
            wake_length: 1.0,
            max_panels: 50,
        }
    }
}

/// Wake panel geometry
///
/// Contains the coordinates and geometric properties of wake panels
/// extending from the trailing edge.
#[derive(Debug, Clone)]
pub struct WakePanels {
    /// x-coordinates of wake nodes (first node is at TE)
    pub x: Vec<f64>,
    /// y-coordinates of wake nodes
    pub y: Vec<f64>,
    /// Arc length from trailing edge
    pub s: Vec<f64>,
    /// x-component of panel normal (pointing in flow direction)
    pub nx: Vec<f64>,
    /// y-component of panel normal
    pub ny: Vec<f64>,
    /// Panel angle (angle of normal vector)
    pub apanel: Vec<f64>,
    /// Number of wake panels (= number of nodes - 1)
    pub n_wake: usize,
}

impl WakePanels {
    /// Create wake panels for a given airfoil
    ///
    /// This generates wake panels following the inviscid flow streamlines
    /// behind the trailing edge, similar to XFOIL's XYWAKE subroutine.
    ///
    /// # Arguments
    /// * `airfoil` - Paneled airfoil geometry
    /// * `gamma` - Vortex strength distribution (from inviscid solution)
    /// * `alpha` - Angle of attack (radians)
    /// * `config` - Wake configuration
    ///
    /// # Returns
    /// Wake panel geometry, or None if wake cannot be generated
    pub fn generate(airfoil: &PaneledAirfoil, gamma: &[f64], alpha: f64, config: &WakePanelConfig) -> Option<Self> {
        let n = airfoil.n;
        if n < 3 || gamma.len() != n {
            return None;
        }

        let chord = airfoil.chord;
        let s = &airfoil.s;

        // XFOIL formula for number of wake points:
        // NW = N/12 + 10*INT(WAKLEN)
        let n_wake = (n / 12 + 10 * config.wake_length as usize).min(config.max_panels);

        // TE coordinates (average of first and last node)
        let x_te = 0.5 * (airfoil.x[0] + airfoil.x[n - 1]);
        let y_te = 0.5 * (airfoil.y[0] + airfoil.y[n - 1]);

        // TE panel direction (from derivatives)
        // XFOIL: SX = 0.5*(YP(N) - YP(1)), SY = 0.5*(XP(1) - XP(N))
        // YP, XP are spline derivatives; approximate from panel endpoints
        let yp_1 = airfoil.y[1] - airfoil.y[0];
        let xp_1 = airfoil.x[1] - airfoil.x[0];
        let yp_n = airfoil.y[n - 1] - airfoil.y[n - 2];
        let xp_n = airfoil.x[n - 1] - airfoil.x[n - 2];

        let sx = 0.5 * (yp_n - yp_1);
        let sy = 0.5 * (xp_1 - xp_n);
        let smod = (sx * sx + sy * sy).sqrt();

        if smod < 1e-10 {
            return None;
        }

        // Initial wake direction (normal to bisector)
        let nx_init = sx / smod;
        let ny_init = sy / smod;

        // Generate wake spacing (exponential stretching like XFOIL's SETEXP)
        // DS1 = average of first and last surface panel sizes
        let ds1 = 0.5 * ((s[1] - s[0]) + (s[n - 1] - s[n - 2]));
        let wake_length_total = config.wake_length * chord;

        let spacing = generate_exponential_spacing(ds1, wake_length_total, n_wake);

        // Allocate arrays
        let mut x = Vec::with_capacity(n_wake + 1);
        let mut y = Vec::with_capacity(n_wake + 1);
        let mut s_wake = Vec::with_capacity(n_wake + 1);
        let mut nx = Vec::with_capacity(n_wake + 1);
        let mut ny = Vec::with_capacity(n_wake + 1);
        let mut apanel = Vec::with_capacity(n_wake);

        // First wake node: tiny distance behind TE
        // XFOIL: X(I) = XTE - 0.0001*NY(I), Y(I) = YTE + 0.0001*NX(I)
        x.push(x_te - 0.0001 * ny_init);
        y.push(y_te + 0.0001 * nx_init);
        s_wake.push(s[n - 1]); // Arc length continues from TE
        nx.push(nx_init);
        ny.push(ny_init);

        // For simplified implementation: follow freestream direction
        // (XFOIL uses PSILIN to follow streamlines, which requires the full
        // inviscid solution. We'll use a simpler approximation initially.)
        let cosa = alpha.cos();
        let sina = alpha.sin();

        // Wake follows freestream direction
        let wake_dx = cosa;
        let wake_dy = sina;

        // Generate remaining wake points
        for i in 0..n_wake {
            let ds = spacing[i];
            let x_new = x[i] + ds * wake_dx;
            let y_new = y[i] + ds * wake_dy;
            let s_new = s_wake[i] + ds;

            x.push(x_new);
            y.push(y_new);
            s_wake.push(s_new);

            // Normal vector (perpendicular to wake direction, pointing "up")
            // For wake following freestream: normal is (-sin α, cos α)
            nx.push(-sina);
            ny.push(cosa);

            // Panel angle
            apanel.push(alpha + PI / 2.0);
        }

        // Add final normal (same as previous)
        nx.push(-sina);
        ny.push(cosa);

        Some(Self {
            x,
            y,
            s: s_wake,
            nx,
            ny,
            apanel,
            n_wake,
        })
    }

    /// Number of wake nodes (n_wake + 1 for the first node at TE)
    pub fn num_nodes(&self) -> usize {
        self.x.len()
    }
}

/// Generate exponentially stretched spacing (like XFOIL's SETEXP)
///
/// Creates a sequence of spacings that start at `ds1` and grow
/// geometrically to cover total length `s_total` in `n` steps.
fn generate_exponential_spacing(ds1: f64, s_total: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }

    // For geometric series: s_total = ds1 * (1 + r + r^2 + ... + r^(n-1))
    //                               = ds1 * (r^n - 1) / (r - 1)
    // Solve for r (growth ratio)

    // Initial guess
    let mut ratio = 1.0 + s_total / (n as f64 * ds1);

    // Newton iteration to find ratio
    for _ in 0..20 {
        let sum = if (ratio - 1.0).abs() < 1e-10 {
            n as f64
        } else {
            (ratio.powi(n as i32) - 1.0) / (ratio - 1.0)
        };
        let target = s_total / ds1;

        if (sum - target).abs() < 1e-10 {
            break;
        }

        // Derivative of sum w.r.t. ratio
        let dsum = if (ratio - 1.0).abs() < 1e-10 {
            0.5 * (n * (n - 1)) as f64
        } else {
            let rn = ratio.powi(n as i32);
            (n as f64 * rn * (ratio - 1.0) - (rn - 1.0)) / ((ratio - 1.0) * (ratio - 1.0))
        };

        if dsum.abs() > 1e-10 {
            ratio -= (sum - target) / dsum;
            ratio = ratio.max(1.001).min(2.0);
        }
    }

    // Generate spacings
    let mut spacing = Vec::with_capacity(n);
    let mut ds = ds1;
    for _ in 0..n {
        spacing.push(ds);
        ds *= ratio;
    }

    spacing
}

/// Compute wake source influence on a control point
///
/// This implements the wake portion of XFOIL's PSWLIN subroutine.
/// Computes dPsi/dSig (DZDM) and dQtan/dSig (DQDM) for wake sources
/// affecting a control point.
///
/// # Arguments
/// * `x_i`, `y_i` - Control point coordinates
/// * `nx_i`, `ny_i` - Control point normal vector
/// * `wake` - Wake panel geometry
/// * `is_wake_point` - True if control point is on wake (changes sign convention)
///
/// # Returns
/// (dzdm, dqdm) - Influence vectors of length n_wake+1
pub fn wake_source_influence(
    x_i: f64,
    y_i: f64,
    nx_i: f64,
    ny_i: f64,
    wake: &WakePanels,
    is_wake_point: bool,
) -> (Vec<f64>, Vec<f64>) {
    let n_nodes = wake.num_nodes();
    let n_panels = wake.n_wake;

    let mut dzdm = vec![0.0; n_nodes];
    let mut dqdm = vec![0.0; n_nodes];

    if n_panels == 0 {
        return (dzdm, dqdm);
    }

    // Loop over wake panels (similar to PSWLIN loop)
    for jo in 0..n_panels {
        let jp = jo + 1;

        // Neighbor indices for interpolation
        let jm = if jo == 0 { jo } else { jo - 1 };
        let jq = if jo == n_panels - 1 { jp } else { jp + 1 };

        // Panel length
        let dx = wake.x[jp] - wake.x[jo];
        let dy = wake.y[jp] - wake.y[jo];
        let ds = (dx * dx + dy * dy).sqrt();
        if ds < 1e-14 {
            continue;
        }
        let ds_inv = 1.0 / ds;

        // Panel tangent
        let sx = dx * ds_inv;
        let sy = dy * ds_inv;

        // Vectors from panel nodes to control point
        let rx1 = x_i - wake.x[jo];
        let ry1 = y_i - wake.y[jo];
        let rx2 = x_i - wake.x[jp];
        let ry2 = y_i - wake.y[jp];

        // Transform to panel-local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        // Distances squared
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // Sign convention (different for wake points)
        let sgn = if is_wake_point { 1.0 } else { yy.signum() };

        // Logarithmic and angular terms
        let (g1, t1) = if rs1 > 1e-20 {
            let g = rs1.ln();
            let t = (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;
            (g, t)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if rs2 > 1e-20 {
            let g = rs2.ln();
            let t = (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;
            (g, t)
        } else {
            (0.0, 0.0)
        };

        // Normal vector components in panel coords
        let x1i = sx * nx_i + sy * ny_i;
        let yyi = sx * ny_i - sy * nx_i;

        // Panel angle
        let apan = wake.apanel[jo];

        // Midpoint quantities
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = rs0.max(1e-20).ln();
        let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * PI;

        // === First half-panel (1-0) ===
        let dx_inv = 1.0 / (x1 - x0).max(1e-10);
        let psum1 = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
        let pdif1 = ((x1 + x0) * psum1 + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dx_inv;

        // Derivatives for velocity
        let psx1 = -(t1 - apan);
        let psx0 = t0 - apan;
        let psyy1 = 0.5 * (g1 - g0);

        let psni1 = psx1 * x1i + psx0 * 0.5 * (x1i + x1i) + psyy1 * yyi;
        let pdx1 = ((x1 + x0) * psx1 + psum1 + 2.0 * x1 * (t1 - apan) - pdif1) * dx_inv;
        let pdx0_1 = ((x1 + x0) * psx0 + psum1 - 2.0 * x0 * (t0 - apan) + pdif1) * dx_inv;
        let pdyy1 = ((x1 + x0) * psyy1 + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dx_inv;
        let pdni1 = pdx1 * x1i + pdx0_1 * 0.5 * (x1i + x1i) + pdyy1 * yyi;

        // Neighbor panel length for interpolation
        let dsm = ((wake.x[jp] - wake.x[jm]).powi(2) + (wake.y[jp] - wake.y[jm]).powi(2)).sqrt();
        let dsm_inv = 1.0 / dsm.max(1e-10);

        // Accumulate dPsi/dm contributions (first half)
        dzdm[jm] += QOPI * (-psum1 * dsm_inv + pdif1 * dsm_inv);
        dzdm[jo] += QOPI * (-psum1 * ds_inv - pdif1 * ds_inv);
        dzdm[jp] += QOPI * (psum1 * (ds_inv + dsm_inv) + pdif1 * (ds_inv - dsm_inv));

        // Accumulate dQtan/dm contributions (first half)
        dqdm[jm] += QOPI * (-psni1 * dsm_inv + pdni1 * dsm_inv);
        dqdm[jo] += QOPI * (-psni1 * ds_inv - pdni1 * ds_inv);
        dqdm[jp] += QOPI * (psni1 * (ds_inv + dsm_inv) + pdni1 * (ds_inv - dsm_inv));

        // === Second half-panel (0-2) ===
        let dx_inv2 = 1.0 / (x0 - x2).max(1e-10);
        let psum2 = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
        let pdif2 = ((x0 + x2) * psum2 + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dx_inv2;

        // Derivatives for velocity
        let psx0_2 = -(t0 - apan);
        let psx2 = t2 - apan;
        let psyy2 = 0.5 * (g0 - g2);

        let psni2 = psx0_2 * 0.5 * (x1i + x1i) + psx2 * x1i + psyy2 * yyi;
        let pdx0_2 = ((x0 + x2) * psx0_2 + psum2 + 2.0 * x0 * (t0 - apan) - pdif2) * dx_inv2;
        let pdx2 = ((x0 + x2) * psx2 + psum2 - 2.0 * x2 * (t2 - apan) + pdif2) * dx_inv2;
        let pdyy2 = ((x0 + x2) * psyy2 + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dx_inv2;
        let pdni2 = pdx0_2 * 0.5 * (x1i + x1i) + pdx2 * x1i + pdyy2 * yyi;

        // Neighbor panel length for interpolation
        let dsp = ((wake.x[jq] - wake.x[jo]).powi(2) + (wake.y[jq] - wake.y[jo]).powi(2)).sqrt();
        let dsp_inv = 1.0 / dsp.max(1e-10);

        // Accumulate dPsi/dm contributions (second half)
        dzdm[jo] += QOPI * (-psum2 * ds_inv + pdif2 * ds_inv);
        dzdm[jp] += QOPI * (-psum2 * dsp_inv - pdif2 * dsp_inv);
        dzdm[jq] += QOPI * (psum2 * (ds_inv + dsp_inv) + pdif2 * (ds_inv - dsp_inv));

        // Accumulate dQtan/dm contributions (second half)
        dqdm[jo] += QOPI * (-psni2 * ds_inv + pdni2 * ds_inv);
        dqdm[jp] += QOPI * (-psni2 * dsp_inv - pdni2 * dsp_inv);
        dqdm[jq] += QOPI * (psni2 * (ds_inv + dsp_inv) + pdni2 * (ds_inv - dsp_inv));
    }

    (dzdm, dqdm)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exponential_spacing() {
        let ds1 = 0.01;
        let s_total = 1.0;
        let n = 20;

        let spacing = generate_exponential_spacing(ds1, s_total, n);

        assert_eq!(spacing.len(), n);

        // First spacing should be close to ds1
        assert!((spacing[0] - ds1).abs() < 1e-10);

        // Sum should be close to s_total
        let sum: f64 = spacing.iter().sum();
        assert!(
            (sum - s_total).abs() / s_total < 0.01,
            "Sum {} should be close to {}",
            sum,
            s_total
        );

        // Spacing should increase
        for i in 1..n {
            assert!(
                spacing[i] >= spacing[i - 1],
                "Spacing should increase: {} < {}",
                spacing[i],
                spacing[i - 1]
            );
        }
    }
}
