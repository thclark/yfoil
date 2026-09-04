//! Panel method solver for inviscid flow
//!
//! This module assembles and solves the linear system for the vortex
//! panel method with linear vorticity distribution.

use nalgebra::{DMatrix, DVector};

use crate::geometry::PaneledAirfoil;

use super::influence::panel_influence;

use std::f64::consts::PI;

/// Coefficient 1/(4*PI) used in influence calculations
const QOPI: f64 = 0.25 / PI;

/// Coefficient 1/(2*PI) used in TE panel calculations
const HOPI: f64 = 0.5 / PI;

/// BWT factor for bisector control point distance (fraction of minimum TE panel length)
const BWT: f64 = 0.1;

/// Bisector condition data for sharp trailing edge
#[derive(Debug, Clone)]
#[allow(dead_code)] // x_bis/y_bis are XFOIL's XBIS/YBIS, kept for instrumentation parity
struct BisectorCondition {
    /// Control point x-coordinate on bisector
    x_bis: f64,
    /// Control point y-coordinate on bisector
    y_bis: f64,
    /// Cosine of bisector direction
    c_bis: f64,
    /// Sine of bisector direction
    s_bis: f64,
    /// DQDG[j] = dQ_bisector/dGamma[j] (velocity influence from vortex)
    dqdg: Vec<f64>,
    /// DQDM[j] = dQ_bisector/dMass[j] (velocity influence from source)
    dqdm: Vec<f64>,
}

/// Compute bisector condition for sharp trailing edge (XFOIL-style)
///
/// For a sharp TE, XFOIL replaces the flow tangency condition at node n-1
/// with a condition that the velocity along the bisector (between upper
/// and lower TE surfaces) is zero at a point just inside the TE.
///
/// This function computes:
/// 1. The bisector angle and control point location
/// 2. DQDG[j] = velocity influence from vortex strength at node j
/// 3. DQDM[j] = velocity influence from source strength at node j
fn compute_bisector_condition(airfoil: &PaneledAirfoil) -> BisectorCondition {
    let n = airfoil.n;

    // TE point (average of nodes 0 and n-1 for closed airfoil)
    let x_te = 0.5 * (airfoil.x[0] + airfoil.x[n - 1]);
    let y_te = 0.5 * (airfoil.y[0] + airfoil.y[n - 1]);

    // Compute panel tangent directions at TE
    // Upper surface: from node 0 to node 1 (reversed for XFOIL convention)
    let dx1 = airfoil.x[1] - airfoil.x[0];
    let dy1 = airfoil.y[1] - airfoil.y[0];
    let ds1 = (dx1 * dx1 + dy1 * dy1).sqrt();

    // Lower surface: from node n-2 to node n-1
    let dxn = airfoil.x[n - 1] - airfoil.x[n - 2];
    let dyn_te = airfoil.y[n - 1] - airfoil.y[n - 2];
    let dsn = (dxn * dxn + dyn_te * dyn_te).sqrt();

    // XFOIL convention:
    // AG1 = ATAN2(-YP(1), -XP(1)) - reversed tangent at upper TE
    // AG2 = ATAN2(YP(N), XP(N)) - tangent at lower TE
    let ag1 = (-dy1 / ds1).atan2(-dx1 / ds1);

    // For AG2, we need to handle the branch cut like XFOIL's ATANC
    // Simple version: use atan2 and adjust if needed
    let mut ag2 = (dyn_te / dsn).atan2(dxn / dsn);

    // Adjust for branch cut - ensure continuity with AG1
    while ag2 - ag1 > PI {
        ag2 -= 2.0 * PI;
    }
    while ag2 - ag1 < -PI {
        ag2 += 2.0 * PI;
    }

    // Bisector angle
    let a_bis = 0.5 * (ag1 + ag2);
    let c_bis = a_bis.cos();
    let s_bis = a_bis.sin();

    // Minimum panel length adjacent to TE
    let ds_min = ds1.min(dsn);

    // Control point on bisector just ahead of TE
    let x_bis = x_te - BWT * ds_min * c_bis;
    let y_bis = y_te - BWT * ds_min * s_bis;

    // Compute DQDG and DQDM at the bisector control point
    // Normal direction for tangential velocity: perpendicular to bisector
    let nx_bis = -s_bis;
    let ny_bis = c_bis;

    let (dqdg, dqdm) = compute_velocity_influence(airfoil, x_bis, y_bis, nx_bis, ny_bis);

    BisectorCondition {
        x_bis,
        y_bis,
        c_bis,
        s_bis,
        dqdg,
        dqdm,
    }
}

/// Compute velocity influence coefficients at an arbitrary point
///
/// Returns (DQDG, DQDM) where:
/// - DQDG[j] = dQ_tangent/dGamma[j] (velocity influence from vortex at node j)
/// - DQDM[j] = dQ_tangent/dMass[j] (velocity influence from source at node j)
///
/// The tangent direction is defined by the normal (nx, ny) such that
/// tangent = (-ny, nx) in XFOIL convention.
fn compute_velocity_influence(airfoil: &PaneledAirfoil, xi: f64, yi: f64, nxi: f64, nyi: f64) -> (Vec<f64>, Vec<f64>) {
    let n = airfoil.n;
    let mut dqdg = vec![0.0; n];
    let mut dqdm = vec![0.0; n];

    // Compute panel angles (APANEL) - same as elsewhere in this file
    let mut apanel = vec![0.0; n];
    for j in 0..n {
        let jp = if j == n - 1 { 0 } else { j + 1 };
        let sx = airfoil.x[jp] - airfoil.x[j];
        let sy = airfoil.y[jp] - airfoil.y[j];

        if sx == 0.0 && sy == 0.0 {
            apanel[j] = (-airfoil.ny[j]).atan2(-airfoil.nx[j]);
        } else if j == n - 1 {
            if airfoil.sharp_te {
                apanel[j] = PI;
            } else {
                apanel[j] = (-sx).atan2(sy) + PI;
            }
        } else {
            apanel[j] = sx.atan2(-sy);
        }
    }

    // Loop over all panels (JO to JP)
    for jo in 0..n {
        let jp = if jo == n - 1 { 0 } else { jo + 1 };

        let dso = ((airfoil.x[jp] - airfoil.x[jo]).powi(2) + (airfoil.y[jp] - airfoil.y[jo]).powi(2)).sqrt();
        if dso < 1e-14 {
            continue;
        }
        let dsio = 1.0 / dso;

        let apan = apanel[jo];

        // Vectors from panel nodes to control point
        let rx1 = xi - airfoil.x[jo];
        let ry1 = yi - airfoil.y[jo];
        let rx2 = xi - airfoil.x[jp];
        let ry2 = yi - airfoil.y[jp];

        // Unit tangent along panel
        let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
        let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

        // Transform to panel-local coordinates
        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        // Squared distances
        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // For off-airfoil point, use sign of yy to avoid branch cut issues
        let sgn = if yy >= 0.0 { 1.0 } else { -1.0 };

        // Log and arctan terms at panel endpoints
        let (g1, t1) = if rs1 > 1e-24 {
            let t = (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
            (rs1.ln(), t)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if rs2 > 1e-24 {
            let t = (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
            (rs2.ln(), t)
        } else {
            (0.0, 0.0)
        };

        // Direction derivatives (X1I, X2I, YYI in XFOIL)
        let x1i = sx * nxi + sy * nyi;
        let x2i = x1i; // Same as X1I in XFOIL
        let yyi = sx * nyi - sy * nxi;

        // Skip TE panel for source influence (XFOIL line 245)
        if jo != n - 1 {
            // Source influence (DQDM) computation
            let jm = if jo == 0 { jo } else { jo - 1 };
            let jq = if jo == n - 2 { n - 1 } else { (jp + 1) % n };

            // First half-panel (1 to 0)
            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = if rs0 > 1e-24 { rs0.ln() } else { 0.0 };
            let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;

            let dxinv_10 = if (x1 - x0).abs() > 1e-14 { 1.0 / (x1 - x0) } else { 0.0 };

            let psx1 = -(t1 - apan);
            let psx0 = t0 - apan;
            let psyy = 0.5 * (g1 - g0);

            let psum_10 = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdx1 = if dxinv_10.abs() > 1e-14 {
                ((x1 + x0) * psx1 + psum_10 + 2.0 * x1 * (t1 - apan)) * dxinv_10
            } else {
                0.0
            };
            let pdx0 = if dxinv_10.abs() > 1e-14 {
                ((x1 + x0) * psx0 + psum_10 - 2.0 * x0 * (t0 - apan)) * dxinv_10
            } else {
                0.0
            };
            let pdyy_10 = if dxinv_10.abs() > 1e-14 {
                ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv_10
            } else {
                0.0
            };

            let psni_10 = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni_10 = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy_10 * yyi;

            let dsm = ((airfoil.x[jp] - airfoil.x[jm]).powi(2) + (airfoil.y[jp] - airfoil.y[jm]).powi(2)).sqrt();
            let dsim = if dsm > 1e-14 { 1.0 / dsm } else { 0.0 };

            dqdm[jm] += QOPI * (-psni_10 * dsim + pdni_10 * dsim);
            dqdm[jo] += QOPI * (-psni_10 * dsio - pdni_10 * dsio);
            dqdm[jp] += QOPI * (psni_10 * (dsio + dsim) + pdni_10 * (dsio - dsim));

            // Second half-panel (0 to 2)
            let dxinv_02 = if (x0 - x2).abs() > 1e-14 { 1.0 / (x0 - x2) } else { 0.0 };

            let psx0_02 = -(t0 - apan);
            let psx2 = t2 - apan;
            let psyy_02 = 0.5 * (g0 - g2);

            let psum_02 = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdx0_02 = if dxinv_02.abs() > 1e-14 {
                ((x0 + x2) * psx0_02 + psum_02 + 2.0 * x0 * (t0 - apan)) * dxinv_02
            } else {
                0.0
            };
            let pdx2 = if dxinv_02.abs() > 1e-14 {
                ((x0 + x2) * psx2 + psum_02 - 2.0 * x2 * (t2 - apan)) * dxinv_02
            } else {
                0.0
            };
            let pdyy_02 = if dxinv_02.abs() > 1e-14 {
                ((x0 + x2) * psyy_02 + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv_02
            } else {
                0.0
            };

            let psni_02 = psx0_02 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy_02 * yyi;
            let pdni_02 = pdx0_02 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy_02 * yyi;

            let dsp = ((airfoil.x[jq] - airfoil.x[jo]).powi(2) + (airfoil.y[jq] - airfoil.y[jo]).powi(2)).sqrt();
            let dsip = if dsp > 1e-14 { 1.0 / dsp } else { 0.0 };

            dqdm[jo] += QOPI * (-psni_02 * (dsip + dsio) - pdni_02 * (dsip - dsio));
            dqdm[jp] += QOPI * (psni_02 * dsio - pdni_02 * dsio);
            dqdm[jq] += QOPI * (psni_02 * dsip + pdni_02 * dsip);
        }

        // Vortex influence (DQDG) computation
        let dxinv = if (x1 - x2).abs() > 1e-14 { 1.0 / (x1 - x2) } else { 0.0 };

        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let _psid = if dxinv.abs() > 1e-14 {
            ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv
        } else {
            0.0
        };

        let psx1_v = 0.5 * g1;
        let psx2_v = -0.5 * g2;
        let psyy_v = t1 - t2;

        let pdx1_v = if dxinv.abs() > 1e-14 {
            ((x1 + x2) * psx1_v + psis - x1 * g1 - _psid) * dxinv
        } else {
            0.0
        };
        let pdx2_v = if dxinv.abs() > 1e-14 {
            ((x1 + x2) * psx2_v + psis + x2 * g2 + _psid) * dxinv
        } else {
            0.0
        };
        let pdyy_v = if dxinv.abs() > 1e-14 {
            ((x1 + x2) * psyy_v - yy * (g1 - g2)) * dxinv
        } else {
            0.0
        };

        let psni_v = psx1_v * x1i + psx2_v * x2i + psyy_v * yyi;
        let pdni_v = pdx1_v * x1i + pdx2_v * x2i + pdyy_v * yyi;

        dqdg[jo] += QOPI * (psni_v - pdni_v);
        dqdg[jp] += QOPI * (psni_v + pdni_v);
    }

    // TE panel contribution to vortex influence (after the main loop, like XFOIL)
    let jo = n - 1;
    let jp = 0;

    let dso = ((airfoil.x[jp] - airfoil.x[jo]).powi(2) + (airfoil.y[jp] - airfoil.y[jo]).powi(2)).sqrt();

    if dso > 1e-14 {
        let dsio = 1.0 / dso;
        let apan = apanel[jo];

        let rx1 = xi - airfoil.x[jo];
        let ry1 = yi - airfoil.y[jo];
        let rx2 = xi - airfoil.x[jp];
        let ry2 = yi - airfoil.y[jp];

        let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
        let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn = if yy >= 0.0 { 1.0 } else { -1.0 };

        let (g1, t1) = if rs1 > 1e-24 {
            let t = (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
            (rs1.ln(), t)
        } else {
            (0.0, 0.0)
        };

        let (g2, t2) = if rs2 > 1e-24 {
            let t = (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
            (rs2.ln(), t)
        } else {
            (0.0, 0.0)
        };

        let x1i = sx * nxi + sy * nyi;
        let yyi = sx * nyi - sy * nxi;

        // TE panel uses SCS and SDS (sine and cosine of TE angle)
        // For sharp TE: SCS = 0, SDS = 1 (perpendicular panels)
        // In general: these depend on the TE gap angle
        let te_gap_x = airfoil.x[0] - airfoil.x[n - 1];
        let te_gap_y = airfoil.y[0] - airfoil.y[n - 1];
        let te_gap = (te_gap_x * te_gap_x + te_gap_y * te_gap_y).sqrt();

        let (scs, sds) = if te_gap < 1e-8 {
            // Sharp TE
            (0.0, 1.0)
        } else {
            // Blunt TE - compute from gap direction
            let chord = airfoil.chord;
            (te_gap_y / chord, te_gap_x / chord)
        };

        // PSIG and PGAM for TE panel
        let _psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        let _pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

        let psigx1 = -(t1 - apan);
        let psigx2 = t2 - apan;
        let psigyy = 0.5 * (g1 - g2);
        let pgamx1 = 0.5 * g1;
        let pgamx2 = -0.5 * g2;
        let pgamyy = t1 - t2;

        let psigni = psigx1 * x1i + psigx2 * x1i + psigyy * yyi;
        let pgamni = pgamx1 * x1i + pgamx2 * x1i + pgamyy * yyi;

        // XFOIL lines 446-447:
        // DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
        // DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
        dqdg[jo] -= HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
        dqdg[jp] += HOPI * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
    }

    (dqdg, dqdm)
}

/// Result of the inviscid panel method solution
#[derive(Debug, Clone)]
pub struct InviscidSolution {
    /// Vortex strength at each node for α=0° unit solution
    pub gam_0: Vec<f64>,
    /// Vortex strength at each node for α=90° unit solution
    pub gam_90: Vec<f64>,
    /// Surface velocity at panel midpoints for α=0° unit solution
    pub qinv_0: Vec<f64>,
    /// Surface velocity at panel midpoints for α=90° unit solution
    pub qinv_90: Vec<f64>,
    /// Internal streamfunction value for α=0°
    pub psi_0: f64,
    /// Internal streamfunction value for α=90°
    pub psi_90: f64,
    /// Number of airfoil panels
    pub n: usize,
    /// Source influence matrix DIJ for viscous-inviscid coupling
    /// For inviscid-only: n×n matrix (airfoil only)
    /// For viscous mode: (n+n_wake)×(n+n_wake) matrix including wake
    /// DIJ[i,j] gives velocity influence at point i from source at point j
    dij: Option<DMatrix<f64>>,
    /// Source influence on streamfunction BIJ = -DZDM (for debugging)
    bij: Option<DMatrix<f64>>,
}

impl InviscidSolution {
    /// Get vortex distribution at arbitrary angle of attack
    ///
    /// γ(α) = cos(α) * γ₀ + sin(α) * γ₉₀
    ///
    /// Note: This returns node-based vortex strengths, not surface velocities.
    /// For surface velocities, use `velocity_at_alpha` instead.
    pub fn gamma_at_alpha(&self, alpha_rad: f64) -> Vec<f64> {
        let cosa = alpha_rad.cos();
        let sina = alpha_rad.sin();
        self.gam_0
            .iter()
            .zip(&self.gam_90)
            .map(|(&g0, &g90)| cosa * g0 + sina * g90)
            .collect()
    }

    /// Get surface velocity at panel midpoints for arbitrary angle of attack
    ///
    /// Q(α) = cos(α) * Q₀ + sin(α) * Q₉₀
    ///
    /// These are computed at panel midpoints (control points) to avoid
    /// the trailing edge singularity that occurs at panel nodes.
    pub fn velocity_at_alpha(&self, alpha_rad: f64) -> Vec<f64> {
        let cosa = alpha_rad.cos();
        let sina = alpha_rad.sin();
        self.qinv_0
            .iter()
            .zip(&self.qinv_90)
            .map(|(&q0, &q90)| cosa * q0 + sina * q90)
            .collect()
    }

    /// Get surface velocity at nodes (original method, has TE singularity)
    ///
    /// This is the raw gamma distribution which for attached flow
    /// approximates the surface velocity.
    pub fn velocity_at_nodes(&self, alpha_rad: f64) -> Vec<f64> {
        self.gamma_at_alpha(alpha_rad)
    }

    /// Get the number of airfoil panels
    pub fn num_panels(&self) -> usize {
        self.n
    }

    /// Check if DIJ matrix has been computed
    pub fn has_dij(&self) -> bool {
        self.dij.is_some()
    }

    /// Get the source influence matrix DIJ
    ///
    /// DIJ[i,j] gives the velocity influence at panel i from a source at panel j.
    /// The source strength is related to mass defect: σ = d(ρ*Ue*δ*)/ds
    ///
    /// Returns None if DIJ hasn't been computed yet. Use `compute_dij()` first.
    pub fn get_dij(&self) -> Option<&DMatrix<f64>> {
        self.dij.as_ref()
    }

    /// Get the source influence on streamfunction BIJ = -DZDM (for debugging)
    ///
    /// BIJ[i,j] gives -dPsi/dSig at control point i from source at node j
    /// before the AIJ inversion.
    pub fn get_bij(&self) -> Option<&DMatrix<f64>> {
        self.bij.as_ref()
    }

    /// Compute velocity correction from mass defect using DIJ matrix
    ///
    /// Given mass defect m = Ue * δ* at each panel, computes the velocity
    /// correction following XFOIL's exact formula:
    ///   dU[i] = Σ_j (-VTI[i] * VTI[j] * DIJ[i,j] * MASS[j])
    ///
    /// where VTI is the sign conversion between panel and BL coordinates:
    /// - VTI = +1 for upper surface
    /// - VTI = -1 for lower surface
    ///
    /// # Arguments
    /// * `mass_defect` - Mass defect (Ue * δ*) at each panel
    /// * `le_index` - Leading edge index to determine upper/lower surface
    ///
    /// # Returns
    /// Velocity correction at each panel, or None if DIJ not computed
    pub fn velocity_from_mass_defect(&self, mass_defect: &[f64], le_index: usize) -> Option<Vec<f64>> {
        let dij = self.dij.as_ref()?;
        if mass_defect.len() != self.n {
            return None;
        }

        let mut dq = vec![0.0; self.n];

        // VTI sign conversion: +1 for upper surface, -1 for lower surface
        // In yfoil panel ordering: 0..le_index is upper, le_index+1..n-1 is lower
        let vti = |idx: usize| -> f64 {
            if idx <= le_index {
                1.0 // Upper surface
            } else {
                -1.0 // Lower surface
            }
        };

        // XFOIL formula: UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
        //                DUI = DUI + UE_M*MASS(JBL,JS)
        for i in 0..self.n {
            let vti_i = vti(i);
            for j in 0..self.n {
                let m = mass_defect[j];
                let d = dij[(i, j)];
                // Skip NaN or infinite values
                if m.is_finite() && d.is_finite() {
                    let vti_j = vti(j);
                    dq[i] += -vti_i * vti_j * d * m;
                }
            }
        }
        Some(dq)
    }

    /// Compute velocity correction from mass defect for UNSIGNED edge velocity
    ///
    /// This version computes the correction for the unsigned edge velocity |Ue|:
    ///   dUe[i] = -Σ_j DIJ[i,j] * MASS[j]
    ///
    /// The negative sign comes from the fact that displacement thickness
    /// (positive mass defect) reduces the effective velocity.
    pub fn velocity_from_mass_defect_unsigned(&self, mass_defect: &[f64]) -> Option<Vec<f64>> {
        let dij = self.dij.as_ref()?;
        if mass_defect.len() != self.n {
            return None;
        }

        let mut dq = vec![0.0; self.n];

        // Simple formula without VTI: dUe[i] = -Σ_j DIJ[i,j] * MASS[j]
        // The negative sign: positive mass defect (displacement) reduces velocity
        for i in 0..self.n {
            for j in 0..self.n {
                let m = mass_defect[j];
                let d = dij[(i, j)];
                if m.is_finite() && d.is_finite() {
                    dq[i] -= d * m;
                }
            }
        }
        Some(dq)
    }
}

/// Build and solve the inviscid panel method system.
///
/// This implements the vortex panel method with:
/// - Linear vorticity distribution on each panel
/// - Kutta condition at trailing edge: γ₁ + γₙ = 0
/// - Flow tangency (constant streamfunction on surface)
/// - Control points at panel MIDPOINTS (not nodes) to avoid singularity
///
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
///
/// # Returns
/// * `InviscidSolution` containing unit solutions for α=0° and α=90°
pub fn solve_inviscid(airfoil: &PaneledAirfoil) -> InviscidSolution {
    let n = airfoil.n;

    // Build influence matrix AIJ
    // Size is (n+1) x (n+1): n equations for Ψ=const, 1 for Kutta condition
    // Unknowns: γ₁...γₙ (vortex strengths) and Ψᵢₙₜ (internal streamfunction)
    let mut aij = DMatrix::<f64>::zeros(n + 1, n + 1);

    // Right-hand sides for α=0° and α=90°
    let mut rhs_0 = DVector::<f64>::zeros(n + 1);
    let mut rhs_90 = DVector::<f64>::zeros(n + 1);

    // Fill influence matrix for each control point (node i)
    for i in 0..n {
        let x_i = airfoil.x[i];
        let y_i = airfoil.y[i];

        // Sum influence from all panels
        for j in 0..n {
            // Panel j goes from node j to node j+1 (with wraparound for last panel)
            let jp1 = if j == n - 1 { 0 } else { j + 1 };

            let x_j = airfoil.x[j];
            let y_j = airfoil.y[j];
            let x_jp1 = airfoil.x[jp1];
            let y_jp1 = airfoil.y[jp1];

            // Check if control point coincides with panel nodes
            let same_j = i == j;
            let same_jp1 = i == jp1;

            // Skip the TE panel (j = n-1) - it uses special PSIG/PGAM formulas
            // computed after the main loop, not the regular PSIS/PSID formulas.
            // XFOIL line 245: IF(JO.EQ.N) GO TO 11 skips regular vortex influence
            if j == n - 1 {
                continue;
            }

            let influence = panel_influence(x_j, y_j, x_jp1, y_jp1, x_i, y_i, same_j, same_jp1);

            // dΨ/dγⱼ contribution
            aij[(i, j)] += influence.coeff_j();

            // dΨ/dγⱼ₊₁ contribution
            aij[(i, jp1)] += influence.coeff_jp1();
        }

        // dΨ/dΨᵢₙₜ = -1 (we want Ψ_surface = Ψ_internal)
        aij[(i, n)] = -1.0;

        // Freestream contributions (RHS)
        // At α=0°: Ψ_∞ = Q∞ * y, so RHS = -y (for unit Q∞)
        // At α=90°: Ψ_∞ = -Q∞ * x, so RHS = x
        rhs_0[i] = -y_i;
        rhs_90[i] = x_i;

        // TE panel contribution to DZDG (XFOIL lines 407-448)
        // This adds the effect of the TE gap panel (from node n-1 to node 0)
        // which uses PSIG/PGAM formulation with SCS/SDS TE angle factors.
        let te_gap_x = airfoil.x[0] - airfoil.x[n - 1];
        let te_gap_y = airfoil.y[0] - airfoil.y[n - 1];
        let te_gap = (te_gap_x * te_gap_x + te_gap_y * te_gap_y).sqrt();

        if te_gap > 1e-8 {
            // Compute SCS and SDS (TE angle factors)
            // XFOIL computes: ANTE = DXS*DYTE - DYS*DXTE (perpendicular to bisector)
            //                 ASTE = DXS*DXTE + DYS*DYTE (parallel to bisector)
            // where DXS, DYS is the bisector direction (average of reversed upper and lower tangents)

            // Compute bisector direction from panel tangents at TE
            // Upper surface tangent (from node 1 to node 0, reversed)
            let dx1 = airfoil.x[1] - airfoil.x[0];
            let dy1 = airfoil.y[1] - airfoil.y[0];
            let ds1 = (dx1 * dx1 + dy1 * dy1).sqrt();

            // Lower surface tangent (from node n-2 to node n-1)
            let dxn = airfoil.x[n - 1] - airfoil.x[n - 2];
            let dyn_te = airfoil.y[n - 1] - airfoil.y[n - 2];
            let dsn = (dxn * dxn + dyn_te * dyn_te).sqrt();

            // Bisector direction: DXS = 0.5*(-XP(1) + XP(N)), DYS = 0.5*(-YP(1) + YP(N))
            // where XP(1) = dx1/ds1 (tangent at upper TE), XP(N) = dxn/dsn (tangent at lower TE)
            let dxs = 0.5 * (-dx1 / ds1 + dxn / dsn);
            let dys = 0.5 * (-dy1 / ds1 + dyn_te / dsn);

            // TE gap vector: DXTE = X(1) - X(N), DYTE = Y(1) - Y(N)
            let dxte = te_gap_x;
            let dyte = te_gap_y;

            // ANTE, ASTE: perpendicular and parallel projections
            let ante = dxs * dyte - dys * dxte;
            let aste = dxs * dxte + dys * dyte;

            // SCS = ANTE/DSTE, SDS = ASTE/DSTE
            let scs = ante / te_gap;
            let sds = aste / te_gap;

            // TE panel geometry (from node n-1 to node 0)
            let jo = n - 1;
            let jp = 0;
            let dso = te_gap;
            let dsio = 1.0 / dso;

            // Vectors from TE panel nodes to control point
            let rx1 = x_i - airfoil.x[jo];
            let ry1 = y_i - airfoil.y[jo];
            let rx2 = x_i - airfoil.x[jp];
            let ry2 = y_i - airfoil.y[jp];

            // Unit tangent along TE panel
            let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
            let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

            // Transform to panel-local coordinates
            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;

            // Squared distances
            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;

            // Log and arctan terms (with coincidence check)
            // TE panel angle - compute from panel direction, not node tangent
            // XFOIL: For TE panel (j=N-1), APANEL = ATAN2(-SX, SY) + PI
            let apan = if airfoil.sharp_te { PI } else { (-sx).atan2(sy) + PI };
            let (g1, t1) = if i == jo || rs1 < 1e-24 {
                (0.0, 0.0)
            } else {
                (rs1.ln(), x1.atan2(yy))
            };
            let (g2, t2) = if i == jp || rs2 < 1e-24 {
                (0.0, 0.0)
            } else {
                (rs2.ln(), x2.atan2(yy))
            };

            // PSIG and PGAM for TE panel (XFOIL lines 408-409)
            let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
            let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

            // TE panel contribution to DZDG (XFOIL lines 434-438)
            // DZDG(JO) += -HOPI*PSIG*SCS*0.5 + HOPI*PGAM*SDS*0.5
            // DZDG(JP) += +HOPI*PSIG*SCS*0.5 - HOPI*PGAM*SDS*0.5
            aij[(i, jo)] += -HOPI * psig * scs * 0.5 + HOPI * pgam * sds * 0.5;
            aij[(i, jp)] += HOPI * psig * scs * 0.5 - HOPI * pgam * sds * 0.5;
        }
    }

    // For sharp trailing edge, replace the equation for node n-1 with
    // a bisector velocity condition (XFOIL-style).
    // This enforces zero internal velocity along the TE bisector at a
    // control point just inside the TE, which is a better-conditioned
    // physical boundary condition than curvature matching.
    //
    // The bisector condition: Q_bisector(XBIS, YBIS) = 0
    // where (XBIS, YBIS) is on the bisector between upper and lower TE surfaces.
    let bisector_condition = if airfoil.sharp_te && n >= 6 {
        let bc = compute_bisector_condition(airfoil);

        // Clear row n-1
        for j in 0..=n {
            aij[(n - 1, j)] = 0.0;
        }

        // Set AIJ row n-1 = DQDG (velocity influence from vortex strength)
        for j in 0..n {
            aij[(n - 1, j)] = bc.dqdg[j];
        }

        // No Psi_internal dependency for this velocity condition
        aij[(n - 1, n)] = 0.0;

        // RHS = -Q_freestream along bisector direction
        // For unit freestream at α=0°: Q = cos(0) * CBIS + sin(0) * SBIS = CBIS
        // For unit freestream at α=90°: Q = cos(90) * CBIS + sin(90) * SBIS = SBIS
        rhs_0[n - 1] = -bc.c_bis;
        rhs_90[n - 1] = -bc.s_bis;

        Some(bc)
    } else {
        None
    };

    // Kutta condition: γ₁ + γₙ = 0
    // Row n+1 (index n): enforce γ[0] + γ[n-1] = 0
    aij[(n, 0)] = 1.0;
    aij[(n, n - 1)] = 1.0;
    // RHS is 0 for Kutta condition
    rhs_0[n] = 0.0;
    rhs_90[n] = 0.0;

    // Solve the system using LU decomposition
    let lu = aij.clone().lu();

    let solution_0 = lu.solve(&rhs_0).expect("Panel system should be solvable for α=0°");
    let solution_90 = lu.solve(&rhs_90).expect("Panel system should be solvable for α=90°");

    // Extract vortex strengths
    let gam_0: Vec<f64> = solution_0.rows(0, n).iter().copied().collect();
    let gam_90: Vec<f64> = solution_90.rows(0, n).iter().copied().collect();
    let psi_0 = solution_0[n];
    let psi_90 = solution_90[n];

    // Compute velocities at panel midpoints
    // For the midpoint-based formulation, the gamma values at nodes
    // give the surface velocity when properly averaged
    let (qinv_0, qinv_90) = compute_midpoint_velocities(airfoil, &gam_0, &gam_90);

    // Compute source influence matrix DIJ for viscous-inviscid coupling
    // This matrix relates source strengths to velocity perturbations
    let (bij, dij) = compute_source_influence_matrix(airfoil, &aij, bisector_condition.as_ref());

    InviscidSolution {
        gam_0,
        gam_90,
        qinv_0,
        qinv_90,
        psi_0,
        psi_90,
        n,
        dij: Some(dij),
        bij: Some(bij),
    }
}

/// Compute the source influence matrices BIJ and DIJ.
///
/// DIJ[i,j] gives the tangential velocity at panel i induced by a unit
/// source at panel j. This is used for viscous-inviscid coupling where
/// the boundary layer mass defect acts as a source distribution.
///
/// The computation follows XFOIL's exact approach (PSILIN + QDCALC):
/// 1. For each control point i, compute DZDM[j] = dΨ/dσⱼ using XFOIL's
///    two-half-panel scheme with neighboring node interpolation
/// 2. Set BIJ[i,j] = -DZDM[j] (note the negation!)
/// 3. Solve DIJ = AIJ⁻¹ * BIJ via backsubstitution
///
/// Returns (BIJ, DIJ) tuple where BIJ is the pre-inversion matrix for debugging.
/// # Arguments
/// * `airfoil` - Paneled airfoil geometry
/// * `aij` - Already-built vortex influence matrix (n+1 x n+1)
/// * `bisector` - Optional bisector condition for sharp TE (provides DQDM for row n-1)
///
/// # Returns
/// Tuple (BIJ, DIJ) - both matrices of size n x n
fn compute_source_influence_matrix(
    airfoil: &PaneledAirfoil,
    aij: &DMatrix<f64>,
    bisector: Option<&BisectorCondition>,
) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = airfoil.n;

    // Compute panel angles APANEL[j] = angle of panel tangent
    // XFOIL's APCALC computes this from actual panel directions, not node normals:
    //   For j < n-1: APANEL(j) = ATAN2(SX, -SY) where SX,SY is panel direction
    //   For j = n-1 (TE panel): APANEL = ATAN2(-SX, SY) + PI (for non-sharp TE)
    let mut apanel = vec![0.0; n];
    for j in 0..n {
        let jp = if j == n - 1 { 0 } else { j + 1 };
        let sx = airfoil.x[jp] - airfoil.x[j];
        let sy = airfoil.y[jp] - airfoil.y[j];

        if sx == 0.0 && sy == 0.0 {
            // Null panel - use node normal as fallback
            apanel[j] = (-airfoil.ny[j]).atan2(-airfoil.nx[j]);
        } else if j == n - 1 {
            // TE panel - special handling
            if airfoil.sharp_te {
                apanel[j] = std::f64::consts::PI;
            } else {
                apanel[j] = (-sx).atan2(sy) + std::f64::consts::PI;
            }
        } else {
            // Normal panel
            apanel[j] = sx.atan2(-sy);
        }
    }

    // Build source influence matrix BIJ
    // BIJ[i,j] = -DZDM[j] at control point i
    let mut bij = DMatrix::<f64>::zeros(n + 1, n);

    // Distance tolerance for coincident points (may be used in future)
    let s_total = airfoil.s[n - 1] - airfoil.s[0];
    let _seps = s_total * 1.0e-5;

    for i in 0..n {
        // Control point location (at panel node i)
        let xi = airfoil.x[i];
        let yi = airfoil.y[i];

        // Compute DZDM for this control point
        let mut dzdm = vec![0.0; n];

        // Loop over all panels (JO to JP)
        // NOTE: XFOIL skips the TE panel (JO=N) for source influence - line 245 in xpanel.f
        // "IF(JO.EQ.N) GO TO 11" - this is critical for correct TE handling
        for jo in 0..n {
            let jp = if jo == n - 1 { 0 } else { jo + 1 };

            // Skip the TE panel for source influence (XFOIL line 245)
            // The TE panel does NOT contribute to DZDM in XFOIL
            if jo == n - 1 {
                continue;
            }

            // Neighboring node indices for source interpolation
            // XFOIL: JM = JO-1, JQ = JP+1, with special cases at boundaries
            let jm = if jo == 0 { jo } else { jo - 1 };
            let jq = if jo == n - 2 {
                n - 1 // XFOIL special case: JQ = JP for JO = N-1
            } else {
                (jp + 1) % n // Normal case with wrap-around
            };

            // Panel length
            let dso = ((airfoil.x[jp] - airfoil.x[jo]).powi(2) + (airfoil.y[jp] - airfoil.y[jo]).powi(2)).sqrt();
            if dso < 1e-14 {
                continue;
            }
            let dsio = 1.0 / dso;

            // Panel angle
            let apan = apanel[jo];

            // Vectors from panel nodes to control point
            let rx1 = xi - airfoil.x[jo];
            let ry1 = yi - airfoil.y[jo];
            let rx2 = xi - airfoil.x[jp];
            let ry2 = yi - airfoil.y[jp];

            // Unit tangent along panel
            let sx = (airfoil.x[jp] - airfoil.x[jo]) * dsio;
            let sy = (airfoil.y[jp] - airfoil.y[jo]) * dsio;

            // Transform to panel-local coordinates
            let x1 = sx * rx1 + sy * ry1;
            let x2 = sx * rx2 + sy * ry2;
            let yy = sx * ry1 - sy * rx1;

            // Squared distances
            let rs1 = rx1 * rx1 + ry1 * ry1;
            let rs2 = rx2 * rx2 + ry2 * ry2;

            // SGN for reflection handling (on airfoil surface, SGN = 1)
            let sgn = 1.0_f64;

            // Log and arctan terms at panel endpoints
            let (g1, t1) = if i == jo || rs1 < 1e-24 {
                (0.0, 0.0)
            } else {
                (rs1.ln(), (sgn * x1).atan2(sgn * yy))
            };

            let (g2, t2) = if i == jp || rs2 < 1e-24 {
                (0.0, 0.0)
            } else {
                (rs2.ln(), (sgn * x2).atan2(sgn * yy))
            };

            // Midpoint quantities
            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = if rs0 > 1e-24 { rs0.ln() } else { 0.0 };
            let t0 = (sgn * x0).atan2(sgn * yy);

            // ============ First half-panel (1 to 0) ============
            let dxinv_10 = if (x1 - x0).abs() > 1e-14 { 1.0 / (x1 - x0) } else { 0.0 };

            let psum_10 = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
            let pdif_10 = if dxinv_10.abs() > 1e-14 {
                ((x1 + x0) * psum_10 + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv_10
            } else {
                0.0
            };

            // Distance from JM to JP (spanning two panels)
            let dsm = ((airfoil.x[jp] - airfoil.x[jm]).powi(2) + (airfoil.y[jp] - airfoil.y[jm]).powi(2)).sqrt();
            let dsim = if dsm > 1e-14 { 1.0 / dsm } else { 0.0 };

            // dPsi/dm contributions for first half-panel
            // XFOIL: DZDM(JM) += QOPI*(-PSUM*DSIM + PDIF*DSIM)
            //        DZDM(JO) += QOPI*(-PSUM*DSIO - PDIF*DSIO)
            //        DZDM(JP) += QOPI*( PSUM*(DSIO+DSIM) + PDIF*(DSIO-DSIM))
            dzdm[jm] += QOPI * (-psum_10 * dsim + pdif_10 * dsim);
            dzdm[jo] += QOPI * (-psum_10 * dsio - pdif_10 * dsio);
            dzdm[jp] += QOPI * (psum_10 * (dsio + dsim) + pdif_10 * (dsio - dsim));

            // ============ Second half-panel (0 to 2) ============
            let dxinv_02 = if (x0 - x2).abs() > 1e-14 { 1.0 / (x0 - x2) } else { 0.0 };

            let psum_02 = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
            let pdif_02 = if dxinv_02.abs() > 1e-14 {
                ((x0 + x2) * psum_02 + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv_02
            } else {
                0.0
            };

            // Distance from JO to JQ (spanning two panels)
            let dsp = ((airfoil.x[jq] - airfoil.x[jo]).powi(2) + (airfoil.y[jq] - airfoil.y[jo]).powi(2)).sqrt();
            let dsip = if dsp > 1e-14 { 1.0 / dsp } else { 0.0 };

            // dPsi/dm contributions for second half-panel
            // XFOIL: DZDM(JO) += QOPI*(-PSUM*(DSIP+DSIO) - PDIF*(DSIP-DSIO))
            //        DZDM(JP) += QOPI*( PSUM*DSIO - PDIF*DSIO)
            //        DZDM(JQ) += QOPI*( PSUM*DSIP + PDIF*DSIP)
            dzdm[jo] += QOPI * (-psum_02 * (dsip + dsio) - pdif_02 * (dsip - dsio));
            dzdm[jp] += QOPI * (psum_02 * dsio - pdif_02 * dsio);
            dzdm[jq] += QOPI * (psum_02 * dsip + pdif_02 * dsip);
        }

        // Set BIJ[i,j] = -DZDM[j] (XFOIL convention)
        for j in 0..n {
            bij[(i, j)] = -dzdm[j];
        }
    }

    // Handle sharp TE condition for source influence
    // For bisector velocity condition, row n-1 gets -DQDM from the bisector point
    // (XFOIL lines 1102-1105: BIJ(N,J) = -DQDM(J))
    if let Some(bc) = bisector {
        for j in 0..n {
            bij[(n - 1, j)] = -bc.dqdm[j];
        }
    } else if airfoil.sharp_te && n >= 6 {
        // Fallback: zero source influence if no bisector condition provided
        for j in 0..n {
            bij[(n - 1, j)] = 0.0;
        }
    }

    // Row n (Kutta condition) has no source influence
    // (already zero from initialization)

    // Solve AIJ \ BIJ for each column to get dγ/dσ
    // This gives DIJ = dQtan/dSig
    let lu = aij.clone().lu();

    let mut dij = DMatrix::<f64>::zeros(n, n);

    for j in 0..n {
        let bij_col = bij.column(j).clone_owned();
        if let Some(dgam) = lu.solve(&bij_col) {
            // Store the vorticity change at each node
            // XFOIL stores DIJ(I,J) = BIJ(I,J) after backsubstitution
            // where BIJ has been overwritten with the solution
            for i in 0..n {
                dij[(i, j)] = dgam[i];
            }
        }
    }

    // Extract the n x n portion of BIJ for debugging output
    let bij_nxn = bij.view((0, 0), (n, n)).clone_owned();

    (bij_nxn, dij)
}

/// Compute surface velocities at panel midpoints.
///
/// For a vortex panel method, the solved vortex strength γ at each node
/// equals the tangential surface velocity there (by the properties of the
/// potential flow solution). At panel midpoints, we use the average of
/// the two node values.
///
/// Special handling for trailing edge:
/// - TE nodes (indices 0 and n-1) have singular γ values due to the Kutta condition
/// - Panels adjacent to the TE use extrapolation from nearby interior nodes
///
/// This approach follows XFOIL's treatment of the trailing edge singularity.
fn compute_midpoint_velocities(airfoil: &PaneledAirfoil, gam_0: &[f64], gam_90: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = airfoil.n;
    let mut qinv_0 = vec![0.0; n];
    let mut qinv_90 = vec![0.0; n];

    // For most panels, use average of node values
    for i in 0..n {
        let ip1 = if i == n - 1 { 0 } else { i + 1 };

        // Check if either node is a TE node (indices 0 or n-1)
        let node_i_is_te = i == 0 || i == n - 1;
        let node_ip1_is_te = ip1 == 0 || ip1 == n - 1;

        if (node_i_is_te || node_ip1_is_te) && airfoil.sharp_te {
            // Panel touches the TE - use only the non-TE node value
            if node_i_is_te && !node_ip1_is_te {
                // Node i is TE, use node ip1
                qinv_0[i] = gam_0[ip1];
                qinv_90[i] = gam_90[ip1];
            } else if !node_i_is_te && node_ip1_is_te {
                // Node ip1 is TE, use node i
                qinv_0[i] = gam_0[i];
                qinv_90[i] = gam_90[i];
            } else {
                // Both nodes are TE (the TE panel itself) - use average of
                // the two nodes adjacent to TE (nodes 1 and n-2)
                qinv_0[i] = 0.5 * (gam_0[1] + gam_0[n - 2]);
                qinv_90[i] = 0.5 * (gam_90[1] + gam_90[n - 2]);
            }
        } else {
            // Normal panels: average of the two node values
            qinv_0[i] = 0.5 * (gam_0[i] + gam_0[ip1]);
            qinv_90[i] = 0.5 * (gam_90[i] + gam_90[ip1]);
        }
    }

    (qinv_0, qinv_90)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{create_paneled_airfoil, naca_4digit};
    use approx::assert_relative_eq;

    #[test]
    fn test_solve_inviscid_naca0012() {
        let geom = naca_4digit("0012", 120).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Check solution dimensions
        assert_eq!(solution.gam_0.len(), airfoil.n);
        assert_eq!(solution.gam_90.len(), airfoil.n);

        // For symmetric airfoil at α=0°, the solution should be antisymmetric
        // Find max velocity magnitude excluding near-TE points (where singularities occur)
        let mut max_vel = 0.0;
        for i in 5..(airfoil.n - 5) {
            let vel = solution.gam_0[i].abs();
            if vel > max_vel {
                max_vel = vel;
            }
        }

        // Max velocity should be reasonable (around 1-2 for thin airfoil)
        assert!(max_vel > 0.5 && max_vel < 3.0, "Max velocity {} out of range", max_vel);
    }

    #[test]
    fn test_kutta_condition() {
        let geom = naca_4digit("0012", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Kutta condition: γ[0] + γ[n-1] ≈ 0
        let kutta_0 = solution.gam_0[0] + solution.gam_0[airfoil.n - 1];
        let kutta_90 = solution.gam_90[0] + solution.gam_90[airfoil.n - 1];

        assert!(kutta_0.abs() < 1e-10, "Kutta condition violated for α=0°: {}", kutta_0);
        assert!(
            kutta_90.abs() < 1e-10,
            "Kutta condition violated for α=90°: {}",
            kutta_90
        );
    }

    #[test]
    fn test_velocity_at_alpha() {
        let geom = naca_4digit("4412", 100).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // Test superposition of midpoint velocities
        let alpha = 0.1; // ~5.7 degrees
        let vel = solution.velocity_at_alpha(alpha);

        // Manual calculation for midpoint velocities
        let cosa = alpha.cos();
        let sina = alpha.sin();
        for i in 0..airfoil.n {
            let expected = cosa * solution.qinv_0[i] + sina * solution.qinv_90[i];
            assert_relative_eq!(vel[i], expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_stagnation_point_at_le() {
        let geom = naca_4digit("0012", 160).unwrap();
        let airfoil = create_paneled_airfoil(&geom);
        let solution = solve_inviscid(&airfoil);

        // At α=0° for symmetric airfoil, stagnation point should be at LE
        let vel = solution.velocity_at_alpha(0.0);
        let le_idx = airfoil.le_index;

        // Velocity at LE should be near zero (stagnation)
        assert!(
            vel[le_idx].abs() < 0.2,
            "LE velocity {} not near stagnation",
            vel[le_idx]
        );
    }
}
