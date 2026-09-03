//! Panel distribution and repaneling
//!
//! Functions for redistributing panel points on an airfoil surface.

use super::airfoil::{Geometry, PaneledAirfoil};
use super::spline::{d2val, deval, seval, spline};

/// Configuration for XFOIL PANE algorithm
#[derive(Debug, Clone, Copy)]
pub struct PaneConfig {
    /// Curvature bunching parameter (default 1.0)
    /// Higher values = more bunching in high-curvature regions
    pub cvpar: f64,
    /// TE/LE panel density ratio (default 0.15)
    /// Ratio of artificial curvature at TE to LE curvature
    pub cterat: f64,
    /// Refinement area panel density ratio (default 0.2)
    pub ctrrat: f64,
    /// Refinement region on top surface (x/c range)
    pub xsref: Option<(f64, f64)>,
    /// Refinement region on bottom surface (x/c range)
    pub xpref: Option<(f64, f64)>,
}

impl Default for PaneConfig {
    fn default() -> Self {
        Self {
            cvpar: 1.0,
            cterat: 0.15,
            ctrrat: 0.2,
            xsref: None,
            xpref: None,
        }
    }
}

/// Repanel an airfoil using XFOIL's PANE algorithm (curvature-based)
///
/// This implements XFOIL's PANGEN subroutine which distributes panels
/// based on local curvature, placing more panels in high-curvature regions
/// (leading edge) and fewer in low-curvature regions (mid-chord).
///
/// # Arguments
/// * `geometry` - Input geometry
/// * `n_panels` - Target number of panels
/// * `config` - PANE configuration parameters
///
/// # Returns
/// New geometry with redistributed points matching XFOIL's PANE output
pub fn repanel_xfoil(geometry: &Geometry, n_panels: usize, config: &PaneConfig) -> Geometry {
    let nb = geometry.x_c.len();
    if nb < 2 {
        return geometry.clone();
    }

    // Calculate arc length along the buffer airfoil
    let sb = calculate_arc_length(&geometry.x_c, &geometry.y_c);

    // Spline the buffer airfoil coordinates
    let xbp = spline(&geometry.x_c, &sb);
    let ybp = spline(&geometry.y_c, &sb);

    // Normalizing length (~ chord)
    let sbref = 0.5 * (sb[nb - 1] - sb[0]);

    // Compute curvature at each buffer point
    let mut w5: Vec<f64> = (0..nb)
        .map(|i| curvature(sb[i], &geometry.x_c, &xbp, &geometry.y_c, &ybp, &sb).abs() * sbref)
        .collect();

    // Find LE point arc length and curvature
    let sble = find_le_arc_length(&geometry.x_c, &xbp, &geometry.y_c, &ybp, &sb);
    let cvle = curvature(sble, &geometry.x_c, &xbp, &geometry.y_c, &ybp, &sb).abs() * sbref;

    // TE coordinates
    let xbte = 0.5 * (geometry.x_c[0] + geometry.x_c[nb - 1]);
    let ybte = 0.5 * (geometry.y_c[0] + geometry.y_c[nb - 1]);

    // LE coordinates
    let xble = seval(sble, &geometry.x_c, &xbp, &sb);
    let yble = seval(sble, &geometry.y_c, &ybp, &sb);
    let chbsq = (xbte - xble).powi(2) + (ybte - yble).powi(2);

    // Set average curvature over region near LE
    let nk = 3;
    let mut cvsum = 0.0;
    for k in -nk..=nk {
        let frac = k as f64 / nk as f64;
        let sbk = sble + frac * sbref / cvle.max(20.0);
        let cvk = curvature(sbk, &geometry.x_c, &xbp, &geometry.y_c, &ybp, &sb).abs() * sbref;
        cvsum += cvk;
    }
    let cvavg = cvsum / (2 * nk + 1) as f64;

    // Curvature attraction coefficient
    let cc = 6.0 * config.cvpar;

    // Set artificial curvature at TE to bunch panels there
    let cvte = cvavg * config.cterat;
    w5[0] = cvte;
    w5[nb - 1] = cvte;

    // Smooth curvature array
    let smool = (1.0 / cvavg.max(20.0)).max(0.25 / (n_panels / 2) as f64);
    let smoosq = (smool * sbref).powi(2);

    // Set up tri-diagonal system for smoothed curvatures
    let mut w1 = vec![0.0; nb];
    let mut w2 = vec![0.0; nb];
    let mut w3 = vec![0.0; nb];

    w2[0] = 1.0;
    w3[0] = 0.0;

    for i in 1..nb - 1 {
        let dsm = sb[i] - sb[i - 1];
        let dsp = sb[i + 1] - sb[i];
        let dso = 0.5 * (sb[i + 1] - sb[i - 1]);

        if dsm == 0.0 || dsp == 0.0 {
            // Leave curvature at corner point unchanged
            w1[i] = 0.0;
            w2[i] = 1.0;
            w3[i] = 0.0;
        } else {
            w1[i] = smoosq * (-1.0 / dsm) / dso;
            w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3[i] = smoosq * (-1.0 / dsp) / dso;
        }
    }

    w1[nb - 1] = 0.0;
    w2[nb - 1] = 1.0;

    // Fix curvature at LE point
    for i in 1..nb - 1 {
        if (sb[i] - sble).abs() < 1e-10 {
            w1[i] = 0.0;
            w2[i] = 1.0;
            w3[i] = 0.0;
            w5[i] = cvle;
        } else if sb[i - 1] < sble && sb[i] > sble {
            // Modify equation at node just after LE point
            let dsm = sb[i] - sble;
            let dsp = sb[i + 1] - sb[i];
            let dso = 0.5 * (sb[i + 1] - sble);
            w1[i] = 0.0;
            w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3[i] = smoosq * (-1.0 / dsp) / dso;
            w5[i] = w5[i] + smoosq * cvle / (dsm * dso);
            break;
        }
    }

    // Set artificial curvature at refinement regions
    if let Some((xsref1, xsref2)) = config.xsref {
        for i in 1..nb - 1 {
            let xoc = ((geometry.x_c[i] - xble) * (xbte - xble) + (geometry.y_c[i] - yble) * (ybte - yble)) / chbsq;
            if sb[i] < sble && xoc > xsref1 && xoc < xsref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * config.ctrrat;
            }
        }
    }

    if let Some((xpref1, xpref2)) = config.xpref {
        for i in 1..nb - 1 {
            let xoc = ((geometry.x_c[i] - xble) * (xbte - xble) + (geometry.y_c[i] - yble) * (ybte - yble)) / chbsq;
            if sb[i] >= sble && xoc > xpref1 && xoc < xpref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * config.ctrrat;
            }
        }
    }

    // Solve tri-diagonal system for smoothed curvature
    w5 = trisol_curvature(&w1, &w2, &w3, &w5);

    // Normalize curvature array
    let cvmax = w5.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    if cvmax > 0.0 {
        for v in &mut w5 {
            *v /= cvmax;
        }
    }

    // Spline the normalized curvature
    let w6 = spline(&w5, &sb);

    // Set initial guess for node positions
    // Use more nodes than specified for more reliable convergence
    let ipfac = 5;
    let nn = ipfac * (n_panels - 1) + 1;

    // Ratio of panel lengths at TE
    let rdste = 0.667;
    let rtf = (rdste - 1.0) * 2.0 + 1.0;

    let dsavg = (sb[nb - 1] - sb[0]) / ((nn - 3) as f64 + 2.0 * rtf);
    let mut snew = vec![0.0; nn];
    snew[0] = sb[0];
    for i in 1..nn - 1 {
        snew[i] = sb[0] + dsavg * ((i - 1) as f64 + rtf);
    }
    snew[nn - 1] = sb[nb - 1];

    // Newton iteration for new node positions
    for _iter in 0..20 {
        let mut w1_n = vec![0.0; nn];
        let mut w2_n = vec![0.0; nn];
        let mut w3_n = vec![0.0; nn];
        let mut w4 = vec![0.0; nn];

        let cv1 = seval(snew[0], &w5, &w6, &sb);
        let mut cv2 = seval(snew[1], &w5, &w6, &sb);
        let cvs1 = deval(snew[0], &w5, &w6, &sb);
        let mut cvs2 = deval(snew[1], &w5, &w6, &sb);

        let mut cavm = (cv1.powi(2) + cv2.powi(2)).sqrt();
        let (mut cavm_s1, mut cavm_s2) = if cavm == 0.0 {
            (0.0, 0.0)
        } else {
            (cvs1 * cv1 / cavm, cvs2 * cv2 / cavm)
        };

        for i in 1..nn - 1 {
            let dsm = snew[i] - snew[i - 1];
            let dsp = snew[i] - snew[i + 1];
            let cv3 = seval(snew[i + 1], &w5, &w6, &sb);
            let cvs3 = deval(snew[i + 1], &w5, &w6, &sb);

            let cavp = (cv3.powi(2) + cv2.powi(2)).sqrt();
            let (cavp_s2, cavp_s3) = if cavp == 0.0 {
                (0.0, 0.0)
            } else {
                (cvs2 * cv2 / cavp, cvs3 * cv3 / cavp)
            };

            let fm = cc * cavm + 1.0;
            let fp = cc * cavp + 1.0;

            let rez = dsp * fp + dsm * fm;

            // Lower, main, and upper diagonals
            w1_n[i] = -fm + cc * dsm * cavm_s1;
            w2_n[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            w3_n[i] = -fp + cc * dsp * cavp_s3;

            // Residual
            w4[i] = -rez;

            cv2 = cv3;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        // Fix endpoints at TE
        w2_n[0] = 1.0;
        w3_n[0] = 0.0;
        w4[0] = 0.0;
        w1_n[nn - 1] = 0.0;
        w2_n[nn - 1] = 1.0;
        w4[nn - 1] = 0.0;

        // Fudge equations adjacent to TE to get TE panel length ratio RTF
        if rtf != 1.0 {
            let i = 1;
            w4[i] = -((snew[i] - snew[i - 1]) + rtf * (snew[i] - snew[i + 1]));
            w1_n[i] = -1.0;
            w2_n[i] = 1.0 + rtf;
            w3_n[i] = -rtf;

            let i = nn - 2;
            w4[i] = -((snew[i] - snew[i + 1]) + rtf * (snew[i] - snew[i - 1]));
            w3_n[i] = -1.0;
            w2_n[i] = 1.0 + rtf;
            w1_n[i] = -rtf;
        }

        // Solve for changes in node positions
        w4 = trisol_curvature(&w1_n, &w2_n, &w3_n, &w4);

        // Find under-relaxation factor
        let mut rlx = 1.0;
        let mut dmax = 0.0_f64;
        for i in 0..nn - 1 {
            let ds = snew[i + 1] - snew[i];
            let dds = w4[i + 1] - w4[i];
            let dsrat = 1.0 + rlx * dds / ds;
            if dsrat > 4.0 {
                rlx = (4.0 - 1.0) * ds / dds;
            }
            if dsrat < 0.2 {
                rlx = (0.2 - 1.0) * ds / dds;
            }
            dmax = dmax.max(w4[i].abs());
        }

        // Update node positions
        for i in 1..nn - 1 {
            snew[i] += rlx * w4[i];
        }

        if dmax.abs() < 1e-3 {
            break;
        }
    }

    // Set new panel node coordinates by sampling every IPFAC-th point
    let mut x_c = Vec::with_capacity(n_panels);
    let mut y_c = Vec::with_capacity(n_panels);

    for i in 0..n_panels {
        let ind = ipfac * i;
        x_c.push(seval(snew[ind], &geometry.x_c, &xbp, &sb));
        y_c.push(seval(snew[ind], &geometry.y_c, &ybp, &sb));
    }

    Geometry {
        reference: geometry.reference,
        x_c,
        y_c,
    }
}

/// Compute curvature of splined curve at arc length parameter ss
///
/// Curvature κ = (x'*y'' - y'*x'') / |r'|³
fn curvature(ss: f64, x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64]) -> f64 {
    let xd = deval(ss, x, xp, s);
    let yd = deval(ss, y, yp, s);
    let xdd = d2val(ss, x, xp, s);
    let ydd = d2val(ss, y, yp, s);

    let sd = (xd.powi(2) + yd.powi(2)).sqrt();
    let sd = sd.max(0.001 * (s[s.len() - 1] - s[0]) / s.len() as f64);

    (xd * ydd - yd * xdd) / sd.powi(3)
}

/// Find leading edge arc length parameter using XFOIL's LEFIND algorithm
fn find_le_arc_length(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();

    // TE coordinates
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);

    // Find first guess for SLE by locating where dot product changes sign
    let mut i_le = n / 2;
    for i in 2..n - 2 {
        let xi = seval(s[i], x, xp, s);
        let yi = seval(s[i], y, yp, s);
        let dxte = xi - xte;
        let dyte = yi - yte;

        let xip = seval(s[i + 1], x, xp, s);
        let yip = seval(s[i + 1], y, yp, s);
        let dx = xip - xi;
        let dy = yip - yi;

        let dotp = dxte * dx + dyte * dy;
        if dotp < 0.0 {
            i_le = i;
            break;
        }
    }

    let mut sle = s[i_le];
    let dseps = (s[n - 1] - s[0]) * 1e-5;

    // Newton iteration
    for _ in 0..50 {
        let xle = seval(sle, x, xp, s);
        let yle = seval(sle, y, yp, s);
        let dxds = deval(sle, x, xp, s);
        let dyds = deval(sle, y, yp, s);
        let dxdd = d2val(sle, x, xp, s);
        let dydd = d2val(sle, y, yp, s);

        let xchord = xle - xte;
        let ychord = yle - yte;

        let res = xchord * dxds + ychord * dyds;
        let ress = dxds.powi(2) + dyds.powi(2) + xchord * dxdd + ychord * dydd;

        let mut dsle = -res / ress;
        let dsle_limit = 0.02 * (xchord + ychord).abs();
        dsle = dsle.clamp(-dsle_limit, dsle_limit);

        sle += dsle;

        if dsle.abs() < dseps {
            break;
        }
    }

    sle
}

/// Solve tri-diagonal system for curvature smoothing
fn trisol_curvature(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut cp = vec![0.0; n];
    let mut dp = vec![0.0; n];
    let mut x = vec![0.0; n];

    // Forward elimination
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for i in 1..n {
        let m = b[i] - a[i] * cp[i - 1];
        if m.abs() < 1e-20 {
            cp[i] = 0.0;
            dp[i] = 0.0;
        } else {
            cp[i] = c[i] / m;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
        }
    }

    // Back substitution
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    x
}

/// Repanel an airfoil with a new number of panels using modified cosine spacing
///
/// Redistributes points along the airfoil surface with higher density
/// near the leading edge. This is a simpler alternative to the XFOIL PANE
/// algorithm (use [`repanel_xfoil`] for exact XFOIL matching).
///
/// # Arguments
/// * `geometry` - Input geometry
/// * `n_panels` - Target number of panels
/// * `te_le_ratio` - Ratio of TE panel density to LE panel density (XFOIL default: 0.15)
///                   Values < 1.0 mean coarser panels at TE, finer at LE
///                   Value of 1.0 gives symmetric cosine spacing
///
/// # Returns
/// New geometry with redistributed points
pub fn repanel_cosine(geometry: &Geometry, n_panels: usize, te_le_ratio: f64) -> Geometry {
    let n = geometry.x_c.len();

    // Calculate arc length along the surface
    let s = calculate_arc_length(&geometry.x_c, &geometry.y_c);

    // Create splines for x and y
    let xp = spline(&geometry.x_c, &s);
    let yp = spline(&geometry.y_c, &s);

    let s_total = s[n - 1];

    // Find LE arc length (approximately midway for a closed airfoil)
    let sle = find_le_arc_length(&geometry.x_c, &xp, &geometry.y_c, &yp, &s);

    // Generate new parameter values using modified cosine spacing
    // with different densities at TE vs LE
    let mut s_new = Vec::with_capacity(n_panels + 1);

    // Clamp te_le_ratio to reasonable range
    let ratio = te_le_ratio.clamp(0.05, 2.0);

    // Split panels roughly evenly between upper and lower surfaces
    let n_half = n_panels / 2;

    // Upper surface: from TE (s=0) to LE (s=sle)
    // Use ratio to control panel density at TE relative to LE
    let s_upper_total = sle - s[0];
    for i in 0..=n_half {
        let t = i as f64 / n_half as f64;
        // Modified cosine with ratio parameter
        // Standard cosine: s = 0.5*(1 - cos(pi*t))
        // Modified: stretch/compress the t parameter based on ratio
        let t_mod = if ratio < 1.0 {
            // Coarser at TE (t=0), finer at LE (t=1)
            // Use power function to redistribute
            t.powf(1.0 / (1.0 + ratio))
        } else if ratio > 1.0 {
            // Finer at TE, coarser at LE
            1.0 - (1.0 - t).powf(ratio)
        } else {
            t
        };
        let theta = std::f64::consts::PI * t_mod;
        let s_frac = 0.5 * (1.0 - theta.cos());
        s_new.push(s[0] + s_upper_total * s_frac);
    }

    // Lower surface: from LE (s=sle) to TE (s=s_total)
    let s_lower_total = s_total - sle;
    for i in 1..=(n_panels - n_half) {
        let t = i as f64 / (n_panels - n_half) as f64;
        let t_mod = if ratio < 1.0 {
            // Finer at LE (t=0), coarser at TE (t=1)
            1.0 - (1.0 - t).powf(1.0 / (1.0 + ratio))
        } else if ratio > 1.0 {
            t.powf(ratio)
        } else {
            t
        };
        let theta = std::f64::consts::PI * t_mod;
        let s_frac = 0.5 * (1.0 - theta.cos());
        s_new.push(sle + s_lower_total * s_frac);
    }

    // Evaluate splines at new parameter values
    let x_c: Vec<f64> = s_new.iter().map(|&si| seval(si, &geometry.x_c, &xp, &s)).collect();
    let y_c: Vec<f64> = s_new.iter().map(|&si| seval(si, &geometry.y_c, &yp, &s)).collect();

    Geometry {
        reference: geometry.reference,
        x_c,
        y_c,
    }
}

/// Calculate arc length along the surface
fn calculate_arc_length(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut s = vec![0.0; n];

    for i in 1..n {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }

    s
}

/// Create a PaneledAirfoil from raw geometry
///
/// Computes all derived quantities needed for aerodynamic analysis:
/// - Arc length parameterization
/// - Spline coefficients
/// - Normal vectors
/// - Panel angles
/// - Leading edge location
pub fn create_paneled_airfoil(geometry: &Geometry) -> PaneledAirfoil {
    let n = geometry.x_c.len();
    let x = geometry.x_c.clone();
    let y = geometry.y_c.clone();

    // Calculate arc length
    let s = calculate_arc_length(&x, &y);

    // Create splines
    let xp = spline(&x, &s);
    let yp = spline(&y, &s);

    // Calculate normal vectors and panel angles
    let (nx, ny, apanel) = calculate_normals_and_angles(&x, &y, &xp, &yp, &s);

    // Find leading edge using XFOIL's chord-perpendicular criterion
    let (sle, le_index) = find_leading_edge(&x, &y, &s, &xp, &yp);

    // Calculate chord length
    let chord = calculate_chord(&x, &y);

    // Check for sharp trailing edge
    let te_gap = ((x[0] - x[n - 1]).powi(2) + (y[0] - y[n - 1]).powi(2)).sqrt();
    let sharp_te = te_gap < 0.0001 * chord;

    PaneledAirfoil {
        x,
        y,
        s,
        xp,
        yp,
        nx,
        ny,
        apanel,
        n,
        sle,
        le_index,
        chord,
        sharp_te,
        reference: geometry.reference,
    }
}

/// Calculate normal vectors and panel angles
fn calculate_normals_and_angles(
    x: &[f64],
    y: &[f64],
    xp: &[f64],
    yp: &[f64],
    s: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len();
    let mut nx = vec![0.0; n];
    let mut ny = vec![0.0; n];
    let mut apanel = vec![0.0; n];

    for i in 0..n {
        // Get tangent vector from spline derivatives
        let dx_ds = deval(s[i], x, xp, s);
        let dy_ds = deval(s[i], y, yp, s);

        // Tangent magnitude
        let ds = (dx_ds * dx_ds + dy_ds * dy_ds).sqrt();

        // Unit tangent
        let tx = dx_ds / ds;
        let ty = dy_ds / ds;

        // Normal is perpendicular to tangent (pointing outward for CCW ordering)
        // For TE->upper->LE->lower->TE ordering, normal points outward with this sign
        nx[i] = ty;
        ny[i] = -tx;

        // Panel angle (angle of tangent from horizontal)
        apanel[i] = ty.atan2(tx);
    }

    (nx, ny, apanel)
}

/// Find leading edge arc length parameter and index
///
/// Uses XFOIL's LEFIND algorithm: finds where the surface tangent is perpendicular
/// to the chord line connecting the LE point to the TE.
///
/// The defining condition is: (X-XTE, Y-YTE) · (X', Y') = 0 at S = SLE
///
/// Returns (sle, le_index)
fn find_leading_edge(x: &[f64], y: &[f64], s: &[f64], xp: &[f64], yp: &[f64]) -> (f64, usize) {
    let n = x.len();

    // Convergence tolerance (matches XFOIL)
    let dseps = (s[n - 1] - s[0]) * 1.0e-5;

    // Trailing edge coordinates
    let x_te = 0.5 * (x[0] + x[n - 1]);
    let y_te = 0.5 * (y[0] + y[n - 1]);

    // Get first guess for SLE by finding where dot product changes sign
    // This matches XFOIL's approach exactly
    let mut i_le = n / 2; // fallback
    for i in 2..n - 2 {
        let dxte = x[i] - x_te;
        let dyte = y[i] - y_te;
        let dx = x[i + 1] - x[i];
        let dy = y[i + 1] - y[i];
        let dotp = dxte * dx + dyte * dy;
        if dotp < 0.0 {
            i_le = i;
            break;
        }
    }

    let mut s_le = s[i_le];

    // Check for sharp LE case (doubled point)
    if i_le > 0 && (s[i_le] - s[i_le - 1]).abs() < 1e-14 {
        return (s_le, i_le);
    }

    // Newton iteration to get exact SLE value (matches XFOIL exactly)
    for _ in 0..50 {
        let x_le = seval(s_le, x, xp, s);
        let y_le = seval(s_le, y, yp, s);
        let dxds = deval(s_le, x, xp, s);
        let dyds = deval(s_le, y, yp, s);
        let dxdd = d2val(s_le, x, xp, s);
        let dydd = d2val(s_le, y, yp, s);

        let xchord = x_le - x_te;
        let ychord = y_le - y_te;

        // Drive dot product between chord line and LE tangent to zero
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        // Newton delta for SLE
        let mut dsle = -res / ress;

        // Limit step size (matches XFOIL exactly: ABS(XCHORD+YCHORD), not ABS(XCHORD)+ABS(YCHORD))
        let dsle_limit = 0.02 * (xchord + ychord).abs();
        dsle = dsle.max(-dsle_limit).min(dsle_limit);

        s_le += dsle;

        if dsle.abs() < dseps {
            break;
        }
    }

    (s_le, i_le)
}

/// Calculate chord length (TE to LE distance)
fn calculate_chord(x: &[f64], y: &[f64]) -> f64 {
    // TE is at first/last points
    let x_te = (x[0] + x[x.len() - 1]) / 2.0;
    let y_te = (y[0] + y[y.len() - 1]) / 2.0;

    // LE is at minimum x
    let (i_le, _) = x
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();

    let x_le = x[i_le];
    let y_le = y[i_le];

    ((x_te - x_le).powi(2) + (y_te - y_le).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_arc_length_calculation() {
        // Simple square path
        let x = vec![0.0, 1.0, 1.0, 0.0, 0.0];
        let y = vec![0.0, 0.0, 1.0, 1.0, 0.0];
        let s = calculate_arc_length(&x, &y);

        assert_eq!(s[0], 0.0);
        assert!((s[1] - 1.0).abs() < 1e-10);
        assert!((s[2] - 2.0).abs() < 1e-10);
        assert!((s[3] - 3.0).abs() < 1e-10);
        assert!((s[4] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_arc_length_circle() {
        // Points on a circle should give arc length = angle * radius
        let n = 36;
        let radius = 1.0;
        let x: Vec<f64> = (0..=n)
            .map(|i| radius * (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos())
            .collect();
        let y: Vec<f64> = (0..=n)
            .map(|i| radius * (2.0 * std::f64::consts::PI * i as f64 / n as f64).sin())
            .collect();

        let s = calculate_arc_length(&x, &y);

        // Total arc length should be approximately 2*pi*r
        let total_arc = s[n];
        assert_relative_eq!(total_arc, 2.0 * std::f64::consts::PI * radius, epsilon = 0.1);
    }

    #[test]
    fn test_chord_calculation_symmetric() {
        // Simple symmetric airfoil shape
        // TE at x=1, LE at x=0
        let x = vec![1.0, 0.8, 0.5, 0.2, 0.0, 0.2, 0.5, 0.8, 1.0];
        let y = vec![0.0, 0.05, 0.08, 0.06, 0.0, -0.06, -0.08, -0.05, 0.0];

        let chord = calculate_chord(&x, &y);
        assert_relative_eq!(chord, 1.0, epsilon = 0.01);
    }

    #[test]
    fn test_paneled_airfoil_from_naca0012() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // Check basic properties - should produce exactly requested panels
        assert_eq!(paneled.n, 100);
        assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

        // Arc length should be monotonically increasing
        for i in 1..paneled.s.len() {
            assert!(paneled.s[i] > paneled.s[i - 1]);
        }

        // Total arc length should be roughly 2x chord for thin airfoil
        let total_arc = paneled.s[paneled.n - 1];
        assert!(total_arc > 1.8 && total_arc < 2.5);

        // Leading edge should be approximately at the midpoint of arc length
        assert!(paneled.sle > total_arc * 0.3 && paneled.sle < total_arc * 0.7);
    }

    #[test]
    fn test_normal_vectors_unit_length() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // All normal vectors should have unit length
        for i in 0..paneled.n {
            let mag = (paneled.nx[i].powi(2) + paneled.ny[i].powi(2)).sqrt();
            assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_normal_vectors_point_outward() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // For a symmetric airfoil centered on y=0:
        // - Upper surface (y > 0) normals should have ny > 0
        // - Lower surface (y < 0) normals should have ny < 0
        // (with some tolerance near LE/TE where y is close to 0)

        let mut upper_count = 0;
        let mut lower_count = 0;

        for i in 0..paneled.n {
            if paneled.y[i] > 0.02 {
                // Upper surface
                assert!(
                    paneled.ny[i] > 0.0,
                    "Upper surface normal should point up at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.ny[i]
                );
                upper_count += 1;
            } else if paneled.y[i] < -0.02 {
                // Lower surface
                assert!(
                    paneled.ny[i] < 0.0,
                    "Lower surface normal should point down at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.ny[i]
                );
                lower_count += 1;
            }
        }

        // Should have checked at least some points on each surface
        assert!(upper_count > 10);
        assert!(lower_count > 10);
    }

    #[test]
    fn test_sharp_trailing_edge_detection() {
        use crate::geometry::naca::naca_4digit;

        // NACA 0012 with XFOIL-compatible blunt TE should be detected as blunt
        let geom = naca_4digit("0012", 100).unwrap();
        let paneled = create_paneled_airfoil(&geom);

        // The NACA generator uses original coefficients for blunt TE (XFOIL-compatible)
        // Gap = 0.00252 which is > 0.0001 * chord, so not sharp
        assert!(!paneled.sharp_te);
    }

    #[test]
    fn test_repanel_preserves_shape() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100).unwrap();
        let repaneled = repanel_cosine(&original, 150, 0.15);

        // Should have approximately the target number of points
        assert!(repaneled.x_c.len() > 140 && repaneled.x_c.len() < 160);

        // Extents should be preserved
        let orig_max_x = original.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_x = repaneled.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_x, new_max_x, epsilon = 0.01);

        let orig_min_x = original.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
        let new_min_x = repaneled.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
        assert_relative_eq!(orig_min_x, new_min_x, epsilon = 0.01);

        // Maximum thickness should be preserved
        let orig_max_y = original.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_y = repaneled.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_y, new_max_y, epsilon = 0.01);
    }

    #[test]
    fn test_repanel_cosine_spacing() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100).unwrap();
        let repaneled = repanel_cosine(&original, 100, 0.15);

        // Cosine spacing should cluster points near LE and TE
        // Points near x=0 and x=1 should be closer together than at mid-chord

        // Find spacing near LE (x close to 0)
        let mut le_spacings = Vec::new();
        let mut mid_spacings = Vec::new();

        for i in 1..repaneled.x_c.len() {
            let dx = (repaneled.x_c[i] - repaneled.x_c[i - 1]).abs();
            let dy = (repaneled.y_c[i] - repaneled.y_c[i - 1]).abs();
            let ds = (dx * dx + dy * dy).sqrt();

            let avg_x = (repaneled.x_c[i] + repaneled.x_c[i - 1]) / 2.0;

            if avg_x < 0.1 || avg_x > 0.9 {
                le_spacings.push(ds);
            } else if avg_x > 0.4 && avg_x < 0.6 {
                mid_spacings.push(ds);
            }
        }

        if !le_spacings.is_empty() && !mid_spacings.is_empty() {
            let avg_le: f64 = le_spacings.iter().sum::<f64>() / le_spacings.len() as f64;
            let avg_mid: f64 = mid_spacings.iter().sum::<f64>() / mid_spacings.len() as f64;

            // LE/TE spacing should be smaller than mid-chord spacing
            assert!(
                avg_le < avg_mid,
                "LE spacing {} should be smaller than mid spacing {}",
                avg_le,
                avg_mid
            );
        }
    }
}
