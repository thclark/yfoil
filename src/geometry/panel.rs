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

/// PANGEN (xfoil.f), line for line: the curvature-based panel distribution XFOIL generates
/// from the buffer airfoil (`PANE`, and the `NACA` command). `n_panels` is NPAN;
/// `config` carries CVPAR/CTERAT/CTRRAT/XSREF/XPREF. Includes the sharp-LE (IBLE) and
/// corner (doubled-point) paths. The returned geometry is the N node coordinates; SCALC,
/// SEGSPL, LEFIND, TECALC, NCALC and APCALC then run in `create_paneled_airfoil` exactly as
/// PANGEN's tail does.
pub fn repanel_xfoil(geometry: &Geometry, n_panels: usize, config: &PaneConfig) -> Geometry {
    let nb = geometry.x_c.len();
    if nb < 2 {
        return geometry.clone();
    }
    let xb = &geometry.x_c;
    let yb = &geometry.y_c;
    let (xsref1, xsref2) = config.xsref.unwrap_or((1.0, 1.0));
    let (xpref1, xpref2) = config.xpref.unwrap_or((1.0, 1.0));

    // Number of temporary nodes for panel distribution calculation exceeds the specified
    // panel number by factor of IPFAC.
    let ipfac = 5;
    // number of airfoil panel points
    let mut n = n_panels;

    // set arc length spline parameter; spline raw airfoil coordinates
    let sb = scalc(xb, yb);
    let xbp = segspl(xb, &sb);
    let ybp = segspl(yb, &sb);

    // normalizing length (~ chord)
    let sbref = 0.5 * (sb[nb - 1] - sb[0]);

    // set up curvature array
    let mut w5: Vec<f64> = (0..nb)
        .map(|i| curv(sb[i], xb, &xbp, yb, &ybp, &sb).abs() * sbref)
        .collect();

    // locate LE point arc length value and the normalized curvature there
    let sble = lefind(xb, &xbp, yb, &ybp, &sb);
    let cvle = curv(sble, xb, &xbp, yb, &ybp, &sb).abs() * sbref;

    // check for doubled point (sharp corner) at LE; IBLE is 1-based like the Fortran (0 = none)
    let mut ible = 0usize;
    for i in 1..nb {
        if sble == sb[i - 1] && sble == sb[i] {
            ible = i;
            // 'Sharp leading edge'
            break;
        }
    }

    // set LE, TE points
    let xble = seval(sble, xb, &xbp, &sb);
    let yble = seval(sble, yb, &ybp, &sb);
    let xbte = 0.5 * (xb[0] + xb[nb - 1]);
    let ybte = 0.5 * (yb[0] + yb[nb - 1]);
    let chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

    // set average curvature over 2*NK+1 points within Rcurv of LE point
    let nk: i32 = 3;
    let mut cvsum = 0.0;
    for k in -nk..=nk {
        let frac = k as f64 / nk as f64;
        let sbk = sble + frac * sbref / cvle.max(20.0);
        let cvk = curv(sbk, xb, &xbp, yb, &ybp, &sb).abs() * sbref;
        cvsum += cvk;
    }
    let mut cvavg = cvsum / (2 * nk + 1) as f64;

    // dummy curvature for sharp LE
    if ible != 0 {
        cvavg = 10.0;
    }

    // set curvature attraction coefficient actually used
    let cc = 6.0 * config.cvpar;

    // set artificial curvature at TE to bunch panels there
    let cvte = cvavg * config.cterat;
    w5[0] = cvte;
    w5[nb - 1] = cvte;

    // set smoothing length = 1 / averaged LE curvature, but no more than 5% of chord and no
    // less than 1/4 average panel spacing
    let smool = (1.0 / cvavg.max(20.0)).max(0.25 / ((n_panels / 2) as f64));
    let smoosq = (smool * sbref) * (smool * sbref);

    // set up tri-diagonal system for smoothed curvatures
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
            // leave curvature at corner point unchanged
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

    // fix curvature at LE point by modifying equations adjacent to LE
    for i in 1..nb - 1 {
        // (I is 1-based in the Fortran: I = i + 1)
        if sb[i] == sble || i + 1 == ible || i + 1 == ible + 1 {
            // if node falls right on LE point, fix curvature there
            w1[i] = 0.0;
            w2[i] = 1.0;
            w3[i] = 0.0;
            w5[i] = cvle;
        } else if sb[i - 1] < sble && sb[i] > sble {
            // modify equation at node just before LE point
            let dsm = sb[i - 1] - sb[i - 2];
            let dsp = sble - sb[i - 1];
            let dso = 0.5 * (sble - sb[i - 2]);
            w1[i - 1] = smoosq * (-1.0 / dsm) / dso;
            w2[i - 1] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3[i - 1] = 0.0;
            w5[i - 1] += smoosq * cvle / (dsp * dso);

            // modify equation at node just after LE point
            let dsm = sb[i] - sble;
            let dsp = sb[i + 1] - sb[i];
            let dso = 0.5 * (sb[i + 1] - sble);
            w1[i] = 0.0;
            w2[i] = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3[i] = smoosq * (-1.0 / dsp) / dso;
            w5[i] += smoosq * cvle / (dsm * dso);
            break;
        }
    }

    // set artificial curvature at bunching points and fix it there
    for i in 1..nb - 1 {
        // chord-based x/c coordinate
        let xoc = ((xb[i] - xble) * (xbte - xble) + (yb[i] - yble) * (ybte - yble)) / chbsq;
        if sb[i] < sble {
            // check if top side point is in refinement area
            if xoc > xsref1 && xoc < xsref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * config.ctrrat;
            }
        } else {
            // check if bottom side point is in refinement area
            if xoc > xpref1 && xoc < xpref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * config.ctrrat;
            }
        }
    }

    // solve for smoothed curvature array W5
    if ible == 0 {
        trisol(&mut w2, &w1, &mut w3, &mut w5);
    } else {
        let i = ible;
        trisol(&mut w2[..i], &w1[..i], &mut w3[..i], &mut w5[..i]);
        trisol(&mut w2[i..], &w1[i..], &mut w3[i..], &mut w5[i..]);
    }

    // find max curvature; normalize curvature array
    let mut cvmax = 0.0_f64;
    for v in &w5 {
        cvmax = cvmax.max(v.abs());
    }
    for v in &mut w5 {
        *v /= cvmax;
    }

    // spline curvature array
    let w6 = segspl(&w5, &sb);

    // Set initial guess for node positions uniform in s. More nodes than specified (by
    // factor of IPFAC) are temporarily used for more reliable convergence.
    let nn = ipfac * (n - 1) + 1;

    // ratio of lengths of panel at TE to one away from the TE
    let rdste = 0.667;
    let rtf = (rdste - 1.0) * 2.0 + 1.0;

    let mut snew = vec![0.0; nn];
    let mut nn1 = 0usize;
    if ible == 0 {
        let dsavg = (sb[nb - 1] - sb[0]) / ((nn - 3) as f64 + 2.0 * rtf);
        snew[0] = sb[0];
        for i in 1..nn - 1 {
            snew[i] = sb[0] + dsavg * ((i - 1) as f64 + rtf);
        }
        snew[nn - 1] = sb[nb - 1];
    } else {
        let nfrac1 = (n * ible) / nb;
        nn1 = ipfac * (nfrac1 - 1) + 1;
        let dsavg1 = (sble - sb[0]) / ((nn1 - 2) as f64 + rtf);
        snew[0] = sb[0];
        for i in 1..nn1 {
            snew[i] = sb[0] + dsavg1 * ((i - 1) as f64 + rtf);
        }
        let nn2 = nn - nn1 + 1;
        let dsavg2 = (sb[nb - 1] - sble) / ((nn2 - 2) as f64 + rtf);
        for i in 1..nn2 - 1 {
            snew[i - 1 + nn1] = sble + dsavg2 * ((i - 1) as f64 + rtf);
        }
        snew[nn - 1] = sb[nb - 1];
    }

    // Newton iteration loop for new node positions
    let mut w1n = vec![0.0; nn];
    let mut w2n = vec![0.0; nn];
    let mut w3n = vec![0.0; nn];
    let mut w4 = vec![0.0; nn];
    for _iter in 1..=20 {
        // set up tri-diagonal system for node position deltas
        let mut cv2 = seval(snew[1], &w5, &w6, &sb);
        let mut cvs2 = deval(snew[1], &w5, &w6, &sb);
        let cv1 = seval(snew[0], &w5, &w6, &sb);
        let cvs1 = deval(snew[0], &w5, &w6, &sb);
        let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
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
            let cavp = (cv3 * cv3 + cv2 * cv2).sqrt();
            let (cavp_s2, cavp_s3) = if cavp == 0.0 {
                (0.0, 0.0)
            } else {
                (cvs2 * cv2 / cavp, cvs3 * cv3 / cavp)
            };
            let fm = cc * cavm + 1.0;
            let fp = cc * cavp + 1.0;
            let rez = dsp * fp + dsm * fm;
            // lower, main, and upper diagonals
            w1n[i] = -fm + cc * dsm * cavm_s1;
            w2n[i] = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            w3n[i] = -fp + cc * dsp * cavp_s3;
            // residual, requiring that (1 + C*curv)*deltaS is equal on both sides of node i
            w4[i] = -rez;
            cv2 = cv3;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        // fix endpoints (at TE)
        w2n[0] = 1.0;
        w3n[0] = 0.0;
        w4[0] = 0.0;
        w1n[nn - 1] = 0.0;
        w2n[nn - 1] = 1.0;
        w4[nn - 1] = 0.0;

        if rtf != 1.0 {
            // fudge equations adjacent to TE to get TE panel length ratio RTF
            let i = 1;
            w4[i] = -((snew[i] - snew[i - 1]) + rtf * (snew[i] - snew[i + 1]));
            w1n[i] = -1.0;
            w2n[i] = 1.0 + rtf;
            w3n[i] = -rtf;
            let i = nn - 2;
            w4[i] = -((snew[i] - snew[i + 1]) + rtf * (snew[i] - snew[i - 1]));
            w3n[i] = -1.0;
            w2n[i] = 1.0 + rtf;
            w1n[i] = -rtf;
        }

        // fix sharp LE point
        if ible != 0 {
            let i = nn1 - 1;
            w1n[i] = 0.0;
            w2n[i] = 1.0;
            w3n[i] = 0.0;
            w4[i] = sble - snew[i];
        }

        // solve for changes W4 in node position arc length values
        trisol(&mut w2n, &w1n, &mut w3n, &mut w4);

        // find under-relaxation factor to keep nodes from changing order
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
            dmax = w4[i].abs().max(dmax);
        }

        // update node position
        for i in 1..nn - 1 {
            snew[i] += rlx * w4[i];
        }

        if dmax.abs() < 1.0e-3 {
            break;
        }
    }
    // 'Paneling convergence failed.  Continuing anyway...' if the loop ran out

    // set new panel node coordinates
    let mut s = Vec::with_capacity(n + 4);
    let mut x = Vec::with_capacity(n + 4);
    let mut y = Vec::with_capacity(n + 4);
    for i in 0..n {
        let ind = ipfac * i;
        s.push(snew[ind]);
        x.push(seval(snew[ind], xb, &xbp, &sb));
        y.push(seval(snew[ind], yb, &ybp, &sb));
    }

    // go over buffer airfoil again, checking for corners (double points)
    for ib in 0..nb - 1 {
        if sb[ib] == sb[ib + 1] {
            // found one !
            let xbcorn = xb[ib];
            let ybcorn = yb[ib];
            let sbcorn = sb[ib];
            // find current-airfoil panel which contains corner (the node count grows on insertion,
            // so this is a while loop rather than a range)
            let mut i = 0;
            while i < n {
                // keep stepping until first node past corner
                if s[i] <= sbcorn {
                    i += 1;
                    continue;
                }
                // move remainder of panel nodes to make room for additional node
                x.insert(i, xbcorn);
                y.insert(i, ybcorn);
                s.insert(i, sbcorn);
                n += 1;
                // shift nodes adjacent to corner to keep panel sizes comparable
                if i >= 2 {
                    s[i - 1] = 0.5 * (s[i] + s[i - 2]);
                    x[i - 1] = seval(s[i - 1], xb, &xbp, &sb);
                    y[i - 1] = seval(s[i - 1], yb, &ybp, &sb);
                }
                if i + 2 < n {
                    s[i + 1] = 0.5 * (s[i] + s[i + 2]);
                    x[i + 1] = seval(s[i + 1], xb, &xbp, &sb);
                    y[i + 1] = seval(s[i + 1], yb, &ybp, &sb);
                }
                // go on to next input geometry point to check for corner
                break;
            }
        }
    }

    Geometry {
        reference: geometry.reference,
        x_c: x,
        y_c: y,
    }
}

/// SCALC: arc length array of a 2-D point array.
pub fn scalc(x: &[f64], y: &[f64]) -> Vec<f64> {
    let mut s = vec![0.0; x.len()];
    for i in 1..x.len() {
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        s[i] = s[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    s
}

/// SEGSPL: splines X(S) like SPLINE but allows derivative discontinuities at segment joints,
/// defined by identical successive S values.
pub fn segspl(x: &[f64], s: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert!(s[0] != s[1], "SEGSPL:  First input point duplicated");
    assert!(s[n - 1] != s[n - 2], "SEGSPL:  Last  input point duplicated");
    let mut xs = vec![0.0; n];
    let mut iseg0 = 0;
    for iseg in 1..n - 2 {
        if s[iseg] == s[iseg + 1] {
            let seg = spline(&x[iseg0..=iseg], &s[iseg0..=iseg]);
            xs[iseg0..=iseg].copy_from_slice(&seg);
            iseg0 = iseg + 1;
        }
    }
    let seg = spline(&x[iseg0..], &s[iseg0..]);
    xs[iseg0..].copy_from_slice(&seg);
    xs
}

/// CURV: curvature of the splined 2-D curve at S = SS, evaluated from the spline's own cubic.
pub fn curv(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
    let n = s.len();
    let mut ilow = 0usize;
    let mut i = n - 1;
    while i - ilow > 1 {
        let imid = (i + ilow) / 2;
        if ss < s[imid] {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    let ds = s[i] - s[i - 1];
    let t = (ss - s[i - 1]) / ds;
    let cx1 = ds * xs[i - 1] - x[i] + x[i - 1];
    let cx2 = ds * xs[i] - x[i] + x[i - 1];
    let xd = x[i] - x[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;
    let xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    let cy1 = ds * ys[i - 1] - y[i] + y[i - 1];
    let cy2 = ds * ys[i] - y[i] + y[i - 1];
    let yd = y[i] - y[i - 1] + (1.0 - 4.0 * t + 3.0 * t * t) * cy1 + t * (3.0 * t - 2.0) * cy2;
    let ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;
    let mut sd = (xd * xd + yd * yd).sqrt();
    sd = sd.max(0.001 * ds);
    (xd * ydd - yd * xdd) / (sd * sd * sd)
}

/// LEFIND: the leading-edge spline parameter SLE where the surface tangent is normal to the
/// chord line from the TE point.
pub fn lefind(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64]) -> f64 {
    let n = x.len();
    // convergence tolerance
    let dseps = (s[n - 1] - s[0]) * 1.0e-5;
    // set trailing edge point coordinates
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);
    // get first guess for SLE (I = 3..N-2 in the Fortran; the loop variable ends at N-1)
    let mut i = n - 2;
    for ii in 2..n - 2 {
        let dxte = x[ii] - xte;
        let dyte = y[ii] - yte;
        let dx = x[ii + 1] - x[ii];
        let dy = y[ii + 1] - y[ii];
        let dotp = dxte * dx + dyte * dy;
        if dotp < 0.0 {
            i = ii;
            break;
        }
    }
    let mut sle = s[i];
    // check for sharp LE case
    if s[i] == s[i - 1] {
        return sle;
    }
    // Newton iteration to get exact SLE value
    for _iter in 1..=50 {
        let xle = seval(sle, x, xp, s);
        let yle = seval(sle, y, yp, s);
        let dxds = deval(sle, x, xp, s);
        let dyds = deval(sle, y, yp, s);
        let dxdd = d2val(sle, x, xp, s);
        let dydd = d2val(sle, y, yp, s);
        let xchord = xle - xte;
        let ychord = yle - yte;
        let res = xchord * dxds + ychord * dyds;
        let ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;
        let mut dsle = -res / ress;
        dsle = dsle.max(-0.02 * (xchord + ychord).abs());
        dsle = dsle.min(0.02 * (xchord + ychord).abs());
        sle += dsle;
        if dsle.abs() < dseps {
            return sle;
        }
    }
    // 'LEFIND:  LE point not found.  Continuing...'
    s[i]
}

/// TRISOL: solves the tri-diagonal system with main diagonal `a`, lower `b`, upper `c` and
/// right-hand side `d`; `d` is replaced by the solution, `a` and `c` are destroyed.
pub fn trisol(a: &mut [f64], b: &[f64], c: &mut [f64], d: &mut [f64]) {
    let kk = a.len();
    for k in 1..kk {
        let km = k - 1;
        c[km] /= a[km];
        d[km] /= a[km];
        a[k] -= b[k] * c[km];
        d[k] -= b[k] * d[km];
    }
    d[kk - 1] /= a[kk - 1];
    for k in (0..kk - 1).rev() {
        d[k] -= c[k] * d[k + 1];
    }
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
    let sle = lefind(&geometry.x_c, &xp, &geometry.y_c, &yp, &s);

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

    // SCALC / SEGSPL
    let s = calculate_arc_length(&x, &y);
    let xp = spline(&x, &s);
    let yp = spline(&y, &s);

    // NCALC: node normals from the spline derivative arrays
    let (nx, ny) = ncalc(&xp, &yp, &s);

    // LEFIND / GEOPAR: leading edge on the spline; chord is the LE–TE distance (XFOIL's
    // definition — for a NACA section whose nodes straddle the LE this is slightly under 1)
    let (sle, le_index) = find_leading_edge(&x, &y, &s, &xp, &yp);
    let xle = seval(sle, &x, &xp, &s);
    let yle = seval(sle, &y, &yp, &s);
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();

    // TECALC: SHARP = DSTE < 0.0001*CHORD
    let dste = ((x[0] - x[n - 1]).powi(2) + (y[0] - y[n - 1]).powi(2)).sqrt();
    let sharp_te = dste < 0.0001 * chord;

    // APCALC: panel angles (needs SHARP for the TE panel)
    let apanel = apcalc(&x, &y, &nx, &ny, sharp_te);

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

/// NCALC (xpanel.f): unit normal vector components at airfoil panel nodes, from the spline
/// derivative arrays (SEGSPL output), with corner-point averaging where S(I) == S(I+1).
fn ncalc(xp: &[f64], yp: &[f64], s: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = xp.len();
    let mut xn = vec![0.0; n];
    let mut yn = vec![0.0; n];
    if n <= 1 {
        return (xn, yn);
    }
    for i in 0..n {
        let sx = yp[i];
        let sy = -xp[i];
        let smod = (sx * sx + sy * sy).sqrt();
        if smod == 0.0 {
            xn[i] = -1.0;
            yn[i] = 0.0;
        } else {
            xn[i] = sx / smod;
            yn[i] = sy / smod;
        }
    }
    // average normal vectors at corner points
    for i in 0..n - 1 {
        if s[i] == s[i + 1] {
            let sx = 0.5 * (xn[i] + xn[i + 1]);
            let sy = 0.5 * (yn[i] + yn[i + 1]);
            let smod = (sx * sx + sy * sy).sqrt();
            if smod == 0.0 {
                xn[i] = -1.0;
                yn[i] = 0.0;
                xn[i + 1] = -1.0;
                yn[i + 1] = 0.0;
            } else {
                xn[i] = sx / smod;
                yn[i] = sy / smod;
                xn[i + 1] = sx / smod;
                yn[i + 1] = sy / smod;
            }
        }
    }
    (xn, yn)
}

/// APCALC (xpanel.f): angle of each airfoil panel (panel `i` runs from node `i` to `i+1`;
/// the TE panel `n-1` closes from node `n-1` back to node 0).
fn apcalc(x: &[f64], y: &[f64], nx: &[f64], ny: &[f64], sharp: bool) -> Vec<f64> {
    let n = x.len();
    let pi = 4.0 * (1.0_f64).atan();
    let mut apanel = vec![0.0; n];
    for i in 0..n - 1 {
        let sx = x[i + 1] - x[i];
        let sy = y[i + 1] - y[i];
        apanel[i] = if sx == 0.0 && sy == 0.0 {
            (-ny[i]).atan2(-nx[i])
        } else {
            sx.atan2(-sy)
        };
    }
    // TE panel
    let (i, ip) = (n - 1, 0);
    apanel[i] = if sharp {
        pi
    } else {
        let sx = x[ip] - x[i];
        let sy = y[ip] - y[i];
        (-sx).atan2(sy) + pi
    };
    apanel
}

/// LEFIND plus the index of the first node past the LE (kept for callers that want it).
fn find_leading_edge(x: &[f64], y: &[f64], s: &[f64], xp: &[f64], yp: &[f64]) -> (f64, usize) {
    let sle = lefind(x, xp, y, yp, s);
    // YFoil convenience only (XFOIL works with SLE): the node nearest the spline LE
    let i_le = (0..x.len())
        .min_by(|&i, &j| (s[i] - sle).abs().partial_cmp(&(s[j] - sle).abs()).unwrap())
        .unwrap_or(0);
    (sle, i_le)
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
