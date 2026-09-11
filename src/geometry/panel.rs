//! Panel distribution and repaneling
//!
//! Functions for redistributing panel points on an airfoil surface.

use super::airfoil::{Geometry, InvalidGeometryError, PanelledFoil};
use super::spline::{spline_derivatives, spline_second_derivative, spline_slope, spline_value};
use serde::{Deserialize, Serialize};

/// XFOIL's `PPAR` parameters of PANGEN (menu keys `P`, `T`, `R`, `XT`, `XB`). The JSON keys of the
/// provenance record and of a `--panelling` file are these field names; the XFOIL variable of
/// each is in its doc comment and in `docs/xfoil-reference.md`.
#[doc(alias = "PPAR")]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct PangenConfig {
    /// Curvature bunching parameter (XFOIL CVPAR, menu `P`): the curvature attraction
    /// coefficient is 6 × this; 0 gives uniform arc-length spacing. Default 1.0
    pub curvature_bunching: f64,
    /// Fictitious trailing-edge curvature as a fraction of the averaged leading-edge curvature
    /// (XFOIL CTERAT, menu `T`, "TE/LE panel density ratio"): bunches nodes at the trailing
    /// edge. Default 0.15
    pub te_curvature_ratio: f64,
    /// Fictitious curvature inside the refinement windows as a fraction of the leading-edge
    /// curvature (XFOIL CTRRAT, menu `R`). Default 0.2
    pub refined_curvature_ratio: f64,
    /// Upper-surface refinement window in x/c (XFOIL XSREF1, XSREF2, menu `XT`);
    /// `None` is XFOIL's `1.0 1.0`, refinement off
    pub refine_upper: Option<[f64; 2]>,
    /// Lower-surface refinement window in x/c (XFOIL XPREF1, XPREF2, menu `XB`)
    pub refine_lower: Option<[f64; 2]>,
}

impl Default for PangenConfig {
    fn default() -> Self {
        Self {
            curvature_bunching: 1.0,
            te_curvature_ratio: 0.15,
            refined_curvature_ratio: 0.2,
            refine_upper: None,
            refine_lower: None,
        }
    }
}

/// yFoil's own arc-length cosine spacing (`repanel_cosine`), no XFOIL equivalent. `te_bias`
/// warps the cosine parameter: 1 is a plain cosine, below 1 coarser at the trailing edge and
/// finer at the leading edge, above 1 finer at the trailing edge; clamped to 0.05…2. The method
/// writes n + 1 nodes (its historic behaviour, frozen by `tests/repanel_cosine_tests.rs`).
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct CosineConfig {
    /// `None` in a request means the default 0.15 when repanelling existing nodes; a generator's
    /// cosine sampling has no bias and records `null`
    pub te_bias: Option<f64>,
}

impl CosineConfig {
    pub const DEFAULT_TE_BIAS: f64 = 0.15;
}

/// The node-distribution method and its parameters
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "method", rename_all = "lowercase")]
pub enum PanelMethod {
    /// XFOIL's PANGEN (`PANE` / `PPAR`), the default
    Pangen(PangenConfig),
    /// yFoil's own cosine spacing
    Cosine(CosineConfig),
}

/// The trailing-edge gap set with XFOIL's `TGAP` after the nodes are distributed: `gap` in
/// chord units, `blend` the blending distance/c (TGAP's second argument, default 1.0)
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TeGap {
    pub gap: f64,
    #[serde(default = "TeGap::default_blend")]
    pub blend: f64,
}

impl TeGap {
    fn default_blend() -> f64 {
        1.0
    }
}

/// How a geometry is panelled: the complete recipe from a section definition or a loaded set of
/// nodes to the output nodes, and the one structure the CLI flags, a `--panelling` JSON file and
/// the `generator.panelling` provenance record share. It serialises flat:
///
/// ```json
/// {"method": "pangen", "n_nodes": 160, "sharp_te": false, "te_gap": null,
///  "curvature_bunching": 1.0, "te_curvature_ratio": 0.15, "refined_curvature_ratio": 0.2,
///  "refine_upper": null, "refine_lower": null}
/// {"method": "cosine", "n_nodes": 160, "sharp_te": false, "te_gap": null, "te_bias": 0.15}
/// ```
///
/// The trailing-edge treatment (`sharp_te`, or `te_gap`; they are exclusive) is applied to the
/// panelled nodes, after the distribution, so that the TGAP blend profile
/// ½Δ·(x/c)·exp(−(1 − x/c)(1/blend − 1)) is evaluated exactly at every output node rather than
/// splined through a possibly coarse input. XFOIL applies TGAP to the buffer airfoil and PANE
/// follows; the two orders differ by PANGEN's spline interpolation of the displaced buffer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PanelConfig {
    /// Number of panel nodes (XFOIL NPAN)
    pub n_nodes: usize,
    /// Move the two trailing-edge nodes to their midpoint (a closed edge; XFOIL's SHARP path)
    #[serde(default)]
    pub sharp_te: bool,
    /// Set the trailing-edge gap with TGAP; `None` leaves the edge as it is
    #[serde(default)]
    pub te_gap: Option<TeGap>,
    #[serde(flatten)]
    pub method: PanelMethod,
    /// Generators only: the number of nodes at which the analytic section is sampled before
    /// PANGEN (default [`PANGEN_BUFFER_NODES`]). Rejected by [`repanel`] of an existing geometry.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub buffer_nodes: Option<usize>,
}

impl Default for PanelConfig {
    /// PANGEN with XFOIL's defaults, 160 nodes, no trailing-edge treatment
    fn default() -> Self {
        Self {
            n_nodes: 160,
            sharp_te: false,
            te_gap: None,
            method: PanelMethod::Pangen(PangenConfig::default()),
            buffer_nodes: None,
        }
    }
}

/// Nodes at which a generator samples its analytic section before PANGEN: 123 per side, the
/// density of XFOIL's own NACA buffer (`XFOIL_NACA_NSIDE`), leading edge straddled
pub const PANGEN_BUFFER_NODES: usize = 246;

/// The keys a panelling record may carry, per method (`method`, `n_nodes`, `sharp_te`, `te_gap`,
/// `buffer_nodes` are common)
const PANEL_KEYS_COMMON: [&str; 5] = ["method", "n_nodes", "sharp_te", "te_gap", "buffer_nodes"];
const PANEL_KEYS_PANGEN: [&str; 5] = [
    "curvature_bunching",
    "te_curvature_ratio",
    "refined_curvature_ratio",
    "refine_upper",
    "refine_lower",
];
const PANEL_KEYS_COSINE: [&str; 1] = ["te_bias"];

#[derive(thiserror::Error, Debug)]
pub enum RepanelError {
    #[error("a sharp trailing edge (sharp_te) and a trailing-edge gap (te_gap) are contradictory")]
    SharpWithGap,
    #[error("refinement window {0}: two x/c values are needed, the first smaller than the second")]
    RefinementWindow(&'static str),
    #[error("the {key} key belongs to the {method} method, not {given}")]
    KeyOfOtherMethod {
        key: String,
        method: &'static str,
        given: &'static str,
    },
    #[error("unknown panelling key {0}")]
    UnknownKey(String),
    #[error("panelling has no method: {0}")]
    Malformed(String),
    #[error("buffer_nodes applies to generated sections only, not to repanelling an existing geometry")]
    BufferNodesOnRepanel,
    #[error(transparent)]
    InvalidGeometry(#[from] InvalidGeometryError),
}

impl PanelConfig {
    /// The method's name as the record writes it
    pub fn method_name(&self) -> &'static str {
        match self.method {
            PanelMethod::Pangen(_) => "pangen",
            PanelMethod::Cosine(_) => "cosine",
        }
    }

    /// The conflicts a set of flags or a file can express: a closed edge with a gap, and a
    /// refinement window that is not two increasing values
    pub fn validate(&self) -> Result<(), RepanelError> {
        if self.sharp_te && self.te_gap.is_some() {
            return Err(RepanelError::SharpWithGap);
        }
        if let PanelMethod::Pangen(p) = &self.method {
            for (name, window) in [("refine_upper", p.refine_upper), ("refine_lower", p.refine_lower)] {
                if let Some([a, b]) = window {
                    // an error unless a < b (so NaN fails too)
                    if a.partial_cmp(&b) != Some(std::cmp::Ordering::Less) {
                        return Err(RepanelError::RefinementWindow(name));
                    }
                }
            }
        }
        Ok(())
    }

    /// Read a record (the shape [`PanelConfig`] serialises to), rejecting keys that belong to
    /// the other method and unknown keys by name — serde's `deny_unknown_fields` cannot do this
    /// through the flattened method — then validating the combination
    pub fn from_json(value: &serde_json::Value) -> Result<Self, RepanelError> {
        let obj = value
            .as_object()
            .ok_or_else(|| RepanelError::Malformed("not a JSON object".into()))?;
        let method = obj
            .get("method")
            .and_then(|m| m.as_str())
            .ok_or_else(|| RepanelError::Malformed("no \"method\" key (pangen or cosine)".into()))?;
        let (allowed, other, other_name): (&[&str], &[&str], &'static str) = match method {
            "pangen" => (&PANEL_KEYS_PANGEN, &PANEL_KEYS_COSINE, "cosine"),
            "cosine" => (&PANEL_KEYS_COSINE, &PANEL_KEYS_PANGEN, "pangen"),
            m => {
                return Err(RepanelError::Malformed(format!(
                    "method {m:?} is neither pangen nor cosine"
                )))
            }
        };
        let given: &'static str = if method == "pangen" { "pangen" } else { "cosine" };
        for key in obj.keys() {
            if PANEL_KEYS_COMMON.contains(&key.as_str()) || allowed.contains(&key.as_str()) {
                continue;
            }
            if other.contains(&key.as_str()) {
                return Err(RepanelError::KeyOfOtherMethod {
                    key: key.clone(),
                    method: other_name,
                    given,
                });
            }
            return Err(RepanelError::UnknownKey(key.clone()));
        }
        let config: PanelConfig =
            serde_json::from_value(value.clone()).map_err(|e| RepanelError::Malformed(e.to_string()))?;
        config.validate()?;
        Ok(config)
    }
}

/// Apply the trailing-edge treatment of `config` to panelled nodes
pub(crate) fn apply_te_treatment(geometry: Geometry, config: &PanelConfig) -> Geometry {
    if config.sharp_te {
        return geometry.sharpen();
    }
    match config.te_gap {
        Some(TeGap { gap, blend }) => set_te_gap(&geometry, gap, blend),
        None => geometry,
    }
}

/// Write `config` into the geometry's provenance record as `panelling` (creating the record
/// when the geometry has none, as a loaded `.dat` has not)
pub(crate) fn record_panelling(geometry: &mut Geometry, config: &PanelConfig) {
    let panelling = serde_json::to_value(config).expect("PanelConfig serialises");
    match geometry.generator.as_mut() {
        Some(rec) => {
            rec["panelling"] = panelling;
        }
        None => {
            geometry.generator = Some(serde_json::json!({
                "yfoil": env!("CARGO_PKG_VERSION"),
                "panelling": panelling,
            }));
        }
    }
}

/// Redistribute the nodes of an existing geometry with `config` — XFOIL's `PANE`/`PPAR` on a
/// LOADed airfoil — then apply its trailing-edge treatment, validate the result and record the
/// panelling in `generator.panelling`.
#[doc(alias = "PANE")]
#[doc(alias = "PPAR")]
pub fn repanel(geometry: &Geometry, config: &PanelConfig) -> Result<Geometry, RepanelError> {
    config.validate()?;
    if config.buffer_nodes.is_some() {
        return Err(RepanelError::BufferNodesOnRepanel);
    }
    // the record states the bias actually used
    let mut used = *config;
    let distributed = match &mut used.method {
        PanelMethod::Pangen(p) => repanel_by_curvature(geometry, config.n_nodes, p),
        PanelMethod::Cosine(c) => {
            let bias = c.te_bias.unwrap_or(CosineConfig::DEFAULT_TE_BIAS);
            c.te_bias = Some(bias);
            repanel_cosine(geometry, config.n_nodes, bias)
        }
    };
    let mut out = apply_te_treatment(distributed, &used);
    out.validate()?;
    record_panelling(&mut out, &used);
    Ok(out)
}

/// PANGEN (xfoil.f), line for line: the curvature-based panel distribution XFOIL generates
/// from the buffer airfoil (`PANE`, and the `NACA` command). `n_panels` is NPAN;
/// `config` carries CVPAR/CTERAT/CTRRAT/XSREF/XPREF. Includes the sharp-LE (IBLE) and
/// corner (doubled-point) paths. The returned geometry is the N node coordinates; SCALC,
/// SEGSPL, LEFIND, TECALC, NCALC and APCALC then run in `panel_foil` exactly as
/// PANGEN's tail does.
#[doc(alias = "PANGEN")]
pub fn repanel_by_curvature(geometry: &Geometry, n_panels: usize, config: &PangenConfig) -> Geometry {
    let nb = geometry.x.len();
    if nb < 2 {
        return geometry.clone();
    }
    let xb = &geometry.x;
    let yb = &geometry.y;
    // XSREF/XPREF: a window of 1.0 1.0 is XFOIL's "refinement off"
    let [xsref1, xsref2] = config.refine_upper.unwrap_or([1.0, 1.0]);
    let [xpref1, xpref2] = config.refine_lower.unwrap_or([1.0, 1.0]);

    // Number of temporary nodes for panel distribution calculation exceeds the specified
    // panel number by factor of IPFAC.
    let ipfac = 5;
    // number of airfoil panel points
    let mut n = n_panels;

    // set arc length spline parameter; spline raw airfoil coordinates
    let sb = arc_coordinate(xb, yb);
    let xbp = spline_segmented(xb, &sb);
    let ybp = spline_segmented(yb, &sb);

    // normalizing length (~ chord)
    let sbref = 0.5 * (sb[nb - 1] - sb[0]);

    // set up curvature array
    let mut w5: Vec<f64> = (0..nb)
        .map(|i| curvature(sb[i], xb, &xbp, yb, &ybp, &sb).abs() * sbref)
        .collect();

    // locate LE point arc length value and the normalized curvature there
    let sble = find_le(xb, &xbp, yb, &ybp, &sb);
    let cvle = curvature(sble, xb, &xbp, yb, &ybp, &sb).abs() * sbref;

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
    let xble = spline_value(sble, xb, &xbp, &sb);
    let yble = spline_value(sble, yb, &ybp, &sb);
    let xbte = 0.5 * (xb[0] + xb[nb - 1]);
    let ybte = 0.5 * (yb[0] + yb[nb - 1]);
    let chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

    // set average curvature over 2*NK+1 points within Rcurv of LE point
    let nk: i32 = 3;
    let mut cvsum = 0.0;
    for k in -nk..=nk {
        let frac = k as f64 / nk as f64;
        let sbk = sble + frac * sbref / cvle.max(20.0);
        let cvk = curvature(sbk, xb, &xbp, yb, &ybp, &sb).abs() * sbref;
        cvsum += cvk;
    }
    let mut cvavg = cvsum / (2 * nk + 1) as f64;

    // dummy curvature for sharp LE
    if ible != 0 {
        cvavg = 10.0;
    }

    // set curvature attraction coefficient actually used
    let cc = 6.0 * config.curvature_bunching;

    // set artificial curvature at TE to bunch panels there
    let cvte = cvavg * config.te_curvature_ratio;
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
                w5[i] = cvle * config.refined_curvature_ratio;
            }
        } else {
            // check if bottom side point is in refinement area
            if xoc > xpref1 && xoc < xpref2 {
                w1[i] = 0.0;
                w2[i] = 1.0;
                w3[i] = 0.0;
                w5[i] = cvle * config.refined_curvature_ratio;
            }
        }
    }

    // solve for smoothed curvature array W5
    if ible == 0 {
        solve_tridiagonal(&mut w2, &w1, &mut w3, &mut w5);
    } else {
        let i = ible;
        solve_tridiagonal(&mut w2[..i], &w1[..i], &mut w3[..i], &mut w5[..i]);
        solve_tridiagonal(&mut w2[i..], &w1[i..], &mut w3[i..], &mut w5[i..]);
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
    let w6 = spline_segmented(&w5, &sb);

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
        let mut cv2 = spline_value(snew[1], &w5, &w6, &sb);
        let mut cvs2 = spline_slope(snew[1], &w5, &w6, &sb);
        let cv1 = spline_value(snew[0], &w5, &w6, &sb);
        let cvs1 = spline_slope(snew[0], &w5, &w6, &sb);
        let mut cavm = (cv1 * cv1 + cv2 * cv2).sqrt();
        let (mut cavm_s1, mut cavm_s2) = if cavm == 0.0 {
            (0.0, 0.0)
        } else {
            (cvs1 * cv1 / cavm, cvs2 * cv2 / cavm)
        };
        for i in 1..nn - 1 {
            let dsm = snew[i] - snew[i - 1];
            let dsp = snew[i] - snew[i + 1];
            let cv3 = spline_value(snew[i + 1], &w5, &w6, &sb);
            let cvs3 = spline_slope(snew[i + 1], &w5, &w6, &sb);
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
        solve_tridiagonal(&mut w2n, &w1n, &mut w3n, &mut w4);

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
        x.push(spline_value(snew[ind], xb, &xbp, &sb));
        y.push(spline_value(snew[ind], yb, &ybp, &sb));
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
                    x[i - 1] = spline_value(s[i - 1], xb, &xbp, &sb);
                    y[i - 1] = spline_value(s[i - 1], yb, &ybp, &sb);
                }
                if i + 2 < n {
                    s[i + 1] = 0.5 * (s[i] + s[i + 2]);
                    x[i + 1] = spline_value(s[i + 1], xb, &xbp, &sb);
                    y[i + 1] = spline_value(s[i + 1], yb, &ybp, &sb);
                }
                // go on to next input geometry point to check for corner
                break;
            }
        }
    }

    Geometry {
        cm_ref: geometry.cm_ref,
        x,
        y,
        generator: geometry.generator.clone(),
    }
}

/// SCALC: arc length array of a 2-D point array.
#[doc(alias = "SCALC")]
pub fn arc_coordinate(x: &[f64], y: &[f64]) -> Vec<f64> {
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
#[doc(alias = "SEGSPL")]
pub fn spline_segmented(x: &[f64], s: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert!(s[0] != s[1], "SEGSPL:  First input point duplicated");
    assert!(s[n - 1] != s[n - 2], "SEGSPL:  Last  input point duplicated");
    let mut xs = vec![0.0; n];
    let mut iseg0 = 0;
    for iseg in 1..n - 2 {
        if s[iseg] == s[iseg + 1] {
            let seg = spline_derivatives(&x[iseg0..=iseg], &s[iseg0..=iseg]);
            xs[iseg0..=iseg].copy_from_slice(&seg);
            iseg0 = iseg + 1;
        }
    }
    let seg = spline_derivatives(&x[iseg0..], &s[iseg0..]);
    xs[iseg0..].copy_from_slice(&seg);
    xs
}

/// CURV: curvature of the splined 2-D curve at S = SS, evaluated from the spline's own cubic.
#[doc(alias = "CURV")]
pub fn curvature(ss: f64, x: &[f64], xs: &[f64], y: &[f64], ys: &[f64], s: &[f64]) -> f64 {
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
#[doc(alias = "LEFIND")]
pub fn find_le(x: &[f64], xp: &[f64], y: &[f64], yp: &[f64], s: &[f64]) -> f64 {
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
        let xle = spline_value(sle, x, xp, s);
        let yle = spline_value(sle, y, yp, s);
        let dxds = spline_slope(sle, x, xp, s);
        let dyds = spline_slope(sle, y, yp, s);
        let dxdd = spline_second_derivative(sle, x, xp, s);
        let dydd = spline_second_derivative(sle, y, yp, s);
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
#[doc(alias = "TRISOL")]
pub fn solve_tridiagonal(a: &mut [f64], b: &[f64], c: &mut [f64], d: &mut [f64]) {
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

/// yFoil's own repanelling, no XFOIL equivalent: cosine spacing in arc length on each surface,
/// the cosine parameter warped by a power law set by `te_bias` (1 is a plain cosine; below 1
/// coarser at the trailing edge and finer at the leading edge; above 1 finer at the trailing
/// edge; clamped to 0.05…2). Writes `n_panels + 1` nodes. The numerics are frozen —
/// `tests/repanel_cosine_tests.rs` gates every node against golden fixtures — because test
/// fixtures were derived with them; new work goes through [`repanel`] with [`PanelMethod::Pangen`].
pub fn repanel_cosine(geometry: &Geometry, n_panels: usize, te_le_ratio: f64) -> Geometry {
    let n = geometry.x.len();

    // Calculate arc length along the surface
    let s = arc_coordinate(&geometry.x, &geometry.y);

    // Create splines for x and y
    let xp = spline_derivatives(&geometry.x, &s);
    let yp = spline_derivatives(&geometry.y, &s);

    let s_total = s[n - 1];

    // Find LE arc length (approximately midway for a closed airfoil)
    let sle = find_le(&geometry.x, &xp, &geometry.y, &yp, &s);

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
    let x_c: Vec<f64> = s_new.iter().map(|&si| spline_value(si, &geometry.x, &xp, &s)).collect();
    let y_c: Vec<f64> = s_new.iter().map(|&si| spline_value(si, &geometry.y, &yp, &s)).collect();

    Geometry {
        cm_ref: geometry.cm_ref,
        x: x_c,
        y: y_c,
        generator: geometry.generator.clone(),
    }
}

/// TGAP (xgdes.f), line for line: set the trailing-edge gap of the buffer airfoil to `gap`
/// (chord units), blending the change into the section over the fraction `blend` of the chord
/// from the trailing edge (XFOIL's "blending distance/c", clamped to 0…1; 0 moves only the two
/// trailing-edge points). The upper and lower surfaces move apart along the unit vector of the
/// existing gap, or along the mean trailing-edge normal when the edge is sharp, by
/// `½ Δgap · (x/c) · exp(−(1 − x/c)(1/blend − 1))` each. XFOIL then re-splines the buffer
/// airfoil; the returned geometry is the moved nodes, with the provenance record kept.
#[doc(alias = "TGAP")]
pub fn set_te_gap(geometry: &Geometry, gap: f64, blend: f64) -> Geometry {
    let nb = geometry.x.len();
    let xb = &geometry.x;
    let yb = &geometry.y;
    let sb = arc_coordinate(xb, yb);
    let xbp = spline_segmented(xb, &sb);
    let ybp = spline_segmented(yb, &sb);

    let sble = find_le(xb, &xbp, yb, &ybp, &sb);
    let xble = spline_value(sble, xb, &xbp, &sb);
    let yble = spline_value(sble, yb, &ybp, &sb);
    let xbte = 0.5 * (xb[0] + xb[nb - 1]);
    let ybte = 0.5 * (yb[0] + yb[nb - 1]);
    let chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

    let dxn = xb[0] - xb[nb - 1];
    let dyn_ = yb[0] - yb[nb - 1];
    let gap_old = (dxn * dxn + dyn_ * dyn_).sqrt();

    // components of unit vector parallel to TE gap
    let (dxu, dyu) = if gap_old > 0.0 {
        (dxn / gap_old, dyn_ / gap_old)
    } else {
        (-0.5 * (ybp[nb - 1] - ybp[0]), 0.5 * (xbp[nb - 1] - xbp[0]))
    };

    let doc = blend.max(0.0).min(1.0);
    let dgap = gap - gap_old;

    let mut x = xb.clone();
    let mut y = yb.clone();
    // go over each point, changing the y-thickness appropriately
    for i in 0..nb {
        // chord-based x/c
        let xoc = ((xb[i] - xble) * (xbte - xble) + (yb[i] - yble) * (ybte - yble)) / chbsq;

        // thickness factor tails off exponentially away from trailing edge
        let tfac = if doc == 0.0 {
            if i == 0 || i == nb - 1 {
                1.0
            } else {
                0.0
            }
        } else {
            let arg = ((1.0 - xoc) * (1.0 / doc - 1.0)).min(15.0);
            (-arg).exp()
        };

        if sb[i] <= sble {
            x[i] = xb[i] + 0.5 * dgap * xoc * tfac * dxu;
            y[i] = yb[i] + 0.5 * dgap * xoc * tfac * dyu;
        } else {
            x[i] = xb[i] - 0.5 * dgap * xoc * tfac * dxu;
            y[i] = yb[i] - 0.5 * dgap * xoc * tfac * dyu;
        }
    }

    Geometry {
        cm_ref: geometry.cm_ref,
        x,
        y,
        generator: geometry.generator.clone(),
    }
}

/// Translates XFOIL's `SCALC`, `SEGSPL`, `LEFIND`, `TECALC`, `NCALC`, `APCALC`.
///
/// Create a PanelledFoil from raw geometry
///
/// Computes all derived quantities needed for aerodynamic analysis:
/// - Arc length parameterization
/// - Spline coefficients
/// - Normal vectors
/// - Panel angles
/// - Leading edge location
#[doc(alias = "SCALC")]
#[doc(alias = "SEGSPL")]
#[doc(alias = "LEFIND")]
#[doc(alias = "TECALC")]
#[doc(alias = "NCALC")]
#[doc(alias = "APCALC")]
pub fn panel_foil(geometry: &Geometry) -> PanelledFoil {
    let n = geometry.x.len();
    let x = geometry.x.clone();
    let y = geometry.y.clone();

    // SCALC / SEGSPL
    let s = arc_coordinate(&x, &y);
    let xp = spline_derivatives(&x, &s);
    let yp = spline_derivatives(&y, &s);

    // NCALC: node normals from the spline derivative arrays
    let (nx, ny) = node_normals(&xp, &yp, &s);

    // LEFIND / GEOPAR: leading edge on the spline; chord is the LE–TE distance (XFOIL's
    // definition — for a NACA section whose nodes straddle the LE this is slightly under 1)
    let sle = find_le(&x, &xp, &y, &yp, &s);
    // yFoil convenience only (XFOIL works with SLE): the node nearest the spline LE
    let le_index = (0..n)
        .min_by(|&i, &j| (s[i] - sle).abs().partial_cmp(&(s[j] - sle).abs()).unwrap())
        .unwrap_or(0);
    let xle = spline_value(sle, &x, &xp, &s);
    let yle = spline_value(sle, &y, &yp, &s);
    let xte = 0.5 * (x[0] + x[n - 1]);
    let yte = 0.5 * (y[0] + y[n - 1]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();

    // TECALC: SHARP = DSTE < 0.0001*CHORD
    let dste = ((x[0] - x[n - 1]).powi(2) + (y[0] - y[n - 1]).powi(2)).sqrt();
    let sharp_te = dste < 0.0001 * chord;

    // APCALC: panel angles (needs SHARP for the TE panel)
    let apanel = panel_angles(&x, &y, &nx, &ny, sharp_te);

    PanelledFoil {
        x,
        y,
        s,
        dxds: xp,
        dyds: yp,
        normal_x: nx,
        normal_y: ny,
        panel_angle: apanel,
        n_foil_nodes: n,
        s_le: sle,
        i_le_node: le_index,
        chord,
        sharp_te,
        cm_ref: geometry.cm_ref,
    }
}

/// NCALC (xpanel.f): unit normal vector components at airfoil panel nodes, from the spline
/// derivative arrays (SEGSPL output), with corner-point averaging where S(I) == S(I+1).
#[doc(alias = "NCALC")]
fn node_normals(xp: &[f64], yp: &[f64], s: &[f64]) -> (Vec<f64>, Vec<f64>) {
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
#[doc(alias = "APCALC")]
fn panel_angles(x: &[f64], y: &[f64], nx: &[f64], ny: &[f64], sharp: bool) -> Vec<f64> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn source() -> Geometry {
        crate::geometry::naca_4digit("2412", 80, crate::geometry::Thickness::Perpendicular).unwrap()
    }

    #[test]
    fn repanel_with_the_default_config_is_pangen_with_xfoil_defaults() {
        let g = source();
        let via_config = repanel(&g, &PanelConfig::default()).unwrap();
        let direct = repanel_by_curvature(&g, 160, &PangenConfig::default());
        assert_eq!(via_config.x.len(), 160);
        for i in 0..160 {
            assert_eq!(via_config.x[i].to_bits(), direct.x[i].to_bits());
            assert_eq!(via_config.y[i].to_bits(), direct.y[i].to_bits());
        }
        let rec = via_config.generator.unwrap();
        assert_eq!(rec["panelling"]["method"], "pangen");
        assert_eq!(rec["panelling"]["n_nodes"], 160);
        assert_eq!(rec["panelling"]["curvature_bunching"], 1.0);
        assert!(rec["panelling"]["refine_upper"].is_null());
        assert!(rec["panelling"].get("buffer_nodes").is_none());
        // the section's own record is kept underneath
        assert_eq!(rec["designation"], "NACA 2412");
    }

    #[test]
    fn repanel_records_on_a_geometry_without_provenance() {
        let mut g = source();
        g.generator = None;
        let out = repanel(
            &g,
            &PanelConfig {
                n_nodes: 100,
                method: PanelMethod::Cosine(CosineConfig { te_bias: None }),
                ..PanelConfig::default()
            },
        )
        .unwrap();
        assert_eq!(out.x.len(), 101, "the cosine method writes n + 1 nodes");
        let rec = out.generator.unwrap();
        assert_eq!(rec["yfoil"], env!("CARGO_PKG_VERSION"));
        assert_eq!(rec["panelling"]["method"], "cosine");
        assert_eq!(
            rec["panelling"]["te_bias"],
            CosineConfig::DEFAULT_TE_BIAS,
            "the bias actually used"
        );
    }

    #[test]
    fn panel_config_round_trips_flat_json_and_rejects_conflicts() {
        let cfg = PanelConfig {
            n_nodes: 120,
            te_gap: Some(TeGap { gap: 0.002, blend: 0.8 }),
            method: PanelMethod::Pangen(PangenConfig {
                refine_upper: Some([0.2, 0.4]),
                ..PangenConfig::default()
            }),
            ..PanelConfig::default()
        };
        let v = serde_json::to_value(cfg).unwrap();
        assert_eq!(v["method"], "pangen");
        assert_eq!(v["te_gap"]["gap"], 0.002);
        assert_eq!(v["refine_upper"], serde_json::json!([0.2, 0.4]));
        assert!(v.get("buffer_nodes").is_none());
        assert_eq!(PanelConfig::from_json(&v).unwrap(), cfg);

        // a file may omit PPAR keys: XFOIL's defaults apply
        let sparse = serde_json::json!({ "method": "pangen", "n_nodes": 90 });
        assert_eq!(
            PanelConfig::from_json(&sparse).unwrap().method,
            PanelMethod::Pangen(PangenConfig::default())
        );
        let err = |v: serde_json::Value| PanelConfig::from_json(&v).unwrap_err().to_string();
        assert!(
            err(serde_json::json!({ "method": "cosine", "n_nodes": 90, "curvature_bunching": 2.0 }))
                .contains("belongs to the pangen method")
        );
        assert!(
            err(serde_json::json!({ "method": "pangen", "n_nodes": 90, "te_bias": 0.3 }))
                .contains("belongs to the cosine method")
        );
        assert!(err(serde_json::json!({ "method": "pangen", "n_nodes": 90, "typo": 1 }))
            .contains("unknown panelling key typo"));
        assert!(err(
            serde_json::json!({ "method": "pangen", "n_nodes": 90, "sharp_te": true, "te_gap": { "gap": 0.002 } })
        )
        .contains("contradictory"));
        assert!(
            err(serde_json::json!({ "method": "pangen", "n_nodes": 90, "refine_upper": [0.4, 0.2] }))
                .contains("refinement window")
        );
        assert!(err(serde_json::json!({ "n_nodes": 90 })).contains("no \"method\" key"));
    }

    #[test]
    fn te_treatment_follows_the_distribution() {
        let g = source();
        let cfg = PanelConfig {
            n_nodes: 100,
            sharp_te: true,
            ..PanelConfig::default()
        };
        let out = repanel(&g, &cfg).unwrap();
        assert_eq!(
            (out.x[0], out.y[0]),
            (out.x[99], out.y[99]),
            "sharpened after panelling"
        );
        let gap = PanelConfig {
            n_nodes: 100,
            te_gap: Some(TeGap { gap: 0.006, blend: 1.0 }),
            ..PanelConfig::default()
        };
        let out = repanel(&g, &gap).unwrap();
        let d = (out.x[0] - out.x[99]).hypot(out.y[0] - out.y[99]);
        assert!((d - 0.006).abs() < 1e-12, "gap {d}");
        assert!(repanel(
            &g,
            &PanelConfig {
                buffer_nodes: Some(246),
                ..PanelConfig::default()
            }
        )
        .is_err());
    }
    use crate::geometry::Thickness;
    use approx::assert_relative_eq;

    #[test]
    fn test_arc_length_calculation() {
        // Simple square path
        let x = vec![0.0, 1.0, 1.0, 0.0, 0.0];
        let y = vec![0.0, 0.0, 1.0, 1.0, 0.0];
        let s = arc_coordinate(&x, &y);

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

        let s = arc_coordinate(&x, &y);

        // Total arc length should be approximately 2*pi*r
        let total_arc = s[n];
        assert_relative_eq!(total_arc, 2.0 * std::f64::consts::PI * radius, epsilon = 0.1);
    }

    #[test]
    fn test_paneled_airfoil_from_naca0012() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let paneled = panel_foil(&geom);

        // Check basic properties - should produce exactly requested panels
        assert_eq!(paneled.n_foil_nodes, 100);
        assert_relative_eq!(paneled.chord, 1.0, epsilon = 0.05);

        // Arc length should be monotonically increasing
        for i in 1..paneled.s.len() {
            assert!(paneled.s[i] > paneled.s[i - 1]);
        }

        // Total arc length should be roughly 2x chord for thin airfoil
        let total_arc = paneled.s[paneled.n_foil_nodes - 1];
        assert!(total_arc > 1.8 && total_arc < 2.5);

        // Leading edge should be approximately at the midpoint of arc length
        assert!(paneled.s_le > total_arc * 0.3 && paneled.s_le < total_arc * 0.7);
    }

    #[test]
    fn test_normal_vectors_unit_length() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let paneled = panel_foil(&geom);

        // All normal vectors should have unit length
        for i in 0..paneled.n_foil_nodes {
            let mag = (paneled.normal_x[i].powi(2) + paneled.normal_y[i].powi(2)).sqrt();
            assert_relative_eq!(mag, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_normal_vectors_point_outward() {
        use crate::geometry::naca::naca_4digit;

        let geom = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let paneled = panel_foil(&geom);

        // For a symmetric airfoil centered on y=0:
        // - Upper surface (y > 0) normals should have ny > 0
        // - Lower surface (y < 0) normals should have ny < 0
        // (with some tolerance near LE/TE where y is close to 0)

        let mut upper_count = 0;
        let mut lower_count = 0;

        for i in 0..paneled.n_foil_nodes {
            if paneled.y[i] > 0.02 {
                // Upper surface
                assert!(
                    paneled.normal_y[i] > 0.0,
                    "Upper surface normal should point up at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.normal_y[i]
                );
                upper_count += 1;
            } else if paneled.y[i] < -0.02 {
                // Lower surface
                assert!(
                    paneled.normal_y[i] < 0.0,
                    "Lower surface normal should point down at index {}, y={}, ny={}",
                    i,
                    paneled.y[i],
                    paneled.normal_y[i]
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
        let geom = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let paneled = panel_foil(&geom);

        // The NACA generator uses original coefficients for blunt TE (XFOIL-compatible)
        // Gap = 0.00252 which is > 0.0001 * chord, so not sharp
        assert!(!paneled.sharp_te);
    }

    #[test]
    fn test_repanel_preserves_shape() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let repaneled = repanel_cosine(&original, 150, 0.15);

        // Should have approximately the target number of points
        assert!(repaneled.x.len() > 140 && repaneled.x.len() < 160);

        // Extents should be preserved
        let orig_max_x = original.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_x = repaneled.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_x, new_max_x, epsilon = 0.01);

        let orig_min_x = original.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let new_min_x = repaneled.x.iter().cloned().fold(f64::INFINITY, f64::min);
        assert_relative_eq!(orig_min_x, new_min_x, epsilon = 0.01);

        // Maximum thickness should be preserved
        let orig_max_y = original.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let new_max_y = repaneled.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(orig_max_y, new_max_y, epsilon = 0.01);
    }

    #[test]
    fn test_repanel_cosine_spacing() {
        use crate::geometry::naca::naca_4digit;

        let original = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
        let repaneled = repanel_cosine(&original, 100, 0.15);

        // Cosine spacing should cluster points near LE and TE
        // Points near x=0 and x=1 should be closer together than at mid-chord

        // Find spacing near LE (x close to 0)
        let mut le_spacings = Vec::new();
        let mut mid_spacings = Vec::new();

        for i in 1..repaneled.x.len() {
            let dx = (repaneled.x[i] - repaneled.x[i - 1]).abs();
            let dy = (repaneled.y[i] - repaneled.y[i - 1]).abs();
            let ds = (dx * dx + dy * dy).sqrt();

            let avg_x = (repaneled.x[i] + repaneled.x[i - 1]) / 2.0;

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
