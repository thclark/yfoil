//! The foil figure: the aerofoil (and wake) with its panel nodes, boundary-layer quantities drawn
//! normal to the surface at a common scale, and the significant locations (stagnation,
//! transition, separation). Any number of design points (files, or alphas of a polar) can be
//! overlaid with any number of quantities: hue is fixed per quantity, tone and dash pattern per
//! design point.
//!
//! What XFOIL does, and is followed here: PANPLT ticks every aerofoil node outward along the
//! spline normals `NX, NY` by `0.01*CHORD` (`xplots.f:75, 172-174`); CPDISP draws the displacement
//! surface at `X + NX*DSTR` on both sides and splits the wake δ* into upper/lower fractions
//! `DSF1/DSF2` from the TE values, drawing the upper wake edge at `X - N*DSTR*DSF1` and the lower
//! at `X + N*DSTR*DSF2` because the wake normal points to the lower side (`xplots.f:699-754`).
//! XFOIL has no scale factor on δ* and marks no transition or separation on the aerofoil plot; the
//! scaling and the markers are yFoil's.
//!
//! The figure is an ironlab figure (`docs/conventions/figures.md`). ironlab has no equal-aspect
//! axes, so the limits are widened in Rust to the aspect of the plot area the layout gives
//! ([`PLOT_AREA_MM`]), which keeps the aerofoil undistorted at the figure's own size.

use std::collections::BTreeMap;
use std::path::Path;

use ironlab::prelude::*;

use crate::geometry::PanelledFoil;
use crate::output::plot::{
    condition_parameters, describe, point_parameters, regime_labels, rgb, save_figure, shared_labels,
    shared_parameters, PlotError, STUDY,
};
use crate::output::{
    AnalysisOutput, BlQuantity, BoundaryLayerOutput, FoilNodes, PointRecord, PolarOutput, SeparationKind, SideStations,
};
use crate::solver::analysis::FlowConditions;

/// How panel nodes are shown
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PanelStyle {
    /// Plain surface line
    #[default]
    None,
    /// Dark grey surface with a short outward tick at every node (PANPLT)
    Notches,
    /// Mid grey surface with a black dot at every node
    Dots,
}

/// How a BL quantity's magnitude is turned into a normal offset (chord units)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OffsetScale {
    /// Per quantity, over all design points: the largest |value| is drawn `max_offset` chords
    /// off the surface
    Auto { max_offset: f64 },
    /// One multiplier for every quantity (1.0 draws δ*, θ and δ at true geometric size)
    Fixed(f64),
}

impl Default for OffsetScale {
    fn default() -> Self {
        OffsetScale::Auto { max_offset: 0.1 }
    }
}

/// Which significant locations to mark
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MarkerSet {
    pub stagnation: bool,
    pub transition: bool,
    pub separation: bool,
}

impl MarkerSet {
    pub const ALL: MarkerSet = MarkerSet {
        stagnation: true,
        transition: true,
        separation: true,
    };
    pub const NONE: MarkerSet = MarkerSet {
        stagnation: false,
        transition: false,
        separation: false,
    };

    pub fn any(&self) -> bool {
        self.stagnation || self.transition || self.separation
    }
}

impl Default for MarkerSet {
    fn default() -> Self {
        MarkerSet::ALL
    }
}

/// The dash patterns cycled over design points: ironlab's four stroked patterns
pub const DASH_CYCLE: [Dash; 4] = [Dash::Solid, Dash::Dotted, Dash::Dashed, Dash::DashDot];

/// Configuration of the foil figure
#[derive(Debug, Clone)]
pub struct FoilPlotConfig {
    pub title: Option<String>,
    pub panels: PanelStyle,
    /// Notch length as a fraction of chord (PANPLT: 0.01)
    pub notch_length: f64,
    /// Draw the wake nodes and, with quantities, the wake band
    pub show_wake: bool,
    pub quantities: Vec<BlQuantity>,
    pub scale: OffsetScale,
    pub markers: MarkerSet,
    /// Hue per quantity (and per design point for geometry-only input), cycled
    pub palette: Vec<(u8, u8, u8)>,
    /// Lightness multipliers from the first to the last design point
    pub tone_range: (f64, f64),
    /// Dash pattern per design point, cycled
    pub patterns: Vec<Dash>,
}

impl Default for FoilPlotConfig {
    fn default() -> Self {
        Self {
            title: None,
            panels: PanelStyle::None,
            notch_length: 0.01,
            show_wake: false,
            quantities: Vec::new(),
            scale: OffsetScale::default(),
            markers: MarkerSet::ALL,
            // kept clear of the marker side colours (yellow for upper, green for lower)
            palette: vec![
                (0, 100, 200),  // blue
                (200, 50, 50),  // red
                (120, 60, 180), // purple
                (0, 150, 160),  // teal
                (140, 90, 40),  // brown
                (190, 0, 140),  // magenta
            ],
            tone_range: (1.0, 0.55),
            patterns: DASH_CYCLE.to_vec(),
        }
    }
}

/// One design point: a geometry alone, or an operating point with its boundary layer
#[derive(Debug, Clone)]
pub struct DesignPoint {
    pub label: String,
    pub geometry: FoilNodes,
    pub boundary_layer: Option<BoundaryLayerOutput>,
    pub converged: bool,
    /// The aerofoil name, for an analysis
    pub foil: Option<String>,
    /// The flow conditions, for an analysis
    pub conditions: Option<FlowConditions>,
    /// The operating point's record, for an analysis
    pub record: Option<PointRecord>,
}

impl DesignPoint {
    /// Panels only
    pub fn from_geometry(label: impl Into<String>, airfoil: &PanelledFoil) -> Self {
        Self {
            label: label.into(),
            geometry: FoilNodes::from_panelled(airfoil),
            boundary_layer: None,
            converged: true,
            foil: None,
            conditions: None,
            record: None,
        }
    }

    /// One `yfoil analyse -o` record
    pub fn from_analysis(label: impl Into<String>, output: AnalysisOutput) -> Self {
        Self {
            label: label.into(),
            converged: output.results.is_valid(),
            geometry: output.geometry,
            boundary_layer: output.boundary_layer,
            foil: Some(output.foil),
            conditions: Some(output.conditions),
            record: Some(output.results),
        }
    }

    /// The embedded distributions of a `yfoil polar --distributions` record, all of them or the
    /// requested alphas (degrees), labelled `"α=…°"`
    pub fn from_polar(stem: &str, polar: &PolarOutput, alphas: Option<&[f64]>) -> Result<Vec<Self>, PlotError> {
        if polar.distributions.is_empty() {
            return Err(PlotError::Config(format!(
                "{stem}: the polar carries no distributions (run `yfoil polar --distributions`)"
            )));
        }
        let available: Vec<f64> = polar.distributions.iter().map(|d| d.results.alpha_deg).collect();
        let selected: Vec<&AnalysisOutput> = match alphas {
            None => polar.distributions.iter().collect(),
            Some(wanted) => {
                let mut out = Vec::with_capacity(wanted.len());
                for &a in wanted {
                    let found = polar
                        .distributions
                        .iter()
                        .find(|d| (d.results.alpha_deg - a).abs() < 1e-6);
                    match found {
                        Some(d) => out.push(d),
                        None => {
                            let list: Vec<String> = available.iter().map(|a| fmt_alpha(*a)).collect();
                            return Err(PlotError::Config(format!(
                                "{stem}: no distribution at alpha = {a}; available: {}",
                                list.join(", ")
                            )));
                        }
                    }
                }
                out
            }
        };
        Ok(selected
            .into_iter()
            .map(|d| Self::from_analysis(format!("α={}°", fmt_alpha(d.results.alpha_deg)), d.clone()))
            .collect())
    }

    /// The label with a note when the point did not converge
    pub fn display_label(&self) -> String {
        if self.converged {
            self.label.clone()
        } else {
            format!("{} (not converged)", self.label)
        }
    }
}

/// An alpha in degrees for a label: the sweep's `a1 + da*i` round-trips through radians, so
/// 3.0000000000000004 reads as 3
fn fmt_alpha(a: f64) -> String {
    let s = format!("{a:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" {
        "0".to_string()
    } else {
        s.to_string()
    }
}

/// A scale factor for the key: one decimal, in scientific notation outside [1, 10)
fn fmt_scale(k: f64) -> String {
    if (1.0..10.0).contains(&k.abs()) {
        format!("{k:.1}")
    } else {
        format!("{k:.1e}")
    }
}

/// A quantity's symbol for the key, typeset as LaTeX
pub fn quantity_symbol(q: BlQuantity) -> &'static str {
    match q {
        BlQuantity::Dstar => "$\\delta^*$",
        BlQuantity::Theta => "$\\theta$",
        BlQuantity::Delta => "$\\delta$",
        BlQuantity::H => "$H$",
        BlQuantity::Hk => "$H_k$",
        BlQuantity::Hstar => "$H^*$",
        BlQuantity::Ue => "$U_e/V_\\infty$",
        BlQuantity::Cf => "$C_f$",
        BlQuantity::Cdiss => "$C_D$",
        BlQuantity::Sqrtctau => "$\\sqrt{C_\\tau}$ / $N$",
        BlQuantity::Sqrtctaueq => "$\\sqrt{C_{\\tau,\\mathrm{eq}}}$",
        BlQuantity::Us => "$1.6/(1+U_s)$",
        BlQuantity::Cp => "$C_p$",
        BlQuantity::MassDefect => "$m$",
    }
}

/// Design points drawn together must share the panels (BL stations coincide only then)
pub fn check_same_geometry(points: &[DesignPoint]) -> Result<(), PlotError> {
    if let Some(first) = points.first() {
        for p in &points[1..] {
            if !first.geometry.same_panels(&p.geometry) {
                return Err(PlotError::Config(format!(
                    "'{}' and '{}' have different panel geometry; design points overlaid with boundary-layer quantities must share the panels",
                    first.label, p.label
                )));
            }
        }
    }
    Ok(())
}

/// A point in data coordinates
pub type Pt = (f64, f64);

/// Stroke of one (quantity, design point) curve
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveStyle {
    pub rgb: (u8, u8, u8),
    pub dash: Dash,
}

fn tone(rgb: (u8, u8, u8), f: f64) -> (u8, u8, u8) {
    let t = |c: u8| ((c as f64) * f).round().clamp(0.0, 255.0) as u8;
    (t(rgb.0), t(rgb.1), t(rgb.2))
}

/// Lightness multiplier of design point `j` of `n`
fn tone_factor(j: usize, n: usize, range: (f64, f64)) -> f64 {
    if n <= 1 {
        range.0
    } else {
        range.0 + (range.1 - range.0) * (j as f64) / ((n - 1) as f64)
    }
}

/// Hue by quantity index, tone and dash by design point index
pub fn point_style(q_index: usize, point_index: usize, n_points: usize, config: &FoilPlotConfig) -> CurveStyle {
    let base = config.palette[q_index % config.palette.len()];
    CurveStyle {
        rgb: tone(base, tone_factor(point_index, n_points, config.tone_range)),
        dash: config.patterns[point_index % config.patterns.len()],
    }
}

/// One scale factor per quantity over all design points
pub fn scale_factors(points: &[DesignPoint], quantities: &[BlQuantity], scale: OffsetScale) -> Vec<Option<f64>> {
    let chord = points.first().map(|p| p.geometry.chord).unwrap_or(1.0);
    quantities
        .iter()
        .map(|&q| match scale {
            OffsetScale::Fixed(k) => Some(k),
            OffsetScale::Auto { max_offset } => {
                let m = points
                    .iter()
                    .filter_map(|p| p.boundary_layer.as_ref())
                    .map(|bl| bl.max_abs(q))
                    .fold(0.0_f64, f64::max);
                if m > 0.0 {
                    Some(max_offset * chord / m)
                } else {
                    None
                }
            }
        })
        .collect()
}

/// Surface stations of one side offset along the node normals: `X + N*k*q` (CPDISP)
pub fn surface_offset(geom: &FoilNodes, side: &SideStations, col: &[f64], k: f64) -> Vec<(f64, f64)> {
    side.i_node
        .iter()
        .zip(col)
        .map(|(&node, &q)| {
            let i = node - 1;
            (
                geom.x[i] + geom.normal_x[i] * k * q,
                geom.y[i] + geom.normal_y[i] * k * q,
            )
        })
        .collect()
}

/// The two edges of the wake band, each starting at the matching TE surface-offset point:
/// upper `X - N*k*q*f1`, lower `X + N*k*q*f2` (CPDISP; the wake normal points to the lower side)
pub fn wake_band(
    geom: &FoilNodes,
    wake: &SideStations,
    col: &[f64],
    k: f64,
    split: [f64; 2],
    te_upper: (f64, f64),
    te_lower: (f64, f64),
) -> (Vec<Pt>, Vec<Pt>) {
    let n = geom.n();
    let mut upper = vec![te_upper];
    let mut lower = vec![te_lower];
    if let Some(w) = &geom.wake {
        for (&node, &q) in wake.i_node.iter().zip(col) {
            let iw = node - n - 1;
            upper.push((
                w.x[iw] - w.normal_x[iw] * k * q * split[0],
                w.y[iw] - w.normal_y[iw] * k * q * split[0],
            ));
            lower.push((
                w.x[iw] + w.normal_x[iw] * k * q * split[1],
                w.y[iw] + w.normal_y[iw] * k * q * split[1],
            ));
        }
    }
    (upper, lower)
}

const SURFACE_NOTCHED: (u8, u8, u8) = (64, 64, 64);
const SURFACE_DOTTED: (u8, u8, u8) = (112, 112, 112);
const SURFACE_PLAIN: (u8, u8, u8) = (64, 64, 64);
/// The wake line: black, dash-coded per design point (the wake moves with alpha)
const WAKE_LINE: (u8, u8, u8) = (0, 0, 0);
/// Marker colours: one per side (the stagnation point belongs to neither), open shapes per kind
const MARK_UPPER: (u8, u8, u8) = (235, 175, 0);
const MARK_LOWER: (u8, u8, u8) = (0, 165, 0);
const MARK_STAGNATION: (u8, u8, u8) = (40, 40, 40);
/// Marker size (pt)
const MARK_SIZE_PT: f64 = 5.0;

/// The kinds of significant location, each with its own open shape
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum MarkKind {
    Stagnation,
    Transition,
    Forced,
    Separation,
    Reattachment,
}

impl MarkKind {
    fn shape(self) -> Marker {
        match self {
            MarkKind::Stagnation => Marker::Circle,
            MarkKind::Transition => Marker::Square,
            MarkKind::Forced => Marker::Diamond,
            MarkKind::Separation => Marker::TriangleUp,
            MarkKind::Reattachment => Marker::TriangleDown,
        }
    }

    fn label(self) -> &'static str {
        match self {
            MarkKind::Stagnation => "stagnation",
            MarkKind::Transition => "transition",
            MarkKind::Forced => "forced transition",
            MarkKind::Separation => "separation",
            MarkKind::Reattachment => "reattachment",
        }
    }
}

/// A polyline, its stroke and its key entry, in data coordinates
#[derive(Debug)]
struct Curve {
    pts: Vec<Pt>,
    rgb: (u8, u8, u8),
    dash: Dash,
    name: Option<String>,
}

/// A marker at a significant location; `side` is 1 (upper), 2 (lower) or 0 (stagnation)
#[derive(Debug)]
struct Mark {
    at: Pt,
    kind: MarkKind,
    side: usize,
    rgb: (u8, u8, u8),
}

impl Mark {
    /// The key entry of a (kind, side): its shape in its side's colour
    fn key(&self) -> String {
        match self.side {
            1 => format!("{} (upper surface)", self.kind.label()),
            2 => format!("{} (lower surface)", self.kind.label()),
            _ => self.kind.label().to_string(),
        }
    }
}

/// Everything laid out in data coordinates before the figure exists
#[derive(Debug)]
struct Layout {
    surfaces: Vec<Curve>,
    /// One wake line per design point, black with the point's dash pattern
    wakes: Vec<Curve>,
    /// Node ticks: one polyline of two-point segments per surface (or wake)
    notches: Vec<Vec<(Pt, Pt)>>,
    dots: Vec<Vec<Pt>>,
    curves: Vec<Curve>,
    markers: Vec<Mark>,
}

impl Layout {
    /// The curve key in drawing order
    fn curve_key(&self) -> Vec<&str> {
        self.surfaces
            .iter()
            .chain(&self.wakes)
            .chain(&self.curves)
            .filter_map(|c| c.name.as_deref())
            .collect()
    }

    /// Every point drawn, for the bounds
    fn points(&self) -> impl Iterator<Item = Pt> + '_ {
        self.surfaces
            .iter()
            .chain(&self.wakes)
            .chain(&self.curves)
            .flat_map(|c| c.pts.iter().copied())
            .chain(self.notches.iter().flatten().flat_map(|&(a, b)| [a, b]))
            .chain(self.markers.iter().map(|m| m.at))
    }
}

fn surface_rgb(j: usize, n_points: usize, geometry_only: bool, config: &FoilPlotConfig) -> (u8, u8, u8) {
    if geometry_only && n_points > 1 {
        config.palette[j % config.palette.len()]
    } else {
        match config.panels {
            PanelStyle::Notches => SURFACE_NOTCHED,
            PanelStyle::Dots => SURFACE_DOTTED,
            PanelStyle::None => SURFACE_PLAIN,
        }
    }
}

fn layout(points: &[DesignPoint], config: &FoilPlotConfig) -> Result<Layout, PlotError> {
    if points.is_empty() {
        return Err(PlotError::Config("nothing to plot: no design points".to_string()));
    }
    let geometry_only = points.iter().all(|p| p.boundary_layer.is_none());
    if !config.quantities.is_empty() {
        if let Some(p) = points.iter().find(|p| p.boundary_layer.is_none()) {
            return Err(PlotError::Config(format!(
                "'{}' has no boundary layer to plot (geometry or inviscid input)",
                p.label
            )));
        }
        check_same_geometry(points)?;
    }
    if config.show_wake {
        if let Some(p) = points.iter().find(|p| p.geometry.wake.is_none()) {
            return Err(PlotError::Config(format!(
                "'{}' has no wake geometry (only a viscous analysis builds the wake)",
                p.label
            )));
        }
    }

    let n_points = points.len();
    let chord = points[0].geometry.chord;
    let mut lay = Layout {
        surfaces: Vec::new(),
        wakes: Vec::new(),
        notches: Vec::new(),
        dots: Vec::new(),
        curves: Vec::new(),
        markers: Vec::new(),
    };

    // Aerofoil outlines: one per geometry-only point; identical panels otherwise, so once
    let notch = config.notch_length * chord;
    let surface_points: Vec<usize> = if geometry_only {
        (0..n_points).collect()
    } else {
        vec![0]
    };
    for &j in &surface_points {
        let g = &points[j].geometry;
        let outline: Vec<Pt> = g.x.iter().copied().zip(g.y.iter().copied()).collect();
        match config.panels {
            PanelStyle::Notches => lay.notches.push(
                (0..g.n())
                    .map(|i| {
                        (
                            (g.x[i], g.y[i]),
                            (g.x[i] + notch * g.normal_x[i], g.y[i] + notch * g.normal_y[i]),
                        )
                    })
                    .collect(),
            ),
            PanelStyle::Dots => lay.dots.push(outline.clone()),
            PanelStyle::None => {}
        }
        lay.surfaces.push(Curve {
            pts: outline,
            rgb: surface_rgb(j, n_points, geometry_only, config),
            dash: Dash::Solid,
            name: (geometry_only && n_points > 1).then(|| points[j].display_label()),
        });
    }

    // Wakes: XYWAKE rebuilds the wake for every alpha, so each design point carries its own;
    // black, with the point's dash pattern. Wake panels are shown only for a single point.
    // Several points with no quantities are told apart by their wakes' dash patterns.
    if config.show_wake {
        let key_wakes = config.quantities.is_empty() && n_points > 1;
        for (j, p) in points.iter().enumerate() {
            let Some(w) = &p.geometry.wake else { continue };
            let line: Vec<Pt> = w.x.iter().copied().zip(w.y.iter().copied()).collect();
            if n_points == 1 {
                match config.panels {
                    PanelStyle::Notches => lay.notches.push(
                        (0..w.x.len())
                            .map(|i| {
                                (
                                    (
                                        w.x[i] - 0.5 * notch * w.normal_x[i],
                                        w.y[i] - 0.5 * notch * w.normal_y[i],
                                    ),
                                    (
                                        w.x[i] + 0.5 * notch * w.normal_x[i],
                                        w.y[i] + 0.5 * notch * w.normal_y[i],
                                    ),
                                )
                            })
                            .collect(),
                    ),
                    PanelStyle::Dots => lay.dots.push(line.clone()),
                    PanelStyle::None => {}
                }
            }
            lay.wakes.push(Curve {
                pts: line,
                rgb: WAKE_LINE,
                dash: config.patterns[j % config.patterns.len()],
                name: key_wakes.then(|| p.display_label()),
            });
        }
    }

    // BL quantities
    let ks = scale_factors(points, &config.quantities, config.scale);
    for (qi, (&q, k)) in config.quantities.iter().zip(&ks).enumerate() {
        let Some(k) = *k else {
            log::warn!("{}: every value is zero, not drawn", q.label());
            continue;
        };
        for (j, p) in points.iter().enumerate() {
            let bl = p.boundary_layer.as_ref().expect("checked above");
            let g = &p.geometry;
            let st = point_style(qi, j, n_points, config);
            let curve = |pts: Vec<Pt>, name: Option<String>| Curve {
                pts,
                rgb: st.rgb,
                dash: st.dash,
                name,
            };
            let upper = surface_offset(g, &bl.upper, bl.upper.column(q), k);
            let lower = surface_offset(g, &bl.lower, bl.lower.column(q), k);
            let (te_u, te_l) = (
                *upper.last().unwrap_or(&(g.x_te, g.y_te)),
                *lower.last().unwrap_or(&(g.x_te, g.y_te)),
            );
            let name = if n_points > 1 {
                format!(
                    "{} at {}, scaled ×{}",
                    quantity_symbol(q),
                    p.display_label(),
                    fmt_scale(k)
                )
            } else {
                format!("{}, scaled ×{}", quantity_symbol(q), fmt_scale(k))
            };
            lay.curves.push(curve(upper, Some(name)));
            lay.curves.push(curve(lower, None));
            if config.show_wake && g.wake.is_some() && !bl.wake.is_empty() {
                let split = if q == BlQuantity::Dstar {
                    bl.wake_split
                } else {
                    [0.5, 0.5]
                };
                let (wu, wl) = wake_band(g, &bl.wake, bl.wake.column(q), k, split, te_u, te_l);
                lay.curves.push(curve(wu, None));
                lay.curves.push(curve(wl, None));
            }
        }
    }

    // Markers: open shape per kind, colour per side (toned per design point)
    if config.markers.any() {
        let side_colour = |is: usize| if is == 1 { MARK_UPPER } else { MARK_LOWER };
        for (j, p) in points.iter().enumerate() {
            let Some(bl) = &p.boundary_layer else { continue };
            let f = tone_factor(j, n_points, config.tone_range);
            let mut mark = |x: f64, y: f64, kind: MarkKind, side: usize, rgb: (u8, u8, u8)| {
                lay.markers.push(Mark {
                    at: (x, y),
                    kind,
                    side,
                    rgb: tone(rgb, f),
                });
            };
            if config.markers.stagnation {
                let m = &bl.stagnation;
                mark(m.x, m.y, MarkKind::Stagnation, 0, MARK_STAGNATION);
            }
            if config.markers.transition {
                for (is, m) in [1, 2].into_iter().zip(&bl.transition) {
                    let kind = if m.forced {
                        MarkKind::Forced
                    } else {
                        MarkKind::Transition
                    };
                    mark(m.x_transition, m.y_transition, kind, is, side_colour(is));
                }
            }
            if config.markers.separation {
                for m in &bl.derived_separation {
                    let kind = match m.kind {
                        SeparationKind::Separation => MarkKind::Separation,
                        SeparationKind::Reattachment => MarkKind::Reattachment,
                    };
                    mark(m.x, m.y, kind, m.side, side_colour(m.side));
                }
            }
        }
    }

    Ok(lay)
}

/// The plot area (width, height in mm) the foil figure's layout gives its axes, without and with
/// the key row: the default 160 × 100 mm figure less the title, the tick and axis labels and the
/// legend row below. ironlab has no equal-aspect axes, so the limits are widened to this aspect
/// instead. The areas were measured on ironlab 0.12.0 renders; the tick labels' widths move them
/// by a millimetre or so, so the aspect is equal to within about one per cent.
pub const PLOT_AREA_MM: [(f64, f64); 2] = [(142.0, 80.5), (142.0, 57.5)];

/// The limits that show `bounds` (x_min, x_max, y_min, y_max) at one scale on both axes in a plot
/// area of `area` (width, height): the short side is widened about its centre
pub fn equal_aspect_limits(bounds: (f64, f64, f64, f64), area: (f64, f64)) -> ((f64, f64), (f64, f64)) {
    let (x0, x1, y0, y1) = bounds;
    let (dx, dy) = (x1 - x0, y1 - y0);
    let aspect = area.1 / area.0;
    if dy / dx < aspect {
        let half = 0.5 * dx * aspect;
        let c = 0.5 * (y0 + y1);
        ((x0, x1), (c - half, c + half))
    } else {
        let half = 0.5 * dy / aspect;
        let c = 0.5 * (x0 + x1);
        ((c - half, c + half), (y0, y1))
    }
}

/// The parameters of the foil figure: provenance, and for one design point its inputs and
/// outputs; for several, the conditions they share
fn foil_parameters(points: &[DesignPoint]) -> (BTreeMap<String, Parameter>, Vec<&'static str>) {
    let mut parameters = BTreeMap::from([
        ("study".to_string(), Parameter::from(STUDY)),
        ("figure".to_string(), Parameter::from("foil")),
        ("points".to_string(), Parameter::from(points.len() as i64)),
    ]);
    let per_point: Vec<BTreeMap<String, Parameter>> = points
        .iter()
        .map(|p| {
            let mut m = BTreeMap::from([("n_nodes".to_string(), Parameter::from(p.geometry.n() as i64))]);
            if let Some(foil) = &p.foil {
                m.insert("foil".to_string(), Parameter::from(foil.as_str()));
            }
            if let Some(c) = &p.conditions {
                m.extend(condition_parameters(c));
            }
            m
        })
        .collect();
    parameters.extend(shared_parameters(&per_point));
    if let [p] = points {
        if let Some(record) = &p.record {
            parameters.extend(point_parameters(record, p.boundary_layer.as_ref()));
        }
    }
    let labels: Vec<Vec<&'static str>> = points
        .iter()
        .map(|p| p.conditions.as_ref().map(regime_labels).unwrap_or_default())
        .collect();
    (parameters, shared_labels(&labels))
}

/// Build the foil figure
pub fn foil_figure(points: &[DesignPoint], config: &FoilPlotConfig) -> Result<Figure, PlotError> {
    let lay = layout(points, config)?;
    let chord = points[0].geometry.chord;

    // bounds: everything drawn, padded by 5% chord, then widened to the plot area's aspect
    let (mut x0, mut x1, mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for (x, y) in lay.points().filter(|(x, y)| x.is_finite() && y.is_finite()) {
        (x0, x1, y0, y1) = (x0.min(x), x1.max(x), y0.min(y), y1.max(y));
    }
    let pad = 0.05 * chord;
    let named = !lay.curve_key().is_empty() || !lay.markers.is_empty();
    let area = PLOT_AREA_MM[usize::from(named)];
    let (xlim, ylim) = equal_aspect_limits((x0 - pad, x1 + pad, y0 - pad, y1 + pad), area);

    let title = config.title.clone().unwrap_or_else(|| {
        if points.len() == 1 {
            points[0].display_label()
        } else {
            "Foil".to_string()
        }
    });
    let mut fig = Figure::new().title(Text::plain(title));
    if named {
        fig = fig.tiles(2, 1).row_weights([3.0, 1.0]);
    }
    let (parameters, labels) = foil_parameters(points);
    fig = describe(fig, parameters, &labels);
    let key = named.then(|| fig.panel_legend(1, 0).id());

    let mut ax = fig.axes(0, 0);
    let draw = |ax: &mut AxesMut, c: &Curve, width: f64| {
        let (x, y): (Vec<f64>, Vec<f64>) = c.pts.iter().copied().unzip();
        let mut line = ax.plot(x, y);
        line.color(rgb(c.rgb)).dash(c.dash).line_width(width);
        if let Some(name) = &c.name {
            line.display_name(name.as_str());
        }
    };
    for c in &lay.surfaces {
        draw(&mut ax, c, 1.0);
    }
    for c in &lay.wakes {
        draw(&mut ax, c, 1.0);
    }
    // node ticks: one line per surface, its segments separated by NaN
    for segs in &lay.notches {
        let mut x = Vec::with_capacity(3 * segs.len());
        let mut y = Vec::with_capacity(3 * segs.len());
        for &(a, b) in segs {
            x.extend([a.0, b.0, f64::NAN]);
            y.extend([a.1, b.1, f64::NAN]);
        }
        ax.plot(x, y).color(rgb(SURFACE_NOTCHED)).line_width(0.5);
    }
    for d in &lay.dots {
        let (x, y): (Vec<f64>, Vec<f64>) = d.iter().copied().unzip();
        ax.plot(x, y)
            .color(Color::BLACK)
            .dash(Dash::None)
            .marker(Marker::Point)
            .marker_size(1.5);
    }
    for c in &lay.curves {
        draw(&mut ax, c, 1.0);
    }
    // markers: one artist each; the first of each (kind, side) carries its key entry
    let mut keyed = Vec::new();
    for m in &lay.markers {
        let mut line = ax.plot([m.at.0], [m.at.1]);
        line.color(rgb(m.rgb))
            .dash(Dash::None)
            .marker(m.kind.shape())
            .marker_size(MARK_SIZE_PT)
            .line_width(1.0);
        if !keyed.contains(&(m.kind, m.side)) {
            keyed.push((m.kind, m.side));
            line.display_name(m.key());
        }
    }
    ax.xlabel("$x/c$")
        .ylabel("$y/c$")
        .xlim(xlim.0, xlim.1)
        .ylim(ylim.0, ylim.1)
        .grid(true);
    if let Some(key) = key {
        ax.legend_in(key);
    }
    Ok(fig)
}

/// Build and write the foil figure where `path` names it (`.fig` or `.fig.json`)
pub fn plot_foil<P: AsRef<Path>>(points: &[DesignPoint], path: P, config: &FoilPlotConfig) -> Result<(), PlotError> {
    save_figure(&foil_figure(points, config)?, path.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{naca_4digit, panel_foil, Thickness};
    use crate::output::{SideStations, WakeNodes};
    use ironlab::ir::Artist;

    fn geometry(n: usize) -> FoilNodes {
        let geom = naca_4digit("0012", n, Thickness::Perpendicular).unwrap();
        FoilNodes::from_panelled(&panel_foil(&geom))
    }

    /// A synthetic boundary layer on a geometry: every side's column is `value` at every node
    fn synthetic_bl(g: &mut FoilNodes, value: f64) -> BoundaryLayerOutput {
        let n = g.n();
        let nw = 5;
        // a straight wake behind the TE, normal pointing down (as XYWAKE's does)
        g.wake = Some(WakeNodes {
            x: (1..=nw).map(|i| g.x_te + 0.2 * i as f64).collect(),
            y: vec![g.y_te; nw],
            s: (1..=nw).map(|i| g.s[n - 1] + 0.2 * i as f64).collect(),
            normal_x: vec![0.0; nw],
            normal_y: vec![-1.0; nw],
        });
        let le = n / 2;
        let (xw, yw) = {
            let w = g.wake.as_ref().unwrap();
            (w.x.clone(), w.y.clone())
        };
        let side = |nodes: Vec<usize>| {
            let mut s = SideStations::default();
            for (k, node) in nodes.iter().enumerate() {
                s.i_station.push(k + 2);
                s.i_node.push(*node);
                if *node <= n {
                    s.x.push(g.x[node - 1]);
                    s.y.push(g.y[node - 1]);
                } else {
                    s.x.push(xw[node - n - 1]);
                    s.y.push(yw[node - n - 1]);
                }
                s.xi.push(k as f64);
                for col in [
                    &mut s.primaries.ue,
                    &mut s.primaries.theta,
                    &mut s.primaries.dstar,
                    &mut s.primaries.sqrtctau,
                    &mut s.primaries.mass_defect,
                    &mut s.closures.ue_compressible,
                    &mut s.closures.h,
                    &mut s.closures.hk,
                    &mut s.closures.hstar,
                    &mut s.closures.cf,
                    &mut s.closures.cdiss,
                    &mut s.closures.delta,
                    &mut s.closures.sqrtctaueq,
                    &mut s.closures.us,
                    &mut s.closures.retheta,
                    &mut s.closures.machsqd_edge,
                    &mut s.cp,
                ] {
                    col.push(value);
                }
            }
            s
        };
        // upper: LE → TE node 1 (XFOIL's side 1 marches from the stagnation point to the TE)
        let upper = side((1..=le + 1).rev().collect());
        let lower = side((le + 1..=n).collect());
        let wake = side((n + 1..=n + nw).collect());
        BoundaryLayerOutput {
            i_te_station: [upper.len() + 1, lower.len() + 1],
            i_transition_station: [3, 3],
            n_wake_nodes: nw,
            qinf: 1.0,
            te_thickness_normal: 0.0,
            stagnation: crate::output::StagnationMarker {
                i_stagnation_node: le,
                s_stagnation: g.s_le,
                x: g.x_le,
                y: g.y_le,
            },
            transition: [
                crate::output::TransitionMarker {
                    i_station: 3,
                    forced: false,
                    x_transition: 0.3,
                    y_transition: 0.05,
                    s_transition: g.s_le - 0.3,
                },
                crate::output::TransitionMarker {
                    i_station: 3,
                    forced: true,
                    x_transition: 0.6,
                    y_transition: -0.04,
                    s_transition: g.s_le + 0.6,
                },
            ],
            derived_separation: Vec::new(),
            wake_split: [0.6, 0.4],
            upper,
            lower,
            wake,
        }
    }

    fn point(label: &str, value: f64, converged: bool) -> DesignPoint {
        let mut g = geometry(60);
        let bl = synthetic_bl(&mut g, value);
        DesignPoint {
            label: label.to_string(),
            geometry: g,
            boundary_layer: Some(bl),
            converged,
            foil: None,
            conditions: None,
            record: None,
        }
    }

    #[test]
    fn style_is_hue_by_quantity_and_tone_pattern_by_point() {
        let c = FoilPlotConfig::default();
        for hue in &c.palette {
            for side in [MARK_UPPER, MARK_LOWER] {
                let d = (hue.0 as i32 - side.0 as i32).abs()
                    + (hue.1 as i32 - side.1 as i32).abs()
                    + (hue.2 as i32 - side.2 as i32).abs();
                assert!(d > 120, "palette colour {hue:?} is too close to marker colour {side:?}");
            }
        }
        let a = point_style(0, 0, 3, &c);
        let b = point_style(1, 0, 3, &c);
        assert_ne!(a.rgb, b.rgb, "different quantities: different hue");
        assert_eq!(a.dash, b.dash, "same design point: same pattern");
        let a2 = point_style(0, 2, 3, &c);
        assert_eq!(a2.dash, Dash::Dashed);
        assert_eq!(a2.rgb, tone(c.palette[0], 0.55), "last point gets the dark tone");
        let s0 = point_style(0, 0, 1, &c);
        assert_eq!(s0.rgb, c.palette[0], "single point: base hue");
        assert_eq!(s0.dash, Dash::Solid);
        assert_eq!(point_style(0, 4, 6, &c).dash, Dash::Solid, "cycle wraps after four");
        assert_eq!(point_style(7, 0, 1, &c).rgb, c.palette[1], "palette cycles");
    }

    #[test]
    fn scale_is_one_factor_per_quantity_over_all_points() {
        let pts = [point("a", 0.01, true), point("b", 0.04, true)];
        let ks = scale_factors(
            &pts,
            &[BlQuantity::Dstar, BlQuantity::H],
            OffsetScale::Auto { max_offset: 0.1 },
        );
        let chord = pts[0].geometry.chord;
        assert!((ks[0].unwrap() - 0.1 * chord / 0.04).abs() < 1e-12);
        assert!((ks[1].unwrap() - 0.1 * chord / 0.04).abs() < 1e-12);
        let ks = scale_factors(&pts, &[BlQuantity::Dstar], OffsetScale::Fixed(1.0));
        assert_eq!(ks, vec![Some(1.0)]);
        let zero = [point("z", 0.0, true)];
        assert_eq!(
            scale_factors(&zero, &[BlQuantity::Cf], OffsetScale::default()),
            vec![None]
        );
    }

    #[test]
    fn surface_offset_follows_node_normals() {
        let p = point("a", 0.02, true);
        let bl = p.boundary_layer.as_ref().unwrap();
        let k = 3.0;
        let up = surface_offset(&p.geometry, &bl.upper, bl.upper.column(BlQuantity::Dstar), k);
        assert_eq!(up.len(), bl.upper.len());
        for (pt, &node) in up.iter().zip(&bl.upper.i_node) {
            let i = node - 1;
            assert!((pt.0 - (p.geometry.x[i] + p.geometry.normal_x[i] * k * 0.02)).abs() < 1e-15);
            assert!((pt.1 - (p.geometry.y[i] + p.geometry.normal_y[i] * k * 0.02)).abs() < 1e-15);
        }
    }

    #[test]
    fn wake_band_starts_at_te_offsets_and_uses_split() {
        let p = point("a", 0.02, true);
        let g = &p.geometry;
        let bl = p.boundary_layer.as_ref().unwrap();
        let te_u = (10.0, 1.0);
        let te_l = (10.0, -1.0);
        let (wu, wl) = wake_band(
            g,
            &bl.wake,
            bl.wake.column(BlQuantity::Dstar),
            2.0,
            [0.6, 0.4],
            te_u,
            te_l,
        );
        assert_eq!(wu[0], te_u);
        assert_eq!(wl[0], te_l);
        assert_eq!(wu.len(), bl.wake.len() + 1);
        let w = g.wake.as_ref().unwrap();
        // wake normal is (0, -1): upper edge = X - N*k*q*f1 lies *above* the wake line
        assert!((wu[1].1 - (w.y[0] + 2.0 * 0.02 * 0.6)).abs() < 1e-15);
        assert!((wl[1].1 - (w.y[0] - 2.0 * 0.02 * 0.4)).abs() < 1e-15);
    }

    #[test]
    fn legend_names_point_and_quantity_in_the_grid_case() {
        let pts = [point("a", 0.01, true), point("b", 0.02, false)];
        let config = FoilPlotConfig {
            quantities: vec![BlQuantity::Dstar, BlQuantity::Theta],
            show_wake: true,
            panels: PanelStyle::Notches,
            ..Default::default()
        };
        let lay = layout(&pts, &config).unwrap();
        let texts = lay.curve_key();
        assert_eq!(texts.len(), 4);
        assert!(texts[0].starts_with("$\\delta^*$ at a, scaled ×"));
        assert!(texts[1].starts_with("$\\delta^*$ at b (not converged), scaled ×"));
        assert!(texts[2].starts_with("$\\theta$ at a, scaled ×"));
        // four distinct (colour, dash) combinations among the named curves
        let mut styles: Vec<String> = lay
            .curves
            .iter()
            .filter(|c| c.name.is_some())
            .map(|c| format!("{:?}{:?}", c.rgb, c.dash))
            .collect();
        styles.sort();
        styles.dedup();
        assert_eq!(styles.len(), 4);
        // 2 points × 2 quantities × (upper, lower, wake upper, wake lower)
        assert_eq!(lay.curves.len(), 16);
        // aerofoil notches once; wake panels only for a single point; one wake line per point,
        // black, dash-coded like its curves
        assert_eq!(lay.notches.len(), 1);
        assert_eq!(lay.notches[0].len(), 60);
        assert_eq!(lay.wakes.len(), 2);
        assert!(lay.wakes.iter().all(|w| w.rgb == WAKE_LINE));
        assert_eq!(lay.wakes[0].dash, Dash::Solid);
        assert_eq!(lay.wakes[1].dash, Dash::Dotted);
        // single point: no label prefix, and the wake gets its panels
        let one = layout(&pts[..1], &config).unwrap();
        assert!(one.curve_key()[0].starts_with("$\\delta^*$, scaled ×"));
        assert_eq!(fmt_scale(3.0071), "3.0");
        assert_eq!(fmt_scale(0.4567), "4.6e-1");
        assert_eq!(fmt_scale(17.1), "1.7e1");
        assert_eq!(fmt_scale(1.0), "1.0");
        assert_eq!(one.notches.len(), 2);
        assert_eq!(one.notches[1].len(), 5);
        assert_eq!(one.wakes.len(), 1);
    }

    #[test]
    fn errors_are_configuration_errors() {
        let g = DesignPoint::from_geometry(
            "g",
            &panel_foil(&naca_4digit("0012", 60, Thickness::Perpendicular).unwrap()),
        );
        let config = FoilPlotConfig {
            quantities: vec![BlQuantity::Dstar],
            ..Default::default()
        };
        assert!(matches!(
            layout(std::slice::from_ref(&g), &config),
            Err(PlotError::Config(_))
        ));
        let wake_cfg = FoilPlotConfig {
            show_wake: true,
            ..Default::default()
        };
        assert!(matches!(
            layout(std::slice::from_ref(&g), &wake_cfg),
            Err(PlotError::Config(_))
        ));
        // geometry-only input renders with no BL, and two panelings get their own colours
        let g2 = DesignPoint::from_geometry(
            "g2",
            &panel_foil(&naca_4digit("0012", 80, Thickness::Perpendicular).unwrap()),
        );
        let lay = layout(&[g.clone(), g2], &FoilPlotConfig::default()).unwrap();
        assert_eq!(lay.surfaces.len(), 2);
        assert_ne!(lay.surfaces[0].rgb, lay.surfaces[1].rgb);
        assert_eq!(lay.curve_key(), vec!["g", "g2"]);
        // mismatched panels with quantities
        let mut b = point("b", 0.01, true);
        b.geometry.x[3] = f64::from_bits(b.geometry.x[3].to_bits() + 1);
        let err = layout(&[point("a", 0.01, true), b], &config).unwrap_err();
        assert!(err.to_string().contains("different panel geometry"));
    }

    #[test]
    fn equal_aspect_widens_the_short_side_about_its_centre() {
        let ((x0, x1), (y0, y1)) = equal_aspect_limits((0.0, 2.0, -0.1, 0.1), (100.0, 50.0));
        assert_eq!((x0, x1), (0.0, 2.0));
        assert!((y0 + 0.5).abs() < 1e-15 && (y1 - 0.5).abs() < 1e-15);
        let ((x0, x1), (y0, y1)) = equal_aspect_limits((0.0, 1.0, 0.0, 2.0), (100.0, 50.0));
        assert_eq!((y0, y1), (0.0, 2.0));
        assert!((x0 + 1.5).abs() < 1e-15 && (x1 - 2.5).abs() < 1e-15);
    }

    #[test]
    fn figure_has_markers_dots_and_one_keyed_entry_per_kind_and_side() {
        let pts = [point("a", 0.01, true)];
        let config = FoilPlotConfig {
            quantities: vec![BlQuantity::Dstar],
            panels: PanelStyle::Dots,
            ..Default::default()
        };
        let fig = foil_figure(&pts, &config).unwrap();
        assert!(fig.validate().errors.is_empty());
        let ir = fig.ir();
        assert_eq!(ir.axes.len(), 1);
        assert_eq!(ir.panel_legends.len(), 1);
        let lines: Vec<_> = ir.axes[0]
            .artists
            .iter()
            .filter_map(|a| match a {
                Artist::Line(l) => Some(l),
                _ => None,
            })
            .collect();
        // surface + dots + (upper, lower) + stagnation, transition, forced transition
        assert_eq!(lines.len(), 1 + 1 + 2 + 3);
        let names: Vec<&str> = lines
            .iter()
            .filter_map(|l| l.display_name.as_ref().map(|t| t.content.as_str()))
            .collect();
        assert_eq!(names.len(), 4);
        assert!(names[0].starts_with("$\\delta^*$, scaled ×"));
        assert_eq!(
            &names[1..],
            [
                "stagnation",
                "transition (upper surface)",
                "forced transition (lower surface)"
            ]
        );
        let shapes: Vec<Marker> = lines.iter().map(|l| l.marker.shape).collect();
        assert_eq!(shapes[1], Marker::Point);
        assert_eq!(&shapes[4..], [Marker::Circle, Marker::Square, Marker::Diamond]);
        assert_eq!(fig.parameters().get("figure"), Some(&Parameter::from("foil")));
        assert_eq!(fmt_alpha(3.0000000000000004), "3");
        assert_eq!(fmt_alpha(-2.5), "-2.5");
        assert_eq!(fmt_alpha(-0.0), "0");
        assert_eq!(fmt_alpha(12.000000000000002), "12");
    }
}
