//! The foil plot: the airfoil (and wake) with its panel nodes, boundary-layer quantities drawn
//! normal to the surface at a common scale, and the significant locations (stagnation,
//! transition, separation). Any number of design points (files, or alphas of a polar) can be
//! overlaid with any number of quantities: hue is fixed per quantity, tone and dash pattern per
//! design point.
//!
//! What XFOIL does, and is followed here: PANPLT ticks every airfoil node outward along the
//! spline normals `NX, NY` by `0.01*CHORD` (`xplots.f:75, 172-174`); CPDISP draws the displacement
//! surface at `X + NX*DSTR` on both sides and splits the wake δ* into upper/lower fractions
//! `DSF1/DSF2` from the TE values, drawing the upper wake edge at `X - N*DSTR*DSF1` and the lower
//! at `X + N*DSTR*DSF2` because the wake normal points to the lower side (`xplots.f:699-754`).
//! XFOIL has no scale factor on δ* and marks no transition or separation on the airfoil plot; the
//! scaling and the markers are yFoil's.

use std::path::Path;

use crate::geometry::PanelledFoil;
use crate::output::plot::PlotError;
use crate::output::svg::{
    write_png, write_svg, Canvas, LegendAnchor, LegendEntry, LegendKey, LineStyle, Margins, MarkerShape,
};
use crate::output::{
    AnalysisOutput, BlQuantity, BoundaryLayerOutput, FoilNodes, PolarOutput, SeparationKind, SideStations,
};

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

/// Dash patterns, cycled over design points
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DashPattern {
    Solid,
    Dot,
    Dash,
    DotDash,
    DotDotDash,
}

impl DashPattern {
    pub const CYCLE: [DashPattern; 5] = [
        DashPattern::Solid,
        DashPattern::Dot,
        DashPattern::Dash,
        DashPattern::DotDash,
        DashPattern::DotDotDash,
    ];

    /// SVG `stroke-dasharray`, `None` for solid
    pub fn dasharray(self) -> Option<&'static str> {
        match self {
            DashPattern::Solid => None,
            DashPattern::Dot => Some("2,3"),
            DashPattern::Dash => Some("8,4"),
            DashPattern::DotDash => Some("8,4,2,4"),
            DashPattern::DotDotDash => Some("8,4,2,4,2,4"),
        }
    }
}

/// Configuration of the foil plot
#[derive(Debug, Clone)]
pub struct FoilPlotConfig {
    pub width: u32,
    pub height: u32,
    pub title: Option<String>,
    pub background: (u8, u8, u8),
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
    pub patterns: Vec<DashPattern>,
}

impl Default for FoilPlotConfig {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 500,
            title: None,
            background: (255, 255, 255),
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
            patterns: DashPattern::CYCLE.to_vec(),
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
}

impl DesignPoint {
    /// Panels only
    pub fn from_geometry(label: impl Into<String>, airfoil: &PanelledFoil) -> Self {
        Self {
            label: label.into(),
            geometry: FoilNodes::from_panelled(airfoil),
            boundary_layer: None,
            converged: true,
        }
    }

    /// One `yfoil analyze -o` record
    pub fn from_analysis(label: impl Into<String>, output: AnalysisOutput) -> Self {
        Self {
            label: label.into(),
            converged: output.results.is_converged(),
            geometry: output.geometry,
            boundary_layer: output.boundary_layer,
        }
    }

    /// The embedded distributions of a `yfoil polar --distributions` record, all of them or the
    /// requested alphas (degrees), labelled `"{stem} α=…°"`
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

/// Output image format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    Svg,
    Png,
}

impl ImageFormat {
    /// From a file extension (`svg` unless the extension is `png`)
    pub fn from_path(path: &Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some(e) if e.eq_ignore_ascii_case("png") => ImageFormat::Png,
            _ => ImageFormat::Svg,
        }
    }
}

/// A point in data coordinates
pub type Pt = (f64, f64);
/// A segment in data coordinates
type Seg = (Pt, Pt);

/// Stroke of one (quantity, design point) curve
#[derive(Debug, Clone, PartialEq)]
pub struct CurveStyle {
    pub rgb: (u8, u8, u8),
    pub dasharray: Option<String>,
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
    let pattern = config.patterns[point_index % config.patterns.len()];
    CurveStyle {
        rgb: tone(base, tone_factor(point_index, n_points, config.tone_range)),
        dasharray: pattern.dasharray().map(str::to_string),
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
/// Neutral stroke for the shape keys of the marker legend
const MARK_KEY: (u8, u8, u8) = (90, 90, 90);
/// Marker half-size in pixels
const MARK_R: f64 = 4.5;

/// The kinds of significant location, each with its own open shape
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MarkKind {
    Stagnation,
    Transition,
    Forced,
    Separation,
    Reattachment,
}

impl MarkKind {
    const ALL: [MarkKind; 5] = [
        MarkKind::Stagnation,
        MarkKind::Transition,
        MarkKind::Forced,
        MarkKind::Separation,
        MarkKind::Reattachment,
    ];

    fn shape(self) -> MarkerShape {
        match self {
            MarkKind::Stagnation => MarkerShape::Circle,
            MarkKind::Transition => MarkerShape::Square,
            MarkKind::Forced => MarkerShape::Diamond,
            MarkKind::Separation => MarkerShape::Triangle,
            MarkKind::Reattachment => MarkerShape::TriangleDown,
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

/// A polyline plus its style, in data coordinates
#[derive(Debug)]
struct Curve {
    pts: Vec<Pt>,
    style: LineStyle,
}

/// A drawn airfoil outline
#[derive(Debug)]
struct Surface {
    outline: Vec<Pt>,
    style: LineStyle,
}

/// A marker at a significant location
#[derive(Debug)]
struct Mark {
    at: Pt,
    kind: MarkKind,
    rgb: (u8, u8, u8),
}

/// Everything laid out in data coordinates before the canvas exists
#[derive(Debug)]
struct Layout {
    surfaces: Vec<Surface>,
    /// One wake line per design point, black with the point's dash pattern
    wakes: Vec<Curve>,
    notches: Vec<Vec<Seg>>,
    dots: Vec<Vec<Pt>>,
    curves: Vec<Curve>,
    markers: Vec<Mark>,
    legend: Vec<LegendEntry>,
    /// The second key, bottom right: marker shapes and side colours
    marker_legend: Vec<LegendEntry>,
}

fn surface_style(points: &[DesignPoint], j: usize, geometry_only: bool, config: &FoilPlotConfig) -> LineStyle {
    let rgb = if geometry_only && points.len() > 1 {
        config.palette[j % config.palette.len()]
    } else {
        match config.panels {
            PanelStyle::Notches => SURFACE_NOTCHED,
            PanelStyle::Dots => SURFACE_DOTTED,
            PanelStyle::None => SURFACE_PLAIN,
        }
    };
    LineStyle {
        rgb,
        width: 1.5,
        dasharray: None,
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
        legend: Vec::new(),
        marker_legend: Vec::new(),
    };

    // Airfoil outlines: one per geometry-only point; identical panels otherwise, so once
    let notch = config.notch_length * chord;
    let surface_points: Vec<usize> = if geometry_only {
        (0..n_points).collect()
    } else {
        vec![0]
    };
    for &j in &surface_points {
        let g = &points[j].geometry;
        let style = surface_style(points, j, geometry_only, config);
        let outline: Vec<Pt> = g.x.iter().copied().zip(g.y.iter().copied()).collect();
        match config.panels {
            PanelStyle::Notches => {
                let segs: Vec<Seg> = (0..g.n())
                    .map(|i| {
                        (
                            (g.x[i], g.y[i]),
                            (g.x[i] + notch * g.normal_x[i], g.y[i] + notch * g.normal_y[i]),
                        )
                    })
                    .collect();
                lay.notches.push(segs);
            }
            PanelStyle::Dots => lay.dots.push(outline.clone()),
            PanelStyle::None => {}
        }
        if geometry_only && n_points > 1 {
            lay.legend.push(LegendEntry {
                text: points[j].display_label(),
                key: LegendKey::Line(style.clone()),
            });
        }
        lay.surfaces.push(Surface { outline, style });
    }

    // Wakes: XYWAKE rebuilds the wake for every alpha, so each design point carries its own;
    // black, with the point's dash pattern. Wake panels are shown only for a single point.
    if config.show_wake {
        for (j, p) in points.iter().enumerate() {
            let Some(w) = &p.geometry.wake else { continue };
            let line: Vec<Pt> = w.x.iter().copied().zip(w.y.iter().copied()).collect();
            let pattern = config.patterns[j % config.patterns.len()];
            let style = LineStyle {
                rgb: WAKE_LINE,
                width: 1.5,
                dasharray: pattern.dasharray().map(str::to_string),
            };
            if n_points == 1 {
                match config.panels {
                    PanelStyle::Notches => {
                        let segs: Vec<Seg> = (0..w.x.len())
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
                            .collect();
                        lay.notches.push(segs);
                    }
                    PanelStyle::Dots => lay.dots.push(line.clone()),
                    PanelStyle::None => {}
                }
            }
            lay.wakes.push(Curve { pts: line, style });
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
            let style = LineStyle {
                rgb: st.rgb,
                width: 1.5,
                dasharray: st.dasharray.clone(),
            };
            let upper = surface_offset(g, &bl.upper, bl.upper.column(q), k);
            let lower = surface_offset(g, &bl.lower, bl.lower.column(q), k);
            if config.show_wake && g.wake.is_some() && !bl.wake.is_empty() {
                let split = if q == BlQuantity::Dstar {
                    bl.wake_split
                } else {
                    [0.5, 0.5]
                };
                let te_u = *upper.last().unwrap_or(&(g.x_te, g.y_te));
                let te_l = *lower.last().unwrap_or(&(g.x_te, g.y_te));
                let (wu, wl) = wake_band(g, &bl.wake, bl.wake.column(q), k, split, te_u, te_l);
                lay.curves.push(Curve {
                    pts: wu,
                    style: style.clone(),
                });
                lay.curves.push(Curve {
                    pts: wl,
                    style: style.clone(),
                });
            }
            lay.curves.push(Curve {
                pts: upper,
                style: style.clone(),
            });
            lay.curves.push(Curve {
                pts: lower,
                style: style.clone(),
            });
            let text = if n_points > 1 {
                format!("{} at {}, scaled ×{}", q.label(), p.display_label(), fmt_scale(k))
            } else {
                format!("{}, scaled ×{}", q.label(), fmt_scale(k))
            };
            lay.legend.push(LegendEntry {
                text,
                key: LegendKey::Line(style),
            });
        }
    }
    if config.quantities.is_empty() && !geometry_only && n_points > 1 {
        // several design points with no quantities: the wake dash patterns tell them apart
        for (j, p) in points.iter().enumerate() {
            let pattern = config.patterns[j % config.patterns.len()];
            lay.legend.push(LegendEntry {
                text: p.display_label(),
                key: LegendKey::Line(LineStyle {
                    rgb: WAKE_LINE,
                    width: 1.5,
                    dasharray: pattern.dasharray().map(str::to_string),
                }),
            });
        }
    }

    // Markers: open shape per kind, colour per side (toned per design point)
    if config.markers.any() {
        let side_colour = |is: usize| if is == 1 { MARK_UPPER } else { MARK_LOWER };
        for (j, p) in points.iter().enumerate() {
            let Some(bl) = &p.boundary_layer else { continue };
            let f = tone_factor(j, n_points, config.tone_range);
            let mut mark = |x: f64, y: f64, kind: MarkKind, rgb: (u8, u8, u8)| {
                lay.markers.push(Mark {
                    at: (x, y),
                    kind,
                    rgb: tone(rgb, f),
                });
            };
            if config.markers.stagnation {
                let m = &bl.stagnation;
                mark(m.x, m.y, MarkKind::Stagnation, MARK_STAGNATION);
            }
            if config.markers.transition {
                for (is, m) in [1, 2].into_iter().zip(&bl.transition) {
                    let kind = if m.forced {
                        MarkKind::Forced
                    } else {
                        MarkKind::Transition
                    };
                    mark(m.x_transition, m.y_transition, kind, side_colour(is));
                }
            }
            if config.markers.separation {
                for m in &bl.derived_separation {
                    let kind = match m.kind {
                        SeparationKind::Separation => MarkKind::Separation,
                        SeparationKind::Reattachment => MarkKind::Reattachment,
                    };
                    mark(m.x, m.y, kind, side_colour(m.side));
                }
            }
        }
        for kind in MarkKind::ALL {
            if lay.markers.iter().any(|m| m.kind == kind) {
                lay.marker_legend.push(LegendEntry {
                    text: kind.label().to_string(),
                    key: LegendKey::Marker(kind.shape(), MARK_KEY),
                });
            }
        }
        if lay.markers.iter().any(|m| m.kind != MarkKind::Stagnation) {
            for (text, rgb) in [("upper surface", MARK_UPPER), ("lower surface", MARK_LOWER)] {
                lay.marker_legend.push(LegendEntry {
                    text: text.to_string(),
                    key: LegendKey::Line(LineStyle {
                        rgb,
                        width: 2.0,
                        dasharray: None,
                    }),
                });
            }
        }
    }

    Ok(lay)
}

/// Render the plot to an SVG document
pub fn render_foil_svg(points: &[DesignPoint], config: &FoilPlotConfig) -> Result<String, PlotError> {
    let lay = layout(points, config)?;
    let chord = points[0].geometry.chord;

    // Bounds: everything drawn, padded by 5% chord
    let mut bx = (f64::INFINITY, f64::NEG_INFINITY);
    let mut by = (f64::INFINITY, f64::NEG_INFINITY);
    let mut take = |x: f64, y: f64| {
        if x.is_finite() && y.is_finite() {
            bx = (bx.0.min(x), bx.1.max(x));
            by = (by.0.min(y), by.1.max(y));
        }
    };
    for sf in &lay.surfaces {
        sf.outline.iter().for_each(|&(x, y)| take(x, y));
    }
    for w in &lay.wakes {
        w.pts.iter().for_each(|&(x, y)| take(x, y));
    }
    for segs in &lay.notches {
        segs.iter().for_each(|&(a, b)| {
            take(a.0, a.1);
            take(b.0, b.1);
        });
    }
    for c in &lay.curves {
        c.pts.iter().for_each(|&(x, y)| take(x, y));
    }
    for m in &lay.markers {
        take(m.at.0, m.at.1);
    }
    let pad = 0.05 * chord;
    let bounds = (bx.0 - pad, bx.1 + pad, by.0 - pad, by.1 + pad);

    let margins = Margins {
        left: 60.0,
        right: 20.0,
        top: 40.0,
        bottom: 40.0,
    };
    let mut c = Canvas::equal_aspect(config.width, config.height, margins, bounds);
    c.header(config.background);
    let title = config.title.clone().unwrap_or_else(|| {
        if points.len() == 1 {
            points[0].display_label()
        } else {
            "Foil".to_string()
        }
    });
    c.title(&title);
    c.grid("x/c", "y/c");

    for sf in &lay.surfaces {
        c.polyline(&sf.outline, &sf.style);
    }
    for w in &lay.wakes {
        c.polyline(&w.pts, &w.style);
    }
    for segs in &lay.notches {
        c.group_start(&LineStyle {
            rgb: SURFACE_NOTCHED,
            width: 1.0,
            dasharray: None,
        });
        for &(a, b) in segs {
            c.group_line(a, b);
        }
        c.group_end();
    }
    for d in &lay.dots {
        for &p in d {
            c.circle(p, 2.0, (0, 0, 0));
        }
    }
    for curve in &lay.curves {
        c.polyline(&curve.pts, &curve.style);
    }
    for m in &lay.markers {
        c.marker(m.at, m.kind.shape(), MARK_R, m.rgb, 1.5);
    }
    c.legend(&lay.legend, LegendAnchor::BottomLeft);
    c.legend(&lay.marker_legend, LegendAnchor::BottomRight);
    Ok(c.finish())
}

/// Render and write the plot
pub fn plot_foil<P: AsRef<Path>>(
    points: &[DesignPoint],
    path: P,
    format: ImageFormat,
    config: &FoilPlotConfig,
) -> Result<(), PlotError> {
    let svg = render_foil_svg(points, config)?;
    match format {
        ImageFormat::Svg => write_svg(path.as_ref(), &svg),
        ImageFormat::Png => write_png(path.as_ref(), &svg),
    }
}

/// [`plot_foil`] as SVG
pub fn plot_foil_svg<P: AsRef<Path>>(
    points: &[DesignPoint],
    path: P,
    config: &FoilPlotConfig,
) -> Result<(), PlotError> {
    plot_foil(points, path, ImageFormat::Svg, config)
}

/// [`plot_foil`] as PNG
pub fn plot_foil_png<P: AsRef<Path>>(
    points: &[DesignPoint],
    path: P,
    config: &FoilPlotConfig,
) -> Result<(), PlotError> {
    plot_foil(points, path, ImageFormat::Png, config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{naca_4digit, panel_foil, Thickness};
    use crate::output::{SideStations, WakeNodes};

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
        assert_eq!(a.dasharray, b.dasharray, "same design point: same pattern");
        let a2 = point_style(0, 2, 3, &c);
        assert_eq!(a2.dasharray, Some("8,4".to_string()));
        assert_eq!(a2.rgb, tone(c.palette[0], 0.55), "last point gets the dark tone");
        let s0 = point_style(0, 0, 1, &c);
        assert_eq!(s0.rgb, c.palette[0], "single point: base hue");
        assert_eq!(s0.dasharray, None);
        assert_eq!(point_style(0, 5, 6, &c).dasharray, None, "cycle wraps after five");
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
        let texts: Vec<&str> = lay.legend.iter().map(|e| e.text.as_str()).collect();
        assert_eq!(texts.len(), 4);
        assert!(texts[0].starts_with("δ* at a, scaled ×"));
        assert!(texts[1].starts_with("δ* at b (not converged), scaled ×"));
        assert!(texts[2].starts_with("θ at a, scaled ×"));
        // four distinct (colour, dash) combinations
        let mut styles: Vec<String> = lay.legend.iter().map(|e| format!("{:?}", e.key)).collect();
        styles.sort();
        styles.dedup();
        assert_eq!(styles.len(), 4);
        // 2 points × 2 quantities × (upper, lower, wake upper, wake lower)
        assert_eq!(lay.curves.len(), 16);
        // airfoil notches once; wake panels only for a single point; one wake line per point,
        // black, dash-coded like its curves
        assert_eq!(lay.notches.len(), 1);
        assert_eq!(lay.notches[0].len(), 60);
        assert_eq!(lay.wakes.len(), 2);
        assert!(lay.wakes.iter().all(|w| w.style.rgb == WAKE_LINE));
        assert_eq!(lay.wakes[0].style.dasharray, None);
        assert_eq!(lay.wakes[1].style.dasharray, Some("2,3".to_string()));
        // single point: no label prefix, and the wake gets its panels
        let one = layout(&pts[..1], &config).unwrap();
        assert!(one.legend[0].text.starts_with("δ*, scaled ×"));
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
        assert_ne!(lay.surfaces[0].style.rgb, lay.surfaces[1].style.rgb);
        assert!(matches!(&lay.legend[0].key, LegendKey::Line(s) if s.rgb == lay.surfaces[0].style.rgb));
        assert_eq!(lay.legend.len(), 2);
        // mismatched panels with quantities
        let mut b = point("b", 0.01, true);
        b.geometry.x[3] = f64::from_bits(b.geometry.x[3].to_bits() + 1);
        let err = layout(&[point("a", 0.01, true), b], &config).unwrap_err();
        assert!(err.to_string().contains("different panel geometry"));
    }

    #[test]
    fn renders_svg_with_markers() {
        let pts = [point("a", 0.01, true)];
        let config = FoilPlotConfig {
            quantities: vec![BlQuantity::Dstar],
            panels: PanelStyle::Dots,
            ..Default::default()
        };
        let svg = render_foil_svg(&pts, &config).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.trim_end().ends_with("</svg>"));
        // stagnation circle (dark), upper transition square (yellow), forced lower transition
        // diamond (green); the marker legend names the three kinds and the two sides
        assert_eq!(svg.matches("fill=\"none\" stroke=\"rgb(40,40,40)\"").count(), 1);
        assert_eq!(
            svg.matches("<rect x=").count(),
            5,
            "square marker + its legend key + border + 2 legend boxes"
        );
        assert!(svg.contains("<polygon points="));
        assert_eq!(svg.matches("stroke=\"rgb(235,175,0)\"").count(), 2, "square + side key");
        assert_eq!(svg.matches("stroke=\"rgb(0,165,0)\"").count(), 2, "diamond + side key");
        for text in [
            "stagnation",
            "transition",
            "forced transition",
            "upper surface",
            "lower surface",
        ] {
            assert!(svg.contains(&format!(">{text}<")), "{text} in the marker legend");
        }
        assert!(!svg.contains(">separation<"));
        // 60 node dots + stagnation circle + the stagnation legend key
        assert_eq!(svg.matches("<circle").count(), 62);
        assert_eq!(ImageFormat::from_path(Path::new("x.PNG")), ImageFormat::Png);
        assert_eq!(fmt_alpha(3.0000000000000004), "3");
        assert_eq!(fmt_alpha(-2.5), "-2.5");
        assert_eq!(fmt_alpha(-0.0), "0");
        assert_eq!(fmt_alpha(12.000000000000002), "12");
        assert_eq!(ImageFormat::from_path(Path::new("x")), ImageFormat::Svg);
    }
}
