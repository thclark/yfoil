//! The `yfoil plot analysis` and `yfoil plot polar` figures, and what the three `yfoil plot`
//! figures share: their error type and the parameters and labels the ironlab viewer sorts,
//! filters and groups by (`docs/development/conventions/figures.md`). The foil figure is in `foil_plot.rs`.
//!
//! Every figure is an ironlab figure, written as `.fig` (or `.fig.json`) through
//! [`figures::save`], which refuses a figure with validation errors.

use std::collections::BTreeMap;
use std::path::Path;

use ironlab::prelude::*;

use crate::output::figures::{self, FigureError};
use crate::output::{AnalysisOutput, BoundaryLayerOutput, PointRecord, PolarOutput, SeparationKind};
use crate::solver::analysis::FlowConditions;

/// The `study` parameter of every `yfoil plot` figure
pub const STUDY: &str = "yfoil-plot";

/// Error type for plotting operations
#[derive(thiserror::Error, Debug)]
pub enum PlotError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Config(String),
    #[error(transparent)]
    Figure(#[from] FigureError),
}

/// Validate and write a figure where `-o` names it: `.fig` or `.fig.json`
pub fn save_figure(fig: &Figure, path: &Path) -> Result<(), PlotError> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    let lower = name.to_ascii_lowercase();
    if !(lower.ends_with(".fig") || lower.ends_with(".fig.json")) {
        return Err(PlotError::Config(format!(
            "{}: a figure is written as .fig or .fig.json",
            path.display()
        )));
    }
    figures::save(fig, path)?;
    Ok(())
}

// ============================================================================
// Parameters and labels
// ============================================================================

/// The flow conditions as figure parameters (`figures::flow_condition_parameters`), as a map
pub(crate) fn condition_parameters(conditions: &FlowConditions) -> BTreeMap<String, Parameter> {
    figures::flow_condition_parameters(Figure::new(), conditions)
        .parameters()
        .clone()
}

/// The regime labels of a set of flow conditions (`docs/development/conventions/figures.md`, "Labels")
pub(crate) fn regime_labels(conditions: &FlowConditions) -> Vec<&'static str> {
    let mut labels = Vec::new();
    match conditions.re {
        None => labels.push("inviscid"),
        Some(re) if re < 1e5 => labels.push("low-re"),
        Some(_) => {}
    }
    if conditions.mach > 0.0 {
        labels.push("compressible");
    }
    if conditions.re.is_some() && conditions.x_trip.iter().any(|&x| x < 1.0) {
        labels.push("forced-transition");
    }
    labels
}

/// The parameters of one operating point: its inputs (`alpha_deg`) and outputs (forces,
/// transition, separation and TE thicknesses), each left out where it does not exist
pub(crate) fn point_parameters(
    record: &PointRecord,
    boundary_layer: Option<&BoundaryLayerOutput>,
) -> BTreeMap<String, Parameter> {
    let mut fig = Figure::new()
        .parameter("alpha_deg", record.alpha_deg)
        .parameter("converged", record.is_valid());
    if let Some(v) = &record.values {
        fig = figures::optional_parameter(fig, "cl", Some(v.cl));
        fig = figures::optional_parameter(fig, "cd", v.cd);
        fig = figures::optional_parameter(fig, "cm", Some(v.cm));
    }
    if let Some(bl) = boundary_layer {
        for (side, name) in [(1, "upper"), (2, "lower")] {
            let stations = if side == 1 { &bl.upper } else { &bl.lower };
            fig = figures::optional_parameter(
                fig,
                &format!("x_transition_{name}"),
                Some(bl.transition[side - 1].x_transition),
            );
            let separation = bl
                .derived_separation
                .iter()
                .find(|m| m.side == side && m.kind == SeparationKind::Separation)
                .map(|m| m.x);
            fig = figures::optional_parameter(fig, &format!("x_separation_{name}"), separation);
            fig = figures::optional_parameter(
                fig,
                &format!("dstar_te_{name}"),
                stations.primaries.dstar.last().copied(),
            );
            fig = figures::optional_parameter(
                fig,
                &format!("theta_te_{name}"),
                stations.primaries.theta.last().copied(),
            );
        }
    }
    fig.parameters().clone()
}

/// The entries every map holds with the same value
pub(crate) fn shared_parameters(sets: &[BTreeMap<String, Parameter>]) -> BTreeMap<String, Parameter> {
    let Some((first, rest)) = sets.split_first() else {
        return BTreeMap::new();
    };
    first
        .iter()
        .filter(|(k, v)| rest.iter().all(|s| s.get(*k) == Some(*v)))
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect()
}

/// The labels every list holds
pub(crate) fn shared_labels(sets: &[Vec<&'static str>]) -> Vec<&'static str> {
    let Some((first, rest)) = sets.split_first() else {
        return Vec::new();
    };
    first
        .iter()
        .copied()
        .filter(|l| rest.iter().all(|s| s.contains(l)))
        .collect()
}

/// Set a figure's parameters and labels
pub(crate) fn describe(mut fig: Figure, parameters: BTreeMap<String, Parameter>, labels: &[&str]) -> Figure {
    for (name, value) in parameters {
        fig = fig.parameter(name, value);
    }
    for label in labels {
        fig = fig.label(*label);
    }
    fig
}

/// An RGB triple as an ironlab colour
pub(crate) fn rgb(c: (u8, u8, u8)) -> Color {
    Color::rgb(c.0 as f32 / 255.0, c.1 as f32 / 255.0, c.2 as f32 / 255.0)
}

/// Padded [min, max] of the finite values; `floor` guards a degenerate range
fn padded_range<'a>(values: impl Iterator<Item = &'a f64>, pad: f64, floor: f64) -> (f64, f64) {
    let (lo, hi) = values
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    if !lo.is_finite() || !hi.is_finite() {
        return (0.0, 1.0);
    }
    let span = (hi - lo).max(floor);
    (lo - pad * span, hi + pad * span)
}

// ============================================================================
// Analysis: Cp and Ue over x/c
// ============================================================================

/// Configuration of the analysis figure
#[derive(Debug, Clone)]
pub struct AnalysisPlotConfig {
    /// Upper surface colour (RGB), blue
    pub upper_color: (u8, u8, u8),
    /// Lower surface colour (RGB), red
    pub lower_color: (u8, u8, u8),
    /// Aerofoil outline colour (RGB), grey
    pub airfoil_color: (u8, u8, u8),
    /// Figure title (default: "<foil> at α = …°")
    pub title: Option<String>,
}

impl Default for AnalysisPlotConfig {
    fn default() -> Self {
        Self {
            upper_color: (0, 100, 200),
            lower_color: (200, 50, 50),
            airfoil_color: (100, 100, 100),
            title: None,
        }
    }
}

/// Height of the analysis figure (mm): two distribution panels over a strip with the aerofoil
const ANALYSIS_HEIGHT_MM: f64 = 120.0;

/// The analysis figure: −Cp (top) and Ue/U∞ (middle) over x/c, upper surface blue and lower
/// red, over a strip with the aerofoil outline; x is linked across the three.
///
/// Cp is drawn as −Cp because ironlab has no reversed axis (the aerodynamic convention puts
/// suction up).
pub fn analysis_figure(analysis: &AnalysisOutput, config: &AnalysisPlotConfig) -> Result<Figure, PlotError> {
    let (Some((x_upper, cp_upper, q_upper)), Some((x_lower, cp_lower, q_lower))) =
        (analysis.upper_surface(), analysis.lower_surface())
    else {
        return Err(PlotError::Config(format!(
            "{}: the analysis carries no surface distributions (its point is not valid; `yfoil analyse --allow-invalid` keeps them)",
            analysis.foil
        )));
    };
    let minus = |v: &[f64]| v.iter().map(|c| -c).collect::<Vec<f64>>();
    let (minus_cp_upper, minus_cp_lower) = (minus(&cp_upper), minus(&cp_lower));

    let g = &analysis.geometry;
    let chord = g.chord;
    let (x_lo, x_hi) = padded_range(g.x.iter(), 0.0, 0.0);
    let xlim = (x_lo - 0.05 * chord, x_hi + 0.1 * chord);
    let cp_lim = padded_range(minus_cp_upper.iter().chain(&minus_cp_lower), 0.1, 1e-3);
    let q_lim = padded_range(q_upper.iter().chain(&q_lower), 0.1, 1e-3);
    let y_lim = padded_range(g.y.iter(), 0.1, 1e-3);

    let title = config
        .title
        .clone()
        .unwrap_or_else(|| format!("{} at α = {:.1}°", analysis.foil, analysis.results.alpha_deg));
    let mut parameters = BTreeMap::from([
        ("study".to_string(), Parameter::from(STUDY)),
        ("figure".to_string(), Parameter::from("analysis")),
        ("foil".to_string(), Parameter::from(analysis.foil.as_str())),
        ("n_nodes".to_string(), Parameter::from(g.n() as i64)),
    ]);
    parameters.extend(condition_parameters(&analysis.conditions));
    parameters.extend(point_parameters(&analysis.results, analysis.boundary_layer.as_ref()));
    let mut fig = Figure::new()
        .size_mm(figures::WIDTH_MM, ANALYSIS_HEIGHT_MM)
        .tiles(3, 1)
        .row_weights([3.0, 3.0, 1.0])
        .title(Text::plain(title));
    fig = describe(fig, parameters, &regime_labels(&analysis.conditions));

    let (upper, lower) = (rgb(config.upper_color), rgb(config.lower_color));
    let mut ax = fig.axes(0, 0);
    ax.plot(&x_upper, &minus_cp_upper)
        .display_name("Upper surface")
        .color(upper)
        .line_width(1.5);
    ax.plot(&x_lower, &minus_cp_lower)
        .display_name("Lower surface")
        .color(lower)
        .line_width(1.5);
    ax.ylabel("$-C_p$")
        .xlim(xlim.0, xlim.1)
        .ylim(cp_lim.0, cp_lim.1)
        .grid(true)
        .legend(LegendLocation::NorthEast);
    let cp_axes = ax.id();

    let mut ax = fig.axes(1, 0);
    ax.plot(&x_upper, &q_upper).color(upper).line_width(1.5);
    ax.plot(&x_lower, &q_lower).color(lower).line_width(1.5);
    ax.ylabel("$U_e/U_\\infty$")
        .xlim(xlim.0, xlim.1)
        .ylim(q_lim.0, q_lim.1)
        .grid(true);
    let ue_axes = ax.id();

    let mut ax = fig.axes(2, 0);
    ax.plot(&g.x, &g.y).color(rgb(config.airfoil_color));
    ax.xlabel("$x/c$")
        .ylabel("$y/c$")
        .xlim(xlim.0, xlim.1)
        .ylim(y_lim.0, y_lim.1);
    let foil_axes = ax.id();

    fig.link(Dim::X, &[cp_axes, ue_axes, foil_axes])
        .map_err(|e| PlotError::Config(format!("linking the analysis axes: {e}")))?;
    Ok(fig)
}

/// Build and write the analysis figure
pub fn plot_analysis<P: AsRef<Path>>(
    analysis: &AnalysisOutput,
    path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    save_figure(&analysis_figure(analysis, config)?, path.as_ref())
}

// ============================================================================
// Polars: CL–α, CL–CD, CM–α, CD–α
// ============================================================================

/// The points of a polar that can be drawn: valid, and therefore carrying their values.
/// An invalid or never-attempted point has nothing to plot and is left out of every series, so a
/// Kármán–Tsien artefact cannot appear as a bump on a curve.
fn valid_points(polar: &PolarOutput) -> Vec<(f64, &crate::output::PointValues)> {
    polar
        .results
        .iter()
        .filter(|p| p.is_valid())
        .filter_map(|p| p.values.as_ref().map(|v| (p.alpha_deg, v)))
        .collect()
}

/// One polar curve with its legend label and what the figure's parameters are made from
#[derive(Debug, Clone)]
pub struct PolarSeries {
    /// Legend label
    pub label: String,
    /// Angle of attack values (degrees)
    pub alpha: Vec<f64>,
    /// Lift coefficient values
    pub cl: Vec<f64>,
    /// Drag coefficient values (NaN for an inviscid point, which has none)
    pub cd: Vec<f64>,
    /// Moment coefficient values (about quarter chord)
    pub cm: Vec<f64>,
    /// The aerofoil name
    pub foil: String,
    /// The polar's flow conditions
    pub conditions: FlowConditions,
    /// The polar's summary: the largest CL over the valid points and its alpha
    pub cl_max: Option<f64>,
    pub alpha_at_cl_max: Option<f64>,
}

impl PolarSeries {
    /// Build from a `PolarOutput`, keeping only the valid points.
    ///
    /// The label is the polar's `label` if set, otherwise its `foil` name.
    pub fn from_polar_output(polar: &PolarOutput) -> Self {
        let pts = valid_points(polar);
        Self {
            label: polar.label.clone().unwrap_or_else(|| polar.foil.clone()),
            alpha: pts.iter().map(|(a, _)| *a).collect(),
            cl: pts.iter().map(|(_, v)| v.cl).collect(),
            cd: pts.iter().map(|(_, v)| v.cd.unwrap_or(f64::NAN)).collect(),
            cm: pts.iter().map(|(_, v)| v.cm).collect(),
            foil: polar.foil.clone(),
            conditions: polar.conditions.clone(),
            cl_max: polar.summary.cl_max,
            alpha_at_cl_max: polar.summary.alpha_at_cl_max,
        }
    }
}

/// Point-marker style, cycled per series. Open markers come first so that coincident series
/// stay visible through each other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeriesMarker {
    pub shape: Marker,
    pub filled: bool,
}

impl SeriesMarker {
    const fn open(shape: Marker) -> Self {
        Self { shape, filled: false }
    }
    const fn filled(shape: Marker) -> Self {
        Self { shape, filled: true }
    }

    /// Default cycle order: every open marker (`x + o □ △ ◇ ▽`), then the filled variants
    pub const CYCLE: [SeriesMarker; 12] = [
        SeriesMarker::open(Marker::Cross),
        SeriesMarker::open(Marker::Plus),
        SeriesMarker::open(Marker::Circle),
        SeriesMarker::open(Marker::Square),
        SeriesMarker::open(Marker::TriangleUp),
        SeriesMarker::open(Marker::Diamond),
        SeriesMarker::open(Marker::TriangleDown),
        SeriesMarker::filled(Marker::Circle),
        SeriesMarker::filled(Marker::Square),
        SeriesMarker::filled(Marker::TriangleUp),
        SeriesMarker::filled(Marker::Diamond),
        SeriesMarker::filled(Marker::TriangleDown),
    ];
}

/// Configuration of the polar figure
#[derive(Debug, Clone)]
pub struct PolarsPlotConfig {
    /// Title; defaults to the single series' label, or "Polar comparison"
    pub title: Option<String>,
    /// Series colours (RGB), cycled when there are more series than colours
    pub palette: Vec<(u8, u8, u8)>,
    /// Series markers, cycled independently of the palette
    pub markers: Vec<SeriesMarker>,
}

impl Default for PolarsPlotConfig {
    fn default() -> Self {
        Self {
            title: None,
            palette: vec![
                (0, 100, 200),
                (200, 50, 50),
                (30, 150, 60),
                (220, 130, 0),
                (120, 60, 180),
                (0, 150, 160),
            ],
            markers: SeriesMarker::CYCLE.to_vec(),
        }
    }
}

/// Height of the polar figure (mm): a 2×2 sheet over a legend row
const POLAR_HEIGHT_MM: f64 = 140.0;

/// The polar figure: CL–α, CL–CD (drag polar), CM–α and CD–α in a 2×2 grid over one legend.
/// α is linked across the three α panels, and CL between the two CL panels.
pub fn polars_figure(series: &[PolarSeries], config: &PolarsPlotConfig) -> Result<Figure, PlotError> {
    if series.is_empty() {
        return Err(PlotError::Config("no polars to plot".to_string()));
    }
    let title = config.title.clone().unwrap_or_else(|| {
        if series.len() == 1 {
            series[0].label.clone()
        } else {
            "Polar comparison".to_string()
        }
    });

    // shared axis ranges across every series so the panels line up
    let alpha_r = padded_range(series.iter().flat_map(|s| s.alpha.iter()), 0.05, 1.0);
    let cl_r = padded_range(series.iter().flat_map(|s| s.cl.iter()), 0.08, 0.1);
    let cd_r = padded_range(series.iter().flat_map(|s| s.cd.iter()), 0.08, 0.001);
    let cm_r = padded_range(series.iter().flat_map(|s| s.cm.iter()), 0.08, 0.01);

    let mut parameters = BTreeMap::from([
        ("study".to_string(), Parameter::from(STUDY)),
        ("figure".to_string(), Parameter::from("polar")),
    ]);
    let per_series: Vec<BTreeMap<String, Parameter>> = series
        .iter()
        .map(|s| {
            let mut p = condition_parameters(&s.conditions);
            p.insert("foil".to_string(), Parameter::from(s.foil.as_str()));
            p
        })
        .collect();
    parameters.extend(shared_parameters(&per_series));
    if let [s] = series {
        parameters.insert("points".to_string(), Parameter::from(s.alpha.len() as i64));
        for (name, value) in [("cl_max", s.cl_max), ("alpha_cl_max_deg", s.alpha_at_cl_max)] {
            if let Some(v) = value.filter(|v| v.is_finite()) {
                parameters.insert(name.to_string(), Parameter::from(v));
            }
        }
    }
    let labels = shared_labels(&series.iter().map(|s| regime_labels(&s.conditions)).collect::<Vec<_>>());

    // the key: up to four series to a row, its row sized to the rows of entries
    let key_columns = series.len().min(4);
    let key_rows = series.len().div_ceil(4);
    let mut fig = Figure::new()
        .size_mm(figures::WIDTH_MM, POLAR_HEIGHT_MM)
        .tiles(3, 2)
        .row_weights([3.0, 3.0, 0.4 + 0.5 * key_rows as f64])
        .title(Text::plain(title));
    fig = describe(fig, parameters, &labels);
    let key = fig
        .panel_legend_span(2, 0, 1, 2)
        .columns(LegendEntryColumns::Count(key_columns as u32))
        .order(LegendEntryOrder::RowMajor)
        .id();

    // (row, column, x label, y label, x limits, y limits, the series' x and y)
    type Extract = fn(&PolarSeries) -> (&[f64], &[f64]);
    type Panel = (u32, u32, &'static str, &'static str, (f64, f64), (f64, f64), Extract);
    let panels: [Panel; 4] = [
        (0, 0, "$\\alpha$ (deg)", "$C_L$", alpha_r, cl_r, |s| (&s.alpha, &s.cl)),
        (0, 1, "$C_D$", "$C_L$", cd_r, cl_r, |s| (&s.cd, &s.cl)),
        (1, 0, "$\\alpha$ (deg)", "$C_M$", alpha_r, cm_r, |s| (&s.alpha, &s.cm)),
        (1, 1, "$\\alpha$ (deg)", "$C_D$", alpha_r, cd_r, |s| (&s.alpha, &s.cd)),
    ];
    let mut ids = Vec::with_capacity(panels.len());
    for (row, col, xl, yl, xr, yr, extract) in panels {
        let mut ax = fig.axes(row, col);
        for (i, s) in series.iter().enumerate() {
            let color = rgb(config.palette[i % config.palette.len()]);
            let marker = config.markers[i % config.markers.len()];
            let (x, y) = extract(s);
            let mut line = ax.plot(x, y);
            line.display_name(Text::plain(s.label.clone()))
                .color(color)
                .line_width(1.0)
                .marker(marker.shape)
                .marker_size(4.0);
            if marker.filled {
                line.marker_face(color);
            }
        }
        ax.xlabel(xl)
            .ylabel(yl)
            .xlim(xr.0, xr.1)
            .ylim(yr.0, yr.1)
            .grid(true)
            .legend_in(key);
        ids.push(ax.id());
    }
    let link = |e: ironlab::Error| PlotError::Config(format!("linking the polar axes: {e}"));
    fig.link(Dim::X, &[ids[0], ids[2], ids[3]]).map_err(link)?;
    fig.link(Dim::Y, &[ids[0], ids[1]]).map_err(link)?;
    Ok(fig)
}

/// Build and write the polar figure
pub fn plot_polars<P: AsRef<Path>>(
    series: &[PolarSeries],
    path: P,
    config: &PolarsPlotConfig,
) -> Result<(), PlotError> {
    save_figure(&polars_figure(series, config)?, path.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padded_range_ignores_non_finite_values() {
        assert_eq!(padded_range([1.0, f64::NAN, 3.0].iter(), 0.5, 0.0), (0.0, 4.0));
        assert_eq!(padded_range([f64::NAN].iter(), 0.5, 0.0), (0.0, 1.0));
        assert_eq!(padded_range([2.0].iter(), 0.5, 1.0), (1.5, 2.5));
    }

    #[test]
    fn shared_parameters_keep_what_every_set_agrees_on() {
        let a = BTreeMap::from([
            ("mach".to_string(), Parameter::from(0.0)),
            ("re".to_string(), Parameter::from(1e6)),
        ]);
        let b = BTreeMap::from([
            ("mach".to_string(), Parameter::from(0.0)),
            ("re".to_string(), Parameter::from(3e6)),
        ]);
        let shared = shared_parameters(&[a.clone(), b]);
        assert_eq!(shared.len(), 1);
        assert!(shared.contains_key("mach"));
        assert_eq!(shared_parameters(std::slice::from_ref(&a)), a);
        assert_eq!(
            shared_labels(&[vec!["inviscid", "compressible"], vec!["compressible"]]),
            vec!["compressible"]
        );
    }

    #[test]
    fn save_refuses_other_extensions() {
        let dir = tempfile::tempdir().unwrap();
        let fig = Figure::new();
        let err = save_figure(&fig, &dir.path().join("a.svg")).unwrap_err();
        assert!(err.to_string().contains(".fig or .fig.json"));
    }
}
