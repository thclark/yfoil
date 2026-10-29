//! The study's figures, from a run's `summary.json`: the spread histogram, the threshold
//! plateau and the spread against α (`docs/development/conventions/figures.md`). The series
//! cases' points are blue circles and the branch-coverage cases' green squares ([`group_style`]);
//! the error threshold grey and dashed, the warning threshold grey and dash-dotted; a spread that
//! is not a number orange. Each figure's key is a panel legend in a tile of its own under the
//! axes, and each figure names its groups in the `groups` parameter.

use crate::analysis::{
    flagged, gap, groups, histogram, measured, plateau_thresholds, points, points_of, thresholds, worst_spreads,
    worst_spreads_of, DECADES,
};
use crate::drive::val;
use crate::report::named_thresholds;
use ironlab::prelude::*;
use serde_json::Value;
use std::path::Path;
use yfoil::output::figures::{self, BLUE, CATEGORY, GREY, ORANGE, WIDTH_MM};

pub const NAMES: [&str; 3] = ["spread-histogram", "threshold-plateau", "spread-against-alpha"];
const AXES_HEIGHT_MM: f64 = 80.0;
/// The legend tile's height relative to the axes'.
const LEGEND_WEIGHT: f64 = 0.5;
/// Where the histogram draws its not-a-number bin: two decades past the last decade.
const NAN_AT: f64 = 1.0e2;

fn figure(name: &str, s: &Value) -> Figure {
    Figure::new()
        .size_mm(WIDTH_MM, AXES_HEIGHT_MM * (1.0 + LEGEND_WEIGHT))
        .tiles(2, 1)
        .fit_row(1)
        .parameter("study", "conditioning")
        .parameter("figure", name)
        .parameter("warning_threshold", thresholds(s).warning)
        .parameter("error_threshold", thresholds(s).error)
        .parameter("cases", s["cases"].as_array().map_or(0, |a| a.len() as i64))
        .parameter("points", points(s).count() as i64)
        .parameter("groups", groups(s).join(","))
        .tag("validation")
}

/// The marker and colour that tell the groups apart: blue circles for the series cases, green
/// squares for the branch-coverage cases (and purple diamonds for any other).
fn group_style(group: &str) -> (Marker, Color) {
    match group {
        "series" => (Marker::Circle, BLUE),
        "branch-coverage" => (Marker::Square, CATEGORY[0]),
        _ => (Marker::Diamond, CATEGORY[1]),
    }
}

/// The colour of a curve over every group together.
const ALL: Color = Color::rgb(0.0, 0.0, 0.0);

/// The figure's key, a panel legend under the axes in `columns` columns.
fn key(fig: &mut Figure, columns: u32) -> NodeId {
    let mut panel = fig.panel(1, 0);
    let mut legend = panel.legend();
    legend
        .columns(LegendEntryColumns::Count(columns))
        .order(LegendEntryOrder::ColumnMajor);
    legend.id()
}

/// Both thresholds as reference lines: vertical across `y` when `vertical`, else horizontal
/// across `x`.
fn threshold_lines(ax: &mut AxesMut<'_>, s: &Value, vertical: bool, span: [f64; 2]) {
    for (name, t) in named_thresholds(s) {
        let (x, y) = if vertical { ([t, t], span) } else { (span, [t, t]) };
        ax.plot(x, y)
            .label(format!("{name} threshold {t:.0e}"))
            .color(GREY)
            .dash(if name == "error" { Dash::Dashed } else { Dash::DashDot })
            .line_width(1.0);
    }
}

pub fn spread_histogram(run_dir: &Path, s: &Value) -> Result<(), String> {
    let (h, nan) = histogram(&worst_spreads(s));
    // each decade at its geometric centre; an empty decade is not drawn on the log count axis
    let centre: Vec<f64> = (0..h.len())
        .map(|i| cr_libm::pow(10.0, (*DECADES.start() + i as i32) as f64 + 0.5))
        .collect();
    let counts = |h: &[usize]| -> Vec<f64> { h.iter().map(|&n| if n == 0 { f64::NAN } else { n as f64 }).collect() };
    let top = h.iter().chain([&nan]).copied().max().unwrap_or(1) as f64 * 2.0;
    let mut fig = figure("spread-histogram", s);
    let k = key(&mut fig, 2);
    {
        let mut ax = fig.axes(0, 0);
        for (i, g) in groups(s).iter().enumerate() {
            let w = worst_spreads_of(s, Some(g));
            let (hg, nan_g) = histogram(&w);
            let (marker, colour) = group_style(g);
            ax.plot(&centre, counts(&hg))
                .label(format!("{g} points measured ({})", w.len()))
                .color(colour)
                .dash(if i == 0 { Dash::Solid } else { Dash::Dotted })
                .line_width(1.0)
                .marker(marker)
                .marker_size(4.0)
                .marker_edge(colour)
                .marker_face(colour);
            if nan_g > 0 {
                ax.plot([NAN_AT], [nan_g as f64])
                    .label(format!("{g}: spread not a number, drawn at 1e2 ({nan_g})"))
                    .color(ORANGE)
                    .dash(Dash::None)
                    .marker(marker)
                    .marker_size(6.0)
                    .marker_edge(ORANGE)
                    .marker_face(ORANGE);
            }
        }
        threshold_lines(&mut ax, s, true, [0.5, top]);
        ax.xscale(Scale::Log)
            .yscale(Scale::Log)
            .xlim(cr_libm::pow(10.0, *DECADES.start() as f64), 10.0 * NAN_AT)
            .ylim(0.5, top)
            .xlabel("largest spread of a reported quantity over five 1-ULP twins")
            .ylabel("points per decade")
            .legend_in(k);
    }
    figures::save_with_pdf(&fig, &run_dir.join("spread-histogram")).map_err(|e| e.to_string())
}

pub fn threshold_plateau(run_dir: &Path, s: &Value) -> Result<(), String> {
    let w = worst_spreads(s);
    let ts = plateau_thresholds();
    let f: Vec<f64> = ts
        .iter()
        .map(|&t| w.iter().filter(|&&p| flagged(p, t)).count() as f64)
        .collect();
    let top = f.iter().copied().fold(1.0, f64::max) * 2.0;
    let bottom = 1.0;
    let mut fig = figure("threshold-plateau", s);
    let k = key(&mut fig, 1);
    {
        let mut ax = fig.axes(0, 0);
        ax.plot(&ts, &f)
            .label(format!("points beyond it, of {} measured in all groups", w.len()))
            .color(ALL)
            .line_width(1.0)
            .marker(Marker::Circle)
            .marker_size(4.0)
            .marker_edge(ALL)
            .marker_face(ALL);
        threshold_lines(&mut ax, s, true, [bottom, top]);
        // the gap around each threshold: the same points are beyond any threshold in it. Every
        // edge as one series, broken between them, so the key has one entry
        let mut edges: Vec<f64> = named_thresholds(s)
            .iter()
            .flat_map(|(_, t)| {
                let (lo, hi) = gap(w.iter().copied(), *t).interval();
                [lo, hi]
            })
            .collect();
        edges.dedup();
        let (mut ex, mut ey) = (vec![], vec![]);
        for edge in edges.into_iter().filter(|e| *e > 0.0 && e.is_finite()) {
            if !ex.is_empty() {
                ex.push(f64::NAN);
                ey.push(f64::NAN);
            }
            ex.extend([edge, edge]);
            ey.extend([bottom, top]);
        }
        if !ex.is_empty() {
            ax.plot(&ex, &ey)
                .label("the gap around each threshold (largest spread below, smallest above)")
                .color(ALL)
                .dash(Dash::Dotted)
                .line_width(1.0);
        }
        ax.xscale(Scale::Log)
            .yscale(Scale::Log)
            .xlim(ts[0], *ts.last().unwrap())
            .ylim(bottom, top)
            .xlabel("threshold on the largest spread")
            .ylabel("points beyond the threshold")
            .legend_in(k);
    }
    figures::save_with_pdf(&fig, &run_dir.join("threshold-plateau")).map_err(|e| e.to_string())
}

/// Every point's largest spread against its α, all cases together: points within the warning
/// threshold open, points between the thresholds half-filled, ill-conditioned points filled, a
/// spread that is not a number at the top, and the points invalid on a single-run reason (whose
/// twins are not compared) as faint marks along the bottom.
pub fn spread_against_alpha(run_dir: &Path, s: &Value) -> Result<(), String> {
    let t = thresholds(s);
    let (floor, top) = (1e-17, 1e2);
    let mut invalid = vec![];
    for p in points(s) {
        let a = val(&p["alpha_deg"]);
        if !measured(p) && a.is_finite() {
            invalid.push(a);
        }
    }
    let mut fig = figure("spread-against-alpha", s).title(
        "Spread of $q_v$ over twins $q_{v,k}$, \
         $v \\in \\{\\alpha, C_L, C_M, C_D, C_{Dp}, C_{Df}, x_\\mathrm{tr}\\}$, $k \\in \\{1, \\dots, 5\\}$",
    );
    let k = key(&mut fig, 1);
    {
        let mut ax = fig.axes(0, 0);
        if !invalid.is_empty() {
            ax.plot(&invalid, vec![2.0 * floor; invalid.len()])
                .label(format!("single-run invalid, not measured ({})", invalid.len()))
                .color(GREY)
                .dash(Dash::None)
                .marker(Marker::Diamond)
                .marker_size(3.0)
                .marker_edge(GREY)
                .marker_face(figures::with_alpha(GREY, 0.5));
        }
        for g in groups(s) {
            let (mut ok, mut warned, mut ill) = ((vec![], vec![]), (vec![], vec![]), (vec![], vec![]));
            let mut nan = vec![];
            for p in points_of(s, Some(&g)).filter(|p| measured(p)) {
                let (a, w) = (val(&p["alpha_deg"]), val(&p["worst_spread"]));
                if w.is_nan() {
                    nan.push(a);
                } else if flagged(w, t.error) {
                    ill.0.push(a);
                    ill.1.push(w);
                } else if flagged(w, t.warning) {
                    warned.0.push(a);
                    warned.1.push(w);
                } else {
                    ok.0.push(a);
                    ok.1.push(w.max(floor));
                }
            }
            let (marker, colour) = group_style(&g);
            ax.plot(&ok.0, &ok.1)
                .label(format!("{g}: within the warning threshold ({})", ok.0.len()))
                .color(colour)
                .dash(Dash::None)
                .marker(marker)
                .marker_size(3.0)
                .marker_edge(colour);
            if !warned.0.is_empty() {
                ax.plot(&warned.0, &warned.1)
                    .label(format!("{g}: may be ill-conditioned ({})", warned.0.len()))
                    .color(colour)
                    .dash(Dash::None)
                    .marker(marker)
                    .marker_size(4.5)
                    .marker_edge(colour)
                    .marker_face(figures::with_alpha(colour, 0.4));
            }
            if !ill.0.is_empty() {
                ax.plot(&ill.0, &ill.1)
                    .label(format!("{g}: ill-conditioned ({})", ill.0.len()))
                    .color(colour)
                    .dash(Dash::None)
                    .marker(marker)
                    .marker_size(6.0)
                    .marker_edge(colour)
                    .marker_face(colour);
            }
            if !nan.is_empty() {
                ax.plot(&nan, vec![top / 2.0; nan.len()])
                    .label(format!("{g}: NaN spread, at the top ({})", nan.len()))
                    .color(ORANGE)
                    .dash(Dash::None)
                    .marker(marker)
                    .marker_size(6.0)
                    .marker_edge(ORANGE)
                    .marker_face(ORANGE);
            }
        }
        threshold_lines(&mut ax, s, false, [-31.0, 31.0]);
        ax.yscale(Scale::Log)
            .xlim(-31.0, 31.0)
            .ylim(floor, top)
            .xlabel("$\\alpha$ (°)")
            .ylabel("$\\max_v \\max_k |q_v - q_{v,k}| / \\max(|q_v|, |q_{v,k}|, 1)$")
            .legend_in(k);
    }
    figures::save_with_pdf(&fig, &run_dir.join("spread-against-alpha")).map_err(|e| e.to_string())
}

pub fn render(run_dir: &Path, s: &Value) {
    for (name, r) in [
        ("spread-histogram", spread_histogram(run_dir, s)),
        ("threshold-plateau", threshold_plateau(run_dir, s)),
        ("spread-against-alpha", spread_against_alpha(run_dir, s)),
    ] {
        if let Err(e) = r {
            eprintln!("{name} not drawn: {e}");
        }
    }
}
