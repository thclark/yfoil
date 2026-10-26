//! The study's figures, drawn from a run folder's `summary.json` as ironlab figure files
//! (`docs/development/conventions/figures.md`): one polar sheet per case (`polars-<case>.fig`) and the
//! transition-position strip (`transition-position.fig`), each with its PDF for `README.tex`.
//!
//! Every number on a figure, axis limits included, is computed here from `summary.json`.

use crate::{Case, Study};
use cr_libm::{log10, pow};
use ironlab::prelude::*;
use ironlab::{LineStyle, MarkerStyle};
use serde_json::Value;
use std::path::Path;
use yfoil::output::figures::{self, FigureError, BLUE, ORANGE, RED, WIDTH_MM};

const PT_TO_MM: f64 = 25.4 / 72.0;
/// XFOIL's α range for every polar of the study.
const ALPHA_RANGE: (f64, f64) = (-30.0, 30.0);
/// Where a NaN tick sits: at its α, this fraction of the panel's height up, whatever the limits.
const NAN_TICK_HEIGHT: f64 = 0.985;

/// Legend names: one entry per look, shared by every panel of a sheet.
const XFOIL_NAME: &str = "XFOIL (reference)";
const YFOIL_NAME: &str = "yFoil";
const FOLLOWED: &str = "yFoil point: every iteration agrees with XFOIL's";
const INFORMATIONAL: &str = "yFoil point: parts from XFOIL within it, or after a point that did (informational)";
const NONFINITE: &str = "went through a non-finite value within the point";
const NAN_TICK: &str = "plotted value is NaN (tick at the top of the panel at its α)";
const ILL: &str = "parts where XFOIL is ill-conditioned (its own 1-ULP twins part there too)";
const MISMATCH: &str = "parts where XFOIL is well-conditioned (its twins agree): a mismatch, replayed step by step";
const STATIONS: &str = "final stagnation or transition station (IST, ITRAN) differs between the codes";

// The strip's key is as wide as its two plot areas, which the case names leave narrow, so its
// entries are the short forms of the sheet's.
const TRANSITION_FOLLOWED: &str = "followed: every iteration agrees with XFOIL's";
const TRANSITION_INFORMATIONAL: &str = "informational: parts from XFOIL";
const TRANSITION_ILL: &str = "divergent: XFOIL ill-conditioned there";
const TRANSITION_MISMATCH: &str = "mismatch: XFOIL well-conditioned there";

/// How one kind of mark is drawn and what the key calls it.
struct Look {
    name: &'static str,
    shape: Marker,
    size: f64,
    edge: Color,
    /// filled with this colour; open without
    face: Option<Color>,
    edge_width: f64,
}

/// Draw marks of one look as a marker-only series, returning its id; nothing when there are none,
/// so the key lists only what the figure shows. A non-positive value on a logarithmic axis is not
/// drawn.
fn draw(ax: &mut AxesMut<'_>, m: &Marks, look: &Look, log: bool) -> Option<NodeId> {
    if m.x.is_empty() {
        return None;
    }
    let y: Vec<f64> =
        m.y.iter()
            .map(|&v| if log && v <= 0.0 { f64::NAN } else { v })
            .collect();
    let mut line = ax.plot(&m.x, &y);
    line.label(look.name)
        .color(look.edge)
        .dash(Dash::None)
        .line_width(look.edge_width)
        .marker(look.shape)
        .marker_size(look.size)
        .marker_edge(look.edge);
    if let Some(f) = look.face {
        line.marker_face(f);
    }
    Some(line.id())
}

/// The key of a sheet: the first artist of each name, in the order first drawn. One item stands
/// for every artist of its name and look, in every panel.
#[derive(Default)]
struct Key {
    items: Vec<(&'static str, NodeId)>,
}

impl Key {
    fn add(&mut self, name: &'static str, id: Option<NodeId>) {
        if let Some(id) = id {
            if !self.items.iter().any(|(n, _)| *n == name) {
                self.items.push((name, id));
            }
        }
    }
}

/// Mark, at the top edge of a panel, the α of each value that is NaN: a vertical tick at data α
/// and an axes-fraction height, so it stays on the edge whatever the limits, pan or zoom, and
/// takes no part in them. It is unnamed: the key explains the ticks once ([`nan_tick_sample`]).
fn nan_ticks(ax: &mut AxesMut<'_>, alpha: &[f64]) {
    if alpha.is_empty() {
        return;
    }
    ax.plot(alpha, vec![NAN_TICK_HEIGHT; alpha.len()])
        .y_coords(Coords::AxesFraction)
        .clip(false)
        .color(ORANGE)
        .dash(Dash::None)
        .marker(Marker::VLine)
        .marker_size(4.0)
        .marker_edge(ORANGE);
}

/// The key's sample for the NaN ticks, which stand for no single artist.
fn nan_tick_sample() -> LegendSample {
    LegendSample::Line {
        line: LineStyle {
            color: ColorSpec::Rgba { color: ORANGE },
            dash: Dash::None,
            ..LineStyle::default()
        },
        marker: MarkerStyle {
            shape: Marker::VLine,
            size_pt: 4.0,
            face: ColorSpec::None,
            edge: ColorSpec::Rgba { color: ORANGE },
        },
    }
}

/// A case name (a code identifier) in monospace, its LaTeX specials escaped.
fn mono(name: &str) -> String {
    let mut out = String::with_capacity(name.len() + 10);
    for c in name.chars() {
        if matches!(c, '_' | '{' | '}' | '$' | '&' | '%' | '#') {
            out.push('\\');
        }
        out.push(c);
    }
    format!("\\texttt{{{out}}}")
}

/// Plain text with the one character that opens LaTeX math outside it escaped.
fn text(s: &str) -> String {
    s.replace('$', "\\$")
}

fn grey_diamond() -> Color {
    figures::with_alpha(Color::BLACK, 0.35)
}

/// The figure parameters and labels of one case (`docs/development/conventions/figures.md`).
fn describe(fig: Figure, study: &Study, figure: &str, case: &Case, run: &Value) -> Figure {
    let mut fig = fig
        .parameter("study", study.slug)
        .parameter("figure", figure)
        .parameter("case", case.name.as_str())
        .parameter("foil", case.foil.as_str())
        .parameter("n_nodes", case.n_nodes as i64);
    fig = figures::flow_condition_parameters(fig, &case.conditions());
    fig = fig
        .parameter("points", run["yfoil_points"].as_i64().unwrap_or(0))
        .parameter("matched", run["matched"].as_i64().unwrap_or(0));
    if let Some([leg, alpha]) = run["parted"].as_array().map(Vec::as_slice) {
        fig = fig.parameter("parted", format!("leg {leg} α {alpha}°"));
    }
    fig = fig.tag("validation").tag("xfoil-comparison");
    for tag in regime_labels(case) {
        fig = fig.tag(tag);
    }
    fig
}

/// The regime tags a case's inputs imply.
fn regime_labels(case: &Case) -> Vec<&'static str> {
    let mut tags = vec![];
    if case.inviscid {
        tags.push("inviscid");
    }
    if case.mach > 0.0 {
        tags.push("compressible");
    }
    if case.xtr.iter().any(|&x| x < 1.0) {
        tags.push("forced-transition");
    }
    if case.alphas.is_empty() && !case.cls.is_empty() {
        tags.push("fixed-cl");
    }
    if !case.inviscid && case.re < 1e5 {
        tags.push("low-re");
    }
    tags
}

/// A panel of the polar sheet: the quantity, its axis label, whether its axis is logarithmic.
struct Panel {
    key: &'static str,
    ylabel: &'static str,
    log: bool,
}

const fn panel(key: &'static str, ylabel: &'static str, log: bool) -> Panel {
    Panel { key, ylabel, log }
}

/// The value a panel plots for one code at one point; `None` where there is nothing to plot
/// (the value is null in `summary.json`, which is how serde writes NaN, or the code has no
/// such value).
fn value(p: &Value, code: &str, key: &str) -> Option<f64> {
    let side = |k: &str| if k.ends_with("upper") { 0 } else { 1 };
    match key {
        // yFoil's state only (the reference agrees in ITRAN wherever the point was followed)
        "trpos_upper" | "trpos_lower" => {
            if code == "xfoil" {
                return None;
            }
            let t = &p["transition_position"][side(key)];
            if t.is_null() || t["at_te"].as_bool().unwrap_or(false) {
                None
            } else {
                t["fraction"].as_f64()
            }
        }
        // the reference has no ratio: its own twin's spread is the floor the tolerance is built on
        "state_ratio" => {
            if code == "yfoil" {
                p["state"][0]["ratio"].as_f64()
            } else {
                None
            }
        }
        "state_diff" => p["state"][0][if code == "yfoil" { "diff" } else { "floor" }].as_f64(),
        "xtr_upper" | "xtr_lower" => p[format!("{code}_transition")][side(key)].as_f64(),
        _ => p[code][key].as_f64(),
    }
    .filter(|v| v.is_finite())
}

/// Collected marks of one look in one panel.
#[derive(Default)]
struct Marks {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl Marks {
    fn push(&mut self, x: f64, y: f64) {
        self.x.push(x);
        self.y.push(y);
    }
}

/// Axis limits from the finite values of a panel: padded by 5 % of the span (in decades on a
/// logarithmic axis), with a floor where the panel fixes one.
fn limits(values: &[f64], log: bool, floor: Option<f64>) -> (f64, f64) {
    let vals: Vec<f64> = values
        .iter()
        .copied()
        .filter(|v| v.is_finite() && (!log || *v > 0.0))
        .map(|v| if log { log10(v) } else { v })
        .collect();
    let (mut lo, mut hi) = if vals.is_empty() {
        if log {
            (-8.0, 0.0)
        } else {
            (0.0, 1.0)
        }
    } else {
        let lo = vals.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let span = hi - lo;
        if span > 0.0 {
            (lo - 0.05 * span, hi + 0.05 * span)
        } else {
            let pad = if lo.abs() > 0.0 && !log { 0.05 * lo.abs() } else { 0.5 };
            (lo - pad, hi + pad)
        }
    };
    if let Some(f) = floor {
        lo = if log { log10(f) } else { f };
        if hi <= lo {
            hi = lo + 1.0;
        }
    }
    if log {
        (pow(10.0, lo), pow(10.0, hi))
    } else {
        (lo, hi)
    }
}

/// One case's polar sheet: CL and CD (top), the last RMSBL and the iterations (second row), the
/// transition location on each side (third), the position of transition within its station
/// interval on each side (fourth) and, where the fixture kept the reference's per-call state,
/// yFoil's final state against the reference's (fifth).
pub fn polar_sheet(study: &Study, case: &Case, run: &Value) -> Figure {
    let points: Vec<&Value> = run["points"].as_array().map(|a| a.iter().collect()).unwrap_or_default();
    let has_state = points.iter().any(|p| p["state"].is_array());
    let mut panels = vec![
        panel("cl", "$C_L$", false),
        panel("cd", "$C_D$", false),
        panel("rmsbl", "RMSBL at the last iteration", true),
        panel("iterations", "iterations", false),
        panel("xtr_upper", "$x_{tr}/c$, upper side", false),
        panel("xtr_lower", "$x_{tr}/c$, lower side", false),
        panel("trpos_upper", "transition position in its interval, upper", false),
        panel("trpos_lower", "transition position in its interval, lower", false),
    ];
    if has_state {
        panels.push(panel("state_ratio", "worst state array: $|\\Delta|$ / tolerance", true));
        panels.push(panel("state_diff", "worst state array: $|\\Delta|$, floor (red)", true));
    }
    let rows = panels.len().div_ceil(2) as u32;
    // the matplotlib sheet's proportions: 118 pt panels 30 pt apart, about 107 pt for the key,
    // whose row is fitted to it
    let (panel_pt, gap_pt, legend_pt, frame_pt) = (118.0, 30.0, 107.0, 42.0);
    let height_pt = frame_pt + rows as f64 * panel_pt + (rows - 1) as f64 * gap_pt + legend_pt;
    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height_pt * PT_TO_MM)
        .tiles(rows + 1, 2)
        .fit_row(rows)
        .title(format!(
            "{}: {}",
            text(run["section"].as_str().unwrap_or("")),
            mono(&case.name)
        ));
    fig = describe(fig, study, "polar-sheet", case, run);
    let key_panel = fig.panel_span(rows, 0, 1, 2).align(PanelLocation::West).legend().id();
    let mut key = Key::default();
    let mut any_nan = false;

    let mut ids = vec![];
    for (i, pn) in panels.iter().enumerate() {
        let (row, col) = ((i / 2) as u32, (i % 2) as u32);
        let mut all: Vec<f64> = vec![];
        let mut nan_ticks_at: Vec<f64> = vec![];
        let mut ax = fig.axes(row, col);
        let usable = |v: Option<f64>| v.filter(|v| !pn.log || *v > 0.0);
        for leg in [1, 2] {
            let mut pts: Vec<&Value> = points.iter().copied().filter(|p| p["leg"] == leg).collect();
            pts.sort_by(|a, b| a["alpha_deg"].as_f64().partial_cmp(&b["alpha_deg"].as_f64()).unwrap());
            let alpha = |p: &Value| p["alpha_deg"].as_f64().unwrap_or(f64::NAN);
            // the reference first, so that yFoil's markers sit on top where they coincide; it has
            // no series in the ratio and transition-position panels
            let xfoil_pts: Vec<&Value> = if matches!(pn.key, "state_ratio" | "trpos_upper" | "trpos_lower") {
                vec![]
            } else {
                pts.iter().copied().filter(|p| p["xfoil"].is_object()).collect()
            };
            let mut xfoil_nonfinite = Marks::default();
            let (xa, xy): (Vec<f64>, Vec<f64>) = xfoil_pts
                .iter()
                .map(|p| {
                    let v = value(p, "xfoil", pn.key);
                    match v {
                        None => nan_ticks_at.push(alpha(p)),
                        Some(y) if p["xfoil_nonfinite"].as_bool().unwrap_or(false) => xfoil_nonfinite.push(alpha(p), y),
                        _ => {}
                    }
                    (alpha(p), usable(v).unwrap_or(f64::NAN))
                })
                .unzip();
            if !xa.is_empty() {
                all.extend(&xy);
                let id = ax
                    .plot(&xa, &xy)
                    .label(XFOIL_NAME)
                    .color(RED)
                    .dash(Dash::Dashed)
                    .line_width(0.9)
                    .marker(Marker::Cross)
                    .marker_size(3.2)
                    .id();
                key.add(XFOIL_NAME, Some(id));
            }
            let mut followed = Marks::default();
            let mut informational = Marks::default();
            let mut yfoil_nonfinite = Marks::default();
            let mut ill = Marks::default();
            let mut mismatch = Marks::default();
            let mut stations = Marks::default();
            let (ya, yy): (Vec<f64>, Vec<f64>) = pts
                .iter()
                .map(|p| {
                    let a = alpha(p);
                    let v = value(p, "yfoil", pn.key);
                    let Some(y) = v else {
                        nan_ticks_at.push(a);
                        return (a, f64::NAN);
                    };
                    let filled = p["followed"].as_bool().unwrap_or(false);
                    if filled {
                        followed.push(a, y);
                    } else {
                        informational.push(a, y);
                    }
                    if p["yfoil_nonfinite"].as_bool().unwrap_or(false) {
                        yfoil_nonfinite.push(a, y);
                    }
                    match p["outcome"].as_str() {
                        Some("divergent") if !filled => ill.push(a, y),
                        Some("mismatch") => mismatch.push(a, y),
                        _ => {}
                    }
                    if pn.key.starts_with("xtr_") && p["stations_differ"].as_bool().unwrap_or(false) {
                        stations.push(a, y);
                    }
                    (a, usable(v).unwrap_or(f64::NAN))
                })
                .unzip();
            all.extend(&yy);
            if !ya.is_empty() {
                let id = ax.plot(&ya, &yy).label(YFOIL_NAME).color(BLUE).line_width(0.8).id();
                key.add(YFOIL_NAME, Some(id));
            }
            let look = |name, shape, size, edge, face| Look {
                name,
                shape,
                size,
                edge,
                face,
                edge_width: 0.6,
            };
            key.add(
                NONFINITE,
                draw(
                    &mut ax,
                    &xfoil_nonfinite,
                    &look(NONFINITE, Marker::Square, 5.2, ORANGE, None),
                    pn.log,
                ),
            );
            key.add(
                FOLLOWED,
                draw(
                    &mut ax,
                    &followed,
                    &look(FOLLOWED, Marker::Circle, 2.4, BLUE, Some(BLUE)),
                    pn.log,
                ),
            );
            key.add(
                INFORMATIONAL,
                draw(
                    &mut ax,
                    &informational,
                    &look(INFORMATIONAL, Marker::Circle, 2.4, BLUE, None),
                    pn.log,
                ),
            );
            key.add(
                NONFINITE,
                draw(
                    &mut ax,
                    &yfoil_nonfinite,
                    &look(NONFINITE, Marker::Square, 5.2, ORANGE, None),
                    pn.log,
                ),
            );
            key.add(
                ILL,
                draw(
                    &mut ax,
                    &ill,
                    &look(ILL, Marker::Diamond, 5.8, Color::BLACK, None),
                    pn.log,
                ),
            );
            let grey = grey_diamond();
            key.add(
                MISMATCH,
                draw(
                    &mut ax,
                    &mismatch,
                    &look(MISMATCH, Marker::Diamond, 5.8, grey, Some(grey)),
                    pn.log,
                ),
            );
            key.add(
                STATIONS,
                draw(
                    &mut ax,
                    &stations,
                    &look(STATIONS, Marker::TriangleDown, 6.5, Color::BLACK, None),
                    pn.log,
                ),
            );
        }
        // reference lines across the panel at any zoom: the tolerance itself, and the ends and
        // middle of the station interval
        let reference_lines: &[f64] = match pn.key {
            "state_ratio" => &[1.0],
            "trpos_upper" | "trpos_lower" => &[0.0, 0.5, 1.0],
            _ => &[],
        };
        for &y in reference_lines {
            ax.hline(y).color(Color::BLACK).line_width(0.5).dash(Dash::Dotted);
        }
        any_nan |= !nan_ticks_at.is_empty();
        nan_ticks(&mut ax, &nan_ticks_at);
        ax.xlabel("$\\alpha$ / deg")
            .ylabel(pn.ylabel)
            .xlim(ALPHA_RANGE.0, ALPHA_RANGE.1)
            .legend_in(key_panel);
        // fixed limits only where they help reading: the interval with its ends clear of the box,
        // and RMSBL from the convergence threshold's decades up; every other panel autoscales
        match pn.key {
            "trpos_upper" | "trpos_lower" => {
                ax.ylim(-0.02, 1.02);
            }
            "rmsbl" => {
                let (lo, hi) = limits(&all, true, Some(1e-8));
                ax.ylim(lo, hi);
            }
            _ => {}
        }
        if pn.log {
            ax.yscale(Scale::Log);
        }
        ids.push((pn.key, ax.id()));
    }
    fig.link_all_x().expect("link α");
    // the upper/lower pairs plot one quantity: one y range for both
    for pair in [["xtr_upper", "xtr_lower"], ["trpos_upper", "trpos_lower"]] {
        let found: Vec<NodeId> = ids
            .iter()
            .filter(|(k, _)| pair.contains(k))
            .map(|(_, id)| *id)
            .collect();
        if found.len() == 2 {
            fig.link(Dim::Y, &found).expect("link the pair");
        }
    }
    let mut panel = fig.panel_span(rows, 0, 1, 2);
    let mut legend = panel.legend();
    legend.boxed(false).columns(LegendEntryColumns::Count(1));
    for (_, id) in &key.items {
        legend.item(*id);
    }
    if any_nan {
        legend.custom_item(nan_tick_sample(), NAN_TICK);
    }
    fig
}

/// Every recorded point of every case: where transition sits within its station interval, one
/// row per case, named by the case and first at the top, one panel per side.
pub fn transition_strip(study: &Study, cases: &[&Case], runs: &[Value]) -> Figure {
    let n = runs.len().max(1);
    let panel_pt = (11.0 * n as f64 + 8.0).max(120.0);
    // about four key rows, in a row fitted to the key
    let legend_pt = 11.0 * 4.0 + 10.0;
    let height_pt = 42.0 + panel_pt + legend_pt;
    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height_pt * PT_TO_MM)
        .tiles(2, 2)
        .fit_row(1)
        .title("Position of transition in its station interval, every recorded point of every case")
        .parameter("study", study.slug)
        .parameter("figure", "transition-strip")
        .parameter("cases", runs.len() as i64)
        .parameter(
            "points",
            runs.iter()
                .map(|r| r["yfoil_points"].as_i64().unwrap_or(0))
                .sum::<i64>(),
        )
        .parameter(
            "matched",
            runs.iter().map(|r| r["matched"].as_i64().unwrap_or(0)).sum::<i64>(),
        )
        .tag("validation")
        .tag("xfoil-comparison");
    let mut tags: Vec<&str> = cases.iter().flat_map(|c| regime_labels(c)).collect();
    tags.sort_unstable();
    tags.dedup();
    for t in tags {
        fig = fig.tag(t);
    }
    let key = {
        let mut panel = fig.panel_span(1, 0, 1, 2);
        panel.align(PanelLocation::NorthWest);
        let mut legend = panel.legend();
        legend.boxed(false).columns(LegendEntryColumns::Count(1));
        legend.id()
    };
    // a case with no transition on the aerofoil (an inviscid one) keeps its row: rows are categories
    let names: Vec<String> = runs.iter().map(|r| mono(r["name"].as_str().unwrap_or(""))).collect();
    let mut ids = vec![];
    for (col, label) in [(0u32, "upper side"), (1, "lower side")] {
        let mut ax = fig.axes(0, col);
        let mut followed = Marks::default();
        let mut informational = Marks::default();
        let mut ill = Marks::default();
        let mut mismatch = Marks::default();
        for (row, run) in runs.iter().enumerate() {
            let y = row as f64;
            for p in run["points"].as_array().into_iter().flatten() {
                let t = &p["transition_position"][col as usize];
                if t.is_null() || t["at_te"].as_bool().unwrap_or(false) {
                    continue;
                }
                let Some(x) = t["fraction"].as_f64().filter(|x| x.is_finite()) else {
                    continue;
                };
                if p["followed"].as_bool().unwrap_or(false) {
                    followed.push(x, y);
                } else if p["outcome"] == "mismatch" {
                    mismatch.push(x, y);
                } else if p["outcome"] == "divergent" && p["parted_iteration"].as_u64().unwrap_or(0) > 0 {
                    ill.push(x, y);
                } else {
                    informational.push(x, y);
                }
            }
        }
        // the ends and middle of the interval, across every row at any zoom
        for x in [0.0, 0.5, 1.0] {
            ax.vline(x).color(Color::BLACK).line_width(0.4).dash(Dash::Dotted);
        }
        let blue = figures::with_alpha(BLUE, 0.6);
        let look = |name, shape, size, edge, face| Look {
            name,
            shape,
            size,
            edge,
            face,
            edge_width: 0.5,
        };
        let grey = grey_diamond();
        draw(
            &mut ax,
            &followed,
            &look(TRANSITION_FOLLOWED, Marker::Circle, 2.0, blue, Some(blue)),
            false,
        );
        draw(
            &mut ax,
            &informational,
            &look(TRANSITION_INFORMATIONAL, Marker::Circle, 2.0, blue, None),
            false,
        );
        draw(
            &mut ax,
            &ill,
            &look(TRANSITION_ILL, Marker::Diamond, 5.0, Color::BLACK, None),
            false,
        );
        draw(
            &mut ax,
            &mismatch,
            &look(TRANSITION_MISMATCH, Marker::Diamond, 5.0, grey, Some(grey)),
            false,
        );
        ax.ycategories(names.iter().map(String::as_str))
            .title(label)
            .xlabel("position in the interval")
            .xlim(-0.02, 1.02)
            .ylim(-0.5, n as f64 - 0.5)
            .legend_in(key);
        if col == 0 {
            ax.ylabel("case");
        } else {
            ax.ytick_labels_hidden();
        }
        ids.push(ax.id());
    }
    fig.link_all_x().expect("link the interval position");
    fig.link(Dim::Y, &ids).expect("link the case rows");
    fig
}

/// Draw every figure of a run folder from its `summary.json`, with each case's inputs from
/// `cases.toml`: `polars-<case>.fig` / `.pdf` and `transition-position.fig` / `.pdf`. Returns the
/// stems written.
pub fn render(study: &Study, cases: &[Case], run_dir: &Path, s: &Value) -> Result<Vec<String>, FigureError> {
    let runs: Vec<Value> = s["polars"].as_array().cloned().unwrap_or_default();
    // a run folder holds `.fig` and `.pdf`: the SVGs the matplotlib figures left are stale
    for old in std::fs::read_dir(run_dir).into_iter().flatten().flatten() {
        if old.file_name().to_string_lossy().ends_with(".svg") {
            let _ = std::fs::remove_file(old.path());
        }
    }
    let mut written = vec![];
    let mut used = vec![];
    for run in &runs {
        let name = run["name"].as_str().unwrap_or("");
        let Some(case) = cases.iter().find(|c| c.name == name) else {
            eprintln!("  {name}: not in cases.toml — no figure");
            continue;
        };
        used.push(case);
        let stem = format!("polars-{name}");
        figures::save_with_pdf(&polar_sheet(study, case, run), &run_dir.join(&stem))?;
        written.push(stem);
    }
    let fig = transition_strip(study, &used, &runs);
    figures::save_with_pdf(&fig, &run_dir.join("transition-position"))?;
    written.push("transition-position".into());
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn limits_pad_and_floor() {
        let (lo, hi) = limits(&[0.0, 1.0, f64::NAN], false, None);
        assert!((lo + 0.05).abs() < 1e-12 && (hi - 1.05).abs() < 1e-12);
        let (lo, hi) = limits(&[1e-12, 1e-4, 0.0], true, Some(1e-8));
        assert!((lo - 1e-8).abs() < 1e-20 && hi > 1e-4);
        let (lo, hi) = limits(&[], true, Some(1e-8));
        assert!(lo < hi);
    }

    #[test]
    fn case_names_are_monospace_with_latex_specials_escaped() {
        assert_eq!(mono("naca0012_n160"), "\\texttt{naca0012\\_n160}");
        assert_eq!(text("a $5"), "a \\$5");
    }
}
