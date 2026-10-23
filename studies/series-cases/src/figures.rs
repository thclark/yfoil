//! The study's figures, drawn from a run folder's `summary.json` as ironlab figure files
//! (`docs/conventions/figures.md`): one polar sheet per case (`polars-<case>.fig`) and the
//! transition-position strip (`transition-position.fig`), each with its PDF for `README.tex`.
//!
//! Every number on a figure, axis limits included, is computed here from `summary.json`.

use crate::{Case, Study};
use cr_libm::{log10, pow};
use ironlab::prelude::*;
use serde_json::Value;
use std::path::Path;
use yfoil::output::figures::{self, FigureError, BLUE, ORANGE, RED, WIDTH_MM};

const PT_TO_MM: f64 = 25.4 / 72.0;
/// XFOIL's α range for every polar of the study.
const ALPHA_RANGE: (f64, f64) = (-30.0, 30.0);

/// Legend names: one entry per look, shared by every panel of a sheet.
const XFOIL_NAME: &str = "XFOIL (reference)";
const YFOIL_NAME: &str = "yFoil";
const FOLLOWED: &str = "yFoil point: every iteration agrees with XFOIL's";
const INFORMATIONAL: &str = "yFoil point: parts from XFOIL within it, or after a point that did (informational)";
const NONFINITE: &str = "went through a non-finite value within the point";
const NAN_TICK: &str = "plotted value is NaN (marker at the top of the panel at its α)";
const ILL: &str = "parts where XFOIL is ill-conditioned (its own 1-ULP twins part there too)";
const MISMATCH: &str = "parts where XFOIL is well-conditioned (its twins agree): a mismatch, replayed step by step";
const STATIONS: &str = "final stagnation or transition station (IST, ITRAN) differs between the codes";

const TRANSITION_FOLLOWED: &str = "every iteration of the point agrees with XFOIL's";
const TRANSITION_INFORMATIONAL: &str = "parts from XFOIL within the point, or after a point that did (informational)";

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

/// Draw marks of one look as a marker-only series; nothing when there are none, so the key lists
/// only what the figure shows. A non-positive value on a logarithmic axis is not drawn.
fn draw(ax: &mut AxesMut<'_>, m: &Marks, look: &Look, log: bool) {
    if m.x.is_empty() {
        return;
    }
    let y: Vec<f64> =
        m.y.iter()
            .map(|&v| if log && v <= 0.0 { f64::NAN } else { v })
            .collect();
    let mut line = ax.plot(&m.x, &y);
    line.display_name(look.name)
        .color(look.edge)
        .dash(Dash::None)
        .line_width(look.edge_width)
        .marker(look.shape)
        .marker_size(look.size)
        .marker_edge(look.edge);
    if let Some(f) = look.face {
        line.marker_face(f);
    }
}

fn grey_diamond() -> Color {
    figures::with_alpha(Color::BLACK, 0.35)
}

/// The figure parameters and labels of one case (`docs/conventions/figures.md`).
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
    fig = fig.label(study.section()).label("xfoil-comparison");
    for tag in regime_labels(case) {
        fig = fig.label(tag);
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

/// Where a marker sits at the top of a panel: 1.5 % of the panel's height below its top edge.
fn near_top(lo: f64, hi: f64, log: bool) -> f64 {
    if log {
        let (l, h) = (log10(lo), log10(hi));
        pow(10.0, h - 0.015 * (h - l))
    } else {
        hi - 0.015 * (hi - lo)
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
    // the matplotlib sheet's proportions: 118 pt panels 30 pt apart, a 107 pt legend strip
    let (panel_pt, gap_pt, legend_pt, frame_pt) = (118.0, 30.0, 107.0, 42.0);
    let height_pt = frame_pt + rows as f64 * panel_pt + (rows - 1) as f64 * gap_pt + legend_pt;
    let row_weights: Vec<f64> = (0..rows)
        .map(|_| 1.0)
        .chain([legend_pt / (panel_pt + gap_pt)])
        .collect();
    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height_pt * PT_TO_MM)
        .tiles(rows + 1, 2)
        .row_weights(row_weights)
        .title(Text::plain(format!(
            "{}: {}",
            run["section"].as_str().unwrap_or(""),
            case.name
        )));
    fig = describe(fig, study, "polar-sheet", case, run);
    let key = {
        let mut legend = fig.panel_legend_span(rows, 0, 1, 2);
        legend
            .boxed(false)
            .columns(LegendEntryColumns::Count(1))
            .align(LegendLocation::West);
        legend.id()
    };

    let mut ids = vec![];
    for (i, pn) in panels.iter().enumerate() {
        let (row, col) = ((i / 2) as u32, (i % 2) as u32);
        let mut all: Vec<f64> = vec![];
        let mut nan_ticks: Vec<f64> = vec![];
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
                        None => nan_ticks.push(alpha(p)),
                        Some(y) if p["xfoil_nonfinite"].as_bool().unwrap_or(false) => xfoil_nonfinite.push(alpha(p), y),
                        _ => {}
                    }
                    (alpha(p), usable(v).unwrap_or(f64::NAN))
                })
                .unzip();
            if !xa.is_empty() {
                all.extend(&xy);
                ax.plot(&xa, &xy)
                    .display_name(XFOIL_NAME)
                    .color(RED)
                    .dash(Dash::Dashed)
                    .line_width(0.9)
                    .marker(Marker::Cross)
                    .marker_size(3.2);
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
                        nan_ticks.push(a);
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
                ax.plot(&ya, &yy).display_name(YFOIL_NAME).color(BLUE).line_width(0.8);
            }
            let look = |name, shape, size, edge, face| Look {
                name,
                shape,
                size,
                edge,
                face,
                edge_width: 0.6,
            };
            draw(
                &mut ax,
                &xfoil_nonfinite,
                &look(NONFINITE, Marker::Square, 5.2, ORANGE, None),
                pn.log,
            );
            draw(
                &mut ax,
                &followed,
                &look(FOLLOWED, Marker::Circle, 2.4, BLUE, Some(BLUE)),
                pn.log,
            );
            draw(
                &mut ax,
                &informational,
                &look(INFORMATIONAL, Marker::Circle, 2.4, BLUE, None),
                pn.log,
            );
            draw(
                &mut ax,
                &yfoil_nonfinite,
                &look(NONFINITE, Marker::Square, 5.2, ORANGE, None),
                pn.log,
            );
            draw(
                &mut ax,
                &ill,
                &look(ILL, Marker::Diamond, 5.8, Color::BLACK, None),
                pn.log,
            );
            let grey = grey_diamond();
            draw(
                &mut ax,
                &mismatch,
                &look(MISMATCH, Marker::Diamond, 5.8, grey, Some(grey)),
                pn.log,
            );
            draw(
                &mut ax,
                &stations,
                &look(STATIONS, Marker::TriangleDown, 6.5, Color::BLACK, None),
                pn.log,
            );
        }
        // reference lines: the tolerance itself, and the ends and middle of the station interval
        let reference_lines: &[f64] = match pn.key {
            "state_ratio" => &[1.0],
            "trpos_upper" | "trpos_lower" => &[0.0, 0.5, 1.0],
            _ => &[],
        };
        all.extend(reference_lines);
        let (lo, hi) = match pn.key {
            "trpos_upper" | "trpos_lower" => (-0.02, 1.02),
            "rmsbl" => limits(&all, true, Some(1e-8)),
            _ => limits(&all, pn.log, None),
        };
        for &y in reference_lines {
            ax.plot([ALPHA_RANGE.0, ALPHA_RANGE.1], [y, y])
                .color(Color::BLACK)
                .line_width(0.5)
                .dash(Dash::Dotted);
        }
        if !nan_ticks.is_empty() {
            // an axes-fraction tick in the matplotlib figure; here a marker at the panel's top
            let top = near_top(lo, hi, pn.log);
            ax.plot(&nan_ticks, vec![top; nan_ticks.len()])
                .display_name(NAN_TICK)
                .color(ORANGE)
                .dash(Dash::None)
                .marker(Marker::Plus)
                .marker_size(4.0);
        }
        ax.xlabel("$\\alpha$ / deg")
            .ylabel(pn.ylabel)
            .xlim(ALPHA_RANGE.0, ALPHA_RANGE.1)
            .ylim(lo, hi)
            .legend_in(key);
        if pn.log {
            ax.yscale(Scale::Log);
        }
        ids.push((pn.key, ax.id(), lo, hi));
    }
    fig.link_all_x().expect("link α");
    // the upper/lower pairs plot one quantity: one y range for both
    for pair in [["xtr_upper", "xtr_lower"], ["trpos_upper", "trpos_lower"]] {
        let found: Vec<_> = ids.iter().filter(|(k, ..)| pair.contains(k)).collect();
        if let [a, b] = found.as_slice() {
            let (lo, hi) = (a.2.min(b.2), a.3.max(b.3));
            fig.axes((ids.iter().position(|i| i.1 == a.1).unwrap() / 2) as u32, 0)
                .ylim(lo, hi);
            fig.link(Dim::Y, &[a.1, b.1]).expect("link the pair");
        }
    }
    fig
}

/// Every recorded point of every case: where transition sits within its station interval, one
/// row per case (numbered, the key names the case), one panel per side.
pub fn transition_strip(study: &Study, cases: &[&Case], runs: &[Value]) -> Figure {
    let n = runs.len().max(1);
    let panel_pt = (11.0 * n as f64 + 8.0).max(120.0);
    let key_rows = (runs.len() + 4) as f64;
    let legend_pt = 11.0 * key_rows + 10.0;
    let height_pt = 42.0 + panel_pt + legend_pt;
    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height_pt * PT_TO_MM)
        .tiles(2, 2)
        .row_weights([1.0, legend_pt / panel_pt])
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
        .label(study.section())
        .label("xfoil-comparison");
    let mut tags: Vec<&str> = cases.iter().flat_map(|c| regime_labels(c)).collect();
    tags.sort_unstable();
    tags.dedup();
    for t in tags {
        fig = fig.label(t);
    }
    let key = {
        let mut legend = fig.panel_legend_span(1, 0, 1, 2);
        legend
            .boxed(false)
            .columns(LegendEntryColumns::Count(1))
            .align(LegendLocation::NorthWest);
        legend.id()
    };
    let (ylo, yhi) = (0.5, n as f64 + 0.5);
    let mut ids = vec![];
    for (col, label) in [(0u32, "upper side"), (1, "lower side")] {
        let mut ax = fig.axes(0, col);
        let mut followed = Marks::default();
        let mut informational = Marks::default();
        let mut ill = Marks::default();
        let mut mismatch = Marks::default();
        for (row, run) in runs.iter().enumerate() {
            let y = (row + 1) as f64;
            let mut case_marks = Marks::default();
            for p in run["points"].as_array().into_iter().flatten() {
                let t = &p["transition_position"][col as usize];
                if t.is_null() || t["at_te"].as_bool().unwrap_or(false) {
                    continue;
                }
                let Some(x) = t["fraction"].as_f64().filter(|x| x.is_finite()) else {
                    continue;
                };
                case_marks.push(x, y);
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
            // a row label in the matplotlib figure; here a transparent series named by the case,
            // so the key and the data tips say which row is which (a case with no transition on
            // the aerofoil, an inviscid one, keeps its row and its key entry with a NaN point)
            if case_marks.x.is_empty() {
                case_marks.push(f64::NAN, y);
            }
            {
                ax.plot(&case_marks.x, &case_marks.y)
                    .display_name(Text::plain(format!(
                        "{}: {}",
                        row + 1,
                        run["name"].as_str().unwrap_or("")
                    )))
                    .color(Color::rgba(0.0, 0.0, 0.0, 0.0))
                    .dash(Dash::None)
                    .marker(Marker::Circle)
                    .marker_size(5.0);
            }
        }
        for x in [0.0, 0.5, 1.0] {
            ax.plot([x, x], [ylo, yhi])
                .color(Color::BLACK)
                .line_width(0.4)
                .dash(Dash::Dotted);
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
            &look(ILL, Marker::Diamond, 5.0, Color::BLACK, None),
            false,
        );
        draw(
            &mut ax,
            &mismatch,
            &look(MISMATCH, Marker::Diamond, 5.0, grey, Some(grey)),
            false,
        );
        ax.xlabel(format!("position of transition in its station interval, {label}"))
            .ylabel(if col == 0 { "case (number in the key)" } else { "" })
            .xlim(-0.02, 1.02)
            .ylim(ylo, yhi)
            .legend_in(key);
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
    fn nan_marker_sits_inside_the_panel() {
        let y = near_top(1e-8, 1.0, true);
        assert!(y < 1.0 && y > 1e-8);
        let y = near_top(0.0, 1.0, false);
        assert!(y < 1.0 && y > 0.98);
    }
}
