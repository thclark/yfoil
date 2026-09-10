//! The seven-row figures: one column per input family, base case in front, perturbation levels
//! stacked behind it coloured by perturbation size (smallest lightest, drawn nearest the base).
//! The node and alpha-step families share a two-column sheet; the panel family has a one-column
//! figure of the same column width. Publication conventions from `figure_style`: physical size,
//! Times New Roman, math-style labels, no in-figure title (the caption carries it).

use crate::run::{Ending, Point, RunResult};
use crate::{perturb, Family, Level};
use figure_style::{self as figure, pxi, pxu, Label, AXIS_LABEL_PT, LEGEND_PT, TEXT_WIDTH_PT, TICK_PT};
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;
use std::path::Path;

/// The families drawn side by side on the sheet; the others get a figure each
pub const SHEET: [Family; 2] = [Family::Geometry, Family::AlphaStep];

/// plotly's YlGnBu colourscale stops (ColorBrewer), light to dark
const YLGNBU: [(u8, u8, u8); 9] = [
    (255, 255, 217),
    (237, 248, 177),
    (199, 233, 180),
    (127, 205, 187),
    (65, 182, 196),
    (29, 145, 192),
    (34, 94, 168),
    (37, 52, 148),
    (8, 29, 88),
];

/// Colour at position t ∈ [0, 1] of the scale, avoiding the near-white end
fn ylgnbu(t: f64) -> RGBColor {
    let t = (0.12 + 0.88 * t.clamp(0.0, 1.0)) * (YLGNBU.len() - 1) as f64;
    let i = (t.floor() as usize).min(YLGNBU.len() - 2);
    let f = t - i as f64;
    let (a, b) = (YLGNBU[i], YLGNBU[i + 1]);
    let mix = |p: u8, q: u8| (p as f64 + (q as f64 - p as f64) * f).round() as u8;
    RGBColor(mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

struct Row {
    label: fn() -> Label,
    log: bool,
    get: fn(&Point) -> f64,
}

/// Rows, top to bottom. "Upper surface" and "trailing edge" are in the caption.
const ROWS: [Row; 7] = [
    Row {
        label: || Label::new().i("C").subi("L"),
        log: false,
        get: |p| p.cl,
    },
    Row {
        label: || Label::new().i("C").subi("D"),
        log: false,
        get: |p| p.cd,
    },
    Row {
        label: || Label::new().i("δ").sup("*").t(" at TE"),
        log: false,
        get: |p| p.dstar_te_upper,
    },
    Row {
        label: || Label::new().i("H").t(" at TE"),
        log: false,
        get: |p| p.h_te_upper,
    },
    Row {
        label: || Label::new().i("x").t("/").i("c").t(" transition"),
        log: false,
        get: |p| p.xtr_upper,
    },
    Row {
        label: || Label::new().i("x").t("/").i("c").t(" separation"),
        log: false,
        get: |p| p.x_sep_upper,
    },
    Row {
        label: || Label::new().t("log").subt("10").t(" RMSBL"),
        log: true,
        get: |p| p.rmsbl,
    },
];

/// Layout, pt. The sheet is `TEXT_WIDTH_PT` wide: two columns and a gutter; a one-column figure
/// is one column wide. Seven rows plus the legend fit an A4 text height (698 pt).
const GUTTER_PT: f64 = 16.0;
const ROW_H_PT: f64 = 82.0;
const X_LABEL_AREA_PT: f64 = 28.0;
const X_STUB_PT: f64 = 4.0;
const Y_LABEL_AREA_PT: f64 = 34.0;
const MARGIN_PT: f64 = 3.0;
const LEGEND_TITLE_PT: f64 = 12.0;
const LEGEND_ROW_PT: f64 = 10.0;
const LEGEND_COLS: usize = 2;
const MARKER_PT: f64 = 1.2;
const OPEN_MARKER_PT: f64 = 2.2;

pub fn column_width_pt() -> f64 {
    (TEXT_WIDTH_PT - GUTTER_PT * (SHEET.len() - 1) as f64) / SHEET.len() as f64
}

/// Legend rows: the base entry on a row of its own, then the perturbation levels in
/// `LEGEND_COLS` columns. The same for every figure (the largest family decides), so that a
/// one-column figure and the sheet have the same height
fn legend_rows() -> usize {
    let n = Family::ALL.iter().map(|f| crate::levels(*f).len()).max().unwrap_or(1);
    1 + (n - 1).div_ceil(LEGEND_COLS)
}

fn legend_height_pt() -> f64 {
    LEGEND_TITLE_PT + LEGEND_ROW_PT * legend_rows() as f64 + 4.0
}

type Area<DB> = DrawingArea<DB, plotters::coord::Shift>;
type Chart<'a, DB> = ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>;

/// The level's legend text, typeset
fn level_label(l: &Level) -> Label {
    let exp = |v: f64| format!("−{}", -(v.log10().round() as i32));
    if l.magnitude == 0.0 {
        return Label::new()
            .t("base: ")
            .i("N")
            .t(&format!(" = {}, Δ", crate::BASE_PANELS))
            .i("α")
            .t(&format!(" = {}°", crate::ALPHA_STEP_DEG));
    }
    match l.family {
        Family::Geometry => {
            if l.raw == perturb::ULP {
                Label::new().t("nodes ± 1 ULP")
            } else {
                Label::new().t("nodes ± 10").sup(&exp(l.raw))
            }
        }
        Family::Panels => {
            let n = l.actual_panels.unwrap_or(l.panels);
            let lab = Label::new().i("N").t(&format!(" = {n}"));
            if n != l.panels {
                lab.t(&format!(" (req. {})", l.panels))
            } else {
                lab
            }
        }
        // the base row states Δα; the levels give only the increment
        Family::AlphaStep => {
            if l.raw == perturb::ULP {
                Label::new().t("Δ").i("α").t(" + 1 ULP")
            } else {
                Label::new().t("Δ").i("α").t(" + 10").sup(&exp(l.raw)).t("°")
            }
        }
    }
}

/// Draw one column per family (a single family gives its figure, the two `SHEET` families the
/// sheet). Every column has the same width whatever the figure.
pub fn draw_families(path: &Path, families: &[(Family, Vec<RunResult>)]) -> Result<(), Box<dyn std::error::Error>> {
    let ncol = families.len();
    let col_w = column_width_pt();
    let width = pxu(col_w * ncol as f64 + GUTTER_PT * (ncol - 1) as f64);
    let legend_h = legend_height_pt();
    let height = pxu(legend_h + ROW_H_PT * ROWS.len() as f64 + X_LABEL_AREA_PT - X_STUB_PT);
    let root = SVGBackend::new(path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    // columns with gutters between them: breakpoints at every column edge, gutters left empty
    let mut xs = Vec::new();
    for c in 1..ncol {
        let x = (col_w + GUTTER_PT) * c as f64;
        xs.push(pxi(x - GUTTER_PT));
        xs.push(pxi(x));
    }
    // rows: legend, then ROWS.len() rows, the last one taller by the x label area
    let mut ys = vec![pxi(legend_h)];
    for r in 1..ROWS.len() {
        ys.push(pxi(legend_h + ROW_H_PT * r as f64));
    }
    let grid = root.split_by_breakpoints(xs, ys);
    let stride = 2 * ncol - 1;
    let cell = |row: usize, col: usize| &grid[row * stride + 2 * col];

    // common y range per row across the columns, so the columns are comparable. The range is
    // set by the converged points of every series (unconverged states can be anything, up to
    // CD = 1e19); values beyond it are drawn clamped onto the axis edge.
    let alpha_max = families
        .iter()
        .flat_map(|(_, rs)| rs.iter().flat_map(|r| r.points.iter().map(|p| p.alpha_deg)))
        .fold(0.0_f64, f64::max);
    let ranges: Vec<(f64, f64)> = ROWS
        .iter()
        .map(|row| {
            let vals = families.iter().flat_map(|(_, rs)| {
                rs.iter().flat_map(|r| {
                    r.points
                        .iter()
                        .filter(|p| p.converged)
                        .map(|p| (row.get)(p))
                        .filter(|v| v.is_finite() && (!row.log || *v > 0.0))
                })
            });
            let (lo, hi) = vals.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
            if row.log {
                (lo / 10.0, hi * 1e4)
            } else {
                let pad = 0.5 * (hi - lo).max(1e-12);
                (lo - pad, hi + pad)
            }
        })
        .collect();

    for (c, (family, results)) in families.iter().enumerate() {
        draw_legend(cell(0, c), family, results)?;
        for (r, row) in ROWS.iter().enumerate() {
            draw_cell(cell(r + 1, c), row, ranges[r], alpha_max, results, r == ROWS.len() - 1)?;
        }
    }
    root.present()?;
    Ok(())
}

/// Stack order and colour: the base (magnitude 0) on top in black; then levels by increasing
/// magnitude, lightest nearest the base
fn ordered(results: &[RunResult]) -> Vec<(&RunResult, RGBColor, bool)> {
    let mut perturbed: Vec<&RunResult> = results.iter().filter(|r| r.level.magnitude > 0.0).collect();
    perturbed.sort_by(|a, b| a.level.magnitude.partial_cmp(&b.level.magnitude).unwrap());
    let n = perturbed.len().max(2) - 1;
    let mut out: Vec<(&RunResult, RGBColor, bool)> = perturbed
        .iter()
        .enumerate()
        .map(|(i, r)| (*r, ylgnbu(i as f64 / n as f64), false))
        .collect();
    if let Some(base) = results.iter().find(|r| r.level.magnitude == 0.0) {
        out.insert(0, (base, BLACK, true));
    }
    out
}

fn draw_legend<DB: DrawingBackend>(
    area: &Area<DB>,
    family: &Family,
    results: &[RunResult],
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    Label::new().t(family.title()).draw(
        area,
        (pxi(Y_LABEL_AREA_PT + MARGIN_PT), pxi(LEGEND_TITLE_PT - 2.0)),
        AXIS_LABEL_PT,
        false,
    )?;
    let entries = ordered(results);
    let rows = legend_rows();
    let (w, _) = area.dim_in_pixel();
    let col_w = (w as i32 - pxi(Y_LABEL_AREA_PT + MARGIN_PT)) / LEGEND_COLS as i32;
    for (i, (r, colour, is_base)) in entries.iter().enumerate() {
        // entry 0 is the base on its own row; the rest fill the columns below it
        let (col, row) = if i == 0 {
            (0, 0)
        } else {
            ((i - 1) / (rows - 1), 1 + (i - 1) % (rows - 1))
        };
        let x = pxi(Y_LABEL_AREA_PT + MARGIN_PT) + col as i32 * col_w;
        let y = pxi(LEGEND_TITLE_PT + LEGEND_ROW_PT * (row as f64 + 0.5) + 2.0);
        let len = pxi(14.0);
        let style = colour.stroke_width(if *is_base { 2 } else { 1 });
        area.draw(&PathElement::new(vec![(x, y), (x + len, y)], style))?;
        if *is_base {
            area.draw(&Circle::new((x + len / 2, y), pxi(MARKER_PT), colour.filled()))?;
        }
        let mut label = level_label(&r.level);
        if r.ending != Ending::Completed {
            label = label.t(" (hung)");
        }
        label.draw(area, (x + len + pxi(3.0), y + pxi(3.0)), LEGEND_PT, false)?;
    }
    Ok(())
}

fn draw_cell<DB: DrawingBackend>(
    area: &Area<DB>,
    row: &Row,
    (lo, hi): (f64, f64),
    alpha_max: f64,
    results: &[RunResult],
    bottom: bool,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let x_range = -0.5..alpha_max + 0.5;
    let (lo_axis, hi_axis) = if row.log { (lo.log10(), hi.log10()) } else { (lo, hi) };
    let mut chart = ChartBuilder::on(area)
        .margin(pxu(MARGIN_PT))
        .margin_right(pxu(MARGIN_PT + 4.0))
        .x_label_area_size(pxu(if bottom { X_LABEL_AREA_PT } else { X_STUB_PT }))
        .y_label_area_size(pxu(Y_LABEL_AREA_PT))
        .build_cartesian_2d(x_range, lo_axis..hi_axis)?;
    let mut mesh = chart.configure_mesh();
    mesh.label_style(figure::text(TICK_PT))
        .axis_style(BLACK.stroke_width(1))
        .light_line_style(TRANSPARENT)
        .bold_line_style(RGBColor(228, 228, 228))
        .x_labels(6)
        .y_labels(4);
    if !bottom {
        mesh.x_label_formatter(&|_| String::new());
    }
    mesh.draw()?;

    // the stack, drawn largest perturbation first so the base lands on top
    let entries = ordered(results);
    for (r, colour, is_base) in entries.iter().rev() {
        let value = |p: &Point| {
            let v = (row.get)(p);
            if row.log {
                v.log10().clamp(lo_axis, hi_axis)
            } else {
                v.clamp(lo_axis, hi_axis)
            }
        };
        let ok = |p: &Point| {
            let v = (row.get)(p);
            if row.log {
                v > 0.0 && v.is_finite()
            } else {
                !v.is_nan()
            }
        };
        let pts: Vec<(f64, f64)> = r
            .points
            .iter()
            .filter(|p| ok(p))
            .map(|p| (p.alpha_deg, value(p)))
            .collect();
        let unc: Vec<(f64, f64)> = r
            .points
            .iter()
            .filter(|p| !p.converged && ok(p))
            .map(|p| (p.alpha_deg, value(p)))
            .collect();
        series(&mut chart, pts, unc, *colour, *is_base)?;
    }

    // axis descriptions, centred on the plotting area
    let (ox, oy) = area.get_base_pixel();
    let (xr, yr) = chart.plotting_area().get_pixel_range();
    let (cx, cy) = ((xr.start + xr.end) / 2 - ox, (yr.start + yr.end) / 2 - oy);
    if bottom {
        let (_, h) = area.dim_in_pixel();
        Label::new().i("α").t(" (°)").draw_centred(
            area,
            (cx, h as i32 - pxi(MARGIN_PT + 1.0)),
            AXIS_LABEL_PT,
            false,
        )?;
    }
    (row.label)().draw_centred(
        area,
        (pxi(MARGIN_PT) + Label::height(AXIS_LABEL_PT), cy),
        AXIS_LABEL_PT,
        true,
    )?;
    Ok(())
}

fn series<DB: DrawingBackend>(
    chart: &mut Chart<DB>,
    points: Vec<(f64, f64)>,
    unconverged: Vec<(f64, f64)>,
    colour: RGBColor,
    is_base: bool,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let style = colour.stroke_width(if is_base { 2 } else { 1 });
    chart.draw_series(LineSeries::new(points.clone(), style))?;
    if is_base {
        chart.draw_series(points.iter().map(|&p| Circle::new(p, pxi(MARKER_PT), colour.filled())))?;
    }
    chart.draw_series(unconverged.iter().map(|&p| {
        Circle::new(
            p,
            pxi(OPEN_MARKER_PT),
            ShapeStyle {
                color: colour.to_rgba(),
                filled: false,
                stroke_width: 1,
            },
        )
    }))?;
    Ok(())
}
