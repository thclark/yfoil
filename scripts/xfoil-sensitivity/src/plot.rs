//! The seven-row plot: one column per input family, base case in front, perturbation levels
//! stacked behind it coloured by perturbation size (smallest lightest, drawn nearest the base).

use crate::run::{Point, RunResult};
use crate::Family;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;
use std::path::Path;

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
    label: &'static str,
    log: bool,
    get: fn(&Point) -> f64,
}

const ROWS: [Row; 7] = [
    Row {
        label: "CL",
        log: false,
        get: |p| p.cl,
    },
    Row {
        label: "CD",
        log: false,
        get: |p| p.cd,
    },
    Row {
        label: "δ* at TE, upper",
        log: false,
        get: |p| p.dstar_te_upper,
    },
    Row {
        label: "H at TE, upper",
        log: false,
        get: |p| p.h_te_upper,
    },
    Row {
        label: "x/c transition, upper",
        log: false,
        get: |p| p.xtr_upper,
    },
    Row {
        label: "x/c TE separation, upper",
        log: false,
        get: |p| p.x_sep_upper,
    },
    Row {
        label: "RMSBL, final iteration",
        log: true,
        get: |p| p.rmsbl,
    },
];

const COL_W: u32 = 520;
const ROW_H: u32 = 260;
const LEGEND_H: u32 = 250;
const TITLE_H: u32 = 40;
const FOOT_H: u32 = 24;

/// Draw one column per family (a single family gives the per-family figure, all three the sheet)
pub fn draw_families(
    path: &Path,
    foil: &str,
    families: &[(Family, Vec<RunResult>)],
) -> Result<(), Box<dyn std::error::Error>> {
    let ncol = families.len() as u32;
    let width = COL_W * ncol + 20;
    let height = TITLE_H + LEGEND_H + ROW_H * ROWS.len() as u32 + FOOT_H;
    let root = SVGBackend::new(path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    let (title_area, rest) = root.split_vertically(TITLE_H);
    title_area.draw(&Text::new(
        format!(
            "XFOIL 6.99 DP input sensitivity — NACA {foil}, Re {:e}, M {}, Ncrit {}, ITER {}, ALFA per point",
            crate::REYNOLDS,
            crate::MACH,
            crate::NCRIT,
            crate::ITER
        ),
        (12, 12),
        ("sans-serif", 18).into_font(),
    ))?;
    let (legend_area, rest) = rest.split_vertically(LEGEND_H);
    let (grid, _foot) = rest.split_vertically(ROW_H * ROWS.len() as u32);
    let legend_cols = legend_area.split_evenly((1, ncol as usize));
    let cells = grid.split_evenly((ROWS.len(), ncol as usize));

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
    root.draw(&Text::new(
        "Axes span the converged points; off-scale values (unconverged states) are drawn on the axis edge. Open circles: RMSBL did not reach 1e-4 in ITER iterations.",
        (12, (height - 16) as i32),
        ("sans-serif", 12).into_font(),
    ))?;

    for (c, (family, results)) in families.iter().enumerate() {
        draw_legend(&legend_cols[c], family, results)?;
        for (r, row) in ROWS.iter().enumerate() {
            let cell = &cells[r * ncol as usize + c];
            draw_cell(cell, row, ranges[r], alpha_max, results, r == ROWS.len() - 1)?;
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
    area: &DrawingArea<DB, plotters::coord::Shift>,
    family: &Family,
    results: &[RunResult],
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    area.draw(&Text::new(
        family.title().to_string(),
        (30, 6),
        ("sans-serif", 16).into_font(),
    ))?;
    let entries = ordered(results);
    let line_h = ((LEGEND_H - 30) as f64 / entries.len().max(1) as f64).min(16.0) as i32;
    for (i, (r, colour, is_base)) in entries.iter().enumerate() {
        let y = 30 + i as i32 * line_h;
        let style = ShapeStyle {
            color: colour.to_rgba(),
            filled: true,
            stroke_width: if *is_base { 2 } else { 1 },
        };
        area.draw(&PathElement::new(
            vec![(30, y + line_h / 2), (60, y + line_h / 2)],
            style,
        ))?;
        if *is_base {
            area.draw(&Circle::new((45, y + line_h / 2), 3, style))?;
        }
        let note = if r.ending == crate::run::Ending::Completed {
            ""
        } else {
            " (XFOIL hung; killed)"
        };
        area.draw(&Text::new(
            format!("{}{note}", r.level.label),
            (68, y + 2),
            ("sans-serif", 12).into_font(),
        ))?;
    }
    Ok(())
}

fn draw_cell<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
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
    // the stack, drawn largest perturbation first so the base lands on top
    let entries = ordered(results);
    let series = |chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
                  points: Vec<(f64, f64)>,
                  unconverged: Vec<(f64, f64)>,
                  colour: RGBColor,
                  is_base: bool|
     -> Result<(), Box<dyn std::error::Error>> {
        let style = colour.stroke_width(if is_base { 2 } else { 1 });
        chart.draw_series(LineSeries::new(points.clone(), style))?;
        if is_base {
            chart.draw_series(points.iter().map(|&p| Circle::new(p, 3, colour.filled())))?;
        }
        chart.draw_series(unconverged.iter().map(|&p| {
            Circle::new(
                p,
                4,
                ShapeStyle {
                    color: colour.to_rgba(),
                    filled: false,
                    stroke_width: 1,
                },
            )
        }))?;
        Ok(())
    };
    let mut builder = ChartBuilder::on(area);
    builder
        .margin(8)
        .margin_left(12)
        .x_label_area_size(if bottom { 32 } else { 18 })
        .y_label_area_size(64);
    if row.log {
        // log axis: transform to log10 and label accordingly
        let mut chart = builder.build_cartesian_2d(x_range, lo.log10()..hi.log10())?;
        chart
            .configure_mesh()
            .x_desc(if bottom { "α (deg)" } else { "" })
            .y_desc(format!("log10 {}", row.label))
            .label_style(("sans-serif", 11))
            .axis_desc_style(("sans-serif", 12))
            .light_line_style(RGBColor(235, 235, 235))
            .draw()?;
        for (r, colour, is_base) in entries.iter().rev() {
            let val = |p: &Point| (row.get)(p).log10().clamp(lo.log10(), hi.log10());
            let ok = |p: &Point| (row.get)(p) > 0.0 && (row.get)(p).is_finite();
            let pts: Vec<(f64, f64)> = r
                .points
                .iter()
                .filter(|p| ok(p))
                .map(|p| (p.alpha_deg, val(p)))
                .collect();
            let unc: Vec<(f64, f64)> = r
                .points
                .iter()
                .filter(|p| !p.converged && ok(p))
                .map(|p| (p.alpha_deg, val(p)))
                .collect();
            series(&mut chart, pts, unc, *colour, *is_base)?;
        }
    } else {
        let mut chart = builder.build_cartesian_2d(x_range, lo..hi)?;
        chart
            .configure_mesh()
            .x_desc(if bottom { "α (deg)" } else { "" })
            .y_desc(row.label)
            .label_style(("sans-serif", 11))
            .axis_desc_style(("sans-serif", 12))
            .light_line_style(RGBColor(235, 235, 235))
            .draw()?;
        for (r, colour, is_base) in entries.iter().rev() {
            let val = |p: &Point| (row.get)(p).clamp(lo, hi);
            let ok = |p: &Point| !(row.get)(p).is_nan();
            let pts: Vec<(f64, f64)> = r
                .points
                .iter()
                .filter(|p| ok(p))
                .map(|p| (p.alpha_deg, val(p)))
                .collect();
            let unc: Vec<(f64, f64)> = r
                .points
                .iter()
                .filter(|p| !p.converged && ok(p))
                .map(|p| (p.alpha_deg, val(p)))
                .collect();
            series(&mut chart, pts, unc, *colour, *is_base)?;
        }
    }
    Ok(())
}
