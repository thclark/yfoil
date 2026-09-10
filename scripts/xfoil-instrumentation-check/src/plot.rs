//! The 2×2 figure: CL, CD and CM against α (pristine red with crosses, instrumented blue with
//! circles) and the upper-surface trailing-edge boundary layer against α (one colour per
//! variable, pristine solid, instrumented dashed — the caption says so; the figure's legends
//! name the builds and the variables). Drawn from `summary.json` only, at the
//! physical size and in the font of the target document (`figure.rs`).

use figure_style::{self as figure, pxi, pxu, Label, AXIS_LABEL_PT, LEGEND_PT, TEXT_WIDTH_PT, TICK_PT};
use crate::run::Point;
use crate::summary::Summary;
use crate::Build;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;
use plotters::series::DashedLineSeries;
use std::path::Path;

const RED: RGBColor = RGBColor(214, 39, 40);
const BLUE: RGBColor = RGBColor(31, 90, 214);

/// Colours of the BL variables (plotly's D3 category set, minus the red and blue used above)
const BL_COLOURS: [RGBColor; 4] = [
    RGBColor(44, 160, 44),
    RGBColor(148, 103, 189),
    RGBColor(255, 127, 14),
    RGBColor(23, 190, 207),
];

/// A plotted quantity: its label and accessor
struct Var {
    label: fn() -> Label,
    get: fn(&Point) -> f64,
}

/// Upper-surface trailing-edge variables of the fourth panel. Cf is omitted: it is ≤ 0 once
/// the trailing edge separates, and the panel is on a logarithmic axis.
const BL_VARS: [Var; 4] = [
    Var {
        label: || Label::new().i("U").subi("e").t("/").i("V").subi("∞"),
        get: |p| p.te_upper.ue,
    },
    Var {
        label: || Label::new().i("δ").sup("*"),
        get: |p| p.te_upper.dstar,
    },
    Var {
        label: || Label::new().i("θ"),
        get: |p| p.te_upper.theta,
    },
    Var {
        label: || Label::new().i("H").subi("k"),
        get: |p| p.te_upper.hk,
    },
];

/// Layout, pt
const LEGEND_ROW_PT: f64 = 11.0;
const CELL_H_PT: f64 = 172.0;
const Y_LABEL_AREA_PT: f64 = 34.0;
/// Gap between the two columns, so a right-hand axis label reads with its own chart
const GUTTER_PT: f64 = 16.0;
/// Room under the plot for the tick values, a gap, and the axis label
const X_LABEL_AREA_PT: f64 = 28.0;
const MARGIN_PT: f64 = 3.0;
const MARKER_PT: f64 = 2.6;

fn colour(b: Build) -> RGBColor {
    match b {
        Build::Pristine => RED,
        Build::Instrumented => BLUE,
    }
}

fn open(c: RGBColor) -> ShapeStyle {
    ShapeStyle {
        color: c.to_rgba(),
        filled: false,
        stroke_width: 1,
    }
}

type Chart<'a, DB> = ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>;
type Area<DB> = DrawingArea<DB, plotters::coord::Shift>;

/// Markers of a build: crosses for the pristine build, open circles for the instrumented one
fn marker<DB: DrawingBackend>(
    chart: &mut Chart<DB>,
    b: Build,
    pts: &[(f64, f64)],
    c: RGBColor,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let r = pxi(MARKER_PT);
    match b {
        Build::Pristine => chart.draw_series(pts.iter().map(|&p| Cross::new(p, r, c.stroke_width(1))))?,
        Build::Instrumented => chart.draw_series(pts.iter().map(|&p| Circle::new(p, r, open(c))))?,
    };
    Ok(())
}

/// A legend sample at `(x, y)` (line centre): the build's line and marker
fn sample<DB: DrawingBackend>(area: &Area<DB>, b: Build, (x, y): (i32, i32)) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let c = colour(b);
    let (len, r) = (pxi(22.0), pxi(MARKER_PT));
    area.draw(&PathElement::new(vec![(x, y), (x + len, y)], c.stroke_width(1)))?;
    match b {
        Build::Pristine => area.draw(&Cross::new((x + len / 2, y), r, c.stroke_width(1)))?,
        Build::Instrumented => area.draw(&Circle::new((x + len / 2, y), r, open(c)))?,
    }
    Ok(())
}

/// Build a cell's chart: axes without descriptions (those are drawn by `axis_labels`)
fn chart<DB: DrawingBackend>(
    cell: &Area<DB>,
    x_range: std::ops::Range<f64>,
    y_range: std::ops::Range<f64>,
) -> Result<Chart<'_, DB>, Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let mut chart = ChartBuilder::on(cell)
        .margin(pxu(MARGIN_PT))
        .margin_right(pxu(MARGIN_PT + 4.0))
        .x_label_area_size(pxu(X_LABEL_AREA_PT))
        .y_label_area_size(pxu(Y_LABEL_AREA_PT))
        .build_cartesian_2d(x_range, y_range)?;
    chart
        .configure_mesh()
        .label_style(figure::text(TICK_PT))
        .axis_style(BLACK.stroke_width(1))
        .light_line_style(TRANSPARENT)
        .bold_line_style(RGBColor(228, 228, 228))
        .x_labels(8)
        .y_labels(6)
        .draw()?;
    Ok(chart)
}

/// The axis descriptions of a cell, centred on its plotting area
fn axis_labels<DB: DrawingBackend>(
    cell: &Area<DB>,
    chart: &Chart<DB>,
    x_label: &Label,
    y_label: &Label,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let (ox, oy) = cell.get_base_pixel();
    let (xr, yr) = chart.plotting_area().get_pixel_range();
    let (cx, cy) = ((xr.start + xr.end) / 2 - ox, (yr.start + yr.end) / 2 - oy);
    let (_, h) = cell.dim_in_pixel();
    x_label.draw_centred(cell, (cx, h as i32 - pxi(MARGIN_PT + 1.0)), AXIS_LABEL_PT, false)?;
    // rotated text hangs to the left of its baseline: put the baseline one text height in
    y_label.draw_centred(
        cell,
        (pxi(MARGIN_PT) + Label::height(AXIS_LABEL_PT), cy),
        AXIS_LABEL_PT,
        true,
    )?;
    Ok(())
}

pub fn draw(path: &Path, s: &Summary) -> Result<(), Box<dyn std::error::Error>> {
    let legend_h = LEGEND_ROW_PT * s.builds.len() as f64 + 4.0;
    let width = pxu(TEXT_WIDTH_PT);
    let height = pxu(legend_h + 2.0 * CELL_H_PT);
    let root = SVGBackend::new(path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    // legend: one row per build
    let (legend, grid) = root.split_vertically(pxu(legend_h));
    for (i, b) in s.builds.iter().enumerate() {
        let y = pxi(LEGEND_ROW_PT * (i as f64 + 0.5) + 2.0);
        sample(&legend, b.build, (pxi(2.0), y))?;
        let text = format!(
            "{}: {} converged points, series to {}°",
            b.label,
            b.points.len(),
            b.series_end_alpha_deg
                .map(|a| format!("{a:.1}"))
                .unwrap_or_else(|| "—".into())
        );
        Label::new()
            .t(&text)
            .draw(&legend, (pxi(30.0), y + pxi(3.0)), LEGEND_PT, false)?;
    }

    // two columns with a gutter between them: the breakpoint grid has three columns, and the
    // middle one is left empty
    let (gw, gh) = grid.dim_in_pixel();
    let col_w = (gw as i32 - pxi(GUTTER_PT)) / 2;
    let all = grid.split_by_breakpoints([col_w, col_w + pxi(GUTTER_PT)], [gh as i32 / 2]);
    let cells = [&all[0], &all[2], &all[3], &all[5]];
    let alpha_max = s
        .builds
        .iter()
        .flat_map(|b| b.series().map(|p| p.alpha_deg))
        .fold(0.0_f64, f64::max);
    let x_range = -0.5..alpha_max + 0.5;
    let x_label = Label::new().i("α").t(" (°)");

    let coefficient: [Var; 3] = [
        Var {
            label: || Label::new().i("C").subi("L"),
            get: |p| p.cl,
        },
        Var {
            label: || Label::new().i("C").subi("D"),
            get: |p| p.cd,
        },
        Var {
            label: || Label::new().i("C").subi("M"),
            get: |p| p.cm,
        },
    ];
    for (i, Var { label, get }) in coefficient.iter().enumerate() {
        let (lo, hi) = range(s.builds.iter().flat_map(|b| b.series().map(get)));
        let mut ch = chart(cells[i], x_range.clone(), lo..hi)?;
        for b in &s.builds {
            let col = colour(b.build);
            let pts: Vec<(f64, f64)> = b.series().map(|p| (p.alpha_deg, get(p))).collect();
            ch.draw_series(LineSeries::new(pts.clone(), col.stroke_width(1)))?;
            marker(&mut ch, b.build, &pts, col)?;
        }
        axis_labels(cells[i], &ch, &x_label, &label())?;
    }

    // trailing-edge boundary layer, log10 axis, with its variable legend in a strip above it
    let (bl_legend, bl_cell) = cells[3].split_vertically(pxu(LEGEND_ROW_PT));
    let vals = s
        .builds
        .iter()
        .flat_map(|b| b.series().flat_map(|p| BL_VARS.iter().map(move |v| (v.get)(p))))
        .filter(|v| *v > 0.0)
        .map(f64::log10);
    let (lo, hi) = range(vals);
    let mut ch = chart(&bl_cell, x_range.clone(), lo..hi)?;
    let mut lx = pxi(Y_LABEL_AREA_PT + MARGIN_PT);
    let ly = pxi(LEGEND_ROW_PT / 2.0 + 1.0);
    for (k, v) in BL_VARS.iter().enumerate() {
        let col = BL_COLOURS[k];
        for b in &s.builds {
            let pts: Vec<(f64, f64)> = b
                .series()
                .map(|p| (p.alpha_deg, (v.get)(p)))
                .filter(|(_, y)| *y > 0.0)
                .map(|(a, y)| (a, y.log10()))
                .collect();
            match b.build {
                Build::Pristine => {
                    ch.draw_series(LineSeries::new(pts.clone(), col.stroke_width(2)))?;
                }
                Build::Instrumented => {
                    ch.draw_series(DashedLineSeries::new(
                        pts.clone(),
                        pxu(3.0),
                        pxu(2.5),
                        col.stroke_width(2),
                    ))?;
                }
            }
            marker(&mut ch, b.build, &pts, col)?;
        }
        let len = pxi(14.0);
        bl_legend.draw(&PathElement::new(vec![(lx, ly), (lx + len, ly)], col.stroke_width(2)))?;
        let label = (v.label)();
        label.draw(&bl_legend, (lx + len + pxi(3.0), ly + pxi(3.0)), LEGEND_PT, false)?;
        lx += len + pxi(3.0) + label.width(LEGEND_PT) + pxi(9.0);
    }
    axis_labels(
        &bl_cell,
        &ch,
        &x_label,
        &Label::new().t("log").subt("10").t(" value, upper-surface TE"),
    )?;
    root.present()?;
    Ok(())
}

fn range(vals: impl Iterator<Item = f64>) -> (f64, f64) {
    let (lo, hi) = vals
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
    if !lo.is_finite() {
        return (0.0, 1.0);
    }
    let pad = 0.06 * (hi - lo).max(1e-9);
    (lo - pad, hi + pad)
}
