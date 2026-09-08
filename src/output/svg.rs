//! A small hand-written SVG canvas with a true equal-aspect data transform and sub-pixel
//! coordinates. plotters rounds to integer pixels in both of its backends, which is too coarse
//! for boundary-layer offsets of a few pixels; every geometry-shaped plot is drawn here instead,
//! and PNG output is the same SVG rasterised (one drawing path).

use std::fmt::Write as _;
use std::path::Path;

use crate::output::plot::PlotError;

/// Pixel margins around the plot area
#[derive(Debug, Clone, Copy)]
pub(crate) struct Margins {
    pub left: f64,
    pub right: f64,
    pub top: f64,
    pub bottom: f64,
}

/// Stroke style of a line
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LineStyle {
    pub rgb: (u8, u8, u8),
    pub width: f64,
    /// SVG `stroke-dasharray`, `None` for solid
    pub dasharray: Option<String>,
}

/// Open marker shapes for point locations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MarkerShape {
    Circle,
    Square,
    Triangle,
    Diamond,
    TriangleDown,
}

impl MarkerShape {
    /// The shape centred on pixel `(cx, cy)` with half-size `r`, unfilled; stroke attributes
    /// are appended by the caller
    fn element(self, cx: f64, cy: f64, r: f64, stroke_attrs: &str) -> String {
        match self {
            MarkerShape::Circle => {
                format!(r#"<circle cx="{cx:.3}" cy="{cy:.3}" r="{r}" fill="none"{stroke_attrs}/>"#)
            }
            MarkerShape::Square => format!(
                r#"<rect x="{:.3}" y="{:.3}" width="{}" height="{}" fill="none"{stroke_attrs}/>"#,
                cx - r,
                cy - r,
                2.0 * r,
                2.0 * r
            ),
            MarkerShape::Triangle => format!(
                r#"<polygon points="{:.3},{:.3} {:.3},{:.3} {:.3},{:.3}" fill="none"{stroke_attrs}/>"#,
                cx,
                cy - 1.15 * r,
                cx - 1.15 * r,
                cy + 0.85 * r,
                cx + 1.15 * r,
                cy + 0.85 * r
            ),
            MarkerShape::TriangleDown => format!(
                r#"<polygon points="{:.3},{:.3} {:.3},{:.3} {:.3},{:.3}" fill="none"{stroke_attrs}/>"#,
                cx,
                cy + 1.15 * r,
                cx - 1.15 * r,
                cy - 0.85 * r,
                cx + 1.15 * r,
                cy - 0.85 * r
            ),
            MarkerShape::Diamond => format!(
                r#"<polygon points="{:.3},{:.3} {:.3},{:.3} {:.3},{:.3} {:.3},{:.3}" fill="none"{stroke_attrs}/>"#,
                cx,
                cy - 1.25 * r,
                cx + 1.25 * r,
                cy,
                cx,
                cy + 1.25 * r,
                cx - 1.25 * r,
                cy
            ),
        }
    }
}

/// The key drawn next to a legend entry
#[derive(Debug, Clone)]
pub(crate) enum LegendKey {
    Line(LineStyle),
    Marker(MarkerShape, (u8, u8, u8)),
}

/// A legend entry: a key and its text
#[derive(Debug, Clone)]
pub(crate) struct LegendEntry {
    pub text: String,
    pub key: LegendKey,
}

/// Where a legend sits inside the plot area
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LegendAnchor {
    BottomLeft,
    BottomRight,
}

/// `rgb(r,g,b)`
pub(crate) fn rgb(c: (u8, u8, u8)) -> String {
    format!("rgb({},{},{})", c.0, c.1, c.2)
}

/// Escape text for an XML text node
pub(crate) fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// A nice step size for axis ticks (1, 2, 5 × decade)
pub(crate) fn nice_step(range: f64, target_ticks: usize) -> f64 {
    let rough_step = range / target_ticks as f64;
    let magnitude = 10.0_f64.powf(rough_step.log10().floor());
    let normalized = rough_step / magnitude;

    let nice = if normalized < 1.5 {
        1.0
    } else if normalized < 3.0 {
        2.0
    } else if normalized < 7.0 {
        5.0
    } else {
        10.0
    };

    nice * magnitude
}

/// Data-to-pixel canvas. Data coordinates are mapped at one scale in x and y (equal aspect): the
/// data box is fitted inside the margin-adjusted plot area and centred.
pub(crate) struct Canvas {
    pub width: f64,
    pub height: f64,
    pub margins: Margins,
    /// Data coordinates at the centre of the plot area
    cx: f64,
    cy: f64,
    /// Pixels per data unit
    scale: f64,
    buf: String,
}

impl Canvas {
    /// Fit the data box `(x_min, x_max, y_min, y_max)` inside the plot area at one scale
    pub fn equal_aspect(width: u32, height: u32, margins: Margins, bounds: (f64, f64, f64, f64)) -> Self {
        let (x_min, x_max, y_min, y_max) = bounds;
        let width = width as f64;
        let height = height as f64;
        let plot_w = (width - margins.left - margins.right).max(1.0);
        let plot_h = (height - margins.top - margins.bottom).max(1.0);
        let dx = (x_max - x_min).max(1e-12);
        let dy = (y_max - y_min).max(1e-12);
        let scale = (plot_w / dx).min(plot_h / dy);
        Self {
            width,
            height,
            margins,
            cx: 0.5 * (x_min + x_max),
            cy: 0.5 * (y_min + y_max),
            scale,
            buf: String::new(),
        }
    }

    /// Plot area in pixels: (left, top, width, height)
    pub fn plot_area(&self) -> (f64, f64, f64, f64) {
        (
            self.margins.left,
            self.margins.top,
            self.width - self.margins.left - self.margins.right,
            self.height - self.margins.top - self.margins.bottom,
        )
    }

    /// Data bounds visible in the plot area: (x_min, x_max, y_min, y_max)
    pub fn visible(&self) -> (f64, f64, f64, f64) {
        let (_, _, w, h) = self.plot_area();
        let hx = 0.5 * w / self.scale;
        let hy = 0.5 * h / self.scale;
        (self.cx - hx, self.cx + hx, self.cy - hy, self.cy + hy)
    }

    /// Data → pixel
    pub fn to_svg(&self, x: f64, y: f64) -> (f64, f64) {
        let (left, top, w, h) = self.plot_area();
        (
            left + 0.5 * w + (x - self.cx) * self.scale,
            top + 0.5 * h - (y - self.cy) * self.scale,
        )
    }

    pub fn header(&mut self, background: (u8, u8, u8)) {
        let _ = writeln!(
            self.buf,
            r#"<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg">"#,
            w = self.width,
            h = self.height
        );
        let _ = writeln!(
            self.buf,
            r#"<rect width="100%" height="100%" fill="{}"/>"#,
            rgb(background)
        );
    }

    pub fn title(&mut self, text: &str) {
        let _ = writeln!(
            self.buf,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="16">{}</text>"#,
            0.5 * self.width,
            0.65 * self.margins.top,
            escape_xml(text)
        );
    }

    /// Grid lines with numeric tick labels, and the axis titles
    pub fn grid(&mut self, x_label: &str, y_label: &str) {
        let (left, top, w, h) = self.plot_area();
        let (vx0, vx1, vy0, vy1) = self.visible();
        let sx = nice_step(vx1 - vx0, 8);
        let sy = nice_step(vy1 - vy0, 5);
        let _ = writeln!(self.buf, r##"<g stroke="#DDDDDD" stroke-width="0.5">"##);
        let mut ticks_x = Vec::new();
        let mut k = (vx0 / sx).ceil() as i64;
        while (k as f64) * sx <= vx1 + 1e-9 * sx {
            let x = k as f64 * sx;
            let (px, _) = self.to_svg(x, 0.0);
            let _ = writeln!(
                self.buf,
                r#"<line x1="{px:.2}" y1="{:.2}" x2="{px:.2}" y2="{:.2}"/>"#,
                top,
                top + h
            );
            ticks_x.push((x, px));
            k += 1;
        }
        let mut ticks_y = Vec::new();
        let mut k = (vy0 / sy).ceil() as i64;
        while (k as f64) * sy <= vy1 + 1e-9 * sy {
            let y = k as f64 * sy;
            let (_, py) = self.to_svg(0.0, y);
            let _ = writeln!(
                self.buf,
                r#"<line x1="{:.2}" y1="{py:.2}" x2="{:.2}" y2="{py:.2}"/>"#,
                left,
                left + w
            );
            ticks_y.push((y, py));
            k += 1;
        }
        let _ = writeln!(self.buf, "</g>");
        let _ = writeln!(
            self.buf,
            r##"<g font-family="sans-serif" font-size="10" fill="#333333">"##
        );
        for (x, px) in ticks_x {
            let _ = writeln!(
                self.buf,
                r#"<text x="{px:.1}" y="{:.1}" text-anchor="middle">{}</text>"#,
                top + h + 14.0,
                tick_label(x, sx)
            );
        }
        for (y, py) in ticks_y {
            let _ = writeln!(
                self.buf,
                r#"<text x="{:.1}" y="{:.1}" text-anchor="end" dominant-baseline="middle">{}</text>"#,
                left - 6.0,
                py,
                tick_label(y, sy)
            );
        }
        let _ = writeln!(self.buf, "</g>");
        let _ = writeln!(
            self.buf,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">{}</text>"#,
            left + 0.5 * w,
            self.height - 8.0,
            escape_xml(x_label)
        );
        let yl = top + 0.5 * h;
        let _ = writeln!(
            self.buf,
            r#"<text x="14" y="{yl:.1}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 14 {yl:.1})">{}</text>"#,
            escape_xml(y_label)
        );
        let _ = writeln!(
            self.buf,
            r##"<rect x="{left:.1}" y="{top:.1}" width="{w:.1}" height="{h:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##
        );
    }

    /// A polyline through data points
    pub fn polyline(&mut self, pts: &[(f64, f64)], style: &LineStyle) {
        if pts.len() < 2 {
            return;
        }
        let _ = write!(
            self.buf,
            r#"<polyline fill="none" stroke="{}" stroke-width="{}"{} points=""#,
            rgb(style.rgb),
            style.width,
            dash_attr(style)
        );
        for (k, &(x, y)) in pts.iter().enumerate() {
            let (px, py) = self.to_svg(x, y);
            if k > 0 {
                self.buf.push(' ');
            }
            let _ = write!(self.buf, "{px:.6},{py:.6}");
        }
        let _ = writeln!(self.buf, r#""/>"#);
    }

    /// A filled circle of pixel radius `r` at a data point
    pub fn circle(&mut self, at: (f64, f64), r: f64, fill: (u8, u8, u8)) {
        let (px, py) = self.to_svg(at.0, at.1);
        let _ = writeln!(
            self.buf,
            r#"<circle cx="{px:.6}" cy="{py:.6}" r="{r}" fill="{}"/>"#,
            rgb(fill)
        );
    }

    /// Open a group with a shared stroke (for many notches)
    pub fn group_start(&mut self, style: &LineStyle) {
        let _ = writeln!(
            self.buf,
            r#"<g stroke="{}" stroke-width="{}"{}>"#,
            rgb(style.rgb),
            style.width,
            dash_attr(style)
        );
    }

    /// A line inside an open group (no per-line stroke attributes)
    pub fn group_line(&mut self, a: (f64, f64), b: (f64, f64)) {
        let (x1, y1) = self.to_svg(a.0, a.1);
        let (x2, y2) = self.to_svg(b.0, b.1);
        let _ = writeln!(
            self.buf,
            r#"<line x1="{x1:.6}" y1="{y1:.6}" x2="{x2:.6}" y2="{y2:.6}"/>"#
        );
    }

    pub fn group_end(&mut self) {
        let _ = writeln!(self.buf, "</g>");
    }

    /// An open marker of pixel half-size `r` at a data point
    pub fn marker(&mut self, at: (f64, f64), shape: MarkerShape, r: f64, rgb_: (u8, u8, u8), width: f64) {
        let (px, py) = self.to_svg(at.0, at.1);
        let attrs = format!(r#" stroke="{}" stroke-width="{width}""#, rgb(rgb_));
        let _ = writeln!(self.buf, "{}", shape.element(px, py, r, &attrs));
    }

    /// Legend in a bottom corner of the plot area, laid out in columns when the entries would
    /// take more than the empty band below the foil
    pub fn legend(&mut self, entries: &[LegendEntry], anchor: LegendAnchor) {
        if entries.is_empty() {
            return;
        }
        let (left, top, w, h) = self.plot_area();
        let row_h = 16.0;
        let key_w = 30.0;
        let pad = 8.0;
        // keep the legend in the empty band below the foil: at most ~30% of the plot height
        let max_rows = (((0.3 * h - pad) / row_h).floor() as usize).max(1);
        let rows = entries.len().min(max_rows);
        let cols = entries.len().div_ceil(rows);
        let col_w: Vec<f64> = (0..cols)
            .map(|c| {
                let widest = entries
                    .iter()
                    .skip(c * rows)
                    .take(rows)
                    .map(|e| e.text.chars().count())
                    .max()
                    .unwrap_or(0);
                key_w + 8.0 + 6.6 * widest as f64 + 12.0
            })
            .collect();
        let box_w: f64 = col_w.iter().sum::<f64>() + pad;
        let box_h = rows as f64 * row_h + pad;
        let x0 = match anchor {
            LegendAnchor::BottomLeft => left + 10.0,
            LegendAnchor::BottomRight => left + w - 10.0 - box_w,
        };
        let y0 = top + h - 10.0 - box_h;
        let _ = writeln!(
            self.buf,
            r##"<rect x="{x0:.1}" y="{y0:.1}" width="{box_w:.1}" height="{box_h:.1}" fill="white" fill-opacity="0.85" stroke="#999999" stroke-width="0.5"/>"##
        );
        let _ = writeln!(
            self.buf,
            r##"<g font-family="sans-serif" font-size="12" fill="#222222">"##
        );
        let mut cx = x0 + pad;
        for (c, cw) in col_w.iter().enumerate() {
            for (r, e) in entries.iter().skip(c * rows).take(rows).enumerate() {
                let cy = y0 + pad + (r as f64 + 0.5) * row_h;
                match &e.key {
                    LegendKey::Line(style) => {
                        let _ = writeln!(
                            self.buf,
                            r#"<line x1="{:.1}" y1="{cy:.1}" x2="{:.1}" y2="{cy:.1}" stroke="{}" stroke-width="{}"{}/>"#,
                            cx,
                            cx + key_w,
                            rgb(style.rgb),
                            style.width,
                            dash_attr(style)
                        );
                    }
                    LegendKey::Marker(shape, colour) => {
                        let attrs = format!(r#" stroke="{}" stroke-width="1.5""#, rgb(*colour));
                        let _ = writeln!(self.buf, "{}", shape.element(cx + 0.5 * key_w, cy, 4.5, &attrs));
                    }
                }
                let _ = writeln!(
                    self.buf,
                    r#"<text x="{:.1}" y="{cy:.1}" dominant-baseline="middle">{}</text>"#,
                    cx + key_w + 8.0,
                    escape_xml(&e.text)
                );
            }
            cx += cw;
        }
        let _ = writeln!(self.buf, "</g>");
    }

    /// Close the document and return it
    pub fn finish(mut self) -> String {
        self.buf.push_str("</svg>\n");
        self.buf
    }
}

fn dash_attr(style: &LineStyle) -> String {
    match &style.dasharray {
        Some(d) => format!(r#" stroke-dasharray="{d}""#),
        None => String::new(),
    }
}

fn tick_label(v: f64, step: f64) -> String {
    let decimals = if step >= 1.0 {
        0
    } else {
        (-step.log10().floor()) as usize
    };
    let s = format!("{:.*}", decimals, v);
    if s == "-0" || s.trim_start_matches('-').chars().all(|c| c == '0' || c == '.') {
        format!("{:.*}", decimals, 0.0)
    } else {
        s
    }
}

/// Write an SVG document to a file
pub(crate) fn write_svg(path: &Path, svg: &str) -> Result<(), PlotError> {
    std::fs::write(path, svg)?;
    Ok(())
}

/// Rasterise an SVG document to a PNG of the document's own size
pub(crate) fn write_png(path: &Path, svg: &str) -> Result<(), PlotError> {
    let mut opt = resvg::usvg::Options::default();
    opt.fontdb_mut().load_system_fonts();
    let tree = resvg::usvg::Tree::from_str(svg, &opt).map_err(|e| PlotError::Drawing(e.to_string()))?;
    let size = tree.size().to_int_size();
    let mut pixmap = resvg::tiny_skia::Pixmap::new(size.width(), size.height())
        .ok_or_else(|| PlotError::Drawing("zero-sized image".to_string()))?;
    resvg::render(&tree, resvg::tiny_skia::Transform::default(), &mut pixmap.as_mut());
    pixmap.save_png(path).map_err(|e| PlotError::Drawing(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_aspect_uses_one_scale_and_centres() {
        let m = Margins {
            left: 50.0,
            right: 20.0,
            top: 30.0,
            bottom: 40.0,
        };
        // a wide, flat data box on a canvas whose plot area is 1000×400: x-limited
        let c = Canvas::equal_aspect(1070, 470, m, (-0.05, 1.05, -0.1, 0.1));
        let px_per_unit = c.to_svg(1.0, 0.0).0 - c.to_svg(0.0, 0.0).0;
        assert!((px_per_unit - 1000.0 / 1.1).abs() < 1e-9);
        let (x0, y0) = c.to_svg(0.5, 0.0);
        assert!((x0 - (50.0 + 500.0)).abs() < 1e-9);
        assert!((y0 - (30.0 + 200.0)).abs() < 1e-9);
        // one data unit is the same number of pixels in x and y
        let (x1, y1) = c.to_svg(0.6, 0.1);
        assert!(((x1 - x0) - (y0 - y1)).abs() < 1e-9);
        // a tall box is y-limited
        let c = Canvas::equal_aspect(1070, 470, m, (0.0, 1.0, -1.0, 1.0));
        let px_per_unit = c.to_svg(1.0, 0.0).0 - c.to_svg(0.0, 0.0).0;
        assert!((px_per_unit - 200.0).abs() < 1e-9);
    }

    #[test]
    fn nice_steps() {
        assert_eq!(nice_step(1.1, 8), 0.1);
        assert_eq!(nice_step(0.3, 5), 0.05);
        assert_eq!(nice_step(20.0, 5), 5.0);
        assert_eq!(tick_label(0.30000000000000004, 0.1), "0.3");
        assert_eq!(tick_label(-0.0, 0.1), "0.0");
        assert_eq!(tick_label(2.0, 1.0), "2");
    }

    #[test]
    fn legend_lays_out_columns() {
        let m = Margins {
            left: 10.0,
            right: 10.0,
            top: 10.0,
            bottom: 10.0,
        };
        let mut c = Canvas::equal_aspect(300, 80, m, (0.0, 1.0, 0.0, 0.2));
        let entries: Vec<LegendEntry> = (0..6)
            .map(|i| LegendEntry {
                text: format!("e{i}"),
                key: LegendKey::Line(LineStyle {
                    rgb: (0, 0, 0),
                    width: 1.0,
                    dasharray: None,
                }),
            })
            .collect();
        c.legend(&entries, LegendAnchor::BottomLeft);
        c.legend(
            &[LegendEntry {
                text: "m".to_string(),
                key: LegendKey::Marker(MarkerShape::Diamond, (1, 2, 3)),
            }],
            LegendAnchor::BottomRight,
        );
        let svg = c.finish();
        assert_eq!(svg.matches("<line").count(), 6);
        assert!(svg.contains("e5"));
        // 6 entries at 16px do not fit a 60px plot area: two or more columns
        assert!(svg.matches("<text").count() == 7);
        assert!(svg.contains(r#"<polygon points="#));
        assert!(svg.contains(r#"stroke="rgb(1,2,3)""#));
    }
}
