//! Plotting utilities for airfoil visualization
//!
//! Uses the plotters library for PNG output and custom high-precision SVG output.

use plotters::prelude::*;
use std::io::Write;
use std::path::Path;

use crate::geometry::{Geometry, PaneledAirfoil};

/// Configuration for geometry plots
#[derive(Debug, Clone)]
pub struct GeometryPlotConfig {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// Show panel node ticks perpendicular to surface
    pub show_nodes: bool,
    /// Length of node ticks (as fraction of chord)
    pub tick_length: f64,
    /// Airfoil line color (RGB)
    pub line_color: (u8, u8, u8),
    /// Node tick color (RGB)
    pub tick_color: (u8, u8, u8),
    /// Background color (RGB)
    pub background: (u8, u8, u8),
    /// Title for the plot
    pub title: Option<String>,
}

impl Default for GeometryPlotConfig {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 400,
            show_nodes: false,
            tick_length: 0.015,
            line_color: (0, 100, 200),      // Blue
            tick_color: (200, 50, 50),      // Red
            background: (255, 255, 255),    // White
            title: None,
        }
    }
}

/// Error type for plotting operations
#[derive(thiserror::Error, Debug)]
pub enum PlotError {
    #[error("Drawing error: {0}")]
    Drawing(String),
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Plot raw airfoil geometry to an SVG file
pub fn plot_geometry_svg<P: AsRef<Path>>(
    geometry: &Geometry,
    output_path: P,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError> {
    let root = SVGBackend::new(&output_path, (config.width, config.height))
        .into_drawing_area();

    plot_geometry_impl(&root, geometry, None, config)
}

/// Plot paneled airfoil geometry to an SVG file with optional node ticks
///
/// Uses high-precision coordinate output (6 decimal places) for smooth curves.
pub fn plot_paneled_svg<P: AsRef<Path>>(
    airfoil: &PaneledAirfoil,
    output_path: P,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError> {
    write_precision_svg(airfoil, output_path, config)
}

/// Plot raw airfoil geometry to a PNG file
pub fn plot_geometry_png<P: AsRef<Path>>(
    geometry: &Geometry,
    output_path: P,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError> {
    let root = BitMapBackend::new(&output_path, (config.width, config.height))
        .into_drawing_area();

    plot_geometry_impl(&root, geometry, None, config)
}

/// Plot paneled airfoil geometry to a PNG file with optional node ticks
pub fn plot_paneled_png<P: AsRef<Path>>(
    airfoil: &PaneledAirfoil,
    output_path: P,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError> {
    let geometry = Geometry {
        reference: airfoil.reference,
        x_c: airfoil.x.clone(),
        y_c: airfoil.y.clone(),
    };

    let root = BitMapBackend::new(&output_path, (config.width, config.height))
        .into_drawing_area();

    plot_geometry_impl(&root, &geometry, Some(airfoil), config)
}

/// Internal implementation for geometry plotting
fn plot_geometry_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    geometry: &Geometry,
    paneled: Option<&PaneledAirfoil>,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError>
where
    DB::ErrorType: 'static,
{
    let bg_color = RGBColor(config.background.0, config.background.1, config.background.2);
    root.fill(&bg_color)
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Calculate data bounds with some padding
    let x_min = geometry.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = geometry.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = geometry.y_c.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = geometry.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Add padding
    let x_range = x_max - x_min;
    let padding = 0.1;

    let x_min_padded = x_min - padding * x_range;
    let x_max_padded = x_max + padding * x_range;

    // For equal aspect ratio, compute y range based on x range and viewport
    let aspect = config.width as f64 / config.height as f64;
    let y_center = (y_min + y_max) / 2.0;
    let y_half_range = (x_max_padded - x_min_padded) / aspect / 2.0;

    let y_min_padded = y_center - y_half_range;
    let y_max_padded = y_center + y_half_range;

    // Build chart
    let title = config.title.clone().unwrap_or_else(|| "Airfoil Geometry".to_string());

    let mut chart = ChartBuilder::on(root)
        .caption(&title, ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min_padded..x_max_padded, y_min_padded..y_max_padded)
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    chart
        .configure_mesh()
        .x_desc("x/c")
        .y_desc("y/c")
        .draw()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw airfoil outline
    let line_color = RGBColor(config.line_color.0, config.line_color.1, config.line_color.2);
    let points: Vec<(f64, f64)> = geometry
        .x_c
        .iter()
        .zip(geometry.y_c.iter())
        .map(|(&x, &y)| (x, y))
        .collect();

    chart
        .draw_series(LineSeries::new(points, line_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw node ticks if requested and we have paneled airfoil data
    if config.show_nodes {
        if let Some(airfoil) = paneled {
            draw_node_ticks(&mut chart, airfoil, config)?;
        }
    }

    root.present()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    Ok(())
}

/// Draw perpendicular tick marks at each panel node
fn draw_node_ticks<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    airfoil: &PaneledAirfoil,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError>
where
    DB::ErrorType: 'static,
{
    let tick_color = RGBColor(config.tick_color.0, config.tick_color.1, config.tick_color.2);
    let tick_len = config.tick_length;

    // Draw a small line segment at each node, perpendicular to the surface
    // Normal vectors point outward, so we draw inward (negative normal direction)
    for i in 0..airfoil.n {
        let x = airfoil.x[i];
        let y = airfoil.y[i];
        let nx = airfoil.nx[i];
        let ny = airfoil.ny[i];

        // Draw tick inward from surface (opposite to normal direction)
        let x_end = x - nx * tick_len;
        let y_end = y - ny * tick_len;

        chart
            .draw_series(LineSeries::new(
                vec![(x, y), (x_end, y_end)],
                tick_color.stroke_width(1),
            ))
            .map_err(|e| PlotError::Drawing(e.to_string()))?;
    }

    Ok(())
}

/// Write high-precision SVG with full coordinate precision
///
/// Unlike the plotters SVGBackend which rounds to integer pixels,
/// this outputs coordinates with 6 decimal places for smooth curves.
pub fn write_precision_svg<P: AsRef<Path>>(
    airfoil: &PaneledAirfoil,
    output_path: P,
    config: &GeometryPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Calculate bounds with padding
    let x_min = airfoil.x.iter().cloned().fold(f64::INFINITY, f64::min);
    let x_max = airfoil.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let y_min = airfoil.y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = airfoil.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let x_range = x_max - x_min;
    let padding = 0.1;
    let x_min_padded = x_min - padding * x_range;
    let x_max_padded = x_max + padding * x_range;

    // Equal aspect ratio
    let aspect = config.width as f64 / config.height as f64;
    let y_center = (y_min + y_max) / 2.0;
    let y_half_range = (x_max_padded - x_min_padded) / aspect / 2.0;
    let y_min_padded = y_center - y_half_range;
    let y_max_padded = y_center + y_half_range;

    let data_width = x_max_padded - x_min_padded;
    let data_height = y_max_padded - y_min_padded;

    // Coordinate transform: data coords to SVG coords
    let margin = 50.0;
    let plot_width = config.width as f64 - 2.0 * margin;
    let plot_height = config.height as f64 - 2.0 * margin;

    let to_svg_x = |x: f64| -> f64 { margin + (x - x_min_padded) / data_width * plot_width };
    let to_svg_y = |y: f64| -> f64 { margin + (y_max_padded - y) / data_height * plot_height };

    // SVG header
    writeln!(
        file,
        r#"<svg width="{}" height="{}" viewBox="0 0 {} {}" xmlns="http://www.w3.org/2000/svg">"#,
        config.width, config.height, config.width, config.height
    )?;

    // Background
    writeln!(
        file,
        r#"<rect width="100%" height="100%" fill="rgb({},{},{})"/>"#,
        config.background.0, config.background.1, config.background.2
    )?;

    // Title
    let title = config.title.clone().unwrap_or_else(|| "Airfoil Geometry".to_string());
    writeln!(
        file,
        r#"<text x="{}" y="25" text-anchor="middle" font-family="sans-serif" font-size="16">{}</text>"#,
        config.width as f64 / 2.0,
        title
    )?;

    // Axis labels
    writeln!(
        file,
        r#"<text x="{}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        config.width as f64 / 2.0,
        config.height as f64 - 10.0
    )?;
    writeln!(
        file,
        r#"<text x="15" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 15 {})">y/c</text>"#,
        config.height as f64 / 2.0,
        config.height as f64 / 2.0
    )?;

    // Grid lines (light gray)
    writeln!(file, r##"<g stroke="#CCCCCC" stroke-width="0.5">"##)?;
    for i in 0..=10 {
        let x = x_min_padded + (i as f64 / 10.0) * data_width;
        let sx = to_svg_x(x);
        writeln!(file, r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}"/>"#,
            sx, margin, sx, config.height as f64 - margin)?;
    }
    for i in 0..=5 {
        let y = y_min_padded + (i as f64 / 5.0) * data_height;
        let sy = to_svg_y(y);
        writeln!(file, r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}"/>"#,
            margin, sy, config.width as f64 - margin, sy)?;
    }
    writeln!(file, "</g>")?;

    // Plot border
    writeln!(
        file,
        r##"<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin, margin, plot_width, plot_height
    )?;

    // Airfoil polyline with high precision coordinates
    write!(
        file,
        r#"<polyline fill="none" stroke="rgb({},{},{})" stroke-width="2" points=""#,
        config.line_color.0, config.line_color.1, config.line_color.2
    )?;
    for i in 0..airfoil.n {
        let sx = to_svg_x(airfoil.x[i]);
        let sy = to_svg_y(airfoil.y[i]);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Node ticks if requested
    if config.show_nodes {
        writeln!(
            file,
            r#"<g stroke="rgb({},{},{})" stroke-width="1">"#,
            config.tick_color.0, config.tick_color.1, config.tick_color.2
        )?;
        for i in 0..airfoil.n {
            let x = airfoil.x[i];
            let y = airfoil.y[i];
            let x_end = x - airfoil.nx[i] * config.tick_length;
            let y_end = y - airfoil.ny[i] * config.tick_length;

            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}"/>"#,
                to_svg_x(x), to_svg_y(y), to_svg_x(x_end), to_svg_y(y_end)
            )?;
        }
        writeln!(file, "</g>")?;
    }

    writeln!(file, "</svg>")?;

    Ok(())
}
