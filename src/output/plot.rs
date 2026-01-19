//! Plotting utilities for airfoil visualization
//!
//! Uses the plotters library for PNG output and custom high-precision SVG output.

use plotters::prelude::*;
use std::io::Write;
use std::path::Path;

use crate::geometry::{Geometry, PaneledAirfoil};
use crate::output::InviscidAnalysisOutput;

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

/// Configuration for analysis distribution plots (Cp and Ue)
#[derive(Debug, Clone)]
pub struct AnalysisPlotConfig {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// Upper surface color (RGB) - default blue
    pub upper_color: (u8, u8, u8),
    /// Lower surface color (RGB) - default red
    pub lower_color: (u8, u8, u8),
    /// Airfoil outline color (RGB) - default gray
    pub airfoil_color: (u8, u8, u8),
    /// Background color (RGB)
    pub background: (u8, u8, u8),
    /// Title for the plot
    pub title: Option<String>,
}

impl Default for AnalysisPlotConfig {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 800,
            upper_color: (0, 100, 200),     // Blue
            lower_color: (200, 50, 50),     // Red
            airfoil_color: (100, 100, 100), // Gray
            background: (255, 255, 255),    // White
            title: None,
        }
    }
}

/// Plot Cp and Ue distributions from inviscid analysis results
///
/// Creates a two-panel plot with:
/// - Top panel: Cp distribution (negative upward, aerodynamic convention)
/// - Bottom panel: Ue (edge velocity) distribution
///
/// Both panels show upper surface in blue and lower surface in red,
/// with the airfoil outline superimposed at ~1/3 of the y-axis height.
///
/// Uses high-precision SVG output (6 decimal places) for smooth curves.
pub fn plot_analysis_svg<P: AsRef<Path>>(
    analysis: &InviscidAnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    write_analysis_precision_svg(analysis, output_path, config)
}

/// Plot Cp and Ue distributions to PNG
pub fn plot_analysis_png<P: AsRef<Path>>(
    analysis: &InviscidAnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    let root = BitMapBackend::new(&output_path, (config.width, config.height))
        .into_drawing_area();

    plot_analysis_impl(&root, analysis, config)
}

/// Internal implementation for analysis distribution plotting (for PNG)
fn plot_analysis_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    analysis: &InviscidAnalysisOutput,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError>
where
    DB::ErrorType: 'static,
{
    let bg_color = RGBColor(config.background.0, config.background.1, config.background.2);
    root.fill(&bg_color)
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Get surface data
    let (x_upper, cp_upper, vel_upper) = analysis.upper_surface();
    let (x_lower, cp_lower, vel_lower) = analysis.lower_surface();

    // Get airfoil outline
    let airfoil_x = &analysis.stations.x;
    let airfoil_y = &analysis.stations.y;

    // Calculate data ranges for Cp (note: inverted y-axis for negative convention)
    let cp_min = cp_upper.iter().chain(cp_lower.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let cp_max = cp_upper.iter().chain(cp_lower.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let cp_range = cp_max - cp_min;
    let cp_padding = 0.1 * cp_range;

    // Calculate data ranges for Ue
    let vel_min = vel_upper.iter().chain(vel_lower.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let vel_max = vel_upper.iter().chain(vel_lower.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let vel_range = vel_max - vel_min;
    let vel_padding = 0.1 * vel_range;

    // Airfoil y range
    let y_min = airfoil_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = airfoil_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Split into upper and lower regions with margin for title
    let title = config.title.clone().unwrap_or_else(|| {
        format!("{} at α = {:.1}°", analysis.airfoil, analysis.alpha_deg)
    });

    let (upper, lower) = root.split_vertically(config.height / 2);

    // Colors
    let upper_color = RGBColor(config.upper_color.0, config.upper_color.1, config.upper_color.2);
    let lower_color = RGBColor(config.lower_color.0, config.lower_color.1, config.lower_color.2);
    let airfoil_color = RGBColor(config.airfoil_color.0, config.airfoil_color.1, config.airfoil_color.2);

    // ===== Cp Plot (upper panel) =====
    // Y-axis is INVERTED: more negative values at top (suction), positive at bottom (pressure)
    let cp_y_min = cp_min - cp_padding;  // This will be at the TOP (most negative = most suction)
    let cp_y_max = cp_max + cp_padding;  // This will be at the BOTTOM (most positive = pressure)

    let mut cp_chart = ChartBuilder::on(&upper)
        .caption(&title, ("sans-serif", 18))
        .margin(10)
        .x_label_area_size(35)
        .y_label_area_size(50)
        .build_cartesian_2d(
            -0.05..1.1,
            cp_y_max..cp_y_min,  // INVERTED: positive at bottom, negative at top
        )
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    cp_chart
        .configure_mesh()
        .x_desc("x/c")
        .y_desc("Cp")
        .x_labels(10)
        .y_labels(8)
        .light_line_style(TRANSPARENT)
        .draw()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw airfoil outline scaled to ~1/3 of plot height at the bottom of Cp range
    // Map airfoil y to Cp coordinates: scale y to fit in bottom 1/3 of Cp range
    let cp_display_range = cp_y_max - cp_y_min;
    let airfoil_scale = (cp_display_range / 3.0) / (y_max - y_min);
    let airfoil_offset = cp_y_max - 0.05 * cp_display_range; // Position near the bottom (high Cp values)

    let airfoil_cp_points: Vec<(f64, f64)> = airfoil_x.iter()
        .zip(airfoil_y.iter())
        .map(|(&x, &y)| (x, airfoil_offset - (y - y_min) * airfoil_scale))
        .collect();

    cp_chart
        .draw_series(LineSeries::new(airfoil_cp_points, airfoil_color.stroke_width(1)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw Cp distributions
    let upper_cp_points: Vec<(f64, f64)> = x_upper.iter()
        .zip(cp_upper.iter())
        .map(|(&x, &cp)| (x, cp))
        .collect();

    let lower_cp_points: Vec<(f64, f64)> = x_lower.iter()
        .zip(cp_lower.iter())
        .map(|(&x, &cp)| (x, cp))
        .collect();

    cp_chart
        .draw_series(LineSeries::new(upper_cp_points, upper_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?
        .label("Upper surface")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], upper_color.stroke_width(2)));

    cp_chart
        .draw_series(LineSeries::new(lower_cp_points, lower_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?
        .label("Lower surface")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], lower_color.stroke_width(2)));

    cp_chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // ===== Ue Plot (lower panel) =====
    let vel_y_min = vel_min - vel_padding;
    let vel_y_max = vel_max + vel_padding;

    let mut vel_chart = ChartBuilder::on(&lower)
        .margin(10)
        .x_label_area_size(35)
        .y_label_area_size(50)
        .build_cartesian_2d(-0.05..1.1, vel_y_min..vel_y_max)
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    vel_chart
        .configure_mesh()
        .x_desc("x/c")
        .y_desc("Ue/U∞")
        .x_labels(10)
        .y_labels(8)
        .light_line_style(TRANSPARENT)
        .draw()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw airfoil outline scaled to ~1/3 of plot height at the bottom of velocity range
    let vel_display_range = vel_y_max - vel_y_min;
    let airfoil_vel_scale = (vel_display_range / 3.0) / (y_max - y_min);
    let airfoil_vel_offset = vel_y_min + 0.05 * vel_display_range;

    let airfoil_vel_points: Vec<(f64, f64)> = airfoil_x.iter()
        .zip(airfoil_y.iter())
        .map(|(&x, &y)| (x, airfoil_vel_offset + (y - y_min) * airfoil_vel_scale))
        .collect();

    vel_chart
        .draw_series(LineSeries::new(airfoil_vel_points, airfoil_color.stroke_width(1)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw Ue distributions
    let upper_vel_points: Vec<(f64, f64)> = x_upper.iter()
        .zip(vel_upper.iter())
        .map(|(&x, &v)| (x, v))
        .collect();

    let lower_vel_points: Vec<(f64, f64)> = x_lower.iter()
        .zip(vel_lower.iter())
        .map(|(&x, &v)| (x, v))
        .collect();

    vel_chart
        .draw_series(LineSeries::new(upper_vel_points, upper_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    vel_chart
        .draw_series(LineSeries::new(lower_vel_points, lower_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    root.present()
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    Ok(())
}

/// Write high-precision SVG for analysis plots
///
/// Unlike the plotters SVGBackend which rounds to integer pixels,
/// this outputs coordinates with 6 decimal places for smooth curves.
fn write_analysis_precision_svg<P: AsRef<Path>>(
    analysis: &InviscidAnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Get surface data
    let (x_upper, cp_upper, vel_upper) = analysis.upper_surface();
    let (x_lower, cp_lower, vel_lower) = analysis.lower_surface();

    // Get airfoil outline
    let airfoil_x = &analysis.stations.x;
    let airfoil_y = &analysis.stations.y;

    // Calculate data ranges for Cp
    let cp_min = cp_upper.iter().chain(cp_lower.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let cp_max = cp_upper.iter().chain(cp_lower.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let cp_range = cp_max - cp_min;
    let cp_padding = 0.1 * cp_range;
    let cp_y_min = cp_min - cp_padding;
    let cp_y_max = cp_max + cp_padding;

    // Calculate data ranges for Ue
    let vel_min = vel_upper.iter().chain(vel_lower.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let vel_max = vel_upper.iter().chain(vel_lower.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let vel_range = vel_max - vel_min;
    let vel_padding = 0.1 * vel_range;
    let vel_y_min = vel_min - vel_padding;
    let vel_y_max = vel_max + vel_padding;

    // Airfoil y range
    let af_y_min = airfoil_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let af_y_max = airfoil_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Layout dimensions
    let width = config.width as f64;
    let height = config.height as f64;
    let margin_left = 60.0;
    let margin_right = 20.0;
    let margin_top = 40.0;
    let margin_bottom = 40.0;
    let gap = 20.0; // Gap between plots

    let plot_width = width - margin_left - margin_right;
    let plot_height = (height - margin_top - margin_bottom - gap) / 2.0;

    // X data range (same for both plots)
    let x_data_min = -0.05;
    let x_data_max = 1.1;
    let x_data_range = x_data_max - x_data_min;

    // Coordinate transforms for Cp plot (top)
    let cp_plot_top = margin_top;
    let cp_plot_bottom = margin_top + plot_height;
    let cp_data_range = cp_y_max - cp_y_min;

    let to_svg_x = |x: f64| -> f64 {
        margin_left + (x - x_data_min) / x_data_range * plot_width
    };
    // Cp: inverted y-axis (negative at top, positive at bottom)
    let to_svg_cp_y = |cp: f64| -> f64 {
        cp_plot_top + (cp - cp_y_min) / cp_data_range * plot_height
    };

    // Coordinate transforms for Ue plot (bottom)
    let vel_plot_top = cp_plot_bottom + gap;
    let vel_plot_bottom = vel_plot_top + plot_height;
    let vel_data_range = vel_y_max - vel_y_min;

    let to_svg_vel_y = |vel: f64| -> f64 {
        vel_plot_bottom - (vel - vel_y_min) / vel_data_range * plot_height
    };

    // Colors
    let upper_color = format!("rgb({},{},{})", config.upper_color.0, config.upper_color.1, config.upper_color.2);
    let lower_color = format!("rgb({},{},{})", config.lower_color.0, config.lower_color.1, config.lower_color.2);
    let airfoil_color = format!("rgb({},{},{})", config.airfoil_color.0, config.airfoil_color.1, config.airfoil_color.2);
    let grid_color = "#DDDDDD"; // Light grey for major gridlines

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
    let title = config.title.clone().unwrap_or_else(|| {
        format!("{} at α = {:.1}°", analysis.airfoil, analysis.alpha_deg)
    });
    writeln!(
        file,
        r#"<text x="{:.1}" y="25" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold">{}</text>"#,
        width / 2.0,
        title
    )?;

    // ===== Cp Plot (top) =====

    // Plot border
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin_left, cp_plot_top, plot_width, plot_height
    )?;

    // Major gridlines for Cp (light grey, no minor gridlines)
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    // Vertical gridlines at x = 0.0, 0.1, 0.2, ... 1.0
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        if sx > margin_left && sx < margin_left + plot_width {
            writeln!(file, r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, cp_plot_top, sx, cp_plot_bottom)?;
        }
    }
    // Horizontal gridlines (Cp)
    let cp_step = nice_step(cp_data_range, 6);
    let cp_start = (cp_y_min / cp_step).floor() * cp_step;
    let mut cp_val = cp_start;
    while cp_val <= cp_y_max {
        let sy = to_svg_cp_y(cp_val);
        if sy > cp_plot_top && sy < cp_plot_bottom {
            writeln!(file, r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                margin_left, sy, margin_left + plot_width, sy)?;
        }
        cp_val += cp_step;
    }
    writeln!(file, "</g>")?;

    // Cp axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        margin_left + plot_width / 2.0,
        cp_plot_bottom + 30.0
    )?;
    writeln!(
        file,
        r#"<text x="15" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 15 {:.1})">Cp</text>"#,
        cp_plot_top + plot_height / 2.0,
        cp_plot_top + plot_height / 2.0
    )?;

    // Cp tick labels
    writeln!(file, r#"<g font-family="sans-serif" font-size="10" text-anchor="end">"#)?;
    cp_val = cp_start;
    while cp_val <= cp_y_max {
        let sy = to_svg_cp_y(cp_val);
        if sy > cp_plot_top + 5.0 && sy < cp_plot_bottom - 5.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
                margin_left - 5.0, sy + 4.0, cp_val)?;
        }
        cp_val += cp_step;
    }
    writeln!(file, "</g>")?;

    // X tick labels for Cp plot
    writeln!(file, r#"<g font-family="sans-serif" font-size="10" text-anchor="middle">"#)?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
            sx, cp_plot_bottom + 15.0, x)?;
    }
    writeln!(file, "</g>")?;

    // Airfoil outline on Cp plot (scaled to ~1/3 height, centered at Cp=0)
    let cp_display_range_val = cp_y_max - cp_y_min;
    let airfoil_scale = (cp_display_range_val / 3.0) / (af_y_max - af_y_min);

    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="1" points=""#, airfoil_color)?;
    for (i, (&x, &y)) in airfoil_x.iter().zip(airfoil_y.iter()).enumerate() {
        // Map airfoil y directly to Cp coordinates (y=0 on airfoil -> Cp=0)
        let cp_y = -y * airfoil_scale;
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp_y);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Upper surface Cp (blue)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, upper_color)?;
    for (i, (&x, &cp)) in x_upper.iter().zip(cp_upper.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Lower surface Cp (red)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, lower_color)?;
    for (i, (&x, &cp)) in x_lower.iter().zip(cp_lower.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Legend for Cp plot
    let legend_x = margin_left + plot_width - 120.0;
    let legend_y = cp_plot_top + 15.0;
    writeln!(file, r#"<rect x="{:.1}" y="{:.1}" width="110" height="45" fill="white" fill-opacity="0.8" stroke="black" stroke-width="0.5"/>"#,
        legend_x, legend_y)?;
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0, legend_y + 15.0, legend_x + 30.0, legend_y + 15.0, upper_color)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">Upper surface</text>"#,
        legend_x + 35.0, legend_y + 19.0)?;
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0, legend_y + 32.0, legend_x + 30.0, legend_y + 32.0, lower_color)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">Lower surface</text>"#,
        legend_x + 35.0, legend_y + 36.0)?;

    // ===== Ue Plot (bottom) =====

    // Plot border
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin_left, vel_plot_top, plot_width, plot_height
    )?;

    // Major gridlines for Ue (light grey)
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    // Vertical gridlines at x = 0.0, 0.1, 0.2, ... 1.0
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        if sx > margin_left && sx < margin_left + plot_width {
            writeln!(file, r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, vel_plot_top, sx, vel_plot_bottom)?;
        }
    }
    // Horizontal gridlines (Ue)
    let vel_step = nice_step(vel_data_range, 6);
    let vel_start = (vel_y_min / vel_step).floor() * vel_step;
    let mut vel_val = vel_start;
    while vel_val <= vel_y_max {
        let sy = to_svg_vel_y(vel_val);
        if sy > vel_plot_top && sy < vel_plot_bottom {
            writeln!(file, r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                margin_left, sy, margin_left + plot_width, sy)?;
        }
        vel_val += vel_step;
    }
    writeln!(file, "</g>")?;

    // Ue axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        margin_left + plot_width / 2.0,
        vel_plot_bottom + 30.0
    )?;
    writeln!(
        file,
        r#"<text x="15" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 15 {:.1})">Ue/U∞</text>"#,
        vel_plot_top + plot_height / 2.0,
        vel_plot_top + plot_height / 2.0
    )?;

    // Ue tick labels
    writeln!(file, r#"<g font-family="sans-serif" font-size="10" text-anchor="end">"#)?;
    vel_val = vel_start;
    while vel_val <= vel_y_max {
        let sy = to_svg_vel_y(vel_val);
        if sy > vel_plot_top + 5.0 && sy < vel_plot_bottom - 5.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.2}</text>"#,
                margin_left - 5.0, sy + 4.0, vel_val)?;
        }
        vel_val += vel_step;
    }
    writeln!(file, "</g>")?;

    // X tick labels for Ue plot
    writeln!(file, r#"<g font-family="sans-serif" font-size="10" text-anchor="middle">"#)?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
            sx, vel_plot_bottom + 15.0, x)?;
    }
    writeln!(file, "</g>")?;

    // Airfoil outline on Ue plot (scaled to ~1/3 height, centered at Ue=0)
    let vel_display_range_val = vel_y_max - vel_y_min;
    let airfoil_vel_scale = (vel_display_range_val / 3.0) / (af_y_max - af_y_min);

    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="1" points=""#, airfoil_color)?;
    for (i, (&x, &y)) in airfoil_x.iter().zip(airfoil_y.iter()).enumerate() {
        // Map airfoil y directly to Ue coordinates (y=0 on airfoil -> Ue=0)
        let vel_y = y * airfoil_vel_scale;
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel_y);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Upper surface Ue (blue)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, upper_color)?;
    for (i, (&x, &vel)) in x_upper.iter().zip(vel_upper.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Lower surface Ue (red)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, lower_color)?;
    for (i, (&x, &vel)) in x_lower.iter().zip(vel_lower.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    writeln!(file, "</svg>")?;

    Ok(())
}

/// Calculate a nice step size for axis ticks
fn nice_step(range: f64, target_ticks: usize) -> f64 {
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

// ============================================================================
// Polar Plotting
// ============================================================================

/// Configuration for polar plots
#[derive(Debug, Clone)]
pub struct PolarPlotConfig {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// First data series color (RGB) - default blue
    pub color_1: (u8, u8, u8),
    /// Second data series color (RGB) - default red (for comparison)
    pub color_2: (u8, u8, u8),
    /// Background color (RGB)
    pub background: (u8, u8, u8),
    /// Title for the plot
    pub title: Option<String>,
    /// Label for first data series
    pub label_1: String,
    /// Label for second data series (for comparison)
    pub label_2: String,
}

impl Default for PolarPlotConfig {
    fn default() -> Self {
        Self {
            width: 1400,
            height: 600,
            color_1: (0, 100, 200),      // Blue
            color_2: (200, 50, 50),      // Red
            background: (255, 255, 255), // White
            title: None,
            label_1: "YFoil".to_string(),
            label_2: "XFOIL".to_string(),
        }
    }
}

/// Data for a single polar curve
#[derive(Debug, Clone)]
pub struct PolarData {
    /// Angle of attack values (degrees)
    pub alpha: Vec<f64>,
    /// Lift coefficient values
    pub cl: Vec<f64>,
    /// Drag coefficient values
    pub cd: Vec<f64>,
}

impl PolarData {
    /// Create from PolarOutput
    pub fn from_polar_output(polar: &crate::output::PolarOutput) -> Self {
        let alpha: Vec<f64> = polar.points.iter().map(|p| p.alpha_deg).collect();
        let cl: Vec<f64> = polar.points.iter().map(|p| p.cl).collect();
        let cd: Vec<f64> = polar.points.iter().map(|p| p.cd).collect();
        Self { alpha, cl, cd }
    }

    /// Create from raw data vectors
    pub fn new(alpha: Vec<f64>, cl: Vec<f64>, cd: Vec<f64>) -> Self {
        Self { alpha, cl, cd }
    }
}

/// Plot a polar comparison to SVG
///
/// Creates a two-panel plot comparing two polars:
/// - Left panel: CL vs Alpha
/// - Right panel: CD vs Alpha
pub fn plot_polar_comparison_svg<P: AsRef<Path>>(
    polar_1: &PolarData,
    polar_2: &PolarData,
    output_path: P,
    config: &PolarPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Calculate data ranges
    let alpha_min = polar_1.alpha.iter().chain(polar_2.alpha.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let alpha_max = polar_1.alpha.iter().chain(polar_2.alpha.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let alpha_range = alpha_max - alpha_min;
    let alpha_padding = 0.05 * alpha_range;

    let cl_min = polar_1.cl.iter().chain(polar_2.cl.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let cl_max = polar_1.cl.iter().chain(polar_2.cl.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let cl_range = cl_max - cl_min;
    let cl_padding = 0.1 * cl_range;

    let cd_min = polar_1.cd.iter().chain(polar_2.cd.iter())
        .cloned().fold(f64::INFINITY, f64::min);
    let cd_max = polar_1.cd.iter().chain(polar_2.cd.iter())
        .cloned().fold(f64::NEG_INFINITY, f64::max);
    let cd_range = cd_max - cd_min;
    let cd_padding = 0.1 * cd_range;

    // Padded ranges
    let alpha_min_p = alpha_min - alpha_padding;
    let alpha_max_p = alpha_max + alpha_padding;
    let cl_min_p = cl_min - cl_padding;
    let cl_max_p = cl_max + cl_padding;
    let cd_min_p = cd_min - cd_padding;
    let cd_max_p = cd_max + cd_padding;

    // Layout
    let width = config.width as f64;
    let height = config.height as f64;
    let margin_left = 70.0;
    let margin_right = 20.0;
    let margin_top = 50.0;
    let margin_bottom = 50.0;
    let gap = 60.0;

    let plot_width = (width - margin_left - margin_right - gap) / 2.0;
    let plot_height = height - margin_top - margin_bottom;

    // Plot regions
    let cl_plot_left = margin_left;
    let cd_plot_left = margin_left + plot_width + gap;

    // Colors
    let color_1 = format!("rgb({},{},{})", config.color_1.0, config.color_1.1, config.color_1.2);
    let color_2 = format!("rgb({},{},{})", config.color_2.0, config.color_2.1, config.color_2.2);
    let grid_color = "#DDDDDD";

    // Coordinate transforms
    let to_svg_alpha_cl = |alpha: f64| -> f64 {
        cl_plot_left + (alpha - alpha_min_p) / (alpha_max_p - alpha_min_p) * plot_width
    };
    let to_svg_cl = |cl: f64| -> f64 {
        margin_top + plot_height - (cl - cl_min_p) / (cl_max_p - cl_min_p) * plot_height
    };
    let to_svg_alpha_cd = |alpha: f64| -> f64 {
        cd_plot_left + (alpha - alpha_min_p) / (alpha_max_p - alpha_min_p) * plot_width
    };
    let to_svg_cd = |cd: f64| -> f64 {
        margin_top + plot_height - (cd - cd_min_p) / (cd_max_p - cd_min_p) * plot_height
    };

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
    let title = config.title.clone().unwrap_or_else(|| "Polar Comparison".to_string());
    writeln!(
        file,
        r#"<text x="{:.1}" y="30" text-anchor="middle" font-family="sans-serif" font-size="18" font-weight="bold">{}</text>"#,
        width / 2.0,
        title
    )?;

    // ===== CL vs Alpha Plot (left) =====

    // Plot border
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        cl_plot_left, margin_top, plot_width, plot_height
    )?;

    // Gridlines
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    let alpha_step = nice_step(alpha_max_p - alpha_min_p, 8);
    let alpha_start = (alpha_min_p / alpha_step).ceil() * alpha_step;
    let mut alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cl(alpha_val);
        if sx > cl_plot_left && sx < cl_plot_left + plot_width {
            writeln!(file, r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, margin_top, sx, margin_top + plot_height)?;
        }
        alpha_val += alpha_step;
    }
    let cl_step = nice_step(cl_max_p - cl_min_p, 8);
    let cl_start = (cl_min_p / cl_step).ceil() * cl_step;
    let mut cl_val = cl_start;
    while cl_val <= cl_max_p {
        let sy = to_svg_cl(cl_val);
        if sy > margin_top && sy < margin_top + plot_height {
            writeln!(file, r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                cl_plot_left, sy, cl_plot_left + plot_width, sy)?;
        }
        cl_val += cl_step;
    }
    writeln!(file, "</g>")?;

    // Axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="14">Alpha (deg)</text>"#,
        cl_plot_left + plot_width / 2.0,
        margin_top + plot_height + 40.0
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 {:.1} {:.1})">CL</text>"#,
        cl_plot_left - 50.0,
        margin_top + plot_height / 2.0,
        cl_plot_left - 50.0,
        margin_top + plot_height / 2.0
    )?;

    // Tick labels for CL plot
    writeln!(file, r#"<g font-family="sans-serif" font-size="11" text-anchor="end">"#)?;
    cl_val = cl_start;
    while cl_val <= cl_max_p {
        let sy = to_svg_cl(cl_val);
        if sy > margin_top + 5.0 && sy < margin_top + plot_height - 5.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.2}</text>"#,
                cl_plot_left - 5.0, sy + 4.0, cl_val)?;
        }
        cl_val += cl_step;
    }
    writeln!(file, "</g>")?;

    writeln!(file, r#"<g font-family="sans-serif" font-size="11" text-anchor="middle">"#)?;
    alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cl(alpha_val);
        if sx > cl_plot_left + 10.0 && sx < cl_plot_left + plot_width - 10.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.0}</text>"#,
                sx, margin_top + plot_height + 18.0, alpha_val)?;
        }
        alpha_val += alpha_step;
    }
    writeln!(file, "</g>")?;

    // Data series 1 (CL)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_1)?;
    for (i, (&alpha, &cl)) in polar_1.alpha.iter().zip(polar_1.cl.iter()).enumerate() {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 1
    for (&alpha, &cl) in polar_1.alpha.iter().zip(polar_1.cl.iter()) {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        writeln!(file, r#"<circle cx="{:.6}" cy="{:.6}" r="3" fill="{}" />"#, sx, sy, color_1)?;
    }

    // Data series 2 (CL)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_2)?;
    for (i, (&alpha, &cl)) in polar_2.alpha.iter().zip(polar_2.cl.iter()).enumerate() {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 2
    for (&alpha, &cl) in polar_2.alpha.iter().zip(polar_2.cl.iter()) {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        writeln!(file, r#"<rect x="{:.6}" y="{:.6}" width="6" height="6" fill="{}" transform="rotate(45 {:.6} {:.6})"/>"#,
            sx - 3.0, sy - 3.0, color_2, sx, sy)?;
    }

    // ===== CD vs Alpha Plot (right) =====

    // Plot border
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        cd_plot_left, margin_top, plot_width, plot_height
    )?;

    // Gridlines
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cd(alpha_val);
        if sx > cd_plot_left && sx < cd_plot_left + plot_width {
            writeln!(file, r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, margin_top, sx, margin_top + plot_height)?;
        }
        alpha_val += alpha_step;
    }
    let cd_step = nice_step(cd_max_p - cd_min_p, 8);
    let cd_start = (cd_min_p / cd_step).ceil() * cd_step;
    let mut cd_val = cd_start;
    while cd_val <= cd_max_p {
        let sy = to_svg_cd(cd_val);
        if sy > margin_top && sy < margin_top + plot_height {
            writeln!(file, r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                cd_plot_left, sy, cd_plot_left + plot_width, sy)?;
        }
        cd_val += cd_step;
    }
    writeln!(file, "</g>")?;

    // Axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="14">Alpha (deg)</text>"#,
        cd_plot_left + plot_width / 2.0,
        margin_top + plot_height + 40.0
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 {:.1} {:.1})">CD</text>"#,
        cd_plot_left - 50.0,
        margin_top + plot_height / 2.0,
        cd_plot_left - 50.0,
        margin_top + plot_height / 2.0
    )?;

    // Tick labels for CD plot
    writeln!(file, r#"<g font-family="sans-serif" font-size="11" text-anchor="end">"#)?;
    cd_val = cd_start;
    while cd_val <= cd_max_p {
        let sy = to_svg_cd(cd_val);
        if sy > margin_top + 5.0 && sy < margin_top + plot_height - 5.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.4}</text>"#,
                cd_plot_left - 5.0, sy + 4.0, cd_val)?;
        }
        cd_val += cd_step;
    }
    writeln!(file, "</g>")?;

    writeln!(file, r#"<g font-family="sans-serif" font-size="11" text-anchor="middle">"#)?;
    alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cd(alpha_val);
        if sx > cd_plot_left + 10.0 && sx < cd_plot_left + plot_width - 10.0 {
            writeln!(file, r#"<text x="{:.1}" y="{:.1}">{:.0}</text>"#,
                sx, margin_top + plot_height + 18.0, alpha_val)?;
        }
        alpha_val += alpha_step;
    }
    writeln!(file, "</g>")?;

    // Data series 1 (CD)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_1)?;
    for (i, (&alpha, &cd)) in polar_1.alpha.iter().zip(polar_1.cd.iter()).enumerate() {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 1
    for (&alpha, &cd) in polar_1.alpha.iter().zip(polar_1.cd.iter()) {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        writeln!(file, r#"<circle cx="{:.6}" cy="{:.6}" r="3" fill="{}" />"#, sx, sy, color_1)?;
    }

    // Data series 2 (CD)
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_2)?;
    for (i, (&alpha, &cd)) in polar_2.alpha.iter().zip(polar_2.cd.iter()).enumerate() {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 2
    for (&alpha, &cd) in polar_2.alpha.iter().zip(polar_2.cd.iter()) {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        writeln!(file, r#"<rect x="{:.6}" y="{:.6}" width="6" height="6" fill="{}" transform="rotate(45 {:.6} {:.6})"/>"#,
            sx - 3.0, sy - 3.0, color_2, sx, sy)?;
    }

    // Legend (centered at bottom)
    let legend_x = width / 2.0 - 100.0;
    let legend_y = height - 25.0;

    // Series 1
    writeln!(file, r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{}" />"#,
        legend_x, legend_y, color_1)?;
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x - 15.0, legend_y, legend_x + 15.0, legend_y, color_1)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x + 25.0, legend_y + 4.0, config.label_1)?;

    // Series 2
    let legend_x2 = legend_x + 120.0;
    writeln!(file, r#"<rect x="{:.1}" y="{:.1}" width="8" height="8" fill="{}" transform="rotate(45 {:.1} {:.1})"/>"#,
        legend_x2 - 4.0, legend_y - 4.0, color_2, legend_x2, legend_y)?;
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x2 - 15.0, legend_y, legend_x2 + 15.0, legend_y, color_2)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x2 + 25.0, legend_y + 4.0, config.label_2)?;

    writeln!(file, "</svg>")?;

    Ok(())
}

/// Parse an XFOIL polar file and extract alpha, CL, CD data
pub fn parse_xfoil_polar_file<P: AsRef<Path>>(path: P) -> Result<PolarData, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let mut alphas = Vec::new();
    let mut cls = Vec::new();
    let mut cds = Vec::new();

    let mut in_data = false;
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("------") {
            in_data = true;
            continue;
        }
        if in_data && !line.is_empty() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                if let (Ok(alpha), Ok(cl), Ok(cd)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                ) {
                    alphas.push(alpha);
                    cls.push(cl);
                    cds.push(cd);
                }
            }
        }
    }

    Ok(PolarData::new(alphas, cls, cds))
}

/// Stitch two polars together in ascending alpha order
///
/// Used for combining positive and negative alpha sweeps into a single polar.
/// Removes duplicate alpha values (keeps first occurrence).
pub fn stitch_polars(polar_pos: &PolarData, polar_neg: &PolarData) -> PolarData {
    // Reverse negative polar to get ascending order
    let mut neg_alpha: Vec<f64> = polar_neg.alpha.iter().copied().rev().collect();
    let mut neg_cl: Vec<f64> = polar_neg.cl.iter().copied().rev().collect();
    let mut neg_cd: Vec<f64> = polar_neg.cd.iter().copied().rev().collect();

    // Concatenate
    neg_alpha.extend(polar_pos.alpha.iter().copied());
    neg_cl.extend(polar_pos.cl.iter().copied());
    neg_cd.extend(polar_pos.cd.iter().copied());

    // Sort by alpha and remove duplicates
    let mut combined: Vec<(f64, f64, f64)> = neg_alpha.into_iter()
        .zip(neg_cl.into_iter())
        .zip(neg_cd.into_iter())
        .map(|((a, c), d)| (a, c, d))
        .collect();

    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Remove duplicates (keep first occurrence of each alpha)
    let mut seen_alphas = std::collections::HashSet::new();
    let combined: Vec<(f64, f64, f64)> = combined.into_iter()
        .filter(|(a, _, _)| {
            let key = (*a * 1000.0).round() as i64; // Round to 0.001 precision
            seen_alphas.insert(key)
        })
        .collect();

    let alphas: Vec<f64> = combined.iter().map(|(a, _, _)| *a).collect();
    let cls: Vec<f64> = combined.iter().map(|(_, c, _)| *c).collect();
    let cds: Vec<f64> = combined.iter().map(|(_, _, d)| *d).collect();

    PolarData::new(alphas, cls, cds)
}

/// Parse a YFoil polar JSON file into PolarData
pub fn parse_yfoil_polar_file<P: AsRef<Path>>(path: P) -> Result<PolarData, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let polar: crate::output::PolarOutput = serde_json::from_str(&content)
        .map_err(|e| PlotError::Drawing(format!("Failed to parse YFoil polar JSON: {}", e)))?;
    Ok(PolarData::from_polar_output(&polar))
}

// ============================================================================
// Geometry Comparison Plotting
// ============================================================================

/// Configuration for geometry comparison plots
#[derive(Debug, Clone)]
pub struct GeometryComparisonPlotConfig {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// First series color (RGB) - XFOIL reference
    pub series1_color: (u8, u8, u8),
    /// Second series color (RGB) - YFoil generated
    pub series2_color: (u8, u8, u8),
    /// Marker size for data points
    pub marker_size: f64,
    /// Background color (RGB)
    pub background: (u8, u8, u8),
    /// Title for the plot
    pub title: Option<String>,
    /// Label for first series
    pub label_1: String,
    /// Label for second series
    pub label_2: String,
}

impl Default for GeometryComparisonPlotConfig {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 400,
            series1_color: (0, 100, 200),     // Blue for XFOIL
            series2_color: (200, 50, 50),     // Red for YFoil
            marker_size: 4.0,
            background: (255, 255, 255),      // White
            title: None,
            label_1: "XFOIL".to_string(),
            label_2: "YFoil".to_string(),
        }
    }
}

/// Plot geometry comparison to SVG
///
/// Creates an overlay plot comparing two sets of coordinates:
/// - Series 1 (XFOIL): blue line with circle markers
/// - Series 2 (YFoil): red line with cross markers
///
/// Uses high-precision SVG output for smooth curves.
pub fn plot_geometry_comparison_svg<P: AsRef<Path>>(
    coords1: &[(f64, f64)],
    coords2: &[(f64, f64)],
    output_path: P,
    config: &GeometryComparisonPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Calculate combined bounds
    let x_min = coords1.iter().chain(coords2.iter())
        .map(|(x, _)| *x)
        .fold(f64::INFINITY, f64::min);
    let x_max = coords1.iter().chain(coords2.iter())
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = coords1.iter().chain(coords2.iter())
        .map(|(_, y)| *y)
        .fold(f64::INFINITY, f64::min);
    let y_max = coords1.iter().chain(coords2.iter())
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);

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

    // Colors
    let color_1 = format!("rgb({},{},{})", config.series1_color.0, config.series1_color.1, config.series1_color.2);
    let color_2 = format!("rgb({},{},{})", config.series2_color.0, config.series2_color.1, config.series2_color.2);
    let grid_color = "#CCCCCC";

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
    let title = config.title.clone().unwrap_or_else(|| "Geometry Comparison".to_string());
    writeln!(
        file,
        r#"<text x="{}" y="25" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold">{}</text>"#,
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
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
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

    // Series 1 (XFOIL) - blue line with circle markers
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_1)?;
    for (i, (x, y)) in coords1.iter().enumerate() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Circle markers for series 1
    for (x, y) in coords1.iter() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        writeln!(file, r#"<circle cx="{:.6}" cy="{:.6}" r="{}" fill="none" stroke="{}" stroke-width="1.5"/>"#,
            sx, sy, config.marker_size, color_1)?;
    }

    // Series 2 (YFoil) - red line with cross markers
    write!(file, r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#, color_2)?;
    for (i, (x, y)) in coords2.iter().enumerate() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        if i > 0 { write!(file, " ")?; }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Cross markers for series 2
    let m = config.marker_size;
    for (x, y) in coords2.iter() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        writeln!(file, r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
            sx - m, sy - m, sx + m, sy + m, color_2)?;
        writeln!(file, r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
            sx - m, sy + m, sx + m, sy - m, color_2)?;
    }

    // Legend
    let legend_x = config.width as f64 - 150.0;
    let legend_y = margin + 15.0;
    writeln!(file, r#"<rect x="{:.1}" y="{:.1}" width="130" height="50" fill="white" fill-opacity="0.9" stroke="black" stroke-width="0.5"/>"#,
        legend_x, legend_y)?;

    // Series 1 legend
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0, legend_y + 15.0, legend_x + 30.0, legend_y + 15.0, color_1)?;
    writeln!(file, r#"<circle cx="{:.1}" cy="{:.1}" r="{}" fill="none" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0, legend_y + 15.0, config.marker_size, color_1)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">{}</text>"#,
        legend_x + 40.0, legend_y + 19.0, config.label_1)?;

    // Series 2 legend
    writeln!(file, r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0, legend_y + 35.0, legend_x + 30.0, legend_y + 35.0, color_2)?;
    writeln!(file, r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0 - m, legend_y + 35.0 - m, legend_x + 20.0 + m, legend_y + 35.0 + m, color_2)?;
    writeln!(file, r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0 - m, legend_y + 35.0 + m, legend_x + 20.0 + m, legend_y + 35.0 - m, color_2)?;
    writeln!(file, r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">{}</text>"#,
        legend_x + 40.0, legend_y + 39.0, config.label_2)?;

    writeln!(file, "</svg>")?;

    Ok(())
}

/// Parse an XFOIL plain coordinate file (.dat format)
///
/// Returns coordinates as Vec<(x, y)>.
pub fn parse_xfoil_dat_file<P: AsRef<Path>>(path: P) -> Result<Vec<(f64, f64)>, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let mut coords = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                coords.push((x, y));
            }
        }
    }

    Ok(coords)
}
