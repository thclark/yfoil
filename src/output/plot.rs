//! Plotting utilities for airfoil visualization
//!
//! Uses the plotters library for PNG output and custom high-precision SVG output. The airfoil/foil
//! plot lives in `foil_plot.rs` on the `svg.rs` canvas.

use plotters::element::ComposedElement;
use plotters::prelude::*;
use std::io::Write;
use std::path::Path;

use crate::output::svg::nice_step;
use crate::output::AnalysisOutput;

/// Error type for plotting operations
#[derive(thiserror::Error, Debug)]
pub enum PlotError {
    #[error("Drawing error: {0}")]
    Drawing(String),
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Config(String),
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
    analysis: &AnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    write_analysis_precision_svg(analysis, output_path, config)
}

/// Plot Cp and Ue distributions to PNG
pub fn plot_analysis_png<P: AsRef<Path>>(
    analysis: &AnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    let root = BitMapBackend::new(&output_path, (config.width, config.height)).into_drawing_area();

    plot_analysis_impl(&root, analysis, config)
}

/// Internal implementation for analysis distribution plotting (for PNG)
fn plot_analysis_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    analysis: &AnalysisOutput,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError>
where
    DB::ErrorType: 'static,
{
    let bg_color = RGBColor(config.background.0, config.background.1, config.background.2);
    root.fill(&bg_color).map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Get surface data
    let (x_upper, cp_upper, vel_upper) = analysis.upper_surface();
    let (x_lower, cp_lower, vel_lower) = analysis.lower_surface();

    // Get airfoil outline
    let airfoil_x = &analysis.geometry.x;
    let airfoil_y = &analysis.geometry.y;

    // Calculate data ranges for Cp (note: inverted y-axis for negative convention)
    let cp_min = cp_upper
        .iter()
        .chain(cp_lower.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cp_max = cp_upper
        .iter()
        .chain(cp_lower.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cp_range = cp_max - cp_min;
    let cp_padding = 0.1 * cp_range;

    // Calculate data ranges for Ue
    let vel_min = vel_upper
        .iter()
        .chain(vel_lower.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let vel_max = vel_upper
        .iter()
        .chain(vel_lower.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let vel_range = vel_max - vel_min;
    let vel_padding = 0.1 * vel_range;

    // Airfoil y range
    let y_min = airfoil_y.iter().cloned().fold(f64::INFINITY, f64::min);
    let y_max = airfoil_y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Split into upper and lower regions with margin for title
    let title = config
        .title
        .clone()
        .unwrap_or_else(|| format!("{} at α = {:.1}°", analysis.foil, analysis.results.alpha_deg));

    let (upper, lower) = root.split_vertically(config.height / 2);

    // Colors
    let upper_color = RGBColor(config.upper_color.0, config.upper_color.1, config.upper_color.2);
    let lower_color = RGBColor(config.lower_color.0, config.lower_color.1, config.lower_color.2);
    let airfoil_color = RGBColor(config.airfoil_color.0, config.airfoil_color.1, config.airfoil_color.2);

    // ===== Cp Plot (upper panel) =====
    // Y-axis is INVERTED: more negative values at top (suction), positive at bottom (pressure)
    let cp_y_min = cp_min - cp_padding; // This will be at the TOP (most negative = most suction)
    let cp_y_max = cp_max + cp_padding; // This will be at the BOTTOM (most positive = pressure)

    let mut cp_chart = ChartBuilder::on(&upper)
        .caption(&title, ("sans-serif", 18))
        .margin(10)
        .x_label_area_size(35)
        .y_label_area_size(50)
        .build_cartesian_2d(
            -0.05..1.1,
            cp_y_max..cp_y_min, // INVERTED: positive at bottom, negative at top
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

    let airfoil_cp_points: Vec<(f64, f64)> = airfoil_x
        .iter()
        .zip(airfoil_y.iter())
        .map(|(&x, &y)| (x, airfoil_offset - (y - y_min) * airfoil_scale))
        .collect();

    cp_chart
        .draw_series(LineSeries::new(airfoil_cp_points, airfoil_color.stroke_width(1)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw Cp distributions
    let upper_cp_points: Vec<(f64, f64)> = x_upper.iter().zip(cp_upper.iter()).map(|(&x, &cp)| (x, cp)).collect();

    let lower_cp_points: Vec<(f64, f64)> = x_lower.iter().zip(cp_lower.iter()).map(|(&x, &cp)| (x, cp)).collect();

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

    let airfoil_vel_points: Vec<(f64, f64)> = airfoil_x
        .iter()
        .zip(airfoil_y.iter())
        .map(|(&x, &y)| (x, airfoil_vel_offset + (y - y_min) * airfoil_vel_scale))
        .collect();

    vel_chart
        .draw_series(LineSeries::new(airfoil_vel_points, airfoil_color.stroke_width(1)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    // Draw Ue distributions
    let upper_vel_points: Vec<(f64, f64)> = x_upper.iter().zip(vel_upper.iter()).map(|(&x, &v)| (x, v)).collect();

    let lower_vel_points: Vec<(f64, f64)> = x_lower.iter().zip(vel_lower.iter()).map(|(&x, &v)| (x, v)).collect();

    vel_chart
        .draw_series(LineSeries::new(upper_vel_points, upper_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    vel_chart
        .draw_series(LineSeries::new(lower_vel_points, lower_color.stroke_width(2)))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    root.present().map_err(|e| PlotError::Drawing(e.to_string()))?;

    Ok(())
}

/// Write high-precision SVG for analysis plots
///
/// Unlike the plotters SVGBackend which rounds to integer pixels,
/// this outputs coordinates with 6 decimal places for smooth curves.
fn write_analysis_precision_svg<P: AsRef<Path>>(
    analysis: &AnalysisOutput,
    output_path: P,
    config: &AnalysisPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Get surface data
    let (x_upper, cp_upper, vel_upper) = analysis.upper_surface();
    let (x_lower, cp_lower, vel_lower) = analysis.lower_surface();

    // Get airfoil outline
    let airfoil_x = &analysis.geometry.x;
    let airfoil_y = &analysis.geometry.y;

    // Calculate data ranges for Cp
    let cp_min = cp_upper
        .iter()
        .chain(cp_lower.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cp_max = cp_upper
        .iter()
        .chain(cp_lower.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cp_range = cp_max - cp_min;
    let cp_padding = 0.1 * cp_range;
    let cp_y_min = cp_min - cp_padding;
    let cp_y_max = cp_max + cp_padding;

    // Calculate data ranges for Ue
    let vel_min = vel_upper
        .iter()
        .chain(vel_lower.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let vel_max = vel_upper
        .iter()
        .chain(vel_lower.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
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

    let to_svg_x = |x: f64| -> f64 { margin_left + (x - x_data_min) / x_data_range * plot_width };
    // Cp: inverted y-axis (negative at top, positive at bottom)
    let to_svg_cp_y = |cp: f64| -> f64 { cp_plot_top + (cp - cp_y_min) / cp_data_range * plot_height };

    // Coordinate transforms for Ue plot (bottom)
    let vel_plot_top = cp_plot_bottom + gap;
    let vel_plot_bottom = vel_plot_top + plot_height;
    let vel_data_range = vel_y_max - vel_y_min;

    let to_svg_vel_y = |vel: f64| -> f64 { vel_plot_bottom - (vel - vel_y_min) / vel_data_range * plot_height };

    // Colors
    let upper_color = format!(
        "rgb({},{},{})",
        config.upper_color.0, config.upper_color.1, config.upper_color.2
    );
    let lower_color = format!(
        "rgb({},{},{})",
        config.lower_color.0, config.lower_color.1, config.lower_color.2
    );
    let airfoil_color = format!(
        "rgb({},{},{})",
        config.airfoil_color.0, config.airfoil_color.1, config.airfoil_color.2
    );
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
    let title = config
        .title
        .clone()
        .unwrap_or_else(|| format!("{} at α = {:.1}°", analysis.foil, analysis.results.alpha_deg));
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
            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, cp_plot_top, sx, cp_plot_bottom
            )?;
        }
    }
    // Horizontal gridlines (Cp)
    let cp_step = nice_step(cp_data_range, 6);
    let cp_start = (cp_y_min / cp_step).floor() * cp_step;
    let mut cp_val = cp_start;
    while cp_val <= cp_y_max {
        let sy = to_svg_cp_y(cp_val);
        if sy > cp_plot_top && sy < cp_plot_bottom {
            writeln!(
                file,
                r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                margin_left,
                sy,
                margin_left + plot_width,
                sy
            )?;
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
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
                margin_left - 5.0,
                sy + 4.0,
                cp_val
            )?;
        }
        cp_val += cp_step;
    }
    writeln!(file, "</g>")?;

    // X tick labels for Cp plot
    writeln!(
        file,
        r#"<g font-family="sans-serif" font-size="10" text-anchor="middle">"#
    )?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(
            file,
            r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
            sx,
            cp_plot_bottom + 15.0,
            x
        )?;
    }
    writeln!(file, "</g>")?;

    // Airfoil outline on Cp plot (scaled to ~1/3 height, centered at Cp=0)
    let cp_display_range_val = cp_y_max - cp_y_min;
    let airfoil_scale = (cp_display_range_val / 3.0) / (af_y_max - af_y_min);

    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="1" points=""#,
        airfoil_color
    )?;
    for (i, (&x, &y)) in airfoil_x.iter().zip(airfoil_y.iter()).enumerate() {
        // Map airfoil y directly to Cp coordinates (y=0 on airfoil -> Cp=0)
        let cp_y = -y * airfoil_scale;
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp_y);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Upper surface Cp (blue)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        upper_color
    )?;
    for (i, (&x, &cp)) in x_upper.iter().zip(cp_upper.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Lower surface Cp (red)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        lower_color
    )?;
    for (i, (&x, &cp)) in x_lower.iter().zip(cp_lower.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_cp_y(cp);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Legend for Cp plot
    let legend_x = margin_left + plot_width - 120.0;
    let legend_y = cp_plot_top + 15.0;
    writeln!(
        file,
        r#"<rect x="{:.1}" y="{:.1}" width="110" height="45" fill="white" fill-opacity="0.8" stroke="black" stroke-width="0.5"/>"#,
        legend_x, legend_y
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0,
        legend_y + 15.0,
        legend_x + 30.0,
        legend_y + 15.0,
        upper_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">Upper surface</text>"#,
        legend_x + 35.0,
        legend_y + 19.0
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0,
        legend_y + 32.0,
        legend_x + 30.0,
        legend_y + 32.0,
        lower_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">Lower surface</text>"#,
        legend_x + 35.0,
        legend_y + 36.0
    )?;

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
            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, vel_plot_top, sx, vel_plot_bottom
            )?;
        }
    }
    // Horizontal gridlines (Ue)
    let vel_step = nice_step(vel_data_range, 6);
    let vel_start = (vel_y_min / vel_step).floor() * vel_step;
    let mut vel_val = vel_start;
    while vel_val <= vel_y_max {
        let sy = to_svg_vel_y(vel_val);
        if sy > vel_plot_top && sy < vel_plot_bottom {
            writeln!(
                file,
                r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                margin_left,
                sy,
                margin_left + plot_width,
                sy
            )?;
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
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.2}</text>"#,
                margin_left - 5.0,
                sy + 4.0,
                vel_val
            )?;
        }
        vel_val += vel_step;
    }
    writeln!(file, "</g>")?;

    // X tick labels for Ue plot
    writeln!(
        file,
        r#"<g font-family="sans-serif" font-size="10" text-anchor="middle">"#
    )?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(
            file,
            r#"<text x="{:.1}" y="{:.1}">{:.1}</text>"#,
            sx,
            vel_plot_bottom + 15.0,
            x
        )?;
    }
    writeln!(file, "</g>")?;

    // Airfoil outline on Ue plot (scaled to ~1/3 height, centered at Ue=0)
    let vel_display_range_val = vel_y_max - vel_y_min;
    let airfoil_vel_scale = (vel_display_range_val / 3.0) / (af_y_max - af_y_min);

    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="1" points=""#,
        airfoil_color
    )?;
    for (i, (&x, &y)) in airfoil_x.iter().zip(airfoil_y.iter()).enumerate() {
        // Map airfoil y directly to Ue coordinates (y=0 on airfoil -> Ue=0)
        let vel_y = y * airfoil_vel_scale;
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel_y);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Upper surface Ue (blue)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        upper_color
    )?;
    for (i, (&x, &vel)) in x_upper.iter().zip(vel_upper.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Lower surface Ue (red)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        lower_color
    )?;
    for (i, (&x, &vel)) in x_lower.iter().zip(vel_lower.iter()).enumerate() {
        let sx = to_svg_x(x);
        let sy = to_svg_vel_y(vel);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    writeln!(file, "</svg>")?;

    Ok(())
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
        let alpha: Vec<f64> = polar.results.iter().map(|p| p.alpha_deg).collect();
        let cl: Vec<f64> = polar.results.iter().map(|p| p.cl).collect();
        let cd: Vec<f64> = polar.results.iter().map(|p| p.cd.unwrap_or(0.0)).collect();
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
    let alpha_min = polar_1
        .alpha
        .iter()
        .chain(polar_2.alpha.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let alpha_max = polar_1
        .alpha
        .iter()
        .chain(polar_2.alpha.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let alpha_range = alpha_max - alpha_min;
    let alpha_padding = 0.05 * alpha_range;

    let cl_min = polar_1
        .cl
        .iter()
        .chain(polar_2.cl.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cl_max = polar_1
        .cl
        .iter()
        .chain(polar_2.cl.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cl_range = cl_max - cl_min;
    let cl_padding = 0.1 * cl_range;

    let cd_min = polar_1
        .cd
        .iter()
        .chain(polar_2.cd.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cd_max = polar_1
        .cd
        .iter()
        .chain(polar_2.cd.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
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
    let to_svg_alpha_cl =
        |alpha: f64| -> f64 { cl_plot_left + (alpha - alpha_min_p) / (alpha_max_p - alpha_min_p) * plot_width };
    let to_svg_cl =
        |cl: f64| -> f64 { margin_top + plot_height - (cl - cl_min_p) / (cl_max_p - cl_min_p) * plot_height };
    let to_svg_alpha_cd =
        |alpha: f64| -> f64 { cd_plot_left + (alpha - alpha_min_p) / (alpha_max_p - alpha_min_p) * plot_width };
    let to_svg_cd =
        |cd: f64| -> f64 { margin_top + plot_height - (cd - cd_min_p) / (cd_max_p - cd_min_p) * plot_height };

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
            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx,
                margin_top,
                sx,
                margin_top + plot_height
            )?;
        }
        alpha_val += alpha_step;
    }
    let cl_step = nice_step(cl_max_p - cl_min_p, 8);
    let cl_start = (cl_min_p / cl_step).ceil() * cl_step;
    let mut cl_val = cl_start;
    while cl_val <= cl_max_p {
        let sy = to_svg_cl(cl_val);
        if sy > margin_top && sy < margin_top + plot_height {
            writeln!(
                file,
                r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                cl_plot_left,
                sy,
                cl_plot_left + plot_width,
                sy
            )?;
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
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.2}</text>"#,
                cl_plot_left - 5.0,
                sy + 4.0,
                cl_val
            )?;
        }
        cl_val += cl_step;
    }
    writeln!(file, "</g>")?;

    writeln!(
        file,
        r#"<g font-family="sans-serif" font-size="11" text-anchor="middle">"#
    )?;
    alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cl(alpha_val);
        if sx > cl_plot_left + 10.0 && sx < cl_plot_left + plot_width - 10.0 {
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.0}</text>"#,
                sx,
                margin_top + plot_height + 18.0,
                alpha_val
            )?;
        }
        alpha_val += alpha_step;
    }
    writeln!(file, "</g>")?;

    // Data series 1 (CL)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_1
    )?;
    for (i, (&alpha, &cl)) in polar_1.alpha.iter().zip(polar_1.cl.iter()).enumerate() {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 1
    for (&alpha, &cl) in polar_1.alpha.iter().zip(polar_1.cl.iter()) {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        writeln!(
            file,
            r#"<circle cx="{:.6}" cy="{:.6}" r="3" fill="{}" />"#,
            sx, sy, color_1
        )?;
    }

    // Data series 2 (CL)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_2
    )?;
    for (i, (&alpha, &cl)) in polar_2.alpha.iter().zip(polar_2.cl.iter()).enumerate() {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 2
    for (&alpha, &cl) in polar_2.alpha.iter().zip(polar_2.cl.iter()) {
        let sx = to_svg_alpha_cl(alpha);
        let sy = to_svg_cl(cl);
        writeln!(
            file,
            r#"<rect x="{:.6}" y="{:.6}" width="6" height="6" fill="{}" transform="rotate(45 {:.6} {:.6})"/>"#,
            sx - 3.0,
            sy - 3.0,
            color_2,
            sx,
            sy
        )?;
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
            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx,
                margin_top,
                sx,
                margin_top + plot_height
            )?;
        }
        alpha_val += alpha_step;
    }
    let cd_step = nice_step(cd_max_p - cd_min_p, 8);
    let cd_start = (cd_min_p / cd_step).ceil() * cd_step;
    let mut cd_val = cd_start;
    while cd_val <= cd_max_p {
        let sy = to_svg_cd(cd_val);
        if sy > margin_top && sy < margin_top + plot_height {
            writeln!(
                file,
                r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                cd_plot_left,
                sy,
                cd_plot_left + plot_width,
                sy
            )?;
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
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.4}</text>"#,
                cd_plot_left - 5.0,
                sy + 4.0,
                cd_val
            )?;
        }
        cd_val += cd_step;
    }
    writeln!(file, "</g>")?;

    writeln!(
        file,
        r#"<g font-family="sans-serif" font-size="11" text-anchor="middle">"#
    )?;
    alpha_val = alpha_start;
    while alpha_val <= alpha_max_p {
        let sx = to_svg_alpha_cd(alpha_val);
        if sx > cd_plot_left + 10.0 && sx < cd_plot_left + plot_width - 10.0 {
            writeln!(
                file,
                r#"<text x="{:.1}" y="{:.1}">{:.0}</text>"#,
                sx,
                margin_top + plot_height + 18.0,
                alpha_val
            )?;
        }
        alpha_val += alpha_step;
    }
    writeln!(file, "</g>")?;

    // Data series 1 (CD)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_1
    )?;
    for (i, (&alpha, &cd)) in polar_1.alpha.iter().zip(polar_1.cd.iter()).enumerate() {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 1
    for (&alpha, &cd) in polar_1.alpha.iter().zip(polar_1.cd.iter()) {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        writeln!(
            file,
            r#"<circle cx="{:.6}" cy="{:.6}" r="3" fill="{}" />"#,
            sx, sy, color_1
        )?;
    }

    // Data series 2 (CD)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_2
    )?;
    for (i, (&alpha, &cd)) in polar_2.alpha.iter().zip(polar_2.cd.iter()).enumerate() {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Data markers for series 2
    for (&alpha, &cd) in polar_2.alpha.iter().zip(polar_2.cd.iter()) {
        let sx = to_svg_alpha_cd(alpha);
        let sy = to_svg_cd(cd);
        writeln!(
            file,
            r#"<rect x="{:.6}" y="{:.6}" width="6" height="6" fill="{}" transform="rotate(45 {:.6} {:.6})"/>"#,
            sx - 3.0,
            sy - 3.0,
            color_2,
            sx,
            sy
        )?;
    }

    // Legend (centered at bottom)
    let legend_x = width / 2.0 - 100.0;
    let legend_y = height - 25.0;

    // Series 1
    writeln!(
        file,
        r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{}" />"#,
        legend_x, legend_y, color_1
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x - 15.0,
        legend_y,
        legend_x + 15.0,
        legend_y,
        color_1
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x + 25.0,
        legend_y + 4.0,
        config.label_1
    )?;

    // Series 2
    let legend_x2 = legend_x + 120.0;
    writeln!(
        file,
        r#"<rect x="{:.1}" y="{:.1}" width="8" height="8" fill="{}" transform="rotate(45 {:.1} {:.1})"/>"#,
        legend_x2 - 4.0,
        legend_y - 4.0,
        color_2,
        legend_x2,
        legend_y
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x2 - 15.0,
        legend_y,
        legend_x2 + 15.0,
        legend_y,
        color_2
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x2 + 25.0,
        legend_y + 4.0,
        config.label_2
    )?;

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
    let mut combined: Vec<(f64, f64, f64)> = neg_alpha
        .into_iter()
        .zip(neg_cl)
        .zip(neg_cd)
        .map(|((a, c), d)| (a, c, d))
        .collect();

    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Remove duplicates (keep first occurrence of each alpha)
    let mut seen_alphas = std::collections::HashSet::new();
    let combined: Vec<(f64, f64, f64)> = combined
        .into_iter()
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
            series1_color: (0, 100, 200), // Blue for XFOIL
            series2_color: (200, 50, 50), // Red for YFoil
            marker_size: 4.0,
            background: (255, 255, 255), // White
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
    let x_min = coords1
        .iter()
        .chain(coords2.iter())
        .map(|(x, _)| *x)
        .fold(f64::INFINITY, f64::min);
    let x_max = coords1
        .iter()
        .chain(coords2.iter())
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = coords1
        .iter()
        .chain(coords2.iter())
        .map(|(_, y)| *y)
        .fold(f64::INFINITY, f64::min);
    let y_max = coords1
        .iter()
        .chain(coords2.iter())
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
    let color_1 = format!(
        "rgb({},{},{})",
        config.series1_color.0, config.series1_color.1, config.series1_color.2
    );
    let color_2 = format!(
        "rgb({},{},{})",
        config.series2_color.0, config.series2_color.1, config.series2_color.2
    );
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
    let title = config
        .title
        .clone()
        .unwrap_or_else(|| "Geometry Comparison".to_string());
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
        writeln!(
            file,
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}"/>"#,
            sx,
            margin,
            sx,
            config.height as f64 - margin
        )?;
    }
    for i in 0..=5 {
        let y = y_min_padded + (i as f64 / 5.0) * data_height;
        let sy = to_svg_y(y);
        writeln!(
            file,
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}"/>"#,
            margin,
            sy,
            config.width as f64 - margin,
            sy
        )?;
    }
    writeln!(file, "</g>")?;

    // Plot border
    writeln!(
        file,
        r##"<rect x="{}" y="{}" width="{}" height="{}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin, margin, plot_width, plot_height
    )?;

    // Series 1 (XFOIL) - blue line with circle markers
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_1
    )?;
    for (i, (x, y)) in coords1.iter().enumerate() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Circle markers for series 1
    for (x, y) in coords1.iter() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        writeln!(
            file,
            r#"<circle cx="{:.6}" cy="{:.6}" r="{}" fill="none" stroke="{}" stroke-width="1.5"/>"#,
            sx, sy, config.marker_size, color_1
        )?;
    }

    // Series 2 (YFoil) - red line with cross markers
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        color_2
    )?;
    for (i, (x, y)) in coords2.iter().enumerate() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", sx, sy)?;
    }
    writeln!(file, r#""/>"#)?;

    // Cross markers for series 2
    let m = config.marker_size;
    for (x, y) in coords2.iter() {
        let sx = to_svg_x(*x);
        let sy = to_svg_y(*y);
        writeln!(
            file,
            r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
            sx - m,
            sy - m,
            sx + m,
            sy + m,
            color_2
        )?;
        writeln!(
            file,
            r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
            sx - m,
            sy + m,
            sx + m,
            sy - m,
            color_2
        )?;
    }

    // Legend
    let legend_x = config.width as f64 - 150.0;
    let legend_y = margin + 15.0;
    writeln!(
        file,
        r#"<rect x="{:.1}" y="{:.1}" width="130" height="50" fill="white" fill-opacity="0.9" stroke="black" stroke-width="0.5"/>"#,
        legend_x, legend_y
    )?;

    // Series 1 legend
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0,
        legend_y + 15.0,
        legend_x + 30.0,
        legend_y + 15.0,
        color_1
    )?;
    writeln!(
        file,
        r#"<circle cx="{:.1}" cy="{:.1}" r="{}" fill="none" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0,
        legend_y + 15.0,
        config.marker_size,
        color_1
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">{}</text>"#,
        legend_x + 40.0,
        legend_y + 19.0,
        config.label_1
    )?;

    // Series 2 legend
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0,
        legend_y + 35.0,
        legend_x + 30.0,
        legend_y + 35.0,
        color_2
    )?;
    writeln!(
        file,
        r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0 - m,
        legend_y + 35.0 - m,
        legend_x + 20.0 + m,
        legend_y + 35.0 + m,
        color_2
    )?;
    writeln!(
        file,
        r#"<line x1="{:.6}" y1="{:.6}" x2="{:.6}" y2="{:.6}" stroke="{}" stroke-width="1.5"/>"#,
        legend_x + 20.0 - m,
        legend_y + 35.0 + m,
        legend_x + 20.0 + m,
        legend_y + 35.0 - m,
        color_2
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">{}</text>"#,
        legend_x + 40.0,
        legend_y + 39.0,
        config.label_2
    )?;

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

// ============================================================================
// Extended Polar Data (with CM)
// ============================================================================

/// Extended polar data including moment coefficient
#[derive(Debug, Clone)]
pub struct PolarDataWithCm {
    /// Angle of attack values (degrees)
    pub alpha: Vec<f64>,
    /// Lift coefficient values
    pub cl: Vec<f64>,
    /// Drag coefficient values
    pub cd: Vec<f64>,
    /// Moment coefficient values (about quarter chord)
    pub cm: Vec<f64>,
}

impl PolarDataWithCm {
    /// Create from PolarOutput
    pub fn from_polar_output(polar: &crate::output::PolarOutput) -> Self {
        Self {
            alpha: polar.results.iter().map(|p| p.alpha_deg).collect(),
            cl: polar.results.iter().map(|p| p.cl).collect(),
            cd: polar.results.iter().map(|p| p.cd.unwrap_or(0.0)).collect(),
            cm: polar.results.iter().map(|p| p.cm).collect(),
        }
    }

    /// Create from raw data vectors
    pub fn new(alpha: Vec<f64>, cl: Vec<f64>, cd: Vec<f64>, cm: Vec<f64>) -> Self {
        Self { alpha, cl, cd, cm }
    }

    /// Convert to basic PolarData (without CM)
    pub fn to_polar_data(&self) -> PolarData {
        PolarData::new(self.alpha.clone(), self.cl.clone(), self.cd.clone())
    }
}

/// Parse an XFOIL polar file and extract alpha, CL, CD, CM data
pub fn parse_xfoil_polar_file_with_cm<P: AsRef<Path>>(path: P) -> Result<PolarDataWithCm, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let mut alphas = Vec::new();
    let mut cls = Vec::new();
    let mut cds = Vec::new();
    let mut cms = Vec::new();

    let mut in_data = false;
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("------") {
            in_data = true;
            continue;
        }
        if in_data && !line.is_empty() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 5 {
                if let (Ok(alpha), Ok(cl), Ok(cd), Ok(cm)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[4].parse::<f64>(), // CM is column 5 (index 4)
                ) {
                    alphas.push(alpha);
                    cls.push(cl);
                    cds.push(cd);
                    cms.push(cm);
                }
            }
        }
    }

    Ok(PolarDataWithCm::new(alphas, cls, cds, cms))
}

/// Parse a YFoil polar JSON file into PolarDataWithCm
pub fn parse_yfoil_polar_file_with_cm<P: AsRef<Path>>(path: P) -> Result<PolarDataWithCm, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let polar: crate::output::PolarOutput = serde_json::from_str(&content)
        .map_err(|e| PlotError::Drawing(format!("Failed to parse YFoil polar JSON: {}", e)))?;
    Ok(PolarDataWithCm::from_polar_output(&polar))
}

/// Stitch two polars (with CM) together in ascending alpha order
pub fn stitch_polars_with_cm(polar_pos: &PolarDataWithCm, polar_neg: &PolarDataWithCm) -> PolarDataWithCm {
    let mut neg_alpha: Vec<f64> = polar_neg.alpha.iter().copied().rev().collect();
    let mut neg_cl: Vec<f64> = polar_neg.cl.iter().copied().rev().collect();
    let mut neg_cd: Vec<f64> = polar_neg.cd.iter().copied().rev().collect();
    let mut neg_cm: Vec<f64> = polar_neg.cm.iter().copied().rev().collect();

    neg_alpha.extend(polar_pos.alpha.iter().copied());
    neg_cl.extend(polar_pos.cl.iter().copied());
    neg_cd.extend(polar_pos.cd.iter().copied());
    neg_cm.extend(polar_pos.cm.iter().copied());

    let mut combined: Vec<(f64, f64, f64, f64)> = neg_alpha
        .into_iter()
        .zip(neg_cl)
        .zip(neg_cd)
        .zip(neg_cm)
        .map(|(((a, c), d), m)| (a, c, d, m))
        .collect();

    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut seen = std::collections::HashSet::new();
    let combined: Vec<_> = combined
        .into_iter()
        .filter(|(a, _, _, _)| seen.insert((*a * 1000.0).round() as i64))
        .collect();

    PolarDataWithCm::new(
        combined.iter().map(|(a, _, _, _)| *a).collect(),
        combined.iter().map(|(_, c, _, _)| *c).collect(),
        combined.iter().map(|(_, _, d, _)| *d).collect(),
        combined.iter().map(|(_, _, _, m)| *m).collect(),
    )
}

// ============================================================================
// XFOIL DUMP Parsing
// ============================================================================

/// Boundary layer data from XFOIL DUMP command
#[derive(Debug, Clone)]
pub struct XfoilDump {
    /// Arc length
    pub s: Vec<f64>,
    /// X coordinate
    pub x: Vec<f64>,
    /// Y coordinate
    pub y: Vec<f64>,
    /// Edge velocity (Ue/Vinf)
    pub ue: Vec<f64>,
    /// Displacement thickness
    pub dstar: Vec<f64>,
    /// Momentum thickness
    pub theta: Vec<f64>,
    /// Skin friction coefficient
    pub cf: Vec<f64>,
    /// Shape factor (H = δ*/θ)
    pub h: Vec<f64>,
    /// Energy shape factor (H*)
    pub hs: Vec<f64>,
}

impl XfoilDump {
    /// Create empty XfoilDump
    pub fn new() -> Self {
        Self {
            s: Vec::new(),
            x: Vec::new(),
            y: Vec::new(),
            ue: Vec::new(),
            dstar: Vec::new(),
            theta: Vec::new(),
            cf: Vec::new(),
            h: Vec::new(),
            hs: Vec::new(),
        }
    }
}

impl Default for XfoilDump {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse an XFOIL DUMP file
///
/// DUMP file format:
/// #    s        x        y     Ue/Vinf    Dstar     Theta      Cf       H       H*   ...
pub fn parse_xfoil_dump_file<P: AsRef<Path>>(path: P) -> Result<XfoilDump, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let mut dump = XfoilDump::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 9 {
            if let (Ok(s), Ok(x), Ok(y), Ok(ue), Ok(dstar), Ok(theta), Ok(cf), Ok(h), Ok(hs)) = (
                parts[0].parse::<f64>(),
                parts[1].parse::<f64>(),
                parts[2].parse::<f64>(),
                parts[3].parse::<f64>(),
                parts[4].parse::<f64>(),
                parts[5].parse::<f64>(),
                parts[6].parse::<f64>(),
                parts[7].parse::<f64>(),
                parts[8].parse::<f64>(),
            ) {
                dump.s.push(s);
                dump.x.push(x);
                dump.y.push(y);
                dump.ue.push(ue);
                dump.dstar.push(dstar);
                dump.theta.push(theta);
                dump.cf.push(cf);
                dump.h.push(h);
                dump.hs.push(hs);
            }
        }
    }

    Ok(dump)
}

// ============================================================================
// XFOIL CP Parsing
// ============================================================================

/// Pressure coefficient data from XFOIL CPWR command
#[derive(Debug, Clone)]
pub struct XfoilCp {
    /// X coordinate
    pub x: Vec<f64>,
    /// Pressure coefficient
    pub cp: Vec<f64>,
}

impl XfoilCp {
    /// Create empty XfoilCp
    pub fn new() -> Self {
        Self {
            x: Vec::new(),
            cp: Vec::new(),
        }
    }
}

impl Default for XfoilCp {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse an XFOIL CPWR file
///
/// CPWR file format:
/// #      x          Cp
pub fn parse_xfoil_cp_file<P: AsRef<Path>>(path: P) -> Result<XfoilCp, PlotError> {
    let content = std::fs::read_to_string(path)?;
    let mut cp_data = XfoilCp::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(x), Ok(cp)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                cp_data.x.push(x);
                cp_data.cp.push(cp);
            }
        }
    }

    Ok(cp_data)
}

// ============================================================================
// 3-Panel Polar Plot (CL, CD, CM)
// ============================================================================

/// Plot a polar comparison with 3 horizontal panels: CL vs α, CD vs α, CM vs α
pub fn plot_polar_3panel_svg<P: AsRef<Path>>(
    polar_1: &PolarDataWithCm,
    polar_2: &PolarDataWithCm,
    output_path: P,
    config: &PolarPlotConfig,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    // Calculate data ranges for alpha (shared x-axis)
    let alpha_min = polar_1
        .alpha
        .iter()
        .chain(polar_2.alpha.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let alpha_max = polar_1
        .alpha
        .iter()
        .chain(polar_2.alpha.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let alpha_range = alpha_max - alpha_min;
    let alpha_padding = 0.05 * alpha_range.max(1.0);
    let alpha_min_p = alpha_min - alpha_padding;
    let alpha_max_p = alpha_max + alpha_padding;

    // CL range
    let cl_min = polar_1
        .cl
        .iter()
        .chain(polar_2.cl.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cl_max = polar_1
        .cl
        .iter()
        .chain(polar_2.cl.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cl_range = (cl_max - cl_min).max(0.1);
    let cl_padding = 0.1 * cl_range;
    let cl_min_p = cl_min - cl_padding;
    let cl_max_p = cl_max + cl_padding;

    // CD range
    let cd_min = polar_1
        .cd
        .iter()
        .chain(polar_2.cd.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cd_max = polar_1
        .cd
        .iter()
        .chain(polar_2.cd.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cd_range = (cd_max - cd_min).max(0.001);
    let cd_padding = 0.1 * cd_range;
    let cd_min_p = cd_min - cd_padding;
    let cd_max_p = cd_max + cd_padding;

    // CM range
    let cm_min = polar_1
        .cm
        .iter()
        .chain(polar_2.cm.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cm_max = polar_1
        .cm
        .iter()
        .chain(polar_2.cm.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cm_range = (cm_max - cm_min).max(0.01);
    let cm_padding = 0.1 * cm_range;
    let cm_min_p = cm_min - cm_padding;
    let cm_max_p = cm_max + cm_padding;

    // Layout: 3 panels side by side
    let width = config.width as f64;
    let height = config.height as f64;
    let margin_left = 70.0;
    let margin_right = 20.0;
    let margin_top = 50.0;
    let margin_bottom = 60.0;
    let gap = 50.0;

    let total_gap = 2.0 * gap;
    let plot_width = (width - margin_left - margin_right - total_gap) / 3.0;
    let plot_height = height - margin_top - margin_bottom;

    // Panel x positions
    let panel1_x = margin_left;
    let panel2_x = margin_left + plot_width + gap;
    let panel3_x = margin_left + 2.0 * (plot_width + gap);

    // Colors
    let color_1 = format!("rgb({},{},{})", config.color_1.0, config.color_1.1, config.color_1.2);
    let color_2 = format!("rgb({},{},{})", config.color_2.0, config.color_2.1, config.color_2.2);
    let grid_color = "#DDDDDD";

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

    // Helper function for coordinate transforms
    let to_svg_x =
        |x: f64, panel_x: f64| -> f64 { panel_x + (x - alpha_min_p) / (alpha_max_p - alpha_min_p) * plot_width };

    // Draw each panel
    let panels = [
        (panel1_x, cl_min_p, cl_max_p, &polar_1.cl, &polar_2.cl, "CL"),
        (panel2_x, cd_min_p, cd_max_p, &polar_1.cd, &polar_2.cd, "CD"),
        (panel3_x, cm_min_p, cm_max_p, &polar_1.cm, &polar_2.cm, "CM"),
    ];

    for (panel_x, y_min, y_max, data_1, data_2, label) in panels {
        let y_range = y_max - y_min;

        let to_svg_y = |y: f64| -> f64 { margin_top + plot_height - (y - y_min) / y_range * plot_height };

        // Plot border
        writeln!(
            file,
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
            panel_x, margin_top, plot_width, plot_height
        )?;

        // Gridlines
        writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;

        // Vertical gridlines (alpha)
        let alpha_step = nice_step(alpha_max_p - alpha_min_p, 6);
        let alpha_start = (alpha_min_p / alpha_step).ceil() * alpha_step;
        let mut alpha_val = alpha_start;
        while alpha_val <= alpha_max_p {
            let sx = to_svg_x(alpha_val, panel_x);
            if sx > panel_x && sx < panel_x + plot_width {
                writeln!(
                    file,
                    r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                    sx,
                    margin_top,
                    sx,
                    margin_top + plot_height
                )?;
            }
            alpha_val += alpha_step;
        }

        // Horizontal gridlines
        let y_step = nice_step(y_range, 6);
        let y_start = (y_min / y_step).ceil() * y_step;
        let mut y_val = y_start;
        while y_val <= y_max {
            let sy = to_svg_y(y_val);
            if sy > margin_top && sy < margin_top + plot_height {
                writeln!(
                    file,
                    r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                    panel_x,
                    sy,
                    panel_x + plot_width,
                    sy
                )?;
            }
            y_val += y_step;
        }
        writeln!(file, "</g>")?;

        // Y-axis label
        writeln!(
            file,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 {:.1} {:.1})">{}</text>"#,
            panel_x - 50.0,
            margin_top + plot_height / 2.0,
            panel_x - 50.0,
            margin_top + plot_height / 2.0,
            label
        )?;

        // X-axis label
        writeln!(
            file,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">Alpha (deg)</text>"#,
            panel_x + plot_width / 2.0,
            margin_top + plot_height + 40.0
        )?;

        // Y tick labels
        writeln!(file, r#"<g font-family="sans-serif" font-size="10" text-anchor="end">"#)?;
        y_val = y_start;
        let decimals = if label == "CD" { 4 } else { 2 };
        while y_val <= y_max {
            let sy = to_svg_y(y_val);
            if sy > margin_top + 5.0 && sy < margin_top + plot_height - 5.0 {
                writeln!(
                    file,
                    r#"<text x="{:.1}" y="{:.1}">{:.decimals$}</text>"#,
                    panel_x - 5.0,
                    sy + 4.0,
                    y_val,
                    decimals = decimals
                )?;
            }
            y_val += y_step;
        }
        writeln!(file, "</g>")?;

        // X tick labels
        writeln!(
            file,
            r#"<g font-family="sans-serif" font-size="10" text-anchor="middle">"#
        )?;
        alpha_val = alpha_start;
        while alpha_val <= alpha_max_p {
            let sx = to_svg_x(alpha_val, panel_x);
            if sx > panel_x + 10.0 && sx < panel_x + plot_width - 10.0 {
                writeln!(
                    file,
                    r#"<text x="{:.1}" y="{:.1}">{:.0}</text>"#,
                    sx,
                    margin_top + plot_height + 18.0,
                    alpha_val
                )?;
            }
            alpha_val += alpha_step;
        }
        writeln!(file, "</g>")?;

        // Data series 1
        write!(
            file,
            r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
            color_1
        )?;
        for (i, (&alpha, &y)) in polar_1.alpha.iter().zip(data_1.iter()).enumerate() {
            let sx = to_svg_x(alpha, panel_x);
            let sy = to_svg_y(y);
            if i > 0 {
                write!(file, " ")?;
            }
            write!(file, "{:.6},{:.6}", sx, sy)?;
        }
        writeln!(file, r#""/>"#)?;

        // Data markers for series 1
        for (&alpha, &y) in polar_1.alpha.iter().zip(data_1.iter()) {
            let sx = to_svg_x(alpha, panel_x);
            let sy = to_svg_y(y);
            writeln!(
                file,
                r#"<circle cx="{:.6}" cy="{:.6}" r="3" fill="{}" />"#,
                sx, sy, color_1
            )?;
        }

        // Data series 2
        write!(
            file,
            r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
            color_2
        )?;
        for (i, (&alpha, &y)) in polar_2.alpha.iter().zip(data_2.iter()).enumerate() {
            let sx = to_svg_x(alpha, panel_x);
            let sy = to_svg_y(y);
            if i > 0 {
                write!(file, " ")?;
            }
            write!(file, "{:.6},{:.6}", sx, sy)?;
        }
        writeln!(file, r#""/>"#)?;

        // Data markers for series 2
        for (&alpha, &y) in polar_2.alpha.iter().zip(data_2.iter()) {
            let sx = to_svg_x(alpha, panel_x);
            let sy = to_svg_y(y);
            writeln!(
                file,
                r#"<rect x="{:.6}" y="{:.6}" width="6" height="6" fill="{}" transform="rotate(45 {:.6} {:.6})"/>"#,
                sx - 3.0,
                sy - 3.0,
                color_2,
                sx,
                sy
            )?;
        }
    }

    // Legend (centered at bottom)
    let legend_x = width / 2.0 - 100.0;
    let legend_y = height - 20.0;

    // Series 1
    writeln!(
        file,
        r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="{}" />"#,
        legend_x, legend_y, color_1
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x - 15.0,
        legend_y,
        legend_x + 15.0,
        legend_y,
        color_1
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x + 25.0,
        legend_y + 4.0,
        config.label_1
    )?;

    // Series 2
    let legend_x2 = legend_x + 120.0;
    writeln!(
        file,
        r#"<rect x="{:.1}" y="{:.1}" width="8" height="8" fill="{}" transform="rotate(45 {:.1} {:.1})"/>"#,
        legend_x2 - 4.0,
        legend_y - 4.0,
        color_2,
        legend_x2,
        legend_y
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x2 - 15.0,
        legend_y,
        legend_x2 + 15.0,
        legend_y,
        color_2
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12">{}</text>"#,
        legend_x2 + 25.0,
        legend_y + 4.0,
        config.label_2
    )?;

    writeln!(file, "</svg>")?;
    Ok(())
}

// ============================================================================
// Cp/Ue Comparison Plot (2 vertical subplots)
// ============================================================================

/// Plot Cp and Ue comparison (2 stacked vertical subplots)
pub fn plot_cp_ue_comparison_svg<P: AsRef<Path>>(
    yfoil_x: &[f64],
    yfoil_cp: &[f64],
    yfoil_ue: &[f64],
    xfoil_cp: &XfoilCp,
    xfoil_dump: &XfoilDump,
    alpha_deg: f64,
    airfoil: &str,
    output_path: P,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    let width = 1200.0;
    let height = 800.0;
    let margin_left = 70.0;
    let margin_right = 30.0;
    let margin_top = 50.0;
    let margin_bottom = 50.0;
    let gap = 30.0;

    let plot_width = width - margin_left - margin_right;
    let plot_height = (height - margin_top - margin_bottom - gap) / 2.0;

    // Cp ranges
    let cp_min = yfoil_cp
        .iter()
        .chain(xfoil_cp.cp.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let cp_max = yfoil_cp
        .iter()
        .chain(xfoil_cp.cp.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let cp_range = (cp_max - cp_min).max(0.1);
    let cp_padding = 0.1 * cp_range;
    let cp_min_p = cp_min - cp_padding;
    let cp_max_p = cp_max + cp_padding;

    // Ue ranges (using XFOIL dump for Ue)
    let ue_min = yfoil_ue
        .iter()
        .chain(xfoil_dump.ue.iter())
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let ue_max = yfoil_ue
        .iter()
        .chain(xfoil_dump.ue.iter())
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let ue_range = (ue_max - ue_min).max(0.1);
    let ue_padding = 0.1 * ue_range;
    let ue_min_p = ue_min - ue_padding;
    let ue_max_p = ue_max + ue_padding;

    // X range
    let x_min = -0.05;
    let x_max = 1.1;

    // Colors
    let yfoil_color = "rgb(0,100,200)";
    let xfoil_color = "rgb(200,50,50)";
    let grid_color = "#DDDDDD";

    // Coordinate transforms
    let cp_plot_top = margin_top;
    let cp_plot_bottom = margin_top + plot_height;
    let ue_plot_top = cp_plot_bottom + gap;
    let ue_plot_bottom = ue_plot_top + plot_height;

    let to_svg_x = |x: f64| -> f64 { margin_left + (x - x_min) / (x_max - x_min) * plot_width };
    // Cp: inverted y-axis (more negative at top)
    let to_svg_cp_y = |cp: f64| -> f64 { cp_plot_top + (cp - cp_min_p) / (cp_max_p - cp_min_p) * plot_height };
    let to_svg_ue_y = |ue: f64| -> f64 { ue_plot_bottom - (ue - ue_min_p) / (ue_max_p - ue_min_p) * plot_height };

    // SVG header
    writeln!(
        file,
        r#"<svg width="{}" height="{}" viewBox="0 0 {} {}" xmlns="http://www.w3.org/2000/svg">"#,
        width as u32, height as u32, width as u32, height as u32
    )?;
    writeln!(file, r#"<rect width="100%" height="100%" fill="white"/>"#)?;

    // Title
    writeln!(
        file,
        r#"<text x="{:.1}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold">{} at α = {:.1}°</text>"#,
        width / 2.0,
        airfoil,
        alpha_deg
    )?;

    // ===== Cp Plot (top) =====
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin_left, cp_plot_top, plot_width, plot_height
    )?;

    // Cp gridlines
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(
            file,
            r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
            sx, cp_plot_top, sx, cp_plot_bottom
        )?;
    }
    let cp_step = nice_step(cp_max_p - cp_min_p, 6);
    let mut cp_val = (cp_min_p / cp_step).ceil() * cp_step;
    while cp_val <= cp_max_p {
        let sy = to_svg_cp_y(cp_val);
        writeln!(
            file,
            r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
            margin_left,
            sy,
            margin_left + plot_width,
            sy
        )?;
        cp_val += cp_step;
    }
    writeln!(file, "</g>")?;

    // Cp axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        margin_left + plot_width / 2.0,
        cp_plot_bottom + 35.0
    )?;
    writeln!(
        file,
        r#"<text x="20" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 20 {:.1})">Cp</text>"#,
        cp_plot_top + plot_height / 2.0,
        cp_plot_top + plot_height / 2.0
    )?;

    // YFoil Cp
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        yfoil_color
    )?;
    for (i, (&x, &cp)) in yfoil_x.iter().zip(yfoil_cp.iter()).enumerate() {
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_cp_y(cp))?;
    }
    writeln!(file, r#""/>"#)?;

    // XFOIL Cp
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" stroke-dasharray="5,3" points=""#,
        xfoil_color
    )?;
    for (i, (&x, &cp)) in xfoil_cp.x.iter().zip(xfoil_cp.cp.iter()).enumerate() {
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_cp_y(cp))?;
    }
    writeln!(file, r#""/>"#)?;

    // ===== Ue Plot (bottom) =====
    writeln!(
        file,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
        margin_left, ue_plot_top, plot_width, plot_height
    )?;

    // Ue gridlines
    writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
    for i in 0..=10 {
        let x = i as f64 / 10.0;
        let sx = to_svg_x(x);
        writeln!(
            file,
            r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
            sx, ue_plot_top, sx, ue_plot_bottom
        )?;
    }
    let ue_step = nice_step(ue_max_p - ue_min_p, 6);
    let mut ue_val = (ue_min_p / ue_step).ceil() * ue_step;
    while ue_val <= ue_max_p {
        let sy = to_svg_ue_y(ue_val);
        writeln!(
            file,
            r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
            margin_left,
            sy,
            margin_left + plot_width,
            sy
        )?;
        ue_val += ue_step;
    }
    writeln!(file, "</g>")?;

    // Ue axis labels
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        margin_left + plot_width / 2.0,
        ue_plot_bottom + 35.0
    )?;
    writeln!(
        file,
        r#"<text x="20" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 20 {:.1})">Ue/U∞</text>"#,
        ue_plot_top + plot_height / 2.0,
        ue_plot_top + plot_height / 2.0
    )?;

    // YFoil Ue
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" points=""#,
        yfoil_color
    )?;
    for (i, (&x, &ue)) in yfoil_x.iter().zip(yfoil_ue.iter()).enumerate() {
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_ue_y(ue))?;
    }
    writeln!(file, r#""/>"#)?;

    // XFOIL Ue (from dump, use x coordinate)
    write!(
        file,
        r#"<polyline fill="none" stroke="{}" stroke-width="2" stroke-dasharray="5,3" points=""#,
        xfoil_color
    )?;
    for (i, (&x, &ue)) in xfoil_dump.x.iter().zip(xfoil_dump.ue.iter()).enumerate() {
        if i > 0 {
            write!(file, " ")?;
        }
        write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_ue_y(ue))?;
    }
    writeln!(file, r#""/>"#)?;

    // Legend
    let legend_x = margin_left + plot_width - 150.0;
    let legend_y = cp_plot_top + 15.0;
    writeln!(
        file,
        r#"<rect x="{:.1}" y="{:.1}" width="140" height="45" fill="white" fill-opacity="0.9" stroke="black" stroke-width="0.5"/>"#,
        legend_x, legend_y
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x + 10.0,
        legend_y + 15.0,
        legend_x + 35.0,
        legend_y + 15.0,
        yfoil_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">YFoil</text>"#,
        legend_x + 45.0,
        legend_y + 19.0
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2" stroke-dasharray="5,3"/>"#,
        legend_x + 10.0,
        legend_y + 32.0,
        legend_x + 35.0,
        legend_y + 32.0,
        xfoil_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">XFOIL</text>"#,
        legend_x + 45.0,
        legend_y + 36.0
    )?;

    writeln!(file, "</svg>")?;
    Ok(())
}

// ============================================================================
// BL Comparison Plot (6 vertical subplots)
// ============================================================================

/// YFoil BL distribution data for plotting
#[derive(Debug, Clone)]
pub struct YfoilBLDist {
    /// X coordinate
    pub x: Vec<f64>,
    /// Arc length
    pub s: Vec<f64>,
    /// Momentum thickness
    pub theta: Vec<f64>,
    /// Displacement thickness
    pub dstar: Vec<f64>,
    /// Shape factor
    pub h: Vec<f64>,
    /// Energy shape factor
    pub hs: Vec<f64>,
    /// Skin friction coefficient
    pub cf: Vec<f64>,
    /// Edge velocity
    pub ue: Vec<f64>,
}

/// Plot BL comparison with 6 stacked vertical subplots
pub fn plot_bl_comparison_svg<P: AsRef<Path>>(
    yfoil_bl: &YfoilBLDist,
    xfoil_dump: &XfoilDump,
    alpha_deg: f64,
    airfoil: &str,
    output_path: P,
) -> Result<(), PlotError> {
    let mut file = std::fs::File::create(output_path)?;

    let width = 1200.0;
    let height = 1000.0; // Taller for 6 panels
    let margin_left = 80.0;
    let margin_right = 30.0;
    let margin_top = 50.0;
    let margin_bottom = 40.0;
    let gap = 10.0;
    let n_panels = 6;

    let plot_width = width - margin_left - margin_right;
    let total_gap = (n_panels - 1) as f64 * gap;
    let plot_height = (height - margin_top - margin_bottom - total_gap) / n_panels as f64;

    // X range (using x/c)
    let x_min = -0.05;
    let x_max = 1.1;

    let to_svg_x = |x: f64| -> f64 { margin_left + (x - x_min) / (x_max - x_min) * plot_width };

    // Colors
    let yfoil_color = "rgb(0,100,200)";
    let xfoil_color = "rgb(200,50,50)";
    let grid_color = "#DDDDDD";

    // Define panels: (yfoil data, xfoil data, label, index)
    let panels: Vec<(&[f64], &[f64], &str)> = vec![
        (&yfoil_bl.theta, &xfoil_dump.theta, "θ"),
        (&yfoil_bl.dstar, &xfoil_dump.dstar, "δ*"),
        (&yfoil_bl.h, &xfoil_dump.h, "H"),
        (&yfoil_bl.hs, &xfoil_dump.hs, "H*"),
        (&yfoil_bl.cf, &xfoil_dump.cf, "Cf"),
        (&yfoil_bl.ue, &xfoil_dump.ue, "Ue/U∞"),
    ];

    // SVG header
    writeln!(
        file,
        r#"<svg width="{}" height="{}" viewBox="0 0 {} {}" xmlns="http://www.w3.org/2000/svg">"#,
        width as u32, height as u32, width as u32, height as u32
    )?;
    writeln!(file, r#"<rect width="100%" height="100%" fill="white"/>"#)?;

    // Title
    writeln!(
        file,
        r#"<text x="{:.1}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold">{} BL Variables at α = {:.1}°</text>"#,
        width / 2.0,
        airfoil,
        alpha_deg
    )?;

    for (panel_idx, (yfoil_data, xfoil_data, label)) in panels.iter().enumerate() {
        let panel_top = margin_top + panel_idx as f64 * (plot_height + gap);
        let panel_bottom = panel_top + plot_height;

        // Calculate y range
        let y_min = yfoil_data
            .iter()
            .chain(xfoil_data.iter())
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let y_max = yfoil_data
            .iter()
            .chain(xfoil_data.iter())
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let y_range = (y_max - y_min).max(1e-6);
        let y_padding = 0.1 * y_range;
        let y_min_p = y_min - y_padding;
        let y_max_p = y_max + y_padding;

        let to_svg_y = |y: f64| -> f64 { panel_bottom - (y - y_min_p) / (y_max_p - y_min_p) * plot_height };

        // Plot border
        writeln!(
            file,
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="none" stroke="#000000" stroke-width="1"/>"##,
            margin_left, panel_top, plot_width, plot_height
        )?;

        // Gridlines
        writeln!(file, r#"<g stroke="{}" stroke-width="0.5">"#, grid_color)?;
        for i in 0..=10 {
            let x = i as f64 / 10.0;
            let sx = to_svg_x(x);
            writeln!(
                file,
                r#"<line x1="{:.6}" y1="{:.1}" x2="{:.6}" y2="{:.1}"/>"#,
                sx, panel_top, sx, panel_bottom
            )?;
        }
        let y_step = nice_step(y_max_p - y_min_p, 4);
        let mut y_val = (y_min_p / y_step).ceil() * y_step;
        while y_val <= y_max_p {
            let sy = to_svg_y(y_val);
            if sy > panel_top && sy < panel_bottom {
                writeln!(
                    file,
                    r#"<line x1="{:.1}" y1="{:.6}" x2="{:.1}" y2="{:.6}"/>"#,
                    margin_left,
                    sy,
                    margin_left + plot_width,
                    sy
                )?;
            }
            y_val += y_step;
        }
        writeln!(file, "</g>")?;

        // Y-axis label
        writeln!(
            file,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="11" transform="rotate(-90 {:.1} {:.1})">{}</text>"#,
            margin_left - 55.0,
            panel_top + plot_height / 2.0,
            margin_left - 55.0,
            panel_top + plot_height / 2.0,
            label
        )?;

        // Y tick labels (only show a couple)
        writeln!(file, r#"<g font-family="sans-serif" font-size="9" text-anchor="end">"#)?;
        y_val = (y_min_p / y_step).ceil() * y_step;
        while y_val <= y_max_p {
            let sy = to_svg_y(y_val);
            if sy > panel_top + 5.0 && sy < panel_bottom - 5.0 {
                // Format based on magnitude
                let formatted = if y_val.abs() < 0.001 && y_val.abs() > 1e-10 {
                    format!("{:.2e}", y_val)
                } else if y_val.abs() < 1.0 {
                    format!("{:.4}", y_val)
                } else {
                    format!("{:.2}", y_val)
                };
                writeln!(
                    file,
                    r#"<text x="{:.1}" y="{:.1}">{}</text>"#,
                    margin_left - 5.0,
                    sy + 3.0,
                    formatted
                )?;
            }
            y_val += y_step;
        }
        writeln!(file, "</g>")?;

        // YFoil data
        write!(
            file,
            r#"<polyline fill="none" stroke="{}" stroke-width="1.5" points=""#,
            yfoil_color
        )?;
        for (i, (&x, &y)) in yfoil_bl.x.iter().zip(yfoil_data.iter()).enumerate() {
            if i > 0 {
                write!(file, " ")?;
            }
            write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_y(y))?;
        }
        writeln!(file, r#""/>"#)?;

        // XFOIL data
        write!(
            file,
            r#"<polyline fill="none" stroke="{}" stroke-width="1.5" stroke-dasharray="4,2" points=""#,
            xfoil_color
        )?;
        for (i, (&x, &y)) in xfoil_dump.x.iter().zip(xfoil_data.iter()).enumerate() {
            if i > 0 {
                write!(file, " ")?;
            }
            write!(file, "{:.6},{:.6}", to_svg_x(x), to_svg_y(y))?;
        }
        writeln!(file, r#""/>"#)?;
    }

    // X-axis label at bottom
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="12">x/c</text>"#,
        margin_left + plot_width / 2.0,
        height - 10.0
    )?;

    // Legend at bottom right
    let legend_x = margin_left + plot_width - 140.0;
    let legend_y = height - 35.0;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2"/>"#,
        legend_x,
        legend_y,
        legend_x + 25.0,
        legend_y,
        yfoil_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">YFoil</text>"#,
        legend_x + 30.0,
        legend_y + 4.0
    )?;
    writeln!(
        file,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2" stroke-dasharray="4,2"/>"#,
        legend_x + 70.0,
        legend_y,
        legend_x + 95.0,
        legend_y,
        xfoil_color
    )?;
    writeln!(
        file,
        r#"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="11">XFOIL</text>"#,
        legend_x + 100.0,
        legend_y + 4.0
    )?;

    writeln!(file, "</svg>")?;
    Ok(())
}

// ============================================================================
// Multi-polar plotting (CLI `yfoil plot polar`)
// ============================================================================

/// One polar curve with its legend label
#[derive(Debug, Clone)]
pub struct PolarSeries {
    /// Legend label
    pub label: String,
    /// Angle of attack values (degrees)
    pub alpha: Vec<f64>,
    /// Lift coefficient values
    pub cl: Vec<f64>,
    /// Drag coefficient values
    pub cd: Vec<f64>,
    /// Moment coefficient values (about quarter chord)
    pub cm: Vec<f64>,
}

impl PolarSeries {
    /// Build from a `PolarOutput`, keeping only converged points.
    ///
    /// The label is the polar's `label` if set, otherwise its `airfoil` name.
    pub fn from_polar_output(polar: &crate::output::PolarOutput) -> Self {
        let pts: Vec<_> = polar.results.iter().filter(|p| p.is_converged()).collect();
        Self {
            label: polar.label.clone().unwrap_or_else(|| polar.foil.clone()),
            alpha: pts.iter().map(|p| p.alpha_deg).collect(),
            cl: pts.iter().map(|p| p.cl).collect(),
            cd: pts.iter().map(|p| p.cd.unwrap_or(0.0)).collect(),
            cm: pts.iter().map(|p| p.cm).collect(),
        }
    }
}

/// Point-marker shape, cycled per series. Open markers come first so that coincident series
/// (YFoil over XFOIL in the validation plots) stay visible through each other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Marker {
    /// Diagonal cross (no filled variant)
    Cross,
    /// Plus sign (no filled variant)
    Plus,
    /// Hollow circle
    CircleHollow,
    /// Hollow square
    SquareHollow,
    /// Hollow upward triangle
    TriangleHollow,
    /// Hollow diamond
    DiamondHollow,
    /// Hollow downward triangle
    TriangleDownHollow,
    /// Filled circle
    Circle,
    /// Filled square
    Square,
    /// Filled upward triangle
    Triangle,
    /// Filled diamond
    Diamond,
    /// Filled downward triangle
    TriangleDown,
}

impl Marker {
    /// Default cycle order: every open marker (`x + o □ △ ◇ ▽`), then the filled variants
    pub const CYCLE: [Marker; 12] = [
        Marker::Cross,
        Marker::Plus,
        Marker::CircleHollow,
        Marker::SquareHollow,
        Marker::TriangleHollow,
        Marker::DiamondHollow,
        Marker::TriangleDownHollow,
        Marker::Circle,
        Marker::Square,
        Marker::Triangle,
        Marker::Diamond,
        Marker::TriangleDown,
    ];
}

/// Pixel offsets from a marker centre
type Verts = Vec<(i32, i32)>;

/// Vertex lists describing a marker: an optional filled polygon and up to two stroked paths,
/// in pixel offsets from the marker centre with half-width `s`.
fn marker_shape(marker: Marker, s: i32) -> (Verts, Verts, Verts) {
    let circle = |r: i32| -> Verts {
        (0..=16)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / 16.0;
                ((r as f64 * t.cos()).round() as i32, (r as f64 * t.sin()).round() as i32)
            })
            .collect()
    };
    let square = vec![(-s, -s), (s, -s), (s, s), (-s, s), (-s, -s)];
    let triangle = vec![(0, -s - 1), (s + 1, s), (-s - 1, s), (0, -s - 1)];
    let triangle_down = vec![(0, s + 1), (s + 1, -s), (-s - 1, -s), (0, s + 1)];
    let diamond = vec![(0, -s - 1), (s + 1, 0), (0, s + 1), (-s - 1, 0), (0, -s - 1)];
    let filled = |v: Verts| (v, vec![], vec![]);
    let hollow = |v: Verts| (vec![], v, vec![]);
    match marker {
        Marker::Cross => (vec![], vec![(-s, -s), (s, s)], vec![(-s, s), (s, -s)]),
        Marker::Plus => (vec![], vec![(-s - 1, 0), (s + 1, 0)], vec![(0, -s - 1), (0, s + 1)]),
        Marker::CircleHollow => hollow(circle(s)),
        Marker::SquareHollow => hollow(square),
        Marker::TriangleHollow => hollow(triangle),
        Marker::DiamondHollow => hollow(diamond),
        Marker::TriangleDownHollow => hollow(triangle_down),
        Marker::Circle => filled(circle(s)),
        Marker::Square => filled(square),
        Marker::Triangle => filled(triangle),
        Marker::Diamond => filled(diamond),
        Marker::TriangleDown => filled(triangle_down),
    }
}

/// The one concrete element type used for every marker kind: legend line, filled polygon and two
/// stroked paths composed on an anchor. No boxed elements, so no lifetime coupling to the backend.
type MarkerElement<C, DB> = ComposedElement<
    C,
    DB,
    PathElement<(i32, i32)>,
    ComposedElement<
        (i32, i32),
        DB,
        ComposedElement<(i32, i32), DB, Polygon<(i32, i32)>, PathElement<(i32, i32)>>,
        PathElement<(i32, i32)>,
    >,
>;

/// A marker at `at` in a colour. With `legend_line` a horizontal stroke is drawn through it so
/// the same element serves as the key entry.
fn marker_element<C: Clone, DB: DrawingBackend>(
    marker: Marker,
    at: C,
    color: RGBColor,
    size: i32,
    legend_line: bool,
) -> MarkerElement<C, DB> {
    let (fill, stroke1, stroke2) = marker_shape(marker, size);
    let hollow: ShapeStyle = color.stroke_width(1);
    // An empty path draws nothing, so the line is present only for key entries.
    let line = PathElement::new(
        if legend_line { vec![(-16, 0), (16, 0)] } else { vec![] },
        color.stroke_width(2),
    );
    EmptyElement::<C, DB>::at(at)
        + line
        + Polygon::new(fill, color.filled())
        + PathElement::new(stroke1, hollow)
        + PathElement::new(stroke2, hollow)
}

/// Configuration for the multi-polar plot
#[derive(Debug, Clone)]
pub struct PolarsPlotConfig {
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
    /// Title; defaults to the single series' label, or "Polar comparison"
    pub title: Option<String>,
    /// Background colour (RGB)
    pub background: (u8, u8, u8),
    /// Series colours (RGB), cycled when there are more series than colours
    pub palette: Vec<(u8, u8, u8)>,
    /// Series markers, cycled independently of the palette
    pub markers: Vec<Marker>,
    /// Marker half-size in pixels
    pub marker_size: i32,
}

impl Default for PolarsPlotConfig {
    fn default() -> Self {
        Self {
            width: 1400,
            height: 1000,
            title: None,
            background: (255, 255, 255),
            palette: vec![
                (0, 100, 200),
                (200, 50, 50),
                (30, 150, 60),
                (220, 130, 0),
                (120, 60, 180),
                (0, 150, 160),
            ],
            markers: Marker::CYCLE.to_vec(),
            marker_size: 3,
        }
    }
}

/// Plot one or more polars to SVG: CL–α, CL–CD (drag polar), CM–α and CD–α in a 2×2 grid.
pub fn plot_polars_svg<P: AsRef<Path>>(
    series: &[PolarSeries],
    output_path: P,
    config: &PolarsPlotConfig,
) -> Result<(), PlotError> {
    let root = SVGBackend::new(&output_path, (config.width, config.height)).into_drawing_area();
    plot_polars_impl(&root, series, config)
}

/// Plot one or more polars to PNG (same layout as [`plot_polars_svg`]).
pub fn plot_polars_png<P: AsRef<Path>>(
    series: &[PolarSeries],
    output_path: P,
    config: &PolarsPlotConfig,
) -> Result<(), PlotError> {
    let root = BitMapBackend::new(&output_path, (config.width, config.height)).into_drawing_area();
    plot_polars_impl(&root, series, config)
}

/// Padded [min, max] of a set of values; `floor` guards a degenerate range.
fn padded_range<'a>(values: impl Iterator<Item = &'a f64>, pad: f64, floor: f64) -> std::ops::Range<f64> {
    let (lo, hi) = values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
        (lo.min(v), hi.max(v))
    });
    if !lo.is_finite() || !hi.is_finite() {
        return 0.0..1.0;
    }
    let span = (hi - lo).max(floor);
    (lo - pad * span)..(hi + pad * span)
}

fn plot_polars_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    series: &[PolarSeries],
    config: &PolarsPlotConfig,
) -> Result<(), PlotError>
where
    DB::ErrorType: 'static,
{
    if series.is_empty() {
        return Err(PlotError::Drawing("no polars to plot".to_string()));
    }
    let bg = RGBColor(config.background.0, config.background.1, config.background.2);
    root.fill(&bg).map_err(|e| PlotError::Drawing(e.to_string()))?;

    let title = config.title.clone().unwrap_or_else(|| {
        if series.len() == 1 {
            series[0].label.clone()
        } else {
            "Polar comparison".to_string()
        }
    });
    let (title_area, body) = root.split_vertically(40);
    title_area
        .titled(&title, ("sans-serif", 20))
        .map_err(|e| PlotError::Drawing(e.to_string()))?;

    let panels = body.split_evenly((2, 2));

    // Shared axis ranges across every series so the panels line up
    let alpha_r = padded_range(series.iter().flat_map(|s| s.alpha.iter()), 0.05, 1.0);
    let cl_r = padded_range(series.iter().flat_map(|s| s.cl.iter()), 0.08, 0.1);
    let cd_r = padded_range(series.iter().flat_map(|s| s.cd.iter()), 0.08, 0.001);
    let cm_r = padded_range(series.iter().flat_map(|s| s.cm.iter()), 0.08, 0.01);

    // (x-label, y-label, x-range, y-range, point extractor)
    type Extract = fn(&PolarSeries) -> Vec<(f64, f64)>;
    type PanelSpec<'a> = (&'a str, &'a str, std::ops::Range<f64>, std::ops::Range<f64>, Extract);
    let panel_specs: [PanelSpec; 4] = [
        ("α (deg)", "CL", alpha_r.clone(), cl_r.clone(), |s| {
            s.alpha.iter().zip(s.cl.iter()).map(|(&x, &y)| (x, y)).collect()
        }),
        ("CD", "CL", cd_r.clone(), cl_r.clone(), |s| {
            s.cd.iter().zip(s.cl.iter()).map(|(&x, &y)| (x, y)).collect()
        }),
        ("α (deg)", "CM", alpha_r.clone(), cm_r.clone(), |s| {
            s.alpha.iter().zip(s.cm.iter()).map(|(&x, &y)| (x, y)).collect()
        }),
        ("α (deg)", "CD", alpha_r.clone(), cd_r.clone(), |s| {
            s.alpha.iter().zip(s.cd.iter()).map(|(&x, &y)| (x, y)).collect()
        }),
    ];
    for (panel, (xl, yl, xr, yr, extract)) in panels.iter().zip(panel_specs.iter()) {
        let mut chart = ChartBuilder::on(panel)
            .margin(12)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(xr.clone(), yr.clone())
            .map_err(|e| PlotError::Drawing(e.to_string()))?;
        chart
            .configure_mesh()
            .x_desc(*xl)
            .y_desc(*yl)
            .x_labels(8)
            .y_labels(8)
            .bold_line_style(RGBColor(205, 205, 205))
            .light_line_style(RGBColor(235, 235, 235))
            .x_max_light_lines(5)
            .y_max_light_lines(5)
            .draw()
            .map_err(|e| PlotError::Drawing(e.to_string()))?;

        for (i, s) in series.iter().enumerate() {
            let c = config.palette[i % config.palette.len()];
            let color = RGBColor(c.0, c.1, c.2);
            let marker = config.markers[i % config.markers.len()];
            let size = config.marker_size;
            let pts = extract(s);
            chart
                .draw_series(LineSeries::new(pts.clone(), color.stroke_width(2)))
                .map_err(|e| PlotError::Drawing(e.to_string()))?
                .label(&s.label)
                .legend(move |(x, y)| marker_element(marker, (x + 14, y), color, size, true));
            chart
                .draw_series(pts.into_iter().map(|p| marker_element(marker, p, color, size, false)))
                .map_err(|e| PlotError::Drawing(e.to_string()))?;
        }

        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperLeft)
            .background_style(WHITE.mix(0.8))
            .border_style(BLACK)
            .draw()
            .map_err(|e| PlotError::Drawing(e.to_string()))?;
    }

    root.present().map_err(|e| PlotError::Drawing(e.to_string()))?;
    Ok(())
}
