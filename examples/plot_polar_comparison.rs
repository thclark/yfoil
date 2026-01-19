//! Generate polar comparison plot from existing data files
//!
//! Run with: cargo run --example plot_polar_comparison --features plotting

use yfoil::output::{
    parse_xfoil_polar_file, parse_yfoil_polar_file, plot_polar_comparison_svg,
    stitch_polars, PolarPlotConfig,
};

fn main() {
    println!("=== Generating Polar Comparison Plot ===\n");

    // Parse YFoil polars
    println!("Loading YFoil polars...");
    let yfoil_pos = parse_yfoil_polar_file(".tmp/yfoil_polar_pos.json")
        .expect("Failed to parse yfoil_polar_pos.json");
    let yfoil_neg = parse_yfoil_polar_file(".tmp/yfoil_polar_neg.json")
        .expect("Failed to parse yfoil_polar_neg.json");

    // Parse XFOIL polars
    println!("Loading XFOIL polars...");
    let xfoil_pos = parse_xfoil_polar_file(".tmp/xfoil_polar_pos.txt")
        .expect("Failed to parse xfoil_polar_pos.txt");
    let xfoil_neg = parse_xfoil_polar_file(".tmp/xfoil_polar_neg.txt")
        .expect("Failed to parse xfoil_polar_neg.txt");

    // Stitch polars
    println!("Stitching polars...");
    let yfoil_full = stitch_polars(&yfoil_pos, &yfoil_neg);
    let xfoil_full = stitch_polars(&xfoil_pos, &xfoil_neg);

    println!("  YFoil: {} points from alpha={:.1}° to {:.1}°",
        yfoil_full.alpha.len(),
        yfoil_full.alpha.first().unwrap_or(&0.0),
        yfoil_full.alpha.last().unwrap_or(&0.0));
    println!("  XFOIL: {} points from alpha={:.1}° to {:.1}°",
        xfoil_full.alpha.len(),
        xfoil_full.alpha.first().unwrap_or(&0.0),
        xfoil_full.alpha.last().unwrap_or(&0.0));

    // Generate comparison plot
    let config = PolarPlotConfig {
        title: Some("NACA 0012 Polar Comparison: YFoil vs XFOIL (Re=1M)".to_string()),
        ..PolarPlotConfig::default()
    };

    let output_path = ".tmp/polar_comparison.svg";
    println!("\nGenerating comparison plot to {}...", output_path);

    plot_polar_comparison_svg(&yfoil_full, &xfoil_full, output_path, &config)
        .expect("Failed to generate comparison plot");

    println!("Done! Open {} to view the comparison.", output_path);
}
