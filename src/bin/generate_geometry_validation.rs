//! Generate geometry validation documentation
//!
//! This binary generates validation documentation comparing YFoil geometry
//! against XFOIL reference data. It creates:
//! - YFoil geometry files (.json)
//! - Comparison SVG plots
//! - Updated markdown report with metrics
//!
//! Run with: cargo run --bin generate_geometry_validation --features plotting

use std::fs;
use std::path::Path;

use yfoil::geometry::{naca_4digit, repanel_xfoil, Geometry, PaneConfig};
use yfoil::output::{
    parse_xfoil_dat_file, plot_geometry_comparison_svg, GeometryComparisonPlotConfig, PlotError,
};

/// Test case definition
struct TestCase {
    airfoil: &'static str,
    n_panels: usize,
    designation: &'static str,
}

/// Metrics computed for each test case
#[derive(Debug, Clone)]
struct GeometryMetrics {
    n_points: usize,
    x_rms_error: f64,
    x_max_error: f64,
    y_rms_error: f64,
    y_max_error: f64,
    total_rms_error: f64,
    total_max_error: f64,
}

const TEST_CASES: &[TestCase] = &[
    TestCase {
        airfoil: "naca0012",
        n_panels: 81,
        designation: "0012",
    },
    TestCase {
        airfoil: "naca0012",
        n_panels: 160,
        designation: "0012",
    },
    TestCase {
        airfoil: "naca4412",
        n_panels: 81,
        designation: "4412",
    },
    TestCase {
        airfoil: "naca4412",
        n_panels: 160,
        designation: "4412",
    },
];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base_dir = Path::new("docs/validation/assets/geometry");

    println!("YFoil Geometry Validation Generator");
    println!("====================================\n");

    let mut all_metrics = Vec::new();

    for case in TEST_CASES {
        println!(
            "Processing {} with {} panels...",
            case.airfoil.to_uppercase(),
            case.n_panels
        );

        // Generate YFoil geometry
        let yfoil_coords = generate_yfoil_geometry(case)?;
        let case_name = format!("{}-{}", case.airfoil, case.n_panels);

        // Save YFoil geometry to JSON
        let yfoil_path = base_dir.join("yfoil").join(format!("{}.json", case_name));
        save_yfoil_geometry(&yfoil_coords, &yfoil_path)?;
        println!("  Saved YFoil geometry to {}", yfoil_path.display());

        // Load XFOIL reference
        let xfoil_path = base_dir.join("xfoil").join(format!("{}.dat", case_name));
        let xfoil_coords = parse_xfoil_dat_file(&xfoil_path)?;
        println!(
            "  Loaded XFOIL reference ({} points)",
            xfoil_coords.len()
        );

        // Compute metrics
        let metrics = compute_metrics(&xfoil_coords, &yfoil_coords);
        println!(
            "  Metrics: total_max_error = {:.2e}",
            metrics.total_max_error
        );

        // Generate comparison plot
        let plot_path = base_dir.join("plots").join(format!("{}-comparison.svg", case_name));
        generate_comparison_plot(case, &xfoil_coords, &yfoil_coords, &plot_path)?;
        println!("  Generated plot: {}", plot_path.display());

        all_metrics.push((case, metrics));
        println!();
    }

    // Update markdown report
    update_markdown_report(&all_metrics)?;
    println!("Updated docs/validation/geometry/README.md");

    // Print summary
    println!("\n=== Summary ===");
    for (case, metrics) in &all_metrics {
        let case_name = format!("{}-{}", case.airfoil, case.n_panels);
        println!(
            "  {}: max error = {:.2e}",
            case_name,
            metrics.total_max_error
        );
    }

    Ok(())
}

/// Generate YFoil geometry for a test case
fn generate_yfoil_geometry(case: &TestCase) -> Result<Vec<(f64, f64)>, Box<dyn std::error::Error>> {
    // Generate base NACA geometry with extra points
    let base_geometry = naca_4digit(case.designation, 300)?;

    // Repanel using XFOIL's algorithm
    let config = PaneConfig::default();
    let repaneled = repanel_xfoil(&base_geometry, case.n_panels, &config);

    // Convert to coordinate pairs
    let coords: Vec<(f64, f64)> = repaneled
        .x_c
        .iter()
        .zip(repaneled.y_c.iter())
        .map(|(&x, &y)| (x, y))
        .collect();

    Ok(coords)
}

/// Save YFoil geometry to JSON file
fn save_yfoil_geometry(
    coords: &[(f64, f64)],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let geometry = Geometry {
        reference: [0.25, 0.0],
        x_c: coords.iter().map(|(x, _)| *x).collect(),
        y_c: coords.iter().map(|(_, y)| *y).collect(),
    };

    let json = serde_json::to_string_pretty(&geometry)?;
    fs::write(path, json)?;
    Ok(())
}

/// Compute comparison metrics between two sets of coordinates
fn compute_metrics(xfoil: &[(f64, f64)], yfoil: &[(f64, f64)]) -> GeometryMetrics {
    let n = xfoil.len().min(yfoil.len());

    let mut x_sum_sq: f64 = 0.0;
    let mut x_max: f64 = 0.0;
    let mut y_sum_sq: f64 = 0.0;
    let mut y_max: f64 = 0.0;
    let mut total_sum_sq: f64 = 0.0;
    let mut total_max: f64 = 0.0;

    for i in 0..n {
        let dx = (xfoil[i].0 - yfoil[i].0).abs();
        let dy = (xfoil[i].1 - yfoil[i].1).abs();
        let dist = (dx * dx + dy * dy).sqrt();

        x_sum_sq += dx * dx;
        x_max = x_max.max(dx);
        y_sum_sq += dy * dy;
        y_max = y_max.max(dy);
        total_sum_sq += dist * dist;
        total_max = total_max.max(dist);
    }

    let n_f = n as f64;
    GeometryMetrics {
        n_points: n,
        x_rms_error: (x_sum_sq / n_f).sqrt(),
        x_max_error: x_max,
        y_rms_error: (y_sum_sq / n_f).sqrt(),
        y_max_error: y_max,
        total_rms_error: (total_sum_sq / n_f).sqrt(),
        total_max_error: total_max,
    }
}

/// Generate comparison plot
fn generate_comparison_plot(
    case: &TestCase,
    xfoil: &[(f64, f64)],
    yfoil: &[(f64, f64)],
    output_path: &Path,
) -> Result<(), PlotError> {
    let config = GeometryComparisonPlotConfig {
        title: Some(format!(
            "{} - {} Panels",
            case.airfoil.to_uppercase(),
            case.n_panels
        )),
        ..Default::default()
    };

    plot_geometry_comparison_svg(xfoil, yfoil, output_path, &config)
}

/// Update the markdown report
fn update_markdown_report(
    metrics: &[(&TestCase, GeometryMetrics)],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut report = String::new();

    report.push_str("# Geometry Validation\n\n");
    report.push_str("This document validates YFoil's geometry generation against XFOIL reference data.\n\n");

    // Methodology
    report.push_str("## Methodology\n\n");
    report.push_str("The validation compares panel coordinates generated by YFoil's implementation of XFOIL's\n");
    report.push_str("PANE algorithm against coordinates generated directly by XFOIL. Both use the same\n");
    report.push_str("algorithm parameters:\n\n");
    report.push_str("- Panel bunching parameter: 1.0\n");
    report.push_str("- TE/LE panel density ratio: 0.15\n");
    report.push_str("- Refined area/LE panel density ratio: 0.2\n\n");

    // Summary table
    report.push_str("## Summary\n\n");
    report.push_str("| Airfoil | Panels | Total RMS Error | Total Max Error |\n");
    report.push_str("|---------|--------|-----------------|----------------|\n");

    for (case, m) in metrics {
        report.push_str(&format!(
            "| {} | {} | {:.2e} | {:.2e} |\n",
            case.airfoil.to_uppercase(),
            case.n_panels,
            m.total_rms_error,
            m.total_max_error
        ));
    }
    report.push_str("\n");

    // Individual case sections
    for (case, m) in metrics {
        let case_name = format!("{}-{}", case.airfoil, case.n_panels);

        report.push_str(&format!(
            "## {} ({} panels)\n\n",
            case.airfoil.to_uppercase(),
            case.n_panels
        ));

        // Commands used
        report.push_str("### Generation Commands\n\n");
        report.push_str("**XFOIL script:**\n");
        report.push_str(&format!(
            "```\nPLOP\nG F\n\nNACA {}\nPANE\nPPAR\nN {}\n\n\nPSAVE {}.dat\n\nQUIT\n```\n\n",
            case.designation, case.n_panels, case_name
        ));
        report.push_str("**YFoil command (equivalent):**\n");
        report.push_str(&format!(
            "```\nyfoil geom naca {} --paneler xfoil -n {} -o {}.json\n```\n\n",
            case.designation, case.n_panels, case_name
        ));

        // Data files
        report.push_str("### Data Files\n\n");
        report.push_str(&format!(
            "- XFOIL reference: [`../assets/geometry/xfoil/{}.dat`](../assets/geometry/xfoil/{}.dat)\n",
            case_name, case_name
        ));
        report.push_str(&format!(
            "- YFoil generated: [`../assets/geometry/yfoil/{}.json`](../assets/geometry/yfoil/{}.json)\n\n",
            case_name, case_name
        ));

        // Metrics table
        report.push_str("### Metrics\n\n");
        report.push_str("| Metric | Value |\n");
        report.push_str("|--------|-------|\n");
        report.push_str(&format!("| Number of points | {} |\n", m.n_points));
        report.push_str(&format!("| X RMS error | {:.2e} |\n", m.x_rms_error));
        report.push_str(&format!("| X max error | {:.2e} |\n", m.x_max_error));
        report.push_str(&format!("| Y RMS error | {:.2e} |\n", m.y_rms_error));
        report.push_str(&format!("| Y max error | {:.2e} |\n", m.y_max_error));
        report.push_str(&format!("| Total RMS error | {:.2e} |\n", m.total_rms_error));
        report.push_str(&format!("| Total max error | {:.2e} |\n", m.total_max_error));
        report.push_str("\n");

        // Comparison plot
        report.push_str("### Comparison Plot\n\n");
        report.push_str(&format!(
            "![{} comparison](../assets/geometry/plots/{}-comparison.svg)\n\n",
            case_name, case_name
        ));
    }

    // Regeneration instructions
    report.push_str("## Regeneration\n\n");
    report.push_str("To regenerate this validation after code changes:\n\n");
    report.push_str("```bash\ncargo run --bin generate_geometry_validation --features plotting\n```\n\n");
    report.push_str("This will:\n");
    report.push_str("1. Regenerate YFoil geometry for all test cases\n");
    report.push_str("2. Recompute comparison metrics against XFOIL baseline\n");
    report.push_str("3. Regenerate comparison SVG plots\n");
    report.push_str("4. Update this markdown report\n\n");
    report.push_str("**Note:** XFOIL reference data is NOT regenerated (it's the committed baseline).\n");

    fs::write("docs/validation/geometry/README.md", report)?;
    Ok(())
}
