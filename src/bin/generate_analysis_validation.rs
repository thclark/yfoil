//! Generate analysis validation documentation
//!
//! This example generates comprehensive validation documentation comparing
//! YFoil and XFOIL analysis results for NACA 0012 and NACA 4412 airfoils.
//!
//! IMPORTANT: Both XFOIL and YFoil run proper polar sweeps, extracting full
//! BL distributions at validation angles (0°, 5°, 10°, 15°, -5°, -10°, -15°).
//! This ensures solutions are properly initialized as if part of a complete
//! polar sweep (e.g., results at 5° are initialized from 4° → 3° → 2° → 1° → 0°).
//!
//! Usage:
//!   cargo run --bin generate_analysis_validation --features plotting
//!
//! To also run XFOIL (generates baseline data):
//!   cargo run --bin generate_analysis_validation --features plotting -- --run-xfoil

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::output::{
    parse_xfoil_cp_file, parse_xfoil_dump_file, parse_xfoil_polar_file_with_cm,
    plot_bl_comparison_svg, plot_cp_ue_comparison_svg, plot_polar_3panel_svg, stitch_polars_with_cm,
    PolarDataWithCm, PolarPlotConfig, YfoilBLDist,
};
use yfoil::solver::{solve_viscous_with_init, ViscalConfig};

const BASE_DIR: &str = "docs/validation/assets/analysis";
const XFOIL_BIN: &str = "xfoil/xfoil6.99/bin/xfoil";
const REYNOLDS: f64 = 1_000_000.0;
const MACH: f64 = 0.0;
const NCRIT: f64 = 9.0;

/// Validation angles where full BL distributions are compared
const VALIDATION_ANGLES: [i32; 7] = [-15, -10, -5, 0, 5, 10, 15];

/// Airfoil configuration
struct AirfoilConfig {
    name: &'static str,
    geometry_file: &'static str,
    title: &'static str,
}

const AIRFOILS: [AirfoilConfig; 2] = [
    AirfoilConfig {
        name: "naca0012",
        geometry_file: "geometry/naca0012.json",
        title: "NACA 0012",
    },
    AirfoilConfig {
        name: "naca4412",
        geometry_file: "geometry/naca4412.json",
        title: "NACA 4412",
    },
];

/// YFoil analysis result for plotting
struct YfoilAnalysisResult {
    x: Vec<f64>,
    s: Vec<f64>,
    cp: Vec<f64>,
    ue: Vec<f64>,
    theta: Vec<f64>,
    dstar: Vec<f64>,
    h: Vec<f64>,
    hs: Vec<f64>,
    cf: Vec<f64>,
}

/// Results from YFoil sweep including polar and BL distributions at validation angles
struct YfoilSweepResults {
    polar: PolarDataWithCm,
    bl_distributions: HashMap<i32, YfoilAnalysisResult>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let run_xfoil = args.iter().any(|a| a == "--run-xfoil");

    println!("=== Analysis Validation Generator ===\n");

    // Step 1: Run XFOIL if requested
    if run_xfoil {
        println!("Running XFOIL to generate baseline data...");
        run_xfoil_scripts()?;
    } else {
        println!("Skipping XFOIL (use --run-xfoil to generate baseline data)");
    }

    // Step 2: Generate YFoil data and comparison plots for each airfoil
    for airfoil in &AIRFOILS {
        println!("\nProcessing {}...", airfoil.title);

        // Check if XFOIL data exists
        let xfoil_polar_pos = format!("{}/xfoil/{}/polar_pos.txt", BASE_DIR, airfoil.name);
        if !Path::new(&xfoil_polar_pos).exists() {
            println!(
                "  Warning: XFOIL data not found for {}. Run with --run-xfoil first.",
                airfoil.name
            );
            continue;
        }

        // Run YFoil sweep - this collects both polar data and BL at validation angles
        println!("  Running YFoil polar sweep (collecting BL at validation angles)...");
        let yfoil_results = run_yfoil_sweep(airfoil)?;

        // Generate polar comparison
        generate_polar_comparison(airfoil, &yfoil_results.polar)?;

        // Generate BL/Cp comparisons at each validation angle
        for &alpha in &VALIDATION_ANGLES {
            if let Some(bl_result) = yfoil_results.bl_distributions.get(&alpha) {
                generate_angle_comparison(airfoil, alpha, bl_result)?;
            } else {
                println!("    Skipping α = {}°: YFoil did not converge", alpha);
            }
        }
    }

    // Step 3: Generate markdown reports
    println!("\nGenerating documentation...");
    generate_readme()?;
    for airfoil in &AIRFOILS {
        generate_airfoil_report(airfoil)?;
    }

    println!("\nValidation documentation generated successfully!");
    println!("  Main report: docs/validation/analysis/README.md");
    println!("  NACA 0012:   docs/validation/analysis/naca0012.md");
    println!("  NACA 4412:   docs/validation/analysis/naca4412.md");

    Ok(())
}

/// Run XFOIL sweep scripts to generate baseline data
///
/// Only runs scripts named *_sweep.xfoil which perform proper polar sweeps
/// and dump BL data at validation angles.
fn run_xfoil_scripts() -> Result<(), Box<dyn std::error::Error>> {
    let cwd = std::env::current_dir()?;
    let scripts_dir = cwd.join(BASE_DIR).join("scripts");
    let xfoil_bin = cwd.join(XFOIL_BIN);

    if !xfoil_bin.exists() {
        return Err(format!("XFOIL binary not found at {:?}", xfoil_bin).into());
    }

    // Find only *_sweep.xfoil scripts (these do proper polar sweeps with BL dumps)
    let mut scripts: Vec<_> = fs::read_dir(&scripts_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.path().file_name().unwrap().to_string_lossy().to_string();
            name.ends_with("_sweep.xfoil")
        })
        .collect();
    scripts.sort_by_key(|e| e.path());

    for entry in scripts {
        let path = entry.path();
        let script_name = path.file_name().unwrap().to_string_lossy();
        println!("  Running {}...", script_name);

        // Read script content
        let script = fs::read_to_string(&path)?;

        // Run XFOIL with script as stdin and wait for completion
        let mut child = Command::new(&xfoil_bin)
            .current_dir(&scripts_dir)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()?;

        // Write script to stdin
        if let Some(mut stdin) = child.stdin.take() {
            stdin.write_all(script.as_bytes())?;
        }

        // Wait for XFOIL to complete
        let status = child.wait()?;
        if !status.success() {
            eprintln!("    Warning: XFOIL exited with status: {}", status);
        }
    }

    Ok(())
}

/// Run YFoil polar sweep and collect BL distributions at validation angles
///
/// This function performs a proper polar sweep:
/// - Positive sweep: 0° → 1° → 2° → ... → 15° (collecting BL at 0°, 5°, 10°, 15°)
/// - Reinitialize
/// - Negative sweep: -1° → -2° → ... → -15° (collecting BL at -5°, -10°, -15°)
///
/// Each angle's solution is initialized from the previous converged solution,
/// ensuring identical initialization to XFOIL's sweep behavior.
fn run_yfoil_sweep(airfoil: &AirfoilConfig) -> Result<YfoilSweepResults, Box<dyn std::error::Error>>
{
    let geom_path = format!("{}/{}", BASE_DIR, airfoil.geometry_file);
    let geometry = read_geometry_from_file(&geom_path)?;
    let paneled = create_paneled_airfoil(&geometry);

    let conditions = FlowConditions::new(REYNOLDS, MACH, NCRIT, paneled.chord);
    let config = ViscalConfig {
        max_iter: 400,
        ..Default::default()
    };

    let mut alphas = Vec::new();
    let mut cls = Vec::new();
    let mut cds = Vec::new();
    let mut cms = Vec::new();
    let mut bl_distributions: HashMap<i32, YfoilAnalysisResult> = HashMap::new();

    // Positive sweep: 0° to 15°
    let mut prev_dq: Option<Vec<f64>> = None;
    for alpha_deg in 0..=15 {
        let alpha_rad = (alpha_deg as f64).to_radians();
        let result =
            solve_viscous_with_init(&paneled, alpha_rad, &conditions, &config, prev_dq.as_deref());
        if result.converged {
            alphas.push(alpha_deg as f64);
            cls.push(result.cl);
            cds.push(result.cd);
            cms.push(result.cm);
            prev_dq = Some(result.dq_source.clone());

            // Collect BL distribution at validation angles
            if VALIDATION_ANGLES.contains(&alpha_deg) {
                bl_distributions.insert(alpha_deg, extract_bl_result(&result));
            }
        } else {
            println!(
                "    Warning: positive sweep stopped at α = {}° (non-convergence)",
                alpha_deg
            );
            break;
        }
    }

    // Negative sweep (reinitialize): -1° to -15°
    prev_dq = None;
    for alpha_deg in (-15..0).rev() {
        let alpha_rad = (alpha_deg as f64).to_radians();
        let result =
            solve_viscous_with_init(&paneled, alpha_rad, &conditions, &config, prev_dq.as_deref());
        if result.converged {
            alphas.insert(0, alpha_deg as f64);
            cls.insert(0, result.cl);
            cds.insert(0, result.cd);
            cms.insert(0, result.cm);
            prev_dq = Some(result.dq_source.clone());

            // Collect BL distribution at validation angles
            if VALIDATION_ANGLES.contains(&alpha_deg) {
                bl_distributions.insert(alpha_deg, extract_bl_result(&result));
            }
        } else {
            println!(
                "    Warning: negative sweep stopped at α = {}° (non-convergence)",
                alpha_deg
            );
            break;
        }
    }

    Ok(YfoilSweepResults {
        polar: PolarDataWithCm::new(alphas, cls, cds, cms),
        bl_distributions,
    })
}

/// Extract BL result from a viscous solution for plotting
fn extract_bl_result(result: &yfoil::solver::ViscousResult) -> YfoilAnalysisResult {
    let mut x = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();
    let mut theta = Vec::new();
    let mut dstar = Vec::new();
    let mut h = Vec::new();
    let mut hs = Vec::new();
    let mut cf = Vec::new();

    // Upper surface
    for (i, res) in result.bl.upper.iter().enumerate() {
        x.push(res.x);
        s.push(result.bl.s_upper[i]);
        ue.push(res.ue);
        theta.push(res.theta);
        dstar.push(res.dstar);
        h.push(res.h);
        hs.push(res.hs);
        cf.push(res.cf);
    }

    // Lower surface
    for (i, res) in result.bl.lower.iter().enumerate() {
        x.push(res.x);
        s.push(result.bl.s_lower[i]);
        ue.push(res.ue);
        theta.push(res.theta);
        dstar.push(res.dstar);
        h.push(res.h);
        hs.push(res.hs);
        cf.push(res.cf);
    }

    // Wake
    for (i, res) in result.bl.wake.iter().enumerate() {
        x.push(result.bl.x_wake[i]);
        s.push(result.bl.s_wake[i]);
        ue.push(res.ue);
        theta.push(res.theta);
        dstar.push(res.dstar);
        h.push(res.h);
        hs.push(res.hs);
        cf.push(res.cf);
    }

    // Cp from Ue: Cp = 1 - Ue^2
    let cp = ue.iter().map(|&u| 1.0 - u * u).collect();

    YfoilAnalysisResult {
        x,
        s,
        cp,
        ue,
        theta,
        dstar,
        h,
        hs,
        cf,
    }
}

/// Generate polar comparison for an airfoil
fn generate_polar_comparison(
    airfoil: &AirfoilConfig,
    yfoil_polar: &PolarDataWithCm,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("  Generating polar comparison...");

    // Parse XFOIL polars
    let xfoil_pos_path = format!("{}/xfoil/{}/polar_pos.txt", BASE_DIR, airfoil.name);
    let xfoil_neg_path = format!("{}/xfoil/{}/polar_neg.txt", BASE_DIR, airfoil.name);

    if !Path::new(&xfoil_pos_path).exists() {
        println!("    Skipping: XFOIL polar data not found");
        return Ok(());
    }

    let xfoil_pos = parse_xfoil_polar_file_with_cm(&xfoil_pos_path)?;
    let xfoil_neg = if Path::new(&xfoil_neg_path).exists() {
        parse_xfoil_polar_file_with_cm(&xfoil_neg_path)?
    } else {
        PolarDataWithCm::new(vec![], vec![], vec![], vec![])
    };
    let xfoil_polar = stitch_polars_with_cm(&xfoil_pos, &xfoil_neg);

    // Generate plot
    let plot_path = format!("{}/plots/{}/polar.svg", BASE_DIR, airfoil.name);
    let mut config = PolarPlotConfig::default();
    config.title = Some(format!("{} Polar Comparison", airfoil.title));
    config.width = 1600;
    config.height = 500;
    plot_polar_3panel_svg(yfoil_polar, &xfoil_polar, &plot_path, &config)?;

    println!("    Saved polar plot to {}", plot_path);

    Ok(())
}

/// Generate comparison at a specific angle using pre-computed BL result
fn generate_angle_comparison(
    airfoil: &AirfoilConfig,
    alpha: i32,
    yfoil_result: &YfoilAnalysisResult,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("  Generating comparison at α = {}°...", alpha);

    // Check for XFOIL data
    let xfoil_dump_path = format!("{}/xfoil/{}/bl_dump_a{}.txt", BASE_DIR, airfoil.name, alpha);
    let xfoil_cp_path = format!("{}/xfoil/{}/cp_a{}.txt", BASE_DIR, airfoil.name, alpha);

    if !Path::new(&xfoil_dump_path).exists() || !Path::new(&xfoil_cp_path).exists() {
        println!("    Skipping: XFOIL data not found for α = {}°", alpha);
        return Ok(());
    }

    // Parse XFOIL data
    let xfoil_dump = parse_xfoil_dump_file(&xfoil_dump_path)?;
    let xfoil_cp = parse_xfoil_cp_file(&xfoil_cp_path)?;

    // Generate Cp/Ue comparison plot
    let cp_ue_path = format!(
        "{}/plots/{}/cp_ue_a{}.svg",
        BASE_DIR, airfoil.name, alpha
    );
    plot_cp_ue_comparison_svg(
        &yfoil_result.x,
        &yfoil_result.cp,
        &yfoil_result.ue,
        &xfoil_cp,
        &xfoil_dump,
        alpha as f64,
        airfoil.title,
        &cp_ue_path,
    )?;

    // Generate BL comparison plot
    let bl_path = format!("{}/plots/{}/bl_a{}.svg", BASE_DIR, airfoil.name, alpha);
    let yfoil_bl = YfoilBLDist {
        x: yfoil_result.x.clone(),
        s: yfoil_result.s.clone(),
        theta: yfoil_result.theta.clone(),
        dstar: yfoil_result.dstar.clone(),
        h: yfoil_result.h.clone(),
        hs: yfoil_result.hs.clone(),
        cf: yfoil_result.cf.clone(),
        ue: yfoil_result.ue.clone(),
    };
    plot_bl_comparison_svg(&yfoil_bl, &xfoil_dump, alpha as f64, airfoil.title, &bl_path)?;

    Ok(())
}

/// Generate main README
fn generate_readme() -> Result<(), Box<dyn std::error::Error>> {
    let readme_path = "docs/validation/analysis/README.md";
    let mut file = fs::File::create(readme_path)?;

    writeln!(file, "# Analysis Validation")?;
    writeln!(file)?;
    writeln!(
        file,
        "This section validates YFoil's viscous analysis against XFOIL reference data."
    )?;
    writeln!(file)?;
    writeln!(file, "## Test Configuration")?;
    writeln!(file)?;
    writeln!(file, "| Parameter | Value |")?;
    writeln!(file, "|-----------|-------|")?;
    writeln!(file, "| Reynolds Number | 1,000,000 |")?;
    writeln!(file, "| Mach Number | 0.0 |")?;
    writeln!(file, "| Ncrit | 9.0 |")?;
    writeln!(file, "| Max Iterations | 400 |")?;
    writeln!(file, "| Panel Count | 160 |")?;
    writeln!(file)?;
    writeln!(file, "## Validation Methodology")?;
    writeln!(file)?;
    writeln!(
        file,
        "Both XFOIL and YFoil run complete polar sweeps to generate validation data:"
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "1. **Positive sweep**: 0° → 1° → 2° → ... → 15° (each angle initialized from previous)"
    )?;
    writeln!(file, "2. **Reinitialize** at α = -1°")?;
    writeln!(
        file,
        "3. **Negative sweep**: -1° → -2° → ... → -15° (each angle initialized from previous)"
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "BL distributions are extracted at validation angles: 0°, ±5°, ±10°, ±15°"
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "This ensures solutions at each validation angle are properly initialized"
    )?;
    writeln!(
        file,
        "(e.g., the solution at 5° is initialized from 4° → 3° → 2° → 1° → 0°)."
    )?;
    writeln!(file)?;
    writeln!(file, "## Generation Scripts")?;
    writeln!(file)?;
    writeln!(file, "### XFOIL Sweep Script")?;
    writeln!(file)?;
    writeln!(
        file,
        "XFOIL sweep scripts perform a polar sweep, dumping BL distributions at validation angles."
    )?;
    writeln!(
        file,
        "Scripts: [`naca0012_sweep.xfoil`](../assets/analysis/scripts/naca0012_sweep.xfoil), "
    )?;
    writeln!(
        file,
        "[`naca4412_sweep.xfoil`](../assets/analysis/scripts/naca4412_sweep.xfoil)"
    )?;
    writeln!(file)?;
    writeln!(file, "Example (abbreviated NACA 0012 positive sweep):")?;
    writeln!(file, "```")?;
    writeln!(file, "PLOP")?;
    writeln!(file, "G F")?;
    writeln!(file)?;
    writeln!(file, "LOAD ../geometry/naca0012.dat")?;
    writeln!(file, "PCOP")?;
    writeln!(file, "OPER")?;
    writeln!(file, "VISC 1000000")?;
    writeln!(file, "ITER 400")?;
    writeln!(file, "PACC")?;
    writeln!(file, "../xfoil/naca0012/polar_pos.txt")?;
    writeln!(file)?;
    writeln!(file, "ALFA 0")?;
    writeln!(file, "DUMP ../xfoil/naca0012/bl_dump_a0.txt")?;
    writeln!(file, "CPWR ../xfoil/naca0012/cp_a0.txt")?;
    writeln!(file)?;
    writeln!(file, "ALFA 1")?;
    writeln!(file, "ALFA 2")?;
    writeln!(file, "ALFA 3")?;
    writeln!(file, "ALFA 4")?;
    writeln!(file, "ALFA 5")?;
    writeln!(file, "DUMP ../xfoil/naca0012/bl_dump_a5.txt")?;
    writeln!(file, "CPWR ../xfoil/naca0012/cp_a5.txt")?;
    writeln!(file, "...")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(
        file,
        "**Key point**: Each ALFA command initializes from the previous converged solution."
    )?;
    writeln!(
        file,
        "BL dumps are extracted only at validation angles (0°, 5°, 10°, 15°, -5°, -10°, -15°),"
    )?;
    writeln!(
        file,
        "but all intermediate angles are computed to ensure proper initialization."
    )?;
    writeln!(file)?;
    writeln!(file, "### YFoil Equivalent")?;
    writeln!(file)?;
    writeln!(
        file,
        "YFoil performs the same sweep internally in `run_yfoil_sweep()`, collecting"
    )?;
    writeln!(
        file,
        "BL distributions at validation angles while sweeping through all intermediate angles:"
    )?;
    writeln!(file)?;
    writeln!(file, "```rust")?;
    writeln!(file, "// Positive sweep: 0° to 15°")?;
    writeln!(file, "for alpha_deg in 0..=15 {{")?;
    writeln!(file, "    let result = solve_viscous_with_init(&paneled, alpha_rad, &conditions, &config, prev_dq.as_deref());")?;
    writeln!(file, "    if VALIDATION_ANGLES.contains(&alpha_deg) {{")?;
    writeln!(file, "        bl_distributions.insert(alpha_deg, extract_bl_result(&result));")?;
    writeln!(file, "    }}")?;
    writeln!(file, "    prev_dq = Some(result.dq_source.clone());")?;
    writeln!(file, "}}")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(file, "## Test Cases")?;
    writeln!(file)?;
    writeln!(file, "### [NACA 0012](naca0012.md)")?;
    writeln!(file)?;
    writeln!(
        file,
        "Symmetric airfoil test case. Tests basic viscous solution at multiple angles of attack."
    )?;
    writeln!(file)?;
    writeln!(file, "### [NACA 4412](naca4412.md)")?;
    writeln!(file)?;
    writeln!(
        file,
        "Cambered airfoil test case. Tests asymmetric flow behavior and higher-lift conditions."
    )?;
    writeln!(file)?;
    writeln!(file, "## Regenerating Validation")?;
    writeln!(file)?;
    writeln!(file, "```bash")?;
    writeln!(
        file,
        "# Generate XFOIL baseline data (only needed once)"
    )?;
    writeln!(
        file,
        "cargo run --bin generate_analysis_validation --features plotting -- --run-xfoil"
    )?;
    writeln!(file)?;
    writeln!(file, "# Regenerate plots and reports")?;
    writeln!(
        file,
        "cargo run --bin generate_analysis_validation --features plotting"
    )?;
    writeln!(file, "```")?;

    Ok(())
}

/// Generate airfoil-specific report
fn generate_airfoil_report(airfoil: &AirfoilConfig) -> Result<(), Box<dyn std::error::Error>> {
    let report_path = format!("docs/validation/analysis/{}.md", airfoil.name);
    let mut file = fs::File::create(&report_path)?;

    writeln!(file, "# {} Validation", airfoil.title)?;
    writeln!(file)?;
    writeln!(file, "## Polar Comparison")?;
    writeln!(file)?;
    writeln!(
        file,
        "![Polar Comparison](../assets/analysis/plots/{}/polar.svg)",
        airfoil.name
    )?;
    writeln!(file)?;
    writeln!(file, "## Distribution Comparisons by Angle of Attack")?;
    writeln!(file)?;

    for &alpha in &VALIDATION_ANGLES {
        writeln!(file, "### α = {}°", alpha)?;
        writeln!(file)?;
        writeln!(file, "#### Cp and Ue Distributions")?;
        writeln!(file)?;
        writeln!(
            file,
            "![Cp/Ue at α={}°](../assets/analysis/plots/{}/cp_ue_a{}.svg)",
            alpha, airfoil.name, alpha
        )?;
        writeln!(file)?;
        writeln!(file, "#### Boundary Layer Variables")?;
        writeln!(file)?;
        writeln!(
            file,
            "![BL at α={}°](../assets/analysis/plots/{}/bl_a{}.svg)",
            alpha, airfoil.name, alpha
        )?;
        writeln!(file)?;
    }

    Ok(())
}
