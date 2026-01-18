//! YFoil CLI - Aerofoil analysis tool

use std::path::PathBuf;

use clap::{Parser, Subcommand};

use yfoil::bl::FlowConditions;
use yfoil::forces::{calculate_cp, integrate_forces};
use yfoil::geometry::{
    create_paneled_airfoil, naca_4digit, naca_5digit, read_dat_file, read_geometry_from_file,
    repanel_cosine, repanel_xfoil, write_dat_file, write_geometry_to_json, Geometry, PaneConfig,
};
use yfoil::panel::solve_inviscid;
use yfoil::output::PolarOutput;
use yfoil::solver::{compute_polar, solve_viscous, PolarConfig, ViscalConfig};

#[cfg(feature = "plotting")]
use yfoil::output::{plot_paneled_svg, plot_paneled_png, GeometryPlotConfig};

/// YFoil - Rust-based aerofoil analysis tool
#[derive(Parser, Debug)]
#[command(name = "yfoil")]
#[command(version, about = "Rust-based aerofoil analysis tool", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Geometry operations (convert, generate, repanel)
    Geom {
        #[command(subcommand)]
        action: GeomAction,
    },

    /// Analyze airfoil at a single operating point (placeholder)
    Analyze {
        /// Path to geometry file (JSON)
        file: PathBuf,

        /// Angle of attack in degrees
        #[arg(short, long, default_value_t = 0.0)]
        alpha: f64,

        /// Reynolds number
        #[arg(short, long, default_value_t = 1_000_000.0)]
        reynolds: f64,

        /// Mach number
        #[arg(short, long, default_value_t = 0.0)]
        mach: f64,

        /// Critical amplification factor for transition
        #[arg(short, long, default_value_t = 9.0)]
        ncrit: f64,

        /// Inviscid analysis only
        #[arg(long)]
        inviscid: bool,
    },

    /// Generate polar sweep
    Polar {
        /// Path to geometry file (JSON)
        file: PathBuf,

        /// Maximum angle of attack
        #[arg(long, default_value_t = 15.0)]
        alpha_max: f64,

        /// Minimum angle of attack
        #[arg(long, default_value_t = -5.0)]
        alpha_min: f64,

        /// Alpha step size
        #[arg(long, default_value_t = 0.5)]
        alpha_step: f64,

        /// Reynolds number
        #[arg(short, long, default_value_t = 1_000_000.0)]
        reynolds: f64,

        /// Mach number
        #[arg(short, long, default_value_t = 0.0)]
        mach: f64,

        /// Critical amplification factor for transition
        #[arg(short, long, default_value_t = 9.0)]
        ncrit: f64,

        /// Output JSON format
        #[arg(long)]
        json: bool,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

#[derive(Subcommand, Debug)]
enum GeomAction {
    /// Convert between geometry formats
    Convert {
        /// Input file path
        input: PathBuf,

        /// Output format: json, dat
        #[arg(long, default_value = "json")]
        to: String,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Airfoil name (for .dat output)
        #[arg(long, default_value = "Airfoil")]
        name: String,
    },

    /// Generate NACA airfoil
    Naca {
        /// NACA designation (e.g., "0012", "4412", "23015")
        spec: String,

        /// Number of panels
        #[arg(short = 'n', long, default_value_t = 160)]
        panels: usize,

        /// Output format: json, dat
        #[arg(long, default_value = "json")]
        to: String,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Repanel airfoil with new point distribution
    Repanel {
        /// Input file path
        input: PathBuf,

        /// Target number of panels
        #[arg(short = 'n', long, default_value_t = 160)]
        panels: usize,

        /// Repaneling method: xfoil (curvature-based PANE) or cosine (modified cosine spacing)
        #[arg(long, default_value = "xfoil")]
        method: String,

        /// LE/TE panel density ratio (for cosine method only)
        #[arg(long, default_value_t = 0.15)]
        le_ratio: f64,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Display geometry information
    Info {
        /// Input file path
        input: PathBuf,

        /// Output file path for JSON (if not specified, prints summary to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Plot geometry (requires 'plotting' feature)
    Plot {
        /// Input file path
        input: PathBuf,

        /// Output file (SVG or PNG based on extension)
        #[arg(short, long, default_value = "geometry.svg")]
        output: PathBuf,

        /// Show panel node ticks perpendicular to surface
        #[arg(long)]
        nodes: bool,

        /// Length of node ticks as fraction of chord
        #[arg(long, default_value_t = 0.015)]
        tick_length: f64,

        /// Plot title
        #[arg(long)]
        title: Option<String>,

        /// Image width in pixels
        #[arg(long, default_value_t = 1200)]
        width: u32,

        /// Image height in pixels
        #[arg(long, default_value_t = 400)]
        height: u32,
    },
}

fn main() {
    // Initialize logger
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Geom { action } => handle_geom(action),
        Commands::Analyze {
            file,
            alpha,
            reynolds,
            mach,
            ncrit,
            inviscid,
        } => {
            // Read geometry
            let geometry = read_geometry_auto(&file);
            let airfoil = create_paneled_airfoil(&geometry);

            // Convert angle to radians
            let alpha_rad = alpha.to_radians();

            if inviscid {
                // Inviscid-only analysis
                let solution = solve_inviscid(&airfoil);
                // Use node-based velocities for consistency with XFOIL
                let velocity = solution.velocity_at_nodes(alpha_rad);
                let cp = calculate_cp(&velocity, mach);
                let coeffs = integrate_forces(&airfoil, &cp, alpha_rad);

                println!("Inviscid Analysis Results");
                println!("========================");
                println!("Airfoil: {}", file.display());
                println!("Alpha:   {:.2}°", alpha);
                println!("Mach:    {:.3}", mach);
                println!();
                println!("CL  = {:+.6}", coeffs.cl);
                println!("CM  = {:+.6}", coeffs.cm);
                println!("CDp = {:+.6} (pressure drag)", coeffs.cdp);
            } else {
                // Viscous analysis
                let cond = FlowConditions::new(reynolds, mach, ncrit, airfoil.chord);
                let config = ViscalConfig::default();

                let result = solve_viscous(&airfoil, alpha_rad, &cond, &config);

                println!("Viscous Analysis Results");
                println!("========================");
                println!("Airfoil: {}", file.display());
                println!("Alpha:   {:.2}°", alpha);
                println!("Re:      {:.2e}", reynolds);
                println!("Mach:    {:.3}", mach);
                println!("Ncrit:   {:.1}", ncrit);
                println!();
                println!("CL  = {:+.6}", result.cl);
                println!("CD  = {:+.6}", result.cd);
                println!("  CDf = {:+.6} (friction)", result.cdf);
                println!("  CDp = {:+.6} (pressure)", result.cdp);
                println!("CM  = {:+.6}", result.cm);
                println!();
                println!("Transition:");
                println!("  Upper: {:.1}% chord", result.xtr_upper * 100.0);
                println!("  Lower: {:.1}% chord", result.xtr_lower * 100.0);
                println!();
                println!("Convergence:");
                println!("  Iterations: {}", result.iterations);
                println!("  Residual:   {:.2e}", result.residual);
                if result.converged {
                    println!("  Status:     Converged");
                } else {
                    println!("  Status:     NOT CONVERGED");
                }
            }
        }
        Commands::Polar {
            file,
            alpha_max,
            alpha_min,
            alpha_step,
            reynolds,
            mach,
            ncrit,
            json,
            output,
        } => {
            // Read geometry
            let geometry = read_geometry_auto(&file);
            let airfoil = create_paneled_airfoil(&geometry);

            // Set up polar configuration
            let conditions = FlowConditions::new(reynolds, mach, ncrit, airfoil.chord);
            let config = PolarConfig {
                alpha_max,
                alpha_min,
                alpha_step,
                conditions: conditions.clone(),
                ..Default::default()
            };

            // Run polar sweep
            let result = compute_polar(&airfoil, &config);

            // Create output struct
            let airfoil_name = file
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("Unknown");
            let polar_output = PolarOutput::from_polar(&result, airfoil_name);

            if json {
                // JSON output
                let json_str = polar_output.to_json().unwrap();
                if let Some(ref path) = output {
                    std::fs::write(path, &json_str).expect("Failed to write output file");
                    println!("Wrote JSON polar to {}", path.display());
                } else {
                    println!("{}", json_str);
                }
            } else {
                // Human-readable output
                println!("Polar Results");
                println!("=============");
                println!("Airfoil: {}", file.display());
                println!("Re:      {:.2e}", reynolds);
                println!("Mach:    {:.3}", mach);
                println!("Ncrit:   {:.1}", ncrit);
                println!();
                println!(
                    "{:>8} {:>10} {:>10} {:>10} {:>8} {:>8} {:>8} {:>5} {:>10} {:>5}",
                    "Alpha", "CL", "CD", "CM", "L/D", "Xtr_U", "Xtr_L", "Iter", "Residual", "Conv"
                );
                println!("{}", "-".repeat(98));

                for point in &polar_output.points {
                    let conv_marker = if point.converged { "Y" } else { "N" };
                    println!(
                        "{:>8.2} {:>10.5} {:>10.6} {:>10.5} {:>8.2} {:>8.3} {:>8.3} {:>5} {:>10.2e} {:>5}",
                        point.alpha_deg,
                        point.cl,
                        point.cd,
                        point.cm,
                        point.ld,
                        point.xtr_upper,
                        point.xtr_lower,
                        point.iterations,
                        point.residual,
                        conv_marker
                    );
                }

                println!();
                println!("Summary:");
                if let Some(cl_max) = polar_output.summary.cl_max {
                    println!(
                        "  CL_max = {:.4} at alpha = {:.2}°",
                        cl_max,
                        polar_output.summary.alpha_cl_max.unwrap_or(0.0)
                    );
                }
                if let Some(ld_max) = polar_output.summary.ld_max {
                    println!(
                        "  L/D_max = {:.2} at CL = {:.4}",
                        ld_max,
                        polar_output.summary.cl_at_ld_max.unwrap_or(0.0)
                    );
                }
                if let Some(cd0) = polar_output.summary.cd0 {
                    println!("  CD0 = {:.6}", cd0);
                }
                println!(
                    "  Converged: {}/{} points",
                    polar_output.summary.num_converged,
                    polar_output.summary.num_converged + polar_output.summary.num_failed
                );

                if !result.completed {
                    println!();
                    println!("Warning: Sweep stopped early due to consecutive failures");
                }

                // Write to file if requested
                if let Some(ref path) = output {
                    let json_str = polar_output.to_json().unwrap();
                    std::fs::write(path, &json_str).expect("Failed to write output file");
                    println!();
                    println!("Wrote JSON polar to {}", path.display());
                }
            }
        }
    }
}

fn handle_geom(action: GeomAction) {
    match action {
        GeomAction::Convert {
            input,
            to,
            output,
            name,
        } => {
            // Read input file (auto-detect format)
            let geometry = read_geometry_auto(&input);

            let output_path = output.unwrap_or_else(|| {
                let mut p = input.clone();
                p.set_extension(&to);
                p
            });

            match to.as_str() {
                "json" => {
                    if let Err(e) = write_geometry_to_json(&geometry, &output_path) {
                        eprintln!("Error writing JSON: {}", e);
                        std::process::exit(1);
                    }
                    println!("Wrote JSON to {}", output_path.display());
                }
                "dat" => {
                    if let Err(e) = write_dat_file(&geometry, &name, &output_path) {
                        eprintln!("Error writing DAT: {}", e);
                        std::process::exit(1);
                    }
                    println!("Wrote DAT to {}", output_path.display());
                }
                _ => {
                    eprintln!("Unknown output format: {}", to);
                    std::process::exit(1);
                }
            }
        }

        GeomAction::Naca {
            spec,
            panels,
            to,
            output,
        } => {
            let geometry = if spec.len() == 4 {
                match naca_4digit(&spec, panels) {
                    Ok(g) => g,
                    Err(e) => {
                        eprintln!("Error generating NACA airfoil: {}", e);
                        std::process::exit(1);
                    }
                }
            } else if spec.len() == 5 {
                match naca_5digit(&spec, panels) {
                    Ok(g) => g,
                    Err(e) => {
                        eprintln!("Error generating NACA airfoil: {}", e);
                        std::process::exit(1);
                    }
                }
            } else {
                eprintln!("Invalid NACA specification: {} (must be 4 or 5 digits)", spec);
                std::process::exit(1);
            };

            let name = format!("NACA {}", spec);
            let output_path = output.unwrap_or_else(|| {
                PathBuf::from(format!("naca{}.{}", spec, to))
            });

            match to.as_str() {
                "json" => {
                    if let Err(e) = write_geometry_to_json(&geometry, &output_path) {
                        eprintln!("Error writing JSON: {}", e);
                        std::process::exit(1);
                    }
                    println!("Generated NACA {} with {} panels", spec, panels);
                    println!("Wrote JSON to {}", output_path.display());
                }
                "dat" => {
                    if let Err(e) = write_dat_file(&geometry, &name, &output_path) {
                        eprintln!("Error writing DAT: {}", e);
                        std::process::exit(1);
                    }
                    println!("Generated NACA {} with {} panels", spec, panels);
                    println!("Wrote DAT to {}", output_path.display());
                }
                _ => {
                    eprintln!("Unknown output format: {}", to);
                    std::process::exit(1);
                }
            }
        }

        GeomAction::Repanel {
            input,
            panels,
            method,
            le_ratio,
            output,
        } => {
            let geometry = read_geometry_auto(&input);

            let repaneled = match method.to_lowercase().as_str() {
                "xfoil" | "pane" => {
                    // Use XFOIL's curvature-based PANE algorithm
                    let config = PaneConfig::default();
                    repanel_xfoil(&geometry, panels, &config)
                }
                "cosine" => {
                    // Use modified cosine spacing
                    repanel_cosine(&geometry, panels, le_ratio)
                }
                _ => {
                    eprintln!("Unknown repaneling method: {}. Use 'xfoil' or 'cosine'", method);
                    std::process::exit(1);
                }
            };

            let output_path = output.unwrap_or_else(|| {
                let mut p = input.clone();
                let stem = p.file_stem().unwrap().to_str().unwrap();
                p.set_file_name(format!("{}_repaneled.json", stem));
                p
            });

            if let Err(e) = write_geometry_to_json(&repaneled, &output_path) {
                eprintln!("Error writing output: {}", e);
                std::process::exit(1);
            }

            println!(
                "Repaneled from {} to {} points using {} method",
                geometry.x_c.len(),
                repaneled.x_c.len(),
                method
            );
            println!("Wrote to {}", output_path.display());
        }

        GeomAction::Info { input, output } => {
            let geometry = read_geometry_auto(&input);
            let airfoil = create_paneled_airfoil(&geometry);
            let info = yfoil::output::GeometryInfo::from_paneled(&airfoil);

            if let Some(ref path) = output {
                // Write full JSON to file
                let json_str = info.to_json().expect("Failed to serialize geometry info");
                std::fs::write(path, &json_str).expect("Failed to write output file");
                println!("Wrote geometry info to {}", path.display());
            } else {
                // Print summary to stdout
                print_geometry_info(&info.summary);
            }
        }

        GeomAction::Plot {
            input,
            output,
            nodes,
            tick_length,
            title,
            width,
            height,
        } => {
            let geometry = read_geometry_auto(&input);
            let airfoil = create_paneled_airfoil(&geometry);

            #[cfg(feature = "plotting")]
            {
                let config = GeometryPlotConfig {
                    width,
                    height,
                    show_nodes: nodes,
                    tick_length,
                    title,
                    ..Default::default()
                };

                // Determine output format from extension
                let ext = output
                    .extension()
                    .and_then(|e| e.to_str())
                    .unwrap_or("svg")
                    .to_lowercase();

                let result = match ext.as_str() {
                    "png" => plot_paneled_png(&airfoil, &output, &config),
                    _ => plot_paneled_svg(&airfoil, &output, &config),
                };

                match result {
                    Ok(()) => println!("Wrote plot to {}", output.display()),
                    Err(e) => {
                        eprintln!("Error creating plot: {}", e);
                        std::process::exit(1);
                    }
                }
            }

            #[cfg(not(feature = "plotting"))]
            {
                let _ = airfoil;
                let _ = output;
                let _ = nodes;
                let _ = tick_length;
                let _ = title;
                let _ = width;
                let _ = height;
                eprintln!(
                    "Plotting feature not enabled. Rebuild with: cargo build --features plotting"
                );
                std::process::exit(1);
            }
        }
    }
}

fn read_geometry_auto(path: &PathBuf) -> Geometry {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "json" => match read_geometry_from_file(path) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("Error reading JSON file: {}", e);
                std::process::exit(1);
            }
        },
        "dat" => match read_dat_file(path) {
            Ok((name, g)) => {
                println!("Loaded: {}", name);
                g
            }
            Err(e) => {
                eprintln!("Error reading DAT file: {}", e);
                std::process::exit(1);
            }
        },
        _ => {
            // Try JSON first, then DAT
            if let Ok(g) = read_geometry_from_file(path) {
                g
            } else if let Ok((name, g)) = read_dat_file(path) {
                println!("Loaded: {}", name);
                g
            } else {
                eprintln!("Could not read file: {}", path.display());
                std::process::exit(1);
            }
        }
    }
}

fn print_geometry_info(summary: &yfoil::output::GeometrySummary) {
    println!("Geometry Information:");
    println!("  Number of points: {}", summary.n_points);
    println!("  Chord: {:.6}", summary.chord);
    println!(
        "  X range: {:.6} to {:.6}",
        summary.x_range[0], summary.x_range[1]
    );
    println!(
        "  Y range: {:.6} to {:.6}",
        summary.y_range[0], summary.y_range[1]
    );
    println!("  Max thickness: {:.4}", summary.max_thickness);
    println!("  TE gap: {:.6}", summary.te_gap);
    println!(
        "  Sharp TE: {}",
        if summary.sharp_te { "yes" } else { "no" }
    );
    println!(
        "  Reference point: ({:.4}, {:.4})",
        summary.reference[0], summary.reference[1]
    );
    println!("  LE index: {}", summary.le_index);
    println!("  LE arc length: {:.6}", summary.sle);
    println!("  Total arc length: {:.6}", summary.total_arc_length);
    println!("  Max curvature: {:.4}", summary.max_curvature);
    println!(
        "  First point: ({:.6}, {:.6})",
        summary.first_point[0], summary.first_point[1]
    );
    println!(
        "  Last point: ({:.6}, {:.6})",
        summary.last_point[0], summary.last_point[1]
    );
}
