//! YFoil CLI - Aerofoil analysis tool

use std::path::PathBuf;

use clap::{Parser, Subcommand};

use yfoil::bl::FlowConditions;
use yfoil::forces::{calculate_cp, integrate_forces};
use yfoil::geometry::{
    create_paneled_airfoil, naca_4digit, naca_5digit, read_dat_file, read_geometry_from_file,
    repanel, write_dat_file, write_geometry_to_json, Geometry,
};
use yfoil::panel::solve_inviscid;
use yfoil::output::PolarOutput;
use yfoil::solver::{compute_polar, solve_viscous, PolarConfig, ViscalConfig};

#[cfg(feature = "plotting")]
use plotly::{Plot, Scatter};

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

        /// LE/TE panel density ratio
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
    },

    /// Plot geometry (requires 'plotting' feature)
    Plot {
        /// Input file path
        input: PathBuf,

        /// Output HTML file
        #[arg(short, long, default_value = "geometry.html")]
        output: PathBuf,
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
                let velocity = solution.velocity_at_alpha(alpha_rad);
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
                    "{:>8} {:>10} {:>10} {:>10} {:>8} {:>8} {:>8} {:>5} {:>10}",
                    "Alpha", "CL", "CD", "CM", "L/D", "Xtr_U", "Xtr_L", "Iter", "Residual"
                );
                println!("{}", "-".repeat(92));

                for point in &polar_output.points {
                    println!(
                        "{:>8.2} {:>10.5} {:>10.6} {:>10.5} {:>8.2} {:>8.3} {:>8.3} {:>5} {:>10.2e}",
                        point.alpha_deg,
                        point.cl,
                        point.cd,
                        point.cm,
                        point.ld,
                        point.xtr_upper,
                        point.xtr_lower,
                        point.iterations,
                        point.residual
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
            le_ratio,
            output,
        } => {
            let geometry = read_geometry_auto(&input);
            let repaneled = repanel(&geometry, panels, le_ratio);

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
                "Repaneled from {} to {} points",
                geometry.x_c.len(),
                repaneled.x_c.len()
            );
            println!("Wrote to {}", output_path.display());
        }

        GeomAction::Info { input } => {
            let geometry = read_geometry_auto(&input);
            print_geometry_info(&geometry);
        }

        GeomAction::Plot { input, output } => {
            let geometry = read_geometry_auto(&input);

            #[cfg(feature = "plotting")]
            {
                let mut plot = Plot::new();
                let trace = Scatter::new(geometry.x_c.clone(), geometry.y_c.clone());
                plot.add_trace(trace);
                plot.write_html(&output);
                println!("Wrote plot to {}", output.display());
            }

            #[cfg(not(feature = "plotting"))]
            {
                let _ = geometry;
                let _ = output;
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

fn print_geometry_info(geometry: &Geometry) {
    let n = geometry.x_c.len();
    let max_x = geometry.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_x = geometry.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_y = geometry.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y = geometry.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

    println!("Geometry Information:");
    println!("  Number of points: {}", n);
    println!("  X range: {:.6} to {:.6}", min_x, max_x);
    println!("  Y range: {:.6} to {:.6}", min_y, max_y);
    println!("  Max thickness: {:.4}", max_y - min_y);
    println!(
        "  Reference point: ({:.4}, {:.4})",
        geometry.reference[0], geometry.reference[1]
    );
    println!("  First point: ({:.6}, {:.6})", geometry.x_c[0], geometry.y_c[0]);
    println!(
        "  Last point: ({:.6}, {:.6})",
        geometry.x_c[n - 1],
        geometry.y_c[n - 1]
    );
}
