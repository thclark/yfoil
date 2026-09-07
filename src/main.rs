//! YFoil CLI - Aerofoil analysis tool

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use yfoil::geometry::{
    naca_4digit, naca_4digit_xfoil, naca_5digit, naca_5digit_xfoil, panel_foil, read_dat_file, read_geometry_from_file,
    repanel_by_curvature, repanel_cosine, write_dat_file, write_geometry_to_json, Geometry, PaneConfig,
};
use yfoil::output::{AnalysisOutput, InviscidAnalysisOutput, PolarOutput};
use yfoil::solver::analysis::{compute_polar, compute_polar_with, FlowConditions, PolarConfig, Session};

#[cfg(feature = "plotting")]
use yfoil::output::{
    plot_analysis_png, plot_analysis_svg, plot_foil, plot_polars_png, plot_polars_svg, AnalysisPlotConfig, BlQuantity,
    DesignPoint, FoilPlotConfig, ImageFormat, MarkerSet, OffsetScale, PanelStyle, PolarSeries, PolarsPlotConfig,
};

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
    #[command(visible_alias = "geom")]
    Geometry {
        #[command(subcommand)]
        action: GeomAction,
    },

    /// Analyze airfoil at a single operating point (placeholder)
    Analyze {
        /// Path to geometry file (JSON)
        file: PathBuf,

        /// Angle of attack in degrees
        #[arg(short, long, default_value_t = 0.0, allow_negative_numbers = true)]
        alpha: f64,
        /// Specified lift coefficient (fixed-CL mode: alpha becomes the unknown; overrides --alpha)
        #[arg(long, allow_negative_numbers = true)]
        cl: Option<f64>,

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
        /// Maximum VISCAL iterations (XFOIL ITER)
        #[arg(long, alias = "iter", default_value_t = 20)]
        iterations: usize,

        /// Output file path for JSON results
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Generate polar sweep
    Polar {
        /// Path to geometry file (JSON)
        file: PathBuf,

        /// Maximum angle of attack
        #[arg(long, default_value_t = 15.0, allow_negative_numbers = true)]
        alpha_max: f64,

        /// Minimum angle of attack
        #[arg(long, default_value_t = -5.0, allow_negative_numbers = true)]
        alpha_min: f64,

        /// Alpha step size
        #[arg(long, default_value_t = 0.5, allow_negative_numbers = true)]
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

        /// Display label for the polar (legend entry in `yfoil plot polar`); defaults to the geometry file stem
        #[arg(long)]
        label: Option<String>,
        /// Maximum VISCAL iterations (XFOIL ITER)
        #[arg(long, alias = "iter", default_value_t = 20)]
        iterations: usize,

        /// Embed the full point record (geometry, wake, boundary layer) of every sweep point in the
        /// JSON, for `yfoil plot foil polar.json --alpha ...`
        #[arg(long)]
        distributions: bool,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Plot results (analysis distributions or polars)
    Plot {
        #[command(subcommand)]
        action: PlotAction,
    },
}

/// Panel node rendering for `yfoil plot foil`
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum PanelStyleArg {
    /// Plain surface line
    None,
    /// Short outward tick at every node along the spline normal (XFOIL PANPLT)
    Notches,
    /// Black dot at every node
    Dots,
}

/// Output image format for `yfoil plot foil`
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum FormatArg {
    Svg,
    Png,
}

#[derive(Subcommand, Debug)]
enum PlotAction {
    /// Plot the foil: panels, wake, boundary-layer quantities drawn normal to the surface, and the
    /// stagnation, transition and separation points. Input: geometry files (panels only, overlaid),
    /// point analysis JSONs (one design point each; same panels), or one polar JSON written with
    /// `--distributions` (one design point per alpha)
    Foil {
        /// Geometry (.json/.dat), analysis JSON (`yfoil analyze -o`) or polar JSON (`yfoil polar --distributions -o`)
        #[arg(required = true)]
        files: Vec<PathBuf>,

        /// Polar input only: the embedded alphas to plot, comma-separated (default: all)
        #[arg(long, allow_hyphen_values = true)]
        alpha: Option<String>,

        /// Panel node rendering (default: notches for geometry input, none otherwise)
        #[arg(long, value_enum)]
        panels: Option<PanelStyleArg>,

        /// Notch length as a fraction of chord (XFOIL PANPLT: 0.01)
        #[arg(long, default_value_t = 0.01)]
        notch_length: f64,

        /// Draw the wake panels and, with --quantity, the wake band
        #[arg(long)]
        wake: bool,

        /// Boundary-layer quantities to draw normal to the surface, comma-separated:
        /// dstar theta delta h hk hs ue cf cdis ctau ctq uslp cp mass
        #[arg(short, long, value_delimiter = ',')]
        quantity: Vec<String>,

        /// Fixed offset multiplier for every quantity (1 draws δ*, θ and δ at true size); default is auto-scaling
        #[arg(long, conflicts_with = "max_offset")]
        scale: Option<f64>,

        /// Auto-scaling: the largest |value| of each quantity over all design points is drawn this
        /// fraction of chord off the surface
        #[arg(long, default_value_t = 0.1)]
        max_offset: f64,

        /// Markers to draw, comma-separated from: stagnation, transition, separation (default: all)
        #[arg(long, value_delimiter = ',', conflicts_with = "no_markers")]
        markers: Option<Vec<String>>,

        /// Draw no markers
        #[arg(long)]
        no_markers: bool,

        /// Design-point labels in input order (default: the file stem, or "stem α=…°" for polar points)
        #[arg(long)]
        label: Vec<String>,

        /// Output file (default: the first input's stem plus _foil.svg)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format (default: from the output extension; svg if ambiguous)
        #[arg(long, value_enum)]
        format: Option<FormatArg>,

        /// Plot title (default: the design-point label, or "Foil")
        #[arg(long)]
        title: Option<String>,

        /// Image width in pixels
        #[arg(long, default_value_t = 1200)]
        width: u32,

        /// Image height in pixels
        #[arg(long, default_value_t = 500)]
        height: u32,
    },

    /// Plot Cp and Ue distributions from a single-point analysis
    Analysis {
        /// Path to analysis results JSON file
        file: PathBuf,

        /// Output file (SVG or PNG based on extension)
        #[arg(short, long, default_value = "analysis.svg")]
        output: PathBuf,

        /// Plot title (auto-generated from file if not specified)
        #[arg(long)]
        title: Option<String>,

        /// Image width in pixels
        #[arg(long, default_value_t = 1200)]
        width: u32,

        /// Image height in pixels
        #[arg(long, default_value_t = 800)]
        height: u32,
    },

    /// Plot one or more polars (CL–α, CL–CD, CM–α, CD–α); several files are overlaid for comparison
    Polar {
        /// Polar JSON files (from `yfoil polar -o`)
        #[arg(required = true)]
        files: Vec<PathBuf>,

        /// Output file (SVG or PNG based on extension); defaults to the polar's stem plus .svg, or polar_comparison.svg for several
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Plot title (defaults to the polar label, or "Polar comparison")
        #[arg(long)]
        title: Option<String>,

        /// Image width in pixels
        #[arg(long, default_value_t = 1400)]
        width: u32,

        /// Image height in pixels
        #[arg(long, default_value_t = 1000)]
        height: u32,
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

        /// Close the trailing edge (zero TE gap; XFOIL's SHARP path)
        #[arg(long)]
        sharp: bool,
        /// Generator model: "exact" (NACA definition, thickness perpendicular to the camber
        /// line, YFoil's spacing) or "xfoil" (XFOIL's NACA4/NACA5: vertical thickness, 245-point
        /// buffer, then PANGEN to the requested panel count)
        #[arg(long, default_value = "exact")]
        naca_model: String,

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
}

fn main() {
    // Initialize logger
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Geometry { action } => handle_geom(action),
        Commands::Analyze {
            file,
            alpha,
            cl,
            reynolds,
            mach,
            ncrit,
            inviscid,
            output,
            iterations,
        } => {
            // Read geometry
            let geometry = read_geometry_auto(&file);
            let airfoil = panel_foil(&geometry);
            let airfoil_name = file.file_stem().and_then(|s| s.to_str()).unwrap_or("Unknown");

            // Convert angle to radians
            let alpha_rad = alpha.to_radians();

            let spec = FlowConditions {
                re: if inviscid { 0.0 } else { reynolds },
                mach,
                ncrit,
                max_iterations: iterations,
                ..FlowConditions::default()
            };
            let mut session = Session::new(&airfoil, spec.clone());
            let point = match cl {
                Some(clspec) => session.cl(clspec),
                None => session.alpha(alpha_rad),
            };
            let alpha = point.alpha.to_degrees();

            if inviscid {
                println!("Inviscid Analysis Results");
                println!("========================");
                println!("Airfoil: {}", file.display());
                println!("Alpha:   {:.2}°", alpha);
                println!("Mach:    {:.3}", mach);
                println!();
                println!("CL  = {:+.6}", point.cl);
                println!("CM  = {:+.6}", point.cm);
                println!("CDp = {:+.6} (pressure drag)", point.cd_pressure);

                // Write JSON output if requested
                if let Some(ref path) = output {
                    let n = airfoil.n;
                    let velocity: Vec<f64> = session.state.q_inviscid[1..=n].to_vec();
                    let cp: Vec<f64> = session.state.cp_inviscid[1..=n].to_vec();
                    let result = InviscidAnalysisOutput::new(
                        &airfoil,
                        &velocity,
                        &cp,
                        (point.cl, point.cm, point.cd_pressure),
                        alpha,
                        mach,
                        airfoil_name,
                    );
                    let json_str = result.to_json().expect("Failed to serialize results");
                    std::fs::write(path, &json_str).expect("Failed to write output file");
                    println!();
                    println!("Wrote JSON to {}", path.display());
                }
            } else {
                println!("Viscous Analysis Results");
                println!("========================");
                println!("Airfoil: {}", file.display());
                println!("Alpha:   {:.2}°", alpha);
                println!("Re:      {:.2e}", reynolds);
                println!("Mach:    {:.3}", mach);
                println!("Ncrit:   {:.1}", ncrit);
                println!();
                println!("CL  = {:+.6}", point.cl);
                println!("CD  = {:+.6}", point.cd);
                println!("  CDf = {:+.6} (friction)", point.cd_friction);
                println!("  CDp = {:+.6} (pressure, CD - CDf)", point.cd - point.cd_friction);
                println!("CM  = {:+.6}", point.cm);
                println!();
                println!("Transition:");
                println!("  Upper: {:.1}% chord", point.transition_upper[0] * 100.0);
                println!("  Lower: {:.1}% chord", point.transition_lower[0] * 100.0);
                println!();
                println!("Convergence:");
                println!("  Iterations: {}", point.iterations);
                println!("  rms:        {:.2e}", point.residual);
                if point.converged {
                    println!("  Status:     Converged");
                } else {
                    println!("  Status:     NOT CONVERGED");
                }

                if let Some(ref path) = output {
                    let result = AnalysisOutput::from_session(&session, &point, airfoil_name, &spec, false);
                    let json_str = result.to_json().expect("Failed to serialize results");
                    std::fs::write(path, &json_str).expect("Failed to write output file");
                    println!();
                    println!("Wrote JSON to {}", path.display());
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
            label,
            output,
            iterations,
            distributions,
        } => {
            // Read geometry
            let geometry = read_geometry_auto(&file);
            let airfoil = panel_foil(&geometry);
            let airfoil_name = file.file_stem().and_then(|s| s.to_str()).unwrap_or("Unknown");

            // Set up polar configuration
            let config = PolarConfig {
                alpha_max,
                alpha_min,
                alpha_step,
                conditions: FlowConditions {
                    re: reynolds,
                    mach,
                    ncrit,
                    max_iterations: iterations,
                    ..FlowConditions::default()
                },
                ..Default::default()
            };

            // Run polar sweep, capturing every visited point's state when asked to
            let mut records: Vec<AnalysisOutput> = Vec::new();
            let result = if distributions {
                let spec = config.conditions.clone();
                compute_polar_with(&airfoil, &config, &mut |session, p| {
                    records.push(AnalysisOutput::from_session(session, p, airfoil_name, &spec, false));
                })
            } else {
                compute_polar(&airfoil, &config)
            };
            records.sort_by(|a, b| a.result.alpha_deg.partial_cmp(&b.result.alpha_deg).unwrap());

            // Create output struct
            let mut polar_output = PolarOutput::from_polar(&result, airfoil_name);
            polar_output.label = label;
            polar_output.distributions = records;

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
        Commands::Plot { action } => handle_plot(action),
    }
}

#[cfg(feature = "plotting")]
fn handle_plot(action: PlotAction) {
    match action {
        PlotAction::Foil {
            files,
            alpha,
            panels,
            notch_length,
            wake,
            quantity,
            scale,
            max_offset,
            markers,
            no_markers,
            label,
            output,
            format,
            title,
            width,
            height,
        } => {
            let inputs: Vec<FoilInput> = files.iter().map(read_foil_input).collect();
            let stem = |f: &PathBuf| f.file_stem().and_then(|s| s.to_str()).unwrap_or("foil").to_string();
            let alphas: Option<Vec<f64>> = alpha.as_deref().map(|a| parse_list(a, "alpha"));

            let geometry_only = inputs.iter().all(|i| matches!(i, FoilInput::Geometry(_)));
            let all_analysis = inputs.iter().all(|i| matches!(i, FoilInput::Analysis(_)));
            let all_polar = inputs.iter().all(|i| matches!(i, FoilInput::Polar(_)));
            if !(geometry_only || all_analysis || all_polar) {
                fail("input files must all be of one kind: geometry, point analysis JSON, or a polar JSON");
            }
            if alphas.is_some() && !all_polar {
                fail("--alpha selects embedded operating points and needs a polar JSON input");
            }
            if geometry_only {
                if !quantity.is_empty() {
                    fail("geometry input has no boundary layer: --quantity needs analysis or polar JSON input");
                }
                if wake {
                    fail("geometry input has no wake: --wake needs analysis or polar JSON input");
                }
                if markers.is_some() {
                    fail("geometry input has no boundary layer: --markers needs analysis or polar JSON input");
                }
            }

            let mut points: Vec<DesignPoint> = Vec::new();
            match inputs.len() {
                n if all_polar && n > 1 => fail("give one polar JSON; several design points come from --alpha"),
                _ => {}
            }
            for (file, input) in files.iter().zip(inputs) {
                match input {
                    FoilInput::Geometry(g) => points.push(DesignPoint::from_geometry(stem(file), &panel_foil(&g))),
                    FoilInput::Analysis(a) => points.push(DesignPoint::from_analysis(stem(file), *a)),
                    FoilInput::Polar(p) => match DesignPoint::from_polar(&stem(file), &p, alphas.as_deref()) {
                        Ok(ps) => points.extend(ps),
                        Err(e) => fail(&e.to_string()),
                    },
                }
            }
            for (p, l) in points.iter_mut().zip(&label) {
                p.label = l.clone();
            }

            let quantities: Vec<BlQuantity> = quantity
                .iter()
                .map(|q| q.parse::<BlQuantity>().unwrap_or_else(|e| fail(&e)))
                .collect();
            let markers = if no_markers || geometry_only {
                MarkerSet::NONE
            } else if let Some(list) = markers {
                let mut m = MarkerSet::NONE;
                for name in list {
                    match name.trim().to_ascii_lowercase().as_str() {
                        "stagnation" => m.stagnation = true,
                        "transition" => m.transition = true,
                        "separation" => m.separation = true,
                        other => fail(&format!(
                            "unknown marker '{other}' (expected stagnation, transition or separation)"
                        )),
                    }
                }
                m
            } else {
                MarkerSet::ALL
            };
            let panels = match panels {
                Some(PanelStyleArg::None) => PanelStyle::None,
                Some(PanelStyleArg::Notches) => PanelStyle::Notches,
                Some(PanelStyleArg::Dots) => PanelStyle::Dots,
                None if geometry_only => PanelStyle::Notches,
                None => PanelStyle::None,
            };
            let title = title.or_else(|| if all_polar { Some(stem(&files[0])) } else { None });
            let config = FoilPlotConfig {
                width,
                height,
                title,
                panels,
                notch_length,
                show_wake: wake,
                quantities,
                scale: match scale {
                    Some(k) => OffsetScale::Fixed(k),
                    None => OffsetScale::Auto { max_offset },
                },
                markers,
                ..Default::default()
            };

            let output = output.unwrap_or_else(|| PathBuf::from(format!("{}_foil.svg", stem(&files[0]))));
            let format = match format {
                Some(FormatArg::Svg) => ImageFormat::Svg,
                Some(FormatArg::Png) => ImageFormat::Png,
                None => ImageFormat::from_path(&output),
            };
            match plot_foil(&points, &output, format, &config) {
                Ok(()) => println!("Wrote foil plot to {}", output.display()),
                Err(e) => fail(&format!("Error creating plot: {e}")),
            }
        }
        PlotAction::Analysis {
            file,
            output,
            title,
            width,
            height,
        } => {
            let json_str = match std::fs::read_to_string(&file) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("Error reading file: {}", e);
                    std::process::exit(1);
                }
            };

            let analysis: InviscidAnalysisOutput = match serde_json::from_str(&json_str) {
                Ok(a) => a,
                Err(e) => {
                    eprintln!("Error parsing JSON: {}", e);
                    eprintln!("Make sure the file is an inviscid analysis output (from 'yfoil analyze --inviscid -o')");
                    std::process::exit(1);
                }
            };

            let config = AnalysisPlotConfig {
                width,
                height,
                title,
                ..Default::default()
            };

            let result = if is_png(&output) {
                plot_analysis_png(&analysis, &output, &config)
            } else {
                plot_analysis_svg(&analysis, &output, &config)
            };

            match result {
                Ok(()) => println!("Wrote analysis plot to {}", output.display()),
                Err(e) => {
                    eprintln!("Error creating plot: {}", e);
                    std::process::exit(1);
                }
            }
        }
        PlotAction::Polar {
            files,
            output,
            title,
            width,
            height,
        } => {
            let mut series = Vec::with_capacity(files.len());
            for file in &files {
                let json_str = match std::fs::read_to_string(file) {
                    Ok(s) => s,
                    Err(e) => {
                        eprintln!("Error reading {}: {}", file.display(), e);
                        std::process::exit(1);
                    }
                };
                let polar: PolarOutput = match serde_json::from_str(&json_str) {
                    Ok(p) => p,
                    Err(e) => {
                        eprintln!("Error parsing {}: {}", file.display(), e);
                        eprintln!("Make sure the file is a polar output (from 'yfoil polar -o')");
                        std::process::exit(1);
                    }
                };
                let s = PolarSeries::from_polar_output(&polar);
                if s.alpha.is_empty() {
                    eprintln!("Warning: {} has no converged points", file.display());
                }
                series.push(s);
            }

            let output = output.unwrap_or_else(|| {
                if files.len() == 1 {
                    files[0].with_extension("svg")
                } else {
                    PathBuf::from("polar_comparison.svg")
                }
            });
            let config = PolarsPlotConfig {
                width,
                height,
                title,
                ..Default::default()
            };

            let result = if is_png(&output) {
                plot_polars_png(&series, &output, &config)
            } else {
                plot_polars_svg(&series, &output, &config)
            };

            match result {
                Ok(()) => println!("Wrote polar plot to {}", output.display()),
                Err(e) => {
                    eprintln!("Error creating plot: {}", e);
                    std::process::exit(1);
                }
            }
        }
    }
}

/// One input file of `yfoil plot foil`
#[cfg(feature = "plotting")]
enum FoilInput {
    Geometry(Geometry),
    Analysis(Box<AnalysisOutput>),
    Polar(Box<PolarOutput>),
}

/// Geometry (.dat, or JSON with `x_c`), analysis JSON (`result` + `geometry`) or polar JSON
/// (`points`), by trying the parsers in that order of specificity
#[cfg(feature = "plotting")]
fn read_foil_input(path: &PathBuf) -> FoilInput {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
    if ext == "dat" {
        return FoilInput::Geometry(read_geometry_auto(path));
    }
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("Error reading {}: {e}", path.display())));
    if let Ok(p) = serde_json::from_str::<PolarOutput>(&text) {
        return FoilInput::Polar(Box::new(p));
    }
    if let Ok(a) = serde_json::from_str::<AnalysisOutput>(&text) {
        return FoilInput::Analysis(Box::new(a));
    }
    if let Ok(g) = serde_json::from_str::<Geometry>(&text) {
        return FoilInput::Geometry(g);
    }
    fail(&format!(
        "{}: not a geometry, analysis or polar JSON (analysis JSON must come from this version's `yfoil analyze -o`)",
        path.display()
    ))
}

/// Comma-separated numbers
#[cfg(feature = "plotting")]
fn parse_list(text: &str, what: &str) -> Vec<f64> {
    text.split(',')
        .map(|t| {
            t.trim()
                .parse::<f64>()
                .unwrap_or_else(|_| fail(&format!("--{what}: '{t}' is not a number")))
        })
        .collect()
}

#[cfg(feature = "plotting")]
fn fail(msg: &str) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

#[cfg(feature = "plotting")]
fn is_png(path: &std::path::Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("png"))
}

#[cfg(not(feature = "plotting"))]
fn handle_plot(_action: PlotAction) {
    eprintln!("Plotting feature not enabled. Rebuild with: cargo build --features plotting");
    std::process::exit(1);
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
            sharp,
            naca_model,
            output,
        } => {
            let geometry = if naca_model == "xfoil" {
                let buffer = if spec.len() == 4 {
                    naca_4digit_xfoil(&spec)
                } else {
                    naca_5digit_xfoil(&spec)
                };
                match buffer {
                    Ok(b) => repanel_by_curvature(&b, panels, &PaneConfig::default()),
                    Err(e) => {
                        eprintln!("Error generating NACA airfoil: {}", e);
                        std::process::exit(1);
                    }
                }
            } else if spec.len() == 4 {
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
            let geometry = if sharp { geometry.sharpen() } else { geometry };

            let name = format!("NACA {}", spec);
            let output_path = output.unwrap_or_else(|| PathBuf::from(format!("naca{}.{}", spec, to)));

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
                    repanel_by_curvature(&geometry, panels, &config)
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
            let airfoil = panel_foil(&geometry);
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
    }
}

fn read_geometry_auto(path: &PathBuf) -> Geometry {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();

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
    println!("  X range: {:.6} to {:.6}", summary.x_range[0], summary.x_range[1]);
    println!("  Y range: {:.6} to {:.6}", summary.y_range[0], summary.y_range[1]);
    println!("  Max thickness: {:.4}", summary.max_thickness);
    println!("  TE gap: {:.6}", summary.te_gap);
    println!("  Sharp TE: {}", if summary.sharp_te { "yes" } else { "no" });
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
