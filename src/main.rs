//! yFoil CLI - Aerofoil analysis tool

use std::path::{Path, PathBuf};
use yfoil::geometry::Thickness;

use clap::{Parser, Subcommand, ValueEnum};

use yfoil::geometry::{
    naca_4digit_vertical, naca_5digit_vertical, panel_foil, read_dat_file, read_geometry_from_file, repanel,
    write_dat_file, write_geometry_to_json, CosineConfig, Geometry, KarmanTrefftz, PanelConfig, PanelMethod,
    PangenConfig, RepanelError, Section, Series, TeGap,
};
use yfoil::output::{AnalysisOutput, PolarOutput};
use yfoil::solver::analysis::{compute_polar, compute_polar_with, FlowConditions, PolarConfig, Session};

#[cfg(feature = "plotting")]
use yfoil::output::{
    plot_analysis_png, plot_analysis_svg, plot_foil, plot_polars_png, plot_polars_svg, AnalysisPlotConfig, BlQuantity,
    DesignPoint, FoilPlotConfig, ImageFormat, MarkerSet, OffsetScale, PanelStyle, PolarSeries, PolarsPlotConfig,
};

/// yFoil - Rust-based aerofoil analysis tool
#[derive(Parser, Debug)]
#[command(name = "yfoil")]
#[command(version, about = "Rust-based aerofoil analysis tool", long_about = None)]
#[command(arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Geometry operations (convert, generate, repanel)
    #[command(visible_alias = "geom", arg_required_else_help = true)]
    Geometry {
        #[command(subcommand)]
        action: GeomAction,
    },

    /// Analyse an aerofoil at a single operating point
    #[command(arg_required_else_help = true)]
    Analyse {
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

        /// Critical amplification factor for transition (Ncrit)
        #[arg(long, default_value_t = 9.0)]
        ncrit: f64,

        /// Inviscid analysis only
        #[arg(long)]
        inviscid: bool,
        /// Maximum VISCAL iterations (XFOIL ITER)
        #[arg(long, default_value_t = 20)]
        max_iterations: usize,

        /// Output file path for JSON results
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Also emit XFOIL's lagged closure arrays under each side's `lagged_closures`
        #[arg(long)]
        include_lagged_closures: bool,
    },

    /// Generate polar sweep
    #[command(arg_required_else_help = true)]
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

        /// Critical amplification factor for transition (Ncrit)
        #[arg(long, default_value_t = 9.0)]
        ncrit: f64,

        /// Output JSON format
        #[arg(long)]
        json: bool,

        /// Display label for the polar (legend entry in `yfoil plot polar`); defaults to the geometry file stem
        #[arg(long)]
        label: Option<String>,
        /// Maximum VISCAL iterations (XFOIL ITER)
        #[arg(long, default_value_t = 20)]
        max_iterations: usize,

        /// Embed the full analysis record (geometry, wake, boundary layer) of every sweep point in
        /// the JSON, for `yfoil plot foil polar.json --alpha ...`
        #[arg(long)]
        distributions: bool,

        /// Also emit XFOIL's lagged closure arrays in the embedded records
        #[arg(long)]
        include_lagged_closures: bool,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Plot results (analysis distributions or polars)
    #[command(arg_required_else_help = true)]
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
    #[command(arg_required_else_help = true)]
    Foil {
        /// Geometry (.json/.dat), analysis JSON (`yfoil analyse -o`) or polar JSON (`yfoil polar --distributions -o`)
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
    #[command(arg_required_else_help = true)]
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
    #[command(arg_required_else_help = true)]
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

/// The node-distribution methods of `yfoil geometry`
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum PanellingMethod {
    /// XFOIL's PANGEN, the algorithm behind its PANE and PPAR commands (the default)
    Pangen,
    /// yFoil's own cosine spacing: no XFOIL equivalent
    Cosine,
}

/// How the nodes are distributed: the same options on `repanel` and on the generators, either as
/// flags or as one JSON file (`--panelling`) in the shape the geometry file records them.
#[derive(clap::Args, Debug, Clone)]
struct PanellingArgs {
    /// Number of panel nodes, trailing edge round to trailing edge (XFOIL's NPAN, "number of
    /// panel nodes"). The cosine repanelling of an existing geometry writes N + 1 nodes, its
    /// historic behaviour. [both methods]
    #[arg(
        short = 'n',
        long = "nodes",
        value_name = "N",
        default_value_t = 160,
        conflicts_with = "panelling",
        help_heading = "Panelling (both methods)"
    )]
    nodes: usize,

    /// Node-distribution method. `pangen` is XFOIL's PANGEN (its PANE / PPAR commands), a
    /// curvature-weighted spacing tuned by the PPAR parameters below. `cosine` is yFoil's own
    /// cosine spacing and has no XFOIL equivalent: on a generator it samples the analytic
    /// section at cosine chord stations; on `repanel` it is an arc-length cosine with the
    /// --cosine-te-bias warp. [both methods]
    #[arg(long, value_enum, default_value_t = PanellingMethod::Pangen, conflicts_with = "panelling", help_heading = "Panelling (both methods)")]
    method: PanellingMethod,

    /// Read the whole panelling from a JSON file instead of flags: the same shape a geometry
    /// file records under generator.panelling, e.g. {"method": "pangen", "n_nodes": 160,
    /// "sharp_te": false, "te_gap": null, "curvature_bunching": 1.0, ...}. Cannot be combined
    /// with any other panelling flag. [both methods]
    #[arg(long, value_name = "FILE", help_heading = "Panelling (both methods)")]
    panelling: Option<PathBuf>,

    /// Curvature bunching parameter (XFOIL's CVPAR, PPAR menu P): the curvature attraction
    /// coefficient is 6 × this, 0 gives uniform arc-length spacing; default 1.0 [pangen]
    #[arg(
        long,
        value_name = "P",
        conflicts_with = "panelling",
        help_heading = "PANGEN parameters (--method pangen only; an error with cosine)"
    )]
    curvature_bunching: Option<f64>,

    /// Fictitious trailing-edge curvature as a fraction of the leading-edge curvature, bunching
    /// nodes at the trailing edge (XFOIL's CTERAT, PPAR menu T, "TE/LE panel density ratio");
    /// default 0.15 [pangen]
    #[arg(
        long,
        value_name = "T",
        conflicts_with = "panelling",
        help_heading = "PANGEN parameters (--method pangen only; an error with cosine)"
    )]
    te_curvature_ratio: Option<f64>,

    /// Fictitious curvature inside the refinement windows as a fraction of the leading-edge
    /// curvature (XFOIL's CTRRAT, PPAR menu R); default 0.2 [pangen]
    #[arg(
        long,
        value_name = "R",
        conflicts_with = "panelling",
        help_heading = "PANGEN parameters (--method pangen only; an error with cosine)"
    )]
    refined_curvature_ratio: Option<f64>,

    /// Upper-surface refinement window as two x/c limits, e.g. 0.2,0.4 (XFOIL's XSREF1 and
    /// XSREF2, PPAR menu XT); off by default [pangen]
    #[arg(long, value_name = "X1,X2", value_delimiter = ',', num_args = 1.., conflicts_with = "panelling", help_heading = "PANGEN parameters (--method pangen only; an error with cosine)")]
    refine_upper: Option<Vec<f64>>,

    /// Lower-surface refinement window as two x/c limits, e.g. 0.3,0.6 (XFOIL's XPREF1 and
    /// XPREF2, PPAR menu XB); off by default [pangen]
    #[arg(long, value_name = "X1,X2", value_delimiter = ',', num_args = 1.., conflicts_with = "panelling", help_heading = "PANGEN parameters (--method pangen only; an error with cosine)")]
    refine_lower: Option<Vec<f64>>,

    /// Warp of the arc-length cosine when repanelling an existing geometry: 1 is a plain cosine,
    /// below 1 coarser at the trailing edge and finer at the leading edge, above 1 finer at the
    /// trailing edge (clamped 0.05…2); default 0.15. A generator's cosine sampling has no bias,
    /// so this is an error there. [cosine]
    #[arg(
        long,
        value_name = "B",
        conflicts_with = "panelling",
        help_heading = "Cosine parameters (--method cosine only; an error with pangen)"
    )]
    cosine_te_bias: Option<f64>,

    /// After panelling, move the two trailing-edge nodes to their midpoint: a closed edge,
    /// XFOIL's SHARP path. Cannot be combined with --te-gap. [both methods]
    #[arg(long, conflicts_with_all = ["te_gap", "te_blend", "panelling"], help_heading = "Trailing edge (both methods, applied after panelling)")]
    sharp: bool,

    /// After panelling, set the trailing-edge gap in chord units with XFOIL's TGAP: the surfaces
    /// move apart along the gap direction by ½Δ·(x/c)·exp(−(1 − x/c)(1/blend − 1)) each. Note
    /// XFOIL applies TGAP to the buffer airfoil before PANE; yFoil applies it to the panelled
    /// nodes so the blend profile is exact at every node. [both methods]
    #[arg(
        long,
        value_name = "GAP",
        conflicts_with = "panelling",
        help_heading = "Trailing edge (both methods, applied after panelling)"
    )]
    te_gap: Option<f64>,

    /// Blending distance/c of --te-gap, 0..1 (TGAP's second argument); default 1.0 [both methods]
    #[arg(
        long,
        value_name = "F",
        requires = "te_gap",
        conflicts_with = "panelling",
        help_heading = "Trailing edge (both methods, applied after panelling)"
    )]
    te_blend: Option<f64>,
}

impl PanellingArgs {
    /// The panelling these flags (or the file) describe. `generator` says whether the caller
    /// samples an analytic section (where a cosine bias does not apply and `n_buffer_nodes` may)
    /// or repanels an existing geometry.
    fn config(&self, generator: bool) -> Result<PanelConfig, String> {
        if let Some(file) = &self.panelling {
            let text = std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?;
            let value: serde_json::Value =
                serde_json::from_str(&text).map_err(|e| format!("{}: {e}", file.display()))?;
            let config = PanelConfig::from_json(&value).map_err(|e| format!("{}: {e}", file.display()))?;
            if !generator && config.n_buffer_nodes.is_some() {
                return Err(format!("{}: {}", file.display(), RepanelError::BufferNodesOnRepanel));
            }
            if generator && matches!(config.method, PanelMethod::Cosine(CosineConfig { te_bias: Some(_) })) {
                return Err(format!(
                    "{}: te_bias applies to repanelling an existing geometry; a generator's cosine sampling has no bias",
                    file.display()
                ));
            }
            return Ok(config);
        }
        let pangen_flags = [
            ("--curvature-bunching", self.curvature_bunching.is_some()),
            ("--te-curvature-ratio", self.te_curvature_ratio.is_some()),
            ("--refined-curvature-ratio", self.refined_curvature_ratio.is_some()),
            ("--refine-upper", self.refine_upper.is_some()),
            ("--refine-lower", self.refine_lower.is_some()),
        ];
        let window = |name: &str, v: &Option<Vec<f64>>| -> Result<Option<[f64; 2]>, String> {
            match v {
                None => Ok(None),
                Some(v) if v.len() == 2 => Ok(Some([v[0], v[1]])),
                Some(v) => Err(format!("{name}: two x/c values are needed (X1,X2), got {}", v.len())),
            }
        };
        let method = match self.method {
            PanellingMethod::Pangen => {
                if self.cosine_te_bias.is_some() {
                    return Err("--cosine-te-bias belongs to --method cosine, not pangen".into());
                }
                let d = PangenConfig::default();
                PanelMethod::Pangen(PangenConfig {
                    curvature_bunching: self.curvature_bunching.unwrap_or(d.curvature_bunching),
                    te_curvature_ratio: self.te_curvature_ratio.unwrap_or(d.te_curvature_ratio),
                    refined_curvature_ratio: self.refined_curvature_ratio.unwrap_or(d.refined_curvature_ratio),
                    refine_upper: window("--refine-upper", &self.refine_upper)?,
                    refine_lower: window("--refine-lower", &self.refine_lower)?,
                })
            }
            PanellingMethod::Cosine => {
                if let Some((flag, _)) = pangen_flags.iter().find(|(_, given)| *given) {
                    return Err(format!("{flag} belongs to --method pangen, not cosine"));
                }
                if generator && self.cosine_te_bias.is_some() {
                    return Err(
                        "--cosine-te-bias applies to repanelling an existing geometry; a generator's cosine sampling has no bias"
                            .into(),
                    );
                }
                PanelMethod::Cosine(CosineConfig {
                    te_bias: self.cosine_te_bias,
                })
            }
        };
        let config = PanelConfig {
            n_nodes: self.nodes,
            sharp_te: self.sharp,
            te_gap: self.te_gap.map(|gap| TeGap {
                gap,
                blend: self.te_blend.unwrap_or(1.0),
            }),
            method,
            n_buffer_nodes: None,
        };
        config.validate().map_err(|e| e.to_string())?;
        Ok(config)
    }
}

#[derive(Subcommand, Debug)]
enum GeomAction {
    /// Convert between geometry formats: yFoil's JSON (with its provenance record) and the
    /// Selig `.dat` XFOIL loads. Coordinates round-trip bit-exactly (17 significant figures).
    #[command(arg_required_else_help = true)]
    Convert {
        /// Input file path
        input: PathBuf,

        /// Output format: json, dat
        #[arg(long, default_value = "json")]
        to: String,

        /// Output file path
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Aerofoil name (for .dat output)
        #[arg(long, default_value = "Aerofoil")]
        name: String,
    },

    /// Generate a NACA section from its designation, panelled with XFOIL's PANGEN by default.
    ///
    /// Every family with a public definition: 4-digit (2412), 4-digit modified (0012-34),
    /// 5-digit (23012; 23112 reflex), 16-series (16-212), 6-series (63-415, --a for the loading
    /// extent) and 6A-series (64A010). The thickness is laid perpendicular to the mean line, the
    /// NACA definition; --thickness vertical is XFOIL's own NACA4/NACA5 model and exists only to
    /// replicate XFOIL's NACA command. The panelling flags are shared with `repanel`; the output
    /// file records the section and the panelling under "generator".
    #[command(arg_required_else_help = true)]
    Naca {
        /// NACA designation (e.g., "0012", "4412", "23015", "16-212", "63-415", "64A010")
        spec: String,

        /// Output format: json, dat
        #[arg(long, default_value = "json")]
        to: String,

        /// How the thickness form is laid on the mean line. perpendicular: the NACA definition,
        /// the default with either panelling method. vertical: XFOIL's NACA4/NACA5 model (XFOIL
        /// adds the thickness vertically, docs/xfoil-known-issues.md §6.1) on XFOIL's own
        /// 245-point buffer, always PANGEN-panelled; 4- and 5-digit only, and only for replicating
        /// XFOIL's NACA command output
        #[arg(long, value_enum, default_value_t = Thickness::Perpendicular)]
        thickness: Thickness,

        /// Extent of uniform loading of a 6-series (or 16-series) mean line, 0..1 (default 1.0)
        #[arg(long)]
        a: Option<f64>,

        #[command(flatten)]
        panelling: PanellingArgs,

        /// Output file path (default naca<designation>.<format>)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Generate a Kármán–Trefftz section, panelled with XFOIL's PANGEN by default.
    ///
    /// The conformal map of a circle through ζ = 1 (Joukowski when the trailing-edge angle is
    /// 0): an analytic section with an exact potential-flow solution and a sharp trailing edge.
    /// The circle centre sets thickness (x) and camber (y). The panelling flags are shared with
    /// `repanel`; the output file records the parameters and the panelling under "generator".
    #[command(arg_required_else_help = true)]
    KarmanTrefftz {
        /// Circle centre x in the ζ-plane (negative; sets the thickness)
        #[arg(long, default_value_t = -0.1)]
        x_centre: f64,

        /// Circle centre y in the ζ-plane (sets the camber)
        #[arg(long, default_value_t = 0.05)]
        y_centre: f64,

        /// Trailing-edge angle in degrees, 0 ≤ τ < 180
        #[arg(long, default_value_t = 10.0)]
        te_angle: f64,

        /// Output format: json, dat
        #[arg(long, default_value = "json")]
        to: String,

        #[command(flatten)]
        panelling: PanellingArgs,

        /// Output file path (default karman-trefftz.<format>)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Redistribute the nodes of an existing geometry: XFOIL's PANE / PPAR on a loaded airfoil.
    ///
    /// The default method is XFOIL's PANGEN, translated line for line: nodes are placed so that
    /// (1 + 6·CVPAR·curvature)·Δs is the same on every panel, with a fictitious curvature added at
    /// the trailing edge (CTERAT) and inside optional refinement windows (CTRRAT, XSREF, XPREF)
    /// and the curvature field smoothed first. The PPAR parameters are the flags below, by
    /// descriptive names with the XFOIL names in their help. `--method cosine` is yFoil's own
    /// arc-length cosine spacing with no XFOIL equivalent. The input is a `.json` or `.dat`
    /// geometry; the output is JSON only, <stem>_repanelled.json beside the input by default,
    /// carrying the panelling under "generator".
    #[command(arg_required_else_help = true)]
    Repanel {
        /// Input geometry file (.json or .dat)
        input: PathBuf,

        #[command(flatten)]
        panelling: PanellingArgs,

        /// Output file path (JSON; default <input stem>_repanelled.json beside the input)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Display geometry information: node count, chord, extents, thickness, trailing-edge gap,
    /// leading edge, arc length and curvature; with -o the full per-node distributions as JSON
    #[command(arg_required_else_help = true)]
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
        Commands::Analyse {
            file,
            alpha,
            cl,
            reynolds,
            mach,
            ncrit,
            inviscid,
            output,
            max_iterations,
            include_lagged_closures,
        } => {
            let geometry = read_geometry_auto(&file);
            let foil = panel_foil(&geometry);
            let foil_name = file.file_stem().and_then(|s| s.to_str()).unwrap_or("Unknown");

            let conditions = FlowConditions {
                re: if inviscid { None } else { Some(reynolds) },
                mach,
                ncrit,
                max_iterations,
                ..FlowConditions::default()
            };
            let mut session = Session::new(&foil, conditions.clone());
            let point = match cl {
                Some(clspec) => session.cl(clspec),
                None => session.alpha(alpha.to_radians()),
            };
            let alpha_deg = point.alpha.to_degrees();

            println!(
                "{}",
                if inviscid {
                    "Inviscid Analysis Results"
                } else {
                    "Viscous Analysis Results"
                }
            );
            println!("========================");
            println!("Foil:    {}", file.display());
            println!("Alpha:   {:.2}°", alpha_deg);
            if let Some(re) = conditions.re {
                println!("Re:      {:.2e}", re);
            }
            println!("Mach:    {:.3}", mach);
            if !inviscid {
                println!("Ncrit:   {:.1}", ncrit);
            }
            println!();
            println!("CL  = {:+.6}", point.cl);
            if inviscid {
                println!("CM  = {:+.6}", point.cm);
                println!("CDp = {:+.6} (pressure drag)", point.cd_pressure);
            } else {
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
                println!(
                    "  Status:     {}",
                    if point.converged { "Converged" } else { "NOT CONVERGED" }
                );
            }

            if let Some(ref path) = output {
                let result = AnalysisOutput::from_session(&session, &point, foil_name, include_lagged_closures);
                let json_str = result.to_json().expect("Failed to serialize results");
                std::fs::write(path, &json_str).expect("Failed to write output file");
                println!();
                println!("Wrote JSON to {}", path.display());
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
            max_iterations,
            distributions,
            include_lagged_closures,
        } => {
            // Read geometry
            let geometry = read_geometry_auto(&file);
            let aerofoil = panel_foil(&geometry);
            let airfoil_name = file.file_stem().and_then(|s| s.to_str()).unwrap_or("Unknown");

            // Set up polar configuration
            let config = PolarConfig {
                alpha_max,
                alpha_min,
                alpha_step,
                conditions: FlowConditions {
                    re: Some(reynolds),
                    mach,
                    ncrit,
                    max_iterations,
                    ..FlowConditions::default()
                },
                ..Default::default()
            };

            // Run polar sweep, capturing every visited point's state when asked to
            let mut records: Vec<AnalysisOutput> = Vec::new();
            let result = if distributions {
                compute_polar_with(&aerofoil, &config, &mut |session, p| {
                    records.push(AnalysisOutput::from_session(
                        session,
                        p,
                        airfoil_name,
                        include_lagged_closures,
                    ));
                })
            } else {
                compute_polar(&aerofoil, &config)
            };
            records.sort_by(|a, b| a.results.alpha_deg.partial_cmp(&b.results.alpha_deg).unwrap());

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
                println!("Aerofoil: {}", file.display());
                println!("Re:      {:.2e}", reynolds);
                println!("Mach:    {:.3}", mach);
                println!("Ncrit:   {:.1}", ncrit);
                println!();
                println!(
                    "{:>8} {:>10} {:>10} {:>10} {:>8} {:>8} {:>8} {:>5} {:>10} {:>5}",
                    "Alpha", "CL", "CD", "CM", "L/D", "Xtr_U", "Xtr_L", "Iter", "Residual", "Conv"
                );
                println!("{}", "-".repeat(98));

                for point in &polar_output.results {
                    let conv_marker = if point.is_converged() { "Y" } else { "N" };
                    println!(
                        "{:>8.2} {:>10.5} {:>10.6} {:>10.5} {:>8.2} {:>8.3} {:>8.3} {:>5} {:>10.2e} {:>5}",
                        point.alpha_deg,
                        point.cl,
                        point.cd.unwrap_or(0.0),
                        point.cm,
                        point.ldratio.unwrap_or(0.0),
                        point.transition_upper.map_or(0.0, |t| t[0]),
                        point.transition_lower.map_or(0.0, |t| t[0]),
                        point.iterations.unwrap_or(0),
                        point.residual.unwrap_or(0.0),
                        conv_marker
                    );
                }

                println!();
                println!("Summary:");
                if let Some(cl_max) = polar_output.summary.cl_max {
                    println!(
                        "  CL_max = {:.4} at alpha = {:.2}°",
                        cl_max,
                        polar_output.summary.alpha_at_cl_max.unwrap_or(0.0)
                    );
                }
                if let Some(ld_max) = polar_output.summary.ldratio_max {
                    println!(
                        "  L/D_max = {:.2} at CL = {:.4}",
                        ld_max,
                        polar_output.summary.cl_at_ldratio_max.unwrap_or(0.0)
                    );
                }
                if let Some(cd0) = polar_output.summary.cd0 {
                    println!("  CD0 = {:.6}", cd0);
                }
                println!(
                    "  Converged: {}/{} points",
                    polar_output.summary.n_converged,
                    polar_output.summary.n_converged + polar_output.summary.n_failed
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

            let analysis: AnalysisOutput = match serde_json::from_str(&json_str) {
                Ok(a) => a,
                Err(e) => {
                    eprintln!("Error parsing JSON: {}", e);
                    eprintln!("Make sure the file is an analysis output (from 'yfoil analyse -o')");
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
        "{}: not a geometry, analysis or polar JSON (analysis JSON must come from this version's `yfoil analyse -o`)",
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

/// Write a generated or repanelled geometry as JSON or `.dat`
fn write_geometry(geometry: &Geometry, name: &str, verb: &str, to: &str, output_path: &Path) {
    match to {
        "json" => {
            if let Err(e) = write_geometry_to_json(geometry, output_path) {
                eprintln!("Error writing JSON: {}", e);
                std::process::exit(1);
            }
            println!("{verb} {} with {} nodes", name, geometry.x.len());
            println!("Wrote JSON to {}", output_path.display());
        }
        "dat" => {
            if let Err(e) = write_dat_file(geometry, name, output_path) {
                eprintln!("Error writing DAT: {}", e);
                std::process::exit(1);
            }
            println!("{verb} {} with {} nodes", name, geometry.x.len());
            println!("Wrote DAT to {}", output_path.display());
        }
        _ => {
            eprintln!("Unknown output format: {}", to);
            std::process::exit(1);
        }
    }
}

fn fail_with(context: &str, e: impl std::fmt::Display) -> ! {
    eprintln!("{context}: {e}");
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
            to,
            thickness,
            a,
            panelling,
            output,
        } => {
            let context = "Error generating NACA aerofoil";
            let section = Section::from_designation(&spec).unwrap_or_else(|e| fail_with(context, e));
            let section = match a {
                Some(a) => section.with_a(a).unwrap_or_else(|e| fail_with(context, e)),
                None => section,
            };
            let config = panelling.config(true).unwrap_or_else(|e| fail_with(context, e));
            let geometry = match thickness {
                Thickness::Perpendicular => section.panelled(&config).unwrap_or_else(|e| fail_with(context, e)),
                Thickness::Vertical => {
                    // XFOIL's own model: its 245-point NACA4/NACA5 buffer, always PANGEN-panelled
                    if matches!(config.method, PanelMethod::Cosine(_)) {
                        fail_with(context, "--thickness vertical is XFOIL's NACA4/NACA5 model and is always panelled with PANGEN (--method pangen)");
                    }
                    let buffer = match section.series {
                        Series::FourDigit => naca_4digit_vertical(&spec),
                        Series::FiveDigit => naca_5digit_vertical(&spec),
                        _ => fail_with(
                            context,
                            "--thickness vertical is XFOIL's NACA4/NACA5 model: 4- and 5-digit sections only",
                        ),
                    }
                    .unwrap_or_else(|e| fail_with(context, e));
                    let mut g = repanel(&buffer, &config).unwrap_or_else(|e| fail_with(context, e));
                    if let Some(rec) = g.generator.as_mut() {
                        rec["panelling"]["n_buffer_nodes"] = serde_json::json!(buffer.x.len());
                    }
                    g
                }
            };

            let name = section.designation.clone();
            let file_stem: String = spec
                .to_ascii_lowercase()
                .chars()
                .filter(|c| c.is_ascii_alphanumeric() || *c == '-')
                .collect();
            let output_path = output.unwrap_or_else(|| PathBuf::from(format!("naca{}.{}", file_stem, to)));
            write_geometry(&geometry, &name, "Generated", &to, &output_path);
        }

        GeomAction::KarmanTrefftz {
            x_centre,
            y_centre,
            te_angle,
            to,
            panelling,
            output,
        } => {
            let context = "Error generating Kármán–Trefftz aerofoil";
            let section = KarmanTrefftz::new(x_centre, y_centre, te_angle).unwrap_or_else(|e| fail_with(context, e));
            let config = panelling.config(true).unwrap_or_else(|e| fail_with(context, e));
            let geometry = section.panelled(&config).unwrap_or_else(|e| fail_with(context, e));
            let output_path = output.unwrap_or_else(|| PathBuf::from(format!("karman-trefftz.{}", to)));
            write_geometry(&geometry, &section.designation(), "Generated", &to, &output_path);
        }

        GeomAction::Repanel {
            input,
            panelling,
            output,
        } => {
            let context = "Error repanelling";
            let config = panelling.config(false).unwrap_or_else(|e| fail_with(context, e));
            let geometry = read_geometry_auto(&input);
            let repanelled = repanel(&geometry, &config).unwrap_or_else(|e| fail_with(context, e));

            let output_path = output.unwrap_or_else(|| {
                let mut p = input.clone();
                let stem = p.file_stem().unwrap().to_str().unwrap().to_string();
                p.set_file_name(format!("{stem}_repanelled.json"));
                p
            });
            if let Err(e) = write_geometry_to_json(&repanelled, &output_path) {
                fail_with("Error writing output", e);
            }
            println!(
                "Repanelled from {} to {} nodes with {}",
                geometry.x.len(),
                repanelled.x.len(),
                config.method_name()
            );
            println!("Wrote JSON to {}", output_path.display());
        }

        GeomAction::Info { input, output } => {
            let geometry = read_geometry_auto(&input);
            let aerofoil = panel_foil(&geometry);
            let info = yfoil::output::GeometryInfo::from_panelled(&aerofoil);

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
    println!("  Number of points: {}", summary.n_foil_nodes);
    println!("  Chord: {:.6}", summary.chord);
    println!("  X range: {:.6} to {:.6}", summary.x_range[0], summary.x_range[1]);
    println!("  Y range: {:.6} to {:.6}", summary.y_range[0], summary.y_range[1]);
    println!("  Max thickness: {:.4}", summary.y_extent);
    println!("  TE gap: {:.6}", summary.te_gap);
    println!("  Sharp TE: {}", if summary.sharp_te { "yes" } else { "no" });
    println!(
        "  Reference point: ({:.4}, {:.4})",
        summary.cm_ref[0], summary.cm_ref[1]
    );
    println!("  LE index: {}", summary.i_le_node);
    println!("  LE arc length: {:.6}", summary.s_le);
    println!("  Total arc length: {:.6}", summary.s_total);
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
