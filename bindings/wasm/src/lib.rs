//! yFoil's WebAssembly binding.
//!
//! Three calls, each taking and returning JSON text in the shapes the CLI reads and writes, so
//! that a result computed in the browser is the file `yfoil` would have written for the same
//! inputs, byte for byte:
//!
//! - [`naca`]: the geometry `yfoil geometry naca <designation> -o x.json` writes;
//! - [`analyse`]: the record `yfoil analyse <geometry> -o x.json` writes;
//! - [`polar`]: the record `yfoil polar <geometry> --json -o x.json` writes.
//!
//! Options are JSON objects whose keys are the CLI's long flags with `_` for `-`; an absent key
//! takes the CLI's default and an unknown key is an error. Each call has a plain-Rust twin
//! (`*_json`) returning `Result<String, String>`, which is what the wasm exports wrap and what
//! native code (the bit-identity check) calls.

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use yfoil::geometry::{panel_foil, Geometry, PanelConfig, Section};
use yfoil::output::{AnalysisOutput, OutputOptions, PointRecord, PolarOutput, Thresholds};
use yfoil::solver::analysis::{
    compute_polar, compute_polar_with, Classifier, FlowConditions, HaltRule, PolarConfig, Session,
};
use yfoil::solver::Twins;

/// The yFoil version this module was built from (the version in every record's point id)
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_owned()
}

/// `yfoil geometry naca`: the panelled section as geometry JSON. `panelling` is a panelling
/// record (the `panelling` object of a generated geometry, `PanelConfig`), or empty for the
/// default (PANGEN, 160 nodes); `a` is the 6-series mean-line loading extent (`--a`), or NaN.
#[wasm_bindgen]
pub fn naca(designation: &str, panelling: &str, a: f64) -> Result<String, JsError> {
    naca_json(designation, panelling, (!a.is_nan()).then_some(a)).map_err(|e| JsError::new(&e))
}

/// `yfoil analyse`: one operating point of `geometry` (geometry JSON), as the analysis record
#[wasm_bindgen]
pub fn analyse(geometry: &str, options: &str) -> Result<String, JsError> {
    analyse_json(geometry, options).map_err(|e| JsError::new(&e))
}

/// `yfoil polar --json`: an alpha sweep of `geometry` (geometry JSON), as the polar record
#[wasm_bindgen]
pub fn polar(geometry: &str, options: &str) -> Result<String, JsError> {
    polar_json(geometry, options).map_err(|e| JsError::new(&e))
}

/// The bit-identity check's probe of the maths library: `function` (a `cr_libm` function name,
/// or `fma`, the fused multiply-add the CORE-MATH sources call) applied elementwise to `x` (and
/// `y`, `z` for the functions that take them). Not an analysis call: it exists so that the
/// module's transcendentals can be compared bit for bit with a native build's.
#[wasm_bindgen]
pub fn evaluate(function: &str, x: &[f64], y: &[f64], z: &[f64]) -> Result<Vec<f64>, JsError> {
    evaluate_native(function, x, y, z).map_err(|e| JsError::new(&e))
}

/// [`evaluate`] without the JavaScript error type
pub fn evaluate_native(function: &str, x: &[f64], y: &[f64], z: &[f64]) -> Result<Vec<f64>, String> {
    let unary = |f: fn(f64) -> f64| Ok(x.iter().map(|&x| f(x)).collect());
    let binary = |f: fn(f64, f64) -> f64| Ok(x.iter().zip(y).map(|(&x, &y)| f(x, y)).collect());
    match function {
        "exp" => unary(cr_libm::exp),
        "ln" => unary(cr_libm::ln),
        "log10" => unary(cr_libm::log10),
        "tanh" => unary(cr_libm::tanh),
        "sin" => unary(cr_libm::sin),
        "cos" => unary(cr_libm::cos),
        "sin_cos.sin" => unary(|x| cr_libm::sin_cos(x).0),
        "sin_cos.cos" => unary(|x| cr_libm::sin_cos(x).1),
        "atan" => unary(cr_libm::atan),
        "asin" => unary(cr_libm::asin),
        "pow" => binary(cr_libm::pow),
        "atan2" => binary(cr_libm::atan2),
        "hypot" => binary(cr_libm::hypot),
        // on wasm32 a call to the `fma` the C sources link against (WebAssembly has no fused
        // multiply-add instruction); on aarch64 the hardware instruction
        "fma" => Ok(x.iter().zip(y).zip(z).map(|((&x, &y), &z)| x.mul_add(y, z)).collect()),
        f => Err(format!("unknown function {f:?}")),
    }
}

/// [`analyse`]'s options: `yfoil analyse`'s flags and defaults
#[derive(Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AnalyseOptions {
    pub alpha: f64,
    pub cl: Option<f64>,
    pub reynolds: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub inviscid: bool,
    pub max_iterations: usize,
    pub include_lagged_closures: bool,
    pub allow_invalid: bool,
    pub conditioning: Conditioning,
    pub spread_warning_threshold: Option<f64>,
    pub spread_error_threshold: Option<f64>,
    /// The foil name in the record (the CLI uses the geometry file's stem)
    pub name: String,
}

impl Default for AnalyseOptions {
    fn default() -> Self {
        Self {
            alpha: 0.0,
            cl: None,
            reynolds: 1_000_000.0,
            mach: 0.0,
            ncrit: 9.0,
            inviscid: false,
            max_iterations: 20,
            include_lagged_closures: false,
            allow_invalid: false,
            conditioning: Conditioning::On,
            spread_warning_threshold: None,
            spread_error_threshold: None,
            name: "Unknown".to_owned(),
        }
    }
}

/// [`polar`]'s options: `yfoil polar`'s flags and defaults
#[derive(Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PolarOptions {
    pub alpha_max: f64,
    pub alpha_min: f64,
    pub alpha_step: f64,
    pub reynolds: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub label: Option<String>,
    pub max_iterations: usize,
    pub distributions: bool,
    pub include_lagged_closures: bool,
    pub allow_invalid: bool,
    pub conditioning: Conditioning,
    pub spread_warning_threshold: Option<f64>,
    pub spread_error_threshold: Option<f64>,
    pub halt_on: String,
    /// The foil name in the record (the CLI uses the geometry file's stem)
    pub name: String,
}

impl Default for PolarOptions {
    fn default() -> Self {
        Self {
            alpha_max: 15.0,
            alpha_min: -5.0,
            alpha_step: 0.5,
            reynolds: 1_000_000.0,
            mach: 0.0,
            ncrit: 9.0,
            label: None,
            max_iterations: 20,
            distributions: false,
            include_lagged_closures: false,
            allow_invalid: false,
            conditioning: Conditioning::On,
            spread_warning_threshold: None,
            spread_error_threshold: None,
            halt_on: "xfoil".to_owned(),
            name: "Unknown".to_owned(),
        }
    }
}

/// `--conditioning`: whether each point is solved with its twins
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Conditioning {
    On,
    Off,
    Study,
}

/// The CLI's [`OutputOptions`] for these flags, and whether the twins run
fn output_options(
    allow_invalid: bool,
    conditioning: Conditioning,
    warning: Option<f64>,
    error: Option<f64>,
) -> Result<(OutputOptions, bool), String> {
    let twins_run = conditioning != Conditioning::Off;
    let thresholds = Thresholds::from_options(warning, error, twins_run)?;
    let options = OutputOptions {
        allow_invalid,
        thresholds,
        study: conditioning == Conditioning::Study,
    };
    Ok((options, twins_run))
}

fn parse_options<T: Default + for<'de> Deserialize<'de>>(options: &str) -> Result<T, String> {
    if options.trim().is_empty() {
        return Ok(T::default());
    }
    serde_json::from_str(options).map_err(|e| format!("options: {e}"))
}

fn parse_geometry(geometry: &str) -> Result<Geometry, String> {
    let geometry: Geometry = serde_json::from_str(geometry).map_err(|e| format!("geometry: {e}"))?;
    geometry.validate().map_err(|e| format!("geometry: {e}"))?;
    Ok(geometry)
}

/// [`naca`] without the JavaScript error type
pub fn naca_json(designation: &str, panelling: &str, a: Option<f64>) -> Result<String, String> {
    let section = Section::from_designation(designation).map_err(|e| e.to_string())?;
    let section = match a {
        Some(a) => section.with_a(a).map_err(|e| e.to_string())?,
        None => section,
    };
    let config = if panelling.trim().is_empty() {
        PanelConfig::default()
    } else {
        let value: serde_json::Value = serde_json::from_str(panelling).map_err(|e| format!("panelling: {e}"))?;
        PanelConfig::from_json(&value).map_err(|e| format!("panelling: {e}"))?
    };
    let geometry = section.panelled(&config).map_err(|e| e.to_string())?;
    // `write_geometry_to_json`'s formatting
    serde_json::to_string_pretty(&geometry).map_err(|e| e.to_string())
}

/// [`analyse`] without the JavaScript error type
pub fn analyse_json(geometry: &str, options: &str) -> Result<String, String> {
    let options: AnalyseOptions = parse_options(options)?;
    let (output_options, twins_run) = output_options(
        options.allow_invalid,
        options.conditioning,
        options.spread_warning_threshold,
        options.spread_error_threshold,
    )?;
    let geometry = parse_geometry(geometry)?;
    let foil = panel_foil(&geometry);
    let twins = twins_run.then(|| Twins::new(&geometry));
    let conditions = FlowConditions {
        re: if options.inviscid { None } else { Some(options.reynolds) },
        mach: options.mach,
        ncrit: options.ncrit,
        max_iterations: options.max_iterations,
        ..FlowConditions::default()
    };
    let mut session = Session::new(&foil, conditions)
        .with_twins(twins)
        .keeping_twin_sessions(output_options.study);
    let point = match options.cl {
        Some(clspec) => session.cl(clspec),
        None => session.alpha(options.alpha.to_radians()),
    };
    AnalysisOutput::from_session(
        &session,
        &point,
        &options.name,
        options.include_lagged_closures,
        &output_options,
    )
    .to_json()
    .map_err(|e| e.to_string())
}

/// [`polar`] without the JavaScript error type
pub fn polar_json(geometry: &str, options: &str) -> Result<String, String> {
    let options: PolarOptions = parse_options(options)?;
    let (output_options, twins_run) = output_options(
        options.allow_invalid,
        options.conditioning,
        options.spread_warning_threshold,
        options.spread_error_threshold,
    )?;
    let halt_on = HaltRule::parse_list(&options.halt_on, twins_run)?;
    // the rules that count a class judge each point as its record will: the same reasons, the
    // same thresholds
    let thresholds = output_options.thresholds;
    let classifier = Classifier(std::sync::Arc::new(move |p| {
        PointRecord::from_point(p, true, false, thresholds).classes
    }));
    let geometry = parse_geometry(geometry)?;
    let foil = panel_foil(&geometry);
    let config = PolarConfig {
        alpha_max: options.alpha_max,
        alpha_min: options.alpha_min,
        alpha_step: options.alpha_step,
        conditions: FlowConditions {
            re: Some(options.reynolds),
            mach: options.mach,
            ncrit: options.ncrit,
            max_iterations: options.max_iterations,
            ..FlowConditions::default()
        },
        twins: twins_run.then(|| Twins::new(&geometry)),
        keep_twin_sessions: output_options.study && options.distributions,
        halt_on,
        classifier: Some(classifier),
    };
    let mut records: Vec<AnalysisOutput> = Vec::new();
    let result = if options.distributions {
        compute_polar_with(&foil, &config, &mut |session, p| {
            records.push(AnalysisOutput::from_session(
                session,
                p,
                &options.name,
                options.include_lagged_closures,
                &output_options,
            ));
        })
    } else {
        compute_polar(&foil, &config)
    };
    // the CLI's ordering (`partial_cmp`): a stable sort that keeps -0 and +0 in sweep order
    records.sort_by(|a, b| {
        a.results
            .alpha_deg
            .partial_cmp(&b.results.alpha_deg)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut output = PolarOutput::from_polar(&result, &options.name, &output_options);
    output.label = options.label;
    output.distributions = records;
    output.to_json().map_err(|e| e.to_string())
}
