//! Generate validation report comparing yFoil vs XFOIL subroutines
//!
//! This generates markdown documentation with detailed comparison tables
//! for BL closure relations (hkin, hs, cf, etc.) and the amplification rate:
//! `docs/fidelity-to-xfoil/subroutine-equivalence/boundary-layer-closures.md` and
//! `transition-and-amplification.md` beside it.
//!
//! Run with: cargo run --bin generate_subroutine_validation

use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use yfoil::bl::system::amplification_rate;
use yfoil::bl::{cdiss_laminar, cf_laminar, cf_turbulent, hk_from_h, hstar_laminar, hstar_turbulent};

/// `TOL_PURE` of `tests/common/utilities/tolerances.rs`, which a binary cannot import
const REL_TOL: f64 = 1e-12;

/// (case, host, libm) of every reference run the fixtures were cut from, from their manifest.
fn reference_provenance(fixture_dir: &Path) -> Vec<(String, String, String)> {
    let text = fs::read_to_string(fixture_dir.join("manifest.json")).unwrap_or_default();
    let doc: serde_json::Value = serde_json::from_str(&text).unwrap_or_default();
    let mut out = vec![];
    for (case, v) in doc["cases"].as_object().into_iter().flatten() {
        let field = |key: &str| -> String {
            v["xfoil_ref"]
                .as_array()
                .into_iter()
                .flatten()
                .filter_map(|l| l.as_str())
                .find_map(|l| l.strip_prefix(&format!("{key}:")).map(|r| r.trim().to_string()))
                .unwrap_or_else(|| "unrecorded".into())
        };
        out.push((case.clone(), field("host"), field("libm")));
    }
    out
}

/// Calculate relative error
fn relative_error(expected: f64, actual: f64) -> f64 {
    if expected.abs() < 1e-15 {
        (expected - actual).abs()
    } else {
        ((expected - actual) / expected).abs()
    }
}

// Fixture structures (same as in tests)
#[derive(Debug, Clone, Serialize, Deserialize)]
struct HkinFixture {
    input: HkinInput,
    output: HkinOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HkinInput {
    h: f64,
    msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HkinOutput {
    hk: f64,
    hk_h: f64,
    hk_msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CflFixture {
    input: CflInput,
    output: CflOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CflInput {
    hk: f64,
    rt: f64,
    msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CflOutput {
    cf: f64,
    cf_hk: f64,
    cf_rt: f64,
    cf_msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HslFixture {
    input: HslInput,
    output: HslOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HslInput {
    hk: f64,
    rt: f64,
    msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HslOutput {
    hs: f64,
    hs_hk: f64,
    hs_rt: f64,
    hs_msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DilFixture {
    input: DilInput,
    output: DilOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DilInput {
    hk: f64,
    rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DilOutput {
    di: f64,
    di_hk: f64,
    di_rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HstFixture {
    input: HstInput,
    output: HstOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HstInput {
    hk: f64,
    rt: f64,
    msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct HstOutput {
    hs: f64,
    hs_hk: f64,
    hs_rt: f64,
    hs_msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CftFixture {
    input: CftInput,
    output: CftOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CftInput {
    hk: f64,
    rt: f64,
    msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CftOutput {
    cf: f64,
    cf_hk: f64,
    cf_rt: f64,
    cf_msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DamplFixture {
    input: DamplInput,
    output: DamplOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DamplInput {
    hk: f64,
    th: f64,
    rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct DamplOutput {
    ax: f64,
    ax_hk: f64,
    ax_th: f64,
    ax_rt: f64,
}

struct ValidationResult {
    name: String,
    test_cases: usize,
    passed: usize,
    failed: usize,
    max_rel_error: f64,
    status: bool,
    details: Vec<CaseDetail>,
}

struct CaseDetail {
    case_num: usize,
    inputs: String,
    outputs: Vec<OutputComparison>,
    #[allow(dead_code)]
    passed: bool,
}

struct OutputComparison {
    name: String,
    xfoil: f64,
    yfoil: f64,
    rel_error: f64,
}

fn load_fixtures<T: for<'de> Deserialize<'de>>(dir: &Path) -> Vec<T> {
    let mut fixtures = Vec::new();
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.filter_map(Result::ok) {
            if entry.path().extension().is_some_and(|e| e == "json") {
                if let Ok(content) = fs::read_to_string(entry.path()) {
                    if let Ok(fixture) = serde_json::from_str::<T>(&content) {
                        fixtures.push(fixture);
                    }
                }
            }
        }
    }
    fixtures
}

fn validate_hkin(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<HkinFixture> = load_fixtures(&fixture_dir.join("hkin"));
    let mut result = ValidationResult {
        name: "HKIN".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let (hk, hk_h, hk_msq) = hk_from_h(f.input.h, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "HK".to_string(),
                xfoil: f.output.hk,
                yfoil: hk,
                rel_error: relative_error(f.output.hk, hk),
            },
            OutputComparison {
                name: "HK_H".to_string(),
                xfoil: f.output.hk_h,
                yfoil: hk_h,
                rel_error: relative_error(f.output.hk_h, hk_h),
            },
            OutputComparison {
                name: "HK_MSQ".to_string(),
                xfoil: f.output.hk_msq,
                yfoil: hk_msq,
                rel_error: relative_error(f.output.hk_msq, hk_msq),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }

        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("H={:.6}, MSQ={:.6}", f.input.h, f.input.msq),
            outputs,
            passed,
        });
    }

    result
}

fn validate_cfl(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<CflFixture> = load_fixtures(&fixture_dir.join("cfl"));
    let mut result = ValidationResult {
        name: "CFL".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = cf_laminar(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "CF".to_string(),
                xfoil: f.output.cf,
                yfoil: r.value,
                rel_error: relative_error(f.output.cf, r.value),
            },
            OutputComparison {
                name: "CF_HK".to_string(),
                xfoil: f.output.cf_hk,
                yfoil: r.value_d_hk,
                rel_error: relative_error(f.output.cf_hk, r.value_d_hk),
            },
            OutputComparison {
                name: "CF_RT".to_string(),
                xfoil: f.output.cf_rt,
                yfoil: r.value_d_retheta,
                rel_error: relative_error(f.output.cf_rt, r.value_d_retheta),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}, RT={:.2}", f.input.hk, f.input.rt),
            outputs,
            passed,
        });
    }

    result
}

fn validate_hsl(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<HslFixture> = load_fixtures(&fixture_dir.join("hsl"));
    let mut result = ValidationResult {
        name: "HSL".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = hstar_laminar(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "HS".to_string(),
                xfoil: f.output.hs,
                yfoil: r.value,
                rel_error: relative_error(f.output.hs, r.value),
            },
            OutputComparison {
                name: "HS_HK".to_string(),
                xfoil: f.output.hs_hk,
                yfoil: r.value_d_hk,
                rel_error: relative_error(f.output.hs_hk, r.value_d_hk),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}", f.input.hk),
            outputs,
            passed,
        });
    }

    result
}

fn validate_dil(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<DilFixture> = load_fixtures(&fixture_dir.join("dil"));
    let mut result = ValidationResult {
        name: "DIL".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = cdiss_laminar(f.input.hk, f.input.rt);

        let outputs = vec![
            OutputComparison {
                name: "DI".to_string(),
                xfoil: f.output.di,
                yfoil: r.value,
                rel_error: relative_error(f.output.di, r.value),
            },
            OutputComparison {
                name: "DI_HK".to_string(),
                xfoil: f.output.di_hk,
                yfoil: r.value_d_hk,
                rel_error: relative_error(f.output.di_hk, r.value_d_hk),
            },
            OutputComparison {
                name: "DI_RT".to_string(),
                xfoil: f.output.di_rt,
                yfoil: r.value_d_retheta,
                rel_error: relative_error(f.output.di_rt, r.value_d_retheta),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}, RT={:.2}", f.input.hk, f.input.rt),
            outputs,
            passed,
        });
    }

    result
}

fn validate_hst(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<HstFixture> = load_fixtures(&fixture_dir.join("hst"));
    let mut result = ValidationResult {
        name: "HST".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = hstar_turbulent(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "HS".to_string(),
                xfoil: f.output.hs,
                yfoil: r.value,
                rel_error: relative_error(f.output.hs, r.value),
            },
            OutputComparison {
                name: "HS_HK".to_string(),
                xfoil: f.output.hs_hk,
                yfoil: r.value_d_hk,
                rel_error: relative_error(f.output.hs_hk, r.value_d_hk),
            },
            OutputComparison {
                name: "HS_RT".to_string(),
                xfoil: f.output.hs_rt,
                yfoil: r.value_d_retheta,
                rel_error: relative_error(f.output.hs_rt, r.value_d_retheta),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}, RT={:.2}", f.input.hk, f.input.rt),
            outputs,
            passed,
        });
    }

    result
}

fn validate_cft(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<CftFixture> = load_fixtures(&fixture_dir.join("cft"));
    let mut result = ValidationResult {
        name: "CFT".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = cf_turbulent(f.input.hk, f.input.rt, f.input.msq, 1.0);

        let outputs = vec![
            OutputComparison {
                name: "CF".to_string(),
                xfoil: f.output.cf,
                yfoil: r.value,
                rel_error: relative_error(f.output.cf, r.value),
            },
            OutputComparison {
                name: "CF_HK".to_string(),
                xfoil: f.output.cf_hk,
                yfoil: r.value_d_hk,
                rel_error: relative_error(f.output.cf_hk, r.value_d_hk),
            },
            OutputComparison {
                name: "CF_RT".to_string(),
                xfoil: f.output.cf_rt,
                yfoil: r.value_d_retheta,
                rel_error: relative_error(f.output.cf_rt, r.value_d_retheta),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}, RT={:.2}", f.input.hk, f.input.rt),
            outputs,
            passed,
        });
    }

    result
}

fn validate_dampl(fixture_dir: &Path) -> ValidationResult {
    let fixtures: Vec<DamplFixture> = load_fixtures(&fixture_dir.join("dampl"));
    let mut result = ValidationResult {
        name: "DAMPL".to_string(),
        test_cases: fixtures.len(),
        passed: 0,
        failed: 0,
        max_rel_error: 0.0,
        status: true,
        details: Vec::new(),
    };

    for (i, f) in fixtures.iter().enumerate() {
        let r = amplification_rate(f.input.hk, f.input.th, f.input.rt);
        let (ax, ax_hk, ax_th, ax_rt) = (r.rate, r.rate_d_hk, r.rate_d_theta, r.rate_d_retheta);

        let outputs = vec![
            OutputComparison {
                name: "AX".to_string(),
                xfoil: f.output.ax,
                yfoil: ax,
                rel_error: relative_error(f.output.ax, ax),
            },
            OutputComparison {
                name: "AX_HK".to_string(),
                xfoil: f.output.ax_hk,
                yfoil: ax_hk,
                rel_error: relative_error(f.output.ax_hk, ax_hk),
            },
            OutputComparison {
                name: "AX_TH".to_string(),
                xfoil: f.output.ax_th,
                yfoil: ax_th,
                rel_error: relative_error(f.output.ax_th, ax_th),
            },
            OutputComparison {
                name: "AX_RT".to_string(),
                xfoil: f.output.ax_rt,
                yfoil: ax_rt,
                rel_error: relative_error(f.output.ax_rt, ax_rt),
            },
        ];

        let max_err = outputs.iter().map(|o| o.rel_error).fold(0.0, f64::max);
        let passed = max_err < REL_TOL;

        if passed {
            result.passed += 1;
        } else {
            result.failed += 1;
            result.status = false;
        }
        result.max_rel_error = result.max_rel_error.max(max_err);

        result.details.push(CaseDetail {
            case_num: i + 1,
            inputs: format!("HK={:.6}, TH={:.2e}, RT={:.2}", f.input.hk, f.input.th, f.input.rt),
            outputs,
            passed,
        });
    }

    result
}

/// The section the reports are written to, from the repository root: the hand-written `index.md`
/// orients the reader, and this generator owns the two pages below
const OUTPUT_DIR: &str = "docs/fidelity-to-xfoil/subroutine-equivalence";
const CLOSURES_PAGE: &str = "boundary-layer-closures.md";
const TRANSITION_PAGE: &str = "transition-and-amplification.md";

/// The fixture directory on GitHub: the site does not serve the repository's files
const FIXTURES_URL: &str = "https://github.com/thclark/yfoil/tree/main/tests/fixtures/subroutines";

fn summary_table(file: &mut File, results: &[&ValidationResult]) -> std::io::Result<()> {
    writeln!(
        file,
        "| Subroutine | Test Cases | Passed | Failed | Max Relative Error | Status |"
    )?;
    writeln!(
        file,
        "|------------|------------|--------|--------|-------------------|--------|"
    )?;
    for r in results {
        let status = if r.status { "✓" } else { "✗" };
        writeln!(
            file,
            "| {} | {} | {} | {} | {:.2e} | {} |",
            r.name, r.test_cases, r.passed, r.failed, r.max_rel_error, status
        )?;
    }
    writeln!(file)
}

fn result_bullets(file: &mut File, r: &ValidationResult) -> std::io::Result<()> {
    let fixture_dir = r.name.to_lowercase();
    writeln!(file, "- **Test Cases:** {}", r.test_cases)?;
    writeln!(file, "- **Passed:** {}", r.passed)?;
    writeln!(file, "- **Max Relative Error:** {:.2e}", r.max_rel_error)?;
    writeln!(file, "- **Status:** {}", if r.status { "✓ PASS" } else { "✗ FAIL" })?;
    writeln!(
        file,
        "- **Fixtures:** [`tests/fixtures/subroutines/{fixture_dir}/`]({FIXTURES_URL}/{fixture_dir}/)"
    )?;
    writeln!(file)
}

fn sample_table(file: &mut File, r: &ValidationResult, n: usize) -> std::io::Result<()> {
    writeln!(file, "### {} sample calls", r.name)?;
    writeln!(file)?;
    writeln!(file, "| Case | Inputs | Output | XFOIL | yFoil | Rel Error |")?;
    writeln!(file, "|------|--------|--------|-------|-------|-----------|")?;
    for detail in r.details.iter().take(n) {
        for (j, output) in detail.outputs.iter().enumerate() {
            if j == 0 {
                writeln!(
                    file,
                    "| {} | {} | {} | {:.10e} | {:.10e} | {:.2e} |",
                    detail.case_num, detail.inputs, output.name, output.xfoil, output.yfoil, output.rel_error
                )?;
            } else {
                writeln!(
                    file,
                    "| | | {} | {:.10e} | {:.10e} | {:.2e} |",
                    output.name, output.xfoil, output.yfoil, output.rel_error
                )?;
            }
        }
    }
    writeln!(file)
}

/// `boundary-layer-closures.md`: the closures at a glance, how they are compared, and each closure's
/// sample calls
fn write_closure_report(results: &[ValidationResult], fixture_dir: &Path, output_dir: &Path) -> std::io::Result<()> {
    let mut file = File::create(output_dir.join(CLOSURES_PAGE))?;
    let closures: Vec<&ValidationResult> = results.iter().filter(|r| r.name != "DAMPL").collect();

    writeln!(file, "# Boundary-layer closures")?;
    writeln!(file)?;
    writeln!(
        file,
        "Generated by `cargo run --bin generate_subroutine_validation` from the tracked fixtures. Do not edit; \
         regenerate."
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "yFoil's boundary-layer closure functions against the calls the instrumented XFOIL reference made to the \
         same functions. The amplification rate, compared the same way, is on \
         [transition and amplification]({TRANSITION_PAGE})."
    )?;
    writeln!(file)?;
    writeln!(file, "## The closures at a glance")?;
    writeln!(file)?;
    summary_table(&mut file, &closures)?;

    writeln!(file, "## What a closure fixture holds")?;
    writeln!(file)?;
    writeln!(
        file,
        "Each fixture is one call of the closure inside a real run of the instrumented double-precision \
         reference ([the reference XFOIL build](../the-reference-xfoil-build/index.md)): its inputs and outputs, \
         dumped at 17 significant digits (`ES24.16`, the minimum that round-trips an f64) while the reference ran \
         the cases below. yFoil's function is called on the same inputs here and its outputs compared with the \
         dumped ones. The error is relative, absolute below 1e-15."
    )?;
    writeln!(file)?;
    writeln!(file, "| Case | Reference build host | libm |")?;
    writeln!(file, "|------|------|------|")?;
    for (case, host, libm) in reference_provenance(fixture_dir) {
        writeln!(file, "| `{case}` | {host} | {libm} |")?;
    }
    writeln!(file)?;
    writeln!(
        file,
        "This report was generated on `{}-{}`. Both codes take their transcendentals (`exp`, `log`, \
         `tanh`, `**`) from the same correctly rounded functions (CORE-MATH, `cr-libm`), whose results do \
         not depend on the host, and the reference is built with `-ffp-contract=off`, so a translation that \
         evaluates each expression in XFOIL's order gives the same bits: an error of exactly 0 means \
         bit-identical, which is what most closures show. A non-zero error means the translation does not \
         evaluate the expression bit for bit as XFOIL does.",
        std::env::consts::ARCH,
        std::env::consts::OS
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "A case passes within **{:.0e}**, the value of `TOL_PURE` in `tests/common/utilities/tolerances.rs`, \
         to which the tests in `tests/subroutine/closures.rs` hold the same fixtures \
         ([tolerances](../tolerances/index.md)).",
        REL_TOL
    )?;
    writeln!(file)?;
    writeln!(file, "## Closure fixtures and how to regenerate them")?;
    writeln!(file)?;
    writeln!(
        file,
        "Fixture data is stored in [`tests/fixtures/subroutines/`]({FIXTURES_URL}/), one directory per closure, \
         one JSON file per call, with the reference's provenance in `manifest.json`. They are cut from the \
         reference's closure log by `cargo xtask fixtures` for the cases with `closures = true` in \
         `xtask/fixtures-config/cases.toml`."
    )?;
    writeln!(file)?;
    writeln!(file, "```bash")?;
    writeln!(file, "# Regenerate this page and the next (uses the tracked fixtures)")?;
    writeln!(file, "cargo run --bin generate_subroutine_validation")?;
    writeln!(file)?;
    writeln!(file, "# To regenerate the fixtures from the instrumented reference:")?;
    writeln!(
        file,
        "cargo xtask fixtures --case closures_naca0012_n60_re1e6 --case closures_naca0012_n60_re1e6_m03"
    )?;
    writeln!(file, "```")?;
    writeln!(file)?;

    for r in closures {
        writeln!(file, "## {}", r.name)?;
        writeln!(file)?;
        result_bullets(&mut file, r)?;
        sample_table(&mut file, r, 10)?;
    }
    Ok(())
}

/// `transition-and-amplification.md`: DAMPL against the reference's calls, and where TRCHEK2's Newton
/// is tested
fn write_transition_report(results: &[ValidationResult], output_dir: &Path) -> std::io::Result<()> {
    let mut file = File::create(output_dir.join(TRANSITION_PAGE))?;

    writeln!(file, "# Transition and amplification")?;
    writeln!(file)?;
    writeln!(
        file,
        "Generated by `cargo run --bin generate_subroutine_validation` from the tracked fixtures. Do not edit; \
         regenerate."
    )?;
    writeln!(file)?;
    writeln!(
        file,
        "Transition in XFOIL is the e^N envelope method: DAMPL gives the spatial amplification rate dN/dx \
         of the most unstable Tollmien–Schlichting wave from the local shape factor and momentum-thickness \
         Reynolds number, and TRCHEK2 integrates it across each interval with a Newton iteration for the point \
         where N reaches Ncrit. DAMPL is compared here call by call, exactly as the closures are \
         ([what a closure fixture holds]({CLOSURES_PAGE}#what-a-closure-fixture-holds)); TRCHEK2 is compared \
         iterate by iterate by `tests/subroutine/transition.rs`, described at the end of the page."
    )?;
    writeln!(file)?;

    for r in results.iter().filter(|r| r.name == "DAMPL") {
        writeln!(file, "## {}", r.name)?;
        writeln!(file)?;
        writeln!(
            file,
            "The DAMPL subroutine computes the spatial amplification rate dN/dx for the eN transition \
             prediction method. Its error, a few ULP, is within `TOL_PURE` and has not been traced to the \
             operation."
        )?;
        writeln!(file)?;
        result_bullets(&mut file, r)?;
        sample_table(&mut file, r, 15)?;
    }

    writeln!(file, "## TRCHEK2's Newton, iterate by iterate")?;
    writeln!(file)?;
    writeln!(
        file,
        "For every dumped SETBL call the instrumented reference writes each TRCHEK2 call's interval inputs and \
         every iterate of its Newton just before the convergence test (`trchek2_<k>.dat`). yFoil's \
         `check_transition_traced` runs on the same inputs: the iterate count must be equal, the initial guess \
         within `TOL_SOLVER` and the final transition point within `TOL_PURE`. The calls are the upper-side \
         transition intervals of the eight iterations of the 7° point of the NACA 0012 polar \
         ([known issues §7.9](../xfoils-known-issues/the-issue-register.md)). The Newton's early exit when \
         there is no amplification over the interval (`xblsys.f:390`) is taken only in a run the reference drives \
         non-finite; every call the reference recorded taking it is replayed and must leave by the same exit at \
         the same iterate ([non-finite probes](../branch-coverage/non-finite-probes.md))."
    )?;
    writeln!(file)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let project_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixture_dir = project_root.join("tests/fixtures/subroutines");
    let output_dir = project_root.join(OUTPUT_DIR);

    fs::create_dir_all(&output_dir)?;

    println!("Running validation...");

    let results = vec![
        validate_hkin(&fixture_dir),
        validate_cfl(&fixture_dir),
        validate_hsl(&fixture_dir),
        validate_dil(&fixture_dir),
        validate_hst(&fixture_dir),
        validate_cft(&fixture_dir),
        validate_dampl(&fixture_dir),
    ];

    println!("\nResults:");
    for r in &results {
        let status = if r.status { "✓" } else { "✗" };
        println!(
            "  {}: {} cases, max error {:.2e} {}",
            r.name, r.test_cases, r.max_rel_error, status
        );
    }

    println!("\nWriting reports to {}", output_dir.display());
    write_closure_report(&results, &fixture_dir, &output_dir)?;
    write_transition_report(&results, &output_dir)?;

    println!("Done!");
    Ok(())
}
