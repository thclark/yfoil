//! Generate validation report comparing YFoil vs XFOIL subroutines
//!
//! This generates markdown documentation with detailed comparison tables
//! for BL closure relations (hkin, hs, cf, etc.).
//!
//! Run with: cargo run --bin generate_subroutine_validation

use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use yfoil::bl::{cf_lam, cf_turb, dampl, di_lam, hkin, hs_lam, hs_turb};

const REL_TOL: f64 = 1e-12;

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
        let (hk, hk_h, hk_msq) = hkin(f.input.h, f.input.msq);

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
        let r = cf_lam(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "CF".to_string(),
                xfoil: f.output.cf,
                yfoil: r.val,
                rel_error: relative_error(f.output.cf, r.val),
            },
            OutputComparison {
                name: "CF_HK".to_string(),
                xfoil: f.output.cf_hk,
                yfoil: r.val_hk,
                rel_error: relative_error(f.output.cf_hk, r.val_hk),
            },
            OutputComparison {
                name: "CF_RT".to_string(),
                xfoil: f.output.cf_rt,
                yfoil: r.val_rt,
                rel_error: relative_error(f.output.cf_rt, r.val_rt),
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
        let r = hs_lam(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "HS".to_string(),
                xfoil: f.output.hs,
                yfoil: r.val,
                rel_error: relative_error(f.output.hs, r.val),
            },
            OutputComparison {
                name: "HS_HK".to_string(),
                xfoil: f.output.hs_hk,
                yfoil: r.val_hk,
                rel_error: relative_error(f.output.hs_hk, r.val_hk),
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
        let r = di_lam(f.input.hk, f.input.rt);

        let outputs = vec![
            OutputComparison {
                name: "DI".to_string(),
                xfoil: f.output.di,
                yfoil: r.val,
                rel_error: relative_error(f.output.di, r.val),
            },
            OutputComparison {
                name: "DI_HK".to_string(),
                xfoil: f.output.di_hk,
                yfoil: r.val_hk,
                rel_error: relative_error(f.output.di_hk, r.val_hk),
            },
            OutputComparison {
                name: "DI_RT".to_string(),
                xfoil: f.output.di_rt,
                yfoil: r.val_rt,
                rel_error: relative_error(f.output.di_rt, r.val_rt),
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
        let r = hs_turb(f.input.hk, f.input.rt, f.input.msq);

        let outputs = vec![
            OutputComparison {
                name: "HS".to_string(),
                xfoil: f.output.hs,
                yfoil: r.val,
                rel_error: relative_error(f.output.hs, r.val),
            },
            OutputComparison {
                name: "HS_HK".to_string(),
                xfoil: f.output.hs_hk,
                yfoil: r.val_hk,
                rel_error: relative_error(f.output.hs_hk, r.val_hk),
            },
            OutputComparison {
                name: "HS_RT".to_string(),
                xfoil: f.output.hs_rt,
                yfoil: r.val_rt,
                rel_error: relative_error(f.output.hs_rt, r.val_rt),
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
        let r = cf_turb(f.input.hk, f.input.rt, f.input.msq, 1.0);

        let outputs = vec![
            OutputComparison {
                name: "CF".to_string(),
                xfoil: f.output.cf,
                yfoil: r.val,
                rel_error: relative_error(f.output.cf, r.val),
            },
            OutputComparison {
                name: "CF_HK".to_string(),
                xfoil: f.output.cf_hk,
                yfoil: r.val_hk,
                rel_error: relative_error(f.output.cf_hk, r.val_hk),
            },
            OutputComparison {
                name: "CF_RT".to_string(),
                xfoil: f.output.cf_rt,
                yfoil: r.val_rt,
                rel_error: relative_error(f.output.cf_rt, r.val_rt),
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
        let (ax, ax_hk, ax_th, ax_rt) = dampl(f.input.hk, f.input.th, f.input.rt);

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

fn write_summary_report(results: &[ValidationResult], output_dir: &Path) -> std::io::Result<()> {
    let mut file = File::create(output_dir.join("README.md"))?;

    writeln!(file, "# Subroutine Validation")?;
    writeln!(file)?;
    writeln!(
        file,
        "This report validates YFoil's boundary layer closure functions against XFOIL reference values."
    )?;
    writeln!(file)?;
    writeln!(file, "## Summary")?;
    writeln!(file)?;
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

    writeln!(file)?;
    writeln!(file, "## Tolerance")?;
    writeln!(file)?;
    writeln!(file, "All tests use a relative tolerance of **{:.0e}**, ensuring at least 12 significant figures match between YFoil and XFOIL.", REL_TOL)?;
    writeln!(file)?;
    writeln!(file, "## Test Fixtures")?;
    writeln!(file)?;
    writeln!(
        file,
        "Fixture data is stored in [`tests/fixtures/subroutines/`](../../../tests/fixtures/subroutines/):"
    )?;
    writeln!(file)?;
    writeln!(file, "```")?;
    writeln!(file, "tests/fixtures/subroutines/")?;
    writeln!(file, "├── hkin/     - Kinematic shape factor (HKIN)")?;
    writeln!(file, "├── hsl/      - Laminar energy shape factor (HSL)")?;
    writeln!(file, "├── hst/      - Turbulent energy shape factor (HST)")?;
    writeln!(file, "├── cfl/      - Laminar skin friction (CFL)")?;
    writeln!(file, "├── cft/      - Turbulent skin friction (CFT)")?;
    writeln!(file, "├── dil/      - Laminar dissipation (DIL)")?;
    writeln!(file, "└── dampl/    - Amplification rate (DAMPL)")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(
        file,
        "Each fixture is a JSON file containing input parameters and expected XFOIL output."
    )?;
    writeln!(file)?;
    writeln!(file, "## XFOIL Instrumentation")?;
    writeln!(file)?;
    writeln!(
        file,
        "Test fixtures were generated by instrumenting XFOIL with WRITE statements to log"
    )?;
    writeln!(
        file,
        "subroutine inputs and outputs. The instrumentation was added to `xfoil/xfoil6.99/src/xbl.f`."
    )?;
    writeln!(file)?;
    writeln!(file, "**Example instrumentation (HKIN subroutine):**")?;
    writeln!(file)?;
    writeln!(file, "```fortran")?;
    writeln!(file, "C---- Log inputs and outputs for validation")?;
    writeln!(file, "      WRITE(96,'(A)') '=== HKIN ==='")?;
    writeln!(file, "      WRITE(96,'(A,E24.16)') 'H=', H")?;
    writeln!(file, "      WRITE(96,'(A,E24.16)') 'MSQ=', MSQ")?;
    writeln!(file, "      WRITE(96,'(A,E24.16)') 'HK=', HK")?;
    writeln!(file, "      WRITE(96,'(A,E24.16)') 'HK_H=', HK_H")?;
    writeln!(file, "      WRITE(96,'(A,E24.16)') 'HK_MSQ=', HK_MSQ")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(
        file,
        "Similar instrumentation was added to each closure subroutine (HSL, HST, CFL, CFT, DIL, DAMPL)."
    )?;
    writeln!(
        file,
        "The instrumented XFOIL writes to unit 96, which outputs `xfoil_subroutine_log.dat`."
    )?;
    writeln!(file)?;
    writeln!(file, "**XFOIL script used to generate fixture data:**")?;
    writeln!(file)?;
    writeln!(file, "```")?;
    writeln!(file, "PLOP")?;
    writeln!(file, "G F")?;
    writeln!(file)?;
    writeln!(file, "NACA 0012")?;
    writeln!(file, "OPER")?;
    writeln!(file, "VISC 1e6")?;
    writeln!(file, "ITER 50")?;
    writeln!(file, "ALFA 0")?;
    writeln!(file, "ALFA 2")?;
    writeln!(file, "ALFA 5")?;
    writeln!(file, "ALFA 8")?;
    writeln!(file, "ALFA 10")?;
    writeln!(file)?;
    writeln!(file, "QUIT")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(
        file,
        "Running this script with instrumented XFOIL produces the log file, which is then"
    )?;
    writeln!(
        file,
        "parsed by `examples/generate_subroutine_fixtures.rs` to create JSON fixtures."
    )?;
    writeln!(file)?;
    writeln!(file, "## Regenerating Validation")?;
    writeln!(file)?;
    writeln!(file, "```bash")?;
    writeln!(file, "# Regenerate validation report (uses existing fixtures)")?;
    writeln!(file, "cargo run --bin generate_subroutine_validation")?;
    writeln!(file)?;
    writeln!(file, "# To regenerate fixtures from instrumented XFOIL:")?;
    writeln!(file, "cargo run --example generate_subroutine_fixtures")?;
    writeln!(file, "```")?;
    writeln!(file)?;
    writeln!(file, "## Detailed Reports")?;
    writeln!(file)?;
    writeln!(
        file,
        "- [Closure Functions (HKIN, CFL, HSL, DIL, CFT, HST)](closure.md)"
    )?;
    writeln!(file, "- [Transition (DAMPL)](transition.md)")?;

    Ok(())
}

fn write_closure_report(results: &[ValidationResult], output_dir: &Path) -> std::io::Result<()> {
    let mut file = File::create(output_dir.join("closure.md"))?;

    writeln!(file, "# Closure Functions Validation")?;
    writeln!(file)?;
    writeln!(
        file,
        "This report details the validation of YFoil's boundary layer closure relations."
    )?;
    writeln!(file)?;

    for r in results {
        if r.name == "DAMPL" {
            continue;
        } // DAMPL goes in transition.md

        let fixture_dir = r.name.to_lowercase();
        writeln!(file, "## {}", r.name)?;
        writeln!(file)?;
        writeln!(file, "- **Test Cases:** {}", r.test_cases)?;
        writeln!(file, "- **Passed:** {}", r.passed)?;
        writeln!(file, "- **Max Relative Error:** {:.2e}", r.max_rel_error)?;
        writeln!(file, "- **Status:** {}", if r.status { "✓ PASS" } else { "✗ FAIL" })?;
        writeln!(
            file,
            "- **Fixtures:** [`tests/fixtures/subroutines/{}/`](../../../tests/fixtures/subroutines/{}/)",
            fixture_dir, fixture_dir
        )?;
        writeln!(file)?;

        // Show sample of cases (first 10)
        writeln!(file, "### Sample Cases")?;
        writeln!(file)?;
        writeln!(file, "| Case | Inputs | Output | XFOIL | YFoil | Rel Error |")?;
        writeln!(file, "|------|--------|--------|-------|-------|-----------|")?;

        for detail in r.details.iter().take(10) {
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
        writeln!(file)?;
    }

    Ok(())
}

fn write_transition_report(results: &[ValidationResult], output_dir: &Path) -> std::io::Result<()> {
    let mut file = File::create(output_dir.join("transition.md"))?;

    writeln!(file, "# Transition (DAMPL) Validation")?;
    writeln!(file)?;
    writeln!(
        file,
        "This report details the validation of YFoil's amplification rate calculation (DAMPL)."
    )?;
    writeln!(file)?;

    for r in results {
        if r.name != "DAMPL" {
            continue;
        }

        writeln!(file, "## {}", r.name)?;
        writeln!(file)?;
        writeln!(
            file,
            "The DAMPL subroutine computes the spatial amplification rate dN/dx for"
        )?;
        writeln!(file, "the eN transition prediction method.")?;
        writeln!(file)?;
        writeln!(file, "- **Test Cases:** {}", r.test_cases)?;
        writeln!(file, "- **Passed:** {}", r.passed)?;
        writeln!(file, "- **Max Relative Error:** {:.2e}", r.max_rel_error)?;
        writeln!(file, "- **Status:** {}", if r.status { "✓ PASS" } else { "✗ FAIL" })?;
        writeln!(
            file,
            "- **Fixtures:** [`tests/fixtures/subroutines/dampl/`](../../../tests/fixtures/subroutines/dampl/)"
        )?;
        writeln!(file)?;

        writeln!(file, "### Sample Cases")?;
        writeln!(file)?;
        writeln!(file, "| Case | Inputs | Output | XFOIL | YFoil | Rel Error |")?;
        writeln!(file, "|------|--------|--------|-------|-------|-----------|")?;

        for detail in r.details.iter().take(15) {
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
        writeln!(file)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let project_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixture_dir = project_root.join("tests/fixtures/subroutines");
    let output_dir = project_root.join("docs/validation/subroutines");

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
    write_summary_report(&results, &output_dir)?;
    write_closure_report(&results, &output_dir)?;
    write_transition_report(&results, &output_dir)?;

    println!("Done!");
    Ok(())
}
