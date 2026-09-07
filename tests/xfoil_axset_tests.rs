//! AXSET function validation tests against XFOIL
//!
//! These tests load fixtures generated from instrumented XFOIL and verify
//! that YFoil's interval_amplification_rate function produces identical output.

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
use yfoil::bl::system::interval_amplification_rate;

const REL_TOL: f64 = 1e-10;
const ABS_TOL: f64 = 1e-14;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct AxsetInput {
    hk1: f64,
    t1: f64,
    rt1: f64,
    a1: f64,
    hk2: f64,
    t2: f64,
    rt2: f64,
    a2: f64,
    acrit: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct AxsetOutput {
    ax: f64,
    ax_hk1: f64,
    ax_t1: f64,
    ax_rt1: f64,
    ax_a1: f64,
    ax_hk2: f64,
    ax_t2: f64,
    ax_rt2: f64,
    ax_a2: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct AxsetFixture {
    input: AxsetInput,
    output: AxsetOutput,
}

fn relative_error(expected: f64, actual: f64) -> f64 {
    if expected.abs() < ABS_TOL {
        (expected - actual).abs()
    } else {
        ((expected - actual) / expected).abs()
    }
}

fn check_value(name: &str, case_num: usize, expected: f64, actual: f64) {
    let rel_err = relative_error(expected, actual);
    let abs_err = (expected - actual).abs();

    if rel_err > REL_TOL && abs_err > ABS_TOL {
        panic!(
            "Case {}: {} mismatch - XFOIL={:.16e}, YFoil={:.16e}, rel_err={:.2e}",
            case_num, name, expected, actual, rel_err
        );
    }
}

fn load_fixtures(dir: &Path) -> Vec<AxsetFixture> {
    let mut fixtures = Vec::new();
    if let Ok(entries) = fs::read_dir(dir) {
        let mut paths: Vec<_> = entries.filter_map(Result::ok).collect();
        paths.sort_by_key(|e| e.path());

        for entry in paths {
            if entry.path().extension().is_some_and(|e| e == "json") {
                if let Ok(content) = fs::read_to_string(entry.path()) {
                    if let Ok(fixture) = serde_json::from_str::<AxsetFixture>(&content) {
                        fixtures.push(fixture);
                    }
                }
            }
        }
    }
    fixtures
}

#[test]
fn test_axset_against_xfoil_fixtures() {
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/subroutines/axset");

    let fixtures = load_fixtures(&fixture_dir);
    assert!(!fixtures.is_empty(), "No AXSET fixtures found in {:?}", fixture_dir);

    let mut passed = 0;
    let mut failed = 0;
    let mut max_rel_error = 0.0f64;

    for (i, f) in fixtures.iter().enumerate() {
        let case_num = i + 1;
        let inp = &f.input;
        let out = &f.output;

        let result = interval_amplification_rate(
            inp.hk1, inp.t1, inp.rt1, inp.a1, inp.hk2, inp.t2, inp.rt2, inp.a2, inp.acrit, 0,
        );

        // Check each output value
        let checks = [
            ("AX", out.ax, result.rate),
            ("AX_HK1", out.ax_hk1, result.rate_d_hk_station1),
            ("AX_T1", out.ax_t1, result.rate_d_theta_station1),
            ("AX_RT1", out.ax_rt1, result.rate_d_retheta_station1),
            ("AX_A1", out.ax_a1, result.rate_d_ampl_station1),
            ("AX_HK2", out.ax_hk2, result.rate_d_hk_station2),
            ("AX_T2", out.ax_t2, result.rate_d_theta_station2),
            ("AX_RT2", out.ax_rt2, result.rate_d_retheta_station2),
            ("AX_A2", out.ax_a2, result.rate_d_ampl_station2),
        ];

        let mut case_passed = true;
        for (name, expected, actual) in checks {
            let rel_err = relative_error(expected, actual);
            let abs_err = (expected - actual).abs();

            max_rel_error = max_rel_error.max(rel_err);

            if rel_err > REL_TOL && abs_err > ABS_TOL {
                case_passed = false;
                eprintln!(
                    "Case {}: {} FAIL - XFOIL={:.16e}, YFoil={:.16e}, rel_err={:.2e}",
                    case_num, name, expected, actual, rel_err
                );
            }
        }

        if case_passed {
            passed += 1;
        } else {
            failed += 1;
        }
    }

    println!(
        "\nAXSET validation: {} passed, {} failed, max_rel_error={:.2e}",
        passed, failed, max_rel_error
    );

    assert_eq!(failed, 0, "{} AXSET test cases failed", failed);
}

#[test]
fn test_axset_sample_cases() {
    // Test a few specific cases with detailed assertions
    let fixture_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/subroutines/axset");

    let fixtures = load_fixtures(&fixture_dir);

    // Test first case
    if let Some(f) = fixtures.first() {
        let inp = &f.input;
        let out = &f.output;

        let result = interval_amplification_rate(
            inp.hk1, inp.t1, inp.rt1, inp.a1, inp.hk2, inp.t2, inp.rt2, inp.a2, inp.acrit, 0,
        );

        check_value("AX", 1, out.ax, result.rate);
        check_value("AX_HK1", 1, out.ax_hk1, result.rate_d_hk_station1);
        check_value("AX_T1", 1, out.ax_t1, result.rate_d_theta_station1);
        check_value("AX_RT1", 1, out.ax_rt1, result.rate_d_retheta_station1);
        check_value("AX_A1", 1, out.ax_a1, result.rate_d_ampl_station1);
        check_value("AX_HK2", 1, out.ax_hk2, result.rate_d_hk_station2);
        check_value("AX_T2", 1, out.ax_t2, result.rate_d_theta_station2);
        check_value("AX_RT2", 1, out.ax_rt2, result.rate_d_retheta_station2);
        check_value("AX_A2", 1, out.ax_a2, result.rate_d_ampl_station2);
    }

    // Test last case to cover different parameter ranges
    if let Some(f) = fixtures.last() {
        let idx = fixtures.len();
        let inp = &f.input;
        let out = &f.output;

        let result = interval_amplification_rate(
            inp.hk1, inp.t1, inp.rt1, inp.a1, inp.hk2, inp.t2, inp.rt2, inp.a2, inp.acrit, 0,
        );

        check_value("AX", idx, out.ax, result.rate);
        check_value("AX_A1", idx, out.ax_a1, result.rate_d_ampl_station1);
        check_value("AX_A2", idx, out.ax_a2, result.rate_d_ampl_station2);
    }
}
