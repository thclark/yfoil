//! (f) Reference integrity: every reference case `xtask/fixtures-config/cases.toml` marks as a run
//! (`run = true`) or as stepped (`step_calls`) has its test in `tests/execution/runs.rs` or
//! `steps.rs`, and every such test names a case marked so — the record of how each case is tested
//! and the tests themselves cannot drift apart.
//!
//! Fixtures: `xtask/fixtures-config/cases.toml` and the two test modules (no reference data).

use std::path::Path;

fn read(rel: &str) -> String {
    std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)).unwrap()
}

/// The case names a module calls `function("<case>", …)` with.
fn called(module: &str, function: &str) -> Vec<String> {
    let pat = format!("{function}(\"");
    read(module)
        .lines()
        .filter_map(|l| {
            l.split_once(&pat)
                .map(|(_, r)| r.split('"').next().unwrap().to_string())
        })
        .collect()
}

#[test]
fn every_run_and_stepped_case_has_its_test() {
    let cases: toml::Value = toml::from_str(&read("xtask/fixtures-config/cases.toml")).unwrap();
    let mut runs = vec![];
    let mut steps = vec![];
    for c in cases["case"].as_array().unwrap() {
        let name = c["name"].as_str().unwrap().to_string();
        if c.get("run").and_then(|v| v.as_bool()) == Some(true) {
            runs.push(name.clone());
        }
        if c.get("step_calls")
            .and_then(|v| v.as_array())
            .is_some_and(|a| !a.is_empty())
        {
            steps.push(name);
        }
    }
    // the polar state machine has its own run test (polar.rs)
    let mut tested_runs = called("tests/execution/runs.rs", "compare_run");
    tested_runs.extend([
        "naca0012_n60_a2_re1e6".to_string(),
        "naca0012_n60_polar_re1e6".to_string(),
    ]);
    let tested_steps = called("tests/execution/steps.rs", "every_iteration");
    runs.sort();
    tested_runs.sort();
    tested_runs.dedup();
    steps.sort();
    let mut ts = tested_steps.clone();
    ts.sort();
    assert_eq!(tested_runs, runs, "run cases (cases.toml) vs tests/execution/runs.rs");
    assert_eq!(ts, steps, "stepped cases (cases.toml) vs tests/execution/steps.rs");
}
