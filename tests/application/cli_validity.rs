//! (c) yFoil functionality: the validity flags of `yfoil analyse` and `yfoil polar` —
//! `--conditioning on|off|study`, `--halt-on` and the thresholds — refused when they contradict
//! each other, and the study written with every twin's record.
//!
//! Fixture: the panels of `tests/fixtures/xfoil/naca0012_n60_polar30_re1e6/`.

use crate::fixtures;
use std::process::Command;

fn yfoil() -> Command {
    Command::new(env!("CARGO_BIN_EXE_yfoil"))
}

fn panels() -> String {
    fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar30_re1e6/panels.json")
        .to_string_lossy()
        .into_owned()
}

/// The command fails with a usage error that says `expected`.
fn refused(args: &[&str], expected: &str) {
    let out = yfoil().args(args).output().expect("run yfoil");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "{args:?} should be refused; stderr: {stderr}"
    );
    assert!(
        stderr.contains(expected),
        "{args:?}: stderr `{stderr}` lacks `{expected}`"
    );
}

#[test]
fn test_contradictory_validity_flags_are_refused() {
    let p = panels();
    refused(
        &["polar", &p, "--conditioning", "off", "--halt-on", "ill-conditioned"],
        "--conditioning off",
    );
    refused(
        &[
            "analyse",
            &p,
            "--conditioning",
            "off",
            "--spread-error-threshold",
            "1e-3",
        ],
        "--conditioning off",
    );
    refused(&["polar", &p, "--halt-on", "none,xfoil"], "cannot be combined");
    refused(&["polar", &p, "--halt-on", "unconverged,unconverged=2"], "given twice");
    refused(&["polar", &p, "--halt-on", "unconverged=0"], "at least 1");
    refused(&["polar", &p, "--halt-on", "stalled"], "unknown criterion");
    refused(&["polar", &p, "--halt-on", "anomalous"], "#12");
    refused(&["polar", &p, "--no-twins"], "unexpected argument");
}

#[test]
fn test_analyse_study_writes_every_twin_with_its_numbers() {
    let p = panels();
    let dir = std::env::current_dir().unwrap().join(".tmp/cli_validity");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("study.json");
    let out = yfoil()
        .args(["analyse", &p, "--alpha", "2", "--conditioning", "study", "-o"])
        .arg(&path)
        .output()
        .expect("run yfoil");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let study = v["study"].as_array().expect("a study array");
    assert_eq!(study.len(), 5);
    let id = v["results"]["id"].as_str().unwrap();
    for (twin, seed) in study.iter().zip(1..) {
        assert_eq!(twin["results"]["role"], "twin");
        assert_eq!(twin["results"]["twin_of"], id);
        assert_eq!(twin["results"]["seed"], seed);
        assert!(twin["results"]["values"].is_object() && twin["boundary_layer"].is_object());
        assert_ne!(
            twin["geometry"], v["geometry"],
            "each twin is written on its own jittered panels"
        );
    }
}
