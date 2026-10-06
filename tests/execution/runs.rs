//! (b) Execution equivalence, whole runs: each case below takes XFOIL's route through its whole
//! OPER script (or through the call named) with every value within the named tolerances
//! (`run.rs`). How each case is tested, and why, is recorded beside it in
//! `xtask/fixtures-config/cases.toml` (`run`, `run_through`); the polar state machine has its own run
//! test (`polar.rs`).
//!
//! Fixtures: `tests/fixtures/xfoil/<case>/` (`viscal_points.dat`, `viscal_iters_all.dat`) — `cargo
//! xtask fixtures --case <case>`.

use crate::run::compare_run;

/// `naca0012_n60_cl03_re1e6`: the whole run.
#[test]
fn naca0012_n60_cl03_re1e6_run() {
    compare_run("naca0012_n60_cl03_re1e6", 0);
}

/// `naca0012_n60_a2_re1e6_type2`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_type2_run() {
    compare_run("naca0012_n60_a2_re1e6_type2", 0);
}

/// `naca0012_n60_a2_re1e6_m03`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_m03_run() {
    compare_run("naca0012_n60_a2_re1e6_m03", 0);
}

/// `naca0012_n60_a4_re1e5`: the whole run.
#[test]
fn naca0012_n60_a4_re1e5_run() {
    compare_run("naca0012_n60_a4_re1e5", 0);
}

/// `naca0012_n60_a2_re1e6_xtr03`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_xtr03_run() {
    compare_run("naca0012_n60_a2_re1e6_xtr03", 0);
}

/// `naca0012_n60_a2_repeat_re1e6`: the whole run.
#[test]
fn naca0012_n60_a2_repeat_re1e6_run() {
    compare_run("naca0012_n60_a2_repeat_re1e6", 0);
}

/// `naca0012_n60_a2_cl03_re1e6`: the whole run.
#[test]
fn naca0012_n60_a2_cl03_re1e6_run() {
    compare_run("naca0012_n60_a2_cl03_re1e6", 0);
}

/// `naca0012_n60_a2_re1e6_type3`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_type3_run() {
    compare_run("naca0012_n60_a2_re1e6_type3", 0);
}

/// `naca0012_n60_a2_re1e6_m03_type2`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_m03_type2_run() {
    compare_run("naca0012_n60_a2_re1e6_m03_type2", 0);
}

/// `naca0012_n60_a2_re1e6_xtr_coinc`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_xtr_coinc_run() {
    compare_run("naca0012_n60_a2_re1e6_xtr_coinc", 0);
}

/// `naca0012_n60_a2_re1e6_damp`: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_damp_run() {
    compare_run("naca0012_n60_a2_re1e6_damp", 0);
}

/// `naca0012_n160_polar_up22_re1e6_iter100`: calls 1–41.
#[test]
fn naca0012_n160_polar_up22_re1e6_iter100_run() {
    compare_run("naca0012_n160_polar_up22_re1e6_iter100", 41);
}

/// `naca4412_n160_polar_down16_re1e6_iter100`: calls 1–30.
#[test]
fn naca4412_n160_polar_down16_re1e6_iter100_run() {
    compare_run("naca4412_n160_polar_down16_re1e6_iter100", 30);
}

/// `naca4412_n60_inviscid_m07_a10`: every inviscid point.
#[test]
fn naca4412_n60_inviscid_m07_a10_run() {
    compare_run("naca4412_n60_inviscid_m07_a10", 0);
}

/// `kt_n60_inviscid_a4_cl05_cl8`: every inviscid point.
#[test]
fn kt_n60_inviscid_a4_cl05_cl8_run() {
    compare_run("kt_n60_inviscid_a4_cl05_cl8", 0);
}

/// `naca64a010_n60_inviscid_type2_m03`: every inviscid point.
#[test]
fn naca64a010_n60_inviscid_type2_m03_run() {
    compare_run("naca64a010_n60_inviscid_type2_m03", 0);
}

/// `naca23012_n60_a4_re1e6_xtr022_0001`: the whole run.
#[test]
fn naca23012_n60_a4_re1e6_xtr022_0001_run() {
    compare_run("naca23012_n60_a4_re1e6_xtr022_0001", 0);
}

/// `naca0012_n60_a8_re1e4_damp`: the whole run.
#[test]
fn naca0012_n60_a8_re1e4_damp_run() {
    compare_run("naca0012_n60_a8_re1e4_damp", 0);
}

/// `naca0012_n60_polar30_re1e6`: calls 1–7.
#[test]
fn naca0012_n60_polar30_re1e6_run() {
    compare_run("naca0012_n60_polar30_re1e6", 7);
}

/// `naca0012_n60_a2_re1e6`, the reference case: the whole run.
#[test]
fn naca0012_n60_a2_re1e6_run() {
    compare_run("naca0012_n60_a2_re1e6", 0);
}
