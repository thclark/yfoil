//! (f) Reference integrity: BLSOLV's sparse-skip thresholds on the reference case are far from their operands, so the bit-identity gated in `subroutine/blsolv.rs` is not a coincidence of rounding.
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_a2_re1e6/` — `cargo xtask fixtures --case naca0012_n60_a2_re1e6`.

use crate::fixtures;

use fixtures::blsolv_fixtures::{parse_blsolv_input, parse_blsolv_output, XfoilNewtonSystem};
use std::path::PathBuf;
use yfoil::bl::blsolv::{solve_newton_system_traced, BlsolvTrace, NewtonDeltas, NewtonSystem};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// Build the yFoil input exactly as XFOIL's BLSOLV sees it: VZ block enabled, IVTE1/IVZ from
/// IBLSYS, and — critically — `S(N)-S(1)` taken from the fixture, never estimated.
fn yfoil_input(xi: &XfoilNewtonSystem) -> NewtonSystem {
    NewtonSystem {
        n_rows: xi.n_rows,
        diagonal: xi.diagonal.clone(),
        subdiagonal: xi.subdiagonal.clone(),
        rhs: xi.vdel_in.clone(),
        mass_influence: xi.mass_influence.clone(),
        te_block: xi.te_block,
        i_te_row_upper: Some(xi.ivte1_0based()),
        i_wake_row: Some(xi.ivz_0based()),
        elimination_threshold: xi.elimination_threshold,
        s_total: Some(xi.s_total.expect("fixture must carry ARC_LENGTH = S(N)-S(1)")),
    }
}

fn solve_call(
    call: usize,
) -> (
    XfoilNewtonSystem,
    fixtures::blsolv_fixtures::BlsolvOutput,
    NewtonDeltas,
    BlsolvTrace,
) {
    let xi = parse_blsolv_input(&fixture_path("blsolv_input.dat"), call).expect("parse BLSOLV input");
    let xo = parse_blsolv_output(&fixture_path("blsolv_output.dat"), call).expect("parse BLSOLV output");
    assert_eq!(
        xi.n_rows, xo.nsys,
        "call {call}: NSYS mismatch between input and output fixture"
    );
    let mut trace = BlsolvTrace::default();
    let sol = solve_newton_system_traced(yfoil_input(&xi), Some(&mut trace));
    (xi, xo, sol, trace)
}

/// Branch-trace check (plan addendum R3): no sparse-skip comparison `|VTMP| > VACC` may sit
/// within 1e-9 (relative) of its threshold on the reference case. If one ever does, the case is
/// threshold-straddling and must be classified as such, not debugged as a translation error.
#[test]
fn test_blsolv_skip_thresholds_are_well_conditioned() {
    for call in 1..=3 {
        let (xi, _, _, trace) = solve_call(call);
        let taken = trace.skips.iter().filter(|s| s.5).count();
        let (iv, kv, k, margin) = trace.tightest_skip_margin().expect("trace has comparisons");
        println!(
            "call {call}: {} comparisons, {} eliminations taken, tightest margin {margin:+.3e} at (iv={iv}, kv={kv}, k={k}), nsys={}",
            trace.skips.len(),
            taken,
            xi.n_rows
        );
        assert!(
            margin.abs() > 1e-9,
            "call {call}: skip comparison at (iv={iv}, kv={kv}, k={k}) is within {margin:e} of VACC — threshold-straddling case"
        );
    }
}
