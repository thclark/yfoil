//! BLSOLV (xsolve.f) against the tracked XFOIL reference fixture.
//!
//! BLSOLV is a pure function of its input file and contains no transcendentals, so with
//! identical inputs and no FMA contraction (the reference is built with -ffp-contract=off)
//! the result is **bit-identical** on any IEEE-754 host. That is the gate here — stronger
//! than the tolerance floor — and it is what the reference fixture shows: 450/450 values
//! identical on each of the three calls, both columns.
//!
//! History: this file once asserted `max_rel_err < 10.0` with a TODO about "~1% errors that
//! accumulate". Those errors came from the test hard-coding `arc_length = 2.0` instead of the
//! fixture's `S(N)-S(1) = 2.0387…`, which shifts the sparse-skip thresholds VACC2/VACC3 by 2%
//! and flips elimination decisions. The solver was never wrong.

mod fixtures;

use fixtures::blsolv_fixtures::{parse_blsolv_input, parse_blsolv_output, BlsolvInput as XfoilBlsolvInput};
use std::path::PathBuf;
use yfoil::bl::blsolv::{blsolv, blsolv_traced, BlsolvInput, BlsolvSolution, BlsolvTrace};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// Build the YFoil input exactly as XFOIL's BLSOLV sees it: VZ block enabled, IVTE1/IVZ from
/// IBLSYS, and — critically — `S(N)-S(1)` taken from the fixture, never estimated.
fn yfoil_input(xi: &XfoilBlsolvInput) -> BlsolvInput {
    BlsolvInput {
        nsys: xi.nsys,
        va: xi.va.clone(),
        vb: xi.vb.clone(),
        vdel: xi.vdel_in.clone(),
        vm: xi.vm.clone(),
        vz: xi.vz,
        ivte1: Some(xi.ivte1_0based()),
        ivz: Some(xi.ivz_0based()),
        vaccel: xi.vaccel,
        arc_length: Some(xi.arc_length.expect("fixture must carry ARC_LENGTH = S(N)-S(1)")),
    }
}

fn solve_call(
    call: usize,
) -> (
    XfoilBlsolvInput,
    fixtures::blsolv_fixtures::BlsolvOutput,
    BlsolvSolution,
    BlsolvTrace,
) {
    let xi = parse_blsolv_input(&fixture_path("blsolv_input.dat"), call).expect("parse BLSOLV input");
    let xo = parse_blsolv_output(&fixture_path("blsolv_output.dat"), call).expect("parse BLSOLV output");
    assert_eq!(
        xi.nsys, xo.nsys,
        "call {call}: NSYS mismatch between input and output fixture"
    );
    let mut trace = BlsolvTrace::default();
    let sol = blsolv_traced(yfoil_input(&xi), Some(&mut trace));
    (xi, xo, sol, trace)
}

/// Per-(equation, column) scale for the scaled error metric: the largest |VDEL| in that slot.
fn scales(xo: &fixtures::blsolv_fixtures::BlsolvOutput) -> [[f64; 2]; 3] {
    let mut sc = [[0.0_f64; 2]; 3];
    for row in &xo.vdel_out {
        for k in 0..3 {
            for c in 0..2 {
                sc[k][c] = sc[k][c].max(row[k][c].abs());
            }
        }
    }
    sc
}

/// The S1 gate: every VDEL value, both columns, all three calls, bit-identical to XFOIL.
#[test]
fn test_blsolv_bit_identical_to_xfoil_calls_1_to_3() {
    for call in 1..=3 {
        let (xi, xo, sol, _) = solve_call(call);
        let sc = scales(&xo);
        let mut mismatches = Vec::new();
        for iv in 0..xi.nsys {
            for k in 0..3 {
                for c in 0..2 {
                    let (a, b) = (sol.vdel[iv][k][c], xo.vdel_out[iv][k][c]);
                    if a.to_bits() != b.to_bits() {
                        let scaled = (a - b).abs() / sc[k][c].max(f64::MIN_POSITIVE);
                        mismatches.push((iv, k, c, a, b, scaled));
                    }
                }
            }
        }
        assert!(
            mismatches.is_empty(),
            "call {call}: {} of {} VDEL values differ from XFOIL; worst scaled error {:.3e} at (iv={}, k={}, col={}): yfoil={:.17e} xfoil={:.17e}",
            mismatches.len(),
            xi.nsys * 6,
            mismatches.iter().map(|m| m.5).fold(0.0, f64::max),
            mismatches[0].0,
            mismatches[0].1,
            mismatches[0].2,
            mismatches[0].3,
            mismatches[0].4
        );
        println!("call {call}: {}/{} values bit-identical", xi.nsys * 6, xi.nsys * 6);
    }
}

/// Branch-trace check (plan addendum R3): no sparse-skip comparison `|VTMP| > VACC` may sit
/// within 1e-9 (relative) of its threshold on the reference case. If one ever does, the case is
/// threshold-straddling and must be classified as such, not debugged as a translation error.
#[test]
fn test_blsolv_no_threshold_straddling() {
    for call in 1..=3 {
        let (xi, _, _, trace) = solve_call(call);
        let taken = trace.skips.iter().filter(|s| s.5).count();
        let (iv, kv, k, margin) = trace.tightest_skip_margin().expect("trace has comparisons");
        println!(
            "call {call}: {} comparisons, {} eliminations taken, tightest margin {margin:+.3e} at (iv={iv}, kv={kv}, k={k}), nsys={}",
            trace.skips.len(),
            taken,
            xi.nsys
        );
        assert!(
            margin.abs() > 1e-9,
            "call {call}: skip comparison at (iv={iv}, kv={kv}, k={k}) is within {margin:e} of VACC — threshold-straddling case"
        );
    }
}

/// Forward-sweep intermediate state (call 1) against XFOIL's `blsolv_trace.dat`
/// "After complete forward sweep" lines: bit-identical before back-substitution too.
#[test]
fn test_blsolv_forward_sweep_matches_xfoil_trace_call_1() {
    let (_, _, _, trace) = solve_call(1);
    let text = std::fs::read_to_string(fixture_path("blsolv_trace.dat")).unwrap();
    let mut checked = 0;
    let mut in_section = false;
    for line in text.lines() {
        if line.contains("After complete forward sweep") {
            in_section = true;
            continue;
        }
        if !in_section {
            continue;
        }
        // VDEL(  iv,3,1)=  value   (1-based iv, mass-defect row, residual column)
        if let Some(rest) = line.strip_prefix("VDEL(") {
            let (idx, val) = rest.split_once(",3,1)=").expect("trace line shape");
            let iv: usize = idx.trim().parse().unwrap();
            let xfoil: f64 = val.trim().parse().unwrap();
            if iv == 0 || iv > trace.vdel_after_forward.len() {
                continue; // stations beyond NSYS are logged as zeros by the fixed-window dump
            }
            let yfoil = trace.vdel_after_forward[iv - 1][2][0];
            assert_eq!(
                yfoil.to_bits(),
                xfoil.to_bits(),
                "forward-sweep VDEL(3,1) at iv={iv}: yfoil={yfoil:.17e} xfoil={xfoil:.17e}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 5, "expected forward-sweep trace values, found {checked}");
    println!("forward sweep: {checked} logged stations bit-identical");
}

/// Sanity: a truncated subsystem solves without NaN/Inf (no reference value; structure only).
#[test]
fn test_blsolv_small_subset_is_finite() {
    let xi = parse_blsolv_input(&fixture_path("blsolv_input.dat"), 1).expect("parse");
    let small = xi.small_subset(10);
    let sol = blsolv(BlsolvInput {
        nsys: small.nsys,
        va: small.va.clone(),
        vb: small.vb.clone(),
        vdel: small.vdel_in.clone(),
        vm: small.vm.clone(),
        vz: [[0.0; 2]; 3],
        ivte1: None,
        ivz: None,
        vaccel: small.vaccel,
        arc_length: small.arc_length,
    });
    for (iv, row) in sol.vdel.iter().enumerate() {
        for (k, r) in row.iter().enumerate() {
            assert!(r[0].is_finite() && r[1].is_finite(), "NaN/Inf at station {iv} row {k}");
        }
    }
}
