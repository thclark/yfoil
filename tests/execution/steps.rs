//! (b) Execution equivalence, single steps of ill-conditioned solutions: for each case below
//! every iteration of the VISCAL calls named in `xtask/fixtures-config/cases.toml` (`step_calls`) is
//! replayed from XFOIL's exact state on entry (`step.rs`), so no step inherits the drift an
//! ill-conditioned solution amplifies over a whole run. Why each case is tested this way is
//! recorded beside it there.
//!
//! Fixtures: `tests/fixtures/xfoil/<case>/` (`mrchdu_input_<k>.dat`, `update_output_<k>.dat`,
//! `viscal_points.dat`) — `cargo xtask fixtures --case <case>`.

use crate::fixtures::cases::load;
use crate::step::replay;

/// Replay every iteration of the case's step calls (`manifest.json`, `case.step_calls`).
fn every_iteration(case: &str) {
    let c = load(case);
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(c.dir.join("manifest.json")).unwrap()).unwrap();
    let calls: Vec<usize> = m["case"]["step_calls"]
        .as_array()
        .expect("step_calls")
        .iter()
        .map(|v| v.as_u64().unwrap() as usize)
        .collect();
    assert!(!calls.is_empty(), "{case}: no step calls");
    let nit: Vec<usize> = c.points.iter().map(|p| p["NITDONE"].parse().unwrap()).collect();
    for call in calls {
        let k0: usize = nit[..call - 1].iter().sum::<usize>() + 1;
        for it in 1..=nit[call - 1] {
            replay(case, call, it, k0 + it - 1);
        }
    }
}

/// `naca0012_n60_sharp_a2_re1e6`: calls [1].
#[test]
fn naca0012_n60_sharp_a2_re1e6_every_iteration() {
    every_iteration("naca0012_n60_sharp_a2_re1e6");
}

/// `naca0012_n60_a12_re1e6`: calls [1].
#[test]
fn naca0012_n60_a12_re1e6_every_iteration() {
    every_iteration("naca0012_n60_a12_re1e6");
}

/// `naca0012_n160_polar_up22_re1e6_iter100`: calls [42].
#[test]
fn naca0012_n160_polar_up22_re1e6_iter100_every_iteration() {
    every_iteration("naca0012_n160_polar_up22_re1e6_iter100");
}

/// `naca4412_n160_polar_down16_re1e6_iter100`: calls [31].
#[test]
fn naca4412_n160_polar_down16_re1e6_iter100_every_iteration() {
    every_iteration("naca4412_n160_polar_down16_re1e6_iter100");
}

/// `naca64a010_n60_type2_a0_re1e6`: calls [1].
#[test]
fn naca64a010_n60_type2_a0_re1e6_every_iteration() {
    every_iteration("naca64a010_n60_type2_a0_re1e6");
}

/// `naca4412_n60_cl1_clm05_re1e6`: calls [1, 2].
#[test]
fn naca4412_n60_cl1_clm05_re1e6_every_iteration() {
    every_iteration("naca4412_n60_cl1_clm05_re1e6");
}

/// `naca4412_n60_a18_re3e6_iter40`: calls [1].
#[test]
fn naca4412_n60_a18_re3e6_iter40_every_iteration() {
    every_iteration("naca4412_n60_a18_re3e6_iter40");
}

/// `naca4412_n60_a20_re1e6_iter50`: calls [1].
#[test]
fn naca4412_n60_a20_re1e6_iter50_every_iteration() {
    every_iteration("naca4412_n60_a20_re1e6_iter50");
}

/// `naca0012_n60_polar30_re1e6`: calls [8, 16, 18].
#[test]
fn naca0012_n60_polar30_re1e6_every_iteration() {
    every_iteration("naca0012_n60_polar30_re1e6");
}

/// `naca4412_n60_polar30_re3e6`: calls [23].
#[test]
fn naca4412_n60_polar30_re3e6_every_iteration() {
    every_iteration("naca4412_n60_polar30_re3e6");
}

/// `naca4412_n160_polar30_re1e6_m05`: calls [37].
#[test]
fn naca4412_n160_polar30_re1e6_m05_every_iteration() {
    every_iteration("naca4412_n160_polar30_re1e6_m05");
}
