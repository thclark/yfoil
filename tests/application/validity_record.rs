//! (c) yFoil functionality: the validity record and the result record built from it: what an unevaluated record
//! asserts, the incompressible margin, and the output contract for an invalid point (status kept,
//! numbers withheld unless `--allow-invalid`). The classification against XFOIL's events is in
//! `known_issues/kt_pole_7_7.rs`.
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_a2_re1e6/` — `cargo xtask fixtures --case naca0012_n60_a2_re1e6`.

use crate::fixtures;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::output::{PointRecord, PointStatus, Reason};
use yfoil::solver::analysis::{FlowConditions, Session};
use yfoil::solver::validity::ValidityRecord;

#[test]
fn test_default_record_asserts_nothing() {
    let v = ValidityRecord::default();
    assert!(
        v.karman_tsien_margin_forces.is_infinite() && v.karman_tsien_margin_pressure.is_infinite(),
        "an unevaluated margin must not read as a violation"
    );
    assert!(!v.out_of_domain());
    assert!(!v.iteration_exhausted());
    assert!(!v.conditions_substituted());
}

#[test]
fn test_incompressible_margin_is_one() {
    // at M = 0, beta = 1 and bfac = 0, so the denominator is 1 at every node whatever the speed:
    // the Kármán–Tsien correction cannot be out of domain in incompressible flow
    let q: Vec<f64> = (0..=40).map(|i| i as f64 * 0.5).collect();
    let (_, margin) = yfoil::solver::clcalc::compute_cp(40, &q, 1.0, 0.0);
    assert_eq!(margin, 1.0);
}

#[test]
fn test_viscous_iteration_exhaustion_withholds_the_numbers() {
    // the reference case converges in 6 iterations; capped at 3 it cannot, and not marginally:
    // XFOIL's own RMSBL after iteration 3 is recorded in viscal_iters_all.dat, far above EPS1
    let dir = fixtures::require_fixture(fixtures::REF_CASE);
    let iters = std::fs::read_to_string(dir.join("viscal_iters_all.dat")).unwrap();
    let rmsbl3: f64 = iters
        .lines()
        .find(|l| l.split_whitespace().take(3).collect::<Vec<_>>() == ["IT", "1", "3"])
        .and_then(|l| l.split_whitespace().nth(3))
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        rmsbl3 > 100.0 * 1e-4,
        "RMSBL after 3 iterations is {rmsbl3}: not near EPS1"
    );

    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let mut session = Session::new(
        &panel_foil(&geometry),
        FlowConditions {
            re: Some(1.0e6),
            max_iterations: 3,
            ..FlowConditions::default()
        },
    );
    let p = session.alpha(2.0_f64.to_radians());
    assert!(!p.converged);
    let r = PointRecord::from_point(&p, true, false);
    assert_eq!(r.status, PointStatus::Invalid);
    assert!(r.reasons.contains(&Reason::ViscousNotConverged), "{:?}", r.reasons);
    assert!(r.values.is_none(), "an invalid point's numbers are withheld");
    let r = PointRecord::from_point(&p, true, true);
    assert_eq!(
        r.status,
        PointStatus::Invalid,
        "--allow-invalid never changes the status"
    );
    assert!(r.values.is_some(), "--allow-invalid populates the numbers");
}
