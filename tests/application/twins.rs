//! (c) yFoil functionality: the solver's twins (`yfoil::solver::twins`). A twin is the same solve
//! from the same initial state on jittered panels; an identity twin (the base's own panels) must
//! therefore reproduce the base bit for bit at every point — fresh march, inherited boundary
//! layer, fixed CL and inviscid — which proves a twin session carries the whole state it is
//! copied from. The jitter is the one `cargo xtask twins` applies to the reference.
//!
//! Fixtures: the panels of `tests/fixtures/xfoil/naca0012_n60_polar30_re1e6/` and
//! `tests/fixtures/xfoil/kt_n60_inviscid_a4_cl05_cl8/` — `cargo xtask fixtures --case NAME`.

use crate::fixtures;
use yfoil::geometry::jitter::{jitter, TWIN_SEEDS};
use yfoil::geometry::{panel_foil, read_geometry_from_file, Geometry};
use yfoil::output::{
    PointRecord, PointStatus, Reason, Thresholds, Warning, CONDITIONING_ERROR_THRESHOLD,
    CONDITIONING_WARNING_THRESHOLD,
};
use yfoil::solver::analysis::{compute_polar, FlowConditions, PointResult, PolarConfig, Session};
use yfoil::solver::Twins;

fn geometry(case: &str) -> Geometry {
    let path = fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}/panels.json"));
    read_geometry_from_file(path.to_str().unwrap()).expect("panels.json")
}

/// Every number a point reports, as bits, and its exact decisions.
fn fingerprint(p: &PointResult) -> Vec<u64> {
    let mut v: Vec<u64> = [
        p.alpha,
        p.cl,
        p.cd,
        p.cd_friction,
        p.cd_pressure,
        p.cm,
        p.cl_d_alpha,
        p.transition_upper[0],
        p.transition_upper[1],
        p.transition_lower[0],
        p.transition_lower[1],
        p.residual,
    ]
    .iter()
    .map(|x| x.to_bits())
    .collect();
    v.extend([
        p.converged as u64,
        p.iterations as u64,
        p.inviscid_cl_iterations as u64,
        p.nonfinite_station_failures as u64,
        p.i_transition_station[1] as u64,
        p.i_transition_station[2] as u64,
    ]);
    v
}

fn assert_identity_twin_reproduces(p: &PointResult) {
    let twins = p.twins.as_ref().expect("twins were run");
    assert_eq!(twins.len(), 1);
    assert_eq!(
        fingerprint(&twins[0].1),
        fingerprint(p),
        "identity twin of the point at alpha {:.3}°",
        p.alpha.to_degrees()
    );
}

#[test]
fn test_identity_twin_reproduces_every_point_of_a_polar() {
    // ±30° through stall: fresh marches at 0° and inherited boundary layers everywhere else
    let g = geometry("naca0012_n60_polar30_re1e6");
    let config = PolarConfig {
        alpha_max: 30.0,
        alpha_min: -30.0,
        alpha_step: 1.0,
        conditions: FlowConditions::default(),
        twins: Some(Twins::identity(&g)),
        ..PolarConfig::default()
    };
    let polar = compute_polar(&panel_foil(&g), &config);
    assert!(polar.results.len() > 30);
    for p in &polar.results {
        assert_identity_twin_reproduces(p);
    }
}

#[test]
fn test_identity_twin_reproduces_fixed_cl_and_inviscid_points() {
    let g = geometry("naca0012_n60_polar30_re1e6");
    let mut session = Session::new(&panel_foil(&g), FlowConditions::default()).with_twins(Some(Twins::identity(&g)));
    assert_identity_twin_reproduces(&session.alpha(2.0_f64.to_radians()));
    assert_identity_twin_reproduces(&session.cl(0.3));
    assert_identity_twin_reproduces(&session.cl(0.5));

    let g = geometry("kt_n60_inviscid_a4_cl05_cl8");
    let inviscid = FlowConditions {
        re: None,
        ..FlowConditions::default()
    };
    let mut session = Session::new(&panel_foil(&g), inviscid).with_twins(Some(Twins::identity(&g)));
    assert_identity_twin_reproduces(&session.alpha(4.0_f64.to_radians()));
    assert_identity_twin_reproduces(&session.cl(0.8));
}

#[test]
fn test_twins_carry_the_five_seeds_and_do_not_steer_the_base() {
    let g = geometry("naca0012_n60_polar30_re1e6");
    let foil = panel_foil(&g);
    let alone = Session::new(&foil, FlowConditions::default()).alpha(2.0_f64.to_radians());
    let with = Session::new(&foil, FlowConditions::default())
        .with_twins(Some(Twins::new(&g)))
        .alpha(2.0_f64.to_radians());
    assert_eq!(fingerprint(&alone), fingerprint(&with), "the base solve is the same with or without twins");
    assert_eq!(with.id, alone.id, "twins are not part of what the point is");
    let seeds: Vec<u64> = with.twins.as_ref().unwrap().iter().map(|(s, _)| *s).collect();
    assert_eq!(seeds, TWIN_SEEDS);
    for (_, t) in with.twins.as_ref().unwrap() {
        assert_ne!(fingerprint(t), fingerprint(&with), "a jittered twin is not the base bit for bit");
    }
}

/// A converged viscous point with these twins attached (the record logic needs nothing else).
fn point_with_twins(twins: Vec<PointResult>) -> PointResult {
    let base = PointResult {
        alpha: 0.1,
        cl: 0.8,
        cd: 0.01,
        cd_friction: 0.006,
        cd_pressure: 0.004,
        cm: -0.05,
        transition_upper: [0.3, 0.05],
        transition_lower: [0.9, -0.01],
        converged: true,
        iterations: 6,
        id: "base".into(),
        ..PointResult::default()
    };
    PointResult {
        twins: Some(twins.into_iter().enumerate().map(|(k, t)| (k as u64 + 1, t)).collect()),
        ..base
    }
}

fn twin_of(p: &PointResult, change: impl Fn(&mut PointResult)) -> PointResult {
    let mut t = PointResult { twins: None, ..p.clone() };
    change(&mut t);
    t
}

#[test]
fn test_twins_that_agree_leave_the_point_valid_with_its_evidence() {
    let base = point_with_twins(vec![]);
    let p = point_with_twins((0..5).map(|_| twin_of(&base, |t| t.cl += 1e-12)).collect());
    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.status, PointStatus::Valid);
    assert!(r.values.is_some());
    let c = r.diagnostics.unwrap().conditioning.expect("measured");
    assert_eq!(c.seeds, [1, 2, 3, 4, 5]);
    assert_eq!(c.worst, "cl");
    assert!(c.worst_spread() > 0.0 && c.worst_spread() < 1e-11);
    assert_eq!(c.twins_not_valid, 0);
    assert!(r.warnings.is_empty());
}

#[test]
fn test_twins_that_disagree_make_the_point_ill_conditioned_and_withhold_it() {
    let base = point_with_twins(vec![]);
    // one twin of five reaches another solution: CD moves by three times the threshold
    let mut twins: Vec<_> = (0..4).map(|_| twin_of(&base, |_| {})).collect();
    twins.push(twin_of(&base, |t| t.cd += 3.0 * CONDITIONING_ERROR_THRESHOLD));
    let p = point_with_twins(twins);

    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.status, PointStatus::Invalid);
    assert_eq!(r.reasons, [Reason::IllConditioned]);
    assert!(r.values.is_none(), "the numbers are withheld");
    let c = r.diagnostics.as_ref().unwrap().conditioning.as_ref().unwrap();
    assert_eq!(c.worst, "cd");
    assert!((c.worst_spread() - 3.0 * CONDITIONING_ERROR_THRESHOLD).abs() < 1e-12);
    assert_eq!(c.thresholds, Thresholds::default());

    // --allow-invalid shows them and changes nothing else
    let shown = PointRecord::from_point(&p, true, true, Thresholds::default());
    assert_eq!(shown.status, PointStatus::Invalid);
    assert_eq!(shown.reasons, [Reason::IllConditioned]);
    assert_eq!(shown.values.unwrap().cl, 0.8);

    // a spread below the error threshold is not ill-conditioned
    let at = point_with_twins(vec![twin_of(&base, |t| t.cl += CONDITIONING_ERROR_THRESHOLD * 0.8)]);
    assert_eq!(PointRecord::from_point(&at, true, false, Thresholds::default()).status, PointStatus::Valid);
}

#[test]
fn test_a_spread_between_the_thresholds_is_a_warning_and_the_thresholds_are_settable() {
    let base = point_with_twins(vec![]);
    let spread = 3.0 * CONDITIONING_WARNING_THRESHOLD;
    let p = point_with_twins(vec![twin_of(&base, |_| {}), twin_of(&base, |t| t.cm += spread)]);

    // by default: valid, numbers shown, and a warning naming the quantity and its spread
    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.status, PointStatus::Valid);
    assert!(r.values.is_some());
    let [Warning::IllConditioned { quantity, spread: s }] = r.warnings.as_slice() else {
        panic!("one IllConditioned warning expected: {:?}", r.warnings);
    };
    assert_eq!(quantity, "cm");
    assert!((s - spread).abs() < 1e-15);
    assert_eq!(
        r.warnings[0].to_string(),
        "solution may be ill-conditioned: its twins spread cm by 3.00e-8"
    );

    // a spread at or below the warning threshold says nothing
    let quiet = point_with_twins(vec![twin_of(&base, |t| t.cm += 0.5 * CONDITIONING_WARNING_THRESHOLD)]);
    assert!(PointRecord::from_point(&quiet, true, false, Thresholds::default()).warnings.is_empty());

    // a lower error threshold withholds the same point, with no warning beside the reason
    let strict = Thresholds { warning: 1e-9, error: 1e-8 };
    let r = PointRecord::from_point(&p, true, false, strict);
    assert_eq!(r.reasons, [Reason::IllConditioned]);
    assert!(r.warnings.is_empty());
    assert_eq!(r.diagnostics.unwrap().conditioning.unwrap().thresholds, strict);

    // a higher warning threshold silences it
    let lax = Thresholds { warning: 1e-6, error: 1e-4 };
    let r = PointRecord::from_point(&p, true, false, lax);
    assert_eq!(r.status, PointStatus::Valid);
    assert!(r.warnings.is_empty());
}

#[test]
fn test_a_non_finite_twin_is_the_largest_disagreement() {
    let base = point_with_twins(vec![]);
    let p = point_with_twins(vec![twin_of(&base, |_| {}), twin_of(&base, |t| t.cm = f64::NAN)]);
    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.reasons, [Reason::IllConditioned]);
    let c = r.diagnostics.as_ref().unwrap().conditioning.as_ref().unwrap();
    assert_eq!(c.worst, "cm");
    assert!(c.worst_spread().is_nan());
    let json = serde_json::to_string(&r).unwrap();
    assert!(json.contains("\"cm\":null"), "a NaN spread serialises as null: {json}");
}

#[test]
fn test_unconverged_twins_warn_but_only_values_decide() {
    let base = point_with_twins(vec![]);
    // two twins ran out of iterations at values that agree with the base's
    let p = point_with_twins(vec![
        twin_of(&base, |t| t.converged = false),
        twin_of(&base, |t| t.converged = false),
        twin_of(&base, |_| {}),
    ]);
    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.status, PointStatus::Valid, "agreeing values: valid");
    assert_eq!(r.warnings, [Warning::TwinsNotValid { count: 2 }]);
    assert_eq!(r.warnings[0].to_string(), "solution converged but 2 twins did not");
    assert!(serde_json::to_string(&r).unwrap().contains("\"warnings\":[{\"TwinsNotValid\":{\"count\":2}}]"));
}

#[test]
fn test_no_twins_and_invalid_points_carry_no_conditioning() {
    // --no-twins: nothing measured, the record is as it was before twins existed
    let alone = PointResult {
        twins: None,
        ..point_with_twins(vec![])
    };
    let r = PointRecord::from_point(&alone, true, false, Thresholds::default());
    assert_eq!(r.status, PointStatus::Valid);
    assert!(r.diagnostics.unwrap().conditioning.is_none());
    assert!(!serde_json::to_string(&PointRecord::from_point(&alone, true, false, Thresholds::default()))
        .unwrap()
        .contains("conditioning"));

    // a point already invalid on its own computation is not judged on its twins
    let base = point_with_twins(vec![]);
    let mut p = point_with_twins(vec![twin_of(&base, |t| t.cl += 1.0)]);
    p.converged = false;
    let r = PointRecord::from_point(&p, true, false, Thresholds::default());
    assert_eq!(r.reasons, [Reason::ViscousNotConverged]);
    assert!(r.diagnostics.unwrap().conditioning.is_none());
}

#[test]
fn test_jitter_is_the_reference_twins_jitter() {
    // the directions read off the reference's twin 1 of kt_n60_inviscid_a4_cl05_cl8 (its
    // ulp1/panels.dat against panels.dat, nodes 0..8), pinned by value: yFoil's twin k and the
    // reference's twin k see the same panels, and the published twin data stays comparable
    let drawn: Vec<(i8, i8)> = (0..8).map(|i| jitter(1, i)).collect();
    assert_eq!(
        drawn,
        [(-1, 0), (-1, 1), (-1, 0), (-1, 0), (1, 1), (0, -1), (1, 0), (-1, 0)]
    );
}
