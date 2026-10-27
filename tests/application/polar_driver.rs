//! (c) yFoil functionality: the polar driver: `compute_polar` runs XFOIL's polar procedure (0 → max, `INIT`, re-solve
//! 0°, 0 → min) through one session, keeps the re-solved 0° out of the points, and stitches the
//! legs ascending; the re-solved 0° is bit-identical to the first.
//!
//! Fixtures: the panels of `tests/fixtures/xfoil/naca0012_n60_polar_re1e6/` — `cargo xtask fixtures --case naca0012_n60_polar_re1e6`.

use crate::fixtures;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{
    compute_polar, FlowConditions, HaltCriterion, HaltRule, PointResult, PolarConfig, Session,
};

fn spec() -> FlowConditions {
    FlowConditions {
        re: Some(1.0e6),
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    }
}

#[test]
fn test_polar_driver_reproduces_the_manual_sequence() {
    let path = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar_re1e6/panels.json");
    let geometry = read_geometry_from_file(path.to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);

    // the manual sequence: ALFA 0 / ASEQ 1 5 1 / INIT / ALFA 0 / ASEQ -1 -5 -1
    let mut session = Session::new(&airfoil, spec());
    let mut seq_points = vec![session.alpha(0.0)];
    for a in 1..=5 {
        seq_points.push(session.sequence_point((a as f64).to_radians()));
    }
    session.init();
    assert!(
        !session.state().bl_initialised && !session.state().pointers_built,
        "INIT clears LBLINI and LIPAN"
    );
    let reseed = session.alpha(0.0);
    for a in 1..=5 {
        seq_points.push(session.sequence_point(-(a as f64).to_radians()));
    }

    // the re-solved 0° is the first 0° again, bit for bit
    let first = &seq_points[0];
    assert_eq!(
        reseed.iterations, first.iterations,
        "yFoil re-solved 0°: iteration count"
    );
    assert_eq!(reseed.converged, first.converged);
    for (name, a, b) in [
        ("CL", reseed.cl, first.cl),
        ("CD", reseed.cd, first.cd),
        ("CM", reseed.cm, first.cm),
        ("CDF", reseed.cd_friction, first.cd_friction),
        ("CDP", reseed.cd_pressure, first.cd_pressure),
        ("CL_ALF", reseed.cl_d_alpha, first.cl_d_alpha),
        ("RMSBL", reseed.residual, first.residual),
        ("XOCTR1", reseed.transition_upper[0], first.transition_upper[0]),
        ("XOCTR2", reseed.transition_lower[0], first.transition_lower[0]),
    ] {
        assert_eq!(a.to_bits(), b.to_bits(), "yFoil re-solved 0° after INIT: {name}");
    }

    // compute_polar drives exactly this sequence, keeps the 0° seed out of the points, and
    // stitches ascending
    let polar = compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 5.0,
            alpha_min: -5.0,
            alpha_step: 1.0,
            conditions: spec(),
            ..PolarConfig::default()
        },
    );
    assert!(polar.completed);
    let converged: Vec<&PointResult> = seq_points.iter().filter(|p| p.converged).collect();
    assert_eq!(polar.results.len(), converged.len(), "compute_polar point count");
    for p in &polar.results {
        let q = converged
            .iter()
            .find(|q| (q.alpha - p.alpha).abs() < 1e-12)
            .expect("compute_polar alpha present in the manual sequence");
        assert_eq!(
            p.cl.to_bits(),
            q.cl.to_bits(),
            "compute_polar CL at {:.1}°",
            p.alpha.to_degrees()
        );
        assert_eq!(
            p.cd.to_bits(),
            q.cd.to_bits(),
            "compute_polar CD at {:.1}°",
            p.alpha.to_degrees()
        );
        assert_eq!(p.iterations, q.iterations);
    }
    assert!(
        polar.results.windows(2).all(|w| w[0].alpha < w[1].alpha),
        "points ascend in alpha"
    );
}

/// XFOIL's NSEQEX rule: an `ASEQ` halts after `after` points in a row fail to converge (4 by
/// default, as XFOIL), and the alphas it never reached are recorded as not attempted. With no
/// iterations at 0° and the rule at one failure, the first point of each leg gets ASEQ's five
/// iterations from an unconverged start, stops at RMSBL ≈ 2.7e-3 (27 × EPS1), and halts its leg.
#[test]
fn test_polar_driver_halts_a_leg_after_consecutive_failures() {
    assert_eq!(PolarConfig::default().halt_on, vec![HaltRule::XFOIL], "XFOIL's NSEQEX");
    assert_eq!(HaltRule::XFOIL.after, 4, "NSEQEX = 4 (xfoil.f:698)");
    let path = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar_re1e6/panels.json");
    let airfoil = panel_foil(&read_geometry_from_file(path.to_str().unwrap()).expect("panels.json"));
    let config = PolarConfig {
        alpha_max: 3.0,
        alpha_min: -3.0,
        alpha_step: 1.0,
        conditions: FlowConditions {
            max_iterations: 0,
            ..spec()
        },
        halt_on: vec![HaltRule {
            criterion: HaltCriterion::Xfoil,
            after: 1,
        }],
        ..PolarConfig::default()
    };
    let polar = compute_polar(&airfoil, &config);
    assert!(!polar.completed, "both legs halt");
    let deg = |v: &[f64]| -> Vec<f64> { v.iter().map(|a| (a.to_degrees() * 2.0).round() / 2.0).collect() };
    let solved: Vec<f64> = deg(&polar.results.iter().map(|p| p.alpha).collect::<Vec<_>>());
    assert_eq!(solved, vec![-1.0, 0.0, 1.0], "0° and the first point of each leg");
    assert!(polar.results.iter().all(|p| !p.converged && p.residual.is_finite()));
    assert_eq!(
        deg(&polar.not_attempted),
        vec![-3.0, -2.0, 2.0, 3.0],
        "the alphas past each halt"
    );
}

/// A classifier that calls the points at these alphas (degrees) ill-conditioned and nothing else.
fn ill_conditioned_at(alphas: &'static [f64]) -> yfoil::solver::analysis::Classifier {
    yfoil::solver::analysis::Classifier(std::sync::Arc::new(move |p: &PointResult| {
        let deg = (p.alpha.to_degrees() * 2.0).round() / 2.0;
        if alphas.contains(&deg) {
            vec![yfoil::output::ValidityClass::IllConditioned]
        } else {
            Vec::new()
        }
    }))
}

fn halting_polar(rules: Vec<HaltRule>, at: &'static [f64]) -> yfoil::solver::analysis::PolarResult {
    let path = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar_re1e6/panels.json");
    let airfoil = panel_foil(&read_geometry_from_file(path.to_str().unwrap()).expect("panels.json"));
    compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 3.0,
            alpha_min: -2.0,
            alpha_step: 1.0,
            conditions: spec(),
            halt_on: rules,
            classifier: Some(ill_conditioned_at(at)),
            ..PolarConfig::default()
        },
    )
}

const ILL: HaltCriterion = HaltCriterion::Class(yfoil::output::ValidityClass::IllConditioned);

/// `--halt-on ill-conditioned`: a leg halts at its first ill-conditioned point, and only that
/// leg; the alphas past it are not attempted, and the halt names its rule and where it stopped.
#[test]
fn test_a_class_rule_halts_the_leg_it_is_met_on() {
    let rule = HaltRule {
        criterion: ILL,
        after: 1,
    };
    let polar = halting_polar(vec![HaltRule::XFOIL, rule], &[2.0]);
    let deg = |v: &[f64]| -> Vec<f64> { v.iter().map(|a| (a.to_degrees() * 2.0).round() / 2.0).collect() };
    assert_eq!(
        deg(&polar.results.iter().map(|p| p.alpha).collect::<Vec<_>>()),
        [-2.0, -1.0, 0.0, 1.0, 2.0]
    );
    assert_eq!(deg(&polar.not_attempted), [3.0]);
    assert!(!polar.completed);
    assert_eq!(polar.halts.len(), 1);
    assert_eq!(polar.halts[0].rule, rule);
    assert_eq!((polar.halts[0].alpha.to_degrees() * 2.0).round() / 2.0, 2.0);
}

/// Each rule counts consecutive points and is reset by one that does not meet it, as XFOIL's
/// ISEQEX is reset by a converged point: two ill-conditioned points apart do not make two in a row.
#[test]
fn test_a_class_rule_counts_consecutive_points_only() {
    let apart = halting_polar(
        vec![HaltRule {
            criterion: ILL,
            after: 2,
        }],
        &[1.0, 3.0],
    );
    assert!(apart.completed && apart.not_attempted.is_empty() && apart.halts.is_empty());
    let together = halting_polar(
        vec![HaltRule {
            criterion: ILL,
            after: 2,
        }],
        &[1.0, 2.0],
    );
    assert_eq!(together.not_attempted.len(), 1, "halted after 1° and 2°, before 3°");
}

/// `--halt-on none`: no rule, so even XFOIL's NSEQEX does not stop a leg.
#[test]
fn test_no_rules_never_halt() {
    let path = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar_re1e6/panels.json");
    let airfoil = panel_foil(&read_geometry_from_file(path.to_str().unwrap()).expect("panels.json"));
    let polar = compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 3.0,
            alpha_min: -3.0,
            alpha_step: 1.0,
            conditions: FlowConditions {
                max_iterations: 0,
                ..spec()
            },
            halt_on: Vec::new(),
            ..PolarConfig::default()
        },
    );
    assert!(polar.completed && polar.not_attempted.is_empty());
    assert_eq!(polar.results.len(), 7, "every alpha, none converged");
}
