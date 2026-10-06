//! (c) yFoil functionality: point identity and the execution chain. A polar is a state machine, not N independent
//! solves: each point starts from the previous point's boundary layer. The result is presented
//! ascending in alpha, which is *not* the order it was solved in, so every record carries a
//! content-addressed `id` and `initialised_from` citing the record that seeded it.
//!
//! Fixtures: the panels of `tests/fixtures/xfoil/naca0012_n60_polar_re1e6/` — `cargo xtask fixtures --case naca0012_n60_polar_re1e6`.
use crate::fixtures;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::output::{PointStatus, PolarOutput};
use yfoil::solver::analysis::{compute_polar, FlowConditions, PolarConfig};

fn case_dir(case: &str) -> std::path::PathBuf {
    fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}"))
}

fn polar() -> PolarOutput {
    let dir = case_dir("naca0012_n60_polar_re1e6");
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let result = compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 4.0,
            alpha_min: -4.0,
            alpha_step: 1.0,
            conditions: FlowConditions {
                re: Some(1.0e6),
                ncrit: 9.0,
                max_iterations: 20,
                ..FlowConditions::default()
            },
            ..PolarConfig::default()
        },
    );
    PolarOutput::from_polar(&result, "naca0012", false)
}

#[test]
fn test_records_are_presented_ascending_in_alpha() {
    let p = polar();
    assert!(
        p.results.windows(2).all(|w| w[0].alpha_deg < w[1].alpha_deg),
        "results are ordered by alpha, not by when they were solved"
    );
}

#[test]
fn test_execution_order_is_recoverable_by_following_the_chain() {
    let p = polar();
    // walk each leg from its root: the chain is the only record of what ran when, since the
    // presentation is ascending in alpha and the ids carry no order of their own
    let next_after = |id: &str| p.results.iter().find(|r| r.initialised_from.as_deref() == Some(id));
    let root = p.results.iter().find(|r| r.initialised_from.is_none()).unwrap();
    assert!(root.alpha_deg.abs() < 1e-9, "the chain roots at the 0° solve");

    // two legs leave the root, so collect both and check they are the sweep XFOIL runs
    let legs: Vec<Vec<f64>> = p
        .results
        .iter()
        .filter(|r| r.initialised_from.as_deref() == root.id.as_deref())
        .map(|first| {
            let mut leg = vec![(first.alpha_deg * 2.0).round() / 2.0];
            let mut node = first;
            while let Some(n) = next_after(node.id.as_deref().unwrap()) {
                leg.push((n.alpha_deg * 2.0).round() / 2.0);
                node = n;
            }
            leg
        })
        .collect();
    assert_eq!(legs.len(), 2, "0° → max and 0° → min");
    let mut sorted = legs.clone();
    sorted.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
    assert_eq!(sorted[0], vec![-1.0, -2.0, -3.0, -4.0]);
    assert_eq!(sorted[1], vec![1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn test_ids_are_unique_and_well_formed() {
    let p = polar();
    let ids: Vec<&str> = p.results.iter().filter_map(|r| r.id.as_deref()).collect();
    assert!(ids.iter().all(|i| i.len() == yfoil::solver::point_id::ID_LENGTH
        && i.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit())));
    let mut unique: Vec<&str> = ids.clone();
    unique.sort_unstable();
    unique.dedup();
    assert_eq!(unique.len(), ids.len(), "ids collide");
}

#[test]
fn test_ids_are_deterministic() {
    // content-addressed, so the same run gives the same ids and a result file is
    // byte-reproducible
    let (a, b) = (polar(), polar());
    let ids = |p: &PolarOutput| -> Vec<Option<String>> { p.results.iter().map(|r| r.id.clone()).collect() };
    assert_eq!(ids(&a), ids(&b));
}

#[test]
fn test_the_id_covers_the_inputs_the_solve_depended_on() {
    let base = polar();
    let at = |p: &PolarOutput, a: f64| {
        p.results
            .iter()
            .find(|r| (r.alpha_deg - a).abs() < 1e-9)
            .unwrap()
            .id
            .clone()
    };
    // a different Reynolds number is a different computation at the same alpha
    let dir = case_dir("naca0012_n60_polar_re1e6");
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let other = PolarOutput::from_polar(
        &compute_polar(
            &airfoil,
            &PolarConfig {
                alpha_max: 4.0,
                alpha_min: -4.0,
                alpha_step: 1.0,
                conditions: FlowConditions {
                    re: Some(3.0e6),
                    ncrit: 9.0,
                    max_iterations: 20,
                    ..FlowConditions::default()
                },
                ..PolarConfig::default()
            },
        ),
        "naca0012",
        false,
    );
    assert_ne!(at(&base, 0.0), at(&other, 0.0), "Re is part of the identity");
    assert_ne!(at(&base, 2.0), at(&other, 2.0));
}

#[test]
fn test_the_predecessor_is_part_of_the_identity() {
    // 0° is solved twice in a sweep — once to open it, once after INIT to seed the downward
    // leg — and those two solves have identical inputs and no predecessor, so they are the
    // same computation and take the same id. That is what lets the downward leg cite the
    // recorded 0° point instead of a phantom: both legs cite one root (asserted in
    // `test_execution_order_is_recoverable_by_following_the_chain`).
    let p = polar();
    let at = |a: f64| p.results.iter().find(|r| (r.alpha_deg - a).abs() < 1e-9).unwrap();
    // +1 and -1 share a predecessor but differ in operating point, so differ in id
    assert_eq!(at(1.0).initialised_from, at(-1.0).initialised_from);
    assert_ne!(at(1.0).id, at(-1.0).id);
}

#[test]
fn test_every_reference_resolves() {
    let p = polar();
    for r in &p.results {
        if let Some(from) = r.initialised_from.as_deref() {
            assert!(
                p.results.iter().any(|q| q.id.as_deref() == Some(from)),
                "alpha {} cites id {from}, which is not in the result",
                r.alpha_deg
            );
        }
        if r.status == PointStatus::NotAttempted {
            assert!(
                r.id.is_none() && r.initialised_from.is_none(),
                "a point that was never solved has no place in the execution chain"
            );
        }
    }
}
