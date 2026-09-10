//! The polar break points (CLAUDE.md Rule 1, the third outcome): a 0.5° `ASEQ` sweep at
//! `ITER 100` driven past CL_max, where the Newton iteration stops converging and both codes
//! wander for 105 iterations per point.
//!
//! | case | script | what the reference does |
//! |---|---|---|
//! | `naca0012_n160_polar_up22_re1e6_iter100` | `ALFA 0 / ASEQ 0.5 22 0.5` | 20.5°–21.5° unconverged; 22° converges on the separated branch (CL ≈ 0.35) in 7 iterations |
//! | `naca4412_n160_polar_down16_re1e6_iter100` | `ALFA 0 / INIT / ALFA 0 / ASEQ -0.5 -16 -0.5` | −14.5°–−15.5° unconverged; −16° converges on the separated branch (CL ≈ +0.05) in 39 iterations |
//!
//! Every converged call up to the break matches the reference point by point and iteration by
//! iteration. The first unconverged call is threshold-straddling: the reference's own +1-ULP
//! twin parts from it (and flips IST/ITRAN) before yFoil does, and the one-step replays from
//! XFOIL's dumped state show the iterations themselves are faithful. Whether the fourth attempt
//! converges is decided at that noise floor — the twin fails where the reference happened to
//! converge — so yFoil's `NSEQEX` halt after four consecutive failures (0012) or its own landing
//! on the same separated branch at −16° (4412) are both consistent with the reference, and
//! neither is a gate difference nor a translation bug.
//!
//! Rule 6: the 105-iteration wanderings also take branches every other case left open —
//! MRCHUE's inverse (prescribed-Hk) marching in the wake, MRCHDU's extrapolation fallback for a
//! Newton failure with residual > 0.1, BLVAR's Us > 0.95 and Hk → 1 clamps, and TRCHEK2's
//! 30-iteration cap (`cargo xtask coverage`, `docs/validation/coverage.md`).

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::parse_bl_dump;
use std::path::PathBuf;
use utilities::records::{check_call, load, transient_tol, Outcome, Records};
use utilities::replay::{replay_iteration, DumpView};
use yfoil::geometry::{panel_foil, read_geometry_from_file, PanelledFoil};
use yfoil::solver::analysis::{compute_polar, FlowConditions, PointResult, PolarConfig, Session};
use yfoil::solver::specal::sequence_command;
use yfoil::solver::viscal::solve_viscous;

const UP_CASE: &str = "naca0012_n160_polar_up22_re1e6_iter100";
const DOWN_CASE: &str = "naca4412_n160_polar_down16_re1e6_iter100";
/// VISCAL call of the first unconverged point in each case (20.5° and −14.5°)
const UP_BREAK_CALL: usize = 42;
const DOWN_BREAK_CALL: usize = 31;
/// The iteration of that call at which the runs are allowed to part (see the tests' output)
const UP_STRADDLE_ITERATION: usize = 32;
const DOWN_STRADDLE_ITERATION: usize = 16;

fn spec() -> FlowConditions {
    FlowConditions {
        max_iterations: 100,
        ..FlowConditions::default()
    }
}

fn case_dir(case: &str) -> PathBuf {
    fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}"))
}

fn foil(dir: &std::path::Path) -> PanelledFoil {
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    panel_foil(&geometry)
}

/// The case's alpha script from its manifest: (`alphas`, `alphas_after_reinit`), degrees.
fn script(dir: &std::path::Path) -> (Vec<f64>, Vec<f64>) {
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    let list = |key: &str| -> Vec<f64> {
        m["case"][key]
            .as_array()
            .map(|a| a.iter().map(|v| v.as_f64().unwrap()).collect())
            .unwrap_or_default()
    };
    (list("alphas"), list("alphas_after_reinit"))
}

/// Drive the polar script through one session as the pipeline wrote it — `ALFA a0`, `ASEQ` for
/// the rest, then `INIT / ALFA a0 / ASEQ …` — stopping after VISCAL call `last` (0 = all).
/// Returns every point with its comparison outcome, in call order.
fn drive(
    session: &mut Session,
    rec: &Records,
    alphas: &[f64],
    after_reinit: &[f64],
    last: usize,
) -> Vec<(PointResult, Outcome)> {
    let mut out = Vec::new();
    let leg = |session: &mut Session, alphas: &[f64], out: &mut Vec<(PointResult, Outcome)>| {
        for (i, a) in alphas.iter().enumerate() {
            let k = out.len() + 1;
            if last != 0 && k > last {
                return;
            }
            let p = if i == 0 {
                session.alpha(a.to_radians())
            } else {
                session.sequence_point(a.to_radians())
            };
            let o = check_call(rec, k, &p, session.state(), transient_tol(k));
            out.push((p, o));
        }
    };
    leg(session, alphas, &mut out);
    if !after_reinit.is_empty() {
        session.init();
        leg(session, after_reinit, &mut out);
    }
    out
}

fn nitdone(rec: &Records, k: usize) -> usize {
    rec.points[k - 1]["NITDONE"].parse().unwrap()
}

fn lvconv(rec: &Records, k: usize) -> bool {
    rec.points[k - 1]["LVCONV"] == "T"
}

/// The sweep up to and including the break call: every call before it matches, the break call
/// straddles at `expected_iteration`, and every call after it starts from a state the reference
/// cannot reproduce (straddling from iteration 0).
fn check_break(case: &str, break_call: usize, expected_iteration: usize) -> (Records, Vec<(PointResult, Outcome)>) {
    let dir = case_dir(case);
    let rec = load(&dir);
    let (alphas, after) = script(&dir);
    assert_eq!(
        rec.points.len(),
        alphas.len() + after.len(),
        "{case}: VISCAL call count"
    );
    let mut session = Session::new(&foil(&dir), spec());
    let run = drive(&mut session, &rec, &alphas, &after, 0);
    for (i, (_, o)) in run.iter().enumerate() {
        let k = i + 1;
        match k.cmp(&break_call) {
            std::cmp::Ordering::Less => assert_eq!(*o, Outcome::Match, "{case} call {k}"),
            std::cmp::Ordering::Equal => match o {
                Outcome::Straddling { at_iteration, why } => {
                    println!("{case} call {k}: straddling at iteration {at_iteration}: {why}");
                    assert_eq!(*at_iteration, expected_iteration, "{case} call {k}: straddle iteration");
                }
                Outcome::Match => panic!("{case} call {k} now matches to the end (the straddle closed): pin it"),
            },
            std::cmp::Ordering::Greater => assert!(
                matches!(o, Outcome::Straddling { at_iteration: 0, .. }),
                "{case} call {k}: {o:?}"
            ),
        }
    }
    // the break call and the two after it ran out of iterations in both codes
    for k in break_call..break_call + 3 {
        assert_eq!(
            nitdone(&rec, k),
            spec().max_iterations + 5,
            "{case} call {k}: reference NITDONE"
        );
        assert!(!lvconv(&rec, k), "{case} call {k}: reference LVCONV");
        assert_eq!(
            run[k - 1].0.iterations,
            spec().max_iterations + 5,
            "{case} call {k}: yFoil iterations"
        );
        assert!(!run[k - 1].0.converged, "{case} call {k}: yFoil converged");
    }
    // the fourth attempt: the reference converged, its +1-ULP twin did not
    let k4 = break_call + 3;
    assert!(lvconv(&rec, k4), "{case} call {k4}: the reference converged");
    let flips = &rec.floor.as_ref().expect("noise_floor.json").flips;
    assert!(
        flips.iter().any(|f| f == &format!("call {k4}: LVCONV T vs F")),
        "{case} call {k4}: the reference's +1-ULP twin did not converge; flips: {flips:?}"
    );
    (rec, run)
}

#[test]
fn test_naca0012_upward_break_matches_then_straddles() {
    let (rec, run) = check_break(UP_CASE, UP_BREAK_CALL, UP_STRADDLE_ITERATION);
    // 22°: the reference lands on the separated branch in 7 iterations; yFoil, like the twin,
    // runs out, which is its fourth consecutive failure
    let k4 = UP_BREAK_CALL + 3;
    assert_eq!(nitdone(&rec, k4), 7);
    assert!(!run[k4 - 1].0.converged);
    // … so the polar halts as XFOIL's ASEQ would (NSEQEX = 4): 41 points, 20.5°–22° failed
    let polar = compute_polar(
        &foil(&case_dir(UP_CASE)),
        &PolarConfig {
            alpha_max: 22.0,
            alpha_min: 0.0,
            alpha_step: 0.5,
            conditions: spec(),
            ..PolarConfig::default()
        },
    );
    assert_eq!(polar.results.len(), UP_BREAK_CALL - 1);
    assert_eq!(
        polar
            .failed_alphas
            .iter()
            .map(|a| (a.to_degrees() * 2.0).round() / 2.0)
            .collect::<Vec<_>>(),
        [20.5, 21.0, 21.5, 22.0]
    );
    assert!(!polar.completed, "halted by NSEQEX consecutive failures");
}

#[test]
fn test_naca4412_downward_break_matches_then_straddles() {
    let (rec, run) = check_break(DOWN_CASE, DOWN_BREAK_CALL, DOWN_STRADDLE_ITERATION);
    // −16°: both codes land on the separated branch (CL ≈ +0.05 against ≈ −0.9 at −14°), by
    // different paths and within the RMSBL < 1e-4 convergence band of each other
    let k4 = DOWN_BREAK_CALL + 3;
    let p = &run[k4 - 1].0;
    assert!(p.converged, "{DOWN_CASE} call {k4}: yFoil converged");
    let cl_ref: f64 = rec.points[k4 - 1]["CL"].parse().unwrap();
    let cl_attached: f64 = rec.points[DOWN_BREAK_CALL - 2]["CL"].parse().unwrap();
    println!(
        "{DOWN_CASE} call {k4} (−16°): yFoil CL {:.8} in {} iterations, reference CL {cl_ref:.8} in {} iterations, attached branch at −14° CL {cl_attached:.8}",
        p.cl,
        p.iterations,
        nitdone(&rec, k4)
    );
    assert!(
        (p.cl - cl_ref).abs() < (cl_attached - cl_ref).abs() / 2.0,
        "{DOWN_CASE} call {k4}: yFoil CL {} is not on the reference's separated branch (CL {cl_ref}, attached {cl_attached})",
        p.cl
    );
    // the polar keeps −16° and carries on: 30 points, three failed, not halted
    let polar = compute_polar(
        &foil(&case_dir(DOWN_CASE)),
        &PolarConfig {
            alpha_max: 0.0,
            alpha_min: -16.0,
            alpha_step: 0.5,
            conditions: spec(),
            ..PolarConfig::default()
        },
    );
    assert_eq!(polar.results.len(), rec.points.len() - 1 - 3, "0°, −0.5°…−14°, −16°");
    assert_eq!(polar.failed_alphas.len(), 3);
    assert!(polar.completed);
}

/// Replay SETBL/UPDATE call `k` of a polar case from XFOIL's dumped state: the session is
/// driven through every call before the break call, then SPECAL and the VISCAL prologue for the
/// break call's alpha, so the wake, DIJ and pointers are the sweep's own.
fn replay_break_iteration(case: &str, break_call: usize, k: usize) {
    let dir = case_dir(case);
    let rec = load(&dir);
    let (alphas, after) = script(&dir);
    let before: usize = (1..break_call).map(|c| nitdone(&rec, c)).sum();
    assert!(
        k > before && k <= before + nitdone(&rec, break_call),
        "{case}: SETBL call {k} is not in VISCAL call {break_call}"
    );
    let iteration = k - before;
    let d = parse_bl_dump(&dir.join(format!("mrchdu_input_{k}.dat")));
    let o = parse_bl_dump(&dir.join(format!("update_output_{k}.dat")));
    let (d, o) = (
        DumpView {
            header: &d.header,
            bl: &d.bl,
        },
        DumpView {
            header: &o.header,
            bl: &o.bl,
        },
    );
    let session = Session::new(&foil(&dir), spec());
    let all: Vec<f64> = alphas.iter().chain(&after).copied().collect();
    let prepare = |session: &mut Session| {
        let run = drive(session, &rec, &alphas, &after, break_call - 1);
        assert!(
            run.iter().all(|(_, o)| *o == Outcome::Match),
            "{case}: every call before the break matches"
        );
        let (st, sys) = session.parts_mut();
        sequence_command(st, sys, all[break_call - 1].to_radians());
        solve_viscous(st, sys.as_mut(), 0, spec().wake_length, None);
    };
    replay_iteration(case, &dir, &d, &o, session, prepare, break_call, iteration, k);
}

#[test]
fn test_naca0012_break_iterations_replay_from_xfoil_state() {
    for k in [237, 238] {
        replay_break_iteration(UP_CASE, UP_BREAK_CALL, k);
    }
}

#[test]
fn test_naca4412_break_iterations_replay_from_xfoil_state() {
    for k in [181, 182] {
        replay_break_iteration(DOWN_CASE, DOWN_BREAK_CALL, k);
    }
}
