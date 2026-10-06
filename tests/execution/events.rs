//! (b) Execution equivalence, event by event: the non-finite reference runs, one rarely-taken
//! branch at a time.
//!
//! The `non-finite` cases of `cases.toml` are the probes of the branch-coverage study whose
//! reference run prints NaN: a whole-run gate is impossible there (the reference cannot
//! reproduce itself, and NaN arithmetic is where two correct translations may legitimately
//! differ). What *can* be gated is each rarely-taken branch they reach, as XFOIL's own event:
//! the instrumented reference writes `events.dat` — one line per event with its site, SETBL call,
//! side, station and the value that decided it — and `dump_calls` keeps the state dumps of those
//! calls. Each test seeds yFoil with XFOIL's exact state entering that call (`step.rs`: the
//! inviscid side from the panels, the BL side from the dump), runs the one
//! subroutine (or the one iteration) in question, asserts that yFoil's branch trace records the
//! same events, and compares the resulting arrays with the reference's dump: within the usual
//! `TOL_SOLVER` where both are finite, and NaN where the reference has NaN (a NaN in the reference
//! and a number in yFoil, or the reverse, is a mismatch).
//!
//! | event | branch | case, SETBL call | what is replayed |
//! |---|---|---|---|
//! | `MRCHUE_GARBAGE` | xbl.f:788–827 | 16-212 M 0.7, 4412 N=160 M 0.5 … (first march) | MRCHUE from the VISCAL prologue, against `mrchdu_input_1.dat` |
//! | `MRCHDU_GARBAGE` | xbl.f:1111–1135 | every non-finite case, call 1 | MRCHDU from `mrchdu_input_1.dat`, against `mrchdu_output_1.dat` |
//! | `UPDATE_RDN4` | xbl.f:1474 | 64A010 TYPE 2 M 0.3, call 7 | SETBL/BLSOLV/UPDATE from `mrchdu_input_7.dat`, against `update_output_7.dat` |
//! | `UPDATE_ISLAND` | xbl.f:1537 | 4412 M 0.6 ITER 60, call 51 | the same at call 51 |
//! | `STMOVE_UEPS` | xpanel.f:1737 | 63-415 M 0.5, call 18 | UPDATE 18 then QVFUE/GAMQV/STMOVE, against `mrchdu_input_19.dat` |
//! | `STFIND_NOT_FOUND` | xpanel.f:1364 | 16-212 M 0.7, call 1 | UPDATE 1 then STMOVE, against `mrchdu_input_2.dat` |
//! | `STFIND_GAM_NEG` | xpanel.f:1365 | 63-415 M 0.5, VISCAL prologue | STFIND on the prologue state |
//! | `TRCHEK2_TRIP_UPSTREAM` | xblsys.f:432 | 16-212 M 0.7, call 2 | MRCHDU at call 2 (the check runs inside it) |
//! | `BLVAR_WAKE_US` | xblsys.f:843 | 4412 N=160 M 0.5, call 8 | MRCHDU and the iteration at call 8 |
//!
//! These replays are the only gates that exercise the fallbacks with non-finite inputs, and when
//! they were first written (2026-09-13) two of them failed on translation defects the finite
//! cases could not see: the garbage test had been written as `if dmax > 0.1`, the complement of
//! XFOIL's `IF(DMAX .LE. 0.1) GO TO 109`, which is the opposite for a NaN residual
//! (docs/xfoil-known-issues.md §7.4); and MRCHUE's failure path made a single closure call at a
//! wake station where XFOIL makes two, `BLVAR(2)` then `BLVAR(3)` (§7.5). Both are fixed; keep
//! these tests in CI for exactly that reason.

//! Fixtures: `tests/fixtures/xfoil/<case>/` of the five non-finite probes (`events.dat`, the dumps of
//! their event calls) — `cargo xtask fixtures --case <case>`.

use crate::fixtures::cases::{dump as step_dump, jog_dump, load, prologue_jogged, seed, Case};
use crate::fixtures::mrchdu_fixtures::BlDump;
use crate::utilities::tolerances::TOL_SOLVER;
use std::path::Path;
use yfoil::bl::blsolv::solve_newton_system;
use yfoil::bl::mrchdu::{march_prescribed_dstar, MrchduTrace};
use yfoil::bl::mrchue::{march_direct, MrchueTrace};
use yfoil::bl::system::FlowParameters;
use yfoil::solver::analysis::Session;
use yfoil::solver::pointers::{find_stagnation, move_stagnation};
use yfoil::solver::setbl::assemble_newton_system;
use yfoil::solver::setbl::set_mach_re_from_cl;
use yfoil::solver::update::{apply_newton_update, UpdateEvent};
use yfoil::solver::velocity::{set_gamma_from_q_viscous, set_q_viscous_from_ue};

/// One line of `events.dat`: `EV <site> <SETBL call> <IS> <IBL> <value>`.
#[derive(Debug, Clone)]
struct Event {
    site: String,
    call: usize,
    side: usize,
    station: usize,
    value: f64,
}

fn events(dir: &Path) -> Vec<Event> {
    std::fs::read_to_string(dir.join("events.dat"))
        .expect("events.dat")
        .lines()
        .filter_map(|l| {
            let t: Vec<&str> = l.split_whitespace().collect();
            (t.len() == 6 && t[0] == "EV").then(|| Event {
                site: t[1].to_string(),
                call: t[2].parse().unwrap(),
                side: t[3].parse().unwrap(),
                station: t[4].parse().unwrap(),
                value: t[5].parse().unwrap(),
            })
        })
        .collect()
}

fn at(ev: &[Event], site: &str, call: usize) -> Vec<(usize, usize)> {
    ev.iter()
        .filter(|e| e.site == site && e.call == call)
        .map(|e| (e.station, e.side))
        .collect()
}

fn case(name: &str) -> Case {
    load(name)
}

/// The seeds of the jogged copies of each step (as `step.rs`).
const JOGS: [u64; 5] = [1, 2, 3, 4, 5];

/// Sessions at the VISCAL prologue of the case's one ALFA point (wake, QINV, pointers, UINV, DIJ;
/// no BL iteration yet): the exact one, then the copies on jogged panels.
fn prologue(c: &Case) -> (Session, Vec<Session>) {
    (
        prologue_jogged(c, 1, None, None),
        JOGS.iter().map(|&s| prologue_jogged(c, 1, None, Some(s))).collect(),
    )
}

fn dump(c: &Case, name: &str) -> BlDump {
    step_dump(c, name)
}

/// Sessions seeded with XFOIL's state entering SETBL call `k`: the exact one, then the copies whose
/// panels and dumped state are jogged by one ULP (`step.rs`).
fn seeded(c: &Case, k: usize) -> (Session, Vec<Session>) {
    let d = dump(c, &format!("mrchdu_input_{k}.dat"));
    let make = |jog: Option<u64>| {
        let dd = match jog {
            Some(s) => jog_dump(&d, s),
            None => d.clone(),
        };
        let mut session = prologue_jogged(c, 1, Some(dd.real("ALFA")), jog);
        seed(&mut session, &dd, TOL_SOLVER);
        session
    };
    (make(None), JOGS.iter().map(|&s| make(Some(s))).collect())
}

/// The BL parameters SETBL derives for this state (MRCL for the current CL, then COMSET…).
fn params_of(session: &mut Session) -> (FlowParameters, [f64; 3]) {
    let st = session.state_mut();
    let clmr = st.cl;
    let _ = set_mach_re_from_cl(st, clmr);
    let mut params = FlowParameters::new(st.mach, st.re, st.gamma_gas);
    params.amplification_model = st.amplification_model;
    (params, st.ncrit)
}

/// Equal where the reference is non-finite (both NaN, or both the same infinity); otherwise within
/// `TOL_SOLVER` on the scale or within four times the step's own sensitivity `spread` (how far the
/// jogged copies moved the value; non-finite spreads count as none), whichever is larger.
fn same_or_within(a: f64, b: f64, scale: f64, spread: f64) -> bool {
    if a.is_nan() || b.is_nan() {
        return a.is_nan() && b.is_nan();
    }
    if a.is_infinite() || b.is_infinite() {
        return a == b;
    }
    let spread = if spread.is_finite() { spread } else { 0.0 };
    (a - b).abs() <= (TOL_SOLVER * a.abs().max(b.abs()).max(scale)).max(4.0 * spread)
}

/// How far the jogged copies moved a value of the exact one.
fn spread(ours: f64, jogged: impl Iterator<Item = f64>) -> f64 {
    jogged
        .map(|j| (j - ours).abs())
        .filter(|d| d.is_finite())
        .fold(0.0_f64, f64::max)
}

/// Compare the state's BL arrays with a dump's `BL` rows (XSSI UEDG THET DSTR CTAU MASS);
fn assert_arrays(session: &Session, jogged: &[Session], o: &BlDump, what: &str) {
    let st = session.state();
    let value = |s: &Session, is: usize, ibl: usize, m: usize| -> f64 {
        let t = s.state();
        [
            t.xi[is][ibl],
            t.ue[is][ibl],
            t.theta[is][ibl],
            t.dstar[is][ibl],
            t.sqrtctau[is][ibl],
            t.mass_defect[is][ibl],
        ][m]
    };
    let names = ["XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS"];
    let mut worst = 0.0_f64;
    let mut nonfinite = 0;
    let mut mismatches: Vec<String> = vec![];
    for is in 1..=2 {
        let nrows = o.bl[is].len() - 1;
        assert!(
            nrows >= st.n_stations[is],
            "{what}: the dump has {nrows} rows on side {is}, the state {} stations",
            st.n_stations[is]
        );
        let nrows = st.n_stations[is];
        for ibl in 2..=nrows {
            let ours = [
                st.xi[is][ibl],
                st.ue[is][ibl],
                st.theta[is][ibl],
                st.dstar[is][ibl],
                st.sqrtctau[is][ibl],
                st.mass_defect[is][ibl],
            ];
            for (m, name) in names.iter().enumerate() {
                if is == 1 && ibl > st.i_te_station[1] && *name == "MASS" {
                    continue; // MASS is never written on side 1's wake slots (known issues §2.8)
                }
                let b = o.bl[is][ibl][m];
                let scale = (2..=nrows)
                    .map(|j| o.bl[is][j][m])
                    .filter(|v| v.is_finite())
                    .fold(0.0_f64, |acc, v| acc.max(v.abs()));
                let sp = spread(ours[m], jogged.iter().map(|j| value(j, is, ibl, m)));
                if !same_or_within(ours[m], b, scale, sp) {
                    mismatches.push(format!(
                        "{name}({ibl},{is}): yfoil={:.17e} xfoil={b:.17e} |diff|={:.3e} (scale {scale:.3e})",
                        ours[m],
                        (ours[m] - b).abs()
                    ));
                }
                if b.is_finite() {
                    worst = worst.max((ours[m] - b).abs() / ours[m].abs().max(b.abs()).max(scale));
                } else {
                    nonfinite += 1;
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{what}: {} entries outside tolerance, the first: {}",
        mismatches.len(),
        mismatches.iter().take(5).cloned().collect::<Vec<_>>().join("; ")
    );
    println!("{what}: arrays match (worst finite relative difference {worst:.2e}, {nonfinite} non-finite entries reproduced)");
}

fn assert_itran(session: &Session, o: &BlDump, what: &str) {
    assert_eq!(
        session.state().i_transition_station[1..],
        [o.int("ITRAN1"), o.int("ITRAN2")],
        "{what}: ITRAN"
    );
}

/// MRCHUE (the first march) with its garbage-extrapolation events, against `mrchdu_input_1`.
fn check_mrchue(name: &str) {
    let c = case(name);
    let ev = events(&c.dir);
    // the first march runs inside SETBL call 1; the instrumentation stamps its events with the
    // counter as it was at MRCHUE's call (0 or 1)
    let mut expected: Vec<(usize, usize)> = at(&ev, "MRCHUE_GARBAGE", 0);
    expected.extend(at(&ev, "MRCHUE_GARBAGE", 1));
    assert!(!expected.is_empty(), "{name}: no MRCHUE_GARBAGE events");
    let (mut session, mut jogged) = prologue(&c);
    let (params, acrit) = params_of(&mut session);
    let mut trace = MrchueTrace::default();
    march_direct(session.state_mut(), &params, acrit, Some(&mut trace));
    for j in jogged.iter_mut() {
        let (p, a) = params_of(j);
        march_direct(j.state_mut(), &p, a, None);
    }
    let ours: Vec<(usize, usize)> = trace.garbage.iter().map(|(ibl, is, _)| (*ibl, *is)).collect();
    assert_eq!(
        ours, expected,
        "{name}: MRCHUE garbage-extrapolation stations (IBL, IS)"
    );
    let o = dump(&c, "mrchdu_input_1.dat");
    assert_itran(&session, &o, &format!("{name} MRCHUE"));
    assert_arrays(&session, &jogged, &o, &format!("{name} MRCHUE"));
    println!("{name}: MRCHUE garbage extrapolation at {ours:?} reproduced");
}

/// MRCHDU at SETBL call `k` with its garbage-extrapolation events, against `mrchdu_output_k`.
fn check_mrchdu(name: &str, k: usize) {
    let c = case(name);
    let ev = events(&c.dir);
    let expected = at(&ev, "MRCHDU_GARBAGE", k);
    let (mut session, mut jogged) = seeded(&c, k);
    let (params, acrit) = params_of(&mut session);
    let mut trace = MrchduTrace::default();
    march_prescribed_dstar(session.state_mut(), &params, acrit, Some(&mut trace));
    for j in jogged.iter_mut() {
        let (p, a) = params_of(j);
        march_prescribed_dstar(j.state_mut(), &p, a, None);
    }
    let ours: Vec<(usize, usize)> = trace.garbage.iter().map(|(ibl, is, _)| (*ibl, *is)).collect();
    assert_eq!(
        ours, expected,
        "{name} call {k}: MRCHDU garbage-extrapolation stations (IBL, IS)"
    );
    let o = dump(&c, &format!("mrchdu_output_{k}.dat"));
    assert_itran(&session, &o, &format!("{name} MRCHDU call {k}"));
    assert_arrays(&session, &jogged, &o, &format!("{name} MRCHDU call {k}"));
    println!("{name} call {k}: MRCHDU garbage extrapolation at {ours:?} reproduced");
}

/// One SETBL/BLSOLV/UPDATE iteration at call `k`; returns the session after UPDATE and its events.
fn iterate(c: &Case, k: usize) -> (Session, Vec<Session>, Vec<UpdateEvent>) {
    let (mut session, mut jogged) = seeded(c, k);
    let one = |s: &mut Session| {
        let st = s.state_mut();
        let r = assemble_newton_system(st);
        let sol = solve_newton_system(r.newton);
        apply_newton_update(st, &sol.deltas, st.mach_d_cl)
    };
    let u = one(&mut session);
    let uj: Vec<_> = jogged.iter_mut().map(one).collect();
    let o = dump(c, &format!("update_output_{k}.dat"));
    for (what, ours, theirs, sp) in [
        (
            "RLX",
            u.relaxation,
            o.real("RLX"),
            spread(u.relaxation, uj.iter().map(|x| x.relaxation)),
        ),
        (
            "RMSBL",
            u.residual,
            o.real("RMSBL"),
            spread(u.residual, uj.iter().map(|x| x.residual)),
        ),
    ] {
        assert!(
            same_or_within(ours, theirs, 1.0, sp),
            "call {k}: {what}: yfoil={ours:.17e} xfoil={theirs:.17e} (step sensitivity {sp:.2e})"
        );
    }
    assert_arrays(&session, &jogged, &o, &format!("UPDATE call {k}"));
    (session, jogged, u.events)
}

#[test]
fn test_mrchue_garbage_extrapolation_events() {
    for name in [
        "naca16-212_n60_a4_re1e6_m07",
        "naca63-415_n60_a14_re1e6_m05",
        "naca4412_n60_a16_re1e6_m06_iter60",
        "naca4412_n160_a16_re1e6_m05",
        "naca64a010_n60_type2_m03_a0",
    ] {
        check_mrchue(name);
    }
}

#[test]
fn test_mrchdu_garbage_extrapolation_events_call_1() {
    for name in [
        "naca16-212_n60_a4_re1e6_m07",
        "naca63-415_n60_a14_re1e6_m05",
        "naca4412_n60_a16_re1e6_m06_iter60",
        "naca4412_n160_a16_re1e6_m05",
        "naca64a010_n60_type2_m03_a0",
    ] {
        check_mrchdu(name, 1);
    }
}

#[test]
fn test_mrchdu_trip_upstream_of_a_laminar_station_call_2() {
    // TRCHEK2's `TRFORC = XIFORC.GT.X1 .AND. XIFORC.LE.X2` with XIFORC <= X1 fires inside this
    // MRCHDU; the gate is the march reproducing XFOIL's transition and arrays through it
    let name = "naca16-212_n60_a4_re1e6_m07";
    let ev = events(&case(name).dir);
    assert!(
        !at(&ev, "TRCHEK2_TRIP_UPSTREAM", 2).is_empty(),
        "{name}: no TRCHEK2_TRIP_UPSTREAM at call 2"
    );
    check_mrchdu(name, 2);
}

#[test]
fn test_update_ue_relaxation_below_low_call_7() {
    let name = "naca64a010_n60_type2_m03_a0";
    let c = case(name);
    let expected: Vec<UpdateEvent> = at(&events(&c.dir), "UPDATE_RDN4", 7)
        .into_iter()
        .map(|(station, side)| UpdateEvent::UeRelaxationBelowLow { side, station })
        .collect();
    assert!(!expected.is_empty(), "{name}: no UPDATE_RDN4 at call 7");
    let (_, _, ours) = iterate(&c, 7);
    let ours: Vec<UpdateEvent> = ours
        .into_iter()
        .filter(|e| matches!(e, UpdateEvent::UeRelaxationBelowLow { .. }))
        .collect();
    assert_eq!(ours, expected, "{name} call 7: UPDATE Ue relaxation-below-low events");
}

#[test]
fn test_update_negative_ue_island_call_51() {
    let name = "naca4412_n60_a16_re1e6_m06_iter60";
    let c = case(name);
    let expected: Vec<UpdateEvent> = at(&events(&c.dir), "UPDATE_ISLAND", 51)
        .into_iter()
        .map(|(station, side)| UpdateEvent::NegativeUeIsland { side, station })
        .collect();
    assert!(!expected.is_empty(), "{name}: no UPDATE_ISLAND at call 51");
    let (_, _, ours) = iterate(&c, 51);
    let ours: Vec<UpdateEvent> = ours
        .into_iter()
        .filter(|e| matches!(e, UpdateEvent::NegativeUeIsland { .. }))
        .collect();
    assert_eq!(ours, expected, "{name} call 51: UPDATE negative-Ue island events");
}

/// UPDATE at call `k`, then QVFUE/GAMQV/STMOVE as VISCAL runs them, against the state entering
/// call `k + 1`.
fn stmove_after(c: &Case, k: usize) -> (Session, yfoil::solver::pointers::StagnationMove) {
    let (mut session, mut jogged, _) = iterate(c, k);
    let tail = |s: &mut Session| {
        let st = s.state_mut();
        set_q_viscous_from_ue(st);
        set_gamma_from_q_viscous(st);
        move_stagnation(st)
    };
    let moved = tail(&mut session);
    for j in jogged.iter_mut() {
        tail(j);
    }
    let st = session.state();
    let next = dump(c, &format!("mrchdu_input_{}.dat", k + 1));
    assert_eq!(st.i_stagnation_node, next.int("IST"), "call {k}: IST after STMOVE");
    assert!(
        same_or_within(
            st.s_stagnation,
            next.real("SST"),
            1.0,
            spread(st.s_stagnation, jogged.iter().map(|j| j.state().s_stagnation))
        ),
        "call {k}: SST after STMOVE: yfoil={:.17e} xfoil={:.17e}",
        st.s_stagnation,
        next.real("SST")
    );
    assert_arrays(&session, &jogged, &next, &format!("STMOVE after call {k}"));
    (session, moved)
}

#[test]
fn test_stmove_ue_floor_after_call_18() {
    let name = "naca63-415_n60_a14_re1e6_m05";
    let c = case(name);
    let expected: Vec<(usize, usize)> = at(&events(&c.dir), "STMOVE_UEPS", 18)
        .into_iter()
        .map(|(station, side)| (side, station))
        .collect();
    assert!(!expected.is_empty(), "{name}: no STMOVE_UEPS at call 18");
    let (_, moved) = stmove_after(&c, 18);
    assert_eq!(
        moved.ue_floored, expected,
        "{name} call 18: stations floored to UEPS (side, station)"
    );
}

#[test]
fn test_stfind_not_found_after_call_1() {
    let name = "naca16-212_n60_a4_re1e6_m07";
    let c = case(name);
    assert!(
        !at(&events(&c.dir), "STFIND_NOT_FOUND", 1).is_empty(),
        "{name}: no STFIND_NOT_FOUND at call 1"
    );
    let (_, moved) = stmove_after(&c, 1);
    assert!(
        !moved.scan.found,
        "{name} call 1: STFIND must not find a sign change (IST = N/2)"
    );
}

#[test]
fn test_stfind_negative_te_vorticity_in_the_prologue() {
    let name = "naca63-415_n60_a14_re1e6_m05";
    let c = case(name);
    let ev = events(&c.dir);
    let first = ev
        .iter()
        .find(|e| e.site == "STFIND_GAM_NEG" && e.call == 0)
        .unwrap_or_else(|| panic!("{name}: no STFIND_GAM_NEG in the prologue"));
    let (mut session, jogged) = prologue(&c);
    let scan = find_stagnation(session.state_mut());
    assert!(scan.found, "{name}: the prologue scan finds the stagnation point");
    assert_eq!(
        scan.first_negative_node,
        Some(first.station),
        "{name}: first node with GAM < 0 in the scan"
    );
    assert!(
        same_or_within(
            session.state().gamma[first.station],
            first.value,
            1.0,
            spread(
                session.state().gamma[first.station],
                jogged.iter().map(|j| j.state().gamma[first.station])
            )
        ),
        "{name}: GAM at that node: yfoil={:.17e} xfoil={:.17e}",
        session.state().gamma[first.station],
        first.value
    );
}

#[test]
fn test_blvar_wake_us_clamp_call_8() {
    // BLVAR's wake Us clamp fires inside MRCHDU and SETBL at this call; the gate is the march
    // and the iteration reproducing XFOIL through it
    let name = "naca4412_n160_a16_re1e6_m05";
    let c = case(name);
    assert!(
        !at(&events(&c.dir), "BLVAR_WAKE_US", 8).is_empty(),
        "{name}: no BLVAR_WAKE_US at call 8"
    );
    check_mrchdu(name, 8);
    let _ = iterate(&c, 8);
}
