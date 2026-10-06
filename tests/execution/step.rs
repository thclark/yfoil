//! The single-step replay used by the (b) branch tests (`branches.rs`): yFoil seeded with XFOIL's
//! complete state on entry to one step of a reference run, required to take XFOIL's decisions and
//! reproduce its state on exit.
//!
//! A step is one VISCAL iteration (iteration `i` of call `c`, global SETBL call `k`). The branch
//! cover names a whole VISCAL call (`iteration = 0`) where the point is part of an `ASEQ`; every
//! iteration of that call is then replayed as its own step. Steps are chosen by `cargo xtask steps`, which admits only
//! decisions are compared exactly. Values are compared within the named tolerance (`TOL_SOLVER`
//! on the host the fixture was generated on, `TOL_CROSS_HOST` elsewhere, `common/utilities/host.rs`)
//! or within four times the step's own sensitivity, whichever is larger: the same step is also
//! replayed from five jogged copies of its inputs (every panel coordinate, from which yFoil builds
//! the inviscid side, and every BL entry of XFOIL's state moved by −1, 0 or +1 ULP), and a value may differ from XFOIL's by
//! up to four times the furthest any jog moved it. That is the step's own conditioning, measured in
//! the test: an ill-conditioned solution amplifies a last-bit difference within one step, and the
//! test allows for exactly that amplification and no more (`docs/conventions/testing.md`).
//!
//! The replay is self-contained. The inviscid side (SPECAL or SPECCL, the wake, QINV, DIJ) is a
//! function of the panels and the call's operating point, so yFoil builds it from the case's
//! `panels.json`; the BL side and the pointers (IST, SST, ITRAN, the arrays, ALFA, CL) come from
//! XFOIL's dump `mrchdu_input_<k>.dat`. The first iteration of a fresh run has no dump: it is
//! replayed from the panels, MRCHUE included.
//!
//! Compared, per replayed iteration (SETBL call `k`): UPDATE's RLX, RMSBL and the magnitude of
//! RMXBL against `update_output_<k>.dat` — not which variable and station UPDATE reports as the
//! largest change, a console label that rounding decides where two changes tie (known issues §7.8;
//! it is compared exactly on the well-conditioned reference case, `update.rs`); then, after the iteration's
//! tail (QVFUE, GAMQV, STMOVE, CLCALC, CDCALC), either the state entering the next iteration
//! (`mrchdu_input_<k+1>.dat`: IST, SST, ITRAN, XSSI, UEDG, THET, DSTR, CTAU, MASS) or, after the
//! call's last iteration, the converged point (`viscal_points.dat`).

use crate::fixtures;
use crate::fixtures::cases::{blocks, dump, jog_dump, load, prologue_jogged, seed, Op};
use crate::fixtures::mrchdu_fixtures::BlDump;
use crate::utilities::host::same_host;
use crate::utilities::tolerances::{assert_within, TOL_CROSS_HOST, TOL_SOLVER};
use std::collections::HashMap;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::Session;
use yfoil::solver::specal::{alpha_command, cl_command};
use yfoil::solver::viscal::{solve_viscous, IterationRecord};

/// The state's BL arrays against a dump's `BL` rows, each within the tolerance on the array's
/// scale (its largest magnitude on that side) or within four times the furthest any jogged
/// session's entry moved.
fn assert_arrays(session: &Session, jogged: &[Session], o: &BlDump, what: &str, tol: f64) {
    let get = |s: &Session, is: usize, ibl: usize, m: usize| -> f64 {
        let st = s.state();
        [
            st.xi[is][ibl],
            st.ue[is][ibl],
            st.theta[is][ibl],
            st.dstar[is][ibl],
            st.sqrtctau[is][ibl],
            st.mass_defect[is][ibl],
        ][m]
    };
    let st = session.state();
    let names = ["XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS"];
    for is in 1..=2 {
        let nrows = st.n_stations[is];
        assert_eq!(o.nbl(is), nrows, "{what}: NBL({is})");
        for (m, name) in names.iter().enumerate() {
            let scale = (2..=nrows)
                .map(|j| o.bl[is][j][m].abs())
                .filter(|v| v.is_finite())
                .fold(0.0_f64, f64::max);
            for ibl in 2..=nrows {
                // MASS is never written on side 1's wake slots (known issues §2.8)
                if is == 1 && ibl > st.i_te_station[1] && *name == "MASS" {
                    continue;
                }
                let ours = get(session, is, ibl, m);
                let spread = jogged
                    .iter()
                    .map(|j| (get(j, is, ibl, m) - ours).abs())
                    .fold(0.0_f64, f64::max);
                assert_close(
                    ours,
                    o.bl[is][ibl][m],
                    tol,
                    scale,
                    spread,
                    &format!("{what}: {name}({ibl},{is})"),
                );
            }
        }
    }
}

/// `|a − b| ≤ max(tol · max(|a|, |b|, scale), 4 · spread)`: the named tolerance, or four times how
/// far the same step moved this value when its input was jogged by one ULP. NaN where XFOIL has NaN.
fn assert_close(a: f64, b: f64, tol: f64, scale: f64, spread: f64, what: &str) {
    // a non-finite reference value (a run XFOIL drove through a NaN) must be matched exactly
    if !b.is_finite() || !a.is_finite() {
        assert!(
            a.is_nan() && b.is_nan() || a == b,
            "{what}: yfoil={a:.17e} xfoil={b:.17e}"
        );
        return;
    }
    let spread = if spread.is_finite() { spread } else { 0.0 };
    let lim = (tol * a.abs().max(b.abs()).max(scale)).max(4.0 * spread);
    assert!(
        (a - b).abs() <= lim,
        "{what}: yfoil={a:.17e} xfoil={b:.17e} |diff|={:.3e} > allowed {lim:.3e} (step sensitivity {spread:.2e})",
        (a - b).abs()
    );
}

/// The seeds of the jogged copies of a step.
const JOGS: [u64; 5] = [1, 2, 3, 4, 5];

/// Replay step (`call`, `iteration`) of `case`, whose first SETBL call is `setbl`. Iteration 0 is
/// a whole call: each of its iterations is replayed as its own step, from XFOIL's state on entry
/// to it, so no iteration inherits the drift of the ones before.
pub fn replay(case: &str, call: usize, iteration: usize, setbl: usize) {
    if iteration == 0 {
        let c = load(case);
        let nit: usize = c.points[call - 1]["NITDONE"].parse().unwrap();
        for it in 1..=nit {
            replay(case, call, it, setbl + it - 1);
        }
        return;
    }
    let c = load(case);
    let tol = if same_host(&c.dir) { TOL_SOLVER } else { TOL_CROSS_HOST };
    let nit: usize = c.points[call - 1]["NITDONE"].parse().unwrap();
    let (first, last) = if iteration == 0 {
        (1, nit)
    } else {
        (iteration, iteration)
    };
    let fresh = call == 1 && first == 1;
    // the exact replay, then the jogged copies of it
    let start = |jog: Option<u64>| -> Session {
        if fresh {
            prologue_jogged(&c, call, None, jog)
        } else {
            let d = dump(&c, &format!("mrchdu_input_{setbl}.dat"));
            let d = match jog {
                Some(s) => jog_dump(&d, s),
                None => d,
            };
            let alpha = (!matches!(c.ops[call - 1], Op::Cl(_))).then(|| d.real("ALFA"));
            let mut s = prologue_jogged(&c, call, alpha, jog);
            seed(&mut s, &d, tol);
            s
        }
    };
    let mut session = start(None);
    let mut jogged: Vec<Session> = JOGS.iter().map(|&s| start(Some(s))).collect();
    let wake_length = c.spec.wake_length;
    let step = |s: &mut Session| -> IterationRecord {
        let mut trace: Vec<IterationRecord> = vec![];
        let (st, sys) = s.parts_mut();
        solve_viscous(st, sys.as_mut(), 1, wake_length, Some(&mut trace));
        trace.remove(0)
    };
    for it in first..=last {
        let k = setbl + (it - first);
        let what = format!("{case} call {call} iteration {it} (SETBL {k})");
        let r = step(&mut session);
        let rj: Vec<IterationRecord> = jogged.iter_mut().map(step).collect();
        let spread =
            |f: &dyn Fn(&IterationRecord) -> f64| rj.iter().map(|j| (f(j) - f(&r)).abs()).fold(0.0_f64, f64::max);
        let o = dump(&c, &format!("update_output_{k}.dat"));
        assert_close(
            r.relaxation,
            o.real("RLX"),
            tol,
            1.0,
            spread(&|x| x.relaxation),
            &format!("{what}: RLX"),
        );
        assert_close(
            r.residual,
            o.real("RMSBL"),
            tol,
            1.0,
            spread(&|x| x.residual),
            &format!("{what}: RMSBL"),
        );
        assert_close(
            r.residual_max.abs(),
            o.real("RMXBL").abs(),
            tol,
            1.0,
            spread(&|x| x.residual_max.abs()),
            &format!("{what}: |RMXBL|"),
        );
        if it < nit {
            let next = dump(&c, &format!("mrchdu_input_{}.dat", k + 1));
            let st = session.state();
            assert_eq!(st.i_stagnation_node, next.int("IST"), "{what}: IST after STMOVE");
            let sst = jogged
                .iter()
                .map(|j| (j.state().s_stagnation - st.s_stagnation).abs())
                .fold(0.0_f64, f64::max);
            assert_close(
                st.s_stagnation,
                next.real("SST"),
                tol,
                1.0,
                sst,
                &format!("{what}: SST"),
            );
            assert_eq!(
                st.i_transition_station[1..],
                [next.int("ITRAN1"), next.int("ITRAN2")],
                "{what}: ITRAN"
            );
            assert_arrays(
                &session,
                &jogged,
                &next,
                &format!("{what}, entering SETBL {}", k + 1),
                tol,
            );
        } else {
            let p = &c.points[call - 1];
            let st = session.state();
            assert_eq!(r.converged, p["LVCONV"] == "T", "{what}: convergence");
            assert_eq!(st.i_stagnation_node, p["IST"].parse::<usize>().unwrap(), "{what}: IST");
            assert_eq!(
                st.i_transition_station[1..],
                [
                    p["ITRAN1"].parse::<usize>().unwrap(),
                    p["ITRAN2"].parse::<usize>().unwrap()
                ],
                "{what}: ITRAN"
            );
            let get = |s: &Session| -> [f64; 7] {
                let t = s.state();
                [
                    t.cl,
                    t.cm,
                    t.cd,
                    t.cd_friction,
                    t.cd_pressure,
                    t.x_transition[1],
                    t.x_transition[2],
                ]
            };
            let ours = get(&session);
            for (n, name) in ["CL", "CM", "CD", "CDF", "CDP", "XOCTR1", "XOCTR2"].iter().enumerate() {
                let sp = jogged
                    .iter()
                    .map(|j| (get(j)[n] - ours[n]).abs())
                    .fold(0.0_f64, f64::max);
                assert_close(
                    ours[n],
                    p[*name].parse().unwrap(),
                    tol,
                    1.0,
                    sp,
                    &format!("{what}: {name}"),
                );
            }
        }
    }
}

/// Replay inviscid operating point `call` of `case` — one SPECAL (`ALFA`) or one SPECCL (`CL`) —
/// in a fresh session, and require XFOIL's record of it (`specal_points.dat` /
/// `speccl_points.dat`): ALFA, CL, CM, CDP, MINF and, per node, GAM, QINV and CPI within
/// `TOL_SOLVER`, and for SPECCL the alpha Newton's exit iteration ITAL exactly. Nothing carries
/// between inviscid points except the starting alpha of SPECCL's Newton, which is taken from the
/// record of the point before.
pub fn replay_inviscid(case: &str, call: usize) {
    let c = load(case);
    let tol = if same_host(&c.dir) { TOL_SOLVER } else { TOL_CROSS_HOST };
    let g = read_geometry_from_file(c.dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let mut session = Session::new(&panel_foil(&g), c.spec.clone());
    let n_before = |pred: fn(&Op) -> bool| c.ops[..call].iter().filter(|o| pred(o)).count();
    let (file, idx) = match c.ops[call - 1] {
        Op::Cl(_) => ("speccl_points.dat", n_before(|o| matches!(o, Op::Cl(_)))),
        _ => ("specal_points.dat", n_before(|o| !matches!(o, Op::Cl(_)))),
    };
    let what = format!("{case} inviscid point {call}");
    let mut ital = None;
    {
        let (st, sys) = session.parts_mut();
        match c.ops[call - 1] {
            Op::Alfa(a) | Op::Aseq(a) => alpha_command(st, sys, a.to_radians()),
            Op::Cl(cl) => {
                // SPECCL starts its Newton from the alpha the previous point left
                if call > 1 {
                    let (pf, pi) = match c.ops[call - 2] {
                        Op::Cl(_) => ("speccl_points.dat", idx - 1),
                        _ => (
                            "specal_points.dat",
                            c.ops[..call - 1].iter().filter(|o| !matches!(o, Op::Cl(_))).count(),
                        ),
                    };
                    let prev =
                        &blocks(&fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}/{pf}")))[pi - 1];
                    alpha_command(st, sys, prev["ALFA"].parse().unwrap());
                }
                ital = Some(cl_command(st, sys, cl));
            }
        }
    }
    let path = fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}/{file}"));
    let text = std::fs::read_to_string(&path).unwrap();
    // one record per call: its scalars and its nodes' GAM QINV CPI
    type Record = (HashMap<String, f64>, Vec<[f64; 3]>);
    let mut recs: Vec<Record> = vec![];
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "CALL" {
            recs.push((HashMap::new(), vec![]));
        }
        let Some(cur) = recs.last_mut() else { continue };
        if k.trim_start().starts_with("NODE") {
            let x: Vec<f64> = v.split_whitespace().map(|t| t.parse().unwrap()).collect();
            cur.1.push([x[0], x[1], x[2]]);
        } else if let Ok(x) = v.trim().parse::<f64>() {
            cur.0.insert(k.trim().to_string(), x);
        }
    }
    let (h, nodes) = &recs[idx - 1];
    let st = session.state();
    if let Some(ital) = ital {
        assert_eq!(ital as f64, h["ITAL"], "{what}: SPECCL's exit iteration ITAL");
    }
    for (name, ours) in [
        ("ALFA", st.alpha),
        ("CL", st.cl),
        ("CM", st.cm),
        ("CDP", st.cd_pressure),
        ("MINF", st.mach),
    ] {
        assert_within(ours, h[name], tol, 1.0, &format!("{what}: {name}"));
    }
    assert_eq!(nodes.len(), st.n_foil_nodes, "{what}: node count");
    for (i, r) in nodes.iter().enumerate() {
        let n = i + 1;
        assert_within(st.gamma[n], r[0], tol, 1.0, &format!("{what}: GAM({n})"));
        assert_within(st.q_inviscid[n], r[1], tol, 1.0, &format!("{what}: QINV({n})"));
        assert_within(st.cp_inviscid[n], r[2], tol, 1.0, &format!("{what}: CPI({n})"));
    }
}
