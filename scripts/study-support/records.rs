//! Comparison of a `Session` run against the VISCAL-level reference records every case keeps:
//! `viscal_points.dat` (one block per VISCAL call), `viscal_iters_all.dat` (one `IT` line per
//! iteration: `k i RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2 ALFA MINF REINF`),
//! `viscal_iter.dat` (call 1 only, with UPDATE's RMXBL/VMXBL/IMXBL/ISMXBL) and, when the pipeline
//! recorded it, `noise_floor.json` — the same run on the case's five seeded 1-ULP twins (every
//! panel coordinate jogged by −1, 0 or +1 ULP), each value's spread the largest over the set.
//!
//! Exact: iteration count, LVCONV, IST, ITRAN. Values: `|a − b| ≤ max(tol·max(|a|,|b|,scale),
//! FLOOR_FACTOR · floor)` where `floor` is the reference's own 1-ULP spread of that very value
//! (CLAUDE.md Rule 1: tolerances are the measured floor times a safety factor).
//!
//! Third outcome — **divergent** — reported, never silently passed or failed:
//! - a 1-ULP twin itself changes the branch trace (iteration count, convergence, IST, ITRAN)
//!   in an earlier call: every later call starts from a state the reference cannot reproduce,
//!   so the call is divergent from iteration 0 (a flip inside the call itself is left to the
//!   per-iteration classification below, and a call that still matches to its end while the
//!   twin flipped inside it is reported as divergent rather than passed);
//! - every earlier iteration matched, and at this iteration the reference's own 1-ULP spread
//!   exceeds `DIVERGENCE_FLOOR`: the runs are allowed to part here — in any gated transient (RMSBL,
//!   RLX, CL, …), in the branch trace (IST/ITRAN), in the reported limiter, or at the convergence
//!   test itself (`RMSBL < EPS1`: the iteration counts differ, every common iteration matched, and
//!   the reference's RMSBL at the last common one moves by more than `DIVERGENCE_FLOOR`) — the
//!   one-step replay from XFOIL's exact state is the evidence that the step itself is faithful;
//! - UPDATE's *reported* limiter (VMXBL/IMXBL: the largest normalised change) differs while
//!   |RMXBL| agrees within the floor: two near-tied changes, RLX itself unaffected — a tie.

#![allow(dead_code)]

use super::tolerances::{TOL_SOLVER, TOL_TRANSIENT};

/// Safety factor applied to a case's own measured 1-ULP floor (`noise_floor.json`, written by
/// `cargo xtask twins` from the case's five seeded 1-ULP twins, the largest spread over them): a
/// value is accepted when it is within the base tolerance *or* within `FLOOR_FACTOR` × the
/// reference's own spread of that value. A study threshold: the tests do not read the floor
/// (`docs/conventions/testing.md`).
pub const FLOOR_FACTOR: f64 = 4.0;
/// Hypersensitivity threshold for the studies' third outcome. When one of the reference's own
/// 1-ULP twins moves a per-iteration value by more than this (absolute), the reference cannot
/// reproduce itself there; a yFoil run that matched every earlier iteration within
/// `FLOOR_FACTOR` × floor and departs at such an iteration is classified *divergent*
/// (`docs/conventions/terminology.md`), reported, not passed or failed.
pub const DIVERGENCE_FLOOR: f64 = 1e-6;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use yfoil::solver::analysis::PointResult;
use yfoil::solver::blstate::SolverState;

#[derive(Debug, Clone, Default)]
pub struct CallFloor {
    pub point: HashMap<String, f64>,
    pub iterations: Vec<HashMap<String, f64>>,
    /// the twin's spread (max over nodes / stations) of every array of the
    /// call's final state (`viscal_state_<k>.dat`), by array name
    pub state: HashMap<String, f64>,
}

#[derive(Debug, Clone, Default)]
pub struct Floor {
    pub branch_identical: bool,
    pub flips: Vec<String>,
    pub calls: Vec<CallFloor>,
    /// inviscid-only cases: per SPECAL / SPECCL call, the twin's spread of the point values and
    /// the per-node maxima (GAM, QINV, CPI)
    pub specal_calls: Vec<HashMap<String, f64>>,
    pub speccl_calls: Vec<HashMap<String, f64>>,
    /// per dumped SETBL/UPDATE call: max |base − twin| per post-UPDATE array (XSSI UEDG THET DSTR CTAU MASS)
    pub update_output: HashMap<usize, HashMap<String, f64>>,
}

pub struct Records {
    pub points: Vec<HashMap<String, String>>,
    /// inviscid-only cases: one block per SPECAL call (`specal_points.dat`) and per SPECCL call
    /// (`speccl_points.dat`), `NODE(i)` entries holding `GAM QINV CPI`
    pub specal: Vec<HashMap<String, String>>,
    pub speccl: Vec<HashMap<String, String>>,
    /// per call (1-based key): rows of the IT line values after the call number
    pub iters: HashMap<usize, Vec<Vec<f64>>>,
    /// VISCAL call 1's per-iteration blocks of viscal_iter.dat (RMXBL/VMXBL/IMXBL/ISMXBL)
    pub iter1: Vec<HashMap<String, String>>,
    pub floor: Option<Floor>,
}

/// What a call comparison concluded.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Outcome {
    #[default]
    Match,
    /// Reported separately: the reason, and the iteration at which the runs were allowed to part
    /// (0 = the whole case, from the twin's own branch flips)
    Divergent { at_iteration: usize, why: String },
    /// Values or branch trace differ beyond tolerance with the reference reproducing itself:
    /// a translation bug (`check_call` panics on it)
    Mismatch { why: String },
}

fn kv_blocks(text: &str, start_key: &str) -> Vec<HashMap<String, String>> {
    let mut out: Vec<HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == start_key {
            out.push(HashMap::new());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

pub fn load(case_dir: &Path) -> Records {
    let blocks = |name: &str| {
        std::fs::read_to_string(case_dir.join(name))
            .map(|t| kv_blocks(&t, "CALL"))
            .unwrap_or_default()
    };
    let points = blocks("viscal_points.dat");
    let specal = blocks("specal_points.dat");
    let speccl = blocks("speccl_points.dat");
    assert!(
        !points.is_empty() || !specal.is_empty(),
        "{}: neither viscal_points.dat nor specal_points.dat",
        case_dir.display()
    );
    let mut iters: HashMap<usize, Vec<Vec<f64>>> = HashMap::new();
    for l in std::fs::read_to_string(case_dir.join("viscal_iters_all.dat"))
        .unwrap_or_default()
        .lines()
    {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(|t| t.parse().unwrap()).collect();
        iters.entry(v[0] as usize).or_default().push(v[1..].to_vec());
    }
    let iter1 = std::fs::read_to_string(case_dir.join("viscal_iter.dat"))
        .map(|t| kv_blocks(&t, "ITER"))
        .unwrap_or_default();
    let floor = std::fs::read_to_string(case_dir.join("noise_floor.json"))
        .ok()
        .map(|t| {
            let j: serde_json::Value = serde_json::from_str(&t).expect("noise_floor.json");
            let to_map = |v: &serde_json::Value| -> HashMap<String, f64> {
                v.as_object()
                    // a non-finite spread (the reference itself diverged to NaN there) is `null`
                    .map(|o| {
                        o.iter()
                            .map(|(k, x)| (k.clone(), x.as_f64().unwrap_or(f64::NAN)))
                            .collect()
                    })
                    .unwrap_or_default()
            };
            Floor {
                branch_identical: j["branch_identical"].as_bool().unwrap_or(false),
                flips: j["flips"]
                    .as_array()
                    .map(|a| a.iter().map(|s| s.as_str().unwrap_or("").to_string()).collect())
                    .unwrap_or_default(),
                calls: j["calls"]
                    .as_array()
                    .map(|a| {
                        a.iter()
                            .map(|c| CallFloor {
                                point: to_map(&c["point"]),
                                iterations: c["iterations"]
                                    .as_array()
                                    .map(|it| it.iter().map(to_map).collect())
                                    .unwrap_or_default(),
                                state: to_map(&c["state"]),
                            })
                            .collect()
                    })
                    .unwrap_or_default(),
                specal_calls: j["specal_calls"]
                    .as_array()
                    .map(|a| a.iter().map(|c| to_map(&c["point"])).collect())
                    .unwrap_or_default(),
                speccl_calls: j["speccl_calls"]
                    .as_array()
                    .map(|a| a.iter().map(|c| to_map(&c["point"])).collect())
                    .unwrap_or_default(),
                update_output: j["update_output"]
                    .as_object()
                    .map(|o| {
                        o.iter()
                            .filter_map(|(k, v)| k.parse().ok().map(|k| (k, to_map(v))))
                            .collect()
                    })
                    .unwrap_or_default(),
            }
        });
    Records {
        points,
        specal,
        speccl,
        iters,
        iter1,
        floor,
    }
}

/// The arrays of a `viscal_state_<k>.dat` record (the final state of VISCAL call `k`), by
/// name: the node arrays `CPI CPV QINV QVIS GAM` (one entry per node `1..=N+NW`) and the
/// station arrays `XSSI X UEDG THET DSTR CTAU MASS TAU` (one entry per `BL(is, ibl)` row, side
/// 1's `1..=NBL(1)` then side 2's `1..=NBL(2)`). `None` when the case did not keep the record.
pub fn read_state_arrays(path: &Path) -> Option<BTreeMap<String, Vec<f64>>> {
    let text = std::fs::read_to_string(path).ok()?;
    let mut arrays: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for l in text.lines() {
        let (names, rest): (&[&str], &str) = if l.starts_with("NODE(") {
            (&STATE_NODE_ARRAYS, l.split_once(")=")?.1)
        } else if l.starts_with("BL(") {
            (&STATE_STATION_ARRAYS, l.split_once(")=")?.1)
        } else {
            continue;
        };
        let v: Vec<f64> = rest.split_whitespace().map(|t| t.parse().unwrap_or(f64::NAN)).collect();
        for (m, name) in names.iter().enumerate() {
            arrays
                .entry(name.to_string())
                .or_default()
                .push(v.get(m).copied().unwrap_or(f64::NAN));
        }
    }
    Some(arrays)
}

pub const STATE_NODE_ARRAYS: [&str; 5] = ["CPI", "CPV", "QINV", "QVIS", "GAM"];
pub const STATE_STATION_ARRAYS: [&str; 8] = ["XSSI", "X", "UEDG", "THET", "DSTR", "CTAU", "MASS", "TAU"];

/// yFoil's counterpart of `read_state_arrays`: the same arrays in the same order from a
/// `SolverState` after a point (station rows `1..=n_stations[is]` per side).
pub fn state_arrays(st: &SolverState) -> BTreeMap<String, Vec<f64>> {
    let np = st.n_foil_nodes + st.n_wake_nodes;
    let mut arrays: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let node = |v: &Vec<f64>| -> Vec<f64> { (1..=np).map(|i| v.get(i).copied().unwrap_or(f64::NAN)).collect() };
    arrays.insert("CPI".into(), node(&st.cp_inviscid));
    arrays.insert("CPV".into(), node(&st.cp_viscous));
    arrays.insert("QINV".into(), node(&st.q_inviscid));
    arrays.insert("QVIS".into(), node(&st.q_viscous));
    arrays.insert("GAM".into(), node(&st.gamma));
    let mut station = |name: &str, f: &dyn Fn(usize, usize) -> f64| {
        let mut out = vec![];
        for is in 1..=2 {
            for ibl in 1..=st.n_stations[is] {
                out.push(f(is, ibl));
            }
        }
        arrays.insert(name.into(), out);
    };
    station("XSSI", &|is, ibl| st.xi[is][ibl]);
    station("X", &|is, ibl| st.x[st.i_node[is][ibl]]);
    station("UEDG", &|is, ibl| st.ue[is][ibl]);
    station("THET", &|is, ibl| st.theta[is][ibl]);
    station("DSTR", &|is, ibl| st.dstar[is][ibl]);
    station("CTAU", &|is, ibl| st.sqrtctau[is][ibl]);
    station("MASS", &|is, ibl| st.mass_defect[is][ibl]);
    station("TAU", &|is, ibl| st.tau[is][ibl]);
    arrays
}

/// One array's comparison of a call's final state: the worst |yfoil − xfoil| over its entries,
/// where (0-based file index; node number for a node array, row index for a station array),
/// the twin's floor for the array, and the ratio of the difference to the gate
/// `allowed(TOL_SOLVER, floor)` — above 1 the array is outside what the reference's own 1-ULP
/// twin allows.
#[derive(Debug, Clone)]
pub struct ArrayDiff {
    pub name: String,
    pub diff: f64,
    pub index: usize,
    pub floor: f64,
    pub ratio: f64,
    /// the two values at the worst entry
    pub yfoil: f64,
    pub xfoil: f64,
}

/// Compare yFoil's state after a point with the reference's final state of VISCAL call `k`
/// (`viscal_state_<k>.dat` in `case_dir`), array by array. `None` when the case has no state
/// record for the call; arrays of different lengths (the codes ended with a different NBL) are
/// compared over the common prefix and flagged by `ratio = +inf`. Sorted worst ratio first.
pub fn compare_state(rec: &Records, case_dir: &Path, k: usize, st: &SolverState) -> Option<Vec<ArrayDiff>> {
    let theirs = read_state_arrays(&case_dir.join(format!("viscal_state_{k}.dat")))?;
    let ours = state_arrays(st);
    let floors = rec.floor.as_ref().and_then(|f| f.calls.get(k - 1)).map(|c| &c.state);
    let mut out = vec![];
    for (name, b) in &theirs {
        let a = ours.get(name)?;
        let floor = floors.and_then(|f| f.get(name)).copied().unwrap_or(0.0);
        let mut worst = ArrayDiff {
            name: name.clone(),
            diff: 0.0,
            index: 0,
            floor,
            ratio: if a.len() == b.len() { 0.0 } else { f64::INFINITY },
            yfoil: f64::NAN,
            xfoil: f64::NAN,
        };
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            let d = (x - y).abs();
            // NaN on both sides is agreement (the non-finite region); on one side it is a
            // difference no gate accepts
            let (d, ratio) = if x.is_nan() && y.is_nan() {
                (0.0, 0.0)
            } else if d.is_nan() {
                (f64::INFINITY, f64::INFINITY)
            } else {
                (d, d / allowed(*x, *y, TOL_SOLVER, 1.0, floor))
            };
            if ratio > worst.ratio || (ratio == worst.ratio && d > worst.diff) {
                worst.diff = d;
                worst.index = i;
                worst.ratio = ratio;
                worst.yfoil = *x;
                worst.xfoil = *y;
            }
        }
        out.push(worst);
    }
    out.sort_by(|p, q| q.ratio.partial_cmp(&p.ratio).unwrap_or(std::cmp::Ordering::Equal));
    Some(out)
}

/// The allowed |a − b| for one value: the base tolerance or the measured floor × safety factor.
pub fn allowed(a: f64, b: f64, tol: f64, scale: f64, floor: f64) -> f64 {
    (tol * a.abs().max(b.abs()).max(scale)).max(FLOOR_FACTOR * floor)
}

fn assert_value(a: f64, b: f64, tol: f64, scale: f64, floor: f64, what: &str) {
    let lim = allowed(a, b, tol, scale, floor);
    let d = (a - b).abs();
    assert!(
        d <= lim,
        "{what}: yfoil={a:.17e} xfoil={b:.17e} |diff|={d:.3e} > allowed {lim:.3e} (tol {tol:.0e}, scale {scale:.1e}, 1-ULP floor {floor:.3e} × {FLOOR_FACTOR})"
    );
}

/// The VISCAL call a twin flip line refers to (`call 42: …`, `call 42 iteration 17: …`); `None`
/// for a line without one (`VISCAL call count …`), which concerns every call.
pub fn flip_call(line: &str) -> Option<usize> {
    line.strip_prefix("call ")
        .and_then(|r| r.split(|c: char| !c.is_ascii_digit()).next())
        .and_then(|n| n.parse().ok())
}

/// Everything a call comparison measured, whether or not it passed. `outcome` is the verdict;
/// the diffs are what the branch-coverage study tabulates (`scripts/branch-coverage`).
#[derive(Debug, Clone, Default)]
pub struct CallReport {
    pub call: usize,
    /// yFoil's alpha at the end of the call, degrees
    pub alpha_deg: f64,
    /// (yfoil, xfoil) iteration counts
    pub iterations: (usize, usize),
    pub converged: (bool, bool),
    /// yFoil's point as computed: CL, CD, CM and the last iteration's RMSBL (0 for inviscid)
    pub yfoil: [f64; 4],
    /// the reference's, from `viscal_points.dat` (NaN when the record lacks it)
    pub xfoil: [f64; 4],
    /// transition x/c on the upper and lower side at the end of the point, yFoil and the
    /// reference (XOCTR1/2; NaN for an inviscid point)
    pub yfoil_transition: [f64; 2],
    pub xfoil_transition: [f64; 2],
    /// IST, ITRAN(1), ITRAN(2) at the end of the point, yFoil and the reference
    pub yfoil_stations: [usize; 3],
    pub xfoil_stations: [usize; 3],
    /// iterations compared within tolerance before the runs were allowed to part (all of them
    /// on a match)
    pub gated_iterations: usize,
    /// worst |diff| / (the reference's own 1-ULP floor) over the gated transients
    pub worst_floor_ratio: f64,
    /// max |yfoil − xfoil| over the gated iterations, per transient (RMSBL, RLX, CL, CD, CM, ALFA)
    pub transient_diff: Vec<(&'static str, f64)>,
    /// |yfoil − xfoil| at the converged point (CL, CD, CM, CDF, CDP, XOCTR1, XOCTR2), only when
    /// the point itself was gated (a match)
    pub point_diff: Vec<(&'static str, f64)>,
    pub outcome: Outcome,
}

fn value_ok(a: f64, b: f64, tol: f64, scale: f64, floor: f64, what: &str) -> Result<(), String> {
    // a non-finite reference value is reproduced only by the same non-finite value (NaN by NaN,
    // an infinity by the same infinity): the branch trace through it is still gated exactly
    if a.is_nan() && b.is_nan() {
        return Ok(());
    }
    if (a.is_infinite() || b.is_infinite()) && a == b {
        return Ok(());
    }
    let lim = allowed(a, b, tol, scale, floor);
    let d = (a - b).abs();
    if d <= lim {
        Ok(())
    } else {
        Err(format!(
            "{what}: yfoil={a:.17e} xfoil={b:.17e} |diff|={d:.3e} > allowed {lim:.3e} (tol {tol:.0e}, scale {scale:.1e}, 1-ULP floor {floor:.3e} × {FLOOR_FACTOR})"
        ))
    }
}

/// Compare VISCAL call `k` (1-based) of the reference with a `Session` result, asserting the
/// gates (`compare_call` with a panic on a mismatch).
pub fn check_call(rec: &Records, k: usize, p: &PointResult, st: &SolverState, transient_tol: f64) -> Outcome {
    let r = compare_call(rec, k, p, st, transient_tol);
    if let Outcome::Mismatch { why } = &r.outcome {
        panic!("{why}");
    }
    r.outcome
}

/// Compare VISCAL call `k` (1-based) of the reference with a `Session` result; never panics.
/// `transient_tol` is the base tolerance for the per-iteration values.
pub fn compare_call(rec: &Records, k: usize, p: &PointResult, st: &SolverState, transient_tol: f64) -> CallReport {
    let x = &rec.points[k - 1];
    let ctx = format!("call {k} (alpha {:.3}°)", p.alpha.to_degrees());
    let xf = |key: &str| x.get(key).and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN);
    let xi = |key: &str| x.get(key).and_then(|v| v.trim().parse::<usize>().ok()).unwrap_or(0);
    let mut report = CallReport {
        call: k,
        alpha_deg: p.alpha.to_degrees(),
        iterations: (p.iterations, x["NITDONE"].parse::<usize>().unwrap()),
        converged: (p.converged, x["LVCONV"] == "T"),
        yfoil: [p.cl, p.cd, p.cm, p.residual],
        xfoil: [xf("CL"), xf("CD"), xf("CM"), xf("RMSBL")],
        yfoil_transition: [p.transition_upper[0], p.transition_lower[0]],
        xfoil_transition: [xf("XOCTR1"), xf("XOCTR2")],
        yfoil_stations: [
            st.i_stagnation_node,
            p.i_transition_station[1],
            p.i_transition_station[2],
        ],
        xfoil_stations: [xi("IST"), xi("ITRAN1"), xi("ITRAN2")],
        ..Default::default()
    };
    let flips_at = |pred: &dyn Fn(Option<usize>) -> bool| -> Vec<String> {
        rec.floor
            .as_ref()
            .map(|f| f.flips.iter().filter(|l| pred(flip_call(l))).cloned().collect())
            .unwrap_or_default()
    };
    // the twin changed its branch trace in an earlier call: this call is divergent
    // from iteration 0 whatever its values do. The comparison still runs below, so that
    // `gated_iterations` and the diffs record how far yFoil follows anyway (the branch-case
    // polars tabulate it); no value can turn the outcome into a mismatch here
    let earlier = flips_at(&|c| c.is_none_or(|c| c < k));
    let ungated = if earlier.is_empty() {
        None
    } else {
        let why = format!(
            "one of the reference's own 1-ULP twins changes its branch trace before this call: {}",
            earlier.join("; ")
        );
        println!(
            "{ctx}: DIVERGENT — {why}. yfoil: {} iterations, converged={}, CL {:.8}; reference: {} iterations, converged={}, CL {}",
            p.iterations, p.converged, p.cl, x["NITDONE"], x["LVCONV"], x["CL"]
        );
        Some(why)
    };
    let ungated_flag = ungated.is_some();
    let mismatch = |report: &mut CallReport, why: String| {
        if !ungated_flag {
            println!("{ctx}: MISMATCH — {why}");
        }
        // for an ungated call this is converted to the divergent outcome after the loop
        report.outcome = Outcome::Mismatch { why };
    };
    assert_eq!(k, x["CALL"].parse::<usize>().unwrap());

    let within_call = flips_at(&|c| c == Some(k));
    let cf = rec.floor.as_ref().and_then(|f| f.calls.get(k - 1));
    let fl = |m: Option<&HashMap<String, f64>>, name: &str| m.and_then(|h| h.get(name)).copied().unwrap_or(0.0);

    // the twin itself changes this call's iteration count or convergence: the count is then not
    // a gate, the common iterations are, and the call is divergent at the end of them
    let count_flip = within_call
        .iter()
        .any(|l| l.contains("NITDONE") || l.contains("iteration count") || l.contains("LVCONV"));
    // the counts or the convergence differ and the twin's did not: the common iterations are
    // still compared below, because the run may have followed the reference to the last common
    // iteration and parted only at the convergence test (RMSBL < EPS1) where the reference
    // itself moves by more than `DIVERGENCE_FLOOR` — the divergent rule applied to that
    // threshold; otherwise the count is a mismatch, settled after the loop
    let count_differs = !count_flip && (p.iterations != report.iterations.1 || p.converged != report.converged.1);
    let its = &rec.iters[&k];
    if report.iterations.1 != its.len() {
        let why = format!(
            "{ctx}: per-iteration record length {} vs NITDONE {}",
            its.len(),
            report.iterations.1
        );
        mismatch(&mut report, why);
        return report;
    }

    let mut worst_ratio = 0.0_f64;
    let names = ["RMSBL", "RLX", "CL", "CD", "CM", "ALFA", "MINF", "REINF"];
    let mut tdiff = [0.0_f64; 8];
    for (y, r) in p.iteration_records.iter().zip(its) {
        let ictx = format!("{ctx} iteration {}", y.iteration);
        let fi = cf.and_then(|c| c.iterations.get(y.iteration - 1));
        let floor_rmsbl = fl(fi, "RMSBL");
        println!(
            "    iteration {:2}: RMSBL yfoil {:.6e} xfoil {:.6e} |diff| {:.2e} (1-ULP floor {:.2e}) | RLX diff {:.2e} (floor {:.2e}) | CL diff {:.2e} (floor {:.2e}) | IST {}",
            y.iteration,
            y.residual,
            r[1],
            (y.residual - r[1]).abs(),
            floor_rmsbl,
            (y.relaxation - r[2]).abs(),
            fl(fi, "RLX"),
            (y.cl - r[3]).abs(),
            fl(fi, "CL"),
            y.i_stagnation_node
        );

        // UPDATE's reported limiter (call 1 only): the variable and station of the largest change
        let mut limiter_flip: Option<String> = None;
        if k == 1 {
            if let Some(b) = rec.iter1.get(y.iteration - 1) {
                let (xv, xi, xs) = (
                    &b["VMXBL"],
                    b["IMXBL"].parse::<usize>().unwrap(),
                    b["ISMXBL"].parse::<usize>().unwrap(),
                );
                let xr: f64 = b["RMXBL"].parse().unwrap();
                println!(
                    "        RLX limiter: yfoil {}@({},{}) RMXBL {:+.6e} | xfoil {}@({},{}) RMXBL {:+.6e}",
                    y.residual_max_variable,
                    y.i_residual_max_station,
                    y.residual_max_side,
                    y.residual_max,
                    xv,
                    xi,
                    xs,
                    xr
                );
                if &y.residual_max_variable.to_string() != xv
                    || (y.i_residual_max_station, y.residual_max_side) != (xi, xs)
                {
                    // a tie: the twin itself reports a different limiter here, or the two
                    // reported largest changes agree within the reference's own RMXBL spread
                    let twin_flipped = fl(fi, "LIMITER_FLIP") > 0.5;
                    let tie = twin_flipped
                        || (y.residual_max.abs() - xr.abs()).abs()
                            <= allowed(y.residual_max.abs(), xr.abs(), transient_tol, 1.0, fl(fi, "RMXBL"));
                    let msg = format!(
                        "reported RLX limiter yfoil {}@({},{}) |RMXBL| {:.6e} vs xfoil {}@({},{}) |RMXBL| {:.6e}",
                        y.residual_max_variable,
                        y.i_residual_max_station,
                        y.residual_max_side,
                        y.residual_max.abs(),
                        xv,
                        xi,
                        xs,
                        xr.abs()
                    );
                    if tie {
                        println!(
                            "        (tie: {msg} — two near-equal changes; RLX unaffected; twin flipped too: {twin_flipped}, RMXBL floor {:.2e})",
                            fl(fi, "RMXBL")
                        );
                    } else {
                        limiter_flip = Some(msg);
                    }
                }
            }
        }

        // every gated transient of this iteration, in the order they are reported
        let gated = [
            ("RMSBL", y.residual, r[1], 1.0),
            ("RLX", y.relaxation, r[2], 1.0),
            ("CL", y.cl, r[3], 1.0),
            ("CD", y.cd, r[4], 1.0),
            ("CM", y.cm, r[5], 1.0),
            ("ALFA", y.alpha, r[10], 1.0),
            ("MINF", y.mach, r[11], 1.0),
            ("REINF", y.re, r[12], r[12]),
        ];
        let parted_value = gated.iter().find(|(name, ours, theirs, scale)| {
            (ours - theirs).abs() > allowed(*ours, *theirs, transient_tol, *scale, fl(fi, name))
        });
        let trace_ok =
            y.i_stagnation_node == r[7] as usize && y.i_transition_station[1..] == [r[8] as usize, r[9] as usize];
        if (parted_value.is_some() || !trace_ok || limiter_flip.is_some()) && floor_rmsbl > DIVERGENCE_FLOOR {
            let parted = if let Some(m) = limiter_flip.clone() {
                m
            } else if !trace_ok {
                format!(
                    "IST/ITRAN yfoil {}/{}/{} vs xfoil {}/{}/{}",
                    y.i_stagnation_node, y.i_transition_station[1], y.i_transition_station[2], r[7], r[8], r[9]
                )
            } else {
                let (name, ours, theirs, _) = parted_value.unwrap();
                format!("{name} yfoil {ours:.6e} vs xfoil {theirs:.6e}")
            };
            let why = format!(
                "every earlier iteration matched within {FLOOR_FACTOR}× the reference's own 1-ULP floor; at iteration {} the reference itself moves by {floor_rmsbl:.2e} (> {DIVERGENCE_FLOOR:.0e}) under 1 ULP and the runs part ({parted})",
                y.iteration
            );
            println!("{ictx}: DIVERGENT — {why}. Values from here on are not gated; the one-step replay from XFOIL's state is the evidence that the step itself is faithful.");
            report.outcome = Outcome::Divergent {
                at_iteration: y.iteration,
                why,
            };
            break;
        }
        for (m, (name, ours, theirs, scale)) in gated.iter().enumerate() {
            let floor = fl(fi, name);
            if let Err(why) = value_ok(*ours, *theirs, transient_tol, *scale, floor, &format!("{ictx}: {name}")) {
                mismatch(&mut report, why);
                break;
            }
            if floor > 0.0 {
                worst_ratio = worst_ratio.max((ours - theirs).abs() / floor);
            }
            tdiff[m] = tdiff[m].max((ours - theirs).abs());
        }
        if !matches!(report.outcome, Outcome::Match) {
            break;
        }
        if y.i_stagnation_node != r[7] as usize {
            mismatch(
                &mut report,
                format!("{ictx}: IST yfoil {} vs xfoil {}", y.i_stagnation_node, r[7]),
            );
            break;
        }
        if y.i_transition_station[1..] != [r[8] as usize, r[9] as usize] {
            mismatch(
                &mut report,
                format!(
                    "{ictx}: ITRAN yfoil {:?} vs xfoil {:?}",
                    &y.i_transition_station[1..],
                    [r[8] as usize, r[9] as usize]
                ),
            );
            break;
        }
        if let Some(flip) = limiter_flip {
            mismatch(
                &mut report,
                format!(
                    "{ictx}: {flip} (RLX values agree but the largest changes are not tied; floor {floor_rmsbl:.2e})"
                ),
            );
            break;
        }
        report.gated_iterations = y.iteration;
    }
    report.worst_floor_ratio = worst_ratio;
    report.transient_diff = names.iter().zip(tdiff).map(|(n, d)| (*n, d)).collect();
    if count_differs && matches!(report.outcome, Outcome::Match) {
        // every common iteration matched: the runs part at the convergence test of the last one
        let common = p.iterations.min(report.iterations.1);
        let floor_rmsbl = fl(cf.and_then(|c| c.iterations.get(common - 1)), "RMSBL");
        let counts = format!(
            "iteration count yfoil {} (converged {}) vs xfoil {} (converged {})",
            p.iterations, p.converged, report.iterations.1, report.converged.1
        );
        if floor_rmsbl > DIVERGENCE_FLOOR {
            let why = format!(
                "the common {common} iterations matched within {FLOOR_FACTOR}× the reference's own 1-ULP floor; at iteration {common} the reference itself moves by {floor_rmsbl:.2e} (> {DIVERGENCE_FLOOR:.0e}) under 1 ULP and the runs part at the convergence test there ({counts})"
            );
            if !ungated_flag {
                println!("{ctx}: DIVERGENT — {why}");
            }
            report.outcome = Outcome::Divergent {
                at_iteration: common,
                why,
            };
        } else {
            mismatch(&mut report, format!("{ctx}: {counts} (the common iterations matched; RMSBL floor at iteration {common} {floor_rmsbl:.2e})"));
        }
    }
    if let Some(why) = ungated {
        // whatever the values did, the call is divergent from iteration 0; what they did is
        // kept in the reason and in `gated_iterations`
        let note = match &report.outcome {
            Outcome::Match => "yfoil followed the reference to the end of the call anyway".to_string(),
            Outcome::Divergent { at_iteration, why: w } => {
                format!("the runs would part at iteration {at_iteration}: {w}")
            }
            Outcome::Mismatch { why: w } => format!("the runs would part: {w}"),
        };
        println!("{ctx}: (ungated) {note}");
        report.outcome = Outcome::Divergent {
            at_iteration: 0,
            why: format!("{why}; {note}"),
        };
        return report;
    }
    if !matches!(report.outcome, Outcome::Match) {
        return report;
    }
    if p.iterations != report.iterations.1 || p.converged != report.converged.1 {
        // count_flip: the common iterations matched; the runs part where one of them stops
        let at = p.iterations.min(report.iterations.1) + 1;
        let why = format!(
            "one of the reference's own 1-ULP twins changes this call's iteration count / convergence ({}); yfoil {} iterations (converged {}), reference {} (converged {}), the common {} iterations match",
            within_call.join("; "),
            p.iterations,
            p.converged,
            report.iterations.1,
            report.converged.1,
            at - 1
        );
        println!("{ctx}: DIVERGENT — {why}");
        report.outcome = Outcome::Divergent { at_iteration: at, why };
        return report;
    }

    if !within_call.is_empty() {
        let why = format!(
            "yfoil matched the reference to the end of the call, but one of the reference's own 1-ULP twins changes its branch trace inside it: {}",
            within_call.join("; ")
        );
        println!("{ctx}: DIVERGENT — {why}");
        report.outcome = Outcome::Divergent { at_iteration: 0, why };
        return report;
    }

    let fp = cf.map(|c| &c.point);
    for (name, ours, scale) in [
        ("ALFA", p.alpha, 1.0),
        ("CL", p.cl, 1.0),
        ("CM", p.cm, 1.0),
        ("CD", p.cd, 1.0),
        ("CDF", p.cd_friction, 1.0),
        ("CDP", p.cd_pressure, 1.0),
        ("XOCTR1", p.transition_upper[0], 1.0),
        ("XOCTR2", p.transition_lower[0], 1.0),
        ("MINF", st.mach, 1.0),
        ("REINF", st.re, st.re),
    ] {
        let theirs: f64 = x[name].parse().unwrap();
        if let Err(why) = value_ok(ours, theirs, TOL_SOLVER, scale, fl(fp, name), &format!("{ctx}: {name}")) {
            mismatch(&mut report, why);
            return report;
        }
        report.point_diff.push((name, (ours - theirs).abs()));
    }
    if st.i_stagnation_node != x["IST"].parse::<usize>().unwrap() {
        mismatch(
            &mut report,
            format!("{ctx}: IST yfoil {} vs xfoil {}", st.i_stagnation_node, x["IST"]),
        );
        return report;
    }
    let itran = [
        x["ITRAN1"].parse::<usize>().unwrap(),
        x["ITRAN2"].parse::<usize>().unwrap(),
    ];
    if p.i_transition_station[1..] != itran {
        mismatch(
            &mut report,
            format!(
                "{ctx}: ITRAN yfoil {:?} vs xfoil {itran:?}",
                &p.i_transition_station[1..]
            ),
        );
        return report;
    }
    println!(
        "{ctx}: {:2} iterations, converged={} CL {:.8} CD {:.8} CM {:+.8} XTR {:.5}/{:.5} Re {:.0} M {:.3} — match (worst diff/floor {:.2})",
        p.iterations, p.converged, p.cl, p.cd, p.cm, p.transition_upper[0], p.transition_lower[0], st.re, st.mach, worst_ratio
    );
    report
}

/// Compare an inviscid-only call — SPECAL (`kind = "specal"`, one per OPER ALFA) or SPECCL
/// (`"speccl"`, one per OPER CL) — with a `Session` result: SPECCL's exit iteration ITAL is the
/// branch trace, ALFA/CL/CM/CDP and the per-node GAM/QINV/CPI the values, gated at `TOL_SOLVER`
/// or the twin's own spread × `FLOOR_FACTOR`.
pub fn compare_inviscid_call(rec: &Records, kind: &str, k: usize, p: &PointResult, st: &SolverState) -> CallReport {
    let (points, floors) = match kind {
        "specal" => (&rec.specal, rec.floor.as_ref().map(|f| &f.specal_calls)),
        "speccl" => (&rec.speccl, rec.floor.as_ref().map(|f| &f.speccl_calls)),
        _ => panic!("kind is specal or speccl"),
    };
    let x = &points[k - 1];
    let ctx = format!("{kind} call {k} (alpha {:.3}°)", p.alpha.to_degrees());
    let fp = floors.and_then(|f| f.get(k - 1));
    let fl = |name: &str| fp.and_then(|h| h.get(name)).copied().unwrap_or(0.0);
    let xf = |key: &str| x.get(key).and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN);
    let mut report = CallReport {
        call: k,
        alpha_deg: p.alpha.to_degrees(),
        iterations: (
            p.inviscid_cl_iterations,
            x.get("ITAL").and_then(|v| v.parse().ok()).unwrap_or(0),
        ),
        converged: (true, true),
        yfoil: [p.cl, 0.0, p.cm, 0.0],
        xfoil: [xf("CL"), 0.0, xf("CM"), 0.0],
        yfoil_transition: [f64::NAN, f64::NAN],
        xfoil_transition: [f64::NAN, f64::NAN],
        ..Default::default()
    };
    let mismatch = |report: &mut CallReport, why: String| {
        println!("{ctx}: MISMATCH — {why}");
        report.outcome = Outcome::Mismatch { why };
    };
    if kind == "speccl" && report.iterations.0 != report.iterations.1 {
        let why = format!(
            "{ctx}: ITAL yfoil {} vs xfoil {}",
            report.iterations.0, report.iterations.1
        );
        mismatch(&mut report, why);
        return report;
    }
    for (name, ours) in [("ALFA", p.alpha), ("CL", p.cl), ("CM", p.cm), ("CDP", p.cd_pressure)] {
        let theirs: f64 = x[name].parse().unwrap();
        if let Err(why) = value_ok(ours, theirs, TOL_SOLVER, 1.0, fl(name), &format!("{ctx}: {name}")) {
            mismatch(&mut report, why);
            return report;
        }
        report.point_diff.push((name, (ours - theirs).abs()));
        report.worst_floor_ratio = report.worst_floor_ratio.max(if fl(name) > 0.0 {
            (ours - theirs).abs() / fl(name)
        } else {
            0.0
        });
    }
    let n: usize = x["N"].parse().unwrap();
    let mut worst = [0.0_f64; 3];
    for i in 1..=n {
        let row: Vec<f64> = x[&format!("NODE({i:5})")]
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        let ours = [st.gamma[i], st.q_inviscid[i], st.cp_inviscid[i]];
        for (m, name) in ["GAM", "QINV", "CPI"].iter().enumerate() {
            if let Err(why) = value_ok(
                ours[m],
                row[m],
                TOL_SOLVER,
                1.0,
                fl(name),
                &format!("{ctx}: {name}({i})"),
            ) {
                mismatch(&mut report, why);
                return report;
            }
            worst[m] = worst[m].max((ours[m] - row[m]).abs());
        }
    }
    for (m, name) in ["GAM", "QINV", "CPI"].iter().enumerate() {
        report.point_diff.push((name, worst[m]));
    }
    println!(
        "{ctx}: ALFA {:.8}° CL {:.8} CM {:+.8} ITAL {} — match (worst node diff GAM {:.2e} QINV {:.2e} CPI {:.2e})",
        p.alpha.to_degrees(),
        p.cl,
        p.cm,
        report.iterations.1,
        worst[0],
        worst[1],
        worst[2]
    );
    report
}

/// Base tolerance for per-iteration transients of call `k`: a fresh single point is at the
/// solver floor; later calls in a sequence carry the propagated wake/DIJ floor.
pub fn transient_tol(k: usize) -> f64 {
    if k == 1 {
        TOL_SOLVER
    } else {
        TOL_TRANSIENT
    }
}
