//! Comparison of a `Session` run against the VISCAL-level reference records every case keeps:
//! `viscal_points.dat` (one block per VISCAL call), `viscal_iters_all.dat` (one `IT` line per
//! iteration: `k i RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2 ALFA MINF REINF`),
//! `viscal_iter.dat` (call 1 only, with UPDATE's RMXBL/VMXBL/IMXBL/ISMXBL) and, when the pipeline
//! recorded it, `noise_floor.json` — the same run on the +1-ULP twin geometry.
//!
//! Exact: iteration count, LVCONV, IST, ITRAN. Values: `|a − b| ≤ max(tol·max(|a|,|b|,scale),
//! FLOOR_FACTOR · floor)` where `floor` is the reference's own 1-ULP spread of that very value
//! (CLAUDE.md Rule 1: tolerances are the measured floor times a safety factor).
//!
//! Third outcome — **threshold-straddling** — reported, never silently passed or failed:
//! - the +1-ULP twin itself changes the branch trace (iteration count, convergence, IST, ITRAN):
//!   the whole case is straddling;
//! - every earlier iteration matched, and at this iteration the reference's own 1-ULP spread
//!   exceeds `STRADDLE_FLOOR`: the runs are allowed to part here (the one-step replay from
//!   XFOIL's exact state is the evidence that the step itself is faithful);
//! - UPDATE's *reported* limiter (VMXBL/IMXBL: the largest normalised change) differs while
//!   |RMXBL| agrees within the floor: two near-tied changes, RLX itself unaffected — a tie.

#![allow(dead_code)]

use super::tolerances::{FLOOR_FACTOR, STRADDLE_FLOOR, TOL_SOLVER, TOL_TRANSIENT};
use std::collections::HashMap;
use std::path::Path;
use yfoil::solver::analysis::PointResult;
use yfoil::solver::blstate::SolverState;

#[derive(Debug, Clone, Default)]
pub struct CallFloor {
    pub point: HashMap<String, f64>,
    pub iterations: Vec<HashMap<String, f64>>,
}

#[derive(Debug, Clone, Default)]
pub struct Floor {
    pub branch_identical: bool,
    pub flips: Vec<String>,
    pub calls: Vec<CallFloor>,
    /// per dumped SETBL/UPDATE call: max |base − twin| per post-UPDATE array (XSSI UEDG THET DSTR CTAU MASS)
    pub update_output: HashMap<usize, HashMap<String, f64>>,
}

pub struct Records {
    pub points: Vec<HashMap<String, String>>,
    /// per call (1-based key): rows of the IT line values after the call number
    pub iters: HashMap<usize, Vec<Vec<f64>>>,
    /// VISCAL call 1's per-iteration blocks of viscal_iter.dat (RMXBL/VMXBL/IMXBL/ISMXBL)
    pub iter1: Vec<HashMap<String, String>>,
    pub floor: Option<Floor>,
}

/// What a call comparison concluded.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Match,
    /// Reported separately: the reason, and the iteration at which the runs were allowed to part
    /// (0 = the whole case, from the twin's own branch flips)
    Straddling {
        at_iteration: usize,
        why: String,
    },
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
    let points = kv_blocks(
        &std::fs::read_to_string(case_dir.join("viscal_points.dat")).unwrap(),
        "CALL",
    );
    let mut iters: HashMap<usize, Vec<Vec<f64>>> = HashMap::new();
    for l in std::fs::read_to_string(case_dir.join("viscal_iters_all.dat"))
        .unwrap()
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
                    .map(|o| o.iter().map(|(k, x)| (k.clone(), x.as_f64().unwrap())).collect())
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
                            })
                            .collect()
                    })
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
        iters,
        iter1,
        floor,
    }
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

/// Compare VISCAL call `k` (1-based) of the reference with a `Session` result.
/// `transient_tol` is the base tolerance for the per-iteration values.
pub fn check_call(rec: &Records, k: usize, p: &PointResult, st: &SolverState, transient_tol: f64) -> Outcome {
    let x = &rec.points[k - 1];
    let ctx = format!("call {k} (alpha {:.3}°)", p.alpha.to_degrees());
    assert_eq!(k, x["CALL"].parse::<usize>().unwrap());

    if let Some(f) = &rec.floor {
        if !f.branch_identical {
            let why = format!(
                "the reference's own +1-ULP twin changes its branch trace: {}",
                f.flips.join("; ")
            );
            println!(
                "{ctx}: THRESHOLD-STRADDLING — {why}. yfoil: {} iterations, converged={}, CL {:.8}; reference: {} iterations, converged={}, CL {}",
                p.iterations, p.converged, p.cl, x["NITDONE"], x["LVCONV"], x["CL"]
            );
            return Outcome::Straddling { at_iteration: 0, why };
        }
    }
    let cf = rec.floor.as_ref().and_then(|f| f.calls.get(k - 1));
    let fl = |m: Option<&HashMap<String, f64>>, name: &str| m.and_then(|h| h.get(name)).copied().unwrap_or(0.0);

    assert_eq!(
        p.iterations,
        x["NITDONE"].parse::<usize>().unwrap(),
        "{ctx}: iteration count"
    );
    assert_eq!(p.converged, x["LVCONV"] == "T", "{ctx}: LVCONV");
    let its = &rec.iters[&k];
    assert_eq!(p.iterations, its.len(), "{ctx}: per-iteration record length");

    let mut worst_ratio = 0.0_f64;
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

        let rmsbl_ok = (y.residual - r[1]).abs() <= allowed(y.residual, r[1], transient_tol, 1.0, floor_rmsbl);
        if (!rmsbl_ok || limiter_flip.is_some()) && floor_rmsbl > STRADDLE_FLOOR {
            let why = format!(
                "every earlier iteration matched within {FLOOR_FACTOR}× the reference's own 1-ULP floor; at iteration {} the reference itself moves by {floor_rmsbl:.2e} (> {STRADDLE_FLOOR:.0e}) under 1 ULP and the runs part ({})",
                y.iteration,
                limiter_flip.clone().unwrap_or_else(|| format!("RMSBL yfoil {:.6e} vs xfoil {:.6e}", y.residual, r[1]))
            );
            println!("{ictx}: THRESHOLD-STRADDLING — {why}. Values from here on are not gated; the one-step replay from XFOIL's state is the evidence that the step itself is faithful.");
            return Outcome::Straddling {
                at_iteration: y.iteration,
                why,
            };
        }
        for (name, ours, theirs, scale) in [
            ("RMSBL", y.residual, r[1], 1.0),
            ("RLX", y.relaxation, r[2], 1.0),
            ("CL", y.cl, r[3], 1.0),
            ("CD", y.cd, r[4], 1.0),
            ("CM", y.cm, r[5], 1.0),
            ("ALFA", y.alpha, r[10], 1.0),
            ("MINF", y.mach, r[11], 1.0),
            ("REINF", y.re, r[12], r[12]),
        ] {
            let floor = fl(fi, name);
            assert_value(ours, theirs, transient_tol, scale, floor, &format!("{ictx}: {name}"));
            if floor > 0.0 {
                worst_ratio = worst_ratio.max((ours - theirs).abs() / floor);
            }
        }
        assert_eq!(y.i_stagnation_node, r[7] as usize, "{ictx}: IST");
        assert_eq!(
            y.i_transition_station[1..],
            [r[8] as usize, r[9] as usize],
            "{ictx}: ITRAN"
        );
        if let Some(flip) = limiter_flip {
            panic!("{ictx}: {flip} (RLX values agree but the largest changes are not tied; floor {floor_rmsbl:.2e})");
        }
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
        assert_value(
            ours,
            x[name].parse().unwrap(),
            TOL_SOLVER,
            scale,
            fl(fp, name),
            &format!("{ctx}: {name}"),
        );
    }
    assert_eq!(st.i_stagnation_node, x["IST"].parse::<usize>().unwrap(), "{ctx}: IST");
    assert_eq!(
        p.i_transition_station[1..],
        [
            x["ITRAN1"].parse::<usize>().unwrap(),
            x["ITRAN2"].parse::<usize>().unwrap()
        ],
        "{ctx}: ITRAN"
    );
    println!(
        "{ctx}: {:2} iterations, converged={} CL {:.8} CD {:.8} CM {:+.8} XTR {:.5}/{:.5} Re {:.0} M {:.3} — match (worst diff/floor {:.2})",
        p.iterations, p.converged, p.cl, p.cd, p.cm, p.transition_upper[0], p.transition_lower[0], st.re, st.mach, worst_ratio
    );
    Outcome::Match
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
