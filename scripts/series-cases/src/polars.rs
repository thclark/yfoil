//! Full ±30° polars of the study's sections, yFoil against the reference point by point.
//!
//! Each `branch-coverage-polar` case is XFOIL's polar procedure (`ALFA 0 / ASEQ 1 30 1 / INIT /
//! ALFA 0 / ASEQ -1 -30 -1`, ITMAX + 5 per ASEQ point, the sequence halting after NSEQEX = 4
//! consecutive unconverged points). yFoil is driven the same way through its `Session`, applying
//! the same halting rule to its own convergence, and every point the reference recorded is
//! compared with `compare_call` (the whole-run gate), matched by the point it computed (leg and
//! alpha), never by running index: a sequence may halt at a different point in the two codes.
//! The reference may have fewer points than the script: its sequence halted, or the fixture
//! driver's watchdog ended a run that had gone non-finite and hung in XFOIL's plot-label loop.
//! yFoil's sweep continues to its own halt either way, so the figure shows where each code ends
//! up.

use crate::utilities::records::{
    compare_call, compare_state, load, reference_ill_conditioned, transient_tol, ArrayDiff, Outcome,
};
use crate::Case;
use serde_json::{json, Value};
use std::path::Path;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::Session;

/// NSEQEX: XFOIL's ASEQ halts once this many consecutive points fail to converge.
pub const NSEQEX: usize = 4;

pub struct PolarPoint {
    pub leg: usize,
    pub alpha_deg: f64,
    /// CL, CD, CM, last RMSBL, iterations, converged
    pub yfoil: (f64, f64, f64, f64, usize, bool),
    /// the same from the reference, when it recorded the call
    pub xfoil: Option<(f64, f64, f64, f64, usize, bool)>,
    /// match | divergent | mismatch | ungated (the leg had already departed), when the call was
    /// compared
    pub outcome: Option<String>,
    /// whether the reference itself was ill-conditioned at this point, from its twins
    /// (`reference_ill_conditioned`); `None` without twins or a matching call
    pub xfoil_ill_conditioned: Option<bool>,
    pub parted_iteration: Option<usize>,
    /// iterations of the call compared within tolerance (all of them when yFoil followed the
    /// reference to the end of the call, whatever the twin did)
    pub gated_iterations: Option<usize>,
    /// a non-finite value (RMSBL, CL, CD or CM) in any iteration of the point, or in the point
    /// itself — yFoil's run, and the reference's records
    pub yfoil_nonfinite: bool,
    pub xfoil_nonfinite: Option<bool>,
    /// yFoil's station Newton failures with a non-finite residual in the point's marches, and
    /// the `NaN` tokens the reference printed during the call
    pub yfoil_nonfinite_stations: usize,
    pub xfoil_nan_tokens: Option<usize>,
    /// transition x/c (upper, lower) and IST/ITRAN(1)/ITRAN(2) at the end of the point
    pub yfoil_transition: [f64; 2],
    pub xfoil_transition: Option<[f64; 2]>,
    pub yfoil_stations: [usize; 3],
    pub xfoil_stations: Option<[usize; 3]>,
    /// yFoil's final state against the reference's `viscal_state_<k>.dat`,
    /// array by array, worst ratio to the gate first
    pub state: Option<Vec<ArrayDiff>>,
    /// where transition sits within its station interval at the end of the point, per side
    /// (index 1 upper, 2 lower; slot 0 unused), from yFoil's state
    pub transition_position: [Option<TransitionPosition>; 3],
}

/// XFOIL's transition convention: ITRAN is the first turbulent station and the transition arc
/// length XSSITR lies in the interval between stations ITRAN − 1 and ITRAN. `fraction` is where
/// in that interval it sits — 0 at station ITRAN − 1, 1 at station ITRAN — and `nearest` the
/// distance to the nearer of the two stations as a fraction of the interval (0 to 0.5).
#[derive(Debug, Clone, Copy)]
pub struct TransitionPosition {
    pub xi_transition: f64,
    pub i_station: usize,
    pub xi_before: f64,
    pub xi_at: f64,
    pub fraction: f64,
    pub nearest: f64,
    /// ITRAN is the trailing-edge station (no transition on the airfoil: it happens in the wake
    /// or at the TE) — the fraction is then not a position between two airfoil stations
    pub at_te: bool,
}

impl TransitionPosition {
    pub fn from_state(st: &yfoil::solver::blstate::SolverState, side: usize) -> Option<Self> {
        let ibl = st.i_transition_station[side];
        if ibl < 2 || ibl > st.n_stations[side] {
            return None;
        }
        let (xi_before, xi_at) = (st.xi[side][ibl - 1], st.xi[side][ibl]);
        let xi_transition = st.xi_transition[side];
        let fraction = (xi_transition - xi_before) / (xi_at - xi_before);
        Some(Self {
            xi_transition,
            i_station: ibl,
            xi_before,
            xi_at,
            fraction,
            nearest: fraction.min(1.0 - fraction),
            at_te: ibl >= st.i_te_station[side],
        })
    }
    pub fn to_json(self) -> Value {
        json!({
            "xi_transition": self.xi_transition, "i_station": self.i_station, "xi_before": self.xi_before,
            "xi_at": self.xi_at, "fraction": self.fraction, "nearest": self.nearest, "at_te": self.at_te,
        })
    }
}

pub struct PolarRun {
    pub name: String,
    pub section: String,
    pub script: String,
    pub description: String,
    pub points: Vec<PolarPoint>,
    /// calls the reference recorded (its sequences' lengths)
    pub xfoil_calls: usize,
    pub truncated_by_watchdog: bool,
}

impl PolarRun {
    /// the section family the run belongs to, as a file-name slug (`NACA 63-415` → `naca-63-415`)
    pub fn family(&self) -> String {
        self.section
            .to_lowercase()
            .replace(" (sharp te)", "")
            .replace(' ', "-")
            .replace(['–', '—'], "-")
    }
    /// yFoil followed the reference to the end of this call: every iteration within tolerance
    /// (a match, or divergent only because the reference's own twin flipped its trace)
    fn followed(p: &PolarPoint) -> bool {
        match (&p.xfoil, p.gated_iterations) {
            (Some(x), Some(g)) => g == x.4,
            _ => false,
        }
    }
    /// the first recorded point yFoil did not follow to the end, if any: (leg, alpha)
    pub fn parted(&self) -> Option<(usize, f64)> {
        self.points
            .iter()
            .find(|p| p.xfoil.is_some() && !Self::followed(p))
            .map(|p| (p.leg, p.alpha_deg))
    }
    /// recorded points yFoil followed to the end
    pub fn matched(&self) -> usize {
        self.points.iter().filter(|p| Self::followed(p)).count()
    }
    pub fn to_json(&self) -> Value {
        let pt = |v: &(f64, f64, f64, f64, usize, bool)| json!({ "cl": v.0, "cd": v.1, "cm": v.2, "rmsbl": v.3, "iterations": v.4, "converged": v.5 });
        json!({
            "name": self.name, "section": self.section, "family": self.family(), "script": self.script, "description": self.description,
            "xfoil_calls": self.xfoil_calls, "yfoil_points": self.points.len(), "matched": self.matched(),
            "truncated_by_watchdog": self.truncated_by_watchdog, "parted": self.parted(),
            "points": self.points.iter().map(|p| json!({
                "leg": p.leg, "alpha_deg": p.alpha_deg, "yfoil": pt(&p.yfoil),
                "xfoil": p.xfoil.as_ref().map(pt), "outcome": p.outcome, "xfoil_ill_conditioned": p.xfoil_ill_conditioned, "parted_iteration": p.parted_iteration,
                "gated_iterations": p.gated_iterations, "followed": Self::followed(p),
                "yfoil_nonfinite": p.yfoil_nonfinite, "xfoil_nonfinite": p.xfoil_nonfinite,
                "yfoil_nonfinite_stations": p.yfoil_nonfinite_stations, "xfoil_nan_tokens": p.xfoil_nan_tokens,
                "yfoil_transition": p.yfoil_transition, "xfoil_transition": p.xfoil_transition,
                "yfoil_stations": p.yfoil_stations, "xfoil_stations": p.xfoil_stations,
                "stations_differ": p.xfoil_stations.is_some_and(|x| x != p.yfoil_stations),
                "transition_position": [p.transition_position[1].map(TransitionPosition::to_json), p.transition_position[2].map(TransitionPosition::to_json)],
                "state": p.state.as_ref().map(|arrays| arrays.iter().map(|d| json!({
                    "array": d.name, "diff": d.diff, "index": d.index, "floor": d.floor, "ratio": d.ratio,
                    "yfoil": d.yfoil, "xfoil": d.xfoil,
                })).collect::<Vec<_>>()),
            })).collect::<Vec<_>>(),
        })
    }
}

/// One operating point of a series case.
enum Op {
    Alpha(f64),
    Cl(f64),
}

/// The reference's SPECAL records (`specal_points.dat`): per call CL, CM, CDP and the per-node
/// GAM, QINV, CPI.
fn specal_records(dir: &Path) -> Vec<(f64, f64, f64, Vec<[f64; 3]>)> {
    let text = std::fs::read_to_string(dir.join("specal_points.dat")).unwrap_or_default();
    let mut out: Vec<(f64, f64, f64, Vec<[f64; 3]>)> = vec![];
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        let k = k.trim();
        if k == "CALL" {
            out.push((f64::NAN, f64::NAN, f64::NAN, vec![]));
            continue;
        }
        let Some(cur) = out.last_mut() else { continue };
        let num = |s: &str| s.trim().parse::<f64>().unwrap_or(f64::NAN);
        match k {
            "CL" => cur.0 = num(v),
            "CM" => cur.1 = num(v),
            "CDP" => cur.2 = num(v),
            _ if k.starts_with("NODE(") => {
                let t: Vec<f64> = v.split_whitespace().map(num).collect();
                if t.len() >= 3 {
                    cur.3.push([t[0], t[1], t[2]]);
                }
            }
            _ => {}
        }
    }
    out
}

/// An inviscid series case: one SPECAL per alpha, from the panels, compared with the reference's
/// SPECAL record of the same call — CL, CM and CDp, and every node's GAM, QINV and CPI — on the
/// Rule 1 metric (scale 1) within `max(TOL_SOLVER · scale, FLOOR_FACTOR · spread)`, with *spread*
/// the reference's own twins' largest spread of that value at that call (`noise_floor.json`,
/// `specal_calls`); the reference is ill-conditioned at a point where a scalar's spread exceeds
/// `DIVERGENCE_FLOOR`. Inviscid runs have no iterations, transition or boundary layer; those fields
/// are left empty.
pub fn inviscid_sweep(root: &Path, case: &Case, description: &str) -> PolarRun {
    use crate::utilities::records::{DIVERGENCE_FLOOR, FLOOR_FACTOR};
    use crate::utilities::tolerances::TOL_SOLVER;
    let dir = root.join("target/fixtures").join(&case.name);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    let recs = specal_records(&dir);
    // the reference's twins: per call, the largest spread of each scalar and array
    let spreads: Vec<serde_json::Map<String, Value>> = std::fs::read_to_string(dir.join("noise_floor.json"))
        .ok()
        .and_then(|t| serde_json::from_str::<Value>(&t).ok())
        .and_then(|v| v["specal_calls"].as_array().cloned())
        .map(|a| {
            a.iter()
                .map(|c| c["point"].as_object().cloned().unwrap_or_default())
                .collect()
        })
        .unwrap_or_default();
    let mut session = Session::new(&airfoil, case.conditions());
    let mut points = vec![];
    for (j, a) in case.alphas.iter().enumerate() {
        let p = session.alpha(a.to_radians());
        let x = recs.get(j);
        let spread = |key: &str| -> f64 {
            spreads
                .get(j)
                .and_then(|m| m.get(key))
                .map_or(0.0, |v| v.as_f64().unwrap_or(f64::INFINITY))
        };
        let within = |a: f64, b: f64, key: &str| -> bool {
            (a.is_nan() && b.is_nan())
                || (a - b).abs() <= (TOL_SOLVER * a.abs().max(b.abs()).max(1.0)).max(FLOOR_FACTOR * spread(key))
        };
        let outcome = x.map(|x| {
            let nodes_ok = {
                let st = session.state();
                x.3.iter().enumerate().all(|(i, n)| {
                    let i = i + 1;
                    within(st.gamma.get(i).copied().unwrap_or(f64::NAN), n[0], "GAM")
                        && within(st.q_inviscid.get(i).copied().unwrap_or(f64::NAN), n[1], "QINV")
                        && within(st.cp_inviscid.get(i).copied().unwrap_or(f64::NAN), n[2], "CPI")
                })
            };
            if within(p.cl, x.0, "CL") && within(p.cm, x.1, "CM") && within(p.cd_pressure, x.2, "CDP") && nodes_ok {
                "match".to_string()
            } else {
                "mismatch".to_string()
            }
        });
        let ill = (!spreads.is_empty()).then(|| ["CL", "CM", "CDP"].iter().any(|k| !(spread(k) <= DIVERGENCE_FLOOR)));
        points.push(PolarPoint {
            leg: 1,
            alpha_deg: *a,
            yfoil: (p.cl, p.cd_pressure, p.cm, f64::NAN, 0, true),
            xfoil: x.map(|x| (x.0, x.2, x.1, f64::NAN, 0, true)),
            gated_iterations: (outcome.as_deref() == Some("match")).then_some(0),
            outcome,
            xfoil_ill_conditioned: x.and(ill),
            parted_iteration: None,
            yfoil_nonfinite: ![p.cl, p.cm, p.cd_pressure].iter().all(|v| v.is_finite()),
            xfoil_nonfinite: x.map(|x| ![x.0, x.1, x.2].iter().all(|v| v.is_finite())),
            yfoil_nonfinite_stations: 0,
            xfoil_nan_tokens: None,
            yfoil_transition: [f64::NAN, f64::NAN],
            xfoil_transition: None,
            yfoil_stations: [session.state().i_stagnation_node, 0, 0],
            xfoil_stations: None,
            state: None,
            transition_position: [None, None, None],
        });
    }
    PolarRun {
        name: case.name.clone(),
        section: case.section(),
        script: case.script(),
        description: description.to_string(),
        points,
        xfoil_calls: recs.len(),
        truncated_by_watchdog: false,
    }
}

pub fn sweep(root: &Path, case: &Case, description: &str) -> PolarRun {
    if case.inviscid {
        return inviscid_sweep(root, case, description);
    }
    let dir = root.join("target/fixtures").join(&case.name);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    let rec = load(&dir);
    let manifest: Value = serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    // `viscal_nonfinite.dat`: NaN tokens the reference printed during each VISCAL call
    let xfoil_nan: Vec<usize> = std::fs::read_to_string(dir.join("viscal_nonfinite.dat"))
        .map(|t| {
            t.lines()
                .filter_map(|l| l.split_whitespace().last().and_then(|n| n.parse().ok()))
                .collect()
        })
        .unwrap_or_default();
    let truncated = manifest["truncated_by_watchdog"].as_bool().unwrap_or(false);
    // the reference's calls by (leg, alpha): a sequence may halt at a different point in the two
    // codes, so calls are matched by the point they computed, never by running index. The
    // reference's second leg starts where its recorded ALFA returns to the first alpha of the
    // script after INIT.
    let ref_alpha = |i: usize| -> f64 {
        rec.points[i]
            .get("ALFA")
            .and_then(|v| v.parse::<f64>().ok())
            .map(|r| r.to_degrees())
            .unwrap_or(f64::NAN)
    };
    let mut ref_calls: std::collections::HashMap<(usize, i64), usize> = std::collections::HashMap::new();
    {
        let mut leg_no = 1usize;
        let first_after = case.alphas_after_reinit.first().copied();
        for i in 0..rec.points.len() {
            let a = ref_alpha(i);
            if i > 0 && leg_no == 1 {
                if let Some(a0) = first_after {
                    // the first point of leg 2 repeats the leg's seed alpha after a sequence
                    // that was moving away from it
                    if (a - a0).abs() < 1e-6 && (ref_alpha(i - 1) - a0).abs() > 0.5 {
                        leg_no = 2;
                    }
                }
            }
            ref_calls.insert((leg_no, (a * 1000.0).round() as i64), i + 1);
        }
    }
    let mut session = Session::new(&airfoil, case.conditions());
    let mut points: Vec<PolarPoint> = vec![];
    // a leg is a list of operating points: the alphas of an ALFA / ASEQ polar leg (the first an
    // ALFA, the rest ASEQ points with its halting rule), or the CL points of a fixed-CL series
    // (each its own OPER `CL` command, matched to the reference's calls in order)
    let leg = |session: &mut Session, ops: &[Op], leg_no: usize, points: &mut Vec<PolarPoint>| {
        let mut failures = 0usize;
        // once a point of this leg was not followed to the end, every later point starts from a
        // state the reference did not compute: its comparison is informational, never a mismatch
        // (a leg starts afresh: INIT and a fresh march)
        let mut departed = false;
        for (j, op) in ops.iter().enumerate() {
            let p = match op {
                Op::Alpha(a) if j == 0 => session.alpha(a.to_radians()),
                Op::Alpha(a) => session.sequence_point(a.to_radians()),
                Op::Cl(c) => session.cl(*c),
            };
            let a = &match op {
                Op::Alpha(a) => *a,
                Op::Cl(_) => p.alpha.to_degrees(),
            };
            let k = match op {
                Op::Alpha(a) => ref_calls.get(&(leg_no, (a * 1000.0).round() as i64)).copied(),
                Op::Cl(_) => (j < rec.points.len()).then_some(j + 1),
            };
            let mut xfoil_transition = None;
            let mut xfoil_stations = None;
            let state = k.and_then(|k| compare_state(&rec, &dir, k, session.state()));
            let (xfoil, outcome, parted, gated) = if let Some(k) = k {
                let r = compare_call(&rec, k, &p, session.state(), transient_tol(k));
                let (o, at) = match &r.outcome {
                    Outcome::Match => ("match".to_string(), None),
                    Outcome::Divergent { at_iteration, .. } => ("divergent".to_string(), Some(*at_iteration)),
                    Outcome::Mismatch { .. } if departed => ("ungated".to_string(), Some(r.gated_iterations + 1)),
                    Outcome::Mismatch { .. } => ("mismatch".to_string(), Some(r.gated_iterations + 1)),
                };
                if r.gated_iterations < r.iterations.1 {
                    departed = true;
                }
                xfoil_transition = Some(r.xfoil_transition);
                xfoil_stations = Some(r.xfoil_stations);
                (
                    Some((
                        r.xfoil[0],
                        r.xfoil[1],
                        r.xfoil[2],
                        r.xfoil[3],
                        r.iterations.1,
                        r.converged.1,
                    )),
                    Some(o),
                    at,
                    Some(r.gated_iterations),
                )
            } else {
                (None, None, None, None)
            };
            let yfoil_nonfinite = p.nonfinite_station_failures > 0
                || ![p.cl, p.cd, p.cm, p.residual].iter().all(|v| v.is_finite())
                || p.iteration_records
                    .iter()
                    .any(|r| ![r.residual, r.cl, r.cd, r.cm].iter().all(|v| v.is_finite()));
            let xfoil_nonfinite = xfoil.map(|x| {
                k.and_then(|k| xfoil_nan.get(k - 1)).is_some_and(|n| *n > 0)
                    || ![x.0, x.1, x.2, x.3].iter().all(|v| v.is_finite())
                    || rec.iters.get(&k.unwrap_or(0)).is_some_and(|rows| {
                        rows.iter()
                            .any(|r| ![r[1], r[3], r[4], r[5]].iter().all(|v| v.is_finite()))
                    })
            });
            let k_for_nan = k;
            points.push(PolarPoint {
                leg: leg_no,
                alpha_deg: *a,
                yfoil: (p.cl, p.cd, p.cm, p.residual, p.iterations, p.converged),
                xfoil,
                outcome,
                xfoil_ill_conditioned: k.and_then(|k| reference_ill_conditioned(&rec, k)),
                parted_iteration: parted,
                gated_iterations: gated,
                yfoil_nonfinite,
                xfoil_nonfinite,
                yfoil_nonfinite_stations: p.nonfinite_station_failures,
                xfoil_nan_tokens: xfoil.and_then(|_| k_for_nan.map(|k| xfoil_nan.get(k - 1).copied().unwrap_or(0))),
                yfoil_transition: [p.transition_upper[0], p.transition_lower[0]],
                xfoil_transition,
                yfoil_stations: [
                    session.state().i_stagnation_node,
                    p.i_transition_station[1],
                    p.i_transition_station[2],
                ],
                xfoil_stations,
                state,
                transition_position: [
                    None,
                    TransitionPosition::from_state(session.state(), 1),
                    TransitionPosition::from_state(session.state(), 2),
                ],
            });
            // ASEQ's halting rule on yFoil's own convergence (the ALFA 0 that seeds the leg is
            // outside the sequence, as in XFOIL); OPER CL commands do not halt
            if j > 0 && matches!(op, Op::Alpha(_)) {
                if p.converged {
                    failures = 0;
                } else {
                    failures += 1;
                    if failures >= NSEQEX {
                        break;
                    }
                }
            }
        }
    };
    if case.alphas.is_empty() {
        let ops: Vec<Op> = case.cls.iter().map(|c| Op::Cl(*c)).collect();
        leg(&mut session, &ops, 1, &mut points);
    } else {
        let ops: Vec<Op> = case.alphas.iter().map(|a| Op::Alpha(*a)).collect();
        leg(&mut session, &ops, 1, &mut points);
    }
    if !case.alphas_after_reinit.is_empty() {
        session.init();
        let ops: Vec<Op> = case.alphas_after_reinit.iter().map(|a| Op::Alpha(*a)).collect();
        leg(&mut session, &ops, 2, &mut points);
    }
    PolarRun {
        name: case.name.clone(),
        section: case.section(),
        script: case.script(),
        description: description.to_string(),
        points,
        xfoil_calls: rec.points.len(),
        truncated_by_watchdog: truncated,
    }
}
