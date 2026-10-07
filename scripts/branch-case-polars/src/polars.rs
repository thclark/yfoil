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

use crate::utilities::records::{compare_call, compare_state, load, transient_tol, ArrayDiff, Outcome};
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
                "xfoil": p.xfoil.as_ref().map(pt), "outcome": p.outcome, "parted_iteration": p.parted_iteration,
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

pub fn sweep(root: &Path, case: &Case, description: &str) -> PolarRun {
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
    let leg = |session: &mut Session, alphas: &[f64], leg_no: usize, points: &mut Vec<PolarPoint>| {
        let mut failures = 0usize;
        // once a point of this leg was not followed to the end, every later point starts from a
        // state the reference did not compute: its comparison is informational, never a mismatch
        // (a leg starts afresh: INIT and a fresh march)
        let mut departed = false;
        for (j, a) in alphas.iter().enumerate() {
            let p = if j == 0 {
                session.alpha(a.to_radians())
            } else {
                session.sequence_point(a.to_radians())
            };
            let k = ref_calls.get(&(leg_no, (a * 1000.0).round() as i64)).copied();
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
            // outside the sequence, as in XFOIL)
            if j > 0 {
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
    leg(&mut session, &case.alphas, 1, &mut points);
    if !case.alphas_after_reinit.is_empty() {
        session.init();
        leg(&mut session, &case.alphas_after_reinit, 2, &mut points);
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
