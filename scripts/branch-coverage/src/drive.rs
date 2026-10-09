//! Run one case through yFoil's `Session` exactly as its OPER script drove XFOIL, and compare
//! every call with the fixture (`scripts/study-support/records.rs`).

use crate::utilities::records::{
    compare_call, compare_inviscid_call, inviscid_ill_conditioned_calls, load, reference_ill_conditioned,
    transient_tol, CallReport, Outcome,
};
use crate::Case;
use std::path::Path;
use yfoil::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};

/// What the gate of one case concluded.
pub struct CaseGate {
    pub name: String,
    pub calls: Vec<CallReport>,
    /// the fixture's host (its manifest) and whether this is it
    pub fixture_host: String,
    pub same_host: bool,
    /// the reference's own 1-ULP spread: worst per-iteration and per-point values over the case
    pub floor_iteration_max: f64,
    pub floor_point_max: f64,
    /// SETBL/UPDATE/… subroutine calls the reference made (from the gcov measurement) — not
    /// known here; kept for the summary by the caller
    pub reference_iterations: usize,
    /// the calls at which the reference itself was ill-conditioned, from its own twins
    /// (`reference_ill_conditioned`); `None` when the case has no twins
    pub ill_conditioned_calls: Option<Vec<usize>>,
}

impl CaseGate {
    /// match | divergent | mismatch, over every call
    pub fn verdict(&self) -> &'static str {
        if self.calls.iter().any(|c| matches!(c.outcome, Outcome::Mismatch { .. })) {
            "mismatch"
        } else if self
            .calls
            .iter()
            .any(|c| matches!(c.outcome, Outcome::Divergent { .. }))
        {
            "divergent"
        } else {
            "match"
        }
    }
    /// the first call and iteration at which the runs were allowed to part, if any
    pub fn parted_at(&self) -> Option<(usize, usize)> {
        self.calls.iter().find_map(|c| match &c.outcome {
            Outcome::Divergent { at_iteration, .. } => Some((c.call, *at_iteration)),
            Outcome::Mismatch { .. } => Some((c.call, c.gated_iterations + 1)),
            Outcome::Match => None,
        })
    }
    pub fn max_transient(&self, name: &str) -> f64 {
        self.calls
            .iter()
            .flat_map(|c| c.transient_diff.iter().filter(|(n, _)| *n == name).map(|(_, d)| *d))
            .fold(0.0, f64::max)
    }
    pub fn max_point(&self, name: &str) -> Option<f64> {
        self.calls
            .iter()
            .flat_map(|c| c.point_diff.iter().filter(|(n, _)| *n == name).map(|(_, d)| *d))
            .reduce(f64::max)
    }
    pub fn worst_floor_ratio(&self) -> f64 {
        self.calls.iter().map(|c| c.worst_floor_ratio).fold(0.0, f64::max)
    }
}

pub fn conditions(case: &Case) -> FlowConditions {
    FlowConditions {
        re: if case.inviscid { None } else { Some(case.re) },
        mach: case.mach,
        ncrit: case.ncrit,
        max_iterations: case.max_iterations,
        x_trip: if case.xtr.is_empty() {
            [1.0, 1.0]
        } else {
            [case.xtr[0], case.xtr[1]]
        },
        mach_cl_dependence: match case.matyp {
            2 => MachClDependence::InverseSqrtCl,
            _ => MachClDependence::Fixed,
        },
        re_cl_dependence: match case.matyp {
            2 => ReClDependence::InverseSqrtCl,
            3 => ReClDependence::InverseCl,
            _ => ReClDependence::Fixed,
        },
        amplification_model: if case.damp {
            AmplificationModel::ModifiedEnvelope
        } else {
            AmplificationModel::Envelope
        },
        ..FlowConditions::default()
    }
}

pub fn gate(root: &Path, case: &Case) -> CaseGate {
    let dir = root.join("target/fixtures").join(&case.name);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    let rec = load(&dir);
    let fixture_host = crate::utilities::host::Host::of_fixture(&dir);
    let same_host = fixture_host == crate::utilities::host::Host::running();
    let mut session = Session::new(&airfoil, conditions(case));
    let mut calls: Vec<CallReport> = vec![];
    if case.inviscid {
        for (i, a) in case.alphas.iter().enumerate() {
            let p = session.alpha(a.to_radians());
            calls.push(compare_inviscid_call(&rec, "specal", i + 1, &p, session.state()));
        }
        for (i, c) in case.cls.iter().enumerate() {
            let p = session.cl(*c);
            calls.push(compare_inviscid_call(&rec, "speccl", i + 1, &p, session.state()));
        }
    } else {
        let mut k = 0;
        let mut leg = |session: &mut Session, alphas: &[f64], calls: &mut Vec<CallReport>| {
            for (j, a) in alphas.iter().enumerate() {
                if k >= rec.points.len() {
                    // the reference's sequence halted (or its run was ended by the watchdog)
                    // before the script's end: nothing further to compare
                    break;
                }
                k += 1;
                let p = if case.polar && j > 0 {
                    session.sequence_point(a.to_radians())
                } else {
                    session.alpha(a.to_radians())
                };
                calls.push(compare_call(&rec, k, &p, session.state(), transient_tol(k)));
            }
        };
        leg(&mut session, &case.alphas, &mut calls);
        if !case.alphas_after_reinit.is_empty() {
            session.init();
            leg(&mut session, &case.alphas_after_reinit, &mut calls);
        }
        for c in &case.cls {
            if k >= rec.points.len() {
                break;
            }
            k += 1;
            let p = session.cl(*c);
            calls.push(compare_call(&rec, k, &p, session.state(), transient_tol(k)));
        }
    }
    let (mut fi, mut fp) = (0.0_f64, 0.0_f64);
    if let Some(f) = &rec.floor {
        for c in &f.calls {
            for v in c.point.values() {
                fp = fp.max(*v);
            }
            for it in &c.iterations {
                for (key, v) in it {
                    if key != "LIMITER_FLIP" {
                        fi = fi.max(*v);
                    }
                }
            }
        }
        for c in f.specal_calls.iter().chain(&f.speccl_calls) {
            for v in c.values() {
                fp = fp.max(*v);
            }
        }
    }
    let reference_iterations = calls.iter().map(|c| c.iterations.1).sum();
    CaseGate {
        name: case.name.clone(),
        calls,
        fixture_host: fixture_host.to_string(),
        same_host,
        floor_iteration_max: fi,
        floor_point_max: fp,
        reference_iterations,
        ill_conditioned_calls: if rec.points.is_empty() {
            inviscid_ill_conditioned_calls(&rec)
        } else {
            rec.floor.as_ref().map(|_| {
                (1..=rec.points.len())
                    .filter(|&k| reference_ill_conditioned(&rec, k) == Some(true))
                    .collect()
            })
        },
    }
}
