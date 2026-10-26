//! The run's `summary.json`, its `README.md` / `README.tex` index, and the generated
//! documentation page (`docs/validation/branch-coverage/index.md`).

use crate::drive::CaseComparison;
use crate::utilities::publish::Bundle;
use crate::utilities::records::Outcome;
use crate::{Branch, Candidate, Universe};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

/// The trailing `# comment` of each case's `name = "…"` line in cases.toml, its description.
fn descriptions(root: &Path) -> BTreeMap<String, String> {
    let text = fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).unwrap_or_default();
    let mut out = BTreeMap::new();
    for l in text.lines() {
        let Some(rest) = l.trim_start().strip_prefix("name = \"") else {
            continue;
        };
        let Some((name, after)) = rest.split_once('"') else {
            continue;
        };
        if let Some((_, c)) = after.split_once('#') {
            out.insert(name.to_string(), c.trim().to_string());
        }
    }
    out
}

fn fmt_e(x: f64) -> String {
    if x == 0.0 {
        "0".into()
    } else if x.is_nan() {
        "NaN".into()
    } else {
        format!("{x:.1e}")
    }
}

/// Everything the tables and figures show, as one JSON document.
pub fn summary(
    universe: &Universe,
    candidates: &[Candidate],
    cover: &[String],
    compared: &[CaseComparison],
    nonfinite_comparisons: &[CaseComparison],
    group_matches: bool,
) -> Value {
    let root = crate::repo_root();
    let desc = descriptions(&root);
    let in_cover = |name: &str| cover.iter().any(|c| c == name);
    let finite: Vec<&Candidate> = candidates
        .iter()
        .filter(|c| !c.non_finite && !c.polar_study())
        .collect();
    let analysis: Vec<&Branch> = universe.branches.iter().filter(|b| b.is_analysis()).collect();
    // who takes each analysis branch, among finite candidates and among all
    let mut takers_finite: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut takers_all: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for c in candidates
        .iter()
        .filter(|c| c.case.group.as_deref() != Some("pathological"))
    {
        for id in &c.taken {
            takers_all.entry(id.as_str()).or_default().push(c.case.name.as_str());
            if !c.non_finite {
                takers_finite.entry(id.as_str()).or_default().push(c.case.name.as_str());
            }
        }
    }
    let cover_taken: BTreeSet<&str> = candidates
        .iter()
        .filter(|c| in_cover(&c.case.name))
        .flat_map(|c| c.taken.iter().map(String::as_str))
        .collect();

    // the branches taken by the measured `pathological` cases (only when their per-case
    // measurement exists: their reference runs hang or run for hours, so they are not measured by
    // default)
    let takers_pathological: BTreeSet<&str> = candidates
        .iter()
        .filter(|c| c.case.group.as_deref() == Some("pathological"))
        .flat_map(|c| c.taken.iter().map(String::as_str))
        .collect();
    // per-branch status over this study: cover (taken by the finite set) | nonfinite (taken by the
    // non-finite set, not the finite one) | pathological (taken only by a pathological case) |
    // unreachable | loop_entry | open (taken by nothing). The finite set is a cover of every branch
    // any finite candidate takes, so a branch taken by a finite candidate and not by the set is a
    // defect of the cover, not a status.
    let study_status = |b: &Branch| -> &'static str {
        if cover_taken.contains(b.id.as_str()) {
            "cover"
        } else if takers_finite.contains_key(b.id.as_str()) {
            panic!("{}: taken by a finite candidate and not by the finite set", b.id)
        } else if takers_all.contains_key(b.id.as_str()) {
            "nonfinite"
        } else if takers_pathological.contains(b.id.as_str()) {
            "pathological"
        } else {
            match b.status.as_str() {
                "open" => "open",
                "loop_entry" => "loop_entry",
                "unreachable" => "unreachable",
                _ => panic!(
                    "{}: taken in the coverage measurement by a case that is not a candidate here (no per-case measurement)",
                    b.id
                ),
            }
        }
    };
    let mut by_status: BTreeMap<&str, usize> = BTreeMap::new();
    for b in &analysis {
        *by_status.entry(study_status(b)).or_default() += 1;
    }
    let mut universe_status: BTreeMap<&str, usize> = BTreeMap::new();
    for b in &analysis {
        *universe_status.entry(b.status.as_str()).or_default() += 1;
    }
    let mut unreachable_by_class: BTreeMap<String, usize> = BTreeMap::new();
    let mut nonfinite_event_tested = 0usize;
    for b in &analysis {
        match study_status(b) {
            "unreachable" => *unreachable_by_class.entry(b.class.clone()).or_default() += 1,
            "nonfinite" if b.test.is_some() => nonfinite_event_tested += 1,
            _ => {}
        }
    }

    // subroutine rows (analysis path), in the coverage.toml order
    let mut sub_order: Vec<String> = vec![];
    for b in &analysis {
        let k = b.key();
        if !sub_order.contains(&k) {
            sub_order.push(k);
        }
    }
    let mut branch_map = vec![];
    let mut subroutines = vec![];
    for key in &sub_order {
        let bs: Vec<&&Branch> = analysis.iter().filter(|b| b.key() == *key).collect();
        let statuses: Vec<&str> = bs.iter().map(|b| study_status(b)).collect();
        let count = |s: &str| statuses.iter().filter(|x| **x == s).count();
        let (file, name) = key.split_once(':').unwrap();
        let calls = universe
            .subroutines
            .iter()
            .find(|r| r.file == file && r.name == name)
            .map(|r| r.calls)
            .unwrap_or(0);
        let mut covered_by: BTreeMap<&str, usize> = BTreeMap::new();
        for c in candidates.iter().filter(|c| in_cover(&c.case.name)) {
            let n = bs.iter().filter(|b| c.taken.contains(&b.id)).count();
            if n > 0 {
                covered_by.insert(c.case.name.as_str(), n);
            }
        }
        subroutines.push(json!({
            "key": key, "file": file, "name": name, "calls": calls, "branches": bs.len(),
            "cover": count("cover"), "nonfinite": count("nonfinite"), "pathological": count("pathological"),
            "open": count("open"), "loop_entry": count("loop_entry"), "unreachable": count("unreachable"),
            "covered_by": covered_by,
        }));
        branch_map.push(json!({
            "key": key, "file": file, "name": name,
            "statuses": statuses,
            "ids": bs.iter().map(|b| b.id.as_str()).collect::<Vec<_>>(),
        }));
    }

    // candidates and the cover's cases
    let candidates_json: Vec<Value> = candidates
        .iter()
        .map(|c| {
            let unique = c
                .taken
                .iter()
                .filter(|id| {
                    takers_finite
                        .get(id.as_str())
                        .is_some_and(|t| t.len() == 1 && t[0] == c.case.name)
                })
                .count();
            json!({
                "name": c.case.name, "section": c.case.section(), "script": c.case.script(),
                "description": desc.get(&c.case.name).cloned().unwrap_or_default(),
                "finite": !c.non_finite, "nan_count": c.nan_count,
                "branches_taken": c.taken.len(), "unique_among_finite": unique,
                "in_cover": in_cover(&c.case.name),
            })
        })
        .collect();
    let comparison_json = |g: &CaseComparison, c: &Candidate| -> Value {
        let calls: Vec<Value> = g
            .calls
            .iter()
            .map(|r| {
                json!({
                    "call": r.call, "alpha_deg": r.alpha_deg,
                    "iterations_yfoil": r.iterations.0, "iterations_xfoil": r.iterations.1,
                    "converged_yfoil": r.converged.0, "converged_xfoil": r.converged.1,
                    "compared_iterations": r.compared_iterations, "worst_floor_ratio": r.worst_floor_ratio,
                    "outcome": match &r.outcome {
                        Outcome::Match => json!("match"),
                        Outcome::Divergent { at_iteration, why } => json!({"divergent": {"at_iteration": at_iteration, "why": why}}),
                        Outcome::Mismatch { why } => json!({"mismatch": why}),
                    },
                })
            })
            .collect();
        let point = |n: &str| g.max_point(n);
        json!({
            "name": g.name, "section": c.case.section(), "script": c.case.script(),
            "description": desc.get(&g.name).cloned().unwrap_or_default(),
            "inviscid": c.case.inviscid, "nan_count": c.nan_count,
            "calls": calls, "reference_iterations": g.reference_iterations,
            "verdict": g.verdict(), "parted_at": g.parted_at(), "ill_conditioned_calls": g.ill_conditioned_calls,
            "compared_iterations": g.calls.iter().map(|r| r.compared_iterations).sum::<usize>(),
            "metrics": {
                "cl": point("CL"), "cd": point("CD"), "cm": point("CM"),
                "rmsbl": g.max_transient("RMSBL"), "rlx": g.max_transient("RLX"),
                "worst_floor_ratio": g.worst_floor_ratio(),
                "floor_iteration_max": g.floor_iteration_max, "floor_point_max": g.floor_point_max,
            },
        })
    };
    let nonfinite_runs: Vec<Value> = nonfinite_comparisons
        .iter()
        .map(|g| comparison_json(g, candidates.iter().find(|c| c.case.name == g.name).unwrap()))
        .collect();
    let mut cases = vec![];
    for name in cover {
        let c = candidates.iter().find(|c| &c.case.name == name).unwrap();
        let g = compared.iter().find(|g| &g.name == name).unwrap();
        let unique: Vec<&str> = c
            .taken
            .iter()
            .filter(|id| takers_finite.get(id.as_str()).is_some_and(|t| t.len() == 1))
            .map(String::as_str)
            .collect();
        let mut unique_subs: BTreeMap<String, usize> = BTreeMap::new();
        for id in &unique {
            let b = analysis.iter().find(|b| b.id == *id).unwrap();
            *unique_subs.entry(b.subroutine.clone()).or_default() += 1;
        }
        let calls: Vec<Value> = g
            .calls
            .iter()
            .map(|r| {
                json!({
                    "call": r.call, "alpha_deg": r.alpha_deg,
                    "iterations_yfoil": r.iterations.0, "iterations_xfoil": r.iterations.1,
                    "converged_yfoil": r.converged.0, "converged_xfoil": r.converged.1,
                    "compared_iterations": r.compared_iterations, "worst_floor_ratio": r.worst_floor_ratio,
                    "outcome": match &r.outcome {
                        Outcome::Match => json!("match"),
                        Outcome::Divergent { at_iteration, why } => json!({"divergent": {"at_iteration": at_iteration, "why": why}}),
                        Outcome::Mismatch { why } => json!({"mismatch": why}),
                    },
                    "transient_diff": r.transient_diff.iter().map(|(k, v)| (k.to_string(), v)).collect::<BTreeMap<_, _>>(),
                    "point_diff": r.point_diff.iter().map(|(k, v)| (k.to_string(), v)).collect::<BTreeMap<_, _>>(),
                })
            })
            .collect();
        let point = |n: &str| g.max_point(n);
        cases.push(json!({
            "name": name, "section": c.case.section(), "script": c.case.script(),
            "description": desc.get(name).cloned().unwrap_or_default(),
            "inviscid": c.case.inviscid, "sharp_te": c.case.foil.ends_with(":sharp") || c.case.foil == "karman-trefftz",
            "n_nodes": c.case.n_nodes,
            "branches_taken": c.taken.len(), "unique_branches": unique, "unique_by_subroutine": unique_subs,
            "calls": calls, "reference_iterations": g.reference_iterations,
            "verdict": g.verdict(), "parted_at": g.parted_at(), "ill_conditioned_calls": g.ill_conditioned_calls,
            "metrics": {
                "cl": point("CL"), "cd": point("CD"), "cm": point("CM"),
                "gam": point("GAM"), "qinv": point("QINV"), "cpi": point("CPI"),
                "rmsbl": g.max_transient("RMSBL"), "rlx": g.max_transient("RLX"),
                "worst_floor_ratio": g.worst_floor_ratio(),
                "floor_iteration_max": g.floor_iteration_max, "floor_point_max": g.floor_point_max,
            },
        }));
    }
    // every DO loop of the analysis path, from the per-case counts: did its body run on the set
    // (first body line executed), and what became of the arm that leaves the loop by its count
    let sum = |names: &dyn Fn(&Candidate) -> bool, key: &str| -> (Vec<u64>, u64) {
        let mut edges: Vec<u64> = vec![];
        let mut body = 0u64;
        for c in candidates.iter().filter(|c| names(c)) {
            if let Some((e, bcount)) = c.do_loops.get(key) {
                if edges.is_empty() {
                    edges = vec![0; e.len()];
                }
                for (x, y) in edges.iter_mut().zip(e) {
                    *x += y;
                }
                body += bcount;
            }
        }
        (edges, body)
    };
    let mut do_loops = vec![];
    let mut seen: BTreeSet<(String, usize)> = BTreeSet::new();
    for b in &analysis {
        if !b.source.trim_start().to_uppercase().starts_with("DO ") || !seen.insert((b.file.clone(), b.line)) {
            continue;
        }
        let key = format!("{}:{}", b.file, b.line);
        let is_cap = analysis
            .iter()
            .any(|x| x.file == b.file && x.line == b.line && x.why.is_some());
        let (set_edges, set_body) = sum(&|c| in_cover(&c.case.name), &key);
        let (fin_edges, fin_body) = sum(&|c| !c.non_finite, &key);
        let (all_edges, all_body) = sum(&|_| true, &key);
        // the iterate arm is the one with the larger count over the finite candidates; the other
        // is the arm that leaves the loop by its count
        let exit_idx = if fin_edges.len() == 2 {
            if fin_edges[0] >= fin_edges[1] {
                1
            } else {
                0
            }
        } else {
            usize::MAX
        };
        let exit_of = |e: &Vec<u64>| e.get(exit_idx).copied().unwrap_or(0);
        let exit = if exit_of(&set_edges) > 0 {
            format!("taken {}× on the set (the loop ran out its count)", exit_of(&set_edges))
        } else if exit_of(&fin_edges) > 0 {
            "taken outside the set".to_string()
        } else if exit_of(&all_edges) > 0 {
            "taken only in non-finite runs".to_string()
        } else if all_body == 0 {
            "never reached".to_string()
        } else {
            "never (the loop is always left early by a GO TO/RETURN, or its exit sits on its last statement)"
                .to_string()
        };
        do_loops.push(json!({
            "subroutine": b.key(), "line": b.line, "source": b.source,
            "iteration_cap": is_cap,
            "body_iterations_set": set_body, "body_iterations_finite": fin_body, "body_iterations_all": all_body,
            "edges_set": set_edges, "edges_finite": fin_edges, "edges_all": all_edges,
            "exit_arm": exit,
        }));
    }
    // the branches not taken by the finite set and not unreachable: the non-finite set's, the
    // pathological cases' and the ones nothing takes
    let mut open = vec![];
    for b in &analysis {
        let st = study_status(b);
        if st == "open" || st == "nonfinite" || st == "pathological" {
            open.push(json!({
                "id": b.id, "subroutine": b.key(), "line": b.line, "branch": b.branch, "source": b.source,
                "status": st, "class": b.class, "why": b.why, "test": b.test,
                "taken_by": takers_all.get(b.id.as_str()).cloned().unwrap_or_else(|| {
                    candidates
                        .iter()
                        .filter(|c| c.taken.contains(&b.id))
                        .map(|c| c.case.name.as_str())
                        .collect()
                }),
            }));
        }
    }
    // the case × subroutine matrix
    let matrix: Vec<Vec<f64>> = cover
        .iter()
        .map(|name| {
            let c = candidates.iter().find(|c| &c.case.name == name).unwrap();
            sub_order
                .iter()
                .map(|key| {
                    let bs: Vec<&&Branch> = analysis.iter().filter(|b| b.key() == *key).collect();
                    if bs.is_empty() {
                        0.0
                    } else {
                        bs.iter().filter(|b| c.taken.contains(&b.id)).count() as f64 / bs.len() as f64
                    }
                })
                .collect()
        })
        .collect();
    json!({
        "universe": {
            "cases_measured": universe.cases,
            "branches_translated_set": universe.branches.len(),
            "subroutines_translated_set": universe.subroutines.len(),
            "branches_analysis": analysis.len(),
            "universe_status": universe_status,
            "study_status": by_status,
            "unreachable_by_class": unreachable_by_class,
            "nonfinite_event_tested": nonfinite_event_tested,
            "pathological_measured": candidates.iter().filter(|c| c.case.group.as_deref() == Some("pathological")).count(),
            "finite_candidates": finite.len(),
            "candidates": candidates.len(),
        },
        "geometry_excluded": { "files": crate::GEOMETRY_FILES, "subroutines": crate::GEOMETRY_SUBROUTINES },
        "subroutines": subroutines,
        "branch_map": branch_map,
        "candidates": candidates_json,
        "cover": cover,
        "group_matches_cover": group_matches,
        "cases": cases,
        "nonfinite_runs": nonfinite_runs,
        "open": open,
        "do_loops": do_loops,
        "matrix": { "cases": cover, "subroutines": sub_order, "fraction": matrix },
    })
}

/// What each status in the tables and figures means. Everything here is measured on the
/// *reference* (gcov on pristine double-precision XFOIL); yFoil is never instrumented — its side
/// of the claim is the whole-run comparison, which asserts that on every run of the finite set
/// yFoil follows XFOIL's branch trace and values.
pub const CLASSIFICATION_MD: &str = "\
### What the statuses mean\n\n\
A *branch* is one outgoing edge of a conditional as gfortran emitted it for the reference build: an `IF` \
gives two edges, every operand of a `.AND.`/`.OR.` chain its own pair, a `DO` line the pair of the \
loop-control test gfortran placed there (iterate again, or leave the loop by its count), and a loop whose \
back-edge test lands on its last statement a pair there. Which arm is *fallthrough* and which is *jump* is \
the compiler's choice, so an edge is named by line and index, never by true/false; the source line and the \
counts of its sibling edges identify the arm. Every count is gcov's on the **reference**: XFOIL 6.99, \
pristine, double precision, compiled with `-fprofile-arcs -ftest-coverage`. yFoil is not instrumented and \
its own control flow is never measured; yFoil's side of the claim is the whole-run comparison below, which \
asserts that on each run of the finite set yFoil reproduces XFOIL's branch trace (iteration counts, IST/ITRAN, \
LVCONV, SPECCL's ITAL) and values. A branch taken by the finite set is therefore one the reference took on \
a run that yFoil matched.\n\n\
Every branch of the analysis path has exactly one of five statuses, drawn in these colours in the branch map:\n\n\
- **taken by the finite set** (green) — the reference took the edge on at least one run of the \
  `branch-coverage` group, and yFoil matches every run of that group. The set is chosen to take every branch \
  that any finite candidate takes, and the study stops with an error if it does not.\n\
- **taken by the non-finite set** (blue) — the edge is taken only on runs of the `non-finite` group, whose \
  reference output contains `NaN`. In every case found so far XFOIL produced the non-finite value in \
  TRCHEK2, whose transition-point Newton lands on the interval end; the MRCHUE/MRCHDU iterate fallback \
  then puts finite values on the station and the run carries on. These runs cannot be compared whole, \
  because the reference's own 1-ULP twins wander by O(10²) on them. Each edge is tested instead as an \
  **event**: the instrumented reference writes `events.dat` (site, SETBL call, side, station, deciding value) \
  and keeps the state dumps of those calls, and a test in `tests/execution/events.rs` seeds yFoil with \
  XFOIL's exact state entering the call, runs the one subroutine or the one iteration, asserts that yFoil's \
  branch trace records the same event, and compares the result with the reference's dump — NaN where the \
  reference has NaN. (TRCHEK2's `AX ≤ 0` exit, taken many times within one call, is recorded \
  call by call instead — `trchek2_ax_exit.dat`, the inputs of each TRCHEK2 call that takes it — and \
  `tests/subroutine/transition.rs` replays each one.) The non-finite set is part of the test suite (`group = \"non-finite\"`, tracked), so \
  these branches are covered by tests one call at a time; what is not claimed is that yFoil follows a whole \
  non-finite run. The table *Branches taken only by the non-finite set* names the test of each edge.\n\
- **taken only in pathological cases** (orange) — the edge is taken only by a case of the `pathological` \
  group: copies of the runs that hang in XFOIL's plot label or take hours (known issues §4), which cannot be \
  tests. Those cases are measured only on request (`cargo xtask coverage --per-case --with-group \
  pathological`); the status table says whether this run measured them.\n\
- **not reachable** (grey) — no execution of XFOIL's analysis path takes the edge. Either a recorded \
  code-reading argument says so (*annotated unreachable*, in one of the five classes below), or the edge is \
  the never-taken arm of a `DO` line's loop-control test.\n\
- **not taken** (red) — no measured run took the edge, and no argument says it cannot be taken. This is a \
  statement about the search, not about the code: no case has yet been constructed that reaches the edge. \
  Each such edge carries a *reach* note — the input the code reading says would reach it.\n\n\
An **annotated unreachable** edge has a recorded claim (`xtask/fixtures-config/coverage.toml`, with a class \
and the reason) that no input reaches it. The claim is a code-reading argument about the reference's own \
source, never \"we could not build a case\" and never a statement about yFoil. Every claim is falsifiable \
by the measurement: a case that takes an annotated edge makes the tool report the annotation as stale (this \
happened to `xoper.f:2784` branch 5, whose claim missed that MRCL's Mach clamp sets `M_CLS = 0`; the claim \
was withdrawn and the edge is now taken). The classes are:\n\n\
- *structural* — impossible by construction of the data the analysis path feeds the routine: coincident \
  consecutive nodes that `LOAD` removes, `IBL ≤ 2` inside an `IF(IBL.GT.3)` block, a quantity kept positive \
  by every update that touches it, a subroutine only ever called under the flag it tests;\n\
- *mode* — guarded by a flag the analysis path never sets: inverse design (`GEOLIN`), the image airfoil \
  (`LIMAGE`, whose only toggle is commented out), flap hinge moments, interactive prompts. These are \
  features outside yFoil's scope, and the edges are unreachable in XFOIL too when driven as the fixture \
  pipeline drives it;\n\
- *guard* — array-bound and illegal-input `STOP`s on Fortran's fixed dimensions (`IVX`, `IWX`, …); yFoil \
  has no fixed dimensions and the cases are within XFOIL's;\n\
- *compiler* — edges gfortran emits that correspond to no decision in the source;\n\
- *numerical* — reachable in principle, but only by a floating-point coincidence (a real landing bitwise \
  on a comparison's other operand: the stagnation point exactly on a node, SETEXP's initial guess exactly \
  1.0) or by exhausting the iteration cap of a Newton iteration that converges on every analysis-path input \
  (SETEXP's 100, SINVRT's 10, MRCHUE's 25 at the similarity station). No input can be designed to reach \
  them; if anything does, it is the parameter/fuzz harness.\n\n\
A **DO-line edge** comes from the control test a Fortran `DO` loop is compiled with: *iterate again* or \
*leave the loop by its count* (the count exhausted, or zero to begin with). gfortran attributes that test to \
the `DO` line as a pair of edges, and the measured counts show two shapes: for a top-tested loop the second \
arm is the loop's normal completion, taken once per execution; for a loop gfortran rotated, the real exit \
sits on the loop's last statement and the `DO` line's second arm is a duplicate that is never taken. Neither \
arm is a body that ran zero times as distinct from a body that ran and finished: the analysis path's \
station, node and wake counts never make a loop empty, and the counts cannot tell an empty pass from a \
completed one. The *iterate* arm is an ordinary edge, taken whenever a run executed the body, and it appears \
in the branch map beside its partner. The grey DO-line edges are the never-taken second arms: a loop always \
left early by a `GO TO` or `RETURN` (its normal completion never happens), a rotated loop's duplicate, or \
both arms of a loop inside a block that is never reached. Whether each loop's body actually ran on the set \
is stated explicitly, from the execution count of the first line of its body, in the table *DO loops of the \
analysis path* below. A `DO` line that is an iteration budget (a Newton loop that could be exhausted: \
SPECCL's twenty attempts, SPECAL's CL(M) loop) is different: its second arm is the exhausted budget, a real \
outcome, so it carries a reach note and both arms are counted like any other conditional.\n\n\
The same classes, with every annotated edge and its reason, are listed in `docs/development/testing/coverage.md`; \
that page is the union over every tracked case, this one the minimal set.\n";

fn verdict_word(v: &str) -> &str {
    match v {
        "match" => "match",
        "divergent" => "divergent",
        _ => "MISMATCH",
    }
}

fn parted(case: &Value) -> String {
    match case["parted_at"].as_array() {
        Some(a) if a.len() == 2 => format!("call {} it. {}", a[0], a[1]),
        _ => "—".into(),
    }
}

fn metric(case: &Value, key: &str) -> String {
    match case["metrics"][key].as_f64() {
        Some(x) => fmt_e(x),
        None => "—".into(),
    }
}

/// The case table (definition and coverage) as Markdown rows.
fn case_table_md(s: &Value) -> String {
    let mut t = String::from(
        "| case | section | OPER | branches | unique | unique branches in |\n|---|---|---|---:|---:|---|\n",
    );
    for c in s["cases"].as_array().unwrap() {
        let subs = c["unique_by_subroutine"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| format!("{k} ({v})"))
            .collect::<Vec<_>>()
            .join(", ");
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["section"].as_str().unwrap(),
            c["script"].as_str().unwrap(),
            c["branches_taken"],
            c["unique_branches"].as_array().unwrap().len(),
            if subs.is_empty() { "—".to_string() } else { subs }
        );
    }
    t
}

/// The whole-run comparison table (difference metrics) as Markdown rows.
fn comparison_table_md(s: &Value) -> String {
    let mut t = String::from(
        "| case | calls | ref. iterations | outcome | parts at | XFOIL ill-conditioned at (its twins) | max \\|ΔCL\\| | max \\|ΔCD\\| | max \\|ΔCM\\| | max \\|ΔRMSBL\\| / \\|ΔGAM\\| | worst diff/floor | ref. 1-ULP spread |\n|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|\n",
    );
    for c in s["cases"].as_array().unwrap() {
        let inviscid = c["inviscid"].as_bool().unwrap_or(false);
        let trans = if inviscid { metric(c, "gam") } else { metric(c, "rmsbl") };
        let floor = if inviscid {
            metric(c, "floor_point_max")
        } else {
            metric(c, "floor_iteration_max")
        };
        let ill = match c["ill_conditioned_calls"].as_array() {
            None => "—".to_string(),
            Some(a) if a.is_empty() => "no call".to_string(),
            Some(a) => format!(
                "call{} {}",
                if a.len() > 1 { "s" } else { "" },
                a.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(", ")
            ),
        };
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {ill} | {} | {} | {} | {}{} | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["calls"].as_array().unwrap().len(),
            c["reference_iterations"],
            verdict_word(c["verdict"].as_str().unwrap()),
            parted(c),
            metric(c, "cl"),
            metric(c, "cd"),
            metric(c, "cm"),
            trans,
            if inviscid { " †" } else { "" },
            metric(c, "worst_floor_ratio"),
            floor
        );
    }
    t
}

/// The whole-run outcome of every non-finite probe.
fn nonfinite_table_md(s: &Value) -> String {
    let mut t = String::from(
        "| case | OPER | `NaN` in reference output | ref. iterations | compared | outcome | parts at | max \\|ΔRMSBL\\| (compared) | worst diff/floor | ref. 1-ULP spread |\n|---|---|---:|---:|---:|---|---|---:|---:|---:|\n",
    );
    for c in s["nonfinite_runs"].as_array().unwrap() {
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["script"].as_str().unwrap(),
            c["nan_count"],
            c["reference_iterations"],
            c["compared_iterations"],
            verdict_word(c["verdict"].as_str().unwrap()),
            parted(c),
            metric(c, "rmsbl"),
            metric(c, "worst_floor_ratio"),
            metric(c, "floor_iteration_max")
        );
    }
    t
}

fn status_table_md(s: &Value) -> String {
    let u = &s["universe"];
    let st = |k: &str| u["study_status"][k].as_u64().unwrap_or(0);
    let cls = |k: &str| u["unreachable_by_class"][k].as_u64().unwrap_or(0);
    let tested = u["nonfinite_event_tested"].as_u64().unwrap_or(0);
    let measured = u["pathological_measured"].as_u64().unwrap_or(0);
    let pathological = if measured == 0 {
        "not measured in this run".to_string()
    } else {
        format!("{} (over {measured} measured cases)", st("pathological"))
    };
    format!(
        "| | branches |\n|---|---:|\n\
         | translated set (all files) | {} |\n\
         | analysis path (panelling generators excluded) | {} |\n\
         | **taken by the finite set** (green) | **{}** |\n\
         | **taken by the non-finite set** (blue) | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;with a test that replays it | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;without one | {} |\n\
         | **taken only in pathological cases** (orange) | {} |\n\
         | **not reachable** (grey) | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;annotated *structural* / *mode* / *guard* / *compiler* / *numerical* | {} / {} / {} / {} / {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;never-taken DO-line edges (loop left early or never reached) | {} |\n\
         | **not taken** (red) | {} |\n",
        u["branches_translated_set"],
        u["branches_analysis"],
        st("cover"),
        st("nonfinite"),
        tested,
        st("nonfinite") - tested,
        pathological,
        st("unreachable") + st("loop_entry"),
        cls("structural"),
        cls("mode"),
        cls("guard"),
        cls("compiler"),
        cls("numerical"),
        st("loop_entry"),
        st("open")
    )
}

/// The rows of `s["open"]` with the given status, as a Markdown table: the source, which cases
/// take the edge, the reach note and the test of its event.
fn branch_list_md(s: &Value, status: &str) -> String {
    let rows: Vec<&Value> = s["open"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|o| o["status"].as_str() == Some(status))
        .collect();
    if rows.is_empty() {
        return String::new();
    }
    let mut t = String::from("| branch | source | taken by | how it is reached | test |\n|---|---|---|---|---|\n");
    for o in rows {
        let taken_by = o["taken_by"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| format!("`{}`", v.as_str().unwrap()))
            .collect::<Vec<_>>()
            .join(", ");
        let test = o["test"]
            .as_str()
            .map(|g| {
                let (file, test) = g.split_once("::").unwrap_or((g, ""));
                format!("`{test}` in `{file}`")
            })
            .unwrap_or_else(|| "none yet".into());
        t += &format!(
            "| `{}` | `{}` | {} | {} | {} |\n",
            o["id"].as_str().unwrap(),
            o["source"].as_str().unwrap().replace('|', "\\|"),
            if taken_by.is_empty() {
                "—".to_string()
            } else {
                taken_by
            },
            o["why"].as_str().unwrap_or(""),
            test
        );
    }
    t
}

/// Every DO loop of the analysis path: whether its body ran on the set (execution count of the
/// first line of the body) and what became of the arm that leaves the loop by its count.
fn do_loop_table_md(s: &Value) -> String {
    let loops = s["do_loops"].as_array().unwrap();
    let ran = loops
        .iter()
        .filter(|l| l["body_iterations_set"].as_u64().unwrap() > 0)
        .count();
    let completed = loops
        .iter()
        .filter(|l| l["exit_arm"].as_str().unwrap().contains("on the set (the loop ran out"))
        .count();
    let never = loops
        .iter()
        .filter(|l| l["exit_arm"].as_str().unwrap().starts_with("never ("))
        .count();
    let unreached = loops
        .iter()
        .filter(|l| l["exit_arm"].as_str().unwrap() == "never reached")
        .count();
    let mut t = format!(
        "One row per `DO` statement of the analysis path. *Body ran* is the number of times the first \
         executable line of the loop body ran over the runs of the set — the explicit statement that the loop \
         was entered, independent of how gfortran attributed its control edges. *Exit arm* is the arm of the \
         `DO` line's control test that leaves the loop by its count — for an iteration budget such as SPECCL's \
         `DO 100 ITAL=1, 20` that is the exhausted budget, for a search loop such as STFIND's scan it is \
         'not found'. Of the {} loops, the set ran the body of {}; the exit arm was taken on the set for {} \
         (the loop ran out its count at least once), {} were never left by their count on any run (always \
         left early by a `GO TO`/`RETURN`, or the exit sits on the loop's last statement), and {} are never \
         reached at all (inside a block the analysis path does not enter, or a subroutine XFOIL never \
         calls).\n\n\
         | subroutine | line | source | body ran (set) | exit arm |\n|---|---:|---|---:|---|\n",
        loops.len(),
        ran,
        completed,
        never,
        unreached
    );
    for l in loops {
        let body = l["body_iterations_set"].as_u64().unwrap();
        t += &format!(
            "| {} | {} | `{}` | {} | {} |\n",
            l["subroutine"].as_str().unwrap(),
            l["line"],
            l["source"].as_str().unwrap().replace('|', "\\|"),
            if body > 0 {
                format!("{body}")
            } else {
                "**0**".to_string()
            },
            l["exit_arm"].as_str().unwrap()
        );
    }
    t + "\n"
}

fn excluded_md(s: &Value) -> String {
    let rows: Vec<&Value> = s["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| !c["finite"].as_bool().unwrap())
        .collect();
    if rows.is_empty() {
        return "None of the tracked candidates printed a non-finite value.\n".into();
    }
    let mut t = String::from("| case | OPER | `NaN` occurrences in the reference output |\n|---|---|---:|\n");
    for c in rows {
        t += &format!(
            "| `{}` | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["script"].as_str().unwrap(),
            c["nan_count"]
        );
    }
    t
}

fn body_md(s: &Value, figures: bool, figure_prefix: &str) -> String {
    let u = &s["universe"];
    let mut m = String::new();
    m += "## Which branches are counted\n\n";
    m += "Before a minimal set can be chosen, the branches it has to take must be fixed: this section \
          defines them, and every count on the page is out of this total.\n\n";
    m += &format!(
        "The translated set is the {} subroutines of `xtask/fixtures-config/coverage.toml`, measured by \
         gcov on the pristine double-precision reference (`docs/development/testing/coverage.md`). ",
        s["universe"]["subroutines_translated_set"]
    );
    m += "The *analysis path* is the translated set without the panelling generators and their helpers \
          (`naca.f`, `xgeom.f`, `spline.f`, and `NACA`/`PANGEN` in `xfoil.f`), which the `pangen_*` cases test \
          bitwise (`docs/development/testing/geometry/README.md`). A branch is one outgoing edge of a conditional as \
          gfortran emitted it; the counts below are edges.\n\n";
    m += &status_table_md(s);
    m += "\n";
    m += CLASSIFICATION_MD;
    m += &format!(
        "\n{} tracked solver cases were measured one at a time; {} have a finite reference run and are the \
         candidates for the finite set. A case whose reference prints `NaN` is not a candidate: a non-finite \
         trajectory cannot be compared as a whole run (NaN comparisons and `MIN`/`MAX` are where two \
         correct codes legitimately differ), so the branches only such a run reaches are tested as events \
         instead.\n\n",
        u["candidates"], u["finite_candidates"]
    );
    m += "## The minimal finite set\n\n";
    m += "Branches taken by exactly one finite candidate make that candidate essential; the remainder is \
          covered greedily by marginal gain and pruned until no case can be removed. The set is recorded as \
          `group = \"branch-coverage\"` in `xtask/fixtures-config/cases.toml`, and this study refuses to \
          generate if the two differ. *Branches* is the number of analysis-path edges the case takes; \
          *unique* the ones no other finite candidate takes.\n\n";
    m += &case_table_md(s);
    m += "\nEvery case: yFoil generates the panels (`--method cosine`), XFOIL `LOAD`s them (bitwise handoff \
          asserted), one OPER script, one instrumented run plus its five 1-ULP twins.\n\n";
    if figures {
        m += &format!(
            "[Branch map]({figure_prefix}branch-map.fig) ([PDF]({figure_prefix}branch-map.pdf)): every branch of \
             the analysis path, one row per subroutine in XFOIL's file order: taken by the finite set (green), \
             taken by the non-finite set (blue), taken only in pathological cases (orange), not reachable — \
             annotated unreachable or a never-taken DO-line edge (grey), not taken (red).\n\n\
             [Case × subroutine]({figure_prefix}case-matrix.fig) ([PDF]({figure_prefix}case-matrix.pdf)): \
             fraction of each subroutine's branches taken by each case of the set; a dot marks a subroutine in \
             which the case takes a branch no other finite candidate takes.\n\n"
        );
    }
    m += "## DO loops of the analysis path\n\n";
    m += &do_loop_table_md(s);
    m += "## Whole-run comparison of the finite set\n\n";
    m += "Each case of the set is run through yFoil's `Session` exactly as its OPER script drove XFOIL and \
          compared with its fixture by the studies' floor comparison (`studies/study-support/records.rs`; the \
          tests compare single steps instead): iteration count, LVCONV, IST and ITRAN exact; every \
          per-iteration transient and the converged point within `max(tol · scale, FLOOR_FACTOR · floor)`, \
          where *floor* is the reference's own 1-ULP spread of that value; **divergent** when every earlier \
          iteration matched and the reference itself moves by more than `DIVERGENCE_FLOOR` there. The \
          differences are absolute, maxima over the compared calls and iterations. *Worst diff/floor* is the \
          largest difference relative to the reference's own spread of the same value (a value below \
          `FLOOR_FACTOR` = 4 is within tolerance); *ref. 1-ULP spread* is the largest spread the twins showed \
          on the case, per-iteration (viscous) or per-point (inviscid) — when it is large the comparison is \
          only as sharp as the reference itself. † inviscid case: the per-node |ΔGAM| maximum (RMSBL does not \
          exist).\n\n";
    m += &comparison_table_md(s);
    m += "\n## Whole-run comparison of the non-finite set\n\n";
    m += "This is a measurement, not a claim — the reference cannot reproduce itself on these runs — of how \
          far yFoil follows each one. The comparison above is applied call by call: *parts at* is the first \
          iteration at which the runs were allowed to part (the reference's own 1-ULP spread exceeds \
          `DIVERGENCE_FLOOR` there and the values differ) or, for a mismatch, the iteration after the last one \
          compared; *compared* counts the iterations matched before that. The NaN count is the number of \
          `NaN` tokens in the reference's own output for the run.\n\n";
    m += &nonfinite_table_md(s);
    m += "\n## Full polars\n\nFull ±30° polars of the aerofoil series and regimes of interest are the \
          series-cases validation study, and the runs that show \
          XFOIL's known issues the known-issues study. Those sweeps are not \
          candidates for the finite set, which is a minimal set of single operating points.\n";
    m += "\n## Branches taken only by the non-finite set\n\n";
    let nonfinite = branch_list_md(s, "nonfinite");
    if nonfinite.is_empty() {
        m += "None.\n";
    } else {
        m += "The non-finite set takes these branches and the finite set does not. The last column names the \
              test that replays each one: as an event, seeding yFoil with XFOIL's state dumped entering the \
              SETBL call at which `events.dat` records the event and requiring the same event and XFOIL's \
              state after it, or, for TRCHEK2's `AX ≤ 0` exit, each TRCHEK2 call the reference recorded \
              taking it.\n\n";
        m += &nonfinite;
    }
    m += "\n## Branches taken only in pathological cases\n\n";
    let pathological = branch_list_md(s, "pathological");
    m += &if u["pathological_measured"].as_u64().unwrap_or(0) == 0 {
        "The `pathological` cases were not measured in this run (`cargo xtask coverage --per-case \
         --with-group pathological`; their reference runs hang or take hours).\n"
            .to_string()
    } else if pathological.is_empty() {
        "None: every branch the measured pathological cases take is taken by one of the sets above.\n".to_string()
    } else {
        pathological
    };
    m += "\n## Branches not taken\n\n";
    let open = branch_list_md(s, "open");
    if open.is_empty() {
        m += "None: every branch of the analysis path is taken by one of the sets above or is not reachable.\n";
    } else {
        m += "No measured run takes these branches and no argument says they cannot be taken; the reach note \
              is the input the code reading says would reach each one.\n\n";
        m += &open;
    }
    m += "\n## The non-finite set\n\n";
    m += "The tracked candidates whose reference run printed a non-finite value, and so are not candidates for \
          the finite set.\n\n";
    m += &excluded_md(s);
    m
}

pub fn write_index(run_dir: &Path, s: &Value) {
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let mut m = String::from("# Branch-coverage study\n\n");
    m += &format!(
        "Generated by `cargo run --release -p branch-coverage` on {} (git {}{}, host {}); provenance in \
         `metadata.json`, every number in `summary.json`, figures drawn from it by `figures.rs` (ironlab).\n\n",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        if meta["git_dirty"].as_bool().unwrap_or(false) {
            ", dirty"
        } else {
            ""
        },
        meta["host"].as_str().unwrap_or("")
    );
    m += &body_md(s, true, "");
    fs::write(run_dir.join("README.md"), m).unwrap();
}

fn tex_escape(s: &str) -> String {
    s.replace('\\', "\\textbackslash{}")
        .replace('&', "\\&")
        .replace('%', "\\%")
        .replace('#', "\\#")
        .replace('_', "\\_")
        .replace('|', "\\textbar{}")
        .replace('×', "$\\times$")
        .replace('Δ', "$\\Delta$")
        .replace('α', "$\\alpha$")
        .replace('—', "---")
        .replace('–', "--")
        .replace('†', "$\\dagger$")
}

fn tex_num(s: &str) -> String {
    // 1.2e-10 -> $1.2\times10^{-10}$
    if let Some((m, e)) = s.split_once('e') {
        format!("${m}\\times10^{{{}}}$", e.trim_start_matches('+'))
    } else {
        tex_escape(s)
    }
}

pub fn write_index_tex(run_dir: &Path, s: &Value) {
    let mut f = fs::File::create(run_dir.join("README.tex")).unwrap();
    let w = |f: &mut fs::File, t: &str| writeln!(f, "{t}").unwrap();
    w(
        &mut f,
        "% Generated by `cargo run --release -p branch-coverage`; tables from summary.json. Needs booktabs.",
    );
    w(&mut f, "\\section*{Branch-coverage study}");
    w(&mut f, "\\noindent The minimal set of cases taking every finite-reachable branch of XFOIL's translated analysis path, each compared with its instrumented double-precision reference. Provenance in \\texttt{metadata.json}, every number in \\texttt{summary.json}.");
    // case table
    w(&mut f, "\\begin{table}[h]\\centering\\small");
    w(&mut f, "\\caption{The minimal case set: section, OPER script, analysis-path branches taken and the branches no other finite candidate takes.}");
    w(&mut f, "\\begin{tabular}{@{}llp{5.2cm}rrp{4cm}@{}}\\toprule");
    w(
        &mut f,
        "Case & Section & OPER & Branches & Unique & Unique branches in \\\\ \\midrule",
    );
    for c in s["cases"].as_array().unwrap() {
        let subs = c["unique_by_subroutine"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| format!("{k} ({v})"))
            .collect::<Vec<_>>()
            .join(", ");
        w(
            &mut f,
            &format!(
                "\\texttt{{{}}} & {} & {} & {} & {} & {} \\\\",
                tex_escape(c["name"].as_str().unwrap()),
                tex_escape(c["section"].as_str().unwrap()),
                tex_escape(c["script"].as_str().unwrap()),
                c["branches_taken"],
                c["unique_branches"].as_array().unwrap().len(),
                if subs.is_empty() {
                    "---".to_string()
                } else {
                    tex_escape(&subs)
                }
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    // whole-run comparison table
    w(&mut f, "\\begin{table}[h]\\centering\\small");
    w(&mut f, "\\caption{Whole-run comparison: every case of the set against its fixture. Differences are absolute maxima over the compared calls and iterations; \\emph{worst diff/floor} is relative to the reference's own 1-ULP spread of the same value; $\\dagger$ inviscid case, per-node $|\\Delta\\mathrm{GAM}|$ in place of RMSBL.}");
    w(&mut f, "\\begin{tabular}{@{}lrrlrrrrrr@{}}\\toprule");
    w(&mut f, "Case & Calls & Iter. & Outcome & $|\\Delta C_L|$ & $|\\Delta C_D|$ & $|\\Delta C_M|$ & $|\\Delta\\mathrm{RMSBL}|$ & diff/floor & 1-ULP spread \\\\ \\midrule");
    for c in s["cases"].as_array().unwrap() {
        let inviscid = c["inviscid"].as_bool().unwrap_or(false);
        let trans = if inviscid { metric(c, "gam") } else { metric(c, "rmsbl") };
        let floor = if inviscid {
            metric(c, "floor_point_max")
        } else {
            metric(c, "floor_iteration_max")
        };
        let outcome = match c["verdict"].as_str().unwrap() {
            "match" => "match".to_string(),
            "divergent" => format!("divergent ({})", tex_escape(&parted(c))),
            _ => format!("\\textbf{{mismatch}} ({})", tex_escape(&parted(c))),
        };
        w(
            &mut f,
            &format!(
                "\\texttt{{{}}}{} & {} & {} & {} & {} & {} & {} & {} & {} & {} \\\\",
                tex_escape(c["name"].as_str().unwrap()),
                if inviscid { "$^\\dagger$" } else { "" },
                c["calls"].as_array().unwrap().len(),
                c["reference_iterations"],
                outcome,
                tex_num(&metric(c, "cl")),
                tex_num(&metric(c, "cd")),
                tex_num(&metric(c, "cm")),
                tex_num(&trans),
                tex_num(&metric(c, "worst_floor_ratio")),
                tex_num(&floor)
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    // status table
    let u = &s["universe"];
    let st = |k: &str| u["study_status"][k].as_u64().unwrap_or(0);
    w(&mut f, "\\begin{table}[h]\\centering\\small\\caption{Branches of the analysis path by status over the study.}\\begin{tabular}{@{}lr@{}}\\toprule & Branches \\\\ \\midrule");
    w(
        &mut f,
        &format!("Translated set (all files) & {} \\\\", u["branches_translated_set"]),
    );
    w(
        &mut f,
        &format!(
            "Analysis path (panelling generators excluded) & {} \\\\",
            u["branches_analysis"]
        ),
    );
    let cls = |k: &str| u["unreachable_by_class"][k].as_u64().unwrap_or(0);
    let tested = u["nonfinite_event_tested"].as_u64().unwrap_or(0);
    let measured = u["pathological_measured"].as_u64().unwrap_or(0);
    w(
        &mut f,
        &format!("Taken by the finite set (green) & \\textbf{{{}}} \\\\", st("cover")),
    );
    w(
        &mut f,
        &format!(
            "Taken by the non-finite set (blue), with a test that replays it / without & {} / {} \\\\",
            tested,
            st("nonfinite") - tested
        ),
    );
    w(
        &mut f,
        &if measured == 0 {
            "Taken only in pathological cases (orange) & not measured \\\\".to_string()
        } else {
            format!(
                "Taken only in pathological cases (orange) & {} \\\\",
                st("pathological")
            )
        },
    );
    w(
        &mut f,
        &format!(
            "Not reachable (grey): structural / mode / guard / compiler / numerical & {} / {} / {} / {} / {} \\\\",
            cls("structural"),
            cls("mode"),
            cls("guard"),
            cls("compiler"),
            cls("numerical")
        ),
    );
    w(
        &mut f,
        &format!(
            "Not reachable (grey): never-taken DO-line edges & {} \\\\",
            st("loop_entry")
        ),
    );
    w(
        &mut f,
        &format!(
            "Not taken (red) & {} \\\\ \\bottomrule\\end{{tabular}}\\end{{table}}",
            st("open")
        ),
    );
    w(&mut f, "\\noindent Every count is gcov's on the pristine double-precision reference; yFoil is never instrumented, and a branch taken by the finite set is one the reference took on a run of the set that yFoil matches as a whole run. \\emph{Taken by the non-finite set}: taken only on runs whose reference output contains NaN --- such a run cannot be compared whole, so each branch is tested as an event, replayed from XFOIL's state at that call, NaN where XFOIL has NaN. \\emph{Taken only in pathological cases}: taken only by the runs that hang or take hours, which cannot be tests. \\emph{Not reachable}: a recorded, falsifiable code-reading claim that no execution of XFOIL's analysis path can take the edge (structural, mode, guard, compiler or numerical: a floating-point coincidence or an exhausted Newton cap), or the never-taken arm of the loop-control test on a DO line. \\emph{Not taken}: no measured run took the edge and no argument says it cannot be taken; each carries a reach note. The Markdown page carries the full definitions.");
    w(&mut f, "\\begin{figure}[h]\\centering\\includegraphics{branch-map.pdf}\\caption{Every branch of the analysis path, one row per subroutine in XFOIL's file order: taken by the finite set (green), taken by the non-finite set (blue), taken only in pathological cases (orange), not reachable (grey), not taken (red).}\\end{figure}");
    w(&mut f, "\\begin{figure}[h]\\centering\\includegraphics{case-matrix.pdf}\\caption{Fraction of each subroutine's branches taken by each case of the set; a dot marks a subroutine in which the case takes a branch no other finite candidate takes; a subroutine the case never entered is light grey.}\\end{figure}");
}

/// The run's docs bundle (`publish/`): one fragment per table the study's pages include, and the
/// two figures.
pub fn write_publish(run_dir: &Path, s: &Value) {
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let u = &s["universe"];
    let bundle = Bundle::create(run_dir);
    bundle.provenance(
        "cargo run --release -p branch-coverage",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        meta["git_dirty"].as_bool().unwrap_or(false),
        meta["host"].as_str().unwrap_or(""),
    );
    bundle.table(
        "measured",
        &format!(
            "| | count |\n|---|---:|\n\
             | translated subroutines (`xtask/fixtures-config/coverage.toml`) | {} |\n\
             | tracked solver cases measured, one at a time | {} |\n\
             | of them with a finite reference run (candidates for the finite set) | {} |\n\
             | cases in the minimal finite set | {} |\n",
            u["subroutines_translated_set"],
            u["candidates"],
            u["finite_candidates"],
            s["cover"].as_array().map(Vec::len).unwrap_or(0),
        ),
    );
    bundle.table("status", &status_table_md(s));
    bundle.table("finite-set", &case_table_md(s));
    bundle.table("do-loops", &do_loop_table_md(s));
    bundle.table("finite-comparison", &comparison_table_md(s));
    bundle.table("non-finite-comparison", &nonfinite_table_md(s));
    let or_none = |t: String, none: &str| if t.is_empty() { format!("{none}\n") } else { t };
    bundle.table("branches-non-finite", &or_none(branch_list_md(s, "nonfinite"), "None."));
    bundle.table(
        "branches-pathological",
        &if u["pathological_measured"].as_u64().unwrap_or(0) == 0 {
            "The `pathological` cases were not measured in this run (`cargo xtask coverage --per-case \
             --with-group pathological`; their reference runs hang or take hours).\n"
                .to_string()
        } else {
            or_none(
                branch_list_md(s, "pathological"),
                "None: every branch the measured pathological cases take is taken by one of the sets.",
            )
        },
    );
    bundle.table(
        "branches-not-taken",
        &or_none(
            branch_list_md(s, "open"),
            "None: every branch of the analysis path is taken by one of the sets or is not reachable.",
        ),
    );
    bundle.table("non-finite-set", &excluded_md(s));
    bundle.table("subroutine-counts", &crate::figures::subroutine_counts_md(s));
    for name in ["branch-map.fig", "case-matrix.fig"] {
        bundle.figure(&run_dir.join(name));
    }
    println!(
        "publish/ written: cargo xtask publish branch-coverage --run {}",
        meta["generated"].as_str().unwrap_or("")
    );
}
