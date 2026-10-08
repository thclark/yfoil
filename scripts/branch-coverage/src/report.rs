//! The run's `summary.json`, its `README.md` / `README.tex` index, and the generated
//! documentation page (`docs/validation/branch-coverage/README.md`).

use crate::drive::CaseGate;
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
    gated: &[CaseGate],
    nonfinite_gates: &[CaseGate],
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
    for c in candidates {
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

    // per-branch status over this study: cover | finite (taken by a finite candidate outside the
    // cover only) | nonfinite (taken only by NaN-tainted candidates) | open | loop_entry | unreachable
    let study_status = |b: &Branch| -> &'static str {
        if cover_taken.contains(b.id.as_str()) {
            "cover"
        } else if takers_finite.contains_key(b.id.as_str()) {
            "finite"
        } else if takers_all.contains_key(b.id.as_str()) {
            "nonfinite"
        } else {
            match b.status.as_str() {
                "open" => "open",
                "loop_entry" => "loop_entry",
                "unreachable" => "unreachable",
                // taken in the universe run by a case that is not a candidate here (e.g. an
                // untracked --big case): reported as finite-reachable, outside the cover
                _ => "finite",
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
    let mut nonfinite_gated = 0usize;
    for b in &analysis {
        match study_status(b) {
            "unreachable" => *unreachable_by_class.entry(b.class.clone()).or_default() += 1,
            "nonfinite" if b.gate.is_some() => nonfinite_gated += 1,
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
            "cover": count("cover"), "finite": count("finite"), "nonfinite": count("nonfinite"),
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
    let gate_json = |g: &CaseGate, c: &Candidate| -> Value {
        let calls: Vec<Value> = g
            .calls
            .iter()
            .map(|r| {
                json!({
                    "call": r.call, "alpha_deg": r.alpha_deg,
                    "iterations_yfoil": r.iterations.0, "iterations_xfoil": r.iterations.1,
                    "converged_yfoil": r.converged.0, "converged_xfoil": r.converged.1,
                    "gated_iterations": r.gated_iterations, "worst_floor_ratio": r.worst_floor_ratio,
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
            "gated_iterations": g.calls.iter().map(|r| r.gated_iterations).sum::<usize>(),
            "metrics": {
                "cl": point("CL"), "cd": point("CD"), "cm": point("CM"),
                "rmsbl": g.max_transient("RMSBL"), "rlx": g.max_transient("RLX"),
                "worst_floor_ratio": g.worst_floor_ratio(),
                "floor_iteration_max": g.floor_iteration_max, "floor_point_max": g.floor_point_max,
            },
        })
    };
    let nonfinite_runs: Vec<Value> = nonfinite_gates
        .iter()
        .map(|g| gate_json(g, candidates.iter().find(|c| c.case.name == g.name).unwrap()))
        .collect();
    let mut cases = vec![];
    for name in cover {
        let c = candidates.iter().find(|c| &c.case.name == name).unwrap();
        let g = gated.iter().find(|g| &g.name == name).unwrap();
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
                    "gated_iterations": r.gated_iterations, "worst_floor_ratio": r.worst_floor_ratio,
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
            "fixture_host": g.fixture_host, "same_host": g.same_host,
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
    // open and non-finite-only branches of the analysis path
    let mut open = vec![];
    for b in &analysis {
        let st = study_status(b);
        if st == "open" || st == "nonfinite" {
            open.push(json!({
                "id": b.id, "subroutine": b.key(), "line": b.line, "branch": b.branch, "source": b.source,
                "status": st, "class": b.class, "why": b.why, "gate": b.gate,
                "taken_by": takers_all.get(b.id.as_str()).cloned().unwrap_or_default(),
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
            "nonfinite_gated": nonfinite_gated,
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
/// of the claim is the gate, which asserts that on every run of the set yFoil follows XFOIL's
/// branch trace and values.
pub const CLASSIFICATION_MD: &str = "\
### What the statuses mean\n\n\
A *branch* is one outgoing edge of a conditional as gfortran emitted it for the reference build: an `IF` \
gives two edges, every operand of a `.AND.`/`.OR.` chain its own pair, a `DO` line the pair of the loop-control test gfortran placed there (iterate again, or leave the \\
loop by its count), and a \
loop whose back-edge test lands on its last statement a pair there. Which arm is *fallthrough* and which is \
*jump* is the compiler's choice, so an edge is named by line and index, never by true/false; the source line \
and the counts of its sibling edges identify the arm. Every count is gcov's on the **reference**: XFOIL \
6.99, pristine, double precision, compiled with `-fprofile-arcs -ftest-coverage`. yFoil is not instrumented \
and its own control flow is never measured; yFoil's side of the claim is the gate, which asserts that on each \
run of the set it reproduces XFOIL's branch trace (iteration counts, IST/ITRAN, LVCONV, SPECCL's ITAL) and \
values. A branch is therefore *covered* when the reference took it on a run that yFoil matched.\n\n\
The never-covered branches fall into two families. **Logically unreachable** edges cannot be taken by any \
execution of XFOIL's analysis path, whatever the input: the argument is about the reference's own source. \
**Practically unreachable** edges can be taken by XFOIL, but not by any finite input we can construct: \
either the reference has to be in its undefined region (a NaN somewhere in the trajectory) to reach them, or \
a floating-point coincidence has to occur. Neither family is a statement about yFoil.\n\n\
- **covered** — the edge counter is > 0 on at least one run of the `branch-coverage` group, and every run \
  of the group is gated.\n\
- **taken only by finite cases outside the set** — taken on some finite candidate that the cover does not \
  need. Zero here means the set covers everything any finite candidate reaches.\n\
- **non-finite only** (practically unreachable) — taken, but only on runs whose reference output contains \
  `NaN`: XFOIL produced a non-finite value somewhere in that trajectory (in every case found so far, \
  TRCHEK2's transition-point Newton landing on the interval end, after which the MRCHUE/MRCHDU fallback \
  extrapolates finite values over the station and the run carries on). The whole run cannot be gated: the \
  reference's own 1-ULP twins wander by O(10²) on it. What is gated is the **event**: the instrumented \
  reference writes `events.dat` (site, SETBL call, side, station, deciding value) and keeps the state dumps \
  of those calls, and `tests/execution/events.rs` seeds yFoil with XFOIL's exact state entering \
  the call, runs the one subroutine or the one iteration, asserts that yFoil's branch trace records the same \
  event, and compares the result with the reference's dump — NaN where the reference has NaN. The table \
  *Not covered* names the test gating each edge; an edge with a gate is faithful at the call where XFOIL \
  took it, and only the claim that yFoil follows the *whole* non-finite run is not made. The probes are \
  tracked cases (`group = \"non-finite\"`), never in the default selection.\n\
- **open** — never taken on any measured run, finite or not, and not annotated. This is a statement about \
  the search, not about the code: no case has been constructed that reaches the edge. Each open edge \
  carries a *reach* note — the input the code reading says would reach it — and stays open until a \
  finite case takes it.\n\
- **annotated unreachable** — a recorded claim (`xtask/fixtures-config/coverage.toml`, with a class and the \
  reason) that no input reaches the edge. It is a code-reading argument about the reference's own source, \
  never \"we could not build a case\" and never a statement about yFoil. Every claim is falsifiable by the \
  measurement: a case that takes an annotated edge makes the tool report the annotation as stale (this \
  happened to `xoper.f:2784` branch 5, whose claim missed that MRCL's Mach clamp sets `M_CLS = 0`; the \
  claim was withdrawn and the edge is now covered). The first four classes are logically unreachable, the \
  fifth practically:\n\
  - *structural* — impossible by construction of the data the analysis path feeds the routine: coincident \
    consecutive nodes that `LOAD` removes, `IBL ≤ 2` inside an `IF(IBL.GT.3)` block, a quantity kept \
    positive by every update that touches it, a subroutine only ever called under the flag it tests;\n\
  - *mode* — guarded by a flag the analysis path never sets: inverse design (`GEOLIN`), the image airfoil \
    (`LIMAGE`, whose only toggle is commented out), flap hinge moments, interactive prompts. These are \
    features outside yFoil's scope, and the edges are unreachable in XFOIL too when driven as the fixture \
    pipeline drives it;\n\
  - *guard* — array-bound and illegal-input `STOP`s on Fortran's fixed dimensions (`IVX`, `IWX`, …); yFoil \
    has no fixed dimensions and the cases are within XFOIL's;\n\
  - *compiler* — edges gfortran emits that correspond to no decision in the source;\n\
  - *numerical* — reachable in principle, but only by a floating-point coincidence (a real landing bitwise \
    on a threshold: the stagnation point exactly on a node, SETEXP's initial guess exactly 1.0) or by \
    exhausting the iteration cap of a Newton iteration that converges on every analysis-path input \
    (SETEXP's 100, SINVRT's 10, MRCHUE's 25 at the similarity station). No input can be designed to reach \
    them; if anything does, it is the \
    parameter/fuzz harness.\n\
- **DO-line edge** (logically unreachable) — a Fortran `DO` loop is compiled with a control test: *iterate again* or *leave the loop \\
  by its count* (the count exhausted, or zero to begin with). gfortran attributes that test to the `DO` line as \\
  a pair of edges, and the measured counts show two shapes: for a top-tested loop the second arm is the loop's \\
  normal completion, taken once per execution; for a loop gfortran rotated, the real exit sits on the loop's \\
  last statement and the `DO` line's second arm is a duplicate that is never taken. Neither arm is a body \\
  that ran zero times as distinct from a body that ran and finished: the analysis path's station, node and \\
  wake counts never make a loop empty, and the counts cannot tell an empty pass from a completed one. The \\
  *iterate* arm is an ordinary edge, covered whenever a run of the set executed the body, and it appears in \\
  the branch map as a coloured cell beside its partner. This bucket holds the never-taken second arms: a \\
  loop always left early by a `GO TO` or `RETURN` (its normal completion never happens), a rotated loop's \\
  duplicate, or both arms of a loop inside a block that is never reached. Whether each loop's body actually \\
  ran on the set is stated explicitly, from the execution count of the first line of its body, in the table \\
  *DO loops of the analysis path* below. A `DO` line that is an iteration budget (a Newton loop that could \\
  be exhausted: SPECCL's twenty attempts, SPECAL's CL(M) loop) is different: its second arm is the exhausted \\
  budget, a real outcome, so it carries a reach note and both arms are counted like any other conditional.\\n\\n\\
The same classes, with every annotated edge and its reason, are listed in [coverage.md](../coverage.md); \
that page is the union over every tracked case, this one the minimal set and its gates.\n";

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

/// The gate table (difference metrics) as Markdown rows.
fn gate_table_md(s: &Value) -> String {
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
        "| case | OPER | `NaN` in reference output | ref. iterations | gated | outcome | parts at | max \\|ΔRMSBL\\| (gated) | worst diff/floor | ref. 1-ULP spread |\n|---|---|---:|---:|---:|---|---|---:|---:|---:|\n",
    );
    for c in s["nonfinite_runs"].as_array().unwrap() {
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["script"].as_str().unwrap(),
            c["nan_count"],
            c["reference_iterations"],
            c["gated_iterations"],
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
    let gated = u["nonfinite_gated"].as_u64().unwrap_or(0);
    let logical = cls("structural") + cls("mode") + cls("guard") + cls("compiler");
    format!(
        "| | branches |\n|---|---:|\n\
         | translated set (all files) | {} |\n\
         | analysis path (panelling generators excluded) | {} |\n\
         | **covered** — taken by the minimal set, gated | **{}** |\n\
         | taken only by finite cases outside the set | {} |\n\
         | **practically unreachable** — no finite input reaches them | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;taken only in non-finite reference runs, gated as events | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;taken only in non-finite reference runs, not gated | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;*numerical*: a floating-point coincidence or an exhausted Newton cap | {} |\n\
         | **logically unreachable** — no execution of XFOIL's analysis path takes them | {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;*structural* / *mode* / *guard* / *compiler* | {} / {} / {} / {} |\n\
         | &nbsp;&nbsp;&nbsp;&nbsp;never-taken DO-line edges (loop left early or never reached) | {} |\n\
         | **open** — never taken, no case constructed, reach note | {} |\n",
        u["branches_translated_set"],
        u["branches_analysis"],
        st("cover"),
        st("finite"),
        st("nonfinite") + cls("numerical"),
        gated,
        st("nonfinite") - gated,
        cls("numerical"),
        logical + st("loop_entry"),
        cls("structural"),
        cls("mode"),
        cls("guard"),
        cls("compiler"),
        st("loop_entry"),
        st("open")
    )
}

fn open_table_md(s: &Value) -> String {
    let mut t = String::from("| branch | source | status | reach | event gate |\n|---|---|---|---|---|\n");
    for o in s["open"].as_array().unwrap() {
        let taken_by = o["taken_by"].as_array().unwrap();
        let status = match o["status"].as_str().unwrap() {
            "nonfinite" => format!(
                "non-finite only ({})",
                taken_by
                    .iter()
                    .map(|v| format!("`{}`", v.as_str().unwrap()))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            s => s.to_string(),
        };
        let gate = o["gate"]
            .as_str()
            .map(|g| {
                let (file, test) = g.split_once("::").unwrap_or((g, ""));
                format!("`{}` in `{file}`", test)
            })
            .unwrap_or_else(|| "—".into());
        t += &format!(
            "| `{}` | `{}` | {} | {} | {} |\n",
            o["id"].as_str().unwrap(),
            o["source"].as_str().unwrap().replace('|', "\\|"),
            status,
            o["why"].as_str().unwrap_or(""),
            gate
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
    m += "## The universe\n\n";
    m += &format!(
        "The translated set is the {} subroutines of `xtask/fixtures-config/coverage.toml`, measured by \
         gcov on the pristine double-precision reference (`docs/validation/coverage.md`). ",
        s["universe"]["subroutines_translated_set"]
    );
    m += "The analysis path \
          leaves out the panelling generators and their helpers (`naca.f`, `xgeom.f`, `spline.f`, and \
          `NACA`/`PANGEN` in `xfoil.f`): those are gated bitwise by the `pangen_*` cases \
          (`docs/validation/geometry/`). A branch is one outgoing edge of a conditional as gfortran emitted \
          it; the counts below are edges.\n\n";
    m += &status_table_md(s);
    m += "\n";
    m += CLASSIFICATION_MD;
    m += &format!(
        "\n{} tracked solver cases were measured one at a time; {} have a finite reference run and are the \
         candidates of the cover. A case whose reference prints `NaN` is not a candidate: a non-finite \
         trajectory is not reproducible by translation (NaN comparisons and `MIN`/`MAX` are where two \
         correct codes legitimately differ), so a branch reached only that way cannot be gated.\n\n",
        u["candidates"], u["finite_candidates"]
    );
    m += "## The minimal set\n\n";
    m += "Branches taken by exactly one finite candidate make that candidate essential; the remainder is \
          covered greedily by marginal gain and pruned until no case can be removed. The set is recorded as \
          `group = \"branch-coverage\"` in `xtask/fixtures-config/cases.toml` and this study refuses to \
          generate if the two differ. *Branches* is the number of analysis-path edges the case takes; \
          *unique* the ones no other finite candidate takes.\n\n";
    m += &case_table_md(s);
    m += "\nEvery case: yFoil generates the panels (`--method cosine`), XFOIL `LOAD`s them (bitwise handoff \
          asserted), one OPER script, one instrumented run plus its five 1-ULP twins.\n\n";
    if figures {
        m += &format!(
            "![Branch map]({figure_prefix}branch-map.svg)\n\n*Every branch of the analysis path, one row per \
             subroutine in XFOIL's file order: taken by the minimal set (dark), taken only by cases outside it \
             (light), taken only in non-finite reference runs (orange), open (red), annotated unreachable \
             (grey), never-taken DO-line edges — the loop was always left early or never reached (light grey).*\n\n\
             ![Case × subroutine]({figure_prefix}case-matrix.svg)\n\n*Fraction of each subroutine's branches \
             taken by each case of the set; a dot marks a subroutine in which the case takes a branch no other \
             finite candidate takes. The three columns on the right count, per subroutine, every branch, \
             the reachable ones (not annotated unreachable and not a never-taken DO-line edge) and those taken \
             on any measured run, the non-finite probes included.*\n\n"
        );
    }
    m += "## DO loops of the analysis path\n\n";
    m += &do_loop_table_md(s);
    m += "## The gates\n\n";
    m += "Each case of the set is run through yFoil's `Session` exactly as its OPER script drove XFOIL and \
          compared with its fixture by the studies' floor comparison (`scripts/study-support/records.rs`; the tests compare single steps instead): \
          iteration count, LVCONV, IST and ITRAN exact; every per-iteration transient and the converged point \
          within `max(tol · scale, FLOOR_FACTOR · floor)`, where *floor* is the reference's own 1-ULP spread \
          of that value; **divergent** when every earlier iteration matched and the reference \
          itself moves by more than `DIVERGENCE_FLOOR` there. The differences are absolute, maxima over the \
          gated calls and iterations. *Worst diff/floor* is the largest gated difference relative to the \
          reference's own spread of the same value (a value below `FLOOR_FACTOR` = 4 is inside the gate); \
          *ref. 1-ULP spread* is the largest spread the twin showed on the case, per-iteration (viscous) or \
          per-point (inviscid) — when it is large the gate is only as sharp as the reference itself. \
          † inviscid case: the per-node |ΔGAM| maximum (RMSBL does not exist).\n\n";
    m += &gate_table_md(s);
    m += "\n## The non-finite probes, whole run\n\n";
    m += "Not a claim — the reference cannot reproduce itself on these runs — but a measurement of how far \
          yFoil follows each one. The same gate as above is applied call by call: *parts at* is the first \
          iteration at which the runs were allowed to part (the reference's own 1-ULP spread exceeds \
          `DIVERGENCE_FLOOR` there and the values differ) or, for a mismatch, the iteration after the last one \
          gated; *gated* counts the iterations matched before that. The NaN count is the number of `NaN` \
          tokens in the reference's own output for the run.\n\n";
    m += &nonfinite_table_md(s);
    m += "\n## Full polars\n\nThe same sections swept as full ±30° polars, through stall and into the non-finite region, are \
          the [series-cases](../series-cases/README.md) study. Those sweeps are measured here like \
          every tracked case (their branches count as taken) but are never candidates for the cover: the cover \
          is the minimal set of single operating points, and the polars are built from its sections.\n";
    m += "\n## Not covered\n\n";
    m += "Analysis-path branches no case of the set takes, with the recorded way to reach them. *Non-finite \
          only* branches were taken during this study's candidate search, but only in runs where the reference \
          had produced a NaN (the stalled compressible and closed-TE stall probes); each is gated as an event \
          by the test named in the last column (`events.dat` of the probe, replayed from the dumped state of \
          that SETBL call).\n\n";
    m += &open_table_md(s);
    m += "\n## Candidates excluded for a non-finite reference run\n\n";
    m += &excluded_md(s);
    m
}

pub fn write_index(run_dir: &Path, s: &Value, _universe: &Universe, _candidates: &[Candidate]) {
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let mut m = String::from("# Branch-coverage study\n\n");
    m += &format!(
        "Generated by `cargo run --release -p branch-coverage` on {} (git {}{}, host {}); provenance in \
         `metadata.json`, every number in `summary.json`, figures drawn by `plot.py`.\n\n",
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
    w(&mut f, "\\noindent The minimal set of cases taking every finite-reachable branch of XFOIL's translated analysis path, each gated against its instrumented double-precision reference. Provenance in \\texttt{metadata.json}, every number in \\texttt{summary.json}.");
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
    // gate table
    w(&mut f, "\\begin{table}[h]\\centering\\small");
    w(&mut f, "\\caption{The gates: every case of the set against its fixture. Differences are absolute maxima over the gated calls and iterations; \\emph{worst diff/floor} is relative to the reference's own 1-ULP spread of the same value; $\\dagger$ inviscid case, per-node $|\\Delta\\mathrm{GAM}|$ in place of RMSBL.}");
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
    let gated = u["nonfinite_gated"].as_u64().unwrap_or(0);
    w(
        &mut f,
        &format!(
            "Covered: taken by the minimal set, gated & \\textbf{{{}}} \\\\",
            st("cover")
        ),
    );
    w(
        &mut f,
        &format!("Taken only by finite cases outside the set & {} \\\\", st("finite")),
    );
    w(
        &mut f,
        &format!(
            "Practically unreachable: non-finite runs only, gated as events / not gated & {} / {} \\\\",
            gated,
            st("nonfinite") - gated
        ),
    );
    w(
        &mut f,
        &format!(
            "Practically unreachable: numerical (coincidence or exhausted cap) & {} \\\\",
            cls("numerical")
        ),
    );
    w(
        &mut f,
        &format!(
            "Logically unreachable: structural / mode / guard / compiler & {} / {} / {} / {} \\\\",
            cls("structural"),
            cls("mode"),
            cls("guard"),
            cls("compiler")
        ),
    );
    w(&mut f, &format!("Open (never taken, reach note) & {} \\\\", st("open")));
    w(&mut f, &format!("Never-taken DO-line edges (loop always left early or never reached) & {} \\\\ \\bottomrule\\end{{tabular}}\\end{{table}}", st("loop_entry")));
    w(&mut f, "\\noindent Every count is gcov's on the pristine double-precision reference; yFoil is never instrumented, and a branch is \\emph{covered} when the reference took it on a run of the set that yFoil's gate matched. \\emph{Non-finite only}: taken only on runs whose reference output contains NaN --- the run is not gateable and not counted, the event is (replayed from XFOIL's state at that call, NaN where XFOIL has NaN). \\emph{Open}: never taken on any measured run and not annotated --- no case has been constructed for it, each carries a reach note. \\emph{Annotated unreachable}: a recorded, falsifiable code-reading claim that no execution of XFOIL's analysis path can take the edge (structural, mode, guard, compiler or numerical: a floating-point coincidence or an exhausted Newton cap). \\emph{DO-line edge}: the never-taken arm of the loop-control test on a DO line (the loop always left early by a GO TO/RETURN, a rotated loop's duplicate test, or a loop never reached); the iterate arm is an ordinary covered edge, and whether every loop's body ran on the set is tabulated from the body's execution count. The Markdown page carries the full definitions.");
    w(&mut f, "\\begin{figure}[h]\\centering\\includegraphics{branch-map.pdf}\\caption{Every branch of the analysis path, one row per subroutine in XFOIL's file order: taken by the minimal set (dark), taken only by cases outside it (light), taken only in non-finite reference runs (orange), open (red), annotated unreachable (grey), never-taken DO-line edges --- the loop always left early or never reached (light grey).}\\end{figure}");
    w(&mut f, "\\begin{figure}[h]\\centering\\includegraphics{case-matrix.pdf}\\caption{Fraction of each subroutine's branches taken by each case of the set; a dot marks a subroutine in which the case takes a branch no other finite candidate takes; the columns on the right count every branch of the subroutine, the reachable ones (not annotated unreachable, not a never-taken DO-line edge) and those taken on any measured run, the non-finite probes included.}\\end{figure}");
}

pub fn write_docs(docs: &Path, run_dir: &Path, s: &Value, universe: &Universe, candidates: &[Candidate], drawn: bool) {
    let _ = (universe, candidates);
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let mut m = String::from("# Branch-coverage study: the minimal case set\n\n");
    m += &format!(
        "Generated by `cargo run --release -p branch-coverage -- --docs` on {} (git {}{}, host {}). Do not edit; \
         regenerate. The run folder (`scripts/branch-coverage/runs/`, gitignored) holds `summary.json` with every \
         number on this page and `README.tex` with the same tables for the paper.\n\n",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        if meta["git_dirty"].as_bool().unwrap_or(false) {
            ", dirty"
        } else {
            ""
        },
        meta["host"].as_str().unwrap_or("")
    );
    m += "**Question.** CLAUDE.md Rule 6: \"numerically exact\" is claimable only over branches that have been \
          exercised. [coverage.md](../coverage.md) measures the union over every tracked case; this page asks \
          for the *minimum* set of cases that takes every reachable branch of the analysis path, and compares every \
          case of that set, run as a whole, with its instrumented reference — so that the coverage claim and the \
          equivalence claim are made over the same runs. This is a study: the tests gate each branch at a single \
          step whose route yFoil is observed to follow ([branch-gating.md](../branch-gating.md), \
          `docs/conventions/testing.md`), and this page adds no evidence they rely on.\n\n";
    m += &body_md(s, drawn, "");
    fs::write(docs.join("README.md"), m).unwrap();
    for name in ["branch-map.svg", "case-matrix.svg"] {
        if run_dir.join(name).exists() {
            fs::copy(run_dir.join(name), docs.join(name)).unwrap();
        }
    }
}
