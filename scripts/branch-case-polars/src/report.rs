//! The run's `summary.json`, its `README.md` / `README.tex` index, and the generated
//! documentation page (`docs/validation/branch-case-polars/README.md`).

use crate::polars::PolarRun;
use serde_json::{json, Value};
use std::fs;
use std::io::Write;
use std::path::Path;

pub fn summary(runs: &[PolarRun]) -> Value {
    // one figure per section family, in order of first appearance
    let mut families: Vec<(String, String)> = vec![];
    for r in runs {
        if !families.iter().any(|(f, _)| *f == r.family()) {
            families.push((r.family(), r.section.replace(" (sharp TE)", "")));
        }
    }
    json!({
        "nseqex": crate::polars::NSEQEX,
        "families": families.iter().map(|(slug, name)| json!({ "slug": slug, "section": name })).collect::<Vec<_>>(),
        "polars": runs.iter().map(|r| r.to_json()).collect::<Vec<_>>(),
    })
}

fn parted(r: &Value) -> String {
    match r["parted"].as_array() {
        Some(a) if a.len() == 2 => format!("leg {} α {}°", a[0], a[1]),
        _ => "—".into(),
    }
}

/// One row per polar: what each code recorded and where they part.
fn table_md(s: &Value) -> String {
    let mut t = String::from(
        "| case | section | OPER | XFOIL points | yFoil points | followed to the end | parts at | reference hung |\n|---|---|---|---:|---:|---:|---|---|\n",
    );
    for r in s["polars"].as_array().unwrap() {
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} |\n",
            r["name"].as_str().unwrap(),
            r["section"].as_str().unwrap(),
            r["script"].as_str().unwrap(),
            r["xfoil_calls"],
            r["yfoil_points"],
            r["matched"],
            parted(r),
            if r["truncated_by_watchdog"].as_bool().unwrap_or(false) {
                "yes"
            } else {
                "no"
            }
        );
    }
    t
}

fn body_md(s: &Value, figure: bool) -> String {
    let mut m = String::new();
    m += "The sections and conditions the branch-coverage study drew its cases from \
          ([branch-coverage](../branch-coverage/README.md)), each swept as a full polar: 0 → +30°, `INIT`, \
          0 → −30° by 1°, XFOIL's polar procedure (`ALFA 0 / ASEQ`, ITMAX + 5 iterations per point, the sequence \
          halting after NSEQEX = 4 consecutive unconverged points). yFoil is driven the same way through its \
          `Session`, applying the same halting rule to its own convergence, and every point the reference \
          recorded is compared with the studies' whole-run floor comparison (`scripts/study-support/records.rs`): \
          iteration count, LVCONV, IST and ITRAN exact, every transient and the converged point within \
          `max(tol · scale, FLOOR_FACTOR · floor)` with *floor* the reference's own 1-ULP spread of the value, \
          a *threshold* (the divergent comparison of `docs/conventions/terminology.md`) where the reference itself moves by more than \
          `DIVERGENCE_FLOOR`.\n\n";
    m += "The point of the study is what happens past convergence: through stall, and into the region where the \
          reference goes non-finite. The reference may record fewer points than the script — its sequence \
          halted, or the fixture driver's watchdog ended a run that had gone non-finite and hung in XFOIL's \
          plot-label loop ([known issues §4](../../xfoil-known-issues.md)); yFoil's sweep continues to its own \
          halt either way. *Followed to the end* counts the recorded points at which every iteration was within \
          tolerance (a match, or a departure allowed only because one of the reference's own 1-ULP twins had \
          already departed earlier in the sweep); *parts at* is the first recorded point (leg, α) yFoil did not \
          follow to the end. A departure where the reference is unstable against its 1-ULP twins (the largest \
          spread of the value over the five twins above `DIVERGENCE_FLOOR` there) is a *threshold*; a departure \
          where it is stable is a \
          *mismatch*, settled by the one-step replay of CLAUDE.md Rule 3 — seed yFoil with XFOIL's exact state \
          entering that SETBL call and require XFOIL's state after it at `TOL_SOLVER`. \
          `tests/execution/steps.rs` replays every iteration of the calls where these sweeps part that way \
          (`step_calls` in `cases.toml`), each value within the named tolerance or four times the step's own \
          measured sensitivity, and `cargo xtask route` observes yFoil take XFOIL's route through each of those \
          steps. Why a sweep parts where the reference is stable against its twins — MRCHDU's station Newton at \
          a transition station amplifying its seed a thousandfold, which both codes follow iterate for iterate \
          from the same state — is recorded with its measurements in \
          [known issues §7.9](../../xfoil-known-issues.md).\n\n";
    m += "One figure per case: CL and CD against α (top), the last iteration's RMSBL and the iterations taken \
          (second row), and the transition location x/c on the upper and lower side (third row, from XOCTR and \
          yFoil's transition point; an open triangle marks a point whose final stagnation or transition station — \
          IST, ITRAN — differs between the codes), and the position of transition within its station interval on \
          each side (fourth row: 0 at station ITRAN − 1, 1 at station ITRAN; see the section on it below). For \
          every case — the reference writes its final state \
          of every VISCAL call (the per-node Cp, q and γ and the per-station XSSI, x, \
          UEDG, THET, DSTR, CTAU, MASS and TAU, `viscal_state_<k>.dat`) a fifth row compares yFoil's final state \
          of the point with the reference's array by array: the worst array's |Δ| as a ratio to its gate \
          `max(TOL_SOLVER · scale, FLOOR_FACTOR · floor)` with the floor the reference's own 1-ULP twin's spread \
          of that array at that call (left; the dotted line is the gate), and the |Δ| itself against that floor \
          (right). That row is the accumulation measurement: a state that stays at the floor from point to point \
          and then leaves it within one point has not accumulated anything. **XFOIL is red** (dashed, crosses), **yFoil blue** (solid, circles — filled where yFoil \
          followed the reference to the end of the point, hollow where it did not or the point was not gated), \
          so any visible separation of the two colours is a difference between the codes. An orange outer square \
          marks a point during which that code went through a non-finite value: a station Newton whose residual \
          was NaN (the garbage extrapolation then carries the march over it — for the reference, any `NaN` it \
          printed during the call; for yFoil, its own count of such stations) or a non-finite RMSBL, CL, CD or CM \
          in any iteration; where the plotted value is itself NaN and cannot be placed, an orange tick at the top \
          of the panel marks its α. A black outer diamond marks a point within which the run departs from the \
          reference where the reference is **unstable against its 1-ULP twins** (one of its own runs with every \
          panel coordinate moved by one ULP departs there too: a *threshold*, the outcome the tests call \
          divergent, CLAUDE.md Rule 1; the one-step replay from XFOIL's state at that iteration is the \
          evidence that the step is faithful). A filled grey diamond marks a point within which the run departs \
          where the reference is **stable against its 1-ULP twins**: a *mismatch*, to be settled by the same \
          replay. A point yFoil \
          followed to the end carries no diamond even when an earlier departure left it ungated; a point after \
          a departure within the same leg that yFoil did not follow carries none either (its comparison is \
          informational: it started from a state the reference did not compute). Where the red line stops short of the blue one the reference's sequence halted or hung \
          there.\n\n";
    if figure {
        let mut last_family = String::new();
        for r in s["polars"].as_array().unwrap() {
            let family = r["section"].as_str().unwrap().replace(" (sharp TE)", "");
            if family != last_family {
                m += &format!("### {family}\n\n");
                last_family = family;
            }
            m += &format!(
                "**`{}`** — {}\n\n![{}](polars-{}.svg)\n\n",
                r["name"].as_str().unwrap(),
                r["description"].as_str().unwrap_or(""),
                r["name"].as_str().unwrap(),
                r["name"].as_str().unwrap()
            );
        }
    }
    m += &table_md(s);
    if figure {
        m += &transition_md(s);
    }
    let st = state_table_md(s);
    if !st.is_empty() {
        m += "\n### The final state, point by point, into each departure\n\n";
        m += "For the cases that keep the reference's per-call state records: every recorded point at which yFoil \
              parts (the first not followed to the end in its leg) and the points of the leg before it, with the \
              worst array of yFoil's final state against the reference's — its worst entry (node number for a \
              node array, `BL(is, ibl)` row index for a station array, side 1 first), the |Δ| there, the twin's \
              1-ULP floor of the array, and the ratio of |Δ| to the gate. A point followed to the end that ends \
              at a ratio well below 1 handed the next point a state the reference's own twin could not tell from \
              the reference's. Where the run parts, `tests/execution/steps.rs` replays every iteration of the call \
              from the reference's exact state (Rule 3), with each value within the named tolerance or four times \
              the step's own measured sensitivity.\n\n";
        m += &st;
    }
    m
}

/// The rows of the state table: for every case that kept the reference's per-call state records,
/// each recorded point at which yFoil parts and the points of its leg before it, with the worst
/// array of the final state.
fn state_rows(s: &Value) -> Vec<Value> {
    let mut rows = vec![];
    for r in s["polars"].as_array().unwrap() {
        let pts: Vec<&Value> = r["points"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|p| p["state"].is_array())
            .collect();
        if pts.is_empty() {
            continue;
        }
        for leg in [1, 2] {
            let leg_pts: Vec<&&Value> = pts.iter().filter(|p| p["leg"] == leg).collect();
            // up to and including the first recorded point of the leg not followed to the end
            let end = leg_pts
                .iter()
                .position(|p| !p["followed"].as_bool().unwrap_or(false))
                .map(|i| i + 1)
                .unwrap_or(leg_pts.len());
            let start = end.saturating_sub(4);
            for p in &leg_pts[start..end] {
                let w = &p["state"][0];
                rows.push(json!({
                    "case": r["name"], "leg": leg, "alpha_deg": p["alpha_deg"], "followed": p["followed"],
                    "outcome": p["outcome"], "iterations": [p["yfoil"]["iterations"], p["xfoil"]["iterations"]],
                    "array": w["array"], "index": w["index"], "diff": w["diff"], "floor": w["floor"], "ratio": w["ratio"],
                }));
            }
        }
    }
    rows
}

fn state_table_md(s: &Value) -> String {
    let rows = state_rows(s);
    if rows.is_empty() {
        return String::new();
    }
    let mut t = String::from(
        "| case | leg | α | iterations yFoil / XFOIL | followed | worst array (entry) | abs Δ | 1-ULP floor | Δ / gate |\n|---|---:|---:|---:|---|---|---:|---:|---:|\n",
    );
    for r in &rows {
        t += &format!(
            "| `{}` | {} | {} | {} / {} | {} | {} ({}) | {:.1e} | {:.1e} | {:.2} |\n",
            r["case"].as_str().unwrap(),
            r["leg"],
            r["alpha_deg"],
            r["iterations"][0],
            r["iterations"][1],
            if r["followed"].as_bool().unwrap_or(false) {
                "yes"
            } else {
                r["outcome"].as_str().unwrap_or("no")
            },
            r["array"].as_str().unwrap(),
            r["index"],
            r["diff"].as_f64().unwrap_or(f64::NAN),
            r["floor"].as_f64().unwrap_or(f64::NAN),
            r["ratio"].as_f64().unwrap_or(f64::NAN),
        );
    }
    t
}

pub fn write_index(run_dir: &Path, s: &Value) {
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let mut m = String::from("# Branch-case polars\n\n");
    m += &format!(
        "Generated by `cargo run --release -p branch-case-polars` on {} (git {}{}, host {}); provenance in \
         `metadata.json`, every number in `summary.json`, the figure drawn by `plot.py`.\n\n",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        if meta["git_dirty"].as_bool().unwrap_or(false) {
            ", dirty"
        } else {
            ""
        },
        meta["host"].as_str().unwrap_or("")
    );
    m += &body_md(s, true);
    fs::write(run_dir.join("README.md"), m).unwrap();
}

fn tex_escape(s: &str) -> String {
    s.replace('\\', "\\textbackslash{}")
        .replace('&', "\\&")
        .replace('%', "\\%")
        .replace('#', "\\#")
        .replace('_', "\\_")
        .replace('α', "$\\alpha$")
        .replace('—', "---")
        .replace('→', "$\\rightarrow$")
        .replace('°', "$^\\circ$")
}

pub fn write_index_tex(run_dir: &Path, s: &Value) {
    let mut f = fs::File::create(run_dir.join("README.tex")).unwrap();
    let w = |f: &mut fs::File, t: &str| writeln!(f, "{t}").unwrap();
    w(
        &mut f,
        "% Generated by `cargo run --release -p branch-case-polars`; table from summary.json. Needs booktabs.",
    );
    w(&mut f, "\\section*{Branch-case polars}");
    w(&mut f, "\\noindent The branch-coverage sections swept as full polars, $0\\rightarrow+30^\\circ$ and after INIT $0\\rightarrow-30^\\circ$ by $1^\\circ$ with XFOIL's polar procedure (ASEQ, ITMAX+5 per point, the sequence halting after NSEQEX = 4 consecutive unconverged points); yFoil driven the same way and compared point by point with the equivalence gate. Provenance in \\texttt{metadata.json}, every number in \\texttt{summary.json}.");
    for r in s["polars"].as_array().unwrap() {
        w(&mut f, &format!("\\begin{{figure}}[h]\\centering\\includegraphics{{polars-{}.pdf}}\\caption{{{} (\\texttt{{{}}}): $C_L$ and $C_D$ against $\\alpha$ (top), the last iteration's RMSBL and the iterations taken (second row), the transition location on each side (third row; an open triangle where the final IST or ITRAN differs between the codes), the position of transition within its station interval on each side (fourth row) and, where the fixture keeps the reference's per-call state, yFoil's final state against the reference's, worst array, as a ratio to its gate and as $|\\Delta|$ against the twin's 1-ULP floor (fifth row). XFOIL red (dashed, crosses), yFoil blue (solid, circles filled where yFoil followed the reference to the end of the point); an orange outer square marks a point in which that code went through a non-finite value, an orange tick a plotted value that is itself NaN; a black outer diamond a point within which the run departs where the reference is unstable against its 1-ULP twins (a threshold), a filled grey diamond one where it is stable (a mismatch). Where the red line stops short of the blue one the reference's sequence halted or hung there.}}\\end{{figure}}", r["name"].as_str().unwrap(), tex_escape(r["section"].as_str().unwrap()), tex_escape(r["name"].as_str().unwrap())));
    }
    w(&mut f, "\\begin{figure}[h]\\centering\\includegraphics{transition-position.pdf}\\caption{Where transition sits within its station interval (0 at station ITRAN$-$1, 1 at station ITRAN) at the end of every recorded point of every case, per side, from yFoil's state; diamonds where the run parts (black: the reference is unstable against its 1-ULP twins there; grey: it is stable).}\\end{figure}");
    w(&mut f, "\\begin{table}[h]\\centering\\small");
    w(&mut f, "\\caption{Points each code recorded, points yFoil followed to the end, and the first recorded point it did not.}");
    w(&mut f, "\\begin{tabular}{@{}llrrrll@{}}\\toprule");
    w(
        &mut f,
        "Case & Section & XFOIL & yFoil & Followed & Parts at & Hung \\\\ \\midrule",
    );
    for r in s["polars"].as_array().unwrap() {
        w(
            &mut f,
            &format!(
                "\\texttt{{{}}} & {} & {} & {} & {} & {} & {} \\\\",
                tex_escape(r["name"].as_str().unwrap()),
                tex_escape(r["section"].as_str().unwrap()),
                r["xfoil_calls"],
                r["yfoil_points"],
                r["matched"],
                tex_escape(&parted(r)),
                if r["truncated_by_watchdog"].as_bool().unwrap_or(false) {
                    "yes"
                } else {
                    "no"
                }
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    let rows = state_rows(s);
    if !rows.is_empty() {
        w(&mut f, "\\begin{table}[h]\\centering\\small");
        w(&mut f, "\\caption{The final state into each departure: the worst array of yFoil's final state against the reference's, its $|\\Delta|$, the twin's 1-ULP floor of the array and the ratio of $|\\Delta|$ to the gate, for every point at which yFoil parts and the points of its leg before it.}");
        w(&mut f, "\\begin{tabular}{@{}lrrrllrrr@{}}\\toprule");
        w(&mut f, "Case & Leg & $\\alpha$ & Iterations & Followed & Array (entry) & $|\\Delta|$ & Floor & $|\\Delta|$/gate \\\\ \\midrule");
        for r in &rows {
            w(
                &mut f,
                &format!(
                    "\\texttt{{{}}} & {} & {} & {} / {} & {} & {} ({}) & {:.1e} & {:.1e} & {:.2} \\\\",
                    tex_escape(r["case"].as_str().unwrap()),
                    r["leg"],
                    r["alpha_deg"],
                    r["iterations"][0],
                    r["iterations"][1],
                    if r["followed"].as_bool().unwrap_or(false) {
                        "yes"
                    } else {
                        r["outcome"].as_str().unwrap_or("no")
                    },
                    r["array"].as_str().unwrap(),
                    r["index"],
                    r["diff"].as_f64().unwrap_or(f64::NAN),
                    r["floor"].as_f64().unwrap_or(f64::NAN),
                    r["ratio"].as_f64().unwrap_or(f64::NAN),
                ),
            );
        }
        w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    }
}

pub fn write_docs(docs: &Path, run_dir: &Path, s: &Value, drawn: bool) {
    let meta: Value = serde_json::from_str(&fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let mut m = String::from("# Branch-case polars: the branch-coverage sections through stall\n\n");
    m +=
        &format!(
        "Generated by `cargo run --release -p branch-case-polars -- --docs` on {} (git {}{}, host {}). Do not edit; \
         regenerate. The run folder (`scripts/branch-case-polars/runs/`, gitignored) holds `summary.json` with \
         every point of every polar and `README.tex` with the same table for the paper.\n\n",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        if meta["git_dirty"].as_bool().unwrap_or(false) { ", dirty" } else { "" },
        meta["host"].as_str().unwrap_or("")
    );
    m += &body_md(s, drawn);
    fs::write(docs.join("README.md"), m).unwrap();
    for old in fs::read_dir(docs).unwrap().flatten() {
        if old.file_name().to_string_lossy().ends_with(".svg") {
            fs::remove_file(old.path()).unwrap();
        }
    }
    for r in s["polars"].as_array().unwrap() {
        let name = format!("polars-{}.svg", r["name"].as_str().unwrap());
        if run_dir.join(&name).exists() {
            fs::copy(run_dir.join(&name), docs.join(&name)).unwrap();
        }
    }
    if run_dir.join("transition-position.svg").exists() {
        fs::copy(
            run_dir.join("transition-position.svg"),
            docs.join("transition-position.svg"),
        )
        .unwrap();
    }
}

/// The transition position of every recorded point at which a run parts: (case, leg, α,
/// outcome, side, ITRAN, fraction within the interval, distance to the nearer station).
fn parted_transition_rows(s: &Value) -> Vec<Value> {
    let mut rows = vec![];
    for r in s["polars"].as_array().unwrap() {
        for p in r["points"].as_array().unwrap() {
            let parts = p["xfoil"].is_object()
                && !p["followed"].as_bool().unwrap_or(false)
                && p["parted_iteration"].as_u64().unwrap_or(0) > 0;
            if !parts {
                continue;
            }
            for (side, name) in [(0usize, "upper"), (1, "lower")] {
                let t = &p["transition_position"][side];
                if t.is_null() || t["at_te"].as_bool().unwrap_or(false) {
                    continue;
                }
                rows.push(json!({
                    "case": r["name"], "leg": p["leg"], "alpha_deg": p["alpha_deg"], "outcome": p["outcome"],
                    "side": name, "i_station": t["i_station"], "fraction": t["fraction"], "nearest": t["nearest"],
                }));
            }
        }
    }
    rows
}

fn transition_md(s: &Value) -> String {
    let mut m = String::from("\n### Where transition sits between the stations\n\n");
    m += "XFOIL places transition at an arc length XSSITR inside the interval between two consecutive BL \
          stations, ITRAN − 1 and ITRAN (the first turbulent station), and TRDIF treats the transition interval \
          with the two states weighted by where XSSITR falls in it. The figure plots, for every recorded point of \
          every case, the position of transition within that interval at the end of the point — 0 at the station \
          before, 1 at the transition station — on each side, from yFoil's state (the reference agrees in ITRAN \
          at every point followed to the end); the points at which the run parts are the diamonds. Departures \
          clustered at 0, 1 or ½ would reveal a sensitivity to where the panelling puts the stations relative to \
          transition; the per-case figures carry the same quantity against α in their fourth row. Points whose \
          ITRAN is the trailing-edge station (no transition on the airfoil) are left out.\n\n\
          ![transition position](transition-position.svg)\n\n";
    let rows = parted_transition_rows(s);
    if !rows.is_empty() {
        m += "| case | leg | α | outcome | side | ITRAN | position in interval | distance to nearer station |\n\
              |---|---:|---:|---|---|---:|---:|---:|\n";
        for r in &rows {
            m += &format!(
                "| `{}` | {} | {} | {} | {} | {} | {:.3} | {:.3} |\n",
                r["case"].as_str().unwrap(),
                r["leg"],
                r["alpha_deg"],
                r["outcome"].as_str().unwrap_or(""),
                r["side"].as_str().unwrap(),
                r["i_station"],
                r["fraction"].as_f64().unwrap_or(f64::NAN),
                r["nearest"].as_f64().unwrap_or(f64::NAN),
            );
        }
    }
    m
}
