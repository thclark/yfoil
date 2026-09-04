//! `cargo xtask coverage [--case NAME]... [--big] [--rebuild]`
//!
//! Rule 6 completeness: build the reference with gcov instrumentation, run every case in
//! `xtask/fixtures-config/cases.toml` through it (from the case's tracked `xfoil.inp`/`panels.dat`), read
//! the accumulated branch counters back with gcov, and report — per translated subroutine —
//! every branch that was never taken. `xtask/fixtures-config/coverage.toml` names the translated set and
//! carries the annotations for branches that are unreachable on the analysis path; anything
//! never taken and not annotated is *open*, and the fixture set is complete when the open
//! list is empty. The report is `docs/validation/coverage.md`.

use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

#[derive(Deserialize)]
struct Spec {
    #[serde(rename = "file")]
    files: Vec<FileSpec>,
    #[serde(default, rename = "unreachable")]
    unreachable: Vec<Unreachable>,
    /// open branches with a recorded way to reach them (kept with the open list, not counted as covered)
    #[serde(default, rename = "note")]
    notes: Vec<Note>,
    /// subroutines never called on the analysis path, with the reason
    #[serde(default, rename = "dead")]
    dead: Vec<DeadSub>,
}

#[derive(Deserialize, Clone)]
struct Note {
    file: String,
    line: usize,
    #[serde(default)]
    branch: Option<usize>,
    reach: String,
}

#[derive(Deserialize, Clone)]
struct DeadSub {
    file: String,
    subroutine: String,
    /// true = the subroutine is reachable on the analysis path and the entry only records why
    /// no case reaches it yet (still counts as unexplained)
    #[serde(default)]
    open: bool,
    why: String,
}

#[derive(Deserialize)]
struct FileSpec {
    name: String,
    subroutines: Vec<String>,
}

#[derive(Deserialize, Clone)]
struct Unreachable {
    file: String,
    line: usize,
    /// branch index on that line (omit = every never-taken branch on the line)
    #[serde(default)]
    branch: Option<usize>,
    /// structural (impossible by construction), mode (a feature outside the analysis path),
    /// guard (array-bound / illegal-input STOP), compiler (no source-level decision)
    class: String,
    why: String,
}

// gcov --json-format (GCC >= 9)
#[derive(Deserialize)]
struct GcovJson {
    gcc_version: String,
    files: Vec<GcovFile>,
}
#[derive(Deserialize)]
struct GcovFile {
    file: String,
    functions: Vec<GcovFn>,
    lines: Vec<GcovLine>,
}
#[derive(Deserialize)]
struct GcovFn {
    name: String,
    start_line: usize,
    end_line: usize,
    execution_count: u64,
}
#[derive(Deserialize)]
struct GcovLine {
    line_number: usize,
    #[serde(default)]
    function_name: Option<String>,
    #[serde(default)]
    branches: Vec<GcovBranch>,
}
#[derive(Deserialize)]
struct GcovBranch {
    count: u64,
    fallthrough: bool,
}

/// One never-taken branch, classified.
struct Dead {
    line: usize,
    branch: usize,
    fallthrough: bool,
    /// counts of the other branches on the same line (what *was* taken there)
    siblings: Vec<u64>,
    source: String,
    kind: Kind,
    /// unreachable: the class; open: the recorded way to reach it
    class: String,
    why: Option<String>,
}

#[derive(PartialEq, Clone, Copy)]
enum Kind {
    /// a DO loop's zero-trip entry edge (the loop body always ran at least once)
    LoopEntry,
    /// annotated in xtask/fixtures-config/coverage.toml
    Unreachable,
    /// never taken, no annotation — the objective for more cases / the fuzz harness
    Open,
}

struct SubReport {
    file: String,
    name: String,
    calls: u64,
    branches: usize,
    taken: usize,
    dead: Vec<Dead>,
}

pub(crate) fn run(flags: &[String]) {
    let root = super::root();
    let big = flags.iter().any(|f| f == "--big");
    let rebuild = flags.iter().any(|f| f == "--rebuild");
    let selected: Vec<&str> = flags
        .windows(2)
        .filter(|w| w[0] == "--case")
        .map(|w| w[1].as_str())
        .collect();
    let gdir = root.join("target/xfoil-ref/gcov");
    let bin = gdir.join("bin");
    if rebuild || !bin.join("xfoil").exists() {
        super::run(
            Command::new(root.join("scripts/xfoil-build.sh")).arg("--gcov"),
            "scripts/xfoil-build.sh --gcov",
        );
    }
    // reset the counters: one measurement = one set of cases
    for e in fs::read_dir(&bin).unwrap().flatten() {
        if e.path().extension().is_some_and(|x| x == "gcda") {
            fs::remove_file(e.path()).unwrap();
        }
    }

    // run the cases from their tracked inputs
    let cases: super::Cases =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("xtask/fixtures-config/cases.toml"))
            .expect("parse cases.toml");
    let mut ran: Vec<String> = vec![];
    let mut failures = 0;
    for case in &cases.cases {
        if !selected.is_empty() && !selected.contains(&case.name.as_str()) {
            continue;
        }
        if !case.track && !big && selected.is_empty() {
            continue;
        }
        let src = if case.track {
            root.join("tests/fixtures/xfoil").join(&case.name)
        } else {
            root.join("target/fixtures").join(&case.name)
        };
        if !src.join("xfoil.inp").exists() {
            eprintln!(
                "  {}: no xfoil.inp under {} — run `cargo xtask fixtures --case {}` first",
                case.name,
                src.display(),
                case.name
            );
            failures += 1;
            continue;
        }
        let work = root.join("target/coverage").join(&case.name);
        let _ = fs::remove_dir_all(&work);
        fs::create_dir_all(&work).unwrap();
        for f in ["xfoil.inp", "panels.dat"] {
            if src.join(f).exists() {
                fs::copy(src.join(f), work.join(f)).unwrap();
            }
        }
        let inp = fs::File::open(work.join("xfoil.inp")).unwrap();
        let out = fs::File::create(work.join("stdout.txt")).unwrap();
        let st = Command::new(bin.join("xfoil"))
            .current_dir(&work)
            .stdin(inp)
            .stdout(out)
            .stderr(Stdio::inherit())
            .status()
            .expect("run gcov xfoil");
        if !st.success() {
            eprintln!("  {}: xfoil exited {st}", case.name);
            failures += 1;
            continue;
        }
        println!("  ran {}", case.name);
        ran.push(case.name.clone());
    }
    assert!(!ran.is_empty(), "no cases ran");

    // read the counters back
    let spec: Spec =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/coverage.toml")).expect("xtask/fixtures-config/coverage.toml"))
            .expect("parse coverage.toml");
    let gcov = gcov_binary();
    let scratch = root.join("target/coverage/_gcov");
    let _ = fs::remove_dir_all(&scratch);
    fs::create_dir_all(&scratch).unwrap();
    let mut gcc_version = String::new();
    let mut reports: Vec<SubReport> = vec![];
    let mut never_called: Vec<(String, Option<(bool, String)>)> = vec![];
    let mut missing: Vec<String> = vec![];
    for fspec in &spec.files {
        let src_path = gdir.join("src").join(&fspec.name);
        let source: Vec<String> = fs::read_to_string(&src_path)
            .unwrap_or_else(|_| panic!("read {}", src_path.display()))
            .lines()
            .map(str::to_string)
            .collect();
        let st = Command::new(&gcov)
            .args(["-j", "-b", "-o"])
            .arg(&bin)
            .arg(&src_path)
            .current_dir(&scratch)
            .stdout(Stdio::null())
            .status()
            .unwrap_or_else(|e| panic!("run {gcov}: {e}"));
        assert!(st.success(), "{gcov} failed on {}", fspec.name);
        // gcov names its output after the object stem (xpanel.gcda -> xpanel.gcov.json.gz)
        let stem = fspec.name.trim_end_matches(".f");
        let gz = scratch.join(format!("{stem}.gcov.json.gz"));
        let json = Command::new("gzip").args(["-dc"]).arg(&gz).output().expect("gzip -dc");
        assert!(json.status.success(), "gzip -dc {} failed", gz.display());
        let parsed: GcovJson = serde_json::from_slice(&json.stdout).expect("parse gcov json");
        gcc_version = parsed.gcc_version.clone();
        let gf = parsed
            .files
            .iter()
            .find(|f| Path::new(&f.file).file_name() == src_path.file_name())
            .unwrap_or_else(|| panic!("no entry for {} in gcov json", fspec.name));
        // lines by function
        let mut by_fn: HashMap<&str, Vec<&GcovLine>> = HashMap::new();
        for l in &gf.lines {
            if let Some(f) = &l.function_name {
                by_fn.entry(f.as_str()).or_default().push(l);
            }
        }
        for sub in &fspec.subroutines {
            let sym = format!("{}_", sub.to_lowercase());
            let Some(func) = gf.functions.iter().find(|f| f.name == sym) else {
                missing.push(format!("{}:{sub}", fspec.name));
                continue;
            };
            if func.execution_count == 0 {
                let why = spec
                    .dead
                    .iter()
                    .find(|d| d.file == fspec.name && &d.subroutine == sub)
                    .map(|d| (d.open, d.why.clone()));
                never_called.push((format!("{}:{sub}", fspec.name), why));
            }
            let lines = by_fn.get(sym.as_str()).cloned().unwrap_or_default();
            let (mut branches, mut taken) = (0usize, 0usize);
            let mut dead = vec![];
            for l in lines {
                if l.line_number < func.start_line || l.line_number > func.end_line {
                    continue;
                }
                for (bi, b) in l.branches.iter().enumerate() {
                    branches += 1;
                    if b.count > 0 {
                        taken += 1;
                        continue;
                    }
                    let text = source.get(l.line_number - 1).cloned().unwrap_or_default();
                    let ann = spec
                        .unreachable
                        .iter()
                        .find(|u| u.file == fspec.name && u.line == l.line_number && u.branch.is_none_or(|k| k == bi));
                    let note = spec
                        .notes
                        .iter()
                        .find(|u| u.file == fspec.name && u.line == l.line_number && u.branch.is_none_or(|k| k == bi));
                    let is_do = text.trim_start().to_uppercase().starts_with("DO ");
                    // a DO line carrying a reach note is an iteration cap (Newton loop exhausted), not a
                    // zero-trip entry: it stays open
                    let kind = if ann.is_some() {
                        Kind::Unreachable
                    } else if is_do && note.is_none() {
                        Kind::LoopEntry
                    } else {
                        Kind::Open
                    };
                    dead.push(Dead {
                        line: l.line_number,
                        branch: bi,
                        fallthrough: b.fallthrough,
                        siblings: l
                            .branches
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| *k != bi)
                            .map(|(_, s)| s.count)
                            .collect(),
                        source: text.trim().to_string(),
                        kind,
                        class: ann.map(|u| u.class.clone()).unwrap_or_default(),
                        why: ann.map(|u| u.why.clone()).or_else(|| note.map(|n| n.reach.clone())),
                    });
                }
            }
            reports.push(SubReport {
                file: fspec.name.clone(),
                name: sub.clone(),
                calls: func.execution_count,
                branches,
                taken,
                dead,
            });
        }
    }
    // annotations that point at nothing dead are stale claims — report them
    let mut stale: Vec<String> = vec![];
    for u in &spec.unreachable {
        let hit = reports.iter().any(|r| {
            r.file == u.file
                && r.dead
                    .iter()
                    .any(|d| d.line == u.line && u.branch.is_none_or(|k| k == d.branch))
        });
        if !hit {
            stale.push(format!(
                "{}:{}{}",
                u.file,
                u.line,
                u.branch.map(|b| format!(" branch {b}")).unwrap_or_default()
            ));
        }
    }
    for n in &spec.notes {
        let hit = reports.iter().any(|r| {
            r.file == n.file
                && r.dead
                    .iter()
                    .any(|d| d.kind == Kind::Open && d.line == n.line && n.branch.is_none_or(|k| k == d.branch))
        });
        if !hit {
            stale.push(format!("note {}:{}", n.file, n.line));
        }
    }
    for d in &spec.dead {
        if !never_called
            .iter()
            .any(|(n, _)| *n == format!("{}:{}", d.file, d.subroutine))
        {
            stale.push(format!("dead {}:{}", d.file, d.subroutine));
        }
    }

    let md = render(
        &reports,
        &ran,
        &never_called,
        &missing,
        &stale,
        &gcc_version,
        &gcov,
        big,
    );
    let out = root.join("docs/validation/coverage.md");
    fs::write(&out, md).unwrap();
    let open: usize = reports
        .iter()
        .map(|r| r.dead.iter().filter(|d| d.kind == Kind::Open).count())
        .sum();
    let total: usize = reports.iter().map(|r| r.branches).sum();
    let taken: usize = reports.iter().map(|r| r.taken).sum();
    println!(
        "coverage: {} subroutines, {taken}/{total} branches taken, {open} open, {} never called ({} unexplained), {} not found, {} stale annotations -> {}",
        reports.len(),
        never_called.len(),
        never_called.iter().filter(|(_, w)| w.as_ref().is_none_or(|(o, _)| *o)).count(),
        missing.len(),
        stale.len(),
        out.display()
    );
    if failures > 0 {
        eprintln!("{failures} case(s) failed to run");
        std::process::exit(1);
    }
}

/// gcov must match the gfortran that produced the .gcno files (Homebrew ships it as gcov-<major>).
fn gcov_binary() -> String {
    if let Ok(g) = std::env::var("GCOV") {
        return g;
    }
    let major = Command::new("gfortran")
        .arg("-dumpversion")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().split('.').next().unwrap_or("").to_string())
        .unwrap_or_default();
    let cand = format!("gcov-{major}");
    if Command::new(&cand)
        .arg("--version")
        .stdout(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
    {
        cand
    } else {
        "gcov".to_string()
    }
}

#[allow(clippy::too_many_arguments)]
fn render(
    reports: &[SubReport],
    ran: &[String],
    never_called: &[(String, Option<(bool, String)>)],
    missing: &[String],
    stale: &[String],
    gcc_version: &str,
    gcov: &str,
    big: bool,
) -> String {
    let date = Command::new("date")
        .arg("+%Y-%m-%d")
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    let count = |k: Kind| -> usize {
        reports
            .iter()
            .map(|r| r.dead.iter().filter(|d| d.kind == k).count())
            .sum()
    };
    let total: usize = reports.iter().map(|r| r.branches).sum();
    let taken: usize = reports.iter().map(|r| r.taken).sum();
    let (open, loops, unreach) = (count(Kind::Open), count(Kind::LoopEntry), count(Kind::Unreachable));
    let mut s = String::new();
    s += "# Branch coverage of the translated XFOIL subroutines\n\n";
    s += &format!(
        "Generated by `cargo xtask coverage{}` on {date}. Do not edit; regenerate.\n\n",
        if big { " --big" } else { "" }
    );
    s += "**Instrument.** The pristine XFOIL 6.99 double-precision reference (build patch series only, no \
          instrumentation dumps) compiled and linked with `-fprofile-arcs -ftest-coverage` \
          (`scripts/xfoil-build.sh --gcov`, `target/xfoil-ref/gcov/`), run on every case below from its \
          tracked `xfoil.inp`/`panels.dat`, counters read back with gcov's JSON output. A *branch* is one \
          outgoing edge of a conditional as the compiler emitted it: an `IF` contributes two edges, a \
          `.AND.`/`.OR.` chain one pair per operand, a `DO` loop its zero-trip entry pair, and a loop whose \
          back-edge test is attributed to its last statement adds a pair there. Which arm an edge is \
          (fallthrough vs jump) is the compiler's choice and is reported as such, not as true/false; the \
          *other edges* column gives the counts of the remaining edges on the same line, so the dead arm \
          can be read off against the source.\n\n";
    s += &format!(
        "**Toolchain.** gfortran/gcc {gcc_version}, `{gcov}`. **Subroutine set.** `xtask/fixtures-config/coverage.toml` \
         ({} subroutines in {} files).\n\n",
        reports.len(),
        reports
            .iter()
            .map(|r| r.file.as_str())
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    );
    s += &format!(
        "**Cases run ({}).** {}\n\n",
        ran.len(),
        ran.iter().map(|c| format!("`{c}`")).collect::<Vec<_>>().join(", ")
    );
    s += "## Completeness (CLAUDE.md Rule 6)\n\n";
    s += "The fixture set is complete when every reachable branch in the translated subroutines has been \
          taken at least once. Never-taken branches fall into three classes:\n\n\
          - **open** — never taken, no annotation. These are the objective for further cases and for the \
          parameter/fuzz harness. \"Numerically exact\" is not yet claimable over them.\n\
          - **loop entry** — the zero-trip edge of a `DO` loop whose body always ran (or a loop inside an \
          unreached block). Listed separately because they are structurally impossible on the analysis \
          path (`DO IBL=2,NBL(IS)` with `NBL ≥ 2`). A `DO` loop that is an iteration cap (a Newton loop \
          that could be exhausted) carries a reach note and stays in the open list.\n\
          - **unreachable** — annotated in `xtask/fixtures-config/coverage.toml` with a class and the reason: \
          *structural* (impossible by construction on the analysis path), *mode* (a feature outside it: \
          inverse design, image airfoil, flap hinge, interactive prompts), *guard* (array-bound or \
          illegal-input STOP; YFoil has no fixed dimensions), *compiler* (no source-level decision).\n\n\
          Open branches with a known way to reach them carry a *reach* note; those are the next cases to \
          add, and the rest are what the parameter/fuzz harness is for.\n\n";
    s += "| | count |\n|---|---|\n";
    s += &format!("| branches in the translated set | {total} |\n| taken | {taken} |\n| never taken — open | **{open}** |\n| never taken — loop entry | {loops} |\n| never taken — annotated unreachable | {unreach} |\n");
    for class in ["structural", "mode", "guard", "compiler"] {
        let n: usize = reports
            .iter()
            .map(|r| {
                r.dead
                    .iter()
                    .filter(|d| d.kind == Kind::Unreachable && d.class == class)
                    .count()
            })
            .sum();
        s += &format!("| &nbsp;&nbsp;&nbsp;&nbsp;{class} | {n} |\n");
    }
    let unexplained = never_called
        .iter()
        .filter(|(_, w)| w.as_ref().is_none_or(|(o, _)| *o))
        .count();
    s += &format!(
        "| subroutines never called | {} ({unexplained} unexplained) |\n\n",
        never_called.len()
    );
    s += &format!(
        "**Verdict: {}**\n\n",
        if open == 0 && unexplained == 0 {
            "complete — every reachable branch of the translated set has been exercised."
        } else {
            "incomplete — the open branches below are not yet covered by any case."
        }
    );
    if !never_called.is_empty() {
        s += "Never called:\n\n";
        for (n, w) in never_called {
            s += &format!(
                "- `{n}` — {}\n",
                match w {
                    None => "**unexplained**".to_string(),
                    Some((true, why)) => format!("**open** — {why}"),
                    Some((false, why)) => why.clone(),
                }
            );
        }
        s += "\n";
    }
    if !missing.is_empty() {
        s += &format!(
            "Not found in the gcov output (check the name in `xtask/fixtures-config/coverage.toml`): {}.\n\n",
            missing.iter().map(|n| format!("`{n}`")).collect::<Vec<_>>().join(", ")
        );
    }
    if !stale.is_empty() {
        s += &format!(
            "Stale annotations (the branch is now taken or does not exist): {}.\n\n",
            stale.iter().map(|n| format!("`{n}`")).collect::<Vec<_>>().join(", ")
        );
    }

    s += "## Per subroutine\n\n| file | subroutine | calls | branches | taken | open | loop entry | unreachable |\n|---|---|---:|---:|---:|---:|---:|---:|\n";
    for r in reports {
        let k = |kind: Kind| r.dead.iter().filter(|d| d.kind == kind).count();
        let open = k(Kind::Open);
        s += &format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} |\n",
            r.file,
            r.name,
            r.calls,
            r.branches,
            r.taken,
            if open > 0 { format!("**{open}**") } else { "0".into() },
            k(Kind::LoopEntry),
            k(Kind::Unreachable)
        );
    }
    s += "\n";

    let section = |s: &mut String, title: &str, kind: Kind, last: Option<&str>| {
        s.push_str(&format!("## {title}\n\n"));
        let mut any = false;
        for r in reports {
            let rows: Vec<&Dead> = r.dead.iter().filter(|d| d.kind == kind).collect();
            if rows.is_empty() {
                continue;
            }
            any = true;
            s.push_str(&format!("### `{}` — {}\n\n", r.file, r.name));
            match last {
                Some(col) => s.push_str(&format!(
                    "| line | edge | other edges | source | {col} |\n|---:|---|---|---|---|\n"
                )),
                None => s.push_str("| line | edge | other edges | source |\n|---:|---|---|---|\n"),
            }
            for d in rows {
                let edge = format!("{} ({})", d.branch, if d.fallthrough { "fallthrough" } else { "jump" });
                let others = d.siblings.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", ");
                let src = d.source.replace('|', "\\|");
                match last {
                    Some(_) => {
                        let mut why = d.why.clone().unwrap_or_default();
                        if kind == Kind::Unreachable {
                            why = format!("*{}* — {why}", d.class);
                        }
                        s.push_str(&format!("| {} | {edge} | {others} | `{src}` | {why} |\n", d.line));
                    }
                    None => s.push_str(&format!("| {} | {edge} | {others} | `{src}` |\n", d.line)),
                }
            }
            s.push('\n');
        }
        if !any {
            s.push_str("None.\n\n");
        }
    };
    section(
        &mut s,
        "Open branches (never taken, not annotated)",
        Kind::Open,
        Some("reach"),
    );
    section(&mut s, "Annotated unreachable branches", Kind::Unreachable, Some("why"));
    section(
        &mut s,
        "DO-loop zero-trip entry edges never taken",
        Kind::LoopEntry,
        None,
    );
    s
}
