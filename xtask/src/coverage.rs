//! `cargo xtask coverage [--case NAME]... [--group NAME] [--with-group NAME] [--big] [--rebuild] [--per-case]`
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
    /// the test that gates the branch as an event when it is reached only in a non-finite
    /// reference run (`tests/xfoil_nonfinite_events_tests.rs::<test>`)
    #[serde(default)]
    gate: Option<String>,
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
    count: u64,
    #[serde(default)]
    function_name: Option<String>,
    #[serde(default)]
    branches: Vec<GcovBranch>,
    #[serde(default)]
    calls: Vec<GcovCall>,
}
#[derive(Deserialize)]
struct GcovCall {
    returned: u64,
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
    /// a never-taken edge of the loop-control test gfortran placed on a DO line (the arm that leaves
    /// the loop by its count), with no reach note
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
    /// every branch of the subroutine, taken or not (the `--per-case` universe)
    all: Vec<BranchRecord>,
}

/// One branch of the translated set as `target/coverage/branches.json` records it.
struct BranchRecord {
    line: usize,
    branch: usize,
    fallthrough: bool,
    source: String,
    status: &'static str,
    class: String,
    why: Option<String>,
    gate: Option<String>,
}

/// gcov counters of every file of the translated set, summable across cases.
struct Counters {
    gcc_version: String,
    files: HashMap<String, GcovFile>,
    sources: HashMap<String, Vec<String>>,
}

impl Counters {
    /// Add another measurement of the same binary (identical structure: same .gcno files).
    fn add(&mut self, other: &Counters) {
        for (name, gf) in &other.files {
            let mine = self.files.get_mut(name).expect("same file set");
            for (a, b) in mine.functions.iter_mut().zip(&gf.functions) {
                assert_eq!(a.name, b.name, "gcov structure differs between cases");
                a.execution_count += b.execution_count;
            }
            for (a, b) in mine.lines.iter_mut().zip(&gf.lines) {
                assert_eq!(a.line_number, b.line_number, "gcov structure differs between cases");
                for (x, y) in a.branches.iter_mut().zip(&b.branches) {
                    x.count += y.count;
                }
            }
        }
    }
}

fn clear_counters(bin: &Path) {
    for e in fs::read_dir(bin).unwrap().flatten() {
        if e.path().extension().is_some_and(|x| x == "gcda") {
            fs::remove_file(e.path()).unwrap();
        }
    }
}

/// Read the accumulated counters of every file in the translated set back with gcov.
fn read_counters(gcov: &str, gdir: &Path, bin: &Path, spec: &Spec, scratch: &Path) -> Counters {
    let _ = fs::remove_dir_all(scratch);
    fs::create_dir_all(scratch).unwrap();
    let mut out = Counters {
        gcc_version: String::new(),
        files: HashMap::new(),
        sources: HashMap::new(),
    };
    for fspec in &spec.files {
        let src_path = gdir.join("src").join(&fspec.name);
        let source: Vec<String> = fs::read_to_string(&src_path)
            .unwrap_or_else(|_| panic!("read {}", src_path.display()))
            .lines()
            .map(str::to_string)
            .collect();
        let st = Command::new(gcov)
            .args(["-j", "-b", "-o"])
            .arg(bin)
            .arg(&src_path)
            .current_dir(scratch)
            .stdout(Stdio::null())
            .status()
            .unwrap_or_else(|e| panic!("run {gcov}: {e}"));
        assert!(st.success(), "{gcov} failed on {}", fspec.name);
        // gcov names its output after the object stem (xpanel.gcda -> xpanel.gcov.json.gz)
        let stem = fspec.name.trim_end_matches(".f");
        let gz = scratch.join(format!("{stem}.gcov.json.gz"));
        let json = Command::new("gzip").args(["-dc"]).arg(&gz).output().expect("gzip -dc");
        assert!(json.status.success(), "gzip -dc {} failed", gz.display());
        let mut parsed: GcovJson = serde_json::from_slice(&json.stdout).expect("parse gcov json");
        out.gcc_version = parsed.gcc_version.clone();
        let idx = parsed
            .files
            .iter()
            .position(|f| Path::new(&f.file).file_name() == src_path.file_name())
            .unwrap_or_else(|| panic!("no entry for {} in gcov json", fspec.name));
        out.files.insert(fspec.name.clone(), parsed.files.swap_remove(idx));
        out.sources.insert(fspec.name.clone(), source);
    }
    out
}

/// `target/coverage/per-case/<case>/branches.json`: the branches of the translated set this one case took,
/// as `file:line:branch` ids, with the per-subroutine call counts.
fn write_case_branches(work: &Path, case: &str, spec: &Spec, c: &Counters) {
    let mut taken: Vec<String> = vec![];
    let mut calls = serde_json::Map::new();
    for fspec in &spec.files {
        let gf = &c.files[&fspec.name];
        for sub in &fspec.subroutines {
            let sym = format!("{}_", sub.to_lowercase());
            let Some(func) = gf.functions.iter().find(|f| f.name == sym) else {
                continue;
            };
            calls.insert(format!("{}:{sub}", fspec.name), serde_json::json!(func.execution_count));
            for l in &gf.lines {
                if l.function_name.as_deref() != Some(sym.as_str())
                    || l.line_number < func.start_line
                    || l.line_number > func.end_line
                {
                    continue;
                }
                for (bi, b) in l.branches.iter().enumerate() {
                    if b.count > 0 {
                        taken.push(format!("{}:{}:{bi}", fspec.name, l.line_number));
                    }
                }
            }
        }
    }
    // every DO line of the set: the counts of its loop-control edges and of the first executable
    // line of its body (`body` > 0 = the body ran on this case)
    let mut do_loops = serde_json::Map::new();
    for fspec in &spec.files {
        let gf = &c.files[&fspec.name];
        let src = &c.sources[&fspec.name];
        for (i, l) in gf.lines.iter().enumerate() {
            let text = src
                .get(l.line_number - 1)
                .map(|t| t.trim_start().to_uppercase())
                .unwrap_or_default();
            if !text.starts_with("DO ") || l.branches.is_empty() {
                continue;
            }
            let Some(sub) = l.function_name.as_deref() else {
                continue;
            };
            if !fspec
                .subroutines
                .iter()
                .any(|s| format!("{}_", s.to_lowercase()) == sub)
            {
                continue;
            }
            let body = gf.lines[i + 1..]
                .iter()
                .find(|n| n.function_name.as_deref() == Some(sub) && n.line_number > l.line_number)
                .map(|n| n.count)
                .unwrap_or(0);
            do_loops.insert(
                format!("{}:{}", fspec.name, l.line_number),
                serde_json::json!({ "edges": l.branches.iter().map(|b| b.count).collect::<Vec<_>>(), "body": body }),
            );
        }
    }
    let j = serde_json::json!({ "case": case, "subroutine_calls": calls, "taken": taken, "do_loops": do_loops });
    fs::write(work.join("branches.json"), serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

/// `target/coverage/branches.json`: every branch of the translated set with its status over the
/// cases run (taken / open / loop_entry / unreachable), the universe of the branch-coverage study.
fn write_universe(path: &Path, reports: &[SubReport], ran: &[String]) {
    let branches: Vec<serde_json::Value> = reports
        .iter()
        .flat_map(|r| {
            r.all.iter().map(move |b| {
                serde_json::json!({
                    "id": format!("{}:{}:{}", r.file, b.line, b.branch),
                    "file": r.file, "subroutine": r.name, "line": b.line, "branch": b.branch,
                    "fallthrough": b.fallthrough, "source": b.source, "status": b.status,
                    "class": b.class, "why": b.why, "gate": b.gate,
                })
            })
        })
        .collect();
    let subroutines: Vec<serde_json::Value> = reports
        .iter()
        .map(|r| serde_json::json!({ "file": r.file, "name": r.name, "calls": r.calls, "branches": r.branches, "taken": r.taken }))
        .collect();
    let j = serde_json::json!({ "cases": ran, "subroutines": subroutines, "branches": branches });
    fs::write(path, serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

pub(crate) fn run(flags: &[String]) {
    let root = super::root();
    let big = flags.iter().any(|f| f == "--big");
    let rebuild = flags.iter().any(|f| f == "--rebuild");
    let per_case = flags.iter().any(|f| f == "--per-case");
    let group: Option<&str> = flags.windows(2).find(|w| w[0] == "--group").map(|w| w[1].as_str());
    // `--with-group NAME`: the default selection plus that group's (untracked) members
    let with_group: Option<&str> = flags.windows(2).find(|w| w[0] == "--with-group").map(|w| w[1].as_str());
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
    clear_counters(&bin);

    // run the cases from their tracked inputs
    let cases: super::Cases = toml::from_str(
        &fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("xtask/fixtures-config/cases.toml"),
    )
    .expect("parse cases.toml");
    let spec: Spec = toml::from_str(
        &fs::read_to_string(root.join("xtask/fixtures-config/coverage.toml"))
            .expect("xtask/fixtures-config/coverage.toml"),
    )
    .expect("parse coverage.toml");
    let gcov = gcov_binary();
    let scratch = root.join("target/coverage/_gcov");
    let mut ran: Vec<String> = vec![];
    let mut failures = 0;
    // `--per-case`: the counters are read back after every case (and reset before it), so each
    // case's own branch set is known; the aggregate is their sum
    let mut acc: Option<Counters> = None;
    for case in &cases.cases {
        let added = with_group.is_some_and(|g| case.group.as_deref() == Some(g));
        if !case.selected(&selected, group, big) && !added {
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
        if per_case {
            clear_counters(&bin);
        }
        // under the fixture runner's watchdog: XFOIL's plot-label hang ends a run (kept, truncated),
        // a run still going after an hour fails
        // under the fixture runner's watchdog; a run it ends keeps its counts (the gcov build
        // writes them on SIGTERM), so a hung or long run is counted up to where it was ended
        match super::run_xfoil(&bin.join("xfoil"), &work) {
            Ok(super::RunEnd::Finished) => {}
            Ok(super::RunEnd::Hung) => println!("  {}: hung in the plot label; counted up to there", case.name),
            Ok(super::RunEnd::TimedOut) => {
                println!("  {}: still running after an hour; counted up to there", case.name)
            }
            Err(e) => {
                eprintln!("  {}: {e}", case.name);
                failures += 1;
                continue;
            }
        }
        println!("  ran {}", case.name);
        ran.push(case.name.clone());
        if per_case {
            // kept apart from the work dirs, which every plain run recreates
            let one = read_counters(&gcov, &gdir, &bin, &spec, &scratch);
            let per = root.join("target/coverage/per-case").join(&case.name);
            let _ = fs::remove_dir_all(&per);
            fs::create_dir_all(&per).unwrap();
            fs::copy(work.join("stdout.txt"), per.join("stdout.txt")).unwrap();
            write_case_branches(&per, &case.name, &spec, &one);
            acc = Some(match acc.take() {
                None => one,
                Some(mut a) => {
                    a.add(&one);
                    a
                }
            });
        }
    }
    assert!(!ran.is_empty(), "no cases ran");
    let counters = match acc {
        Some(a) => a,
        None => read_counters(&gcov, &gdir, &bin, &spec, &scratch),
    };
    let gcc_version = counters.gcc_version.clone();

    // classify every branch of the translated set
    let mut reports: Vec<SubReport> = vec![];
    let mut never_called: Vec<(String, Option<(bool, String)>)> = vec![];
    let mut missing: Vec<String> = vec![];
    for fspec in &spec.files {
        let source = &counters.sources[&fspec.name];
        let gf = &counters.files[&fspec.name];
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
            let mut all = vec![];
            for l in lines {
                if l.line_number < func.start_line || l.line_number > func.end_line {
                    continue;
                }
                for (bi, b) in l.branches.iter().enumerate() {
                    branches += 1;
                    let text = source.get(l.line_number - 1).cloned().unwrap_or_default();
                    if b.count > 0 {
                        taken += 1;
                        // a taken branch keeps its reach note, if any: the branch-coverage study
                        // reports branches taken only in non-finite runs with it
                        let note = spec.notes.iter().find(|u| {
                            u.file == fspec.name && u.line == l.line_number && u.branch.is_none_or(|k| k == bi)
                        });
                        all.push(BranchRecord {
                            line: l.line_number,
                            branch: bi,
                            fallthrough: b.fallthrough,
                            source: text.trim().to_string(),
                            status: "taken",
                            class: String::new(),
                            why: note.map(|n| n.reach.clone()),
                            gate: note.and_then(|n| n.gate.clone()),
                        });
                        continue;
                    }
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
                    // loop-control edge: it stays open
                    let kind = if ann.is_some() {
                        Kind::Unreachable
                    } else if is_do && note.is_none() {
                        Kind::LoopEntry
                    } else {
                        Kind::Open
                    };
                    let why = ann.map(|u| u.why.clone()).or_else(|| note.map(|n| n.reach.clone()));
                    all.push(BranchRecord {
                        line: l.line_number,
                        branch: bi,
                        fallthrough: b.fallthrough,
                        source: text.trim().to_string(),
                        status: match kind {
                            Kind::Unreachable => "unreachable",
                            Kind::LoopEntry => "loop_entry",
                            Kind::Open => "open",
                        },
                        class: ann.map(|u| u.class.clone()).unwrap_or_default(),
                        why: why.clone(),
                        gate: note.and_then(|n| n.gate.clone()),
                    });
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
                        why,
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
                all,
            });
        }
    }
    // the universe (read by `cargo xtask steps`' report and the branch-coverage study)
    write_universe(&root.join("target/coverage/branches.json"), &reports, &ran);
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
          `.AND.`/`.OR.` chain one pair per operand, a `DO` line the pair of its loop-control test (iterate again, or \
          leave the loop by its count; measured, the second arm is the normal completion of a top-tested loop \
          and a never-taken duplicate for a loop gfortran rotated), and a loop whose \
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
          taken at least once. Every count is the reference's (gcov on pristine XFOIL). yFoil's side of \
          the claim is measured separately: `cargo xtask route` observes yFoil take the reference's route through \
          each step that gates a branch ([branch-gating.md](branch-gating.md)). Never-taken branches fall into three classes:\n\n\
          - **open** — never taken on any case run here, and not annotated: a statement about the case set, \
          not about the code — no case has been constructed that reaches the edge. These are the objective \
          for further cases and for the parameter/fuzz harness. \"Numerically exact\" is not yet claimable \
          over them.\n\
          - **DO-line edge** — a never-taken arm of the loop-control test on a `DO` line, for a loop whose body \
          ran (its *iterate* arm is counted as taken like any other edge): the arm that leaves the loop by its \
          count was never taken because the loop is always left early by a `GO TO`/`RETURN`, or because \
          gfortran rotated the loop and placed the real exit on its last statement; or both arms of a loop \
          inside a block never reached. None of them is a body that ran zero times: the analysis path's \
          station, node and wake counts never make a loop empty. A `DO` line that is an iteration cap (a \
          Newton loop that could be exhausted) carries a reach note and both of its edges stay in the open \
          list until taken.\n\
          - **unreachable** — a recorded claim, annotated in `xtask/fixtures-config/coverage.toml` with a \
          class and the reason, that no execution of XFOIL's analysis path can take the edge: a code-reading \
          argument about the reference's own source, never \"no case was found\" and never about yFoil, and \
          falsifiable — a case that takes an annotated edge is reported below as a stale annotation. Classes: \
          *structural* (impossible by construction on the analysis path), *mode* (a feature outside it: \
          inverse design, image airfoil, flap hinge, interactive prompts), *guard* (array-bound or \
          illegal-input STOP; yFoil has no fixed dimensions), *compiler* (no source-level decision), \
          *numerical* (reachable only by a floating-point coincidence — a real landing bitwise on a \
          threshold — or by the iteration cap of a Newton iteration that converges on the analysis path; \
          no input can be designed to reach it, and the parameter/fuzz harness is the instrument).\n\n\
          Open branches with a known way to reach them carry a *reach* note; those are the next cases to \
          add, and the rest are what the parameter/fuzz harness is for.\n\n";
    s += "| | count |\n|---|---|\n";
    s += &format!("| branches in the translated set | {total} |\n| taken | {taken} |\n| never taken — open | **{open}** |\n| never taken — DO-line edge | {loops} |\n| never taken — annotated unreachable | {unreach} |\n");
    for class in ["structural", "mode", "guard", "compiler", "numerical"] {
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

    s += "## Per subroutine\n\n| file | subroutine | calls | branches | taken | open | DO-line | unreachable |\n|---|---|---:|---:|---:|---:|---:|---:|\n";
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
    section(&mut s, "DO-line loop-control edges never taken", Kind::LoopEntry, None);
    s
}

/// The gcov build's directory and binary, building it if needed.
pub(crate) fn gcov_build(root: &Path, rebuild: bool) -> std::path::PathBuf {
    let gdir = root.join("target/xfoil-ref/gcov");
    if rebuild || !gdir.join("bin/xfoil").exists() {
        super::run(
            Command::new(root.join("scripts/xfoil-build.sh")).arg("--gcov"),
            "scripts/xfoil-build.sh --gcov",
        );
    }
    gdir
}

/// Delete the gcov build's counters, so the next run measures only itself.
pub(crate) fn reset_counters(gdir: &Path) {
    for e in fs::read_dir(gdir.join("bin")).unwrap().flatten() {
        if e.path().extension().is_some_and(|x| x == "gcda") {
            fs::remove_file(e.path()).unwrap();
        }
    }
}

/// The accumulated counts of the translated subroutines (`coverage.toml`).
pub(crate) struct Counts {
    /// every branch, keyed `file:line:index` — the identifiers the coverage report and the step
    /// cover use
    pub branches: HashMap<String, u64>,
    /// every subroutine's call count, keyed by its name
    pub calls: HashMap<String, u64>,
    /// the calls made from every line that makes one, keyed `file:line`
    pub sites: HashMap<String, u64>,
}

pub(crate) fn branch_counts(root: &Path, gdir: &Path, scratch: &Path) -> Counts {
    let spec: Spec = toml::from_str(
        &fs::read_to_string(root.join("xtask/fixtures-config/coverage.toml"))
            .expect("xtask/fixtures-config/coverage.toml"),
    )
    .expect("parse coverage.toml");
    let gcov = gcov_binary();
    let bin = gdir.join("bin");
    let _ = fs::remove_dir_all(scratch);
    fs::create_dir_all(scratch).unwrap();
    let mut out = HashMap::new();
    let mut calls = HashMap::new();
    let mut sites = HashMap::new();
    for fspec in &spec.files {
        let src_path = gdir.join("src").join(&fspec.name);
        let st = Command::new(&gcov)
            .args(["-j", "-b", "-o"])
            .arg(&bin)
            .arg(&src_path)
            .current_dir(scratch)
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .unwrap_or_else(|e| panic!("run {gcov}: {e}"));
        assert!(st.success(), "{gcov} failed on {}", fspec.name);
        let stem = fspec.name.trim_end_matches(".f");
        let gz = scratch.join(format!("{stem}.gcov.json.gz"));
        let json = Command::new("gzip").args(["-dc"]).arg(&gz).output().expect("gzip -dc");
        assert!(json.status.success(), "gzip -dc {} failed", gz.display());
        let parsed: GcovJson = serde_json::from_slice(&json.stdout).expect("parse gcov json");
        let gf = parsed
            .files
            .iter()
            .find(|f| Path::new(&f.file).file_name() == src_path.file_name())
            .unwrap_or_else(|| panic!("no entry for {} in gcov json", fspec.name));
        for sub in &fspec.subroutines {
            let sym = format!("{}_", sub.to_lowercase());
            let Some(func) = gf.functions.iter().find(|f| f.name == sym) else {
                continue;
            };
            calls.insert(sub.clone(), func.execution_count);
            for l in &gf.lines {
                if l.function_name.as_deref() != Some(sym.as_str())
                    || l.line_number < func.start_line
                    || l.line_number > func.end_line
                {
                    continue;
                }
                for (bi, b) in l.branches.iter().enumerate() {
                    out.insert(format!("{}:{}:{bi}", fspec.name, l.line_number), b.count);
                }
                if !l.calls.is_empty() {
                    sites.insert(
                        format!("{}:{}", fspec.name, l.line_number),
                        l.calls.iter().map(|c| c.returned).sum(),
                    );
                }
            }
        }
    }
    Counts {
        branches: out,
        calls,
        sites,
    }
}
