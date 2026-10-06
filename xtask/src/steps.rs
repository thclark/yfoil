//! The step cover: which single steps of which reference runs the (b) branch tests replay.
//!
//! `cargo xtask twins --case NAME...` reruns a case's reference (`cargo xtask fixtures` must have
//! made its work directory, `target/fixtures/<case>/`) on five seeded 1-ULP twins of its panels,
//! `ulp1` … `ulp5` beside it. A twin jogs every panel coordinate, x and y independently, by −1, 0
//! or +1 ULP drawn from splitmix64 of (seed, node), so it is reproducible from its seed
//! (`docs/conventions/terminology.md`, *twin*).
//!
//! `cargo xtask steps [--case NAME]... [--rebuild]` then, for every candidate case:
//!
//! 1. **attributes branches to steps.** A step is one VISCAL iteration (iteration `i` of call
//!    `c`, global SETBL call `k`) or, for an inviscid case, one SPECAL/SPECCL call. The gcov build
//!    of the reference reruns the case's script truncated after step (c, i) — the first `c`
//!    operating points, the last with `ITER i` — and the branches whose counts rose over the
//!    previous step's prefix are the ones that step takes. A step's prefix also runs the
//!    end-of-call code once, which cancels between consecutive prefixes except where the
//!    convergence decision differs: that is attributed to the step that made it.
//! 2. **decides whether each step is well-conditioned** from the twins: every twin reproduces the
//!    reference's iteration counts, IST, ISTB and ITRAN through the step, and the twins' spread of
//!    RMSBL, RLX, CL, CD, CM and ALFA at the step is at most a quarter of `TOL_SOLVER`.
//! 3. **chooses the cover**: the fewest well-conditioned steps that take every branch any
//!    candidate step takes (essential steps first, then greedily, then pruned), and lists the
//!    branches that only ill-conditioned steps take — those get no test
//!    (`docs/conventions/testing.md`, rule 3).
//!
//! Output: `target/coverage/steps.json` (every step: its branches and its conditioning) and
//! `xtask/fixtures-config/step-cover.toml` (the chosen steps, tracked: it is the input of the
//! fixture cases and the list of tests in `tests/execution/branches.rs`).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

/// The named tolerance a step's twin spread is measured against (`TOL_SOLVER`,
/// `tests/common/utilities/tolerances.rs`), and the fraction of it the spread may use.
const TOL_SOLVER: f64 = 1e-10;
const CONDITIONING_FRACTION: f64 = 0.25;
const SEEDS: [u64; 5] = [1, 2, 3, 4, 5];

pub fn twins(flags: &[String]) {
    let root = super::root();
    let xfoil = root.join("target/xfoil-ref/instrumented/bin/xfoil");
    assert!(
        xfoil.exists(),
        "no instrumented reference: run `cargo xtask xfoil-build`"
    );
    let names: Vec<&str> = flags
        .windows(2)
        .filter(|w| w[0] == "--case")
        .map(|w| w[1].as_str())
        .collect();
    assert!(!names.is_empty(), "cargo xtask twins --case NAME...");
    let mut failures = 0;
    for name in names {
        let work = root.join("target/fixtures").join(name);
        if !work.join("xfoil.inp").exists() {
            eprintln!("  {name}: no work directory — run `cargo xtask fixtures --case {name}` first");
            failures += 1;
            continue;
        }
        let results: Vec<Result<(), String>> = std::thread::scope(|scope| {
            let handles: Vec<_> = SEEDS
                .iter()
                .map(|&seed| {
                    let (work, xfoil) = (&work, &xfoil);
                    scope.spawn(move || twin(xfoil, work, seed))
                })
                .collect();
            handles.into_iter().map(|h| h.join().expect("twin thread")).collect()
        });
        let errs: Vec<String> = results.into_iter().filter_map(Result::err).collect();
        if errs.is_empty() {
            println!("  {name}: {} twins", SEEDS.len());
        } else {
            eprintln!("  {name}: {}", errs.join("; "));
            failures += 1;
        }
    }
    if failures > 0 {
        std::process::exit(1);
    }
}

/// The jog of node `i`'s coordinates in twin `seed`: −1, 0 or +1 ULP for x and for y.
fn jog(seed: u64, i: usize) -> (i8, i8) {
    let mut z = (seed << 32) ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let three = |v: u64| -> i8 { (v % 3) as i8 - 1 };
    (three(z & 0xFFFF_FFFF), three(z >> 32))
}

fn twin(xfoil: &Path, work: &Path, seed: u64) -> Result<(), String> {
    let dir = work.join(format!("ulp{seed}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let src = fs::read_to_string(work.join("panels.dat")).map_err(|e| e.to_string())?;
    let mut lines = src.lines();
    let mut out = format!("{}\n", lines.next().unwrap_or(""));
    let step = |v: f64, by: i8| match by {
        1 => v.next_up(),
        -1 => v.next_down(),
        _ => v,
    };
    for (i, l) in lines.filter(|l| l.split_whitespace().count() >= 2).enumerate() {
        let t: Vec<&str> = l.split_whitespace().collect();
        let x: f64 = t[0].parse().map_err(|_| format!("bad coordinate {l}"))?;
        let y: f64 = t[1].parse().map_err(|_| format!("bad coordinate {l}"))?;
        let (jx, jy) = jog(seed, i);
        out.push_str(&format!(" {:.17e}  {:.17e}\n", step(x, jx), step(y, jy)));
    }
    fs::write(dir.join("panels.dat"), out).map_err(|e| e.to_string())?;
    for f in ["xfoil.inp", "dump_calls.txt"] {
        if work.join(f).exists() {
            fs::copy(work.join(f), dir.join(f)).map_err(|e| e.to_string())?;
        }
    }
    let st = Command::new(xfoil)
        .current_dir(&dir)
        .stdin(fs::File::open(dir.join("xfoil.inp")).map_err(|e| e.to_string())?)
        .stdout(fs::File::create(dir.join("stdout.txt")).map_err(|e| e.to_string())?)
        .stderr(Stdio::null())
        .status()
        .map_err(|e| e.to_string())?;
    if st.success() {
        Ok(())
    } else {
        Err(format!("twin {seed}: xfoil exited {st}"))
    }
}

/// One step of a reference run.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
struct Step {
    case: String,
    /// VISCAL call (viscous) or SPECAL/SPECCL call (inviscid), 1-based
    call: usize,
    /// iteration within the call (0 for an inviscid call)
    iteration: usize,
    /// global SETBL call number (0 for an inviscid call)
    setbl: usize,
    /// NITDONE of the step's VISCAL call (0 for an inviscid call)
    call_iterations: usize,
    taken: Vec<String>,
    /// how many times the step called each translated subroutine (`cargo xtask route` compares
    /// yFoil's calls with these)
    #[serde(default)]
    calls: BTreeMap<String, i64>,
    /// how many calls the step made from each call site `route-map.toml` names (`file:line`)
    #[serde(default)]
    sites: BTreeMap<String, i64>,
    conditioned: Option<bool>,
    /// why the step is not well-conditioned, or the worst twin spread relative to the tolerance
    conditioning: String,
}

pub fn steps(flags: &[String]) {
    let root = super::root();
    let rebuild = flags.iter().any(|f| f == "--rebuild");
    let selected: Vec<&str> = flags
        .windows(2)
        .filter(|w| w[0] == "--case")
        .map(|w| w[1].as_str())
        .collect();
    let cases: super::Cases =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("parse cases.toml");
    if flags.iter().any(|f| f == "--from-json") {
        let all: Vec<Step> = serde_json::from_str(
            &fs::read_to_string(root.join("target/coverage/steps.json")).expect("target/coverage/steps.json"),
        )
        .expect("parse steps.json");
        let cover = choose_cover(&all);
        fs::write(root.join("xtask/fixtures-config/step-cover.toml"), &cover).unwrap();
        write_report(&root, &all, &cover);
        println!("cover recomputed from target/coverage/steps.json");
        return;
    }
    let gdir = super::coverage::gcov_build(&root, rebuild);
    let mut all: Vec<Step> = vec![];
    for case in &cases.cases {
        let candidate = if selected.is_empty() {
            case.track
                && !case.geometry_only
                && case.tgap.is_empty()
                && (!case.alphas.is_empty() || !case.cls.is_empty())
        } else {
            selected.contains(&case.name.as_str())
        };
        if !candidate {
            continue;
        }
        match case_steps(&root, &gdir, case) {
            Ok(s) => {
                println!(
                    "  {}: {} steps, {} well-conditioned",
                    case.name,
                    s.len(),
                    s.iter().filter(|x| x.conditioned == Some(true)).count()
                );
                all.extend(s);
            }
            Err(e) => eprintln!("  {}: {e}", case.name),
        }
    }
    // with --case, merge with the attribution already measured for the cases not rerun now
    let out = root.join("target/coverage/steps.json");
    fs::create_dir_all(out.parent().unwrap()).unwrap();
    if let (false, Ok(text)) = (selected.is_empty(), fs::read_to_string(&out)) {
        if let Ok(previous) = serde_json::from_str::<Vec<Step>>(&text) {
            let rerun: BTreeSet<String> = all.iter().map(|s| s.case.clone()).collect();
            let mut kept: Vec<Step> = previous.into_iter().filter(|s| !rerun.contains(&s.case)).collect();
            kept.extend(all);
            all = kept;
        }
    }
    fs::write(&out, serde_json::to_string_pretty(&all).unwrap()).unwrap();
    let cover = choose_cover(&all);
    fs::write(root.join("xtask/fixtures-config/step-cover.toml"), &cover).unwrap();
    write_report(&root, &all, &cover);
    println!(
        "steps -> {}; cover -> xtask/fixtures-config/step-cover.toml",
        out.display()
    );
}

/// The call sites `xtask/fixtures-config/route-map.toml` names, whose per-step calls are recorded.
fn mapped_sites(root: &Path) -> Vec<String> {
    let Ok(text) = fs::read_to_string(root.join("xtask/fixtures-config/route-map.toml")) else {
        return vec![];
    };
    let map: toml::Value = toml::from_str(&text).expect("parse route-map.toml");
    map.get("site")
        .and_then(|s| s.as_array())
        .into_iter()
        .flatten()
        .map(|s| s["site"].as_str().expect("route-map.toml: site").to_string())
        .collect()
}

/// The operating points of a script, in order: `(line index, Some(point within its ASEQ))` for
/// the points of an `ASEQ a1 aN da` (alpha formed as XFOIL forms it, `AA1 + DAA·(p−1)`), `None`
/// for an `ALFA`/`CL` command.
fn operating_points(script: &str) -> Vec<(usize, Option<usize>)> {
    let mut out = vec![];
    for (n, l) in script.lines().enumerate() {
        let t: Vec<&str> = l.split_whitespace().collect();
        match t.first().copied() {
            Some("ALFA") | Some("CL") => out.push((n, None)),
            Some("ASEQ") => {
                let v: Vec<f64> = t[1..4].iter().map(|x| x.parse().unwrap()).collect();
                let np = if v[2] != 0.0 {
                    ((v[1] - v[0]) / v[2] + 0.5).floor() as usize + 1
                } else {
                    1
                };
                out.extend((1..=np).map(|p| (n, Some(p))));
            }
            _ => {}
        }
    }
    out
}

/// For each VISCAL call in order, the index (into `operating_points`) of the point it ran: an `ASEQ`
/// halts after NSEQEX = 4 consecutive points that do not converge, skipping the rest of that
/// sequence (`converged`, per call, from `viscal_points.dat`).
fn executed_points(script: &str, converged: &[bool]) -> Vec<usize> {
    const NSEQEX: usize = 4;
    let ops = operating_points(script);
    let mut out = vec![];
    let mut i = 0;
    while i < ops.len() && out.len() < converged.len() {
        let (line, seq) = ops[i];
        if seq.is_none() {
            out.push(i);
            i += 1;
            continue;
        }
        // one ASEQ: its points are ops[i..] with the same line
        let mut failures = 0;
        let mut j = i;
        while j < ops.len() && ops[j].0 == line {
            if out.len() >= converged.len() {
                break;
            }
            out.push(j);
            failures = if converged[out.len() - 1] { 0 } else { failures + 1 };
            j += 1;
            if failures >= NSEQEX {
                break;
            }
        }
        while i < ops.len() && ops[i].0 == line {
            i += 1;
        }
    }
    out
}

/// `LVCONV` of every VISCAL call, from the work directory's `viscal_points.dat`.
fn converged_per_call(work: &Path) -> Vec<bool> {
    fs::read_to_string(work.join("viscal_points.dat"))
        .unwrap_or_default()
        .lines()
        .filter_map(|l| l.split_once('='))
        .filter(|(k, _)| k.trim() == "LVCONV")
        .map(|(_, v)| v.trim() == "T")
        .collect()
}

/// The script truncated after step (c, i): everything before operating point `c`, then point
/// `c` (an index into `operating_points`) itself — given `ITER i` when `i > 0` (only for an `ALFA`/`CL` command: inside an `ASEQ`
/// an iteration cap would cap every point of the sequence, so an `ASEQ` point is a whole-call
/// step, `i = 0`, truncated as `ASEQ a1 ap da`, which forms the same alphas).
fn truncated_script(script: &str, c: usize, i: usize) -> String {
    let lines: Vec<&str> = script.lines().collect();
    let (line, seq) = operating_points(script)[c];
    let mut out = String::new();
    for l in &lines[..line] {
        if l.trim() == "QUIT" {
            break;
        }
        // the output commands (CPWR, DUMP) evaluate closures for what they write: they are not
        // part of any step, and would otherwise be counted in the next point's first step
        if l.starts_with("CPWR ") || l.starts_with("DUMP ") {
            continue;
        }
        out.push_str(l);
        out.push('\n');
    }
    match seq {
        None => {
            if i > 0 {
                out.push_str(&format!("ITER {i}\n"));
            }
            out.push_str(lines[line]);
            out.push('\n');
        }
        Some(p) => {
            assert_eq!(i, 0, "an ASEQ point is a whole-call step");
            let t: Vec<&str> = lines[line].split_whitespace().collect();
            let v: Vec<f64> = t[1..4].iter().map(|x| x.parse().unwrap()).collect();
            let ap = v[0] + v[2] * (p - 1) as f64;
            out.push_str(&format!("ASEQ {} {} {}\n", t[1], ap, t[3]));
        }
    }
    out.push_str("\nQUIT\n");
    out
}

/// `NITDONE` of every VISCAL call, from the work directory's `viscal_points.dat`.
fn iterations_per_call(work: &Path) -> Result<Vec<usize>, String> {
    let text = fs::read_to_string(work.join("viscal_points.dat")).map_err(|e| format!("viscal_points.dat: {e}"))?;
    Ok(text
        .lines()
        .filter_map(|l| l.split_once('='))
        .filter(|(k, _)| k.trim() == "NITDONE")
        .map(|(_, v)| v.trim().parse().unwrap())
        .collect())
}

/// `IT` rows of `viscal_iters_all.dat` keyed (call, iteration): RMSBL RLX CL CD CM ISTB IST
/// ITRAN1 ITRAN2 ALFA.
fn iteration_rows(dir: &Path) -> Result<BTreeMap<(usize, usize), Vec<f64>>, String> {
    let text = fs::read_to_string(dir.join("viscal_iters_all.dat")).map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut out = BTreeMap::new();
    for l in text.lines() {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(|t| t.parse().unwrap_or(f64::NAN)).collect();
        out.insert((v[0] as usize, v[1] as usize), v[2..12].to_vec());
    }
    Ok(out)
}

/// Whether step (c, i) is well-conditioned: every twin's branch trace equals the reference's
/// through it, and the twins' spread at it is within the fraction of the tolerance.
fn conditioning(work: &Path, nit: &[usize], c: usize, i: usize) -> (Option<bool>, String) {
    if i == 0 {
        // a whole-call step: every iteration of the call must be well-conditioned
        let mut worst = String::new();
        for ii in 1..=nit[c - 1] {
            let (ok, why) = conditioning(work, nit, c, ii);
            if ok != Some(true) {
                return (ok, format!("iteration {ii}: {why}"));
            }
            worst = why;
        }
        return (Some(true), format!("every iteration well-conditioned (last: {worst})"));
    }
    let Ok(base) = iteration_rows(work) else {
        return (None, "no reference iteration records".into());
    };
    let mut worst = 0.0_f64;
    for seed in SEEDS {
        let dir = work.join(format!("ulp{seed}"));
        if !dir.exists() {
            return (None, format!("no twin {seed}: run `cargo xtask twins --case <case>`"));
        }
        let Ok(tw) = iteration_rows(&dir) else {
            return (Some(false), format!("twin {seed} has no iteration records"));
        };
        let Ok(tnit) = iterations_per_call(&dir) else {
            return (Some(false), format!("twin {seed} has no points"));
        };
        for (cc, &n) in nit.iter().enumerate().take(c - 1) {
            if tnit.get(cc) != Some(&n) {
                return (
                    Some(false),
                    format!(
                        "twin {seed}: call {} takes {:?} iterations, not {n}",
                        cc + 1,
                        tnit.get(cc)
                    ),
                );
            }
        }
        let tn = tnit.get(c - 1).copied().unwrap_or(0);
        if tn < i || (i == nit[c - 1] && tn != i) {
            return (
                Some(false),
                format!("twin {seed}: call {c} takes {tn} iterations, not {}", nit[c - 1]),
            );
        }
        for (&(cc, ii), r) in base.range((1, 1)..=(c, i)) {
            let Some(t) = tw.get(&(cc, ii)) else {
                return (Some(false), format!("twin {seed}: no iteration ({cc}, {ii})"));
            };
            // ISTB IST ITRAN1 ITRAN2 are columns 5..9
            if r[5..9] != t[5..9] {
                return (
                    Some(false),
                    format!("twin {seed}: branch trace differs at ({cc}, {ii})"),
                );
            }
        }
        let (r, t) = (&base[&(c, i)], &tw[&(c, i)]);
        for col in [0usize, 1, 2, 3, 4, 9] {
            let (a, b) = (r[col], t[col]);
            let rel = (a - b).abs() / (TOL_SOLVER * a.abs().max(b.abs()).max(1.0));
            worst = if rel.is_nan() { f64::INFINITY } else { worst.max(rel) };
        }
    }
    let ok = worst <= CONDITIONING_FRACTION;
    (Some(ok), format!("twin spread {worst:.3} × TOL_SOLVER"))
}

/// One SPECAL (`specal_points.dat`) or SPECCL (`speccl_points.dat`) record: its scalars and the
/// per-node values, in call order.
fn inviscid_records(dir: &Path, file: &str) -> Vec<(HashMap<String, f64>, Vec<f64>)> {
    let Ok(text) = fs::read_to_string(dir.join(file)) else {
        return vec![];
    };
    let mut out: Vec<(HashMap<String, f64>, Vec<f64>)> = vec![];
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        let k = k.trim();
        if k == "CALL" {
            out.push((HashMap::new(), vec![]));
        }
        let Some(cur) = out.last_mut() else { continue };
        if k.starts_with("NODE") {
            cur.1
                .extend(v.split_whitespace().map(|t| t.parse::<f64>().unwrap_or(f64::NAN)));
        } else if let Ok(x) = v.trim().parse::<f64>() {
            cur.0.insert(k.to_string(), x);
        }
    }
    out
}

/// The branch events of the `events.dat` span (between `SPECAL_ENTER` markers) that contains
/// operating point `c` of the script; a CL point has no marker of its own, so it shares the span
/// of the SPECAL before it (or the span before the first marker).
fn event_span(dir: &Path, script: &str, c: usize) -> Vec<String> {
    let ops: Vec<bool> = script
        .lines()
        .map(str::trim_start)
        .filter(|l| l.starts_with("ALFA ") || l.starts_with("CL "))
        .map(|l| l.starts_with("ALFA "))
        .collect();
    let span = ops[..c].iter().filter(|&&a| a).count();
    let text = fs::read_to_string(dir.join("events.dat")).unwrap_or_default();
    let mut k = 0;
    let mut out = vec![];
    for l in text.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() < 2 {
            continue;
        }
        if t[1] == "SPECAL_ENTER" {
            k += 1;
        }
        if k == span {
            out.push(t[1..5].join(" "));
        }
    }
    out
}

/// Whether inviscid operating point `c` is well-conditioned: every twin's branch events in its
/// span and its SPECCL exit iteration equal the reference's, and the twins' spread of ALFA, CL,
/// CM, CDP, MINF and every node's values is within the fraction of the tolerance.
fn conditioning_inviscid(work: &Path, script: &str, c: usize) -> (Option<bool>, String) {
    let ops: Vec<bool> = script
        .lines()
        .map(str::trim_start)
        .filter(|l| l.starts_with("ALFA ") || l.starts_with("CL "))
        .map(|l| l.starts_with("ALFA "))
        .collect();
    let alfa = ops[c - 1];
    let (file, idx) = if alfa {
        ("specal_points.dat", ops[..c].iter().filter(|&&a| a).count())
    } else {
        ("speccl_points.dat", ops[..c].iter().filter(|&&a| !a).count())
    };
    let base = inviscid_records(work, file);
    let Some((bs, bn)) = base.get(idx - 1) else {
        return (None, format!("no reference record {idx} in {file}"));
    };
    let bev = event_span(work, script, c);
    let mut worst = 0.0_f64;
    for seed in SEEDS {
        let dir = work.join(format!("ulp{seed}"));
        if !dir.exists() {
            return (None, format!("no twin {seed}: run `cargo xtask twins --case <case>`"));
        }
        if event_span(&dir, script, c) != bev {
            return (Some(false), format!("twin {seed}: branch events differ"));
        }
        let tw = inviscid_records(&dir, file);
        let Some((ts, tn)) = tw.get(idx - 1) else {
            return (Some(false), format!("twin {seed}: no record {idx}"));
        };
        if bs.get("ITAL") != ts.get("ITAL") {
            return (
                Some(false),
                format!("twin {seed}: ITAL {:?} vs {:?}", ts.get("ITAL"), bs.get("ITAL")),
            );
        }
        let mut pairs: Vec<(f64, f64)> = ["ALFA", "CL", "CM", "CDP", "MINF"]
            .iter()
            .filter_map(|k| Some((*bs.get(*k)?, *ts.get(*k)?)))
            .collect();
        if bn.len() != tn.len() {
            return (Some(false), format!("twin {seed}: node count differs"));
        }
        pairs.extend(bn.iter().copied().zip(tn.iter().copied()));
        for (a, b) in pairs {
            let rel = (a - b).abs() / (TOL_SOLVER * a.abs().max(b.abs()).max(1.0));
            worst = if rel.is_nan() {
                if a.is_nan() && b.is_nan() {
                    worst
                } else {
                    f64::INFINITY
                }
            } else {
                worst.max(rel)
            };
        }
    }
    (
        Some(worst <= CONDITIONING_FRACTION),
        format!("twin spread {worst:.3} × TOL_SOLVER"),
    )
}

fn case_steps(root: &Path, gdir: &Path, case: &super::Case) -> Result<Vec<Step>, String> {
    let work = root.join("target/fixtures").join(&case.name);
    if !work.join("xfoil.inp").exists() {
        return Err(format!(
            "no work directory — run `cargo xtask fixtures --case {}`",
            case.name
        ));
    }
    let script = fs::read_to_string(work.join("xfoil.inp")).unwrap();
    let ops = operating_points(&script);
    // the point each call ran (an ASEQ that halted skips the rest of its points)
    let point_of: Vec<usize> = if case.inviscid {
        (0..ops.len()).collect()
    } else {
        executed_points(&script, &converged_per_call(&work))
    };
    let plan: Vec<(usize, usize)> = if case.inviscid {
        (1..=ops.len()).map(|c| (c, 0)).collect()
    } else {
        let nit = iterations_per_call(&work)?;
        if nit.len() != point_of.len() {
            return Err(format!(
                "{} VISCAL calls for {} executed points",
                nit.len(),
                point_of.len()
            ));
        }
        nit.iter()
            .enumerate()
            .flat_map(|(c, &n)| {
                let whole = ops[point_of[c]].1.is_some();
                (if whole { 0..=0 } else { 1..=n }).map(move |i| (c + 1, i))
            })
            .collect()
    };
    let nit = if case.inviscid {
        vec![]
    } else {
        iterations_per_call(&work)?
    };
    let run_dir = root.join("target/coverage/steps").join(&case.name);
    let _ = fs::remove_dir_all(&run_dir);
    fs::create_dir_all(&run_dir).unwrap();
    fs::copy(work.join("panels.dat"), run_dir.join("panels.dat")).unwrap();
    let scratch = run_dir.join("_gcov");
    let mut previous: HashMap<String, u64> = HashMap::new();
    // the calls of the script before its first operating point (LOAD, OPER's set-up), which are
    // not part of the first step's route
    let first_line = operating_points(&script)[0].0;
    let head: String = script.lines().take(first_line).map(|l| format!("{l}\n")).collect();
    fs::write(run_dir.join("xfoil.inp"), head + "\nQUIT\n").unwrap();
    super::coverage::reset_counters(gdir);
    let st = Command::new(gdir.join("bin/xfoil"))
        .current_dir(&run_dir)
        .stdin(fs::File::open(run_dir.join("xfoil.inp")).unwrap())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map_err(|e| e.to_string())?;
    if !st.success() {
        return Err(format!("the script's head: gcov xfoil exited {st}"));
    }
    let head_counts = super::coverage::branch_counts(root, gdir, &scratch);
    let mut previous_calls: HashMap<String, u64> = head_counts.calls;
    let mut previous_sites: HashMap<String, u64> = head_counts.sites;
    let mapped = mapped_sites(root);
    let mut out = vec![];
    for (n, &(c, i)) in plan.iter().enumerate() {
        fs::write(run_dir.join("xfoil.inp"), truncated_script(&script, point_of[c - 1], i)).unwrap();
        super::coverage::reset_counters(gdir);
        let st = Command::new(gdir.join("bin/xfoil"))
            .current_dir(&run_dir)
            .stdin(fs::File::open(run_dir.join("xfoil.inp")).unwrap())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map_err(|e| e.to_string())?;
        if !st.success() {
            return Err(format!("step ({c}, {i}): gcov xfoil exited {st}"));
        }
        let now = super::coverage::branch_counts(root, gdir, &scratch);
        let (counts, calls_now, sites_now) = (now.branches, now.calls, now.sites);
        let calls: BTreeMap<String, i64> = calls_now
            .iter()
            .map(|(s, &v)| (s.clone(), v as i64 - previous_calls.get(s).copied().unwrap_or(0) as i64))
            .filter(|(_, d)| *d != 0)
            .collect();
        let sites: BTreeMap<String, i64> = mapped
            .iter()
            .map(|s| {
                let d =
                    sites_now.get(s).copied().unwrap_or(0) as i64 - previous_sites.get(s).copied().unwrap_or(0) as i64;
                (s.clone(), d)
            })
            .filter(|(_, d)| *d != 0)
            .collect();
        let taken: Vec<String> = counts
            .iter()
            .filter(|(id, &v)| v > previous.get(*id).copied().unwrap_or(0))
            .map(|(id, _)| id.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        // a whole-call step is numbered by its first SETBL call
        let setbl = if case.inviscid {
            0
        } else {
            nit[..c - 1].iter().sum::<usize>() + i.max(1)
        };
        let (conditioned, why) = if case.inviscid {
            conditioning_inviscid(&work, &script, c)
        } else {
            conditioning(&work, &nit, c, i)
        };
        out.push(Step {
            case: case.name.clone(),
            call: c,
            iteration: i,
            setbl,
            call_iterations: if case.inviscid { 0 } else { nit[c - 1] },
            taken,
            calls,
            sites,
            conditioned,
            conditioning: why,
        });
        previous = counts;
        previous_calls = calls_now;
        previous_sites = sites_now;
        if n % 10 == 9 {
            println!("    {}: {}/{} steps", case.name, n + 1, plan.len());
        }
    }
    Ok(out)
}

/// The fewest steps taking every branch that any step takes: essential steps (the only step
/// taking some branch) first, then greedily by the number of branches still uncovered (ties: fewer
/// SETBL calls into the run, so the replay's fixture comes from early in it), then pruned of any step
/// the rest make redundant. Every step is eligible, ill-conditioned solutions included
/// (`docs/conventions/testing.md`, rule 4); the twins' conditioning is recorded beside each step for
/// information. A step listed in `xtask/fixtures-config/route.toml` as divergent (its replay takes
/// a different route from XFOIL's) or unmeasured (its replay could not be run) is excluded
/// (`cargo xtask route`).
fn choose_cover(all: &[Step]) -> String {
    let divergent = route_divergent();
    let good: Vec<&Step> = all
        .iter()
        .filter(|s| !divergent.contains(&(s.case.clone(), s.call, s.iteration)))
        .collect();
    let universe: BTreeSet<&String> = all.iter().flat_map(|s| &s.taken).collect();
    let coverable: BTreeSet<&String> = good.iter().flat_map(|s| &s.taken).collect();
    let mut chosen: Vec<&Step> = vec![];
    let mut left: BTreeSet<&String> = coverable.clone();
    for b in &coverable {
        let takers: Vec<&&Step> = good.iter().filter(|s| s.taken.contains(b)).collect();
        if takers.len() == 1 && !chosen.iter().any(|c| std::ptr::eq(*c, *takers[0])) {
            chosen.push(takers[0]);
        }
    }
    for s in &chosen {
        for b in &s.taken {
            left.remove(b);
        }
    }
    while !left.is_empty() {
        let best = good
            .iter()
            .max_by_key(|s| {
                let gain = s.taken.iter().filter(|b| left.contains(b)).count();
                (gain, std::cmp::Reverse(s.setbl))
            })
            .unwrap();
        for b in &best.taken {
            left.remove(b);
        }
        chosen.push(best);
    }
    // prune
    let mut k = 0;
    while k < chosen.len() {
        let rest: BTreeSet<&String> = chosen
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != k)
            .flat_map(|(_, s)| &s.taken)
            .collect();
        if coverable.iter().all(|b| rest.contains(b)) {
            chosen.remove(k);
        } else {
            k += 1;
        }
    }
    chosen.sort_by(|a, b| (&a.case, a.call, a.iteration).cmp(&(&b.case, b.call, b.iteration)));
    let ill: Vec<&&String> = universe.difference(&coverable).collect();
    let mut s = String::from(
        "# Generated by `cargo xtask steps` — do not edit by hand.\n#\n# The fewest steps of the candidate reference runs that take every branch any candidate step\n# takes, excluding the steps route.toml lists as divergent or unmeasured (`cargo xtask route`). Each\n# [[step]] is replayed by one test in tests/execution/branches.rs.\n",
    );
    s.push_str(&format!(
        "# {} branches taken by some step; {} by a step whose route agrees; {} steps chosen.\n",
        universe.len(),
        coverable.len(),
        chosen.len()
    ));
    s.push_str("\n# Branches taken only at steps whose replay diverges from XFOIL's route: no test gates these.\nroute_divergent_only = [\n");
    for b in ill {
        s.push_str(&format!("  \"{b}\",\n"));
    }
    s.push_str("]\n");
    for st in &chosen {
        s.push_str(&format!(
            "\n[[step]]\ncase = \"{}\"\ncall = {}\niteration = {}\nsetbl = {}\ncall_iterations = {}\n# {}\nbranches = [{}]\n",
            st.case,
            st.call,
            st.iteration,
            st.setbl,
            st.call_iterations,
            st.conditioning,
            st.taken
                .iter()
                .filter(|b| coverable.contains(b))
                .map(|b| format!("\"{b}\""))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    s
}

#[derive(serde::Deserialize)]
struct Cover {
    #[serde(default)]
    step: Vec<CoverStep>,
}

#[derive(serde::Deserialize)]
struct CoverStep {
    case: String,
    call: usize,
    iteration: usize,
    setbl: usize,
    call_iterations: usize,
}

/// Fold `xtask/fixtures-config/step-cover.toml` into the cases: every case with a chosen step is
/// tracked, dumps the SETBL calls the step replays, and keeps the files its test reads
/// (`tests/execution/step.rs`) in addition to its own `keep` list.
pub(crate) fn apply_cover(root: &Path, cases: &mut [super::Case]) {
    let Ok(text) = fs::read_to_string(root.join("xtask/fixtures-config/step-cover.toml")) else {
        return;
    };
    let cover: Cover = toml::from_str(&text).expect("parse step-cover.toml");
    for st in &cover.step {
        let case = cases
            .iter_mut()
            .find(|c| c.name == st.case)
            .unwrap_or_else(|| panic!("step-cover.toml names {}, which cases.toml does not define", st.case));
        let mut files: Vec<String> = vec![];
        let mut calls: Vec<usize> = vec![];
        if case.inviscid {
            let n_alpha = case.alphas.len();
            files.push(
                if st.call <= n_alpha {
                    "specal_points.dat"
                } else {
                    "speccl_points.dat"
                }
                .into(),
            );
            if st.call > n_alpha && st.call > 1 {
                files.push(
                    if st.call - 1 <= n_alpha {
                        "specal_points.dat"
                    } else {
                        "speccl_points.dat"
                    }
                    .into(),
                );
            }
        } else {
            let (first, last) = if st.iteration == 0 {
                (1, st.call_iterations)
            } else {
                (st.iteration, st.iteration)
            };
            let k0 = st.setbl;
            if !(st.call == 1 && first == 1) {
                files.push(format!("mrchdu_input_{k0}.dat"));
            }
            for it in first..=last {
                let k = k0 + (it - first);
                calls.push(k);
                files.push(format!("update_output_{k}.dat"));
                if it < st.call_iterations {
                    calls.push(k + 1);
                    files.push(format!("mrchdu_input_{}.dat", k + 1));
                }
            }
            files.push("viscal_points.dat".into());
        }
        case.track = true;
        for k in calls {
            if !case.dump_calls.contains(&k) {
                case.dump_calls.push(k);
            }
        }
        case.dump_calls.sort_unstable();
        for f in files {
            if !case.keep.contains(&f) {
                case.keep.push(f);
            }
        }
    }
}

/// For the VISCAL calls in `step_calls`, from the work directory's `viscal_points.dat`: the SETBL
/// calls to dump (every iteration of each call) and the files the single-step tests read
/// (`tests/execution/step.rs`: the entering dump of each iteration — the first iteration of a fresh
/// run excepted — UPDATE's output, and the call's converged point).
pub(crate) fn step_call_files(work: &Path, step_calls: &[usize]) -> (Vec<usize>, Vec<String>) {
    let nit = iterations_per_call(work).expect("viscal_points.dat");
    let mut calls = vec![];
    let mut files = vec!["viscal_points.dat".to_string()];
    for &c in step_calls {
        let k0: usize = nit[..c - 1].iter().sum::<usize>() + 1;
        for it in 1..=nit[c - 1] {
            let k = k0 + it - 1;
            calls.push(k);
            files.push(format!("update_output_{k}.dat"));
            if !(c == 1 && it == 1) {
                files.push(format!("mrchdu_input_{k}.dat"));
            }
        }
    }
    (calls, files)
}

/// The steps recorded as divergent in `xtask/fixtures-config/route.toml`, as
/// (case, call, iteration); empty before the route has been measured.
fn route_divergent() -> BTreeSet<(String, usize, usize)> {
    #[derive(serde::Deserialize)]
    struct Route {
        #[serde(default)]
        divergent: Vec<RouteStep>,
        #[serde(default)]
        unmeasured: Vec<RouteStep>,
    }
    #[derive(serde::Deserialize)]
    struct RouteStep {
        case: String,
        call: usize,
        iteration: usize,
    }
    let root = super::root();
    let Ok(text) = fs::read_to_string(root.join("xtask/fixtures-config/route.toml")) else {
        return BTreeSet::new();
    };
    let r: Route = toml::from_str(&text).expect("parse route.toml");
    r.divergent
        .into_iter()
        .chain(r.unmeasured)
        .map(|s| (s.case, s.call, s.iteration))
        .collect()
}

/// `docs/validation/branch-gating.md`: how full branch coverage is established — every candidate
/// case and how it is tested (from `cases.toml`), its steps and their conditioning, the cover and
/// the test of each chosen step, and the branches no step whose route agrees takes.
fn write_report(root: &Path, all: &[Step], cover: &str) {
    #[derive(serde::Deserialize)]
    struct Chosen {
        case: String,
        call: usize,
        iteration: usize,
        setbl: usize,
        branches: Vec<String>,
    }
    #[derive(serde::Deserialize)]
    struct CoverFile {
        #[serde(default)]
        step: Vec<Chosen>,
        #[serde(default)]
        route_divergent_only: Vec<String>,
    }
    let c: CoverFile = toml::from_str(cover).expect("parse the cover");
    let cases: super::Cases =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).unwrap()).unwrap();
    let how = |name: &str| -> String {
        let Some(k) = cases.cases.iter().find(|k| k.name == name) else {
            return "—".into();
        };
        let mut v = vec![];
        if k.run {
            v.push(if k.run_through > 0 {
                format!("run, calls 1–{}", k.run_through)
            } else {
                "run".into()
            });
        }
        if !k.step_calls.is_empty() {
            v.push(format!("every iteration of calls {:?} as steps", k.step_calls));
        }
        if k.events {
            v.push("events".into());
        }
        if v.is_empty() {
            "branch cover only".into()
        } else {
            v.join("; ")
        }
    };
    let universe: BTreeSet<&String> = all.iter().flat_map(|s| &s.taken).collect();
    let mut md = String::from(
        "# Branch gating\n\nGenerated by `cargo xtask steps` — do not edit. The method is in\n\
         [`docs/conventions/testing.md`](../conventions/testing.md), \"How branch coverage is established\"\n\
         and \"What gated means\".\n\n",
    );
    md.push_str(&format!(
        "{} branches of the translated subroutines are taken by some step of the {} candidate cases below; \
         {} steps are chosen to gate them all{}.\n\n",
        universe.len(),
        all.iter().map(|s| &s.case).collect::<BTreeSet<_>>().len(),
        c.step.len(),
        if c.route_divergent_only.is_empty() {
            String::new()
        } else {
            format!(
                ", except {} taken only at steps whose replay diverges from XFOIL's route",
                c.route_divergent_only.len()
            )
        }
    ));
    md.push_str("## Candidate cases\n\n| Case | How it is tested | Steps | Well-conditioned (twins) | Branches |\n|---|---|---:|---:|---:|\n");
    let names: BTreeSet<&String> = all.iter().map(|s| &s.case).collect();
    for n in names {
        let st: Vec<&Step> = all.iter().filter(|s| &s.case == n).collect();
        let taken: BTreeSet<&String> = st.iter().flat_map(|s| &s.taken).collect();
        md.push_str(&format!(
            "| `{n}` | {} | {} | {} | {} |\n",
            how(n),
            st.len(),
            st.iter().filter(|s| s.conditioned == Some(true)).count(),
            taken.len()
        ));
    }
    md.push_str("\n## The cover\n\nEach step is replayed by a test of the same name in `tests/execution/branches.rs`.\n\n| Case | Call | Iteration | SETBL | Branches |\n|---|---:|---:|---:|---:|\n");
    for s in &c.step {
        let it = if s.iteration == 0 {
            "all".to_string()
        } else {
            s.iteration.to_string()
        };
        md.push_str(&format!(
            "| `{}` | {} | {} | {} | {} |\n",
            s.case,
            s.call,
            it,
            s.setbl,
            s.branches.len()
        ));
    }
    // the route measurement (`cargo xtask route`): the steps excluded from the cover
    if let Ok(text) = fs::read_to_string(root.join("xtask/fixtures-config/route.toml")) {
        md.push_str("\n## The route\n\nA step is eligible for the cover only if yFoil is observed to take XFOIL's route through it: the same call count of every translated subroutine (`cargo xtask route`, [`testing.md`](../conventions/testing.md), \"What gated means\").");
        if let Some(summary) = text
            .lines()
            .find(|l| l.starts_with("# ") && l.contains("steps measured"))
        {
            md.push_str(&format!(" Measured: {}", summary.trim_start_matches("# ")));
        }
        let r: toml::Value = toml::from_str(&text).expect("parse route.toml");
        let excluded: Vec<&toml::Value> = ["divergent", "unmeasured"]
            .iter()
            .filter_map(|k| r.get(*k).and_then(|v| v.as_array()))
            .flatten()
            .collect();
        if excluded.is_empty() {
            md.push_str("\n\nNo step is excluded.\n");
        } else {
            md.push_str("\n\n| Case | Call | Iteration | Found by | Why |\n|---|---:|---:|---|---|\n");
            for d in excluded {
                md.push_str(&format!(
                    "| `{}` | {} | {} | {} | {} |\n",
                    d["case"].as_str().unwrap_or(""),
                    d["call"].as_integer().unwrap_or(0),
                    d["iteration"].as_integer().unwrap_or(0),
                    d.get("found").and_then(|f| f.as_str()).unwrap_or("not measured"),
                    d.get("why").and_then(|f| f.as_str()).unwrap_or("")
                ));
            }
        }
    }
    md.push_str("\n## Branches no test gates\n\n");
    if c.route_divergent_only.is_empty() {
        md.push_str("None among the branches the candidates take.\n");
    } else {
        for b in &c.route_divergent_only {
            md.push_str(&format!("- `{b}`\n"));
        }
    }
    // the branches the tracked cases take outside any solver step (`cargo xtask coverage`)
    if let Ok(text) = fs::read_to_string(root.join("target/coverage/branches.json")) {
        let u: serde_json::Value = serde_json::from_str(&text).unwrap_or_default();
        let mut by_sub: BTreeMap<String, usize> = BTreeMap::new();
        for b in u["branches"].as_array().into_iter().flatten() {
            let id = b["id"].as_str().unwrap_or("").to_string();
            if b["status"] == "taken" && !universe.contains(&id) {
                *by_sub
                    .entry(b["subroutine"].as_str().unwrap_or("?").to_string())
                    .or_default() += 1;
            }
        }
        if !by_sub.is_empty() {
            md.push_str("\n## Branches taken outside any solver step\n\nTaken by the tracked cases (`cargo xtask coverage`) but in no step of a candidate run: the geometry routines, which only the geometry-only cases run and `tests/subroutine/pangen.rs` gates against XFOIL's own PANGEN dumps.\n\n| Subroutine | Branches |\n|---|---:|\n");
            for (sub, n) in by_sub {
                md.push_str(&format!("| {sub} | {n} |\n"));
            }
        }
    }
    md.push_str("\nBranches of the translated subroutines that no case takes at all are listed, with how each could be reached, in [coverage.md](coverage.md).\n");
    fs::write(root.join("docs/validation/branch-gating.md"), md).unwrap();
}
