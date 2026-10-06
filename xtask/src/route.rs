//! `cargo xtask route [--case NAME]... [--reuse]`: is yFoil's route through each candidate step
//! XFOIL's? The strong meaning of "gated" (`docs/conventions/testing.md`, "What gated means").
//!
//! XFOIL's route through a step is its call count of every translated subroutine in the step,
//! measured by `cargo xtask steps` from the gcov build (`target/coverage/steps.json`, `calls`).
//! yFoil's is measured by replaying the same step from XFOIL's dumped state in a build of the
//! `execution` tests with LLVM source coverage (`-C instrument-coverage --cfg yfoil_route`,
//! `tests/execution/route.rs`, which describes how its profile windows line up with XFOIL's
//! prefix-run differences), and reading the call count of every function carrying the
//! subroutine's `#[doc(alias)]`.
//!
//! 1. **Dumps.** Every candidate case is rerun by the instrumented reference in
//!    `target/route/fixtures/<case>/` with every SETBL call dumped (only the `mrchdu_input_<k>`
//!    files the replays read are kept). `--reuse` keeps an existing run.
//! 2. **Replays.** The coverage build replays every step; its profiles go to
//!    `target/route/profiles/`.
//! 3. **Comparison.** Per step and subroutine: a subroutine translated by one function is compared
//!    by call count; one split over several functions (BLDIF's equations, say), by whether it was
//!    called at all. Subroutines with no alias are not compared, and are listed. Where yFoil
//!    makes no call at an XFOIL call site — work done inline, or a call whose result nothing
//!    reads — `xtask/fixtures-config/route-map.toml` names the site with the reason, and the
//!    calls XFOIL made from it in the step (gcov's per-line call counts) are taken off XFOIL's
//!    counts before comparing.
//!
//! Output: `target/route/route.json` (every step, every subroutine's two counts) and
//! `xtask/fixtures-config/route.toml`, the steps whose route differs (`[[divergent]]`, read by
//! the step cover, which then excludes them) and those that could not be replayed
//! (`[[unmeasured]]`, excluded too). Entries found by other means (`found` other than "route")
//! are kept.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

#[derive(serde::Deserialize, Clone)]
struct Step {
    case: String,
    call: usize,
    iteration: usize,
    #[serde(default)]
    calls: BTreeMap<String, i64>,
    #[serde(default)]
    sites: BTreeMap<String, i64>,
}

#[derive(serde::Serialize)]
struct Compared {
    case: String,
    call: usize,
    iteration: usize,
    measured: bool,
    /// subroutine → (XFOIL's calls, yFoil's calls, how compared, agrees)
    subroutines: BTreeMap<String, (i64, i64, &'static str, bool)>,
}

/// One translated function: the subroutine its alias names, and how it appears in a profile.
struct Translation {
    alias: String,
    /// `yfoil::solver::setbl`, from the source file
    module: String,
    name: String,
}

/// Every `#[doc(alias = "X")]` on a function under `src/`.
fn translations(root: &Path) -> Vec<Translation> {
    fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
        for e in fs::read_dir(dir).unwrap().flatten() {
            let p = e.path();
            if p.is_dir() {
                walk(&p, out);
            } else if p.extension().is_some_and(|x| x == "rs") {
                out.push(p);
            }
        }
    }
    let mut files = vec![];
    walk(&root.join("src"), &mut files);
    let mut out = vec![];
    for f in files {
        let rel = f.strip_prefix(root.join("src")).unwrap().with_extension("");
        let mut parts: Vec<String> = rel.iter().map(|s| s.to_string_lossy().to_string()).collect();
        if parts.last().is_some_and(|l| l == "mod" || l == "lib" || l == "main") {
            parts.pop();
        }
        if parts.first().is_some_and(|p| p == "bin") {
            continue;
        }
        let module = std::iter::once("yfoil".to_string())
            .chain(parts)
            .collect::<Vec<_>>()
            .join("::");
        let text = fs::read_to_string(&f).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        for (n, l) in lines.iter().enumerate() {
            let Some(rest) = l.trim().strip_prefix("#[doc(alias = \"") else {
                continue;
            };
            let alias = rest.split('"').next().unwrap().to_string();
            // the next `fn` within the attribute block
            for m in &lines[n + 1..(n + 8).min(lines.len())] {
                let t = m.trim();
                if t.starts_with("#[") || t.starts_with("///") {
                    continue;
                }
                if let Some(i) = t.find("fn ") {
                    let name: String = t[i + 3..]
                        .chars()
                        .take_while(|c| c.is_alphanumeric() || *c == '_')
                        .collect();
                    out.push(Translation {
                        alias: alias.clone(),
                        module: module.clone(),
                        name,
                    });
                }
                break;
            }
        }
    }
    out
}

fn llvm_bin(tool: &str) -> PathBuf {
    let sysroot = Command::new("rustc")
        .args(["--print", "sysroot"])
        .output()
        .expect("rustc");
    let sysroot = String::from_utf8(sysroot.stdout).unwrap().trim().to_string();
    for e in fs::read_dir(Path::new(&sysroot).join("lib/rustlib")).unwrap().flatten() {
        let p = e.path().join("bin").join(tool);
        if p.exists() {
            return p;
        }
    }
    panic!("{tool} not found: rustup component add llvm-tools");
}

/// Every yFoil function's entry count in one raw profile, by demangled path (hash stripped).
fn profile_counts(profdata: &Path, profile: &Path) -> BTreeMap<String, i64> {
    let out = Command::new(profdata)
        .args(["show", "--all-functions"])
        .arg(profile)
        .output()
        .expect("llvm-profdata show");
    assert!(out.status.success(), "llvm-profdata show {}", profile.display());
    let mut counts = BTreeMap::new();
    let mut current: Option<String> = None;
    for l in String::from_utf8_lossy(&out.stdout).lines() {
        let t = l.trim();
        if l.starts_with("  ") && !l.starts_with("    ") && t.ends_with(':') {
            let sym = &t[..t.len() - 1];
            current = Some(format!("{:#}", rustc_demangle::demangle(sym)));
        } else if let (Some(v), Some(name)) = (t.strip_prefix("Function count: "), current.as_ref()) {
            *counts.entry(name.clone()).or_insert(0) += v.parse::<i64>().unwrap();
        }
    }
    counts
}

/// Does the demangled path `path` name the function `t`? A free function is `module::name`; a
/// method `<module::Type>::name` or `<module::Type as Trait>::name`.
fn names(path: &str, t: &Translation) -> bool {
    let Some(head) = path.strip_suffix(&format!("::{}", t.name)) else {
        return false;
    };
    let head = head.trim_start_matches('<');
    head == t.module
        || head
            .strip_prefix(&format!("{}::", t.module))
            .is_some_and(|ty| !ty.contains("::") || ty.contains(" as "))
}

pub fn route(flags: &[String]) {
    let root = super::root();
    let reuse = flags.iter().any(|f| f == "--reuse");
    let selected: Vec<&str> = flags
        .windows(2)
        .filter(|w| w[0] == "--case")
        .map(|w| w[1].as_str())
        .collect();
    let steps: Vec<Step> = serde_json::from_str(
        &fs::read_to_string(root.join("target/coverage/steps.json"))
            .expect("target/coverage/steps.json: run `cargo xtask steps` first"),
    )
    .expect("parse steps.json");
    let steps: Vec<Step> = steps
        .into_iter()
        .filter(|s| selected.is_empty() || selected.contains(&s.case.as_str()))
        .collect();
    let cases: BTreeSet<String> = steps.iter().map(|s| s.case.clone()).collect();
    let base = root.join("target/route");

    // 1. the reference runs, every SETBL call dumped
    let xfoil = root.join("target/xfoil-ref/instrumented/bin/xfoil");
    assert!(xfoil.exists(), "build the reference first: cargo xtask xfoil-build");
    for case in &cases {
        let work = root.join("target/fixtures").join(case);
        let dir = base.join("fixtures").join(case);
        if reuse && dir.join("manifest.json").exists() {
            continue;
        }
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for f in ["manifest.json", "panels.json", "panels.dat", "xfoil.inp"] {
            if work.join(f).exists() {
                fs::copy(work.join(f), dir.join(f)).unwrap();
            }
        }
        let total: usize = fs::read_to_string(work.join("viscal_points.dat"))
            .map(|t| {
                t.lines()
                    .filter_map(|l| l.split_once('='))
                    .filter(|(k, _)| k.trim() == "NITDONE")
                    .map(|(_, v)| v.trim().parse::<usize>().unwrap())
                    .sum()
            })
            .unwrap_or(0);
        assert!(total <= 1000, "{case}: {total} SETBL calls, the dump list holds 1000");
        if total > 0 {
            let list: Vec<String> = (1..=total).map(|k| k.to_string()).collect();
            fs::write(dir.join("dump_calls.txt"), list.join("\n") + "\n").unwrap();
        }
        let st = Command::new(&xfoil)
            .current_dir(&dir)
            .stdin(fs::File::open(dir.join("xfoil.inp")).unwrap())
            .stdout(fs::File::create(dir.join("stdout.txt")).unwrap())
            .stderr(Stdio::null())
            .status()
            .expect("run xfoil");
        assert!(st.success(), "{case}: xfoil exited {st}");
        // keep what the replays read
        for e in fs::read_dir(&dir).unwrap().flatten() {
            let n = e.file_name().to_string_lossy().to_string();
            let keep = n.starts_with("mrchdu_input_")
                || [
                    "manifest.json",
                    "panels.json",
                    "panels.dat",
                    "xfoil.inp",
                    "viscal_points.dat",
                    "dump_calls.txt",
                ]
                .contains(&n.as_str());
            if !keep {
                let _ = fs::remove_file(e.path());
            }
        }
        println!("  {case}: reference rerun, {total} SETBL calls dumped");
    }

    // 2. the replays, in the coverage build
    let build = Command::new("cargo")
        .args([
            "test",
            "--release",
            "--test",
            "execution",
            "--no-run",
            "--message-format=json",
        ])
        .env("RUSTFLAGS", "-C instrument-coverage --cfg yfoil_route")
        .env("CARGO_TARGET_DIR", base.join("build"))
        // build scripts are instrumented too, and write a profile when they run
        .env("LLVM_PROFILE_FILE", base.join("build/build-script-%p.profraw"))
        .current_dir(&root)
        .stderr(Stdio::inherit())
        .output()
        .expect("cargo test --no-run");
    assert!(build.status.success(), "coverage build failed");
    let exe = String::from_utf8_lossy(&build.stdout)
        .lines()
        .filter_map(|l| serde_json::from_str::<serde_json::Value>(l).ok())
        .filter_map(|v| v["executable"].as_str().map(PathBuf::from))
        .next_back()
        .expect("the execution test binary");
    let profiles = base.join("profiles");
    let _ = fs::remove_dir_all(&profiles);
    fs::create_dir_all(&profiles).unwrap();
    let list = base.join("steps.txt");
    let lines: Vec<String> = steps
        .iter()
        .map(|s| format!("{} {} {}", s.case, s.call, s.iteration))
        .collect();
    fs::write(&list, lines.join("\n") + "\n").unwrap();
    let run = Command::new(&exe)
        .args(["route::measure", "--exact", "--nocapture"])
        .env("YFOIL_ROUTE_STEPS", &list)
        .env("YFOIL_ROUTE_OUT", &profiles)
        .env("YFOIL_ROUTE_FIXTURES", base.join("fixtures"))
        .env("LLVM_PROFILE_FILE", base.join("default_%p.profraw"))
        .current_dir(&root)
        .stderr(Stdio::null())
        .output()
        .expect("run the route measurement");
    let stdout = String::from_utf8_lossy(&run.stdout);
    let unmeasured: BTreeSet<String> = stdout
        .lines()
        .filter_map(|l| l.strip_prefix("ROUTE-UNMEASURED "))
        .map(|s| s.trim().to_string())
        .collect();

    // 3. the comparison
    let profdata = llvm_bin("llvm-profdata");
    let all = translations(&root);
    let mut by_alias: BTreeMap<String, Vec<&Translation>> = BTreeMap::new();
    for t in &all {
        by_alias.entry(t.alias.clone()).or_default().push(t);
    }
    let spec: toml::Value =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/coverage.toml")).expect("coverage.toml"))
            .unwrap();
    let subroutines: Vec<String> = spec["file"]
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|f| f["subroutines"].as_array().unwrap().iter())
        .map(|s| s.as_str().unwrap().to_string())
        .collect();
    let untranslated: Vec<&String> = subroutines.iter().filter(|s| !by_alias.contains_key(*s)).collect();
    // The call sites where yFoil's translation makes no call although XFOIL does
    // (`xtask/fixtures-config/route-map.toml`, each with its reason): work yFoil does inline
    // (SETBL's COMSET) and calls whose result nothing reads, which yFoil does not translate
    // (SETBL's BLMID(3) at the trailing edge, the BLMID calls of the MRCHUE/MRCHDU fallbacks). The
    // calls XFOIL made from them in a step are measured by gcov (`steps.json`, `sites`) and taken
    // off XFOIL's counts before comparing, with the calls the dropped routine makes in turn (BLMID(1)
    // calls CFL; BLMID(2) CFT and CFL), so yFoil is compared with what XFOIL did less what it
    // translates without a call.
    let map: toml::Value =
        toml::from_str(&fs::read_to_string(root.join("xtask/fixtures-config/route-map.toml")).expect("route-map.toml"))
            .expect("parse route-map.toml");
    // (site, subroutine, calls of the subroutine per call made at the site)
    let uncalled: Vec<(String, String, i64)> = map
        .get("site")
        .and_then(|o| o.as_array())
        .into_iter()
        .flatten()
        .flat_map(|o| {
            let site = o["site"].as_str().unwrap().to_string();
            o["calls"]
                .as_table()
                .unwrap()
                .iter()
                .map(move |(sub, n)| (site.clone(), sub.clone(), n.as_integer().unwrap()))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut compared = vec![];
    for s in &steps {
        let key = format!("{} {} {}", s.case, s.call, s.iteration);
        let tag = format!("{}_{}_{}", s.case, s.call, s.iteration);
        let mut parts: Vec<(i64, PathBuf)> = vec![];
        for e in fs::read_dir(&profiles).unwrap().flatten() {
            let n = e.file_name().to_string_lossy().to_string();
            let Some(part) = n
                .strip_prefix(&format!("{tag}_"))
                .and_then(|p| p.strip_suffix(".profraw"))
            else {
                continue;
            };
            let sign = if part.starts_with('B') { -1 } else { 1 };
            parts.push((sign, e.path()));
        }
        let measured = !unmeasured.contains(&key) && !parts.is_empty();
        let mut yfoil: BTreeMap<String, i64> = BTreeMap::new();
        if measured {
            for (sign, p) in &parts {
                for (f, n) in profile_counts(&profdata, p) {
                    *yfoil.entry(f).or_insert(0) += sign * n;
                }
            }
        }
        let mut subs = BTreeMap::new();
        if measured {
            for sub in &subroutines {
                let Some(ts) = by_alias.get(sub) else { continue };
                let counts: Vec<i64> = ts
                    .iter()
                    .map(|t| yfoil.iter().filter(|(p, _)| names(p, t)).map(|(_, n)| *n).sum())
                    .collect();
                // XFOIL's calls less those made from the sites yFoil translates without a call
                let x = s.calls.get(sub).copied().unwrap_or(0)
                    - uncalled
                        .iter()
                        .filter(|(_, u, _)| u == sub)
                        .map(|(site, _, n)| n * s.sites.get(site).copied().unwrap_or(0))
                        .sum::<i64>();
                let (y, how, agrees) = if ts.len() == 1 {
                    (counts[0], "count", counts[0] == x)
                } else {
                    let y = counts.iter().copied().max().unwrap_or(0);
                    (y, "called", (y > 0) == (x > 0))
                };
                if x != 0 || y != 0 {
                    subs.insert(sub.clone(), (x, y, how, agrees));
                }
            }
        }
        compared.push(Compared {
            case: s.case.clone(),
            call: s.call,
            iteration: s.iteration,
            measured,
            subroutines: subs,
        });
    }
    fs::write(
        base.join("route.json"),
        serde_json::to_string_pretty(&compared).unwrap(),
    )
    .unwrap();

    // route.toml: what was found by other means is kept
    let path = root.join("xtask/fixtures-config/route.toml");
    let kept: Vec<toml::Value> = fs::read_to_string(&path)
        .ok()
        .and_then(|t| toml::from_str::<toml::Value>(&t).ok())
        .and_then(|v| v.get("divergent").and_then(|d| d.as_array().cloned()))
        .unwrap_or_default()
        .into_iter()
        .filter(|d| d.get("found").and_then(|f| f.as_str()) != Some("route"))
        .collect();
    let mut s = String::from(
        "# Steps whose replay does not take XFOIL's route: divergent comparisons\n# (docs/conventions/terminology.md). The step cover excludes them (`cargo xtask steps`).\n# Entries with found = \"route\" are written by `cargo xtask route` — do not edit them by hand;\n# the others record a divergence found by other means, with the evidence.\n",
    );
    let diverged: Vec<&Compared> = compared
        .iter()
        .filter(|c| c.measured && c.subroutines.values().any(|v| !v.3))
        .collect();
    let not_measured: Vec<&Compared> = compared.iter().filter(|c| !c.measured).collect();
    s.push_str(&format!(
        "#\n# {} steps measured, {} divergent, {} not measured; subroutines with no translation alias (not compared): {}.\n",
        compared.len() - not_measured.len(),
        diverged.len(),
        not_measured.len(),
        untranslated.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
    ));
    for d in &kept {
        s.push_str("\n[[divergent]]\n");
        s.push_str(&toml::to_string(d).unwrap());
    }
    for c in &diverged {
        let why: Vec<String> = c
            .subroutines
            .iter()
            .filter(|(_, v)| !v.3)
            .map(|(k, v)| format!("{k} XFOIL {} yFoil {}", v.0, v.1))
            .collect();
        if kept.iter().any(|d| {
            d.get("case").and_then(|x| x.as_str()) == Some(c.case.as_str())
                && d.get("call").and_then(|x| x.as_integer()) == Some(c.call as i64)
                && d.get("iteration").and_then(|x| x.as_integer()) == Some(c.iteration as i64)
        }) {
            continue;
        }
        s.push_str(&format!(
            "\n[[divergent]]\ncase = \"{}\"\ncall = {}\niteration = {}\nfound = \"route\"\nwhy = \"{}\"\n",
            c.case,
            c.call,
            c.iteration,
            why.join("; ")
        ));
    }
    for c in &not_measured {
        s.push_str(&format!(
            "\n[[unmeasured]]\ncase = \"{}\"\ncall = {}\niteration = {}\n",
            c.case, c.call, c.iteration
        ));
    }
    fs::write(&path, s).unwrap();
    println!(
        "route: {} steps measured, {} divergent, {} not measured -> {}",
        compared.len() - not_measured.len(),
        diverged.len(),
        not_measured.len(),
        path.display()
    );
}
