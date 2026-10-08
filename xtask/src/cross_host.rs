//! `cargo xtask cross-host`: the reference's own spread between two hosts' maths libraries, per
//! whole-run case and per compared value — the floor of a cross-host gate (CLAUDE.md Rule 1).
//!
//! Both codes take `exp`, `log`, `pow`, `tanh`… from the host's libm, and Apple's and glibc's differ
//! by one ULP on a fraction of inputs, so a fixture generated on one host cannot be met more closely
//! on another than the reference itself moves between them. `TOL_CROSS_HOST` covers one replayed
//! step; a whole run of an ill-conditioned case (an exhausted Newton, a sweep into a stall break)
//! moves further. This command, run on the *other* host (`scripts/cross-host.sh` does so in a Linux
//! container), reruns the instrumented reference in a work directory on each tracked run case's own
//! `panels.dat` and `xfoil.inp` — the same geometry bit for bit — and compares its records with the
//! tracked ones over exactly what the run tests compare: every iteration and converged point of the
//! VISCAL calls through `run_through`, and every inviscid point with its nodes. The metric is the
//! tests' `|a − b| / max(|a|, |b|, scale)`. It writes `tests/fixtures/cross-host/spread.json`:
//! both hosts and their libraries, per case the largest spread of each value, whether the
//! reference's route (iteration counts, IST, ITRAN) differs between the hosts, and the length and
//! FNV-1a hash of every tracked record it was measured against, so a test can tell a stale
//! measurement from a current one. The tests gate a run cross-host at
//! `max(tolerance, TOL_CROSS_HOST, CROSS_HOST_FACTOR × spread)` (`tests/common/utilities/`).

use crate::{run_xfoil, Case, Cases, RunEnd};
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::process::Command;

/// FNV-1a (64-bit) of a file's bytes, written out so the test side computes the same
/// (`tests/common/utilities/cross_host.rs`).
pub(crate) fn fnv64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// A real the reference wrote with `ES24.16`, three-digit exponents (`1.96+116`) included.
fn real(s: &str) -> f64 {
    let s = s.trim();
    if let Ok(v) = s.parse() {
        return v;
    }
    match s.get(1..).and_then(|t| t.rfind(['+', '-'])) {
        Some(j) => format!("{}e{}", &s[..=j], &s[j + 1..]).parse().unwrap_or(f64::NAN),
        None => f64::NAN,
    }
}

fn err(a: f64, b: f64, scale: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(scale)
}

/// `IT k i …` rows by call (1-based), as `tests/execution/run.rs` reads them.
fn iterations(text: &str) -> BTreeMap<usize, Vec<Vec<f64>>> {
    let mut out: BTreeMap<usize, Vec<Vec<f64>>> = BTreeMap::new();
    for l in text.lines() {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(real).collect();
        out.entry(v[0] as usize).or_default().push(v[1..].to_vec());
    }
    out
}

/// One record: its `KEY = value` scalars and its `NODE` rows.
type Record = (BTreeMap<String, f64>, Vec<Vec<f64>>);

/// `KEY = value` blocks opened by `CALL`, with the `NODE` rows of each block kept apart.
fn blocks(text: &str) -> Vec<Record> {
    let mut out: Vec<Record> = vec![];
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        let k = k.trim();
        if k == "CALL" {
            out.push((BTreeMap::new(), vec![]));
        }
        let Some(cur) = out.last_mut() else { continue };
        if k.starts_with("NODE") {
            cur.1.push(v.split_whitespace().map(real).collect());
        } else {
            cur.0.insert(k.to_string(), real(v));
        }
    }
    out
}

fn widen(m: &mut Map<String, Value>, name: &str, e: f64) {
    let cur = m.get(name).and_then(Value::as_f64).unwrap_or(0.0);
    // a NaN spread (a value non-finite on one host only) is recorded as such: no floor covers it
    if e.is_nan() || cur.is_nan() {
        m.insert(name.to_string(), json!(f64::NAN.to_string()));
    } else if e > cur || !m.contains_key(name) {
        m.insert(name.to_string(), json!(e));
    }
}

/// The spread of one case: (per-iteration, per-point and per-node maps, route note).
fn compare(case: &Case, ours: &Path, theirs: &Path) -> (Value, Option<String>) {
    let (mut it_m, mut pt_m, mut node_m) = (Map::new(), Map::new(), Map::new());
    let mut route = None;
    let read = |d: &Path, f: &str| fs::read_to_string(d.join(f)).unwrap_or_default();
    if case.inviscid {
        for f in ["specal_points.dat", "speccl_points.dat"] {
            for (a, b) in blocks(&read(theirs, f)).iter().zip(blocks(&read(ours, f)).iter()) {
                for k in ["ALFA", "CL", "CM", "CDP", "MINF"] {
                    if let (Some(x), Some(y)) = (a.0.get(k), b.0.get(k)) {
                        widen(&mut pt_m, k, err(*x, *y, 1.0));
                    }
                }
                if a.0.get("ITAL") != b.0.get("ITAL") && route.is_none() {
                    route = Some(format!("{f}: ITAL differs"));
                }
                for (r, s) in a.1.iter().zip(b.1.iter()) {
                    for (j, k) in ["GAM", "QINV", "CPI"].iter().enumerate() {
                        widen(&mut node_m, k, err(r[j], s[j], 1.0));
                    }
                }
            }
        }
    } else {
        let through = if case.run_through == 0 {
            usize::MAX
        } else {
            case.run_through
        };
        let (a, b) = (
            iterations(&read(theirs, "viscal_iters_all.dat")),
            iterations(&read(ours, "viscal_iters_all.dat")),
        );
        // columns after the call number: it RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2 ALFA MINF REINF
        let cols = [
            ("RMSBL", 1),
            ("RLX", 2),
            ("CL", 3),
            ("CD", 4),
            ("CM", 5),
            ("ALFA", 10),
            ("MINF", 11),
            ("REINF", 12),
        ];
        for (k, rows) in a.iter().filter(|(k, _)| **k <= through) {
            let other = b.get(k).cloned().unwrap_or_default();
            if rows.len() != other.len() && route.is_none() {
                route = Some(format!(
                    "call {k}: iteration count {} here, {} there",
                    other.len(),
                    rows.len()
                ));
            }
            for (r, s) in rows.iter().zip(other.iter()) {
                if r[7..10] != s[7..10] && route.is_none() {
                    route = Some(format!("call {k} iteration {}: IST/ITRAN differ", r[0]));
                }
                for (name, j) in cols {
                    let scale = if name == "REINF" { r[j].abs() } else { 1.0 };
                    widen(&mut it_m, name, err(r[j], s[j], scale));
                }
            }
        }
        let (pa, pb) = (
            blocks(&read(theirs, "viscal_points.dat")),
            blocks(&read(ours, "viscal_points.dat")),
        );
        for (n, (x, y)) in pa.iter().zip(pb.iter()).enumerate() {
            if n + 1 > through {
                break;
            }
            for k in [
                "ALFA", "CL", "CM", "CD", "CDF", "CDP", "XOCTR1", "XOCTR2", "MINF", "REINF",
            ] {
                if let (Some(p), Some(q)) = (x.0.get(k), y.0.get(k)) {
                    let scale = if k == "REINF" { p.abs() } else { 1.0 };
                    widen(&mut pt_m, k, err(*p, *q, scale));
                }
            }
        }
    }
    let mut v = json!({ "iteration": it_m, "point": pt_m });
    if !node_m.is_empty() {
        v["node"] = Value::Object(node_m);
    }
    (v, route)
}

fn host() -> (String, String) {
    let u = |a: &str| {
        Command::new("uname")
            .arg(a)
            .output()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .unwrap_or_default()
    };
    let libm = Command::new("ldd")
        .arg("--version")
        .output()
        .map(|o| {
            String::from_utf8_lossy(&o.stdout)
                .lines()
                .next()
                .unwrap_or("")
                .to_string()
        })
        .unwrap_or_default();
    (format!("{}-{}", u("-m"), u("-s")), libm)
}

pub(crate) fn run(flags: &[String]) {
    let root = crate::root();
    let out_path = flags
        .windows(2)
        .find(|w| w[0] == "--out")
        .map(|w| root.join(&w[1]))
        .unwrap_or_else(|| root.join("tests/fixtures/cross-host/spread.json"));
    let cases: Cases = toml::from_str(
        &fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("xtask/fixtures-config/cases.toml"),
    )
    .expect("parse cases.toml");
    let xfoil = root.join("target/xfoil-ref/instrumented/bin/xfoil");
    if !xfoil.exists() {
        crate::xfoil_build(&[]);
    }
    let (here, libm) = host();
    let mut fixture_host: Option<(String, String)> = None;
    let mut out = Map::new();
    for case in cases.cases.iter().filter(|c| c.track && c.run) {
        let fixture = root.join("tests/fixtures/xfoil").join(&case.name);
        let manifest: Value =
            serde_json::from_str(&fs::read_to_string(fixture.join("manifest.json")).unwrap()).unwrap();
        let line = |key: &str| -> String {
            manifest["xfoil_ref"]
                .as_array()
                .and_then(|a| a.iter().filter_map(Value::as_str).find(|s| s.starts_with(key)))
                .map(|s| s[key.len()..].trim().to_string())
                .unwrap_or_default()
        };
        let fh = (line("host:"), line("libm:"));
        match &fixture_host {
            None => fixture_host = Some(fh.clone()),
            Some(h) if *h != fh => panic!("{}: fixtures from two hosts ({:?}, {:?})", case.name, h, fh),
            _ => {}
        }
        let work = root.join("target/cross-host").join(&case.name);
        let _ = fs::remove_dir_all(&work);
        fs::create_dir_all(&work).unwrap();
        for f in ["panels.dat", "xfoil.inp"] {
            fs::copy(fixture.join(f), work.join(f)).unwrap_or_else(|e| panic!("{}: {f}: {e}", case.name));
        }
        print!("{} ... ", case.name);
        match run_xfoil(&xfoil, &work) {
            Ok(RunEnd::Finished) | Ok(RunEnd::Hung) => {}
            other => panic!("{}: the reference did not run to its end: {other:?}", case.name),
        }
        let (spread, route) = compare(case, &work, &fixture);
        let mut files = Map::new();
        for f in [
            "viscal_iters_all.dat",
            "viscal_points.dat",
            "specal_points.dat",
            "speccl_points.dat",
        ] {
            if let Ok(b) = fs::read(fixture.join(f)) {
                files.insert(
                    f.to_string(),
                    json!({ "bytes": b.len(), "fnv64": format!("{:016x}", fnv64(&b)) }),
                );
            }
        }
        println!("{}", route.as_deref().unwrap_or("same route"));
        out.insert(
            case.name.clone(),
            json!({ "run_through": case.run_through, "route_differs": route, "records": files, "spread": spread }),
        );
    }
    let (fh, fl) = fixture_host.expect("no tracked run cases");
    assert_ne!(
        fh.split('-').nth(1).map(str::to_ascii_lowercase),
        here.split('-').nth(1).map(str::to_ascii_lowercase),
        "cross-host measures the fixture host against another: run it on a host whose libm differs from the fixtures' ({fh})"
    );
    let doc = json!({
        "about": "The reference's own spread between the fixtures' host and another host's maths library, per tracked \
                  run case and per compared value, in the tests' metric |a - b| / max(|a|, |b|, scale). On a host whose \
                  libm is `measured_host.libm`, the run tests gate each value at max(tolerance, TOL_CROSS_HOST, \
                  CROSS_HOST_FACTOR x spread). Written by `cargo xtask cross-host` (scripts/cross-host.sh); never edited.",
        "fixture_host": { "host": fh, "libm": fl },
        "measured_host": { "host": here, "libm": libm },
        "cases": out,
    });
    fs::create_dir_all(out_path.parent().unwrap()).unwrap();
    fs::write(&out_path, serde_json::to_string_pretty(&doc).unwrap() + "\n").unwrap();
    println!("spread -> {}", out_path.display());
}
