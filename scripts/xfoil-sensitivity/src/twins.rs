//! The twins mode: the reference against its own seeded 1-ULP twins, case by case.
//!
//! `cargo xtask twins` reruns every case on five seeded twins of its panels
//! (every coordinate jogged by −1, 0 or +1 ULP; CLAUDE.md Rule 1) and leaves the runs under
//! `target/fixtures/<case>/` (the reference) and `target/fixtures/<case>/ulp<seed>/`. This mode
//! reads those runs back, copies their per-call records into the run folder for reprocessing,
//! and derives, per recorded point (one VISCAL call: one alpha of one leg of the sweep):
//!
//! - the converged-point scalars — CL, CD, CM, CDF, CDP, the final RMSBL, the transition x/c
//!   of each side — and, from the final state record (`viscal_state_<k>.dat`), the
//!   boundary-layer scalars of each side: δ*, θ and H at the trailing-edge station, the largest
//!   H, and the trailing-edge separation x/c;
//! - how the call finished: converged (LVCONV) or stopped at its iteration limit (NITDONE =
//!   NITER1 without convergence);
//! - for every twin, the point that computed the same (leg, alpha), the absolute difference of
//!   every scalar and the largest absolute difference of every state array (per node: CPI, CPV,
//!   QINV, QVIS, GAM; per station: XSSI, X, UEDG, THET, DSTR, CTAU, MASS, TAU), and whether it
//!   finished differently from the reference;
//! - the envelope: the largest of those differences over the twins, per scalar and per array.
//!
//! Everything goes into `twins.json`; `plot.py` draws one sheet per case from it.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde_json::{json, Value};

/// The scalars recorded per point, in the order they are tabled
pub const SCALARS: [&str; 19] = [
    "cl",
    "cd",
    "cm",
    "cdf",
    "cdp",
    "rmsbl",
    "xtr_upper",
    "xtr_lower",
    "dstar_te_upper",
    "dstar_te_lower",
    "theta_te_upper",
    "theta_te_lower",
    "h_te_upper",
    "h_te_lower",
    "h_max_upper",
    "h_max_lower",
    "x_sep_upper",
    "x_sep_lower",
    "iterations",
];

const NODE_ARRAYS: [&str; 5] = ["CPI", "CPV", "QINV", "QVIS", "GAM"];
const STATION_ARRAYS: [&str; 8] = ["XSSI", "X", "UEDG", "THET", "DSTR", "CTAU", "MASS", "TAU"];

/// One recorded point of one run
#[derive(Debug, Clone)]
pub struct RunPoint {
    pub call: usize,
    pub leg: usize,
    pub alpha_deg: f64,
    pub converged: bool,
    pub iterations: usize,
    pub iteration_limit: usize,
    /// stopped at the iteration limit without converging
    pub capped: bool,
    pub scalars: BTreeMap<&'static str, f64>,
    /// per side (1 upper, 2 lower): the station arrays as recorded, by name
    pub arrays: BTreeMap<String, Vec<f64>>,
}

/// One run: the reference or one twin
#[derive(Debug, Clone)]
pub struct Run {
    pub label: String,
    pub seed: Option<u64>,
    pub points: Vec<RunPoint>,
    pub truncated_by_watchdog: bool,
}

fn real(s: &str) -> f64 {
    s.trim().parse().unwrap_or(f64::NAN)
}

/// `viscal_points.dat`: one block per VISCAL call, `KEY= value` lines
fn point_blocks(text: &str) -> Vec<BTreeMap<String, String>> {
    let mut out = vec![];
    for l in text.lines() {
        if let Some((k, v)) = l.split_once('=') {
            if k.trim() == "CALL" {
                out.push(BTreeMap::new());
            }
            if let Some(b) = out.last_mut() {
                b.insert(k.trim().to_string(), v.trim().to_string());
            }
        }
    }
    out
}

/// A state record's header (`KEY= value`) and its arrays by name
type StateRecord = (BTreeMap<String, String>, BTreeMap<String, Vec<f64>>);

/// The arrays of `viscal_state_<k>.dat` and its header: node arrays by name, station arrays
/// by name per side (`"UEDG/1"`, `"UEDG/2"`), and the header values
fn state(path: &Path) -> Option<StateRecord> {
    let text = fs::read_to_string(path).ok()?;
    let mut header = BTreeMap::new();
    let mut arrays: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for l in text.lines() {
        if let Some(rest) = l.strip_prefix("NODE(") {
            let (_, vals) = rest.split_once(")=")?;
            let v: Vec<f64> = vals.split_whitespace().map(real).collect();
            for (m, name) in NODE_ARRAYS.iter().enumerate() {
                arrays
                    .entry(name.to_string())
                    .or_default()
                    .push(v.get(m).copied().unwrap_or(f64::NAN));
            }
        } else if let Some(rest) = l.strip_prefix("BL(") {
            let (idx, vals) = rest.split_once(")=")?;
            let (is, _) = idx.split_once(',')?;
            let is = is.trim();
            let v: Vec<f64> = vals.split_whitespace().map(real).collect();
            for (m, name) in STATION_ARRAYS.iter().enumerate() {
                arrays
                    .entry(format!("{name}/{is}"))
                    .or_default()
                    .push(v.get(m).copied().unwrap_or(f64::NAN));
            }
        } else if let Some((k, v)) = l.split_once('=') {
            header.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    Some((header, arrays))
}

/// The boundary-layer scalars of one side from its station arrays and the header: (δ*, θ, H at
/// the trailing-edge station, the largest H on the airfoil, the trailing-edge separation x/c)
fn side_scalars(header: &BTreeMap<String, String>, arrays: &BTreeMap<String, Vec<f64>>, is: usize) -> [f64; 5] {
    let get = |name: &str| arrays.get(&format!("{name}/{is}")).cloned().unwrap_or_default();
    let (x, theta, dstar, tau) = (get("X"), get("THET"), get("DSTR"), get("TAU"));
    let te: usize = header
        .get(&format!("IBLTE{is}"))
        .and_then(|v| v.parse().ok())
        .unwrap_or(theta.len());
    // rows are stations 1..=NBL (index ibl − 1); the trailing-edge station is IBLTE
    let row = |v: &Vec<f64>, ibl: usize| v.get(ibl - 1).copied().unwrap_or(f64::NAN);
    let (d, t) = (row(&dstar, te), row(&theta, te));
    let h_te = d / t;
    let h_max = (2..=te)
        .map(|ibl| row(&dstar, ibl) / row(&theta, ibl))
        .filter(|h| h.is_finite())
        .fold(f64::NAN, f64::max);
    // walk back from the trailing edge while τ ≤ 0: the separation point is the first station of
    // that run; attached at the TE when the run is empty
    let mut sep = te + 1;
    while sep > 2 && row(&tau, sep - 1) <= 0.0 {
        sep -= 1;
    }
    let x_sep = if sep > te { 1.0 } else { row(&x, sep) };
    [d, t, h_te, h_max, x_sep]
}

/// Read an inviscid run back: one SPECAL call per point (`specal_points.dat`: ALFA, CL, CM, CDP
/// and the per-node GAM, QINV, CPI). SPECAL is direct: every call converges, none is capped.
pub fn load_inviscid_run(dir: &Path, label: &str, seed: Option<u64>) -> Option<Run> {
    let text = fs::read_to_string(dir.join("specal_points.dat")).ok()?;
    let mut points = vec![];
    let mut block: Option<StateRecord> = None;
    let mut blocks = vec![];
    for l in text.lines() {
        if let Some(rest) = l.strip_prefix("NODE(") {
            if let (Some((_, arrays)), Some((_, vals))) = (block.as_mut(), rest.split_once(")=")) {
                let v: Vec<f64> = vals.split_whitespace().map(real).collect();
                for (m, name) in ["GAM", "QINV", "CPI"].iter().enumerate() {
                    arrays
                        .entry(name.to_string())
                        .or_default()
                        .push(v.get(m).copied().unwrap_or(f64::NAN));
                }
            }
        } else if let Some((k, v)) = l.split_once('=') {
            if k.trim() == "CALL" {
                if let Some(b) = block.take() {
                    blocks.push(b);
                }
                block = Some((BTreeMap::new(), BTreeMap::new()));
            }
            if let Some((h, _)) = block.as_mut() {
                h.insert(k.trim().to_string(), v.trim().to_string());
            }
        }
    }
    if let Some(b) = block.take() {
        blocks.push(b);
    }
    for (i, (h, arrays)) in blocks.into_iter().enumerate() {
        let g = |k: &str| real(h.get(k).map(String::as_str).unwrap_or("NaN"));
        let mut scalars: BTreeMap<&'static str, f64> = BTreeMap::new();
        for name in SCALARS {
            scalars.insert(name, f64::NAN);
        }
        scalars.insert("cl", g("CL"));
        scalars.insert("cm", g("CM"));
        scalars.insert("cdp", g("CDP"));
        scalars.insert("iterations", 0.0);
        points.push(RunPoint {
            call: i + 1,
            leg: 1,
            alpha_deg: g("ALFA").to_degrees(),
            converged: true,
            iterations: 0,
            iteration_limit: 0,
            capped: false,
            scalars,
            arrays,
        });
    }
    Some(Run {
        label: label.to_string(),
        seed,
        points,
        truncated_by_watchdog: false,
    })
}

/// Read a run directory (the reference or a twin) back
pub fn load_run(dir: &Path, label: &str, seed: Option<u64>, first_after: Option<f64>) -> Option<Run> {
    let text = fs::read_to_string(dir.join("viscal_points.dat")).ok()?;
    let blocks = point_blocks(&text);
    let truncated = fs::read_to_string(dir.join("manifest.json"))
        .ok()
        .and_then(|t| serde_json::from_str::<Value>(&t).ok())
        .and_then(|m| m["truncated_by_watchdog"].as_bool())
        .unwrap_or(false);
    let mut points = vec![];
    let mut leg = 1usize;
    let mut prev_alpha = f64::NAN;
    for (i, b) in blocks.iter().enumerate() {
        let call = i + 1;
        let alpha_deg = real(b.get("ALFA").map(String::as_str).unwrap_or("NaN")).to_degrees();
        // leg 2 starts where the recorded alpha returns to the seed alpha of the second script
        // after a sequence that was moving away from it
        if let Some(a0) = first_after {
            if leg == 1 && i > 0 && (alpha_deg - a0).abs() < 1e-6 && (prev_alpha - a0).abs() > 0.5 {
                leg = 2;
            }
        }
        prev_alpha = alpha_deg;
        let converged = b.get("LVCONV").map(|v| v == "T").unwrap_or(false);
        let iterations: usize = b.get("NITDONE").and_then(|v| v.parse().ok()).unwrap_or(0);
        let limit: usize = b.get("NITER1").and_then(|v| v.parse().ok()).unwrap_or(0);
        let mut scalars: BTreeMap<&'static str, f64> = BTreeMap::new();
        let g = |k: &str| real(b.get(k).map(String::as_str).unwrap_or("NaN"));
        scalars.insert("cl", g("CL"));
        scalars.insert("cd", g("CD"));
        scalars.insert("cm", g("CM"));
        scalars.insert("cdf", g("CDF"));
        scalars.insert("cdp", g("CDP"));
        scalars.insert("rmsbl", g("RMSBL"));
        scalars.insert("xtr_upper", g("XOCTR1"));
        scalars.insert("xtr_lower", g("XOCTR2"));
        scalars.insert("iterations", iterations as f64);
        let (header, arrays) = state(&dir.join(format!("viscal_state_{call}.dat"))).unwrap_or_default();
        for (is, side) in [(1usize, "upper"), (2, "lower")] {
            let [d, t, h, hmax, xsep] = side_scalars(&header, &arrays, is);
            scalars.insert(
                if side == "upper" {
                    "dstar_te_upper"
                } else {
                    "dstar_te_lower"
                },
                d,
            );
            scalars.insert(
                if side == "upper" {
                    "theta_te_upper"
                } else {
                    "theta_te_lower"
                },
                t,
            );
            scalars.insert(if side == "upper" { "h_te_upper" } else { "h_te_lower" }, h);
            scalars.insert(if side == "upper" { "h_max_upper" } else { "h_max_lower" }, hmax);
            scalars.insert(if side == "upper" { "x_sep_upper" } else { "x_sep_lower" }, xsep);
        }
        points.push(RunPoint {
            call,
            leg,
            alpha_deg,
            converged,
            iterations,
            iteration_limit: limit,
            capped: !converged && limit > 0 && iterations >= limit,
            scalars,
            arrays,
        });
    }
    Some(Run {
        label: label.to_string(),
        seed,
        points,
        truncated_by_watchdog: truncated,
    })
}

/// How a twin's point is matched to the reference's: by the (leg, alpha) it computed, or — for
/// a sweep of prescribed CL points, whose alpha is what the solve finds — by call number
fn key(p: &RunPoint, by_call: bool) -> (usize, i64) {
    if by_call {
        (p.leg, p.call as i64)
    } else {
        (p.leg, (p.alpha_deg * 1000.0).round() as i64)
    }
}

/// Largest |a − b| over the entries of two arrays; NaN when the lengths differ or an entry is
/// non-finite on one side only (both NaN counts as agreement)
fn array_diff(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() {
        return f64::NAN;
    }
    let mut worst = 0.0_f64;
    for (x, y) in a.iter().zip(b) {
        if x.is_nan() && y.is_nan() {
            continue;
        }
        let d = (x - y).abs();
        if d.is_nan() {
            return f64::NAN;
        }
        worst = worst.max(d);
    }
    worst
}

/// A case's twins record for `twins.json`
pub fn case_json(name: &str, section: &str, script: &str, reference: &Run, twins: &[Run], by_call: bool) -> Value {
    let point_scalars = |p: &RunPoint| -> Value {
        let mut o = serde_json::Map::new();
        o.insert("call".into(), json!(p.call));
        o.insert("leg".into(), json!(p.leg));
        o.insert("alpha_deg".into(), json!(p.alpha_deg));
        o.insert("converged".into(), json!(p.converged));
        o.insert("iterations".into(), json!(p.iterations));
        o.insert("iteration_limit".into(), json!(p.iteration_limit));
        o.insert("capped".into(), json!(p.capped));
        for name in SCALARS {
            o.insert(name.into(), json!(p.scalars.get(name).copied().unwrap_or(f64::NAN)));
        }
        Value::Object(o)
    };
    let ref_points: Vec<Value> = reference.points.iter().map(point_scalars).collect();
    let mut twins_json = vec![];
    // envelope accumulators per reference point: per scalar and per array
    let mut env_scalar: Vec<BTreeMap<&'static str, f64>> = vec![BTreeMap::new(); reference.points.len()];
    let mut env_array: Vec<BTreeMap<String, f64>> = vec![BTreeMap::new(); reference.points.len()];
    let mut env_completion: Vec<bool> = vec![false; reference.points.len()];
    let mut env_present: Vec<usize> = vec![0; reference.points.len()];
    let fold_max = |slot: &mut f64, d: f64| {
        // a NaN difference (the twin is non-finite where the reference is not) is the worst
        if slot.is_nan() {
            return;
        }
        *slot = if d.is_nan() { f64::NAN } else { slot.max(d) };
    };
    for t in twins {
        let by_key: BTreeMap<(usize, i64), &RunPoint> = t.points.iter().map(|p| (key(p, by_call), p)).collect();
        let mut pts = vec![];
        for (i, r) in reference.points.iter().enumerate() {
            let Some(p) = by_key.get(&key(r, by_call)) else {
                continue;
            };
            env_present[i] += 1;
            let mut o = point_scalars(p);
            let obj = o.as_object_mut().unwrap();
            let mut diff = serde_json::Map::new();
            for name in SCALARS {
                let (a, b) = (
                    p.scalars.get(name).copied().unwrap_or(f64::NAN),
                    r.scalars.get(name).copied().unwrap_or(f64::NAN),
                );
                let d = if a.is_nan() && b.is_nan() { 0.0 } else { (a - b).abs() };
                diff.insert(name.into(), json!(d));
                fold_max(env_scalar[i].entry(name).or_insert(0.0), d);
            }
            let mut adiff = serde_json::Map::new();
            for (aname, ra) in &r.arrays {
                let d = p.arrays.get(aname).map(|pa| array_diff(pa, ra)).unwrap_or(f64::NAN);
                adiff.insert(aname.clone(), json!(d));
                fold_max(env_array[i].entry(aname.clone()).or_insert(0.0), d);
            }
            let completion_differs = p.converged != r.converged || p.capped != r.capped;
            env_completion[i] |= completion_differs;
            obj.insert("diff".into(), Value::Object(diff));
            obj.insert("array_diff".into(), Value::Object(adiff));
            obj.insert("completion_differs".into(), json!(completion_differs));
            pts.push(o);
        }
        twins_json.push(json!({
            "label": t.label, "seed": t.seed, "recorded_points": t.points.len(),
            "truncated_by_watchdog": t.truncated_by_watchdog,
            "completion_differs": pts.iter().filter(|p| p["completion_differs"].as_bool().unwrap_or(false)).count(),
            "points": pts,
        }));
    }
    let envelope: Vec<Value> = reference
        .points
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut o = serde_json::Map::new();
            o.insert("call".into(), json!(r.call));
            o.insert("leg".into(), json!(r.leg));
            o.insert("alpha_deg".into(), json!(r.alpha_deg));
            o.insert("twins_present".into(), json!(env_present[i]));
            o.insert("completion_differs".into(), json!(env_completion[i]));
            let mut sc = serde_json::Map::new();
            for name in SCALARS {
                sc.insert(name.into(), json!(env_scalar[i].get(name).copied().unwrap_or(f64::NAN)));
            }
            o.insert("scalars".into(), Value::Object(sc));
            let mut ar = serde_json::Map::new();
            for (k, v) in &env_array[i] {
                ar.insert(k.clone(), json!(v));
            }
            o.insert("arrays".into(), Value::Object(ar));
            Value::Object(o)
        })
        .collect();
    // the extents of every plotted quantity over the reference's converged points, and of the
    // finite differences over every twin: what the axes are set from
    let mut extents = serde_json::Map::new();
    for name in SCALARS {
        let extent = |converged_only: bool| {
            reference
                .points
                .iter()
                .filter(|p| p.converged || !converged_only)
                .map(|p| p.scalars.get(name).copied().unwrap_or(f64::NAN))
                .filter(|v| v.is_finite())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)))
        };
        // over the converged points; over every finite point when none converged
        let (lo, hi) = match extent(true) {
            (lo, hi) if lo.is_finite() => (lo, hi),
            _ => extent(false),
        };
        let (dlo, dhi) = env_scalar
            .iter()
            .filter_map(|m| m.get(name))
            .copied()
            .filter(|v| v.is_finite() && *v > 0.0)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
        extents.insert(
            name.into(),
            json!({ "min": lo, "max": hi, "diff_min": dlo, "diff_max": dhi }),
        );
    }
    let (alo, ahi) = reference
        .points
        .iter()
        .map(|p| p.alpha_deg)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
    extents.insert("alpha_deg".into(), json!({ "min": alo, "max": ahi }));
    json!({
        "name": name, "section": section, "script": script,
        "seeds": twins.iter().filter_map(|t| t.seed).collect::<Vec<_>>(),
        "reference": {
            "recorded_points": reference.points.len(),
            "converged": reference.points.iter().filter(|p| p.converged).count(),
            "capped": reference.points.iter().filter(|p| p.capped).count(),
            "truncated_by_watchdog": reference.truncated_by_watchdog,
            "points": ref_points,
        },
        "twins": twins_json,
        "envelope": envelope,
        "extents": extents,
    })
}

/// Copy a run's per-call records into the run folder for reprocessing
pub fn archive(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap();
    for entry in fs::read_dir(src).unwrap().flatten() {
        let n = entry.file_name().to_string_lossy().to_string();
        let keep = n.starts_with("viscal_")
            || n.starts_with("specal_")
            || n.starts_with("speccl_")
            || n == "noise_floor.json"
            || n == "manifest.json"
            || n == "xfoil.inp"
            || n == "panels.dat"
            || n == "panels.json"
            || n == "dump_calls.txt";
        if keep && entry.path().is_file() {
            fs::copy(entry.path(), dst.join(&n)).unwrap();
        }
    }
}

/// The work directory of a case's reference run (`cargo xtask fixtures`) and of its twins (`cargo xtask twins`)
pub fn work_dirs(root: &Path, case: &str) -> (PathBuf, Vec<(u64, PathBuf)>) {
    let work = root.join("target/fixtures").join(case);
    let mut twins = vec![];
    for entry in fs::read_dir(&work).into_iter().flatten().flatten() {
        let n = entry.file_name().to_string_lossy().to_string();
        if let Some(seed) = n.strip_prefix("ulp").and_then(|s| s.parse::<u64>().ok()) {
            if entry.path().join("viscal_points.dat").exists() || entry.path().join("specal_points.dat").exists() {
                twins.push((seed, entry.path()));
            }
        }
    }
    twins.sort();
    (work, twins)
}

/// `README.md` of the run folder for the twins mode
pub fn write_index(run_dir: &Path, j: &Value) {
    fs::write(
        run_dir.join("README.md"),
        format!("# XFOIL against its 1-ULP twins\n\n{}", body_md(j, run_dir, true)),
    )
    .unwrap();
}

fn table_md(j: &Value) -> String {
    let mut t = String::from(
        "| case | section | N | ITER | OPER | reference points | converged | at the iteration limit | twins | twin points (min–max) | finished differently (any twin) | worst envelope CL | worst envelope CD | worst envelope $x_{tr}$ upper | worst envelope δ* TE upper |\n|---|---|---:|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|\n",
    );
    for c in j["cases"].as_array().unwrap() {
        let r = &c["reference"];
        let twins = c["twins"].as_array().unwrap();
        let counts: Vec<usize> = twins
            .iter()
            .map(|t| t["recorded_points"].as_u64().unwrap_or(0) as usize)
            .collect();
        let differs = c["envelope"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|e| e["completion_differs"].as_bool().unwrap_or(false))
            .count();
        // the worst envelope over the reference's converged points
        let worst = |name: &str| -> String {
            // n/a for a scalar the run does not have (the boundary-layer scalars of an
            // inviscid case)
            if r["points"]
                .as_array()
                .unwrap()
                .iter()
                .all(|p| p[name].as_f64().is_none())
            {
                return "n/a".into();
            }
            let v = c["envelope"]
                .as_array()
                .unwrap()
                .iter()
                .zip(r["points"].as_array().unwrap())
                .filter(|(_, p)| p["converged"].as_bool().unwrap_or(false))
                .map(|(e, _)| e["scalars"][name].as_f64().unwrap_or(f64::NAN))
                .fold(f64::NEG_INFINITY, |a, b| {
                    if b.is_nan() || a.is_nan() {
                        f64::NAN
                    } else {
                        a.max(b)
                    }
                });
            if v.is_nan() {
                "NaN".into()
            } else if v == f64::NEG_INFINITY {
                "—".into()
            } else {
                format!("{v:.1e}")
            }
        };
        t += &format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | {} | {}–{} | {} | {} | {} | {} | {} |\n",
            c["name"].as_str().unwrap(),
            c["section"].as_str().unwrap(),
            c["n_nodes"],
            c["max_iterations"],
            c["script"].as_str().unwrap(),
            r["recorded_points"],
            r["converged"],
            r["capped"],
            twins.len(),
            counts.iter().min().copied().unwrap_or(0),
            counts.iter().max().copied().unwrap_or(0),
            differs,
            worst("cl"),
            worst("cd"),
            worst("xtr_upper"),
            worst("dstar_te_upper"),
        );
    }
    t
}

fn body_md(j: &Value, run_dir: &Path, figures: bool) -> String {
    let mut m = String::new();
    m += &format!(
        "Generated by `cargo run --release -p xfoil-sensitivity -- --twins --group {}` on {} (git {}{}). \
         Provenance in `metadata.json`; every number in `twins.json`; the per-call records of the reference and \
         of every twin under `twins/<case>/`, copied from the `cargo xtask fixtures` and `cargo xtask twins` runs they were read from.\n\n",
        j["group"].as_str().unwrap_or(""),
        j["created_utc"].as_str().unwrap_or(""),
        j["commit"].as_str().map(|s| &s[..s.len().min(9)]).unwrap_or(""),
        if j["dirty"].as_bool().unwrap_or(false) {
            ", dirty"
        } else {
            ""
        }
    );
    m += "The cases are the `twins-baseline` group of `xtask/fixtures-config/cases.toml`: the series-cases \
          sections and conditions at 240 panels and ITER 200, swept 0 → +30°, `INIT`, 0 → −30° by 1° with XFOIL's \
          polar procedure (`ALFA 0 / ASEQ`, ITMAX + 5 per sequence point, the sequence halting after NSEQEX = 4 \
          consecutive unconverged points), so that the reference is given every chance to converge and the \
          panelling is fine enough not to be the story; each case's N and ITER are recorded with its results. \
          The instrumented double-precision XFOIL 6.99 reference run on each case as generated, and on five \
          seeded 1-ULP twins of its panels — every coordinate jogged by −1, 0 or +1 ULP, x and y independently, \
          the jog of each coordinate drawn from splitmix64 of the seed and the node, so a twin is reproducible from \
          its seed (`cargo xtask twins`; `docs/conventions/terminology.md`, *twin*). Per recorded point — one VISCAL call, one alpha of \
          one leg of the sweep — the converged-point scalars (CL, CD, CM, CDF, CDP, the final RMSBL, the \
          transition x/c of each side, the iterations taken), the boundary-layer scalars of each side from the \
          final state record (δ*, θ and H at the trailing-edge station, the largest H, the trailing-edge \
          separation x/c), and how the call finished: converged, or stopped at its iteration limit. A twin's point \
          is matched to the reference's by (leg, α); its absolute difference in every scalar and the largest \
          absolute difference in every final-state array (per node CPI, CPV, QINV, QVIS, GAM; per station and side \
          XSSI, X, UEDG, THET, DSTR, CTAU, MASS, TAU) are recorded, and the **envelope** is the largest of those \
          over the twins, per scalar and per array, per point. These envelopes are the reference's own \
          sensitivity to 1 ULP of its input at every point of every case: the floor against which yFoil is gated.\n\n";
    m += "One sheet per case. **Left column**, the baseline: CL, CD, the upper-side transition x/c and δ* at the \
          upper-side trailing-edge station against α, the reference in black on top, the twins beneath it in \
          YlGnBu (one colour per seed, darker for a higher seed), open circles where the call did not converge, \
          an orange tick at the top of the panel where the reference stopped at its iteration limit. **Right \
          column**: log |twin − reference| of the same quantity, each twin in its colour, the envelope in red, \
          and an orange diamond on a twin's line where that twin finished differently from the reference (one \
          converged and the other stopped at its limit). Both legs of a sweep are drawn as one curve, ascending \
          in α. Where a twin has no point at a reference α (its sequence halted earlier) the line has a gap.\n\n";
    if figures {
        for c in j["cases"].as_array().unwrap() {
            let name = c["name"].as_str().unwrap();
            m += &format!(
                "### {}: `{}`\n\nN = {} nodes, ITER {}; {}\n\n![{name}](twins-{name}.svg)\n\n",
                c["section"].as_str().unwrap(),
                name,
                c["n_nodes"],
                c["max_iterations"],
                c["script"].as_str().unwrap()
            );
        }
    }
    m += "## Summary\n\nWorst envelope: the largest envelope over the reference's converged points — an O(1) value there \
          means a twin left the reference's solution branch at a point the reference converged on, which is the \
          reference being unstable against one ULP at that point; NaN where a twin was non-finite at a point \
          the reference was not; — where the reference converged nowhere.\n\n";
    m += &table_md(j);
    let _ = run_dir;
    m
}

fn tex_escape(s: &str) -> String {
    s.replace('\\', "\\textbackslash{}")
        .replace('&', "\\&")
        .replace('%', "\\%")
        .replace('#', "\\#")
        .replace('_', "\\_")
        .replace('α', "$\\alpha$")
        .replace('—', "---")
        .replace('–', "--")
        .replace('°', "$^\\circ$")
}

/// `README.tex`: the sheets and the summary table for the paper
pub fn write_index_tex(run_dir: &Path, j: &Value) {
    use std::io::Write;
    let mut f = fs::File::create(run_dir.join("README.tex")).unwrap();
    let w = |f: &mut fs::File, s: &str| writeln!(f, "{s}").unwrap();
    w(
        &mut f,
        "% Generated by `cargo run --release -p xfoil-sensitivity -- --twins`; table from twins.json. Needs booktabs.",
    );
    w(&mut f, "\\section*{XFOIL against its 1-ULP twins}");
    w(&mut f, "\\noindent The reference run on each case and on five seeded 1-ULP twins of its panels (every coordinate jogged by $-1$, $0$ or $+1$ ULP from the seed); per recorded point the converged-point and boundary-layer scalars, each twin's differences and the envelope (the largest over the twins). Provenance in \\texttt{metadata.json}, every number in \\texttt{twins.json}.");
    for c in j["cases"].as_array().unwrap() {
        let name = c["name"].as_str().unwrap();
        w(&mut f, &format!("\\begin{{figure}}[p]\\centering\\includegraphics{{twins-{name}.pdf}}\\caption{{{} (\\texttt{{{}}}), {}: left, $C_L$, $C_D$, upper-side transition $x/c$ and $\\delta^*$ at the upper-side trailing edge against $\\alpha$ — the reference black on top, the twins in YlGnBu beneath, open circles unconverged, an orange tick where the reference stopped at its iteration limit; right, $|$twin $-$ reference$|$ per twin, the envelope in red, an orange diamond where a twin finished differently from the reference.}}\\end{{figure}}", tex_escape(c["section"].as_str().unwrap()), tex_escape(name), tex_escape(c["script"].as_str().unwrap())));
    }
    w(&mut f, "\\begin{table}[h]\\centering\\small");
    w(&mut f, "\\caption{Per case: the reference's recorded points, how many converged and how many stopped at the iteration limit, the twins, and the worst envelope of $C_L$, $C_D$, the upper transition $x/c$ and the upper trailing-edge $\\delta^*$ over the converged points.}");
    w(&mut f, "\\begin{tabular}{@{}lrrrrrrrrrrr@{}}\\toprule");
    w(
        &mut f,
        "Case & $N$ & ITER & Points & Conv. & Limit & Twins & Differ & $C_L$ & $C_D$ & $x_{tr}$ & $\\delta^*_{TE}$ \\\\ \\midrule",
    );
    for c in j["cases"].as_array().unwrap() {
        let r = &c["reference"];
        let differs = c["envelope"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|e| e["completion_differs"].as_bool().unwrap_or(false))
            .count();
        let worst = |name: &str| -> String {
            // n/a for a scalar the run does not have (the boundary-layer scalars of an
            // inviscid case)
            if r["points"]
                .as_array()
                .unwrap()
                .iter()
                .all(|p| p[name].as_f64().is_none())
            {
                return "n/a".into();
            }
            let v = c["envelope"]
                .as_array()
                .unwrap()
                .iter()
                .zip(r["points"].as_array().unwrap())
                .filter(|(_, p)| p["converged"].as_bool().unwrap_or(false))
                .map(|(e, _)| e["scalars"][name].as_f64().unwrap_or(f64::NAN))
                .fold(f64::NEG_INFINITY, |a, b| {
                    if b.is_nan() || a.is_nan() {
                        f64::NAN
                    } else {
                        a.max(b)
                    }
                });
            if v.is_nan() {
                "NaN".into()
            } else if v == f64::NEG_INFINITY {
                "---".into()
            } else {
                format!("{v:.1e}")
            }
        };
        w(
            &mut f,
            &format!(
                "\\texttt{{{}}} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} & {} \\\\",
                tex_escape(c["name"].as_str().unwrap()),
                c["n_nodes"],
                c["max_iterations"],
                r["recorded_points"],
                r["converged"],
                r["capped"],
                c["twins"].as_array().unwrap().len(),
                differs,
                worst("cl"),
                worst("cd"),
                worst("xtr_upper"),
                worst("dstar_te_upper")
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
}

/// The documentation page: `docs/validation/xfoil-sensitivity/README.md` and the sheets (SVG)
pub fn write_docs(docs: &Path, run_dir: &Path, j: &Value) {
    fs::create_dir_all(docs).unwrap();
    let mut m = String::from("# XFOIL against its 1-ULP twins: the reference's own noise floor, case by case\n\n");
    m +=
        "Generated by `cargo run --release -p xfoil-sensitivity -- --twins --docs`; do not edit — regenerate. The run \
          folder (`scripts/xfoil-sensitivity/runs/`, gitignored) holds `twins.json` with every number and the per-call \
          records of the reference and every twin for reprocessing.\n\n";
    m += &body_md(j, run_dir, true);
    fs::write(docs.join("README.md"), m).unwrap();
    for old in fs::read_dir(docs).unwrap().flatten() {
        if old.file_name().to_string_lossy().ends_with(".svg") {
            fs::remove_file(old.path()).unwrap();
        }
    }
    for c in j["cases"].as_array().unwrap() {
        let name = format!("twins-{}.svg", c["name"].as_str().unwrap());
        if run_dir.join(&name).exists() {
            fs::copy(run_dir.join(&name), docs.join(&name)).unwrap();
        }
    }
}
