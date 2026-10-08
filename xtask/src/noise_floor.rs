//! The reference's noise floor from its seeded 1-ULP twins (study data, never read by a test):
//! `cargo xtask twins` runs the twins (`steps.rs`) and then writes `noise_floor.json` beside the
//! case's work directory from the base run and the `ulp<seed>` directories — per VISCAL call the
//! largest absolute spread of every per-point and per-iteration value over the twins, the per-array
//! spread of the dumped post-UPDATE state and of the final state, and every branch flip. The
//! branch-coverage and series-cases studies read it (`docs/conventions/terminology.md`,
//! *noise floor*). Ported unchanged from the fixture pipeline that wrote it before 2026-10-06.

use std::fs;
use std::path::Path;

/// Write `noise_floor.json` for the case in `work` from its twins `ulp<seed>`, `seeds`; returns a
/// one-line summary.
pub(crate) fn write(work: &Path, inviscid: bool, seeds: &[u64]) -> Result<String, String> {
    let mut merged: Option<serde_json::Value> = None;
    for &t in seeds {
        let ulp = work.join(format!("ulp{t}"));
        let j = if inviscid {
            inviscid_twin_summary(work, &ulp)?
        } else {
            twin_record(work, &ulp)?
        };
        merged = Some(match merged {
            Some(m) => merge_max(&m, &j),
            None => j,
        });
    }
    let mut merged = merged.ok_or("no twins")?;
    merged["method"] = serde_json::json!(format!(
        "{} seeded twins, each jogging every panel coordinate by -1, 0 or +1 ULP (x and y independently, splitmix64 of (seed, node)); same OPER script; the largest absolute |base − twin| per value over the set, a flip in any twin",
        seeds.len()
    ));
    merged["twins"] = serde_json::json!(seeds.len());
    merged["seeds"] = serde_json::json!(seeds);
    fs::write(
        work.join("noise_floor.json"),
        serde_json::to_string_pretty(&merged).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    let flips = merged["flips"].as_array().map_or(0, |a| a.len());
    Ok(if flips == 0 {
        format!("noise floor written, branch trace identical in {} twins", seeds.len())
    } else {
        format!("noise floor written, {flips} branch flips over {} twins", seeds.len())
    })
}

/// One twin's spread against the base run (the record `merge_max` combines).
fn twin_record(work: &Path, ulp: &Path) -> Result<serde_json::Value, String> {
    let pa = parse_points(&work.join("viscal_points.dat"))?;
    let pb = parse_points(&ulp.join("viscal_points.dat"))?;
    let ia = parse_iters(&work.join("viscal_iters_all.dat"))?;
    let ib = parse_iters(&ulp.join("viscal_iters_all.dat"))?;

    // call 1's viscal_iter.dat carries UPDATE's RMXBL and the reported limiter VMXBL/IMXBL/ISMXBL
    let qa = parse_blocks(&work.join("viscal_iter.dat"), "ITER").unwrap_or_default();
    let qb = parse_blocks(&ulp.join("viscal_iter.dat"), "ITER").unwrap_or_default();

    let mut flips: Vec<String> = Vec::new();
    if pa.len() != pb.len() {
        flips.push(format!("VISCAL call count {} vs {}", pa.len(), pb.len()));
    }
    let point_keys = [
        "ALFA", "CL", "CM", "CD", "CDF", "CDP", "XOCTR1", "XOCTR2", "MINF", "REINF", "RMSBL",
    ];
    let iter_keys = ["RMSBL", "RLX", "CL", "CD", "CM", "ALFA", "MINF", "REINF"];
    let iter_cols = [1usize, 2, 3, 4, 5, 10, 11, 12];
    let mut calls = Vec::new();
    for (k, (a, b)) in pa.iter().zip(&pb).enumerate() {
        let call = k + 1;
        for key in ["NITDONE", "LVCONV", "IST", "ITRAN1", "ITRAN2"] {
            if a.get(key) != b.get(key) {
                flips.push(format!(
                    "call {call}: {key} {} vs {}",
                    a.get(key).cloned().unwrap_or_default(),
                    b.get(key).cloned().unwrap_or_default()
                ));
            }
        }
        let mut point = serde_json::Map::new();
        for key in point_keys {
            if let (Some(x), Some(y)) = (a.get(key), b.get(key)) {
                let (x, y): (f64, f64) = (x.parse().unwrap_or(0.0), y.parse().unwrap_or(0.0));
                point.insert(key.to_string(), serde_json::json!((x - y).abs()));
            }
        }
        let (ra, rb) = (
            ia.get(&call).cloned().unwrap_or_default(),
            ib.get(&call).cloned().unwrap_or_default(),
        );
        if ra.len() != rb.len() {
            flips.push(format!("call {call}: iteration count {} vs {}", ra.len(), rb.len()));
        }
        let mut iterations = Vec::new();
        for (i, (x, y)) in ra.iter().zip(&rb).enumerate() {
            for (name, col) in [("IST", 7usize), ("ITRAN1", 8), ("ITRAN2", 9)] {
                if x[col] != y[col] {
                    flips.push(format!(
                        "call {call} iteration {}: {name} {} vs {}",
                        i + 1,
                        x[col],
                        y[col]
                    ));
                }
            }
            let mut m = serde_json::Map::new();
            for (key, col) in iter_keys.iter().zip(iter_cols) {
                m.insert(key.to_string(), serde_json::json!((x[col] - y[col]).abs()));
            }
            if call == 1 {
                if let (Some(ba), Some(bb)) = (qa.get(i), qb.get(i)) {
                    let (ra, rb): (f64, f64) = (
                        ba.get("RMXBL").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                        bb.get("RMXBL").and_then(|v| v.parse().ok()).unwrap_or(0.0),
                    );
                    m.insert("RMXBL".to_string(), serde_json::json!((ra - rb).abs()));
                    let same = ba.get("VMXBL") == bb.get("VMXBL")
                        && ba.get("IMXBL") == bb.get("IMXBL")
                        && ba.get("ISMXBL") == bb.get("ISMXBL");
                    m.insert(
                        "LIMITER_FLIP".to_string(),
                        serde_json::json!(if same { 0.0 } else { 1.0 }),
                    );
                }
            }
            iterations.push(serde_json::Value::Object(m));
        }
        let mut record = serde_json::json!({ "call": call, "point": point, "iterations": iterations });
        // the final state's per-array spread (max over nodes / stations of |a - b|), when the
        // case keeps the state records; missing rows count as NaN, which `json_f64` reads back
        if let (Some(sa), Some(sb)) = (
            read_state_arrays(&work.join(format!("viscal_state_{call}.dat"))),
            read_state_arrays(&ulp.join(format!("viscal_state_{call}.dat"))),
        ) {
            let mut spread = serde_json::Map::new();
            for (name, a) in &sa {
                let b = sb.get(name);
                let worst = a
                    .iter()
                    .enumerate()
                    .map(|(i, x)| b.and_then(|b| b.get(i)).map_or(f64::NAN, |y| (x - y).abs()))
                    .fold(0.0_f64, |acc, d| {
                        if d.is_nan() || acc.is_nan() {
                            f64::NAN
                        } else {
                            acc.max(d)
                        }
                    });
                spread.insert(name.clone(), serde_json::json!(worst));
            }
            record["state"] = serde_json::Value::Object(spread);
        }
        calls.push(record);
    }
    // per-station spreads of the dumped post-UPDATE arrays (replay harness): max |base - twin|
    // over stations, per array, keyed by SETBL call number; with UPDATE's own header scalars
    // (RLX, RMSBL, RMXBL, CL, DAC). A SETBL call number counts iterations from the start of
    // the run, so base and twin dumps of the same number are the same (VISCAL call, iteration)
    // only while their iteration counts agree: a dump after a count flip has no twin at the
    // same iteration and records no floor
    let nitdone = |pts: &[std::collections::HashMap<String, String>]| -> Vec<usize> {
        pts.iter()
            .map(|p| p.get("NITDONE").and_then(|v| v.trim().parse().ok()).unwrap_or(0))
            .collect()
    };
    let (na, nb) = (nitdone(&pa), nitdone(&pb));
    let aligned = |k: usize| -> bool {
        // the VISCAL call containing SETBL k in the base run, and the count before it in both
        let mut before = 0usize;
        for (c, n) in na.iter().enumerate() {
            if k <= before + n {
                let twin_before: usize = nb.iter().take(c).sum();
                return twin_before == before;
            }
            before += n;
        }
        false
    };
    let mut arrays = serde_json::Map::new();
    for entry in fs::read_dir(work).map_err(|e| e.to_string())?.flatten() {
        let n = entry.file_name().to_string_lossy().to_string();
        let Some(k) = n.strip_prefix("update_output_").and_then(|r| r.strip_suffix(".dat")) else {
            continue;
        };
        let twin = ulp.join(&n);
        if !twin.exists() || !k.parse::<usize>().is_ok_and(aligned) {
            continue;
        }
        let names = ["XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS"];
        let mut worst = [0.0_f64; 6];
        let rows = |p: &Path| -> Vec<Vec<f64>> {
            fs::read_to_string(p)
                .unwrap_or_default()
                .lines()
                .filter(|l| l.starts_with("BL("))
                .map(|l| {
                    l.split_once(")=")
                        .unwrap()
                        .1
                        .split_whitespace()
                        .map(|t| t.parse().unwrap_or(f64::NAN))
                        .collect()
                })
                .collect()
        };
        for (a, b) in rows(&entry.path()).iter().zip(rows(&twin)) {
            for m in 0..6 {
                worst[m] = worst[m].max((a[m] - b[m]).abs());
            }
        }
        let mut per = serde_json::Map::new();
        for (m, name) in names.iter().enumerate() {
            per.insert(name.to_string(), serde_json::json!(worst[m]));
        }
        let header = |p: &Path| -> std::collections::HashMap<String, f64> {
            fs::read_to_string(p)
                .unwrap_or_default()
                .lines()
                .filter(|l| !l.starts_with("BL(") && !l.starts_with("BLX("))
                .filter_map(|l| l.split_once('='))
                .filter_map(|(key, v)| v.trim().parse::<f64>().ok().map(|x| (key.trim().to_string(), x)))
                .collect()
        };
        let (ha, hb) = (header(&entry.path()), header(&twin));
        for key in ["RLX", "RMSBL", "RMXBL", "CL", "DAC"] {
            if let (Some(x), Some(y)) = (ha.get(key), hb.get(key)) {
                per.insert(key.to_string(), serde_json::json!((x - y).abs()));
            }
        }
        // the transition arc length MRCHDU found in the same SETBL call (`mrchdu_output_<k>.dat`)
        let m = format!("mrchdu_output_{k}.dat");
        if ulp.join(&m).exists() {
            let (ma, mb) = (header(&work.join(&m)), header(&ulp.join(&m)));
            for key in ["XSSITR1", "XSSITR2"] {
                if let (Some(x), Some(y)) = (ma.get(key), mb.get(key)) {
                    per.insert(key.to_string(), serde_json::json!((x - y).abs()));
                }
            }
        }
        arrays.insert(k.to_string(), serde_json::Value::Object(per));
    }
    let branch_identical = flips.is_empty();
    Ok(serde_json::json!({
        "branch_identical": branch_identical,
        "flips": flips,
        "calls": calls,
        "update_output": arrays,
    }))
}

/// Merge two twins' records: numbers by their maximum (a `null` — a non-finite spread — wins),
/// booleans by `and`, the `flips` lists by union in order, objects by key, arrays element by
/// element (an element only one side has is kept).
fn merge_max(a: &serde_json::Value, b: &serde_json::Value) -> serde_json::Value {
    use serde_json::Value;
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            let mut out = serde_json::Map::new();
            for (k, v) in x {
                if k == "flips" {
                    let mut list: Vec<Value> = v.as_array().cloned().unwrap_or_default();
                    for f in y.get(k).and_then(|f| f.as_array()).into_iter().flatten() {
                        if !list.contains(f) {
                            list.push(f.clone());
                        }
                    }
                    out.insert(k.clone(), Value::Array(list));
                } else if let Some(w) = y.get(k) {
                    out.insert(k.clone(), merge_max(v, w));
                } else {
                    out.insert(k.clone(), v.clone());
                }
            }
            for (k, w) in y {
                if !x.contains_key(k) {
                    out.insert(k.clone(), w.clone());
                }
            }
            Value::Object(out)
        }
        (Value::Array(x), Value::Array(y)) => {
            let n = x.len().max(y.len());
            Value::Array(
                (0..n)
                    .map(|i| match (x.get(i), y.get(i)) {
                        (Some(v), Some(w)) => merge_max(v, w),
                        (Some(v), None) => v.clone(),
                        (None, Some(w)) => w.clone(),
                        (None, None) => Value::Null,
                    })
                    .collect(),
            )
        }
        (Value::Number(_), Value::Number(_)) => {
            let (p, q) = (json_f64(a), json_f64(b));
            if p.is_nan() || q.is_nan() {
                Value::Null
            } else {
                serde_json::json!(p.max(q))
            }
        }
        (Value::Bool(p), Value::Bool(q)) => Value::Bool(*p && *q),
        (Value::Null, _) | (_, Value::Null) => Value::Null,
        _ => a.clone(),
    }
}

/// The +1-ULP twin of an inviscid-only case: SPECAL (`specal_points.dat`) and SPECCL
/// (`speccl_points.dat`) blocks compared call by call — the SPECCL exit iteration ITAL is the
/// branch trace, the point values and the per-node GAM/QINV/CPI the spreads.
fn inviscid_twin_summary(work: &Path, ulp: &Path) -> Result<serde_json::Value, String> {
    let mut flips: Vec<String> = Vec::new();
    let mut worst = 0.0_f64;
    let mut out = serde_json::Map::new();
    for (file, kind) in [("specal_points.dat", "specal"), ("speccl_points.dat", "speccl")] {
        let (a, b) = (work.join(file), ulp.join(file));
        if !a.exists() && !b.exists() {
            continue;
        }
        let pa = parse_points(&a)?;
        let pb = parse_points(&b)?;
        if pa.len() != pb.len() {
            flips.push(format!("{kind} call count {} vs {}", pa.len(), pb.len()));
        }
        let mut calls = Vec::new();
        for (k, (x, y)) in pa.iter().zip(&pb).enumerate() {
            let call = k + 1;
            if x.get("ITAL") != y.get("ITAL") {
                flips.push(format!(
                    "{kind} call {call}: ITAL {} vs {}",
                    x.get("ITAL").cloned().unwrap_or_default(),
                    y.get("ITAL").cloned().unwrap_or_default()
                ));
            }
            let mut point = serde_json::Map::new();
            for key in ["ALFA", "CL", "CM", "CDP", "MINF", "REINF"] {
                if let (Some(u), Some(v)) = (x.get(key), y.get(key)) {
                    let d = (u.parse::<f64>().unwrap_or(0.0) - v.parse::<f64>().unwrap_or(0.0)).abs();
                    worst = worst.max(d);
                    point.insert(key.to_string(), serde_json::json!(d));
                }
            }
            // per-node spreads (`NODE(i)= GAM QINV CPI`)
            let mut nodes = [0.0_f64; 3];
            for (key, u) in x.iter().filter(|(k, _)| k.starts_with("NODE")) {
                let Some(v) = y.get(key) else { continue };
                let ru: Vec<f64> = u.split_whitespace().filter_map(|t| t.parse().ok()).collect();
                let rv: Vec<f64> = v.split_whitespace().filter_map(|t| t.parse().ok()).collect();
                for m in 0..3 {
                    if let (Some(p), Some(q)) = (ru.get(m), rv.get(m)) {
                        nodes[m] = nodes[m].max((p - q).abs());
                    }
                }
            }
            for (m, name) in ["GAM", "QINV", "CPI"].iter().enumerate() {
                worst = worst.max(nodes[m]);
                point.insert(name.to_string(), serde_json::json!(nodes[m]));
            }
            calls.push(serde_json::json!({ "call": call, "point": point }));
        }
        out.insert(format!("{kind}_calls"), serde_json::Value::Array(calls));
    }
    let branch_identical = flips.is_empty();
    let mut j = serde_json::json!({
        "branch_identical": branch_identical,
        "flips": flips,
        "calls": [],
        "update_output": {},
    });
    for (k, v) in out {
        j[k] = v;
    }
    let _ = worst;
    Ok(j)
}

/// The arrays of a `viscal_state_<k>.dat` record, by name: the node arrays (`CPI`, `CPV`,
/// `QINV`, `QVIS`, `GAM`, one entry per node in file order) and the station arrays
/// (`XSSI`, `X`, `UEDG`, `THET`, `DSTR`, `CTAU`, `MASS`, `TAU`, one entry per `BL(is, ibl)` row
/// in file order — side 1 then side 2 with its wake). `None` when the file is absent.
fn read_state_arrays(path: &Path) -> Option<std::collections::BTreeMap<String, Vec<f64>>> {
    const NODE: [&str; 5] = ["CPI", "CPV", "QINV", "QVIS", "GAM"];
    const BL: [&str; 8] = ["XSSI", "X", "UEDG", "THET", "DSTR", "CTAU", "MASS", "TAU"];
    let text = fs::read_to_string(path).ok()?;
    let mut arrays: std::collections::BTreeMap<String, Vec<f64>> = Default::default();
    for l in text.lines() {
        let (names, rest): (&[&str], &str) = if l.starts_with("NODE(") {
            (&NODE, l.split_once(")=")?.1)
        } else if l.starts_with("BL(") {
            (&BL, l.split_once(")=")?.1)
        } else {
            continue;
        };
        let v: Vec<f64> = rest.split_whitespace().map(|t| t.parse().unwrap_or(f64::NAN)).collect();
        for (m, name) in names.iter().enumerate() {
            arrays
                .entry(name.to_string())
                .or_default()
                .push(v.get(m).copied().unwrap_or(f64::NAN));
        }
    }
    Some(arrays)
}

fn parse_blocks(p: &Path, start_key: &str) -> Result<Vec<std::collections::HashMap<String, String>>, String> {
    let text = fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    let mut out: Vec<std::collections::HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == start_key {
            out.push(Default::default());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    Ok(out)
}

fn parse_points(p: &Path) -> Result<Vec<std::collections::HashMap<String, String>>, String> {
    let text = fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    let mut out: Vec<std::collections::HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "CALL" {
            out.push(Default::default());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    Ok(out)
}

fn parse_iters(p: &Path) -> Result<std::collections::HashMap<usize, Vec<Vec<f64>>>, String> {
    let text = fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    let mut out: std::collections::HashMap<usize, Vec<Vec<f64>>> = Default::default();
    for l in text.lines() {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(|t| t.parse().unwrap_or(f64::NAN)).collect();
        out.entry(v[0] as usize).or_default().push(v[1..].to_vec());
    }
    Ok(out)
}

/// A spread from `noise_floor.json`: serde_json writes a non-finite f64 (a NaN reference value,
/// as a diverged run leaves them) as `null`, read back as NaN.
fn json_f64(v: &serde_json::Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}
