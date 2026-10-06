//! Closure fixtures: `tests/fixtures/subroutines/<closure>/<case>_NNN.json`.
//!
//! The instrumented reference logs the inputs and outputs of every call of the BL closure
//! functions into `xfoil_subroutine_log.dat` (`=== HKIN ===` blocks of `KEY= value` lines; patch
//! 16). A case with `closures = true` is run like any other — yFoil's panels, `LOAD`ed — and its
//! log is cut into one JSON file per distinct input set (`{"input": {...}, "output": {...}}`, keys
//! the lower-cased Fortran names), at most `CAP` per closure and case. A case's files carry its
//! name as their prefix, so several cases (incompressible and compressible) share a closure's
//! directory; `tests/fixtures/subroutines/manifest.json` records, per case, the reference build
//! that produced them. The gates are `tests/subroutine/closures.rs` and `axset.rs`.

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

/// At most this many distinct calls per closure and case.
const CAP: usize = 100;

/// One closure: its log section, the directory its fixtures go to, and its input and output
/// keys as the log writes them.
struct Closure {
    section: &'static str,
    dir: &'static str,
    inputs: &'static [&'static str],
    outputs: &'static [&'static str],
}

const CLOSURES: &[Closure] = &[
    Closure {
        section: "HKIN",
        dir: "hkin",
        inputs: &["H", "MSQ"],
        outputs: &["HK", "HK_H", "HK_MSQ"],
    },
    Closure {
        section: "CFL",
        dir: "cfl",
        inputs: &["HK", "RT", "MSQ"],
        outputs: &["CF", "CF_HK", "CF_RT", "CF_MSQ"],
    },
    Closure {
        section: "HSL",
        dir: "hsl",
        inputs: &["HK", "RT", "MSQ"],
        outputs: &["HS", "HS_HK", "HS_RT", "HS_MSQ"],
    },
    Closure {
        section: "DIL",
        dir: "dil",
        inputs: &["HK", "RT"],
        outputs: &["DI", "DI_HK", "DI_RT"],
    },
    Closure {
        section: "HST",
        dir: "hst",
        inputs: &["HK", "RT", "MSQ"],
        outputs: &["HS", "HS_HK", "HS_RT", "HS_MSQ"],
    },
    Closure {
        section: "CFT",
        dir: "cft",
        inputs: &["HK", "RT", "MSQ"],
        outputs: &["CF", "CF_HK", "CF_RT", "CF_MSQ"],
    },
    Closure {
        section: "DAMPL",
        dir: "dampl",
        inputs: &["HK", "TH", "RT"],
        outputs: &["AX", "AX_HK", "AX_TH", "AX_RT"],
    },
    Closure {
        section: "AXSET",
        dir: "axset",
        inputs: &["HK1", "T1", "RT1", "A1", "HK2", "T2", "RT2", "A2", "ACRIT"],
        outputs: &[
            "AX", "AX_HK1", "AX_T1", "AX_RT1", "AX_A1", "AX_HK2", "AX_T2", "AX_RT2", "AX_A2",
        ],
    },
    Closure {
        section: "BLKIN",
        dir: "blkin",
        inputs: &["T2", "D2", "U2", "HSTINV", "GM1BL", "RSTBL", "HVRAT", "REYBL"],
        outputs: &[
            "M2", "H2", "HK2", "RT2", "HK2_T2", "HK2_D2", "HK2_U2", "RT2_T2", "RT2_U2",
        ],
    },
    Closure {
        section: "BLVAR",
        dir: "blvar",
        inputs: &["ITYP", "HK2", "RT2", "M2", "T2", "D2", "S2"],
        outputs: &["HS2", "CF2", "DI2", "US2", "CQ2", "DE2"],
    },
];

/// Cut case `case`'s log (`work/xfoil_subroutine_log.dat`) into fixtures under `out`
/// (`tests/fixtures/subroutines`), replacing that case's earlier files, and record the case in
/// `out/manifest.json` with the reference build's manifest lines.
pub fn write(case: &str, work: &Path, out: &Path, reference: &str) -> Result<String, String> {
    let log = fs::read_to_string(work.join("xfoil_subroutine_log.dat"))
        .map_err(|e| format!("xfoil_subroutine_log.dat: {e}"))?;
    let mut blocks: Vec<(String, BTreeMap<String, String>)> = vec![];
    for line in log.lines() {
        let t = line.trim();
        if let Some(name) = t.strip_prefix("=== ").and_then(|r| r.strip_suffix(" ===")) {
            blocks.push((name.to_string(), BTreeMap::new()));
        } else if let (Some((_, cur)), Some((k, v))) = (blocks.last_mut(), t.split_once('=')) {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    let mut counts = vec![];
    for c in CLOSURES {
        let dir = out.join(c.dir);
        fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
        for e in fs::read_dir(&dir).map_err(|e| e.to_string())?.flatten() {
            if e.file_name().to_string_lossy().starts_with(&format!("{case}_")) {
                fs::remove_file(e.path()).map_err(|e| e.to_string())?;
            }
        }
        let mut seen: HashSet<Vec<String>> = HashSet::new();
        let mut n = 0;
        for (name, kv) in &blocks {
            if name != c.section || n >= CAP {
                continue;
            }
            let Some(key) = c.inputs.iter().map(|k| kv.get(*k).cloned()).collect::<Option<Vec<_>>>() else {
                continue;
            };
            if !c.outputs.iter().all(|k| kv.contains_key(*k)) || !seen.insert(key) {
                continue;
            }
            let num = |k: &str| -> serde_json::Value {
                let v = fortran_real(&kv[k]);
                if k == "ITYP" {
                    serde_json::json!(v as i64)
                } else {
                    serde_json::json!(v)
                }
            };
            let obj = |keys: &[&str]| -> serde_json::Value {
                serde_json::Value::Object(keys.iter().map(|k| (k.to_lowercase(), num(k))).collect())
            };
            n += 1;
            let j = serde_json::json!({ "input": obj(c.inputs), "output": obj(c.outputs) });
            fs::write(
                dir.join(format!("{case}_{n:03}.json")),
                serde_json::to_string_pretty(&j).unwrap(),
            )
            .map_err(|e| e.to_string())?;
        }
        counts.push(format!("{} {n}", c.dir));
    }
    let mpath = out.join("manifest.json");
    let mut manifest: serde_json::Value = fs::read_to_string(&mpath)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_else(|| serde_json::json!({ "generated_by": "cargo xtask fixtures", "cases": {} }));
    manifest["cases"][case] = serde_json::json!({ "xfoil_ref": reference.lines().collect::<Vec<_>>() });
    fs::write(&mpath, serde_json::to_string_pretty(&manifest).unwrap()).map_err(|e| e.to_string())?;
    Ok(counts.join(", "))
}

/// A Fortran `ES24.16` real. With a three-digit exponent the format has no room for the `E` and
/// writes `3.0480272680724994-314`; that form is read too.
fn fortran_real(s: &str) -> f64 {
    let s = s.trim();
    if let Ok(v) = s.parse() {
        return v;
    }
    match s[1..].rfind(['+', '-']) {
        // the sign at s[j + 1] starts the exponent
        Some(j) => format!("{}e{}", &s[..=j], &s[j + 1..]).parse().unwrap_or(f64::NAN),
        None => f64::NAN,
    }
}
