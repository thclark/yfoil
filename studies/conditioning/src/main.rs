//! Conditioning: the thresholds for yFoil's own twins (`yfoil::output::Thresholds`).
//!
//! The series cases (`group = "series"`: the aerofoil series and flow regimes of interest under
//! normal conditions) and the branch-coverage cases (`group = "branch-coverage"`: the operating
//! points that take XFOIL's rarely-taken branches) of `xtask/fixtures-config/cases.toml` are used
//! as a set of yFoil operating points: each case's panels and its OPER procedure (`ALFA` points,
//! `ASEQ` with NSEQEX halting, `INIT`, `CL` points; viscous or inviscid; `TYPE`, `XTR`, `DAMP`),
//! driven through yFoil with the solver's per-point twins. A case whose procedure uses a setting
//! the driver does not express ([`drive::UNEXPRESSED`]) is not driven, and is listed as such. Nothing the reference
//! computed is read. For each point valid on single-run evidence the study takes the largest
//! spread over the reported quantities and the five twins (`|a − b| / max(|a|, |b|, 1)`), and
//! shows the two populations — the noise floor and the points whose twins reach another
//! solution — and the gap between them, inside which the error threshold sits and across which
//! the set of withheld points does not change; and where in the noise floor the warning
//! threshold sits.
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo run --release -p conditioning
//! cargo run --release -p conditioning -- --render studies/conditioning/runs/<UTC>
//! cargo xtask publish conditioning [--run <UTC>]
//! ```
//!
//! A run writes `runs/<UTC>/` (`metadata.json`, `summary.json`, `README.md`, the figures with
//! their PDFs, `publish/`), never under `docs/`.

mod analysis;
mod drive;
mod figures;
mod report;
mod utilities;

use serde_json::{json, Value};
use std::path::{Path, PathBuf};

fn rerender(run_dir: &Path) {
    let s: Value = serde_json::from_str(&std::fs::read_to_string(run_dir.join("summary.json")).expect("summary.json"))
        .expect("parse summary.json");
    figures::render(run_dir, &s);
    report::write_readme(run_dir, &s);
    report::write_publish(run_dir, &s);
}

/// The groups of cases the study drives, in the order they are reported.
pub const GROUPS: [&str; 2] = ["series", "branch-coverage"];

/// The raw `[[case]]` tables of cases.toml, by name: for the settings `series_cases::Case` does
/// not carry.
fn raw_cases(root: &Path) -> std::collections::HashMap<String, toml::Value> {
    let t: toml::Value =
        toml::from_str(&std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("parse cases.toml");
    t["case"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| (c["name"].as_str().unwrap().to_string(), c.clone()))
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Some(i) = args.iter().position(|a| a == "--render") {
        rerender(&args.get(i + 1).map(PathBuf::from).expect("--render <run_dir>"));
        return;
    }
    let root = series_cases::repo_root();
    let cases = series_cases::load_cases(&root);
    let raw = raw_cases(&root);
    let mut records = vec![];
    // cases of the study's groups that were not driven, with why: reported, never silently dropped
    let mut not_driven: Vec<Value> = vec![];
    for case in GROUPS
        .iter()
        .flat_map(|g| cases.iter().filter(move |c| c.group.as_deref() == Some(*g)))
    {
        let name = &case.name;
        let group = case.group.clone().unwrap_or_default();
        let unexpressed: Vec<&str> = drive::UNEXPRESSED
            .iter()
            .copied()
            .filter(|k| raw[name].get(*k).is_some())
            .collect();
        if !unexpressed.is_empty() {
            let why = format!(
                "its procedure sets {}, which the driver does not express",
                unexpressed.join(", ")
            );
            eprintln!("  {name}: NOT DRIVEN — {why}");
            not_driven.push(json!({"name": name, "group": group, "why": why}));
            continue;
        }
        let Some(panels) = drive::panels(&root, name) else {
            let why = format!("no panels.json (run `cargo xtask fixtures --case {name}`)");
            eprintln!("  {name}: NOT DRIVEN — {why}");
            not_driven.push(json!({"name": name, "group": group, "why": why}));
            continue;
        };
        let start = std::time::Instant::now();
        let r = drive::case(&panels, case);
        println!(
            "{name}: {} points ({:.1} s)",
            r["points"].as_array().unwrap().len(),
            start.elapsed().as_secs_f64()
        );
        records.push(r);
    }
    assert!(!records.is_empty(), "no case of the study's groups could be driven");

    let stamp = {
        let out = std::process::Command::new("date")
            .args(["-u", "+%Y-%m-%dT%H-%M-%SZ"])
            .output()
            .expect("date");
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };
    let run_dir = root.join("studies/conditioning/runs").join(&stamp);
    std::fs::create_dir_all(&run_dir).unwrap();
    let metadata = json!({
        "study": "conditioning",
        "generated": stamp,
        "git_describe": series_cases::git(&root, &["describe", "--tags", "--always"]).unwrap_or_else(|| "untagged".into()),
        "git_commit": series_cases::git(&root, &["rev-parse", "HEAD"]).unwrap_or_default(),
        "git_dirty": series_cases::git(&root, &["status", "--porcelain"]).map(|s| !s.is_empty()).unwrap_or(true),
        "host": format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        "nseqex": drive::NSEQEX,
    });
    std::fs::write(
        run_dir.join("metadata.json"),
        serde_json::to_string_pretty(&metadata).unwrap(),
    )
    .unwrap();
    let summary = json!({
        "thresholds": yfoil::output::Thresholds::default(),
        "seeds": yfoil::geometry::jitter::TWIN_SEEDS,
        "groups": GROUPS,
        "cases": records,
        "not_driven": not_driven,
    });
    std::fs::write(
        run_dir.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .unwrap();
    rerender(&run_dir);
    println!("run folder: {}", run_dir.display());
    println!("publish it with: cargo xtask publish conditioning --run {stamp}");
}
