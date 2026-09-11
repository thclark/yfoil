//! Aerofoil series study (`docs/validation/aerofoil-series/`).
//!
//! Draws every generator family yFoil offers — NACA 4-digit, 4-digit modified, 5-digit,
//! 16-series, 6-series, 6A-series and the Kármán–Trefftz analytic section — as one figure per
//! family with one panel per parameter: the family default (the validation set's member) in
//! black on top, the members that vary the parameter behind it coloured by value (YlGnBu), all
//! panels on the same axes (x −0.2…1.2, y −0.5…0.5, equal scale). It also compares the
//! generators with the public-domain NASA/PDAS `naca456` reference fixtures
//! (`tests/fixtures/naca456/`) and tabulates the worst difference per case.
//!
//! Every invocation creates a dated folder `runs/<UTC datetime>/` next to this crate holding
//! `metadata.json` (git version tag or `untagged`, commit sha, dirty flag, the naca456 fixture
//! manifest), one geometry JSON per drawn section under `<series>/<panel>/<slug>.json` (each with
//! its `generator` provenance record), `summary.json` (every section, panel and comparison
//! number) and the `README.md` / `README.tex` index; the figures (SVG and PDF) are then drawn from
//! those files by `plot.py` through `scripts/figures/render.sh` (matplotlib; presentation only,
//! every number comes from the JSON). `runs/` is gitignored. With `--docs` the figures and a
//! generated Markdown page are also written to `docs/validation/aerofoil-series/`, the only
//! place a run's products are tracked.
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo run --release -p aerofoil-series -- [--docs] [--series naca-6]
//! ```

mod naca456;
mod report;
mod study;

use std::path::{Path, PathBuf};

/// Axis extents of every panel, chord units
pub const X_RANGE: (f64, f64) = (-0.2, 1.2);
pub const Y_RANGE: (f64, f64) = (-0.5, 0.5);

/// Draw a run folder's figures with `scripts/figures/render.sh <study> <run_dir>` (matplotlib
/// through uv); on failure the data are complete and the command is printed for a later run.
pub fn render(study: &str, run_dir: &Path) -> bool {
    let root = repo_root();
    let status = std::process::Command::new(root.join("scripts/figures/render.sh"))
        .arg(study)
        .arg(run_dir)
        .status();
    match status {
        Ok(st) if st.success() => true,
        other => {
            eprintln!(
                "figures not drawn ({}): run `scripts/figures/render.sh {study} {}`",
                other.map(|s| s.to_string()).unwrap_or_else(|e| e.to_string()),
                run_dir.display()
            );
            false
        }
    }
}

/// Repository root: this crate lives at scripts/aerofoil-series
pub fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("repository root")
        .to_path_buf()
}

/// `runs/` next to this crate
pub fn runs_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("runs")
}

fn git(root: &Path, args: &[&str]) -> Option<String> {
    let out = std::process::Command::new("git")
        .args(args)
        .current_dir(root)
        .output()
        .ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// UTC timestamp `YYYY-MM-DDTHH-MM-SSZ` (filesystem-safe), from the Unix time
fn utc_stamp() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs() as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    // civil-from-days (Howard Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    format!("{y:04}-{mo:02}-{d:02}T{h:02}-{m:02}-{s:02}Z")
}

/// `metadata.json` of a run folder
pub fn write_metadata(run_dir: &Path, series: &[&str]) {
    let root = repo_root();
    let sha = git(&root, &["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".into());
    let version = git(&root, &["describe", "--exact-match", "--tags", "HEAD"]).unwrap_or_else(|| "untagged".into());
    let dirty = git(&root, &["status", "--porcelain"])
        .map(|s| !s.is_empty())
        .unwrap_or(true);
    let naca456_manifest: serde_json::Value =
        std::fs::read_to_string(root.join("tests/fixtures/naca456/manifest.json"))
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or(serde_json::Value::Null);
    let j = serde_json::json!({
        "study": "aerofoil-series",
        "created_utc": run_dir.file_name().unwrap().to_string_lossy(),
        "version": version,
        "commit": sha,
        "dirty": dirty,
        "yfoil": env!("CARGO_PKG_VERSION"),
        "series": series,
        "n_panels": study::N_PANELS,
        "axes": { "x": X_RANGE, "y": Y_RANGE },
        "naca456_fixtures": naca456_manifest,
    });
    std::fs::write(run_dir.join("metadata.json"), serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

struct Args {
    docs: bool,
    only: Option<String>,
}

fn parse_args() -> Args {
    let mut a = Args {
        docs: false,
        only: None,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--docs" => a.docs = true,
            "--series" => {
                a.only = Some(argv[i + 1].clone());
                i += 1;
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    a
}

fn main() {
    let args = parse_args();
    let run_dir = runs_root().join(utc_stamp());
    std::fs::create_dir_all(&run_dir).unwrap();

    let all = study::series();
    let selected: Vec<&study::SeriesStudy> = all
        .iter()
        .filter(|s| args.only.as_deref().is_none_or(|o| o == s.slug))
        .collect();
    assert!(!selected.is_empty(), "no series named {:?}", args.only);
    let slugs: Vec<&str> = selected.iter().map(|s| s.slug).collect();
    write_metadata(&run_dir, &slugs);

    for s in &selected {
        // the geometries, one file per drawn section
        for p in &s.panels {
            let dir = run_dir.join(s.slug).join(p.slug);
            std::fs::create_dir_all(&dir).unwrap();
            for m in &p.members {
                yfoil::geometry::write_geometry_to_json(&m.geometry, dir.join(format!("{}.json", m.slug))).unwrap();
            }
        }
        yfoil::geometry::write_geometry_to_json(
            &s.default.geometry,
            run_dir.join(s.slug).join(format!("{}.json", s.default.slug)),
        )
        .unwrap();
    }

    let comparison = naca456::compare(&repo_root());
    report::write_summary_json(&run_dir, &selected, &comparison);
    report::write_index(&run_dir, &selected, &comparison);
    report::write_index_tex(&run_dir, &selected, &comparison);
    let drawn = render("aerofoil-series", &run_dir);
    if args.docs {
        assert!(drawn, "--docs needs the figures");
        let docs = repo_root().join("docs/validation/aerofoil-series");
        report::write_docs(&docs, &run_dir, &all, &comparison);
        println!("wrote {}", docs.display());
    }
    println!("run folder: {}", run_dir.display());
}
