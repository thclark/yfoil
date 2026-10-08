//! Branch-coverage study (`docs/validation/branch-coverage/`).
//!
//! CLAUDE.md Rule 6: "numerically exact" is claimable only over branches that have been
//! exercised. `cargo xtask coverage` measures which branches of the translated XFOIL subroutines
//! the reference takes over *all* tracked cases; this study asks the sharper question — what is
//! the **minimum set of cases** that takes every branch of the *analysis path* (the panelling
//! generators and their spline/geometry helpers are excluded: they are gated by the `pangen_*`
//! cases, `docs/validation/geometry/`) — and then gates every case of that set against its
//! instrumented reference, so that the coverage claim and the equivalence claim are made over
//! the same runs.
//!
//! What it does, in order:
//!
//! 1. `cargo xtask coverage --per-case --with-group non-finite` (unless `--no-measure`): the gcov
//!    reference is run on every tracked case, and on the untracked `non-finite` probes if their
//!    fixtures exist (`cargo xtask fixtures --group non-finite`), separately, writing `target/coverage/branches.json` (the universe: every
//!    branch with its status over all cases) and `target/coverage/per-case/<case>/branches.json`
//!    (the branches that case took, with the reference's stdout beside it).
//! 2. Candidates are the tracked solver cases whose reference run is finite — a case whose
//!    reference prints `NaN` is excluded (a non-finite trajectory is not reproducible by
//!    translation: NaN comparisons and MIN/MAX are where two correct codes may legitimately
//!    differ) and reported.
//! 3. The minimum cover of the finite-reachable analysis branches: branches taken by exactly
//!    one candidate make that candidate *essential*; the rest is a greedy cover by marginal gain,
//!    pruned until no case can be removed. The result must equal the `group = "branch-coverage"`
//!    set of `xtask/fixtures-config/cases.toml`, which is how the group is kept honest.
//! 4. Every case of the cover — and, as a measurement rather than a claim, every non-finite
//!    probe — is run through yFoil's `Session` and compared with its fixture
//!    (`scripts/study-support/records.rs`, the studies' floor comparison; the tests compare single steps instead, `docs/conventions/testing.md`): iteration counts and
//!    branch trace exact, values within tolerance or the reference's own 1-ULP spread ×
//!    `FLOOR_FACTOR`, the third outcome divergent.
//!
//! The same sections swept as full ±30° polars, through stall and into the non-finite region, are
//! the `series-cases` study (`scripts/series-cases`).
//!
//! Every invocation creates `runs/<UTC datetime>/` next to this crate with `metadata.json`,
//! `summary.json` (every number of the tables and figures), `README.md` / `README.tex`; the
//! figures are drawn from `summary.json` by `plot.py` through `scripts/figures/render.sh`
//! (matplotlib, presentation only). With `--docs` the Markdown page and the SVG figures are also
//! written to `docs/validation/branch-coverage/`, the only tracked products.
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo run --release -p branch-coverage -- [--docs] [--no-measure]
//! ```

mod cover;
mod drive;
mod report;

mod utilities;

use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// Files and subroutines of the translated set that belong to panelling / geometry generation,
/// not to the analysis path (covered by the `pangen_*` and `naca456` gates).
pub const GEOMETRY_FILES: &[&str] = &["naca.f", "xgeom.f", "spline.f"];
pub const GEOMETRY_SUBROUTINES: &[&str] = &["xfoil.f:NACA", "xfoil.f:PANGEN"];

/// One case as `xtask/fixtures-config/cases.toml` records it (the fields the study needs).
#[derive(Debug, Deserialize, Clone)]
pub struct Case {
    pub name: String,
    pub foil: String,
    pub n_nodes: usize,
    #[serde(default)]
    pub alphas: Vec<f64>,
    #[serde(default)]
    pub alphas_after_reinit: Vec<f64>,
    #[serde(default)]
    pub polar: bool,
    #[serde(default)]
    pub cls: Vec<f64>,
    #[serde(default)]
    pub matyp: usize,
    #[serde(default)]
    pub xtr: Vec<f64>,
    #[serde(default)]
    pub damp: bool,
    #[serde(default)]
    pub geometry_only: bool,
    #[serde(default)]
    pub tgap: Vec<f64>,
    #[serde(default)]
    pub inviscid: bool,
    #[serde(default)]
    pub group: Option<String>,
    #[serde(default)]
    pub re: f64,
    #[serde(default)]
    pub mach: f64,
    #[serde(default = "default_ncrit")]
    pub ncrit: f64,
    #[serde(default = "default_iter")]
    pub max_iterations: usize,
    #[serde(default)]
    pub track: bool,
}
fn default_ncrit() -> f64 {
    9.0
}
fn default_iter() -> usize {
    20
}

#[derive(Deserialize)]
struct Cases {
    #[serde(rename = "case")]
    cases: Vec<Case>,
}

impl Case {
    /// A solver case with a tracked fixture (not a geometry-only or TGAP case). The untracked
    /// `non-finite` probes are measured too, but never candidates of the cover.
    pub fn is_solver_case(&self) -> bool {
        self.track && !self.geometry_only && self.tgap.is_empty()
    }
    /// The OPER script in one line, for the tables.
    pub fn script(&self) -> String {
        let mut parts: Vec<String> = vec![];
        if self.inviscid {
            parts.push(format!("M {}", self.mach));
        } else {
            parts.push(format!("Re {:.0e}", self.re));
            parts.push(format!("M {}", self.mach));
            parts.push(format!("Ncrit {}", self.ncrit));
            parts.push(format!("ITER {}", self.max_iterations));
        }
        if self.matyp != 0 {
            parts.push(format!("TYPE {}", self.matyp));
        }
        if !self.xtr.is_empty() {
            parts.push(format!("XTR {} {}", self.xtr[0], self.xtr[1]));
        }
        if self.damp {
            parts.push("DAMP".into());
        }
        let seq = |a: &[f64]| -> String {
            if self.polar && a.len() > 1 {
                format!("ALFA {} / ASEQ {} {} {}", a[0], a[1], a[a.len() - 1], a[1] - a[0])
            } else {
                a.iter().map(|x| format!("ALFA {x}")).collect::<Vec<_>>().join(" / ")
            }
        };
        if !self.alphas.is_empty() {
            parts.push(seq(&self.alphas));
        }
        if !self.alphas_after_reinit.is_empty() {
            parts.push(format!("INIT / {}", seq(&self.alphas_after_reinit)));
        }
        for c in &self.cls {
            parts.push(format!("CL {c}"));
        }
        parts.join("; ")
    }
    /// A short section name for the tables (`naca4:0012:sharp` → `NACA 0012 (sharp TE)`).
    pub fn section(&self) -> String {
        let mut it = self.foil.split(':');
        let kind = it.next().unwrap_or("");
        let spec = it.next().unwrap_or("");
        let sharp = it.next() == Some("sharp");
        let name = match kind {
            "karman-trefftz" => "Kármán–Trefftz".to_string(),
            _ => format!("NACA {spec}"),
        };
        if sharp {
            format!("{name} (sharp TE)")
        } else {
            name
        }
    }
}

/// The universe file `target/coverage/branches.json`.
#[derive(Deserialize, Clone)]
pub struct Universe {
    pub cases: Vec<String>,
    pub subroutines: Vec<SubroutineRow>,
    pub branches: Vec<Branch>,
}
#[derive(Deserialize, Clone)]
pub struct SubroutineRow {
    pub file: String,
    pub name: String,
    pub calls: u64,
    pub branches: usize,
    pub taken: usize,
}
#[derive(Deserialize, Clone)]
pub struct Branch {
    pub id: String,
    pub file: String,
    pub subroutine: String,
    pub line: usize,
    pub branch: usize,
    pub source: String,
    /// taken | open | loop_entry | unreachable (over all cases of the measurement)
    pub status: String,
    #[serde(default)]
    pub class: String,
    #[serde(default)]
    pub why: Option<String>,
    /// the test gating the branch as an event (`tests/…::test_…`), for non-finite-only branches
    #[serde(default)]
    pub gate: Option<String>,
}
impl Branch {
    pub fn key(&self) -> String {
        format!("{}:{}", self.file, self.subroutine)
    }
    pub fn is_analysis(&self) -> bool {
        !GEOMETRY_FILES.contains(&self.file.as_str()) && !GEOMETRY_SUBROUTINES.contains(&self.key().as_str())
    }
}

/// One candidate: a tracked solver case with its measured branch set.
impl Candidate {
    /// A sweep of another study (`series-cases`) or a `pathological` copy
    /// sets): measured, never a cover candidate.
    pub fn polar_study(&self) -> bool {
        matches!(self.case.group.as_deref(), Some("series") | Some("pathological"))
    }
}

pub struct Candidate {
    pub case: Case,
    /// branch ids of the analysis path this case took
    pub taken: BTreeSet<String>,
    /// per-subroutine call counts (`file:SUB` → calls)
    pub calls: BTreeMap<String, u64>,
    /// per DO line (`file:line`): the counts of its loop-control edges and of the first
    /// executable line of its body on this case
    pub do_loops: BTreeMap<String, (Vec<u64>, u64)>,
    /// the reference printed NaN on this case (its stdout under target/coverage/per-case/<case>/)
    pub non_finite: bool,
    /// number of `NaN` occurrences in that stdout
    pub nan_count: usize,
}

pub fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("repository root")
        .to_path_buf()
}

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

/// Draw the run folder's figures with `scripts/figures/render.sh branch-coverage <run_dir>`.
pub fn render(run_dir: &Path) -> bool {
    let root = repo_root();
    let status = std::process::Command::new(root.join("scripts/figures/render.sh"))
        .arg("branch-coverage")
        .arg(run_dir)
        .status();
    match status {
        Ok(st) if st.success() => true,
        other => {
            eprintln!(
                "figures not drawn ({}): run `scripts/figures/render.sh branch-coverage {}`",
                other.map(|s| s.to_string()).unwrap_or_else(|e| e.to_string()),
                run_dir.display()
            );
            false
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let docs = args.iter().any(|a| a == "--docs");
    let measure = !args.iter().any(|a| a == "--no-measure");
    let root = repo_root();

    // 1. the per-case measurement
    if measure {
        // every tracked case plus the untracked non-finite probes (`cargo xtask fixtures --group
        // non-finite` first), so that the branches only they reach are reported as such
        println!("measuring: cargo xtask coverage --per-case --with-group non-finite");
        let st = std::process::Command::new("cargo")
            .args(["xtask", "coverage", "--per-case", "--with-group", "non-finite"])
            .current_dir(&root)
            .status()
            .expect("cargo xtask coverage");
        assert!(st.success(), "cargo xtask coverage --per-case failed");
    }
    let universe: Universe = serde_json::from_str(
        &std::fs::read_to_string(root.join("target/coverage/branches.json"))
            .expect("target/coverage/branches.json — run without --no-measure"),
    )
    .expect("parse branches.json");
    let cases: Cases =
        toml::from_str(&std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("parse cases.toml");

    // 2. candidates
    let analysis_ids: BTreeSet<&str> = universe
        .branches
        .iter()
        .filter(|b| b.is_analysis())
        .map(|b| b.id.as_str())
        .collect();
    let mut candidates: Vec<Candidate> = vec![];
    for case in cases
        .cases
        .iter()
        .filter(|c| c.is_solver_case() || c.group.as_deref() == Some("non-finite"))
    {
        let dir = root.join("target/coverage/per-case").join(&case.name);
        let Ok(text) = std::fs::read_to_string(dir.join("branches.json")) else {
            eprintln!(
                "  {}: no per-case measurement (not in the universe run) — skipped",
                case.name
            );
            continue;
        };
        let j: serde_json::Value = serde_json::from_str(&text).unwrap();
        let taken: BTreeSet<String> = j["taken"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .filter(|id| analysis_ids.contains(id.as_str()))
            .collect();
        let calls: BTreeMap<String, u64> = j["subroutine_calls"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, v)| (k.clone(), v.as_u64().unwrap_or(0)))
            .collect();
        let do_loops: BTreeMap<String, (Vec<u64>, u64)> = j["do_loops"]
            .as_object()
            .map(|o| {
                o.iter()
                    .map(|(k, v)| {
                        (
                            k.clone(),
                            (
                                v["edges"]
                                    .as_array()
                                    .unwrap()
                                    .iter()
                                    .map(|x| x.as_u64().unwrap_or(0))
                                    .collect(),
                                v["body"].as_u64().unwrap_or(0),
                            ),
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        let stdout = std::fs::read(dir.join("stdout.txt")).unwrap_or_default();
        let nan_count = String::from_utf8_lossy(&stdout).matches("NaN").count();
        candidates.push(Candidate {
            case: case.clone(),
            taken,
            calls,
            do_loops,
            non_finite: nan_count > 0,
            nan_count,
        });
    }
    // the cover is a set of single operating points (and short sequences); the ±30° polars are
    // the series-cases study's cases, built from this cover's sections, and are measured
    // here but never chosen for it
    let finite: Vec<&Candidate> = candidates
        .iter()
        .filter(|c| !c.non_finite && !c.polar_study())
        .collect();
    println!(
        "{} candidate cases, {} eligible for the cover ({} excluded for NaN, {} polar-study sweeps)",
        candidates.len(),
        finite.len(),
        candidates.iter().filter(|c| c.non_finite).count(),
        candidates.iter().filter(|c| c.polar_study() && !c.non_finite).count()
    );

    // 3. the minimum cover
    let target: BTreeSet<String> = finite.iter().flat_map(|c| c.taken.iter().cloned()).collect();
    let cover = cover::minimum_cover(&finite, &target);
    println!("cover: {} cases — {}", cover.len(), cover.join(", "));
    let group: BTreeSet<String> = cases
        .cases
        .iter()
        .filter(|c| c.group.as_deref() == Some("branch-coverage"))
        .map(|c| c.name.clone())
        .collect();
    let computed: BTreeSet<String> = cover.iter().cloned().collect();
    let group_matches = group == computed;
    if !group_matches {
        eprintln!(
            "WARNING: the `branch-coverage` group in cases.toml ({}) differs from the computed cover ({})",
            group.iter().cloned().collect::<Vec<_>>().join(", "),
            computed.iter().cloned().collect::<Vec<_>>().join(", ")
        );
    }

    // 4. gate every case of the cover against its fixture
    let mut gated: Vec<drive::CaseGate> = vec![];
    for name in &cover {
        let cand = candidates.iter().find(|c| &c.case.name == name).unwrap();
        println!("== {name} ==");
        gated.push(drive::gate(&root, &cand.case));
    }

    // 4b. the whole-run gate on the non-finite probes too — not a claim (their reference cannot
    // reproduce itself), a measurement of how far yFoil follows each one before the runs part
    let mut nonfinite_gated: Vec<drive::CaseGate> = vec![];
    for cand in candidates.iter().filter(|c| c.non_finite) {
        println!("== {} (non-finite, whole run) ==", cand.case.name);
        nonfinite_gated.push(drive::gate(&root, &cand.case));
    }

    // 5. outputs
    let stamp = chrono_utc_now();
    let run_dir = runs_root().join(&stamp);
    std::fs::create_dir_all(&run_dir).unwrap();
    let metadata = serde_json::json!({
        "study": "branch-coverage",
        "generated": stamp,
        "git_describe": git(&root, &["describe", "--tags", "--always"]).unwrap_or_else(|| "untagged".into()),
        "git_commit": git(&root, &["rev-parse", "HEAD"]).unwrap_or_default(),
        "git_dirty": git(&root, &["status", "--porcelain"]).map(|s| !s.is_empty()).unwrap_or(true),
        "host": format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        "universe_cases": universe.cases,
        "group_matches_cover": group_matches,
    });
    std::fs::write(
        run_dir.join("metadata.json"),
        serde_json::to_string_pretty(&metadata).unwrap(),
    )
    .unwrap();
    let summary = report::summary(&universe, &candidates, &cover, &gated, &nonfinite_gated, group_matches);
    std::fs::write(
        run_dir.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .unwrap();
    report::write_index(&run_dir, &summary, &universe, &candidates);
    report::write_index_tex(&run_dir, &summary);
    let drawn = render(&run_dir);
    if docs {
        let docs_dir = root.join("docs/validation/branch-coverage");
        std::fs::create_dir_all(&docs_dir).unwrap();
        report::write_docs(&docs_dir, &run_dir, &summary, &universe, &candidates, drawn);
        println!("docs written to {}", docs_dir.display());
    }
    println!("run folder: {}", run_dir.display());
    if !group_matches {
        eprintln!(
            "update `group = \"branch-coverage\"` in xtask/fixtures-config/cases.toml to the computed cover and re-run"
        );
        std::process::exit(1);
    }
}

/// UTC timestamp `YYYY-MM-DDTHH-MM-SSZ` without a chrono dependency.
fn chrono_utc_now() -> String {
    let out = std::process::Command::new("date")
        .args(["-u", "+%Y-%m-%dT%H-%M-%SZ"])
        .output()
        .expect("date");
    String::from_utf8_lossy(&out.stdout).trim().to_string()
}
