//! Branch-case polars (`docs/validation/branch-case-polars/`).
//!
//! The sections and conditions the branch-coverage study drew its cases from
//! (`scripts/branch-coverage`, `docs/validation/branch-coverage/`), each swept as a full polar —
//! 0 → +30°, `INIT`, 0 → −30° by 1°, XFOIL's polar procedure (`ALFA 0 / ASEQ`, ITMAX + 5 per
//! point, the sequence halting after NSEQEX = 4 consecutive unconverged points) — and yFoil
//! driven the same way through its `Session`, compared point by point with the same gate the
//! studies' floor comparison (`scripts/study-support/records.rs`). The point of the study is what
//! happens *past* convergence: through stall, and into the region where the reference goes
//! non-finite. There the reference cannot reproduce itself (its 1-ULP twins flip its branch
//! trace), so a departure is reported as divergent at the point where it happens, and
//! the figure shows where each code ends up.
//!
//! The cases are `group = "branch-case-polars"` in `xtask/fixtures-config/cases.toml`; their
//! fixtures come from `cargo xtask fixtures --group branch-case-polars`, whose watchdog ends a
//! reference run that has gone non-finite and hung in XFOIL's plot-label loop
//! (docs/xfoil-known-issues.md §4) and keeps what it wrote.
//!
//! Every invocation creates `runs/<UTC datetime>/` next to this crate with `metadata.json`,
//! `summary.json` (every point of every polar, both codes, with the outcome of each comparison),
//! `README.md` / `README.tex`; the figure (`polars.svg`, `polars.pdf`) is drawn from
//! `summary.json` by `plot.py` through `scripts/figures/render.sh`. With `--docs` the Markdown
//! page and the SVG are also written to `docs/validation/branch-case-polars/`.
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo run --release -p branch-case-polars -- [--docs]
//! ```

mod polars;
mod report;
mod utilities;

use serde::Deserialize;
use std::path::{Path, PathBuf};
use yfoil::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
use yfoil::solver::analysis::FlowConditions;

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
    pub matyp: usize,
    #[serde(default)]
    pub xtr: Vec<f64>,
    #[serde(default)]
    pub damp: bool,
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
    /// The OPER script in one line, for the tables.
    pub fn script(&self) -> String {
        let mut parts: Vec<String> = vec![
            format!("Re {:.0e}", self.re),
            format!("M {}", self.mach),
            format!("Ncrit {}", self.ncrit),
            format!("ITER {}", self.max_iterations),
        ];
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
            if a.len() > 1 {
                format!("ALFA {} / ASEQ {} {} {}", a[0], a[1], a[a.len() - 1], a[1] - a[0])
            } else {
                a.iter().map(|x| format!("ALFA {x}")).collect::<Vec<_>>().join(" / ")
            }
        };
        parts.push(seq(&self.alphas));
        if !self.alphas_after_reinit.is_empty() {
            parts.push(format!("INIT / {}", seq(&self.alphas_after_reinit)));
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
    /// The flow conditions the case's OPER script set.
    pub fn conditions(&self) -> FlowConditions {
        FlowConditions {
            re: Some(self.re),
            mach: self.mach,
            ncrit: self.ncrit,
            max_iterations: self.max_iterations,
            x_trip: if self.xtr.is_empty() {
                [1.0, 1.0]
            } else {
                [self.xtr[0], self.xtr[1]]
            },
            mach_cl_dependence: match self.matyp {
                2 => MachClDependence::InverseSqrtCl,
                _ => MachClDependence::Fixed,
            },
            re_cl_dependence: match self.matyp {
                2 => ReClDependence::InverseSqrtCl,
                3 => ReClDependence::InverseCl,
                _ => ReClDependence::Fixed,
            },
            amplification_model: if self.damp {
                AmplificationModel::ModifiedEnvelope
            } else {
                AmplificationModel::Envelope
            },
            ..FlowConditions::default()
        }
    }
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

/// Draw the run folder's figure with `scripts/figures/render.sh branch-case-polars <run_dir>`.
pub fn render(run_dir: &Path) -> bool {
    let root = repo_root();
    let status = std::process::Command::new(root.join("scripts/figures/render.sh"))
        .arg("branch-case-polars")
        .arg(run_dir)
        .status();
    match status {
        Ok(st) if st.success() => true,
        other => {
            eprintln!(
                "figure not drawn ({}): run `scripts/figures/render.sh branch-case-polars {}`",
                other.map(|s| s.to_string()).unwrap_or_else(|e| e.to_string()),
                run_dir.display()
            );
            false
        }
    }
}

/// The trailing `# comment` of each case's `name = "…"` line in cases.toml, its description.
fn descriptions(root: &Path) -> std::collections::BTreeMap<String, String> {
    let text = std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).unwrap_or_default();
    let mut out = std::collections::BTreeMap::new();
    for l in text.lines() {
        let Some(rest) = l.trim_start().strip_prefix("name = \"") else {
            continue;
        };
        let Some((name, after)) = rest.split_once('"') else {
            continue;
        };
        if let Some((_, c)) = after.split_once('#') {
            out.insert(name.to_string(), c.trim().to_string());
        }
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let docs = args.iter().any(|a| a == "--docs");
    let root = repo_root();
    let cases: Cases =
        toml::from_str(&std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("parse cases.toml");
    let desc = descriptions(&root);
    let mut runs: Vec<polars::PolarRun> = vec![];
    for case in cases
        .cases
        .iter()
        .filter(|c| c.group.as_deref() == Some("branch-case-polars"))
    {
        let dir = root.join("target/fixtures").join(&case.name);
        if !dir.join("viscal_points.dat").exists() {
            eprintln!(
                "  {}: no fixture under {} (run `cargo xtask fixtures --group branch-case-polars`) — skipped",
                case.name,
                dir.display()
            );
            continue;
        }
        println!("== {} ==", case.name);
        runs.push(polars::sweep(
            &root,
            case,
            desc.get(&case.name).map(String::as_str).unwrap_or(""),
        ));
    }
    assert!(!runs.is_empty(), "no polar fixtures found");

    let stamp = {
        let out = std::process::Command::new("date")
            .args(["-u", "+%Y-%m-%dT%H-%M-%SZ"])
            .output()
            .expect("date");
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };
    let run_dir = runs_root().join(&stamp);
    std::fs::create_dir_all(&run_dir).unwrap();
    let metadata = serde_json::json!({
        "study": "branch-case-polars",
        "generated": stamp,
        "git_describe": git(&root, &["describe", "--tags", "--always"]).unwrap_or_else(|| "untagged".into()),
        "git_commit": git(&root, &["rev-parse", "HEAD"]).unwrap_or_default(),
        "git_dirty": git(&root, &["status", "--porcelain"]).map(|s| !s.is_empty()).unwrap_or(true),
        "host": format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        "nseqex": polars::NSEQEX,
    });
    std::fs::write(
        run_dir.join("metadata.json"),
        serde_json::to_string_pretty(&metadata).unwrap(),
    )
    .unwrap();
    let summary = report::summary(&runs);
    std::fs::write(
        run_dir.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .unwrap();
    report::write_index(&run_dir, &summary);
    report::write_index_tex(&run_dir, &summary);
    let drawn = render(&run_dir);
    if docs {
        let docs_dir = root.join("docs/validation/branch-case-polars");
        std::fs::create_dir_all(&docs_dir).unwrap();
        report::write_docs(&docs_dir, &run_dir, &summary, drawn);
        println!("docs written to {}", docs_dir.display());
    }
    println!("run folder: {}", run_dir.display());
}
