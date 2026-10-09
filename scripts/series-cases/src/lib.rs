//! Series cases (`docs/validation/series-cases/`).
//!
//! Runs of the aerofoil series and flow conditions of interest, shown to an external reader:
//! full polars (XFOIL's polar procedure, 0 → +30°, `INIT`, 0 → −30° by 1°), inviscid alpha sweeps
//! (one SPECAL per alpha) and fixed-CL sweeps (one OPER `CL` per point), each with yFoil driven
//! the same way through its `Session` and compared point by point with the reference
//! (`scripts/study-support/records.rs`), and each run with its five seeded 1-ULP twins so the
//! study says where XFOIL itself was ill-conditioned. The point is what happens *past*
//! convergence: through stall, and into the region where the reference goes non-finite.
//!
//! The cases are `group = "series"` in `xtask/fixtures-config/cases.toml`; their
//! fixtures come from `cargo xtask fixtures --group series`, whose watchdog ends a
//! reference run that has gone non-finite and hung in XFOIL's plot-label loop
//! (docs/xfoil-known-issues.md §4) and keeps what it wrote.
//!
//! Every invocation creates `runs/<UTC datetime>/` next to this crate with `metadata.json`,
//! `summary.json` (every point of every polar, both codes, with the outcome of each comparison),
//! `README.md` / `README.tex`; the figure (`polars.svg`, `polars.pdf`) is drawn from
//! `summary.json` by `plot.py` through `scripts/figures/render.sh`. With `--docs` the Markdown
//! page and the SVG are also written to `docs/validation/series-cases/`.
//!
//! The library runs any [`Study`]: the series cases (`cargo run --release -p series-cases -- [--docs]`)
//! and the known-issues cases (`cargo run --release -p known-issues -- [--docs]`) share every
//! comparison, table and figure and differ only in their [`Study`] descriptor.

pub mod polars;
pub mod report;
pub mod utilities;

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
    /// OPER `CL` points (a fixed-CL series case)
    #[serde(default)]
    pub cls: Vec<f64>,
    /// no `VISC`: each ALFA is one SPECAL
    #[serde(default)]
    pub inviscid: bool,
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
        let mut parts: Vec<String> = if self.inviscid {
            vec!["inviscid".into(), format!("M {}", self.mach)]
        } else {
            vec![
                format!("Re {:.0e}", self.re),
                format!("M {}", self.mach),
                format!("Ncrit {}", self.ncrit),
                format!("ITER {}", self.max_iterations),
            ]
        };
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
        if self.inviscid {
            if let (Some(a), Some(b)) = (self.alphas.first(), self.alphas.last()) {
                parts.push(format!("ALFA {a} … {b} (one SPECAL each)"));
            }
        } else if self.alphas.is_empty() {
            if let (Some(a), Some(b)) = (self.cls.first(), self.cls.last()) {
                parts.push(format!("CL {a} … {b} ({} points)", self.cls.len()));
            }
        } else {
            parts.push(seq(&self.alphas));
        }
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
            re: (!self.inviscid).then_some(self.re),
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

pub fn git(root: &Path, args: &[&str]) -> Option<String> {
    let out = std::process::Command::new("git")
        .args(args)
        .current_dir(root)
        .output()
        .ok()?;
    out.status
        .success()
        .then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Draw the run folder's figures with `scripts/figures/render.sh <study> <run_dir>`.
pub fn render(study: &Study, run_dir: &Path) -> bool {
    let root = repo_root();
    let status = std::process::Command::new(root.join("scripts/figures/render.sh"))
        .arg(study.slug)
        .arg(run_dir)
        .status();
    match status {
        Ok(st) if st.success() => true,
        other => {
            eprintln!(
                "figure not drawn ({}): run `scripts/figures/render.sh {} {}`",
                other.map(|s| s.to_string()).unwrap_or_else(|e| e.to_string()),
                study.slug,
                run_dir.display()
            );
            false
        }
    }
}

/// The trailing `# comment` of each case's `name = "…"` line in cases.toml, its description.
pub fn descriptions(root: &Path) -> std::collections::BTreeMap<String, String> {
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

/// What distinguishes one case study from another: the cases it runs, where it writes, and how
/// its pages introduce themselves. The series cases and the known-issues cases share every
/// comparison, table and figure.
pub struct Study {
    /// the crate and the docs directory (`scripts/<slug>/`, `docs/validation/<slug>/`)
    pub slug: &'static str,
    /// the `group` of `xtask/fixtures-config/cases.toml` it runs
    pub group: &'static str,
    /// the page title
    pub title: &'static str,
    /// the first paragraph of the page: what the study is for (Markdown)
    pub purpose: &'static str,
    /// the same for `README.tex`
    pub purpose_tex: &'static str,
}

impl Study {
    pub fn runs_root(&self) -> PathBuf {
        repo_root().join("scripts").join(self.slug).join("runs")
    }
}

/// The series cases (`docs/validation/series-cases/`).
pub const SERIES: Study = Study {
    slug: "series-cases",
    group: "series",
    title: "Series cases",
    purpose: "Runs of the aerofoil series and flow regimes of interest under normal conditions, shown to an \
              external reader: that yFoil follows the reference through each of them. A member is here \
              because it corresponds, more loosely, to a member of the branch-coverage cases, or because it \
              shows a generator family or a physical regime that yFoil can be used for; a run whose purpose is \
              to show one of XFOIL's known issues belongs to the known-issues cases instead.",
    purpose_tex: "Runs of the aerofoil series and flow regimes of interest under normal conditions: that yFoil follows the reference through each of them.",
};

/// The known-issues cases (`docs/validation/known-issues/`).
pub const KNOWN_ISSUES: Study = Study {
    slug: "known-issues",
    group: "known-issues",
    title: "Known-issues cases",
    purpose: "Runs that show one of XFOIL's known issues (`docs/xfoil-known-issues.md`) arising, and whether \
              yFoil does the same: the reference hanging in its plot-label loop, the Kármán–Tsien pole, and the \
              solution's and convergence's sensitivity to the panelling (a case run at N = 60, 160 and 240 \
              panels, each with the iteration limit it was studied with). Each case's description names the \
              section of the known issues it shows. The study is to grow until it shows every documented issue \
              one for one.",
    purpose_tex: "Runs that show one of XFOIL's known issues arising, and whether yFoil does the same.",
};

/// Run a study over its cases' fixtures; `--docs` also writes its page under `docs/validation/`.
pub fn run(study: &Study) {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let docs = args.iter().any(|a| a == "--docs");
    let root = repo_root();
    let cases: Cases =
        toml::from_str(&std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("parse cases.toml");
    let desc = descriptions(&root);
    let mut runs: Vec<polars::PolarRun> = vec![];
    for case in cases.cases.iter().filter(|c| c.group.as_deref() == Some(study.group)) {
        let dir = root.join("target/fixtures").join(&case.name);
        if !dir.join("viscal_points.dat").exists() && !dir.join("specal_points.dat").exists() {
            eprintln!(
                "  {}: no fixture under {} (run `cargo xtask fixtures --group {}`) — skipped",
                case.name,
                dir.display(),
                study.group
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
    assert!(!runs.is_empty(), "no {} fixtures found", study.group);

    let stamp = {
        let out = std::process::Command::new("date")
            .args(["-u", "+%Y-%m-%dT%H-%M-%SZ"])
            .output()
            .expect("date");
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };
    let run_dir = study.runs_root().join(&stamp);
    std::fs::create_dir_all(&run_dir).unwrap();
    let metadata = serde_json::json!({
        "study": study.slug,
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
    report::write_index(study, &run_dir, &summary);
    report::write_index_tex(study, &run_dir, &summary);
    let drawn = render(study, &run_dir);
    if docs {
        let docs_dir = root.join("docs/validation").join(study.slug);
        std::fs::create_dir_all(&docs_dir).unwrap();
        report::write_docs(study, &docs_dir, &run_dir, &summary, drawn);
        println!("docs written to {}", docs_dir.display());
    }
    println!("run folder: {}", run_dir.display());
}
