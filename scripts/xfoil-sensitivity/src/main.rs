//! XFOIL input-sensitivity study (`docs/validation/xfoil-sensitivity/`).
//!
//! Runs the instrumented double-precision XFOIL 6.99 reference on one aerofoil over an upward
//! alpha sweep, once for the base case and once per perturbation level of three input families —
//! node coordinates (1 ULP up to 1e-7 chord), panel count (160 down to 96) and alpha step
//! (0.5° plus 1 ULP up to 1e-5°) — and records seven per-alpha quantities per level. `plot.py`
//! (matplotlib, through `scripts/figures/render.sh`) draws them from `metrics.json` alone: base
//! case in front, perturbations stacked behind it coloured by perturbation size (YlGnBu), the node
//! and alpha-step families side by side as one sheet, every family as a figure of the same column
//! width. Every number in a figure — including the axis extents — is written by this program;
//! the script only presents them.
//!
//! Every invocation creates a dated folder `runs/<UTC datetime>/` next to this crate holding
//! `metadata.json` (git version tag or `untagged`, commit sha, dirty flag, the reference build's
//! manifest and the study constants), one directory per XFOIL run under `<foil>/<family>/<level>/`
//! with the panels (`panels.dat`, written at 17 significant digits and asserted bitwise against
//! XFOIL's post-LOAD dump), the OPER script, XFOIL's stdout and the instrumented per-call records
//! (`viscal_points.dat`, `viscal_state_<k>.dat`), and the post-processed products: the figures (SVG and PDF),
//! `summary.json` (the per-level table), `metrics.json` (every per-alpha quantity) and the
//! `README.md` / `README.tex` index. `runs/` is gitignored.
//!
//! The **twins mode** (`--twins`) is the second study of this crate: the reference against its
//! own seeded 1-ULP twins, case by case (`src/twins.rs`). It reads the runs `cargo xtask
//! fixtures` and `cargo xtask twins` leave under `target/fixtures/<case>/` (the reference and
//! `ulp<seed>/`), archives
//! their per-call records under `twins/<case>/` in the run folder, derives every converged-point
//! scalar and boundary-layer scalar per point, the differences twin by twin, the envelope (the
//! largest over the twins) of every scalar and every state array, and how each call finished
//! (converged, or stopped at its iteration limit) — all in `twins.json` — and draws one sheet per
//! case: the baseline column (reference black on top, twins in YlGnBu beneath) and the error
//! column (log |twin − reference| per twin, the envelope in red, an orange diamond where a twin
//! finished differently from the reference). `--docs` publishes the sheets to
//! `docs/validation/xfoil-sensitivity/`. The default case set is the `twins-baseline` group of
//! `cases.toml`: the series-cases sections and conditions at 240 panels and ITER 200, swept
//! 0 → ±30° by 1°, so the reference is given every chance to converge; each case's `n_nodes` and
//! `max_iterations` are recorded in `twins.json` and tabled.
//!
//! Usage (from the repository root, reference built with `cargo xtask xfoil-build`):
//!
//! ```text
//! cargo run --release -p xfoil-sensitivity -- [--foil 0012] [--family geometry|panels|alpha-step]
//!                                            [--jobs N] [--resume RUN_DIR] [--plot-only RUN_DIR]
//! cargo run --release -p xfoil-sensitivity -- --twins [--group series] [--case NAME]...
//!                                            [--docs] [--plot-only RUN_DIR]
//! ```

mod perturb;
mod run;
mod twins;

use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------------------------
// Study definition. These ranges are the ones to change; everything else follows from them.
// ---------------------------------------------------------------------------------------------

/// Base operating point
pub const BASE_FOIL: &str = "0012";
pub const BASE_PANELS: usize = 160;
pub const REYNOLDS: f64 = 1.0e6;
pub const MACH: f64 = 0.0;
pub const NCRIT: f64 = 9.0;
pub const ITER: usize = 100;
pub const ALPHA_STEP_DEG: f64 = 0.5;
/// Points in the upward sweep: 0°, step, 2·step … (NPOINTS − 1)·step
pub const NPOINTS: usize = 51;

/// Geometry perturbation amplitudes in chord units, after the 1-ULP level (`perturb::ULP`).
pub const GEOMETRY_EPS: &[f64] = &[perturb::ULP, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7];
/// Panel counts of the perturbed cases (160 − 1, 2, 4, 8, 16, 32, 64)
pub const PANEL_COUNTS: &[usize] = &[159, 158, 156, 152, 144, 128, 96];
/// Alpha-step perturbations in degrees, after the 1-ULP level
pub const ALPHA_STEP_DELTAS: &[f64] = &[
    perturb::ULP,
    1e-15,
    1e-14,
    1e-13,
    1e-12,
    1e-11,
    1e-10,
    1e-9,
    1e-8,
    1e-7,
    1e-6,
    1e-5,
];
/// Seed of the node-perturbation pattern (`perturb::pattern`)
pub const SEED: u64 = 0x5eed_f011_0012_2026;
/// Kill an XFOIL run after this long, or once its stdout has not grown for `STALL_SECS`: the
/// reference spins forever in its plot-label formatter once CL or CD is non-finite
/// (docs/xfoil-known-issues.md §4), and it prints a line per iteration until then
pub const RUN_TIMEOUT_SECS: u64 = 1800;
pub const STALL_SECS: u64 = 90;

/// One input family of the study
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Geometry,
    Panels,
    AlphaStep,
}

impl Family {
    pub const ALL: [Family; 3] = [Family::Geometry, Family::Panels, Family::AlphaStep];
    pub fn slug(self) -> &'static str {
        match self {
            Family::Geometry => "geometry",
            Family::Panels => "panels",
            Family::AlphaStep => "alpha-step",
        }
    }
    pub fn title(self) -> &'static str {
        match self {
            Family::Geometry => "Node coordinates perturbed",
            Family::Panels => "Panel count reduced",
            Family::AlphaStep => "Alpha step perturbed",
        }
    }
}

/// One run: the base case (`magnitude` 0) or a perturbation level of a family
#[derive(Debug, Clone)]
pub struct Level {
    pub family: Family,
    /// Directory name under the family
    pub slug: String,
    /// Legend text
    pub label: String,
    /// Size of the perturbation, used only to order and colour the stack (0 = base)
    pub magnitude: f64,
    /// The level's defining value as tabled: ε or δ (`perturb::ULP` for the 1-ULP level), or N
    pub raw: f64,
    pub panels: usize,
    pub geometry_eps: Option<f64>,
    pub alpha_step: f64,
    /// Node count actually generated (yFoil's generator returns an even count), set when the
    /// run is read back
    pub actual_panels: Option<usize>,
}

pub fn levels(family: Family) -> Vec<Level> {
    let base = Level {
        family,
        slug: "base".into(),
        label: format!("base: N = {BASE_PANELS}, Δα = {ALPHA_STEP_DEG}°"),
        magnitude: 0.0,
        raw: 0.0,
        panels: BASE_PANELS,
        geometry_eps: None,
        alpha_step: ALPHA_STEP_DEG,
        actual_panels: None,
    };
    let mut out = vec![base.clone()];
    match family {
        Family::Geometry => {
            for &eps in GEOMETRY_EPS {
                out.push(Level {
                    slug: perturb::eps_slug(eps),
                    label: format!("nodes ± {}", perturb::eps_label(eps)),
                    magnitude: if eps == perturb::ULP { 1e-16 } else { eps },
                    raw: eps,
                    geometry_eps: Some(eps),
                    ..base.clone()
                });
            }
        }
        Family::Panels => {
            for &n in PANEL_COUNTS {
                out.push(Level {
                    slug: format!("n{n}"),
                    label: format!("N = {n}"),
                    magnitude: (BASE_PANELS - n) as f64,
                    raw: n as f64,
                    panels: n,
                    ..base.clone()
                });
            }
        }
        Family::AlphaStep => {
            for &d in ALPHA_STEP_DELTAS {
                let step = if d == perturb::ULP {
                    perturb::next_up(ALPHA_STEP_DEG)
                } else {
                    ALPHA_STEP_DEG + d
                };
                out.push(Level {
                    slug: perturb::eps_slug(d),
                    label: format!("Δα = {ALPHA_STEP_DEG}° + {}", perturb::eps_label(d)),
                    magnitude: if d == perturb::ULP { step - ALPHA_STEP_DEG } else { d },
                    raw: d,
                    alpha_step: step,
                    ..base.clone()
                });
            }
        }
    }
    out
}

/// The families drawn side by side on the sheet (`plot.py`); the others get a figure each
pub const SHEET: [Family; 2] = [Family::Geometry, Family::AlphaStep];

/// Draw the run folder's figures with `scripts/figures/render.sh xfoil-sensitivity <run_dir>`
/// (matplotlib through uv); on failure the data are complete and the command is printed.
pub fn render(run_dir: &Path) {
    let status = std::process::Command::new(repo_root().join("scripts/figures/render.sh"))
        .arg("xfoil-sensitivity")
        .arg(run_dir)
        .status();
    if !matches!(status, Ok(st) if st.success()) {
        eprintln!(
            "figures not drawn: run `scripts/figures/render.sh xfoil-sensitivity {}`",
            run_dir.display()
        );
    }
}

/// Repository root: this crate lives at scripts/xfoil-sensitivity
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
pub fn write_metadata(run_dir: &Path, foil: &str, families: &[Family]) {
    let root = repo_root();
    let sha = git(&root, &["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".into());
    let version = git(&root, &["describe", "--exact-match", "--tags", "HEAD"]).unwrap_or_else(|| "untagged".into());
    let dirty = git(&root, &["status", "--porcelain"])
        .map(|s| !s.is_empty())
        .unwrap_or(true);
    let reference: Vec<String> = std::fs::read_to_string(root.join("target/xfoil-ref/manifest.txt"))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    let j = serde_json::json!({
        "study": "xfoil-sensitivity",
        "created_utc": run_dir.file_name().unwrap().to_string_lossy(),
        "version": version,
        "commit": sha,
        "dirty": dirty,
        "reference": reference,
        "foil": format!("naca{foil}"),
        "families": families.iter().map(|f| f.slug()).collect::<Vec<_>>(),
        "base": { "panels": BASE_PANELS, "re": REYNOLDS, "mach": MACH, "ncrit": NCRIT, "iter": ITER,
                  "alpha_step_deg": ALPHA_STEP_DEG, "npoints": NPOINTS, "seed": format!("{SEED:#x}") },
        "levels": { "geometry_eps": GEOMETRY_EPS, "panel_counts": PANEL_COUNTS, "alpha_step_deltas": ALPHA_STEP_DELTAS },
        "ulp_marker": perturb::ULP,
        "stall_secs": STALL_SECS,
        "timeout_secs": RUN_TIMEOUT_SECS,
    });
    std::fs::write(run_dir.join("metadata.json"), serde_json::to_string_pretty(&j).unwrap()).unwrap();
}

struct Args {
    foil: String,
    families: Vec<Family>,
    jobs: usize,
    /// An existing run folder to continue (levels already marked done are kept)
    resume: Option<PathBuf>,
    /// An existing run folder to post-process only
    plot_only: Option<PathBuf>,
    /// The twins mode: the reference against its seeded 1-ULP twins, case by case
    twins: bool,
    /// `--group NAME` / `--case NAME` selection of the twins mode (default: the
    /// series-cases group)
    group: String,
    cases: Vec<String>,
    /// publish the twins sheets to docs/validation/xfoil-sensitivity/
    docs: bool,
}

fn parse_args() -> Args {
    let mut a = Args {
        foil: BASE_FOIL.to_string(),
        families: Family::ALL.to_vec(),
        jobs: 4,
        resume: None,
        plot_only: None,
        twins: false,
        group: "twins-baseline".into(),
        cases: vec![],
        docs: false,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--foil" => {
                a.foil = argv[i + 1].clone();
                i += 1;
            }
            "--family" => {
                a.families = vec![match argv[i + 1].as_str() {
                    "geometry" => Family::Geometry,
                    "panels" => Family::Panels,
                    "alpha-step" => Family::AlphaStep,
                    other => panic!("unknown family {other}"),
                }];
                i += 1;
            }
            "--jobs" => {
                a.jobs = argv[i + 1].parse().expect("--jobs N");
                i += 1;
            }
            "--resume" => {
                a.resume = Some(PathBuf::from(&argv[i + 1]));
                i += 1;
            }
            "--plot-only" => {
                a.plot_only = Some(PathBuf::from(&argv[i + 1]));
                i += 1;
            }
            "--twins" => a.twins = true,
            "--docs" => a.docs = true,
            "--group" => {
                a.group = argv[i + 1].clone();
                i += 1;
            }
            "--case" => {
                a.cases.push(argv[i + 1].clone());
                i += 1;
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    a
}

/// A case of `xtask/fixtures-config/cases.toml`, as the twins mode needs it
#[derive(Debug, Clone, serde::Deserialize)]
pub struct Case {
    pub name: String,
    pub foil: String,
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
    pub track: bool,
    #[serde(default)]
    pub re: f64,
    #[serde(default)]
    pub mach: f64,
    #[serde(default)]
    pub ncrit: f64,
    #[serde(default)]
    pub max_iterations: usize,
    #[serde(default)]
    pub n_nodes: usize,
    #[serde(default)]
    pub inviscid: bool,
    #[serde(default)]
    pub cls: Vec<f64>,
}

#[derive(Debug, serde::Deserialize)]
struct Cases {
    case: Vec<Case>,
}

impl Case {
    /// The OPER conditions and script, as the polar study tables them
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
            if a.is_empty() {
                return String::new();
            }
            if self.polar && a.len() > 1 {
                format!("ALFA {} / ASEQ {} {} {}", a[0], a[1], a[a.len() - 1], a[1] - a[0])
            } else {
                a.iter().map(|v| format!("ALFA {v}")).collect::<Vec<_>>().join(" / ")
            }
        };
        parts.push(seq(&self.alphas));
        if !self.cls.is_empty() {
            parts.push(format!(
                "CL {} … {} ({} points)",
                self.cls[0],
                self.cls[self.cls.len() - 1],
                self.cls.len()
            ));
        }
        if !self.alphas_after_reinit.is_empty() {
            parts.push(format!("INIT / {}", seq(&self.alphas_after_reinit)));
        }
        parts.join("; ")
    }
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

/// The twins mode, start to finish
fn twins_main(args: &Args) {
    let root = repo_root();
    let cases: Cases =
        toml::from_str(&std::fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("cases.toml"))
            .expect("cases.toml");
    let selected: Vec<Case> = cases
        .case
        .iter()
        .filter(|c| {
            if !args.cases.is_empty() {
                args.cases.contains(&c.name)
            } else {
                c.group.as_deref() == Some(args.group.as_str())
            }
        })
        .cloned()
        .collect();
    assert!(
        !selected.is_empty(),
        "no case selected (group {}, cases {:?})",
        args.group,
        args.cases
    );
    let run_dir = match &args.plot_only {
        Some(d) => {
            assert!(
                d.join("twins.json").exists(),
                "{} is not a twins run folder",
                d.display()
            );
            d.clone()
        }
        None => {
            let d = runs_root().join(utc_stamp());
            std::fs::create_dir_all(&d).unwrap();
            d
        }
    };
    if args.plot_only.is_none() {
        let sha = git(&root, &["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".into());
        let version = git(&root, &["describe", "--exact-match", "--tags", "HEAD"]).unwrap_or_else(|| "untagged".into());
        let dirty = git(&root, &["status", "--porcelain"])
            .map(|s| !s.is_empty())
            .unwrap_or(true);
        let reference: Vec<String> = std::fs::read_to_string(root.join("target/xfoil-ref/manifest.txt"))
            .unwrap_or_default()
            .lines()
            .map(str::to_string)
            .collect();
        let mut cases_json = vec![];
        for case in &selected {
            let (work, twin_dirs) = twins::work_dirs(&root, &case.name);
            let records = if case.inviscid {
                "specal_points.dat"
            } else {
                "viscal_points.dat"
            };
            assert!(
                work.join(records).exists() && !twin_dirs.is_empty(),
                "{}: no reference run with twins under {} — run `cargo xtask fixtures --case {2}` and `cargo xtask twins --case {2}` first",
                case.name,
                work.display(),
                case.name
            );
            let first_after = case.alphas_after_reinit.first().copied();
            let load = |dir: &Path, label: &str, seed: Option<u64>| {
                if case.inviscid {
                    twins::load_inviscid_run(dir, label, seed)
                } else {
                    twins::load_run(dir, label, seed, first_after)
                }
            };
            let reference = load(&work, "reference", None).expect("reference records");
            let mut tw = vec![];
            for (seed, dir) in &twin_dirs {
                if let Some(r) = load(dir, &format!("twin {seed}"), Some(*seed)) {
                    tw.push(r);
                }
            }
            let case_dir = run_dir.join("twins").join(&case.name);
            twins::archive(&work, &case_dir.join("reference"));
            for (seed, dir) in &twin_dirs {
                twins::archive(dir, &case_dir.join(format!("twin{seed}")));
            }
            println!(
                "== {} == reference {} points ({} converged, {} at the iteration limit), {} twins",
                case.name,
                reference.points.len(),
                reference.points.iter().filter(|p| p.converged).count(),
                reference.points.iter().filter(|p| p.capped).count(),
                tw.len()
            );
            let by_call = case.alphas.is_empty() && !case.cls.is_empty();
            let mut cj = twins::case_json(&case.name, &case.section(), &case.script(), &reference, &tw, by_call);
            // the basic parameters, in the record in case the case definition is lost
            cj["n_nodes"] = serde_json::json!(case.n_nodes);
            cj["max_iterations"] = serde_json::json!(case.max_iterations);
            cj["re"] = serde_json::json!(case.re);
            cj["mach"] = serde_json::json!(case.mach);
            cj["ncrit"] = serde_json::json!(case.ncrit);
            cj["alphas"] = serde_json::json!(case.alphas);
            cj["alphas_after_reinit"] = serde_json::json!(case.alphas_after_reinit);
            cj["inviscid"] = serde_json::json!(case.inviscid);
            // the rows of the case's sheet (presentation reads them; an inviscid case has no BL)
            cj["plot_rows"] = if case.inviscid {
                serde_json::json!([["cl", "$C_L$"], ["cm", "$C_M$"], ["cdp", "$C_{D,p}$"]])
            } else {
                serde_json::json!([
                    ["cl", "$C_L$"],
                    ["cd", "$C_D$"],
                    ["xtr_upper", "$x/c$ transition, upper"],
                    ["dstar_te_upper", "$\\delta^*$ at TE, upper"]
                ])
            };
            cases_json.push(cj);
        }
        let j = serde_json::json!({
            "study": "xfoil-sensitivity/twins",
            "created_utc": run_dir.file_name().unwrap().to_string_lossy(),
            "version": version, "commit": sha, "dirty": dirty, "reference": reference,
            "group": args.group, "scalars": twins::SCALARS,
            "cases": cases_json,
        });
        std::fs::write(run_dir.join("twins.json"), serde_json::to_string_pretty(&j).unwrap()).unwrap();
        twins::write_index(&run_dir, &j);
        twins::write_index_tex(&run_dir, &j);
    }
    render(&run_dir);
    if args.docs {
        let j: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(run_dir.join("twins.json")).unwrap()).unwrap();
        twins::write_docs(&root.join("docs/validation/xfoil-sensitivity"), &run_dir, &j);
    }
    println!("run folder: {}", run_dir.display());
}

fn main() {
    let args = parse_args();
    if args.twins {
        twins_main(&args);
        return;
    }
    let root = repo_root();
    let xfoil = root.join("target/xfoil-ref/instrumented/bin/xfoil");
    assert!(
        xfoil.exists(),
        "instrumented reference missing at {}: run `cargo xtask xfoil-build`",
        xfoil.display()
    );
    let run_dir = match (&args.plot_only, &args.resume) {
        (Some(d), _) => {
            assert!(d.join("metadata.json").exists(), "{} is not a run folder", d.display());
            d.clone()
        }
        (None, Some(d)) => {
            std::fs::create_dir_all(d).unwrap();
            d.clone()
        }
        (None, None) => {
            let d = runs_root().join(utc_stamp());
            std::fs::create_dir_all(&d).unwrap();
            d
        }
    };
    if args.plot_only.is_none() {
        write_metadata(&run_dir, &args.foil, &args.families);
    }
    let foil_dir = run_dir.join(format!("naca{}", args.foil));

    let mut families = Vec::new();
    for family in &args.families {
        let levels = levels(*family);
        if args.plot_only.is_none() {
            run::run_family(&xfoil, &foil_dir, &args.foil, &levels, args.jobs);
        }
        let results: Vec<run::RunResult> = levels
            .iter()
            .map(|l| run::load(&foil_dir.join(family.slug()).join(&l.slug), l))
            .collect();
        families.push((*family, results));
    }

    run::write_summary_json(&run_dir, &args.foil, &families);
    run::write_index(&run_dir, &args.foil, &families);
    run::write_index_tex(&run_dir, &args.foil, &families);
    render(&run_dir);
    println!("run folder: {}", run_dir.display());
}
