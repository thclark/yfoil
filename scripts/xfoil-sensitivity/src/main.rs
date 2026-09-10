//! XFOIL input-sensitivity study (`docs/validation/xfoil-sensitivity/`).
//!
//! Runs the instrumented double-precision XFOIL 6.99 reference on one aerofoil over an upward
//! alpha sweep, once for the base case and once per perturbation level of three input families —
//! node coordinates (1 ULP up to 1e-7 chord), panel count (160 down to 96) and alpha step
//! (0.5° plus 1 ULP up to 1e-5°) — and plots seven per-alpha quantities per family, base case in
//! front, perturbations stacked behind it coloured by perturbation size (YlGnBu).
//!
//! Every invocation creates a dated folder `runs/<UTC datetime>/` next to this crate holding
//! `metadata.json` (git version tag or `untagged`, commit sha, dirty flag, the reference build's
//! manifest and the study constants), one directory per XFOIL run under `<foil>/<family>/<level>/`
//! with the panels (`panels.dat`, written at 17 significant digits and asserted bitwise against
//! XFOIL's post-LOAD dump), the OPER script, XFOIL's stdout and the instrumented per-call records
//! (`viscal_points.dat`, `viscal_state_<k>.dat`), and the post-processed products: the SVG plots,
//! `summary.json` (the per-level table), `metrics.json` (every per-alpha quantity) and the
//! `README.md` / `README.tex` index. `runs/` is gitignored.
//!
//! Usage (from the repository root, reference built with `cargo xtask xfoil-build`):
//!
//! ```text
//! cargo run --release -p xfoil-sensitivity -- [--foil 0012] [--family geometry|panels|alpha-step]
//!                                            [--jobs N] [--resume RUN_DIR] [--plot-only RUN_DIR]
//! ```

mod perturb;
mod plot;
mod run;

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
}

fn parse_args() -> Args {
    let mut a = Args {
        foil: BASE_FOIL.to_string(),
        families: Family::ALL.to_vec(),
        jobs: 4,
        resume: None,
        plot_only: None,
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
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    a
}

fn main() {
    let args = parse_args();
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

    // per-family columns and the sheet
    for (family, results) in &families {
        let path = run_dir.join(format!("naca{}_{}.svg", args.foil, family.slug()));
        plot::draw_families(&path, &args.foil, std::slice::from_ref(&(*family, results.clone()))).unwrap();
        println!("wrote {}", path.display());
    }
    if families.len() == Family::ALL.len() {
        let path = run_dir.join(format!("naca{}_sheet.svg", args.foil));
        plot::draw_families(&path, &args.foil, &families).unwrap();
        println!("wrote {}", path.display());
    }
    run::write_summary_json(&run_dir, &args.foil, &families);
    run::write_index(&run_dir, &args.foil, &families);
    run::write_index_tex(&run_dir, &args.foil, &families);
    println!("run folder: {}", run_dir.display());
}
