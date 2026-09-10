//! Driving one reference build and reading XFOIL's standard output files back.

use crate::{Build, ITER, MACH, NCRIT, PANELS, REYNOLDS, RUN_TIMEOUT_SECS, STALL_SECS};
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// Name of the `DUMP` file written after the `ALFA` at this angle
pub fn dump_name(alpha_deg: f64) -> String {
    format!("alpha_{alpha_deg:.3}.bl")
}

pub const POLAR_NAME: &str = "polar.pol";

/// Write the run directory's OPER script. Geometry comes from XFOIL's own `NACA` generator and
/// `PPAR` paneling, so both builds start from the same code at the same precision; the polar
/// accumulates every converged point in `polar.pol` (`PLRADD` appends as it goes), and each
/// `ALFA` is followed by a `DUMP` of the boundary layer at XFOIL's standard precision.
fn prepare(dir: &Path, foil: &str, alphas: &[f64]) {
    fs::create_dir_all(dir).unwrap();
    let mut s = String::from("PLOP\nG F\n\n");
    // PPAR: after a change, the first blank line regenerates the panels and returns to the
    // menu; the second leaves it
    s += &format!("NACA {foil}\nPPAR\nN {PANELS}\n\n\n");
    s += &format!("OPER\nVISC {REYNOLDS}\nMACH {MACH}\nVPAR\nN {NCRIT}\n\nITER {ITER}\n");
    s += &format!("PACC\n{POLAR_NAME}\n\n");
    for &a in alphas {
        s += &format!("ALFA {a}\nDUMP {}\n", dump_name(a));
    }
    s += "PACC\n\nQUIT\n";
    fs::write(dir.join("xfoil.inp"), s).unwrap();
}

/// How a run ended
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Ending {
    Completed,
    /// stdout stopped growing: XFOIL hung on a non-finite state and was killed
    Stalled,
    /// still running at RUN_TIMEOUT_SECS and was killed
    TimedOut,
}

impl Ending {
    pub fn text(self) -> &'static str {
        match self {
            Ending::Completed => "completed",
            Ending::Stalled => "XFOIL hung on a non-finite state (known issues §4); killed",
            Ending::TimedOut => "killed after the timeout",
        }
    }
}

fn run_xfoil(xfoil: &Path, dir: &Path) -> Ending {
    for f in fs::read_dir(dir).unwrap().flatten() {
        if f.file_name() != "xfoil.inp" {
            let _ = fs::remove_file(f.path());
        }
    }
    let inp = fs::File::open(dir.join("xfoil.inp")).unwrap();
    let out = fs::File::create(dir.join("stdout.txt")).unwrap();
    let mut child = Command::new(xfoil)
        .current_dir(dir)
        .stdin(inp)
        .stdout(out)
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn xfoil");
    let start = Instant::now();
    let (mut last_size, mut last_growth) = (0u64, Instant::now());
    loop {
        if let Some(_status) = child.try_wait().unwrap() {
            return Ending::Completed;
        }
        let size = fs::metadata(dir.join("stdout.txt")).map(|m| m.len()).unwrap_or(0);
        if size != last_size {
            last_size = size;
            last_growth = Instant::now();
        }
        let ending = if last_growth.elapsed() > Duration::from_secs(STALL_SECS) {
            Some(Ending::Stalled)
        } else if start.elapsed() > Duration::from_secs(RUN_TIMEOUT_SECS) {
            Some(Ending::TimedOut)
        } else {
            None
        };
        if let Some(e) = ending {
            let _ = child.kill();
            let _ = child.wait();
            return e;
        }
        std::thread::sleep(Duration::from_millis(500));
    }
}

/// Keep the script and XFOIL's standard outputs only. The instrumented build also writes its
/// subroutine log, per-call state and Newton dumps — hundreds of MB that this study, by
/// design, never reads.
fn prune(dir: &Path) {
    for f in fs::read_dir(dir).unwrap().flatten() {
        let name = f.file_name().to_string_lossy().to_string();
        let keep = name == "xfoil.inp"
            || name == "stdout.txt"
            || name == POLAR_NAME
            || (name.starts_with("alpha_") && name.ends_with(".bl"))
            || name == "done";
        if !keep {
            let _ = fs::remove_file(f.path());
        }
    }
}

/// Run one build unless its directory is already marked done
pub fn run_build(xfoil: &Path, dir: &Path, build: Build, foil: &str, alphas: &[f64]) {
    if dir.join("done").exists() {
        println!("{}: already done, kept", build.slug());
        return;
    }
    prepare(dir, foil, alphas);
    let t = Instant::now();
    let ending = run_xfoil(xfoil, dir);
    prune(dir);
    fs::write(dir.join("done"), format!("{ending:?}\n")).unwrap();
    println!(
        "{}: {:.0} s, {}",
        build.slug(),
        t.elapsed().as_secs_f64(),
        ending.text()
    );
}

/// Boundary-layer state at one trailing-edge station, as `DUMP` writes it
/// (`Ue/Vinf F9.5`, `Dstar F10.6`, `Theta F10.6`, `Cf F10.6`, `H F10.4`)
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TeState {
    pub ue: f64,
    pub dstar: f64,
    pub theta: f64,
    pub cf: f64,
    pub hk: f64,
}

/// One converged polar point with its trailing-edge BL state
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Point {
    pub alpha_deg: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdp: f64,
    pub cm: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    /// Node 1 of the dump: the upper-surface trailing-edge station (IBLTE(1))
    pub te_upper: TeState,
    /// Node N of the dump: the lower-surface trailing-edge station (IBLTE(2))
    pub te_lower: TeState,
}

/// One build's run, read back from its standard output files
#[derive(Debug, Clone)]
pub struct BuildResult {
    pub build: Build,
    pub ending: Ending,
    /// Converged points, in sweep order (the polar file records only LVCONV points)
    pub points: Vec<Point>,
    /// Requested angles that have no polar entry: unconverged, or after a kill
    pub missing_alpha_deg: Vec<f64>,
    /// Raw file bytes for the like-for-like comparison, keyed by file name
    pub files: Vec<(String, Vec<u8>)>,
    pub stdout: Vec<u8>,
}

/// The polar file's data rows, keyed by column name from its header line
fn read_polar(path: &Path) -> Vec<std::collections::HashMap<String, f64>> {
    let text = fs::read_to_string(path).unwrap_or_default();
    let lines: Vec<&str> = text.lines().collect();
    let Some(dash) = lines.iter().position(|l| l.trim().starts_with("----")) else {
        return Vec::new();
    };
    let names: Vec<String> = lines[dash - 1].split_whitespace().map(str::to_string).collect();
    lines[dash + 1..]
        .iter()
        .filter_map(|l| {
            let vals: Vec<f64> = l.split_whitespace().filter_map(|t| t.parse().ok()).collect();
            (vals.len() >= 7 && vals.len() <= names.len()).then(|| names.iter().cloned().zip(vals).collect())
        })
        .collect()
}

/// The two trailing-edge rows of a `DUMP` file: node 1 (upper TE) and node N (lower TE)
fn read_dump(path: &Path) -> Option<(TeState, TeState)> {
    let text = fs::read_to_string(path).ok()?;
    let rows: Vec<Vec<f64>> = text
        .lines()
        .filter(|l| !l.trim_start().starts_with('#'))
        .map(|l| l.split_whitespace().filter_map(|t| t.parse().ok()).collect())
        .filter(|v: &Vec<f64>| v.len() >= 8)
        .collect();
    let te = |r: &Vec<f64>| TeState {
        ue: r[3],
        dstar: r[4],
        theta: r[5],
        cf: r[6],
        hk: r[7],
    };
    let (first, last) = (rows.first()?, rows.get(PANELS - 1)?);
    // both trailing-edge nodes sit at x = 1 (PPAR N nodes; the wake rows follow)
    assert!(
        (first[1] - last[1]).abs() < 1e-3,
        "{}: rows 1 and {PANELS} are not both at the trailing edge",
        path.display()
    );
    Some((te(first), te(last)))
}

/// Read a completed run back
pub fn load(dir: &Path, build: Build, alphas: &[f64]) -> BuildResult {
    let ending = match fs::read_to_string(dir.join("done")).unwrap_or_default().trim() {
        "Stalled" => Ending::Stalled,
        "TimedOut" => Ending::TimedOut,
        _ => Ending::Completed,
    };
    let rows = read_polar(&dir.join(POLAR_NAME));
    let mut points = Vec::new();
    let mut missing = Vec::new();
    for &a in alphas {
        // F7.3 in the polar file: match to the requested angle within half a unit of it
        let Some(row) = rows.iter().find(|r| (r["alpha"] - a).abs() < 5e-4) else {
            missing.push(a);
            continue;
        };
        let Some((te_upper, te_lower)) = read_dump(&dir.join(dump_name(a))) else {
            missing.push(a);
            continue;
        };
        points.push(Point {
            alpha_deg: a,
            cl: row["CL"],
            cd: row["CD"],
            cdp: row["CDp"],
            cm: row["CM"],
            xtr_top: row["Top_Xtr"],
            xtr_bot: row["Bot_Xtr"],
            te_upper,
            te_lower,
        });
    }
    let mut files = vec![(
        POLAR_NAME.to_string(),
        fs::read(dir.join(POLAR_NAME)).unwrap_or_default(),
    )];
    for &a in alphas {
        if let Ok(b) = fs::read(dir.join(dump_name(a))) {
            files.push((dump_name(a), b));
        }
    }
    BuildResult {
        build,
        ending,
        points,
        missing_alpha_deg: missing,
        files,
        stdout: fs::read(dir.join("stdout.txt")).unwrap_or_default(),
    }
}
