//! Driving the instrumented reference and reading its per-call records back.

use crate::{perturb, Family, Level, ITER, MACH, NCRIT, NPOINTS, REYNOLDS, RUN_TIMEOUT_SECS, SEED, STALL_SECS};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};
use yfoil::geometry::{naca_4digit, write_dat_file, write_geometry_to_json, Geometry, Thickness};

/// One alpha point of one run, from `viscal_points.dat` and `viscal_state_<k>.dat`
#[derive(Debug, Clone)]
#[allow(dead_code)] // every recorded quantity is kept; later stages of the study plot more of them
pub struct Point {
    pub alpha_deg: f64,
    pub converged: bool,
    pub iterations: usize,
    /// RMSBL of the final iteration
    pub rmsbl: f64,
    pub cl: f64,
    pub cd: f64,
    pub cm: f64,
    /// Transition x/c on the upper surface (XOCTR(1))
    pub xtr_upper: f64,
    /// δ* at the upper-surface trailing-edge station, DSTR(IBLTE(1), 1)
    pub dstar_te_upper: f64,
    /// Shape factor δ*/θ at that station
    pub h_te_upper: f64,
    /// x/c of the upper-surface trailing-edge separation point: the first station of the run of
    /// τ ≤ 0 stations that reaches the trailing edge; 1.0 when the flow is attached at the TE
    pub x_sep_upper: f64,
}

/// One run: its level and the points XFOIL recorded (fewer than NPOINTS if it was killed)
#[derive(Debug, Clone)]
pub struct RunResult {
    pub level: Level,
    pub points: Vec<Point>,
    pub ending: Ending,
}

fn geometry(foil: &str, n: usize) -> Geometry {
    naca_4digit(foil, n, Thickness::Perpendicular).expect("NACA 4-digit designation")
}

/// Write the run directory's inputs: `panels.json`, `panels.dat`, `xfoil.inp`
fn prepare(dir: &Path, foil: &str, level: &Level, base_pattern: &[(f64, f64)]) -> Geometry {
    fs::create_dir_all(dir).unwrap();
    let mut g = geometry(foil, level.panels);
    if let Some(eps) = level.geometry_eps {
        assert_eq!(
            level.panels,
            crate::BASE_PANELS,
            "geometry levels perturb the base node set"
        );
        g = perturb::perturb_geometry(&g, base_pattern, eps);
    }
    write_geometry_to_json(&g, dir.join("panels.json")).unwrap();
    write_dat_file(&g, &format!("NACA {foil} {}", level.slug), dir.join("panels.dat")).unwrap();
    let mut s = String::from("PLOP\nG F\n\nLOAD panels.dat\nOPER\n");
    s += &format!("VISC {REYNOLDS}\nMACH {MACH}\nVPAR\nN {NCRIT}\n\nITER {ITER}\n");
    for i in 0..NPOINTS {
        s += &format!("ALFA {:.17e}\n", i as f64 * level.alpha_step);
    }
    s += "\nQUIT\n";
    fs::write(dir.join("xfoil.inp"), s).unwrap();
    g
}

/// Rule 4: XFOIL's post-LOAD dump of the nodes must be bitwise what we wrote
fn check_handoff(dir: &Path, g: &Geometry) -> Result<usize, String> {
    let dump = fs::read_to_string(dir.join("xfoil_panels.dat")).map_err(|_| "xfoil_panels.dat missing")?;
    let mut n = 0;
    for line in dump.lines() {
        let p: Vec<&str> = line.split_whitespace().collect();
        if p.len() != 7 || p[0].parse::<usize>().is_err() {
            continue;
        }
        let i: usize = p[0].parse().unwrap();
        let (x, y): (f64, f64) = (p[1].parse().map_err(|_| "bad x")?, p[2].parse().map_err(|_| "bad y")?);
        if i == 0 || i > g.x.len() || x.to_bits() != g.x[i - 1].to_bits() || y.to_bits() != g.y[i - 1].to_bits() {
            return Err(format!("node {i} differs"));
        }
        n += 1;
    }
    if n != g.x.len() {
        return Err(format!("{n} of {} nodes dumped", g.x.len()));
    }
    Ok(n)
}

/// How a run ended
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ending {
    Completed,
    /// stdout stopped growing: XFOIL hung on a non-finite state and was killed
    Stalled,
    /// still running at RUN_TIMEOUT_SECS and was killed
    TimedOut,
}

fn run_xfoil(xfoil: &Path, dir: &Path) -> Ending {
    for f in fs::read_dir(dir).unwrap().flatten() {
        let name = f.file_name().to_string_lossy().to_string();
        if name.starts_with("viscal_") || name.starts_with("xfoil_") || name == "stdout.txt" {
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

/// Keep the study's inputs and per-call records; the instrumented build also writes a
/// subroutine log, pointer/Ue/DIJ and Newton dumps of a few hundred MB per run that nothing here
/// reads
fn prune(dir: &Path) {
    for f in fs::read_dir(dir).unwrap().flatten() {
        let name = f.file_name().to_string_lossy().to_string();
        let keep = name.starts_with("viscal_")
            || name.starts_with("panels.")
            || name == "xfoil.inp"
            || name == "xfoil_panels.dat"
            || name == "stdout.txt"
            || name == "done";
        if !keep {
            let _ = fs::remove_file(f.path());
        }
    }
}

fn complete(dir: &Path) -> bool {
    dir.join("done").exists()
}

/// Run every level of a family that has not completed yet, `jobs` at a time
pub fn run_family(xfoil: &Path, runs_root: &Path, foil: &str, levels: &[Level], jobs: usize) {
    let family = levels[0].family;
    // one pattern for the base geometry's node count (XFOIL's N: yFoil's generator returns N nodes)
    let base_nodes = geometry(foil, crate::BASE_PANELS).x.len();
    let base_pattern = perturb::pattern(SEED, base_nodes);
    let todo: Vec<&Level> = levels
        .iter()
        .filter(|l| !complete(&runs_root.join(family.slug()).join(&l.slug)))
        .collect();
    println!("{}: {} of {} runs to do", family.slug(), todo.len(), levels.len());
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..jobs.max(1) {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                let Some(level) = todo.get(i) else { break };
                let dir = runs_root.join(family.slug()).join(&level.slug);
                let _ = fs::remove_file(dir.join("done"));
                let g = prepare(&dir, foil, level, &base_pattern);
                let t = Instant::now();
                let ending = run_xfoil(xfoil, &dir);
                prune(&dir);
                match check_handoff(&dir, &g) {
                    Ok(n) => println!(
                        "  {}/{}: {} nodes bitwise identical, {:.0} s{}",
                        family.slug(),
                        level.slug,
                        n,
                        t.elapsed().as_secs_f64(),
                        match ending {
                            Ending::Completed => "",
                            Ending::Stalled => " (STALLED on a non-finite state, killed)",
                            Ending::TimedOut => " (TIMED OUT, killed)",
                        }
                    ),
                    Err(e) => panic!("{}/{}: HANDOFF FAILED: {e}", family.slug(), level.slug),
                }
                fs::write(dir.join("done"), format!("{ending:?}\n")).unwrap();
            });
        }
    });
}

fn kv_blocks(text: &str, start_key: &str) -> Vec<HashMap<String, String>> {
    let mut out: Vec<HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == start_key {
            out.push(HashMap::new());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

/// Fortran ES24.16 output, including gfortran's `Infinity` / `NaN` spellings
fn real(s: &str) -> f64 {
    match s.trim() {
        "Infinity" | "+Infinity" => f64::INFINITY,
        "-Infinity" => f64::NEG_INFINITY,
        t if t.contains("NaN") => f64::NAN,
        t => t.parse().unwrap_or(f64::NAN),
    }
}

/// `viscal_state_<k>.dat`: header and the side-1 station rows
/// (`XSSI X UEDG THET DSTR CTAU MASS TAU`)
struct State {
    header: HashMap<String, String>,
    side1: Vec<[f64; 8]>,
}

fn read_state(path: &Path) -> Option<State> {
    let text = fs::read_to_string(path).ok()?;
    let mut header = HashMap::new();
    let mut side1 = Vec::new();
    for l in text.lines() {
        if let Some(rest) = l.strip_prefix("BL(") {
            let (idx, vals) = rest.split_once(")=")?;
            let (is, _ibl) = idx.split_once(',')?;
            if is.trim() == "1" {
                let v: Vec<f64> = vals.split_whitespace().map(real).collect();
                side1.push(std::array::from_fn(|m| v.get(m).copied().unwrap_or(f64::NAN)));
            }
        } else if l.starts_with("NODE(") {
            continue;
        } else if let Some((k, v)) = l.split_once('=') {
            header.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    Some(State { header, side1 })
}

/// Read a completed run back
pub fn load(dir: &Path, level: &Level) -> RunResult {
    let ending = match fs::read_to_string(dir.join("done")).unwrap_or_default().trim() {
        "Stalled" => Ending::Stalled,
        "TimedOut" => Ending::TimedOut,
        _ => Ending::Completed,
    };
    let text = fs::read_to_string(dir.join("viscal_points.dat")).unwrap_or_default();
    let mut points = Vec::new();
    for (i, rec) in kv_blocks(&text, "CALL").iter().enumerate() {
        let k = i + 1;
        let st = read_state(&dir.join(format!("viscal_state_{k}.dat")));
        let (dstar, h, xsep) = match &st {
            Some(s) => {
                // rows are 1-based stations with slot 0 = the stagnation point; side 1 ends at IBLTE(1)
                let te: usize = s
                    .header
                    .get("IBLTE1")
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(s.side1.len());
                let row = |ibl: usize| s.side1.get(ibl - 1).copied().unwrap_or([f64::NAN; 8]);
                let te_row = row(te);
                let dstar = te_row[4];
                let h = te_row[4] / te_row[3];
                // walk back from the TE while τ ≤ 0
                let mut sep = te + 1;
                while sep > 2 && row(sep - 1)[7] <= 0.0 {
                    sep -= 1;
                }
                let xsep = if sep > te { 1.0 } else { row(sep)[1] };
                (dstar, h, xsep)
            }
            None => (f64::NAN, f64::NAN, f64::NAN),
        };
        points.push(Point {
            alpha_deg: real(&rec["ALFA"]).to_degrees(),
            converged: rec.get("LVCONV").map(|v| v == "T").unwrap_or(false),
            iterations: rec.get("NITDONE").and_then(|v| v.parse().ok()).unwrap_or(0),
            rmsbl: real(rec.get("RMSBL").map(String::as_str).unwrap_or("NaN")),
            cl: real(&rec["CL"]),
            cd: real(&rec["CD"]),
            cm: real(&rec["CM"]),
            xtr_upper: real(&rec["XOCTR1"]),
            dstar_te_upper: dstar,
            h_te_upper: h,
            x_sep_upper: xsep,
        });
    }
    // the panel family is labelled by the node count actually generated (yFoil's generator
    // returns an even count, so a request for 159 yields 158)
    let mut level = level.clone();
    if level.family == Family::Panels && level.magnitude > 0.0 {
        if let Ok(g) = fs::read_to_string(dir.join("panels.json")) {
            if let Some(n) = serde_json::from_str::<serde_json::Value>(&g)
                .ok()
                .and_then(|v| v["x"].as_array().map(|a| a.len()))
            {
                level.label = if n == level.panels {
                    format!("N = {n}")
                } else {
                    format!("N = {n} (requested {})", level.panels)
                };
                level.actual_panels = Some(n);
            }
        }
    }
    RunResult { level, points, ending }
}

fn ending_text(e: Ending) -> &'static str {
    match e {
        Ending::Completed => "completed",
        Ending::Stalled => "XFOIL hung on a non-finite state (known issues §4); killed",
        Ending::TimedOut => "killed after the timeout",
    }
}

/// `summary.json` (the per-level table) and `metrics.json` (every per-alpha quantity)
pub fn write_summary_json(run_dir: &Path, foil: &str, families: &[(Family, Vec<RunResult>)]) {
    let level_json = |r: &RunResult| {
        serde_json::json!({
            "family": r.level.family.slug(),
            "level": r.level.slug,
            "label": r.level.label,
            "magnitude": r.level.magnitude,
            "panels": r.level.panels,
            "geometry_eps": r.level.geometry_eps,
            "alpha_step_deg": r.level.alpha_step,
            "ending": format!("{:?}", r.ending),
        })
    };
    let summary: Vec<serde_json::Value> = families
        .iter()
        .flat_map(|(_, rs)| rs.iter())
        .map(|r| {
            let mut j = level_json(r);
            let o = j.as_object_mut().unwrap();
            o.insert("points".into(), r.points.len().into());
            o.insert(
                "converged".into(),
                r.points.iter().filter(|p| p.converged).count().into(),
            );
            o.insert(
                "first_unconverged_alpha_deg".into(),
                r.points.iter().find(|p| !p.converged).map(|p| p.alpha_deg).into(),
            );
            o.insert(
                "last_converged_alpha_deg".into(),
                r.points.iter().rev().find(|p| p.converged).map(|p| p.alpha_deg).into(),
            );
            j
        })
        .collect();
    fs::write(
        run_dir.join("summary.json"),
        serde_json::to_string_pretty(&serde_json::json!({ "foil": format!("naca{foil}"), "levels": summary })).unwrap(),
    )
    .unwrap();
    let metrics: Vec<serde_json::Value> = families
        .iter()
        .flat_map(|(_, rs)| rs.iter())
        .map(|r| {
            let mut j = level_json(r);
            let pts: Vec<serde_json::Value> = r
                .points
                .iter()
                .map(|p| {
                    serde_json::json!({
                        "alpha_deg": p.alpha_deg, "converged": p.converged, "iterations": p.iterations,
                        "rmsbl": p.rmsbl, "cl": p.cl, "cd": p.cd, "cm": p.cm, "xtr_upper": p.xtr_upper,
                        "dstar_te_upper": p.dstar_te_upper, "h_te_upper": p.h_te_upper, "x_sep_upper": p.x_sep_upper,
                    })
                })
                .collect();
            j.as_object_mut().unwrap().insert("points".into(), pts.into());
            j
        })
        .collect();
    fs::write(
        run_dir.join("metrics.json"),
        serde_json::to_string_pretty(&serde_json::json!({ "foil": format!("naca{foil}"), "runs": metrics })).unwrap(),
    )
    .unwrap();
}

/// `README.md` of the run folder: what was run, and how each level ended
pub fn write_index(run_dir: &Path, foil: &str, families: &[(Family, Vec<RunResult>)]) {
    let mut f = fs::File::create(run_dir.join("README.md")).unwrap();
    writeln!(f, "# XFOIL input sensitivity — NACA {foil}\n").unwrap();
    writeln!(
        f,
        "Generated by `cargo run --release -p xfoil-sensitivity -- --foil {foil}` (provenance in `metadata.json`, \
         table data in `summary.json`, every per-alpha quantity in `metrics.json`, XFOIL's raw records under \
         `naca{foil}/<family>/<level>/`); the instrumented double-precision XFOIL 6.99 reference, \
         `ALFA` per point (ITER {ITER}), Re {REYNOLDS:e}, M {MACH}, Ncrit {NCRIT}, {NPOINTS} points from 0° in \
         steps of {}°. Node perturbations: `src/perturb.rs`.\n",
        crate::ALPHA_STEP_DEG
    )
    .unwrap();
    writeln!(
        f,
        "Figures are drawn at the document's physical size ({} pt sheet width, {} at {}–{} pt) and come as SVG and \
         PDF; include the PDF at natural size (`\\includegraphics{{...pdf}}`, no `width=`, no `\\resizebox`).\n",
        figure_style::TEXT_WIDTH_PT,
        figure_style::FONT,
        figure_style::TICK_PT,
        figure_style::AXIS_LABEL_PT
    )
    .unwrap();
    writeln!(f, "| Figure | Columns |\n|---|---|").unwrap();
    writeln!(
        f,
        "| [naca{foil}_sheet.svg](naca{foil}_sheet.svg) | {} |",
        crate::plot::SHEET
            .iter()
            .map(|fm| format!("[{}](naca{foil}_{}.svg)", fm.title(), fm.slug()))
            .collect::<Vec<_>>()
            .join(", ")
    )
    .unwrap();
    writeln!(
        f,
        "| [naca{foil}_panels.svg](naca{foil}_panels.svg) | {} (same column width as the sheet) |",
        Family::Panels.title()
    )
    .unwrap();
    for (family, results) in families {
        writeln!(f, "\n## {}\n", family.title()).unwrap();
        writeln!(
            f,
            "| Level | Points | Converged | First unconverged α | Last converged α | Run |\n|---|---|---|---|---|---|"
        )
        .unwrap();
        for r in results {
            let conv = r.points.iter().filter(|p| p.converged).count();
            let first_fail = r
                .points
                .iter()
                .find(|p| !p.converged)
                .map(|p| format!("{:.1}°", p.alpha_deg));
            let last_conv = r
                .points
                .iter()
                .rev()
                .find(|p| p.converged)
                .map(|p| format!("{:.1}°", p.alpha_deg));
            writeln!(
                f,
                "| {} | {} | {} | {} | {} | {} |",
                r.level.label,
                r.points.len(),
                conv,
                first_fail.unwrap_or_else(|| "—".into()),
                last_conv.unwrap_or_else(|| "—".into()),
                ending_text(r.ending)
            )
            .unwrap();
        }
    }
}

/// The level's legend text in LaTeX
fn label_tex(level: &Level, actual_nodes: Option<usize>) -> String {
    let pow = |v: f64| -> String {
        if v == perturb::ULP {
            "1 ULP".to_string()
        } else {
            format!("$10^{{{}}}$", v.log10().round() as i32)
        }
    };
    if level.magnitude == 0.0 {
        return format!(
            "base: $N = {}$, $\\Delta\\alpha = {}^\\circ$",
            crate::BASE_PANELS,
            crate::ALPHA_STEP_DEG
        );
    }
    match level.family {
        Family::Geometry => format!("nodes $\\pm$ {}", pow(level.raw)),
        Family::Panels => match actual_nodes {
            Some(n) if n != level.panels => format!("$N = {n}$ (requested {})", level.panels),
            _ => format!("$N = {}$", level.panels),
        },
        Family::AlphaStep => {
            if level.raw == perturb::ULP {
                format!("$\\Delta\\alpha = {}^\\circ + $ 1 ULP", crate::ALPHA_STEP_DEG)
            } else {
                format!(
                    "$\\Delta\\alpha = {}^\\circ + 10^{{{}}}{{}}^\\circ$",
                    crate::ALPHA_STEP_DEG,
                    level.raw.log10().round() as i32
                )
            }
        }
    }
}

/// `README.tex` of the run folder: the same index with LaTeX tables and the figures included
/// through the `svg` package (needs Inkscape at build time; drop the figures otherwise)
pub fn write_index_tex(run_dir: &Path, foil: &str, families: &[(Family, Vec<RunResult>)]) {
    let mut f = fs::File::create(run_dir.join("README.tex")).unwrap();
    let w = |f: &mut fs::File, s: &str| writeln!(f, "{s}").unwrap();
    w(&mut f, "\\documentclass[11pt,a4paper]{article}");
    w(&mut f, "\\usepackage[margin=20mm]{geometry}");
    w(&mut f, "\\usepackage{amsmath,booktabs,longtable,siunitx,graphicx,url}");
    w(&mut f, &format!("\\title{{XFOIL input sensitivity --- NACA {foil}}}"));
    w(
        &mut f,
        "\\author{yFoil project}\n\\date{\\today}\n\\begin{document}\n\\maketitle",
    );
    w(&mut f, &format!(
        "\\noindent Generated by \\texttt{{cargo run --release -p xfoil-sensitivity -- --foil {foil}}} (provenance in \
         \\texttt{{metadata.json}}, table data in \\texttt{{summary.json}}, per-alpha quantities in \\texttt{{metrics.json}}, \
         XFOIL's raw records under \\texttt{{naca{foil}/}}): the instrumented double-precision XFOIL~6.99 reference, one \
         \\texttt{{ALFA}} per point (\\texttt{{ITER {ITER}}}), $Re = \\num{{{REYNOLDS:e}}}$, $M = {MACH}$, \
         $N_{{\\mathrm{{crit}}}} = {NCRIT}$, {NPOINTS} points from $0^\\circ$ in steps of ${}^\\circ$. \
         The node-perturbation method is documented in \\texttt{{src/perturb.rs}}.",
        crate::ALPHA_STEP_DEG
    ));
    w(&mut f, &format!(
        "\n\\section*{{Figures}}\n\\noindent Figures are drawn at the document's physical size ({} pt sheet width, {} at \
         {}--{} pt), so they are included at natural size: no \\texttt{{width=}} and no \\texttt{{\\resizebox}}, or the text \
         no longer matches the document.\n\\begin{{table}}[h]\\centering\\begin{{tabular}}{{@{{}}ll@{{}}}}\\toprule\nFigure & Content \\\\ \\midrule",
        figure_style::TEXT_WIDTH_PT,
        figure_style::FONT,
        figure_style::TICK_PT,
        figure_style::AXIS_LABEL_PT
    ));
    w(
        &mut f,
        &format!("\\texttt{{naca{foil}\\_sheet.pdf}} & node coordinates and alpha step, one column each \\\\"),
    );
    for fm in Family::ALL {
        w(
            &mut f,
            &format!(
                "\\texttt{{naca{foil}\\_{}.pdf}} & {} (one column) \\\\",
                fm.slug().replace('-', "\\-"),
                fm.title()
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    for (family, results) in families {
        w(&mut f, &format!("\n\\section*{{{}}}", family.title()));
        w(
            &mut f,
            "{\\renewcommand\\arraystretch{1.1}\n\\begin{longtable}{@{}lrrrrl@{}}\\toprule",
        );
        w(&mut f, "Level & Points & Converged & First unconverged $\\alpha$ & Last converged $\\alpha$ & Run \\\\ \\midrule\\endhead");
        for r in results {
            let conv = r.points.iter().filter(|p| p.converged).count();
            let first_fail = r
                .points
                .iter()
                .find(|p| !p.converged)
                .map(|p| format!("${:.1}^\\circ$", p.alpha_deg));
            let last_conv = r
                .points
                .iter()
                .rev()
                .find(|p| p.converged)
                .map(|p| format!("${:.1}^\\circ$", p.alpha_deg));
            let actual = r
                .level
                .label
                .strip_prefix("N = ")
                .and_then(|t| t.split_whitespace().next())
                .and_then(|t| t.parse().ok());
            w(
                &mut f,
                &format!(
                    "{} & {} & {} & {} & {} & {} \\\\",
                    label_tex(&r.level, actual),
                    r.points.len(),
                    conv,
                    first_fail.unwrap_or_else(|| "---".into()),
                    last_conv.unwrap_or_else(|| "---".into()),
                    match r.ending {
                        Ending::Completed => "completed",
                        Ending::Stalled => "XFOIL hung on a non-finite state; killed",
                        Ending::TimedOut => "killed after the timeout",
                    }
                ),
            );
        }
        w(&mut f, "\\bottomrule\\end{longtable}}");
    }
    w(&mut f, "\n\\section*{Figures}\n\\noindent Axes span the converged points; off-scale values (unconverged states) are drawn on the axis edge. Open circles: RMSBL did not reach $10^{-4}$ in \\texttt{ITER} iterations. Rows, top to bottom: $C_L$, $C_D$, $\\delta^*$ and $H$ at the upper-surface trailing edge, upper-surface transition and trailing-edge separation $x/c$, and $\\log_{10}$ of the final RMSBL.\n");
    w(&mut f, &format!("\\begin{{figure}}[p]\\centering\\includegraphics{{naca{foil}_sheet.pdf}}\\caption{{XFOIL input sensitivity, NACA {foil}: node coordinates (left) and alpha step (right) perturbed; base case black, perturbations behind it coloured by size.}}\\end{{figure}}"));
    w(&mut f, &format!("\\begin{{figure}}[p]\\centering\\includegraphics{{naca{foil}_panels.pdf}}\\caption{{XFOIL input sensitivity, NACA {foil}: panel count reduced.}}\\end{{figure}}"));
    w(&mut f, "\\end{document}");
}
