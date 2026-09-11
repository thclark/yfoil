//! `summary.json`: every value the plot and the index use, and the file-level comparison.
//! The plot is drawn from this file alone, so the chain XFOIL files → summary → figure is
//! reproducible from the run folder.

use crate::run::{BuildResult, Ending, Point};
use crate::{Build, ALPHA_STEP_DEG, ITER, MACH, NCRIT, PANELS, REYNOLDS};
use std::fs;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BuildSeries {
    pub build: Build,
    pub label: String,
    pub ending: Ending,
    /// Every converged point XFOIL recorded, in sweep order
    pub points: Vec<Point>,
    /// Requested angles with no converged record
    pub missing_alpha_deg: Vec<f64>,
    /// Last angle of the unbroken converged run from the start of the sweep: the plotted series
    /// ends here (one unconverged point ends it; nothing beyond is plotted)
    pub series_end_alpha_deg: Option<f64>,
    /// The plotted series: the points of the unbroken converged run, in sweep order
    pub series: Vec<Point>,
}

impl BuildSeries {
    /// The points that are plotted: the unbroken converged run
    pub fn series(&self) -> impl Iterator<Item = &Point> {
        let end = self.series_end_alpha_deg.unwrap_or(f64::NEG_INFINITY);
        self.points.iter().filter(move |p| p.alpha_deg <= end)
    }
}

/// Largest absolute difference between the builds over the points both recorded, per metric
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MetricDiff {
    pub metric: String,
    pub max_abs_diff: f64,
    pub at_alpha_deg: Option<f64>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Comparison {
    /// `polar.pol` byte-identical between the builds
    pub polar_identical: bool,
    /// `DUMP` files present in both builds, and how many of those are byte-identical
    pub dumps_compared: usize,
    pub dumps_identical: usize,
    pub dumps_differing: Vec<String>,
    /// stdout byte-identical (informative: the instrumented build may print more)
    pub stdout_identical: bool,
    /// Angles converged in both builds
    pub common_points: usize,
    pub metric_diffs: Vec<MetricDiff>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Conditions {
    pub panels: usize,
    pub re: f64,
    pub mach: f64,
    pub ncrit: f64,
    pub iter: usize,
    pub alpha_step_deg: f64,
    pub npoints: usize,
}

/// Extent of a quantity over the plotted series of every build
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Extent {
    pub min: f64,
    pub max: f64,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Summary {
    pub foil: String,
    pub conditions: Conditions,
    pub builds: Vec<BuildSeries>,
    pub comparison: Comparison,
    /// Extents of `alpha_deg` and every metric over the plotted series of both builds (positive
    /// values only for the trailing-edge BL variables, which are plotted on a log axis): what a
    /// plot's axis range is set from
    pub extents: std::collections::BTreeMap<String, Extent>,
}

/// A named accessor on a point
pub type Metric = (&'static str, fn(&Point) -> f64);

/// Metrics compared and plotted
pub const METRICS: &[Metric] = &[
    ("cl", |p| p.cl),
    ("cd", |p| p.cd),
    ("cdp", |p| p.cdp),
    ("cm", |p| p.cm),
    ("xtr_top", |p| p.xtr_top),
    ("xtr_bot", |p| p.xtr_bot),
    ("te_upper.ue", |p| p.te_upper.ue),
    ("te_upper.dstar", |p| p.te_upper.dstar),
    ("te_upper.theta", |p| p.te_upper.theta),
    ("te_upper.cf", |p| p.te_upper.cf),
    ("te_upper.hk", |p| p.te_upper.hk),
    ("te_lower.ue", |p| p.te_lower.ue),
    ("te_lower.dstar", |p| p.te_lower.dstar),
    ("te_lower.theta", |p| p.te_lower.theta),
    ("te_lower.cf", |p| p.te_lower.cf),
    ("te_lower.hk", |p| p.te_lower.hk),
];

fn series_end(points: &[Point], missing: &[f64]) -> Option<f64> {
    let first_missing = missing.iter().copied().fold(f64::INFINITY, f64::min);
    points
        .iter()
        .filter(|p| p.alpha_deg < first_missing)
        .map(|p| p.alpha_deg)
        .fold(None, |m, a| Some(m.map_or(a, |m: f64| m.max(a))))
}

pub fn build(foil: &str, npoints: usize, results: &[BuildResult]) -> Summary {
    let builds: Vec<BuildSeries> = results
        .iter()
        .map(|r| {
            let end = series_end(&r.points, &r.missing_alpha_deg);
            BuildSeries {
                build: r.build,
                label: r.build.label().to_string(),
                ending: r.ending,
                points: r.points.clone(),
                missing_alpha_deg: r.missing_alpha_deg.clone(),
                series_end_alpha_deg: end,
                series: r
                    .points
                    .iter()
                    .filter(|p| p.alpha_deg <= end.unwrap_or(f64::NEG_INFINITY))
                    .cloned()
                    .collect(),
            }
        })
        .collect();
    let mut extents = std::collections::BTreeMap::new();
    let extent = |name: &str, get: &dyn Fn(&Point) -> f64, positive: bool| -> Extent {
        let (lo, hi) = builds
            .iter()
            .flat_map(|b| b.series.iter())
            .map(get)
            .filter(|v| v.is_finite() && (!positive || *v > 0.0))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
        let _ = name;
        Extent { min: lo, max: hi }
    };
    // a quantity with no admissible value (a BL variable that is never positive) has no extent
    let mut insert = |name: &str, e: Extent| {
        if e.min.is_finite() && e.max.is_finite() {
            extents.insert(name.to_string(), e);
        }
    };
    insert("alpha_deg", extent("alpha_deg", &|p| p.alpha_deg, false));
    for (name, get) in METRICS {
        let positive = name.starts_with("te_");
        insert(name, extent(name, get, positive));
    }
    let (a, b) = (&results[0], &results[1]);
    fn file<'r>(r: &'r BuildResult, name: &str) -> Option<&'r Vec<u8>> {
        r.files.iter().find(|(n, _)| n == name).map(|(_, b)| b)
    }
    let polar_identical = file(a, crate::run::POLAR_NAME) == file(b, crate::run::POLAR_NAME);
    let mut dumps_compared = 0;
    let mut dumps_differing = Vec::new();
    for (name, bytes) in &a.files {
        if name == crate::run::POLAR_NAME {
            continue;
        }
        if let Some(other) = file(b, name) {
            dumps_compared += 1;
            if other != bytes {
                dumps_differing.push(name.clone());
            }
        }
    }
    let mut metric_diffs = Vec::new();
    let common: Vec<(&Point, &Point)> = a
        .points
        .iter()
        .filter_map(|p| b.points.iter().find(|q| q.alpha_deg == p.alpha_deg).map(|q| (p, q)))
        .collect();
    for (name, get) in METRICS {
        let mut best = (0.0_f64, None);
        for (p, q) in &common {
            let d = (get(p) - get(q)).abs();
            if d > best.0 {
                best = (d, Some(p.alpha_deg));
            }
        }
        metric_diffs.push(MetricDiff {
            metric: name.to_string(),
            max_abs_diff: best.0,
            at_alpha_deg: best.1,
        });
    }
    Summary {
        extents,
        foil: foil.to_string(),
        conditions: Conditions {
            panels: PANELS,
            re: REYNOLDS,
            mach: MACH,
            ncrit: NCRIT,
            iter: ITER,
            alpha_step_deg: ALPHA_STEP_DEG,
            npoints,
        },
        builds,
        comparison: Comparison {
            polar_identical,
            dumps_compared,
            dumps_identical: dumps_compared - dumps_differing.len(),
            dumps_differing,
            stdout_identical: a.stdout == b.stdout,
            common_points: common.len(),
            metric_diffs,
        },
    }
}

fn fmt_alpha(a: Option<f64>) -> String {
    a.map(|a| format!("{a:.1}°")).unwrap_or_else(|| "—".into())
}

/// `README.md` of the run folder
pub fn write_index(run_dir: &Path, s: &Summary) {
    let mut f = fs::File::create(run_dir.join("README.md")).unwrap();
    let c = &s.conditions;
    writeln!(f, "# XFOIL instrumentation cross-check — NACA {}\n", s.foil).unwrap();
    writeln!(
        f,
        "Generated by `cargo run --release -p xfoil-instrumentation-check` (provenance in `metadata.json`; every \
         plotted value and the file comparison in `summary.json`; XFOIL's raw files under `pristine/` and \
         `instrumented/`). Both double-precision XFOIL 6.99 builds ran the same OPER script: `NACA {}` with \
         `PPAR N {}`, Re {:e}, M {}, Ncrit {}, `ITER {}`, `PACC` on, one `ALFA` per point followed by a `DUMP`, \
         {} points from {}° in steps of {}°. Only XFOIL's standard output files are compared — the polar file \
         (`F7.3/F9.4/F10.5`) and the BL dumps (`F9.5/F10.6`) — never the instrumented build's own records.\n",
        s.foil,
        c.panels,
        c.re,
        c.mach,
        c.ncrit,
        c.iter,
        c.npoints,
        crate::ALPHA_START_DEG,
        c.alpha_step_deg
    )
    .unwrap();
    writeln!(
        f,
        "Figure: [naca{0}_check.svg](naca{0}_check.svg), and `naca{0}_check.pdf` for LaTeX, drawn by `plot.py` \
         (matplotlib, via `scripts/figures/render.sh`) from `summary.json` alone at the document's physical size in \
         its font (`scripts/figures/style.py`); include it at natural size: `\\includegraphics{{naca{0}_check.pdf}}`, \
         no `width=`, no `\\resizebox`.\n",
        s.foil,
    )
    .unwrap();
    writeln!(f, "## Files\n\n| File | Result |\n|---|---|").unwrap();
    let cmp = &s.comparison;
    writeln!(
        f,
        "| `polar.pol` | {} |",
        if cmp.polar_identical {
            "byte-identical"
        } else {
            "**differs**"
        }
    )
    .unwrap();
    writeln!(
        f,
        "| `alpha_*.bl` ({} compared) | {} byte-identical{} |",
        cmp.dumps_compared,
        cmp.dumps_identical,
        if cmp.dumps_differing.is_empty() {
            String::new()
        } else {
            format!("; **differ**: {}", cmp.dumps_differing.join(", "))
        }
    )
    .unwrap();
    writeln!(
        f,
        "| `stdout.txt` | {} (informative only) |\n",
        if cmp.stdout_identical {
            "byte-identical"
        } else {
            "differs"
        }
    )
    .unwrap();
    writeln!(
        f,
        "## Runs\n\n| Build | Converged points | First missing α | Series end α | Run |\n|---|---|---|---|---|"
    )
    .unwrap();
    for b in &s.builds {
        writeln!(
            f,
            "| {} | {} | {} | {} | {} |",
            b.label,
            b.points.len(),
            fmt_alpha(b.missing_alpha_deg.first().copied()),
            fmt_alpha(b.series_end_alpha_deg),
            b.ending.text()
        )
        .unwrap();
    }
    writeln!(
        f,
        "\n## Largest difference per metric over the {} angles converged in both builds\n\n| Metric | max abs diff | at α |\n|---|---|---|",
        cmp.common_points
    )
    .unwrap();
    for d in &cmp.metric_diffs {
        writeln!(
            f,
            "| `{}` | {:e} | {} |",
            d.metric,
            d.max_abs_diff,
            fmt_alpha(d.at_alpha_deg)
        )
        .unwrap();
    }
    writeln!(f, "\n## Plotted values\n").unwrap();
    writeln!(
        f,
        "The unbroken converged run of each build, as drawn (values as XFOIL wrote them).\n"
    )
    .unwrap();
    writeln!(
        f,
        "| α | Build | CL | CD | CM | Ue/V∞ TE up | δ* TE up | θ TE up | Hk TE up |\n|---|---|---|---|---|---|---|---|---|"
    )
    .unwrap();
    for b in &s.builds {
        for p in b.series() {
            writeln!(
                f,
                "| {:.1} | {} | {:.4} | {:.5} | {:.4} | {:.5} | {:.6} | {:.6} | {:.4} |",
                p.alpha_deg,
                b.build.slug(),
                p.cl,
                p.cd,
                p.cm,
                p.te_upper.ue,
                p.te_upper.dstar,
                p.te_upper.theta,
                p.te_upper.hk
            )
            .unwrap();
        }
    }
}

fn tex_alpha(a: Option<f64>) -> String {
    a.map(|a| format!("${a:.1}^\\circ$")).unwrap_or_else(|| "---".into())
}

/// `README.tex` of the run folder: the same index with LaTeX tables and the figure included
/// through the `svg` package (needs Inkscape at build time; drop the figure otherwise)
pub fn write_index_tex(run_dir: &Path, s: &Summary) {
    let mut f = fs::File::create(run_dir.join("README.tex")).unwrap();
    let w = |f: &mut fs::File, t: &str| writeln!(f, "{t}").unwrap();
    let c = &s.conditions;
    let cmp = &s.comparison;
    w(&mut f, "\\documentclass[11pt,a4paper]{article}");
    w(&mut f, "\\usepackage[margin=20mm]{geometry}");
    w(&mut f, "\\usepackage{amsmath,booktabs,longtable,siunitx,graphicx,url}");
    w(
        &mut f,
        &format!("\\title{{XFOIL instrumentation cross-check --- NACA {}}}", s.foil),
    );
    w(
        &mut f,
        "\\author{yFoil project}\n\\date{\\today}\n\\begin{document}\n\\maketitle",
    );
    w(&mut f, &format!(
        "\\noindent Generated by \\texttt{{cargo run --release -p xfoil-instrumentation-check}} (provenance in \
         \\texttt{{metadata.json}}; every plotted value and the file comparison in \\texttt{{summary.json}}; XFOIL's raw \
         files under \\texttt{{pristine/}} and \\texttt{{instrumented/}}). Both double-precision XFOIL~6.99 builds ran \
         the same OPER script: \\texttt{{NACA {}}} with \\texttt{{PPAR N {}}}, $Re = \\num{{{:e}}}$, $M = {}$, \
         $N_{{\\mathrm{{crit}}}} = {}$, \\texttt{{ITER {}}}, \\texttt{{PACC}} on, one \\texttt{{ALFA}} per point followed \
         by a \\texttt{{DUMP}}, {} points from ${}^\\circ$ in steps of ${}^\\circ$. Only XFOIL's standard output files \
         are compared --- the polar file (\\texttt{{F7.3/F9.4/F10.5}}) and the BL dumps (\\texttt{{F9.5/F10.6}}) --- \
         never the instrumented build's own records.",
        s.foil, c.panels, c.re, c.mach, c.ncrit, c.iter, c.npoints, crate::ALPHA_START_DEG, c.alpha_step_deg
    ));
    w(&mut f, "\n\\section*{Files}\n\\begin{table}[h]\\centering\\begin{tabular}{@{}ll@{}}\\toprule\nFile & Result \\\\ \\midrule");
    w(
        &mut f,
        &format!(
            "\\texttt{{polar.pol}} & {} \\\\",
            if cmp.polar_identical {
                "byte-identical"
            } else {
                "\\textbf{differs}"
            }
        ),
    );
    w(
        &mut f,
        &format!(
            "\\texttt{{alpha\\_*.bl}} ({} compared) & {} byte-identical{} \\\\",
            cmp.dumps_compared,
            cmp.dumps_identical,
            if cmp.dumps_differing.is_empty() {
                String::new()
            } else {
                format!(
                    "; \\textbf{{differ}}: {}",
                    cmp.dumps_differing.join(", ").replace('_', "\\_")
                )
            }
        ),
    );
    w(
        &mut f,
        &format!(
            "\\texttt{{stdout.txt}} & {} (informative only) \\\\",
            if cmp.stdout_identical {
                "byte-identical"
            } else {
                "differs"
            }
        ),
    );
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    w(
        &mut f,
        "\n\\section*{Runs}\n\\begin{table}[h]\\centering\\begin{tabular}{@{}lrrrl@{}}\\toprule",
    );
    w(
        &mut f,
        "Build & Converged points & First missing $\\alpha$ & Series end $\\alpha$ & Run \\\\ \\midrule",
    );
    for b in &s.builds {
        w(
            &mut f,
            &format!(
                "{} & {} & {} & {} & {} \\\\",
                b.label,
                b.points.len(),
                tex_alpha(b.missing_alpha_deg.first().copied()),
                tex_alpha(b.series_end_alpha_deg),
                match b.ending {
                    Ending::Completed => "completed",
                    Ending::Stalled => "XFOIL hung on a non-finite state; killed",
                    Ending::TimedOut => "killed after the timeout",
                }
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    w(&mut f, &format!(
        "\n\\section*{{Largest difference per metric}}\n\\noindent Over the {} angles converged in both builds.\n\
         \\begin{{table}}[h]\\centering\\begin{{tabular}}{{@{{}}lrr@{{}}}}\\toprule\nMetric & max $|\\Delta|$ & at $\\alpha$ \\\\ \\midrule",
        cmp.common_points
    ));
    for d in &cmp.metric_diffs {
        w(
            &mut f,
            &format!(
                "\\texttt{{{}}} & \\num{{{:e}}} & {} \\\\",
                d.metric.replace('_', "\\_"),
                d.max_abs_diff,
                tex_alpha(d.at_alpha_deg)
            ),
        );
    }
    w(&mut f, "\\bottomrule\\end{tabular}\\end{table}");
    w(
        &mut f,
        "\n\\section*{Figure}\n\\noindent The figure is drawn by \\texttt{plot.py} (matplotlib) from \
         \\texttt{summary.json} at the document's physical size in its font (\\texttt{scripts/figures/style.py}), \
         so it is included at natural size: no \\texttt{width=} and no \\texttt{\\resizebox}, or the text no \
         longer matches the document.\n",
    );
    w(&mut f, &format!(
        "\\begin{{figure}}[h]\\centering\\includegraphics{{naca{}_check.pdf}}\\caption{{XFOIL 6.99 DP, pristine (red, \
         crosses) against instrumented (blue, circles): NACA {}, $Re = \\num{{{:e}}}$, $M = {}$, $N_{{\\mathrm{{crit}}}} = {}$, \
         \\texttt{{ITER {}}}, one \\texttt{{ALFA}} per point. Lower right: upper-surface trailing-edge boundary layer, one \
         colour per variable, pristine solid, instrumented dashed. Each series ends at the first angle the build did not \
         converge.}}\\end{{figure}}",
        s.foil, s.foil, c.re, c.mach, c.ncrit, c.iter
    ));
    w(&mut f, "\n\\section*{Plotted values}\n\\noindent The unbroken converged run of each build, as drawn (values as XFOIL wrote them).\n");
    w(&mut f, "{\\small\\begin{longtable}{@{}rlrrrrrrr@{}}\\toprule");
    w(&mut f, "$\\alpha$ & Build & $C_L$ & $C_D$ & $C_M$ & $U_e/V_\\infty$ & $\\delta^*$ & $\\theta$ & $H_k$ \\\\ \\midrule\\endhead");
    for b in &s.builds {
        for p in b.series() {
            w(
                &mut f,
                &format!(
                    "{:.1} & {} & {:.4} & {:.5} & {:.4} & {:.5} & {:.6} & {:.6} & {:.4} \\\\",
                    p.alpha_deg,
                    b.build.slug(),
                    p.cl,
                    p.cd,
                    p.cm,
                    p.te_upper.ue,
                    p.te_upper.dstar,
                    p.te_upper.theta,
                    p.te_upper.hk
                ),
            );
        }
    }
    w(&mut f, "\\bottomrule\\end{longtable}}");
    w(&mut f, "\\end{document}");
}
