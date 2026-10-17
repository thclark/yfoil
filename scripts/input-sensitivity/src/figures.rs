//! The study's figures, built from a run folder's `metadata.json` and `metrics.json` alone, so a
//! run folder can be redrawn without rerunning the reference (`--render RUN_DIR`).
//!
//! Seven stacked rows share α — CL, CD, δ* and H at the upper-surface trailing edge, the
//! upper-surface transition and trailing-edge separation x/c, and the final RMSBL on a log axis —
//! one column per input family. The base case is black, 1.5 pt, with filled markers, drawn last
//! so it lies on top; the perturbation levels are drawn behind it, largest first, coloured by
//! perturbation size (`figures::ylgnbu`, lightest = smallest). Unconverged points carry open
//! circles. The axes span the converged extents `metrics.json` records; a value beyond them (an
//! unconverged state can be anything, up to CD = 1e19) is drawn on the axis edge, and a value
//! that cannot be drawn (non-finite, or non-positive on the log axis) is left out. Each column
//! has a title and its own legend panel below it. The node and alpha-step families share a
//! two-column sheet; every family also gets a one-column figure.
//!
//! Each figure is saved as `<run>/naca<foil>_<name>.fig` with its PDF beside it (`README.tex`
//! includes the PDFs).

use std::collections::BTreeMap;
use std::path::Path;

use ironlab::prelude::*;
use serde::Deserialize;
use yfoil::output::figures;

use crate::{perturb, Family, SHEET};

#[derive(Deserialize)]
struct Metadata {
    foil: String,
    base: Base,
    ulp_marker: f64,
}

#[derive(Deserialize)]
struct Base {
    panels: usize,
    re: f64,
    mach: f64,
    ncrit: f64,
    iter: usize,
    alpha_step_deg: f64,
}

#[derive(Deserialize)]
struct Metrics {
    levels_by_family: BTreeMap<String, Vec<LevelRecord>>,
    extents: BTreeMap<String, Extent>,
}

#[derive(Deserialize, Clone, Copy)]
struct Extent {
    min: f64,
    max: f64,
}

#[derive(Deserialize)]
struct LevelRecord {
    family: String,
    magnitude: f64,
    raw: f64,
    panels: usize,
    actual_panels: usize,
    geometry_eps: Option<f64>,
    ending: String,
    points: Vec<PointRecord>,
}

/// A point as `metrics.json` writes it; a non-finite value is written as `null`
#[derive(Deserialize)]
struct PointRecord {
    alpha_deg: f64,
    converged: bool,
    rmsbl: Option<f64>,
    cl: Option<f64>,
    cd: Option<f64>,
    xtr_upper: Option<f64>,
    dstar_te_upper: Option<f64>,
    h_te_upper: Option<f64>,
    x_sep_upper: Option<f64>,
}

type Get = fn(&PointRecord) -> Option<f64>;

/// The rows, top to bottom: `metrics.json` key, value, axis label, log axis
const ROWS: [(&str, Get, &str, bool); 7] = [
    ("cl", |p| p.cl, "$C_L$", false),
    ("cd", |p| p.cd, "$C_D$", false),
    ("dstar_te_upper", |p| p.dstar_te_upper, "$\\delta^*$ at TE", false),
    ("h_te_upper", |p| p.h_te_upper, "$H$ at TE", false),
    ("xtr_upper", |p| p.xtr_upper, "$x/c$ transition", false),
    ("x_sep_upper", |p| p.x_sep_upper, "$x/c$ separation", false),
    ("rmsbl", |p| p.rmsbl, "RMSBL", true),
];

/// Height (mm) of a unit row weight: a data row has weight 1, so the figure's height (at the
/// default width) follows from its rows and its legend
const ROW_HEIGHT_MM: f64 = 29.0;
/// Legend columns: two under a full-width column, one under each column of the sheet (two
/// would not fit its width at the default font)
const LEGEND_COLUMNS_FULL_WIDTH: u32 = 2;
const LEGEND_COLUMNS_SHEET: u32 = 1;
const BASE_WIDTH_PT: f64 = 1.5;
const LEVEL_WIDTH_PT: f64 = 0.75;
const MARKER_PT: f64 = 2.4;
const OPEN_MARKER_PT: f64 = 4.4;
/// Height of one legend entry, and of the legend's insets, relative to a data row
const LEGEND_ENTRY_WEIGHT: f64 = 0.14;
const LEGEND_INSET_WEIGHT: f64 = 0.36;

fn exponent(v: f64) -> i32 {
    -(cr_libm::log10(v).round() as i32)
}

/// Legend text from the level's fields
fn level_label(level: &LevelRecord, meta: &Metadata) -> String {
    if level.magnitude == 0.0 {
        return format!(
            "base: $N$ = {}, $\\Delta\\alpha$ = {}°",
            meta.base.panels, meta.base.alpha_step_deg
        );
    }
    let ulp = meta.ulp_marker;
    let mut text = match level.family.as_str() {
        "geometry" => {
            let eps = level.geometry_eps.unwrap_or(level.raw);
            if eps == ulp {
                "nodes ± 1 ULP".to_string()
            } else {
                format!("nodes ± $10^{{-{}}}$", exponent(eps))
            }
        }
        "panels" => {
            let (n, req) = (level.actual_panels, level.panels);
            let mut t = format!("$N$ = {n}");
            if n != req {
                t += &format!(" (req. {req})");
            }
            t
        }
        _ => {
            if level.raw == ulp {
                "$\\Delta\\alpha$ + 1 ULP".to_string()
            } else {
                format!("$\\Delta\\alpha + 10^{{-{}}}{{}}^\\circ$", exponent(level.raw))
            }
        }
    };
    if level.ending != "Completed" {
        text += " (hung)";
    }
    text
}

/// The axis range of a row: the converged extent, padded by half of it on a linear axis, or
/// from a tenth of the smallest to 1e4 times the largest on the log axis
fn y_limits(e: Extent, log: bool) -> (f64, f64) {
    if log {
        (e.min / 10.0, e.max * 1e4)
    } else {
        let pad = 0.5 * (e.max - e.min).max(1e-12);
        (e.min - pad, e.max + pad)
    }
}

/// Levels in legend order: the base first, then by increasing magnitude
fn ordered(levels: &[LevelRecord]) -> Vec<&LevelRecord> {
    let mut out: Vec<&LevelRecord> = levels.iter().filter(|l| l.magnitude == 0.0).collect();
    let mut rest: Vec<&LevelRecord> = levels.iter().filter(|l| l.magnitude > 0.0).collect();
    rest.sort_by(|a, b| a.magnitude.total_cmp(&b.magnitude));
    out.extend(rest);
    out
}

/// One figure of the given families, one column each, saved as `naca<foil>_<name>`
fn draw(run_dir: &Path, meta: &Metadata, metrics: &Metrics, families: &[Family], name: &str) -> Result<(), String> {
    let foil = meta.foil.trim_start_matches("naca");
    let extent = |key: &str| {
        metrics
            .extents
            .get(key)
            .copied()
            .ok_or_else(|| format!("metrics.json has no extent for {key}"))
    };
    let alpha = extent("alpha_deg")?;
    let ncol = families.len() as u32;
    let nrow = ROWS.len() as u32;
    let legend_columns = if ncol > 1 {
        LEGEND_COLUMNS_SHEET
    } else {
        LEGEND_COLUMNS_FULL_WIDTH
    };
    let legend_rows = families
        .iter()
        .map(|f| metrics.levels_by_family.get(f.slug()).map_or(0, Vec::len))
        .max()
        .unwrap_or(0)
        .div_ceil(legend_columns as usize);
    let legend_weight = LEGEND_ENTRY_WEIGHT * legend_rows as f64 + LEGEND_INSET_WEIGHT;
    let mut weights = vec![1.0; ROWS.len()];
    weights.push(legend_weight);
    let height_mm = ROW_HEIGHT_MM * (ROWS.len() as f64 + legend_weight);

    let slugs: Vec<&str> = families.iter().map(|f| f.slug()).collect();
    let mut fig = Figure::new()
        .size_mm(figures::WIDTH_MM, height_mm)
        .tiles(nrow + 1, ncol)
        .row_weights(weights)
        .parameter("study", "input-sensitivity")
        .parameter("figure", "sensitivity-sheet")
        .parameter("foil", meta.foil.as_str())
        .parameter("family", slugs.join(","))
        .parameter("n_nodes", meta.base.panels as i64)
        .parameter("viscous", true)
        .parameter("re", meta.base.re)
        .parameter("mach", meta.base.mach)
        .parameter("ncrit", meta.base.ncrit)
        .parameter("max_iterations", meta.base.iter as i64)
        .parameter("alpha_step_deg", meta.base.alpha_step_deg)
        .label("testing");
    if meta.base.mach > 0.0 {
        fig = fig.label("compressible");
    }
    if meta.base.re < 1e5 {
        fig = fig.label("low-re");
    }

    let mut row_axes: Vec<Vec<NodeId>> = vec![Vec::new(); ROWS.len()];
    for (c, family) in families.iter().enumerate() {
        let c = c as u32;
        let levels = metrics
            .levels_by_family
            .get(family.slug())
            .ok_or_else(|| format!("metrics.json has no family {}", family.slug()))?;
        let levels = ordered(levels);
        let n = (levels.len().max(2) - 2).max(1) as f64;
        let key = {
            let mut panel = fig.panel_legend(nrow, c);
            panel
                .columns(LegendEntryColumns::Count(legend_columns))
                .order(LegendEntryOrder::RowMajor)
                .align(LegendLocation::North);
            panel.id()
        };
        for (r, (quantity, get, label, log)) in ROWS.iter().enumerate() {
            let (lo, hi) = y_limits(extent(quantity)?, *log);
            let mut ax = fig.axes(r as u32, c);
            // the stack, largest perturbation first so the base lands on top
            for (i, level) in levels.iter().enumerate().rev() {
                let is_base = level.magnitude == 0.0;
                let colour = if is_base {
                    Color::BLACK
                } else {
                    figures::ylgnbu((i - 1) as f64 / n)
                };
                // drawable points, clipped to the axis range
                let drawn: Vec<(f64, f64, bool)> = level
                    .points
                    .iter()
                    .filter_map(|p| {
                        let v = get(p)?;
                        (v.is_finite() && (!log || v > 0.0)).then(|| (p.alpha_deg, v.clamp(lo, hi), p.converged))
                    })
                    .collect();
                let a: Vec<f64> = drawn.iter().map(|d| d.0).collect();
                let v: Vec<f64> = drawn.iter().map(|d| d.1).collect();
                {
                    let mut line = ax.plot(&a, &v);
                    line.color(colour)
                        .line_width(if is_base { BASE_WIDTH_PT } else { LEVEL_WIDTH_PT })
                        .display_name(level_label(level, meta));
                    if is_base {
                        line.marker(Marker::Circle).marker_size(MARKER_PT).marker_face(colour);
                    }
                }
                let (ua, uv): (Vec<f64>, Vec<f64>) = drawn.iter().filter(|d| !d.2).map(|d| (d.0, d.1)).unzip();
                if !ua.is_empty() {
                    ax.plot(&ua, &uv)
                        .color(colour)
                        .dash(Dash::None)
                        .marker(Marker::Circle)
                        .marker_size(OPEN_MARKER_PT);
                }
            }
            ax.xlim(-0.5, alpha.max + 0.5)
                .ylim(lo, hi)
                .ylabel(*label)
                .legend_in(key);
            if *log {
                ax.yscale(Scale::Log);
            }
            if r == 0 {
                ax.title(family.title());
            }
            if r + 1 == ROWS.len() {
                ax.xlabel("$\\alpha$ (°)");
            }
            row_axes[r].push(ax.id());
        }
    }
    fig.link_all_x().map_err(|e| format!("{e}"))?;
    if ncol > 1 {
        for ids in &row_axes {
            fig.link(Dim::Y, ids).map_err(|e| format!("{e}"))?;
        }
    }
    figures::save_with_pdf(&fig, &run_dir.join(format!("naca{foil}_{name}"))).map_err(|e| e.to_string())
}

/// The figure names a run folder gets: one per family it holds, and the sheet if it holds both
/// of the sheet's families
pub fn names(families: &[String]) -> Vec<String> {
    let mut out = families.to_vec();
    if SHEET.iter().all(|f| families.iter().any(|s| s == f.slug())) {
        out.push("sheet".into());
    }
    out
}

/// Draw every figure of a run folder from its `metadata.json` and `metrics.json`; returns the
/// figure files written
pub fn render(run_dir: &Path) -> Result<Vec<String>, String> {
    let read =
        |f: &str| std::fs::read_to_string(run_dir.join(f)).map_err(|e| format!("{}: {e}", run_dir.join(f).display()));
    let meta: Metadata = serde_json::from_str(&read("metadata.json")?).map_err(|e| format!("metadata.json: {e}"))?;
    let metrics: Metrics = serde_json::from_str(&read("metrics.json")?).map_err(|e| format!("metrics.json: {e}"))?;
    debug_assert_eq!(meta.ulp_marker, perturb::ULP);
    let foil = meta.foil.trim_start_matches("naca").to_string();
    let present: Vec<String> = Family::ALL
        .iter()
        .map(|f| f.slug().to_string())
        .filter(|s| metrics.levels_by_family.contains_key(s))
        .collect();
    let mut written = Vec::new();
    for name in names(&present) {
        let families: Vec<Family> = if name == "sheet" {
            SHEET.to_vec()
        } else {
            Family::ALL.into_iter().filter(|f| f.slug() == name).collect()
        };
        draw(run_dir, &meta, &metrics, &families, &name)?;
        // the matplotlib SVG this figure replaces
        let stale = run_dir.join(format!("naca{foil}_{name}.svg"));
        if stale.exists() {
            std::fs::remove_file(&stale).map_err(|e| format!("{}: {e}", stale.display()))?;
        }
        written.push(format!("naca{foil}_{name}.fig"));
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sheet_needs_both_of_its_families() {
        let all: Vec<String> = ["geometry", "panels", "alpha-step"].map(String::from).to_vec();
        assert_eq!(names(&all).last().map(String::as_str), Some("sheet"));
        assert!(!names(&["geometry".into()]).contains(&"sheet".to_string()));
    }

    /// A run folder of two levels per family, with a non-finite and an off-scale value: every
    /// figure is written and loads back
    #[test]
    fn renders_a_run_folder() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../.tmp/input-sensitivity-render-test");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let meta = serde_json::json!({ "foil": "naca0012", "ulp_marker": 0.0, "base": {
            "panels": 160, "re": 1e6, "mach": 0.0, "ncrit": 9.0, "iter": 100, "alpha_step_deg": 0.5 } });
        let point = |a: f64, converged: bool, cd: Option<f64>| {
            serde_json::json!({ "alpha_deg": a, "converged": converged, "rmsbl": 1e-5, "cl": 0.1 * a,
                "cd": cd, "xtr_upper": 0.5, "dstar_te_upper": 0.01, "h_te_upper": 2.0, "x_sep_upper": 1.0 })
        };
        let level = |family: &str, magnitude: f64, raw: f64| {
            serde_json::json!({ "family": family, "magnitude": magnitude, "raw": raw, "panels": 160,
                "actual_panels": 160, "geometry_eps": null, "ending": "Completed",
                "points": [point(0.0, true, Some(0.01)), point(0.5, false, None), point(1.0, false, Some(1e19))] })
        };
        let extent = |min: f64, max: f64| serde_json::json!({ "min": min, "max": max });
        let metrics = serde_json::json!({
            "levels_by_family": {
                "geometry": [level("geometry", 0.0, 0.0), level("geometry", 1e-16, 0.0)],
                "panels": [level("panels", 0.0, 0.0), level("panels", 2.0, 158.0)],
                "alpha-step": [level("alpha-step", 0.0, 0.0), level("alpha-step", 1e-9, 1e-9)],
            },
            "extents": {
                "alpha_deg": extent(0.0, 1.0), "cl": extent(0.0, 0.1), "cd": extent(0.01, 0.01),
                "dstar_te_upper": extent(0.01, 0.01), "h_te_upper": extent(2.0, 2.0),
                "xtr_upper": extent(0.5, 0.5), "x_sep_upper": extent(1.0, 1.0), "rmsbl": extent(1e-5, 1e-5),
            },
        });
        std::fs::write(dir.join("metadata.json"), meta.to_string()).unwrap();
        std::fs::write(dir.join("metrics.json"), metrics.to_string()).unwrap();
        let written = render(&dir).unwrap();
        assert_eq!(written.len(), 4);
        for f in written {
            let fig = Figure::load(dir.join(&f)).unwrap();
            assert!(fig.validate().errors.is_empty(), "{f}");
            assert!(dir.join(f.replace(".fig", ".pdf")).exists());
        }
    }

    #[test]
    fn log_limits_span_four_decades_above() {
        let (lo, hi) = y_limits(Extent { min: 1e-6, max: 1e-4 }, true);
        assert_eq!((lo, hi), (1e-7, 1.0));
        let (lo, hi) = y_limits(Extent { min: 0.0, max: 1.0 }, false);
        assert_eq!((lo, hi), (-0.5, 1.5));
    }
}
