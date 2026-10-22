//! The instrumentation-check figure, built from `summary.json` alone as an ironlab figure
//! (`docs/contributing/conventions/figures.md`): a legend panel across the top with one entry per build (its
//! converged count and where its plotted series ends), then a 2×2 grid against α — CL, CD and CM
//! (pristine red with crosses, instrumented blue with open circles) and the upper-surface
//! trailing-edge boundary layer on a logarithmic axis (one colour per variable, pristine solid,
//! instrumented dashed). α is linked across all four panels; no two panels plot the same
//! quantity, so no y axis is linked. Every number, axis extents included, is computed here.

use std::path::{Path, PathBuf};

use ironlab::prelude::*;
use yfoil::output::figures;
use yfoil::solver::analysis::FlowConditions;

use crate::run::Point;
use crate::summary::{Extent, Summary};
use crate::Build;

/// A plotted quantity: its `summary.json` metric (the key of its extent), accessor and label
type Quantity = (&'static str, fn(&Point) -> f64, &'static str);

/// The coefficient panels
const COEFFICIENTS: [Quantity; 3] = [
    ("cl", |p| p.cl, "$C_L$"),
    ("cd", |p| p.cd, "$C_D$"),
    ("cm", |p| p.cm, "$C_M$"),
];

/// The upper-surface trailing-edge variables of the fourth panel (Cf omitted: ≤ 0 once separated)
const BL_VARIABLES: [Quantity; 4] = [
    ("te_upper.ue", |p| p.te_upper.ue, r"$U_e/V_\infty$"),
    ("te_upper.dstar", |p| p.te_upper.dstar, r"$\delta^*$"),
    ("te_upper.theta", |p| p.te_upper.theta, r"$\theta$"),
    ("te_upper.hk", |p| p.te_upper.hk, "$H_k$"),
];

/// Height of the sheet (mm): a legend row above two rows of plots, taller than the default aspect
const HEIGHT_MM: f64 = 130.0;
const MARKER_PT: f64 = 4.0;

/// Colour, marker and dash of a build
fn style(build: Build) -> (Color, Marker, Dash) {
    match build {
        Build::Pristine => (figures::RED, Marker::Cross, Dash::Solid),
        Build::Instrumented => (figures::BLUE, Marker::Circle, Dash::Dashed),
    }
}

/// The figure's file stem in the run folder: `naca<foil>_check`
pub fn stem(run_dir: &Path, s: &Summary) -> PathBuf {
    run_dir.join(format!("naca{}_check", s.foil))
}

fn padded(e: &Extent) -> (f64, f64) {
    let pad = 0.06 * (e.max - e.min).max(1e-9);
    (e.min - pad, e.max + pad)
}

/// The flow conditions of the study's OPER script, as yFoil's `conditions` record names them:
/// XFOIL's defaults (free transition, TYPE 1, no DAMP) except where the script sets them
fn conditions(s: &Summary) -> FlowConditions {
    FlowConditions {
        re: Some(s.conditions.re),
        mach: s.conditions.mach,
        ncrit: s.conditions.ncrit,
        max_iterations: s.conditions.iter,
        ..FlowConditions::default()
    }
}

pub fn build(s: &Summary) -> Figure {
    let fig = Figure::new()
        .size_mm(figures::WIDTH_MM, HEIGHT_MM)
        .tiles(3, 2)
        .row_weights([0.5, 2.0, 2.0])
        .label("testing")
        .label("xfoil-comparison")
        .parameter("study", "xfoil-instrumentation-check")
        .parameter("figure", "instrumentation-check")
        .parameter("foil", format!("naca{}", s.foil))
        .parameter("n_nodes", s.conditions.panels as i64)
        .parameter("points", s.conditions.npoints as i64)
        .parameter("polar_identical", s.comparison.polar_identical)
        .parameter("dumps_compared", s.comparison.dumps_compared as i64)
        .parameter("dumps_identical", s.comparison.dumps_identical as i64);
    let mut fig = figures::flow_condition_parameters(fig, &conditions(s));
    let key = fig
        .panel_legend_span(0, 0, 1, 2)
        .boxed(false)
        .columns(LegendEntryColumns::Count(1))
        .id();

    let alpha_max = s.extents.get("alpha_deg").map_or(0.0, |e| e.max);
    let (x_lo, x_hi) = (-0.5, alpha_max + 0.5);
    let alphas = |series: &[Point]| series.iter().map(|p| p.alpha_deg).collect::<Vec<_>>();
    let mut ids = Vec::new();

    for (i, (metric, get, label)) in COEFFICIENTS.iter().enumerate() {
        let mut ax = fig.axes(1 + i as u32 / 2, i as u32 % 2);
        for b in &s.builds {
            let (colour, marker, _) = style(b.build);
            let end = b.series_end_alpha_deg.map_or("—".to_string(), |a| format!("{a:.1}°"));
            let values: Vec<f64> = b.series.iter().map(get).collect();
            ax.plot(alphas(&b.series), values)
                .display_name(format!(
                    "{}: {} converged points, series to {end}",
                    b.label,
                    b.points.len()
                ))
                .color(colour)
                .line_width(0.75)
                .marker(marker)
                .marker_size(MARKER_PT);
        }
        ax.xlabel(r"$\alpha$ (°)")
            .ylabel(*label)
            .xlim(x_lo, x_hi)
            .legend_in(key);
        if let Some(e) = s.extents.get(*metric) {
            let (lo, hi) = padded(e);
            ax.ylim(lo, hi);
        }
        ids.push(ax.id());
    }

    // the log panel: the extent over every variable's positive values, widened by a factor 1.3
    let (lo, hi) = BL_VARIABLES
        .iter()
        .filter_map(|(metric, _, _)| s.extents.get(*metric))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), e| {
            (lo.min(e.min), hi.max(e.max))
        });
    let mut ax = fig.axes(2, 1);
    for (k, (_, get, label)) in BL_VARIABLES.iter().enumerate() {
        for b in &s.builds {
            let (_, marker, dash) = style(b.build);
            // a non-positive value has no place on a log axis: it breaks the line
            let values: Vec<f64> = b
                .series
                .iter()
                .map(get)
                .map(|v| if v > 0.0 { v } else { f64::NAN })
                .collect();
            let mut line = ax.semilogy(alphas(&b.series), values);
            line.color(figures::CATEGORY[k])
                .line_width(1.5)
                .dash(dash)
                .marker(marker)
                .marker_size(MARKER_PT);
            // one entry per variable, drawn from the pristine (solid) line
            if b.build == Build::Pristine {
                line.display_name(*label);
            }
        }
    }
    ax.xlabel(r"$\alpha$ (°)")
        .ylabel("value, upper-surface TE")
        .xlim(x_lo, x_hi)
        .legend(LegendLocation::West)
        .legend_columns(LegendEntryColumns::Count(2));
    if lo.is_finite() && hi.is_finite() {
        ax.ylim(lo / 1.3, hi * 1.3);
    }
    ids.push(ax.id());
    fig.link(Dim::X, &ids).expect("the four axes belong to this figure");
    fig
}

/// Build the figure from `summary.json` and write `naca<foil>_check.fig` and `.pdf` into the
/// run folder (the PDF for its LaTeX index)
pub fn render(run_dir: &Path, s: &Summary) -> Result<PathBuf, figures::FigureError> {
    let stem = stem(run_dir, s);
    figures::save_with_pdf(&build(s), &stem)?;
    Ok(stem.with_extension("fig"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A summary of two builds whose series end at different angles, with a non-positive Ue
    fn summary() -> Summary {
        let point = |a: f64| {
            serde_json::json!({ "alpha_deg": a, "cl": 0.1 * a, "cd": 0.005 + 1e-4 * a * a, "cdp": 0.001,
                "cm": -0.01 * a, "xtr_top": 0.5, "xtr_bot": 0.6,
                "te_upper": { "ue": 0.9 - 0.1 * a, "dstar": 0.003 * (1.0 + a), "theta": 0.002 * (1.0 + a), "cf": 0.001, "hk": 1.6 + a },
                "te_lower": { "ue": -0.9, "dstar": 0.003, "theta": 0.002, "cf": 0.001, "hk": 1.6 } })
        };
        let build = |name: &str, label: &str, n: usize| {
            let pts: Vec<_> = (0..n).map(|i| point(i as f64)).collect();
            serde_json::json!({ "build": name, "label": label, "ending": "Completed", "points": pts,
                "missing_alpha_deg": [], "series_end_alpha_deg": (n - 1) as f64, "series": pts })
        };
        let extent = |lo: f64, hi: f64| serde_json::json!({ "min": lo, "max": hi });
        serde_json::from_value(serde_json::json!({
            "foil": "0012",
            "conditions": { "panels": 160, "re": 1e6, "mach": 0.0, "ncrit": 9.0, "iter": 100, "alpha_step_deg": 1.0, "npoints": 12 },
            "builds": [build("pristine", "XFOIL 6.99 DP, pristine", 12), build("instrumented", "XFOIL 6.99 DP, instrumented", 11)],
            "comparison": { "polar_identical": false, "dumps_compared": 11, "dumps_identical": 11, "dumps_differing": [],
                "stdout_identical": false, "common_points": 11, "metric_diffs": [] },
            "extents": { "alpha_deg": extent(0.0, 11.0), "cl": extent(0.0, 1.1), "cd": extent(0.005, 0.0171),
                "cm": extent(-0.11, 0.0), "te_upper.ue": extent(0.1, 0.9), "te_upper.dstar": extent(0.003, 0.036),
                "te_upper.theta": extent(0.002, 0.024), "te_upper.hk": extent(1.6, 12.6) }
        }))
        .unwrap()
    }

    #[test]
    fn figure_validates_and_exports_without_legend_warnings() {
        let fig = build(&summary());
        let report = fig.validate();
        assert!(report.errors.is_empty(), "{:?}", report.errors);
        // temporary files go in the repository's .tmp/, never /tmp
        let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../.tmp")
            .join(format!("xfoil-instrumentation-check-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let export = fig.export_pdf(dir.join("check.pdf")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(export.scene.is_empty(), "{:?}", export.scene);
        assert_eq!(fig.parameters().get("re"), Some(&1e6.into()));
        assert_eq!(fig.labels(), ["testing", "xfoil-comparison"]);
    }
}
