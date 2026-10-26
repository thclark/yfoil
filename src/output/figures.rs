//! Shared conventions for yFoil's ironlab figures: the colours, the ordered-family colours that
//! replace a colourmap, the parameters and tags the ironlab viewer sorts, filters and groups
//! by, and saving. The rules are in `docs/development/conventions/figures.md`.
//!
//! Figure size and fonts are ironlab's defaults (set up for AMS journals). A figure sets
//! `size_mm` only to change its aspect ratio, keeping the default width.

use std::path::Path;

use ironlab::prelude::*;

use crate::solver::analysis::FlowConditions;

/// ironlab's default figure width (mm), which a tall figure keeps while setting its height.
pub const WIDTH_MM: f64 = 160.0;

/// The reference (XFOIL) in every comparison.
pub const RED: Color = Color::rgb(0xd6 as f32 / 255.0, 0x27 as f32 / 255.0, 0x28 as f32 / 255.0);
/// yFoil in every comparison.
pub const BLUE: Color = Color::rgb(0x1f as f32 / 255.0, 0x5a as f32 / 255.0, 0xd6 as f32 / 255.0);
/// Non-finite behaviour.
pub const ORANGE: Color = Color::rgb(0xff as f32 / 255.0, 0x7f as f32 / 255.0, 0x0e as f32 / 255.0);
/// Secondary marks (grey diamonds, reference lines).
pub const GREY: Color = Color::rgb(0.5, 0.5, 0.5);
/// D3's category colours after red and blue, for unordered families of four.
pub const CATEGORY: [Color; 4] = [
    Color::rgb(0x2c as f32 / 255.0, 0xa0 as f32 / 255.0, 0x2c as f32 / 255.0),
    Color::rgb(0x94 as f32 / 255.0, 0x67 as f32 / 255.0, 0xbd as f32 / 255.0),
    Color::rgb(0xff as f32 / 255.0, 0x7f as f32 / 255.0, 0x0e as f32 / 255.0),
    Color::rgb(0x17 as f32 / 255.0, 0xbe as f32 / 255.0, 0xcf as f32 / 255.0),
];

/// ColorBrewer's YlGnBu, the nine stops matplotlib interpolates.
const YLGNBU: [[u8; 3]; 9] = [
    [0xff, 0xff, 0xd9],
    [0xed, 0xf8, 0xb1],
    [0xc7, 0xe9, 0xb4],
    [0x7f, 0xcd, 0xbb],
    [0x41, 0xb6, 0xc4],
    [0x1d, 0x91, 0xc0],
    [0x22, 0x5e, 0xa8],
    [0x25, 0x34, 0x94],
    [0x08, 0x1d, 0x58],
];

/// The colour of member `t` ∈ [0, 1] of an ordered family (lightest = smallest), on YlGnBu
/// from 0.2 so the near-white end is never used. This replaces the colourmap the matplotlib
/// figures used: ironlab lines take explicit colours, and these are the colours matplotlib drew
/// (its 256-entry lookup table, sampled the way `Colormap.__call__` samples it).
pub fn ylgnbu(t: f64) -> Color {
    let x = 0.2 + 0.8 * t.clamp(0.0, 1.0);
    let index = ((x * 256.0) as usize).min(255);
    let s = index as f64 / 255.0 * (YLGNBU.len() - 1) as f64;
    let lo = (s.floor() as usize).min(YLGNBU.len() - 2);
    let f = s - lo as f64;
    let c = |k: usize| (YLGNBU[lo][k] as f64 * (1.0 - f) + YLGNBU[lo + 1][k] as f64 * f) / 255.0;
    Color::rgb(c(0) as f32, c(1) as f32, c(2) as f32)
}

/// The same colour made translucent.
pub fn with_alpha(c: Color, alpha: f32) -> Color {
    Color::rgba(c.r, c.g, c.b, alpha)
}

/// The flow conditions as figure parameters, named as in yFoil's `conditions` record.
pub fn flow_condition_parameters(fig: Figure, c: &FlowConditions) -> Figure {
    let mut fig = fig
        .parameter("viscous", c.re.is_some())
        .parameter("mach", c.mach)
        .parameter("ncrit", c.ncrit)
        .parameter("max_iterations", c.max_iterations as i64)
        .parameter("x_trip_upper", c.x_trip[0])
        .parameter("x_trip_lower", c.x_trip[1])
        .parameter("mach_cl_dependence", serde_name(c.mach_cl_dependence))
        .parameter("re_cl_dependence", serde_name(c.re_cl_dependence))
        .parameter("amplification_model", serde_name(c.amplification_model));
    if let Some(re) = c.re {
        fig = fig.parameter("re", re);
    }
    fig
}

/// A conditions enum under its serde name, as the `conditions` record writes it.
fn serde_name(v: impl serde::Serialize) -> String {
    serde_json::to_value(v)
        .ok()
        .and_then(|v| v.as_str().map(str::to_owned))
        .unwrap_or_default()
}

/// Add a parameter only when it has a finite value: a missing or non-finite output (an inviscid
/// point's TE δ*, a side that never separates) is left out rather than written as a fake number.
pub fn optional_parameter(fig: Figure, name: &str, value: Option<f64>) -> Figure {
    match value {
        Some(v) if v.is_finite() => fig.parameter(name, v),
        _ => fig,
    }
}

/// Figure-saving errors.
#[derive(Debug, thiserror::Error)]
pub enum FigureError {
    #[error("figure {path} does not validate: {errors}")]
    Invalid { path: String, errors: String },
    #[error("figure {path}: {source}")]
    Ironlab { path: String, source: ironlab::Error },
}

/// Refuse a figure with validation errors, or whose layout leaves something off the page or
/// out of its legend: entries that do not fit, two entries of one legend with the same name and
/// different samples, or text beyond the page edge. The layout warnings are scene warnings, not
/// validation errors, so a figure with one would otherwise be written silently.
fn check(fig: &Figure, path: &Path) -> Result<(), FigureError> {
    let report = fig.validate();
    let mut errors: Vec<String> = report.errors.iter().map(|e| format!("{e:?}")).collect();
    if errors.is_empty() {
        errors.extend(
            fig.layout_warnings()
                .into_iter()
                .filter(|w| {
                    matches!(
                        w.kind,
                        ironlab::SceneWarningKind::LegendOverflow { .. }
                            | ironlab::SceneWarningKind::LegendNameClash { .. }
                            | ironlab::SceneWarningKind::TextOutsidePage { .. }
                    )
                })
                .map(|w| w.message),
        );
    }
    if errors.is_empty() {
        return Ok(());
    }
    Err(FigureError::Invalid {
        path: path.display().to_string(),
        errors: errors.join("; "),
    })
}

/// Validate and save a figure; the encoding follows the extension (`.fig` or `.fig.json`).
/// A figure that does not validate or does not lay out ([`check`]) is never written.
pub fn save(fig: &Figure, path: &Path) -> Result<(), FigureError> {
    check(fig, path)?;
    fig.save(path).map_err(|source| FigureError::Ironlab {
        path: path.display().to_string(),
        source,
    })
}

/// Validate and save `<stem>.fig` and, for the run folders' LaTeX reports, `<stem>.pdf`.
pub fn save_with_pdf(fig: &Figure, stem: &Path) -> Result<(), FigureError> {
    let fig_path = stem.with_extension("fig");
    check(fig, &fig_path)?;
    let pdf = stem.with_extension("pdf");
    fig.export_pdf(&pdf).map_err(|source| FigureError::Ironlab {
        path: pdf.display().to_string(),
        source,
    })?;
    save(fig, &fig_path)
}

/// The documentation markup that embeds a figure file, `file_name` relative to the page (a study
/// section's pages embed `figures/<name>.fig`; `docs/development/conventions/figures.md`,
/// "Embedding"). The figure's name is the file's stem. The docs embed `.fig` only: `.fig.json`
/// is too slow to load in the site.
pub fn embed_html(file_name: &str, alt: &str) -> String {
    let stem = file_name.rsplit('/').next().unwrap();
    let name = stem.trim_end_matches(".fig");
    let alt = alt.replace('&', "&amp;").replace('"', "&quot;").replace('<', "&lt;");
    // `format` is explicit: the site's cache-busting appends `?v=<hash>` to `src`
    format!("<ironlab-figure src=\"{file_name}\" format=\"fig\" name=\"{name}\" alt=\"{alt}\">{alt}</ironlab-figure>")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ylgnbu_matches_matplotlibs_lookup_table() {
        // matplotlib: colormaps["YlGnBu"](0.2) and (1.0)
        let lo = ylgnbu(0.0);
        let hi = ylgnbu(1.0);
        // matplotlib 3.11.1: colormaps["YlGnBu"](0.2 + 0.8 * 0.37)[:3] * 255
        let mid = ylgnbu(0.37);
        for (got, want) in [(mid.r, 67.918), (mid.g, 183.082), (mid.b, 195.576)] {
            assert!((got * 255.0 - want).abs() < 1e-3, "{got} vs {want}");
        }
        assert!((hi.r - 8.0 / 255.0).abs() < 1e-6 && (hi.b - 88.0 / 255.0).abs() < 1e-6);
        // index 51 of 256: 0.6 of the way from the second stop to the third
        assert!((lo.r - ((0xed as f32) * 0.4 + (0xc7 as f32) * 0.6) / 255.0).abs() < 1e-5);
    }

    #[test]
    fn embed_escapes_alt_and_names_from_the_file() {
        let s = embed_html("polars-a.fig", "CL \"and\" CD");
        assert!(s.contains("name=\"polars-a\""));
        assert!(s.contains("&quot;and&quot;"));
        let s = embed_html("series-cases/polars-a.fig", "CL");
        assert!(s.contains("src=\"series-cases/polars-a.fig\"") && s.contains("name=\"polars-a\""));
    }
}
