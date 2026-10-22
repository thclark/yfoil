//! The family sheets: one ironlab figure per generator family, drawn from a run folder's
//! `summary.json` and the geometry files it names, so a run's figures can be redrawn without
//! rerunning the study (`--render <run_dir>`).
//!
//! Layout: two columns of panels, one per parameter, each panel followed in the row below by its
//! own legend panel (row-major; the default's entry is last, because ironlab lists entries in
//! drawing order and the default is drawn on top). Every panel has the same linked axes
//! (`X_RANGE` × `Y_RANGE`); the members are drawn 0.75 pt in `figures::ylgnbu(i / (n − 1))`
//! (already ordered by parameter value), the family default black 1.5 pt on top. ironlab has no
//! equal-aspect axes, so the figure height and the row weights are chosen here so that the plot
//! areas come out close to equal scale at the default width.
//!
//! A legend has three columns of entries where they fit the panel's width and fewer where they do
//! not (a long default label, such as the Kármán–Trefftz one): the sheet is laid out, compiled by
//! ironlab, and any legend ironlab reports as overflowing is given one column fewer and a taller
//! row, until none overflows.

use std::path::{Path, PathBuf};

use ironlab::prelude::*;
use ironlab::SceneWarningKind;
use serde_json::Value;
use yfoil::output::figures;

use crate::{X_RANGE, Y_RANGE};

/// Columns of panels
const COLUMNS: usize = 2;
/// Columns of entries in a legend, where they fit
const LEGEND_COLUMNS: usize = 3;

// Heights (mm) from which the figure height and the row weights are set. The plot area of a panel
// is the column width less the space ironlab reserves for the tick labels and the y label
// (PLOT_AREA_WIDTH_MM); its height follows from the equal-scale aspect, and the plot row adds the
// title above and the x tick labels and label below. Measured on an export at ironlab's defaults
// (160 mm wide, 9 pt).
const PLOT_AREA_WIDTH_MM: f64 = 63.5;
const PLOT_ROW_EXTRA_MM: f64 = 16.6;
const LEGEND_ENTRY_ROW_MM: f64 = 4.0;
const LEGEND_PAD_MM: f64 = 2.5;
const TITLE_MM: f64 = 8.0;

/// A geometry file's coordinates
fn coordinates(path: &Path) -> (Vec<f64>, Vec<f64>) {
    let g = yfoil::geometry::read_geometry_from_file(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    (g.x, g.y)
}

fn str_of<'a>(v: &'a Value, key: &str) -> &'a str {
    v[key]
        .as_str()
        .unwrap_or_else(|| panic!("summary.json: missing string {key}"))
}

fn f64_of(v: &Value, key: &str) -> f64 {
    v[key]
        .as_f64()
        .unwrap_or_else(|| panic!("summary.json: missing number {key}"))
}

fn members(panel: &Value) -> &Vec<Value> {
    panel["members"].as_array().expect("summary.json: members")
}

/// The family sheet of one `summary.json` series entry, with `legend_columns[k]` columns of
/// entries in panel k's legend; returns the figure and the legend panels' ids
pub fn family_sheet(run_dir: &Path, series: &Value, n_nodes: i64, legend_columns: &[usize]) -> (Figure, Vec<NodeId>) {
    let slug = str_of(series, "slug");
    let panels = series["panels"].as_array().expect("summary.json: panels");
    let default = &series["default"];
    let rows = panels.len().div_ceil(COLUMNS);

    // row weights: a plot row, then the legend row its panels' legends need
    let plot_area_h = PLOT_AREA_WIDTH_MM * (Y_RANGE.1 - Y_RANGE.0) / (X_RANGE.1 - X_RANGE.0);
    let plot_row = plot_area_h + PLOT_ROW_EXTRA_MM;
    let mut weights = Vec::with_capacity(2 * rows);
    for r in 0..rows {
        let entry_rows = (r * COLUMNS..((r + 1) * COLUMNS).min(panels.len()))
            .map(|k| (members(&panels[k]).len() + 1).div_ceil(legend_columns[k]))
            .max()
            .unwrap_or(1);
        weights.push(plot_row);
        weights.push(entry_rows as f64 * LEGEND_ENTRY_ROW_MM + 2.0 * LEGEND_PAD_MM);
    }
    let height = weights.iter().sum::<f64>() + TITLE_MM;

    let mut fig = Figure::new()
        .size_mm(figures::WIDTH_MM, height)
        .tiles((2 * rows) as u32, COLUMNS as u32)
        .row_weights(weights)
        .title(Text::plain(str_of(series, "title")))
        .label("validation")
        .parameter("study", "geometry-validation")
        .parameter("figure", "family-sheet")
        .parameter("series", slug)
        .parameter("n_nodes", n_nodes)
        .parameter("designation", str_of(default, "designation"))
        .parameter("max_thickness", f64_of(default, "max_thickness"))
        .parameter("te_gap", f64_of(default, "te_gap"));

    let (default_x, default_y) = coordinates(&run_dir.join(slug).join(str_of(default, "file")));
    let default_label = format!("{} (default)", str_of(default, "label"));
    let mut axes = Vec::with_capacity(panels.len());
    let mut legends = Vec::with_capacity(panels.len());
    for (k, panel) in panels.iter().enumerate() {
        let (r, c) = ((k / COLUMNS) as u32, (k % COLUMNS) as u32);
        let key = {
            let mut legend = fig.panel_legend(2 * r + 1, c);
            legend
                .columns(LegendEntryColumns::Count(legend_columns[k] as u32))
                .order(LegendEntryOrder::RowMajor)
                .align(LegendLocation::North)
                .boxed(false);
            legend.id()
        };
        let members = members(panel);
        let span = (members.len().max(2) - 1) as f64;
        let mut ax = fig.axes(2 * r, c);
        for (i, m) in members.iter().enumerate() {
            let (x, y) = coordinates(&run_dir.join(slug).join(str_of(panel, "slug")).join(str_of(m, "file")));
            ax.plot(&x, &y)
                .display_name(Text::plain(str_of(m, "label")))
                .color(figures::ylgnbu(i as f64 / span))
                .line_width(0.75);
        }
        // the default last, so that it is drawn on top
        ax.plot(&default_x, &default_y)
            .display_name(Text::plain(default_label.clone()))
            .color(Color::rgb(0.0, 0.0, 0.0))
            .line_width(1.5);
        ax.title(Text::plain(str_of(panel, "title")))
            .xlabel("$x/c$")
            .ylabel("$y/c$")
            .xlim(X_RANGE.0, X_RANGE.1)
            .ylim(Y_RANGE.0, Y_RANGE.1)
            .legend_in(key);
        axes.push(ax.id());
        legends.push(key);
    }
    // the same x/c and y/c on every panel
    fig.link(Dim::X, &axes).expect("the panels belong to this figure");
    fig.link(Dim::Y, &axes).expect("the panels belong to this figure");
    (fig, legends)
}

/// The legends of a figure that ironlab reports as overflowing, compiling it through a PDF export
/// to `probe`
fn overflowing(fig: &Figure, probe: &Path) -> Vec<NodeId> {
    let report = fig
        .export_pdf(probe)
        .unwrap_or_else(|e| panic!("{}: {e}", probe.display()));
    report
        .scene
        .iter()
        .filter(|w| matches!(w.kind, SceneWarningKind::LegendOverflow { .. }))
        .filter_map(|w| w.node)
        .collect()
}

/// Draw one series' sheet to `<stem>.fig` and `<stem>.pdf`, with as many legend columns (up to
/// three) as fit
fn draw(run_dir: &Path, series: &Value, n_nodes: i64, stem: &Path) {
    let n_panels = series["panels"].as_array().expect("summary.json: panels").len();
    let mut columns = vec![LEGEND_COLUMNS; n_panels];
    let probe = stem.with_extension("pdf");
    loop {
        let (fig, legends) = family_sheet(run_dir, series, n_nodes, &columns);
        let over = overflowing(&fig, &probe);
        if over.is_empty() {
            figures::save_with_pdf(&fig, stem).unwrap_or_else(|e| panic!("{e}"));
            return;
        }
        for (k, id) in legends.iter().enumerate() {
            if over.contains(id) {
                assert!(
                    columns[k] > 1,
                    "{}: panel {k}'s legend overflows in one column",
                    stem.display()
                );
                columns[k] -= 1;
            }
        }
    }
}

/// Draw every family sheet of a run folder (`<slug>.fig` and `<slug>.pdf`); returns the stems
pub fn render(run_dir: &Path) -> Vec<PathBuf> {
    let summary_path = run_dir.join("summary.json");
    let summary: Value = serde_json::from_str(
        &std::fs::read_to_string(&summary_path).unwrap_or_else(|e| panic!("{}: {e}", summary_path.display())),
    )
    .expect("summary.json");
    let n_nodes = summary["n_nodes"].as_i64().expect("summary.json: n_nodes");
    summary["series"]
        .as_array()
        .expect("summary.json: series")
        .iter()
        .map(|s| {
            let stem = run_dir.join(str_of(s, "slug"));
            draw(run_dir, s, n_nodes, &stem);
            stem
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use ironlab::ir::Parameter;
    use ironlab::prelude::*;

    /// Every family has its page in the documentation, named the kebab-case of its title and
    /// linked from the section's index, and the page embeds its sheet from the directory beside
    /// it named after it; the sheet loads and names its series
    #[test]
    fn documentation_holds_every_family_sheet() {
        let docs = crate::repo_root().join(crate::DOCS_DIR);
        let index = std::fs::read_to_string(docs.join("index.md")).unwrap();
        for s in crate::study::series() {
            let kebab: String = s
                .title
                .to_lowercase()
                .replace('á', "a")
                .chars()
                .map(|c| if c.is_alphanumeric() { c } else { '-' })
                .collect();
            assert_eq!(s.slug, kebab, "{}: the slug names the page", s.title);
            assert!(
                index.contains(&format!("]({}.md)", s.slug)),
                "{} not linked from the index",
                s.slug
            );
            let page = std::fs::read_to_string(docs.join(format!("{}.md", s.slug))).unwrap();
            let file = format!("{0}/{0}.fig", s.slug);
            let fig = Figure::load(docs.join(&file)).unwrap_or_else(|e| panic!("{file}: {e}"));
            let p = fig.parameters();
            assert_eq!(p.get("series"), Some(&Parameter::String(s.slug.into())), "{file}");
            assert_eq!(
                p.get("figure"),
                Some(&Parameter::String("family-sheet".into())),
                "{file}"
            );
            assert_eq!(
                p.get("study"),
                Some(&Parameter::String("geometry-validation".into())),
                "{file}"
            );
            assert_eq!(
                p.get("designation"),
                Some(&Parameter::String(s.default.designation.clone())),
                "{file}"
            );
            assert!(fig.labels().iter().any(|l| l == "validation"), "{file}");
            assert!(fig.validate().errors.is_empty(), "{file}");
            assert_eq!(fig.ir().axes.len(), s.panels.len(), "{file}: one panel per parameter");
            assert!(
                page.contains(&format!("<ironlab-figure src=\"{file}\"")),
                "{file} not embedded"
            );
        }
    }
}
