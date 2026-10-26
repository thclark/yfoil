//! The family sheets: one ironlab figure per generator family, drawn from a run folder's
//! `summary.json` and the geometry files it names, so a run's figures can be redrawn without
//! rerunning the study (`--render <run_dir>`).
//!
//! Layout: two columns of panels, one per parameter, each panel followed in the row below by its
//! own legend panel. Every panel has the same linked axes (`X_RANGE` × `Y_RANGE`) at equal aspect
//! with the same ticks (`X_TICKS`, `Y_TICKS`), titled by the parameter it varies, aligned with its
//! plot area's left edge. The parameters it holds fixed title its legend: a panel title in full
//! ("position of maximum camber p (m = 0.04, t = 0.12)") is wider than a column at ironlab's
//! default font, and ironlab sets text on one line.
//! The members are drawn 0.75 pt in `figures::ylgnbu(i / (n − 1))` (already ordered by parameter
//! value), the family default black 1.5 pt on top. Each legend lists the default first, on a row of
//! its own, then the members in as many columns as fit the panel's width, in a row fitted to them.

use std::path::{Path, PathBuf};

use ironlab::prelude::*;
use serde_json::Value;
use yfoil::output::figures;

use crate::{X_RANGE, Y_RANGE};

/// Columns of panels
const COLUMNS: usize = 2;
/// The ticks of every panel
const X_TICKS: [f64; 8] = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2];
const Y_TICKS: [f64; 5] = [-0.4, -0.2, 0.0, 0.2, 0.4];
/// The nominal height (mm) of a row of panels with its legend. The legend row is fitted to its
/// entries and the plot takes the rest; equal aspect then sizes each plot area within it.
const ROW_MM: f64 = 95.0;
/// The figure title (mm)
const TITLE_MM: f64 = 10.0;

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

/// A panel title split into the parameter it varies and the ones it holds fixed:
/// "maximum camber m (p = 0.4, t = 0.12)" → ("maximum camber m", "p = 0.4, t = 0.12")
fn split_title(title: &str) -> (&str, Option<&str>) {
    match title.rsplit_once(" (") {
        Some((varied, held)) if held.ends_with(')') => (varied, Some(&held[..held.len() - 1])),
        _ => (title, None),
    }
}

fn members(panel: &Value) -> &Vec<Value> {
    panel["members"].as_array().expect("summary.json: members")
}

/// The family sheet of one `summary.json` series entry
pub fn family_sheet(run_dir: &Path, series: &Value, n_nodes: i64) -> Figure {
    let slug = str_of(series, "slug");
    let panels = series["panels"].as_array().expect("summary.json: panels");
    let default = &series["default"];
    let rows = panels.len().div_ceil(COLUMNS);

    let mut fig = Figure::new()
        .size_mm(figures::WIDTH_MM, rows as f64 * ROW_MM + TITLE_MM)
        .tiles((2 * rows) as u32, COLUMNS as u32)
        .title(Text::plain(str_of(series, "title")))
        .tag("validation")
        .parameter("study", "aerofoil-series")
        .parameter("figure", "family-sheet")
        .parameter("series", slug)
        .parameter("n_nodes", n_nodes)
        .parameter("designation", str_of(default, "designation"))
        .parameter("max_thickness", f64_of(default, "max_thickness"))
        .parameter("te_gap", f64_of(default, "te_gap"));
    for r in 0..rows {
        fig = fig.fit_row((2 * r + 1) as u32);
    }

    let (default_x, default_y) = coordinates(&run_dir.join(slug).join(str_of(default, "file")));
    let default_label = format!("{} (default)", str_of(default, "label"));
    let mut axes = Vec::with_capacity(panels.len());
    for (k, panel) in panels.iter().enumerate() {
        let (r, c) = ((k / COLUMNS) as u32, (k % COLUMNS) as u32);
        let members = members(panel);
        let span = (members.len().max(2) - 1) as f64;
        let mut ax = fig.axes(2 * r, c);
        let mut ids = Vec::with_capacity(members.len());
        for (i, m) in members.iter().enumerate() {
            let (x, y) = coordinates(&run_dir.join(slug).join(str_of(panel, "slug")).join(str_of(m, "file")));
            let id = ax
                .plot(&x, &y)
                .label(Text::plain(str_of(m, "label")))
                .color(figures::ylgnbu(i as f64 / span))
                .line_width(0.75)
                .id();
            ids.push(id);
        }
        // the default last, so that it is drawn on top
        let default_id = ax
            .plot(&default_x, &default_y)
            .label(Text::plain(default_label.clone()))
            .color(Color::rgb(0.0, 0.0, 0.0))
            .line_width(1.5)
            .id();
        let (varied, held) = split_title(str_of(panel, "title"));
        ax.title(Text::plain(varied))
            .title_align(TitleAlign::Left)
            .xlabel("$x/c$")
            .ylabel("$y/c$")
            .xlim(X_RANGE.0, X_RANGE.1)
            .ylim(Y_RANGE.0, Y_RANGE.1)
            .xticks(X_TICKS)
            .yticks(Y_TICKS)
            .aspect(Aspect::Equal)
            .aspect_adjust(AspectAdjust::Box);
        axes.push(ax.id());
        // listed first, on a row of its own, then the members
        let key = {
            let mut panel = fig.panel(2 * r + 1, c);
            panel.align(PanelLocation::NorthWest);
            if let Some(held) = held {
                panel.title(Text::plain(held)).title_align(TitleAlign::Left);
            }
            let mut legend = panel.legend();
            legend.full_row(default_id);
            for id in ids {
                legend.item(id);
            }
            legend
                .columns(LegendEntryColumns::Auto)
                .order(LegendEntryOrder::RowMajor)
                .boxed(false);
            legend.id()
        };
        fig.axes(2 * r, c).legend_in(key);
    }
    // the same x/c and y/c on every panel
    fig.link(Dim::X, &axes).expect("the panels belong to this figure");
    fig.link(Dim::Y, &axes).expect("the panels belong to this figure");
    fig
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
            figures::save_with_pdf(&family_sheet(run_dir, s, n_nodes), &stem).unwrap_or_else(|e| panic!("{e}"));
            stem
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn titles_split_into_the_varied_and_the_held_parameters() {
        assert_eq!(
            super::split_title("position of maximum camber p (m = 0.04, t = 0.12)"),
            ("position of maximum camber p", Some("m = 0.04, t = 0.12"))
        );
        assert_eq!(super::split_title("family"), ("family", None));
    }

    use ironlab::ir::Parameter;
    use ironlab::prelude::*;

    /// Every family's sheet is in the documentation, loads, names its series, and is embedded
    #[test]
    fn documentation_holds_every_family_sheet() {
        let docs = crate::repo_root().join("docs/validation/aerofoil-series/figures");
        let page =
            std::fs::read_to_string(crate::repo_root().join("docs/validation/aerofoil-series/generator-families.md"))
                .unwrap();
        for s in crate::study::series() {
            let file = format!("{}.fig", s.slug);
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
                Some(&Parameter::String("aerofoil-series".into())),
                "{file}"
            );
            assert_eq!(
                p.get("designation"),
                Some(&Parameter::String(s.default.designation.clone())),
                "{file}"
            );
            assert!(fig.tags().iter().any(|l| l == "validation"), "{file}");
            assert!(fig.validate().errors.is_empty(), "{file}");
            assert_eq!(fig.ir().axes.len(), s.panels.len(), "{file}: one panel per parameter");
            assert!(
                page.contains(&format!("<ironlab-figure src=\"figures/{file}\"")),
                "{file} not embedded"
            );
        }
    }
}
