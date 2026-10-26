//! The study's two occupancy figures, built from `summary.json` as ironlab figures
//! (`docs/development/conventions/figures.md`): the branch map (every branch of the analysis path,
//! one row per subroutine grouped by source file, coloured by status) and the case × subroutine
//! matrix (the fraction of each subroutine's branches each case of the minimal set takes).
//!
//! The matplotlib case matrix carried three columns of counts beside its rows (branches,
//! reachable, taken over all cases); ironlab has no margin table yet (thclark/ironlab#121), so the
//! counts are the `subroutine-counts` table the pages include under the figure.

use std::collections::BTreeSet;
use std::path::Path;

use ironlab::prelude::*;
use serde_json::Value;
use yfoil::output::figures::{self, WIDTH_MM};

/// The branch statuses as the figures key them, in the categorical scale's order. The DO-line
/// edges that are never taken (`loop_entry`) are not reachable either, so they share that entry.
const STATUSES: [(&str, (u8, u8, u8), &str); 5] = [
    ("cover", (0x2c, 0xa0, 0x2c), "taken by the finite set"),
    ("nonfinite", (0x1f, 0x77, 0xb4), "taken by the non-finite set"),
    ("pathological", (0xff, 0x7f, 0x0e), "taken only in pathological cases"),
    ("unreachable", (0x9a, 0x9a, 0x9a), "not reachable"),
    ("open", (0xd6, 0x27, 0x28), "not taken"),
];

/// The colour of a cell whose subroutine the case never entered.
const NOT_ENTERED: Color = Color::rgb(0xf4 as f32 / 255.0, 0xf4 as f32 / 255.0, 0xf4 as f32 / 255.0);

/// The height of a row of either figure (mm), set by its tick label.
const ROW_MM: f64 = 3.4;

/// The number of colours of the listed YlGnBu colormap of the case matrix.
const RAMP: usize = 64;

fn rgb((r, g, b): (u8, u8, u8)) -> Color {
    Color::rgb(r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
}

/// A name of code in monospace: LaTeX's reserved `_` escaped.
fn tt(name: &str) -> String {
    format!(r"\texttt{{{}}}", name.replace('_', r"\_"))
}

/// The index of a status in [`STATUSES`].
fn status_index(status: &str) -> f64 {
    let status = if status == "loop_entry" { "unreachable" } else { status };
    STATUSES
        .iter()
        .position(|(key, _, _)| *key == status)
        .unwrap_or_else(|| panic!("unknown branch status {status}")) as f64
}

/// The rows of `branch_map` that have branches, in the study's (coverage.toml) order.
fn rows(s: &Value) -> Vec<&Value> {
    s["branch_map"]
        .as_array()
        .expect("branch_map")
        .iter()
        .filter(|r| r["statuses"].as_array().is_some_and(|a| !a.is_empty()))
        .collect()
}

fn str_of<'a>(v: &'a Value, key: &str) -> &'a str {
    v[key].as_str().unwrap_or_else(|| panic!("{key} in summary.json"))
}

/// The branch map: one row of cells per subroutine of the analysis path, under a bold header per
/// source file, each cell one branch coloured by its status; the key is the status scale.
pub fn branch_map(s: &Value) -> Figure {
    let rows = rows(s);
    let mut files: Vec<(&str, Vec<String>)> = Vec::new();
    for r in &rows {
        let (file, name) = (str_of(r, "file"), str_of(r, "name"));
        match files.last_mut() {
            Some((f, members)) if *f == file => members.push(tt(name)),
            _ => files.push((file, vec![tt(name)])),
        }
    }
    let statuses: Vec<Vec<f64>> = rows
        .iter()
        .map(|r| {
            r["statuses"]
                .as_array()
                .unwrap()
                .iter()
                .map(|st| status_index(st.as_str().unwrap()))
                .collect()
        })
        .collect();
    let branches: usize = statuses.iter().map(Vec::len).sum();
    let lines = rows.len() + 2 * files.len();
    // a row is as tall as its label needs and the widest row spans the width the labels leave, so
    // the cells are not square: 110 square cells to a row would make the rows too thin to label
    let height = 40.0 + ROW_MM * lines as f64;

    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height)
        .tiles(2, 1)
        .fit_row(1)
        .title("Branches of the analysis path by status")
        .parameter("study", "branch-coverage")
        .parameter("figure", "branch-map")
        .parameter("subroutines", rows.len() as i64)
        .parameter("branches", branches as i64)
        .tag("validation");
    let status = fig
        .categorical_scale(STATUSES.map(|(_, c, label)| (rgb(c), label)))
        .id();
    let mut ax = fig.axes(0, 0);
    let placed = ax.ycategory_groups(files.iter().map(|(f, m)| (*f, m.iter())));
    ax.occupancy(placed.place(&statuses)).color_scale(status);
    ax.box_on(false);
    fig.panel(1, 0)
        .legend()
        .item(status)
        .columns(LegendEntryColumns::Count(3));
    fig
}

/// The case × subroutine matrix: the fraction of each subroutine's branches each case of the
/// minimal set takes; a dot marks a subroutine in which the case takes a branch no other finite
/// candidate takes, and a subroutine the case never entered is keyed as such.
pub fn case_matrix(s: &Value) -> Figure {
    let m = &s["matrix"];
    let cases: Vec<&str> = m["cases"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    let all_subs: Vec<&str> = m["subroutines"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    let keep: Vec<usize> = (0..all_subs.len())
        .filter(|&i| s["branch_map"][i]["statuses"].as_array().is_some_and(|a| !a.is_empty()))
        .collect();
    let subs: Vec<&str> = keep.iter().map(|&i| all_subs[i]).collect();
    let name = |key: &str| key.split_once(':').map_or(key, |(_, n)| n).to_string();
    // a case that never entered a subroutine took none of its branches: every subroutine with a
    // conditional takes at least one edge once entered, so a fraction of zero is "not entered"
    let fraction = Matrix::from_fn(subs.len(), cases.len(), |row, col| {
        let f = m["fraction"][col][keep[row]].as_f64().unwrap();
        if f > 0.0 {
            f
        } else {
            f64::NAN
        }
    });
    let mut unique = (Vec::new(), Vec::new());
    for (col, case) in cases.iter().enumerate() {
        let c = s["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["name"].as_str() == Some(case))
            .unwrap_or_else(|| panic!("case {case} in summary.json"));
        let subs_unique: BTreeSet<&str> = c["unique_by_subroutine"]
            .as_object()
            .map(|o| o.keys().map(String::as_str).collect())
            .unwrap_or_default();
        for (row, key) in subs.iter().enumerate() {
            if subs_unique.contains(name(key).as_str()) {
                unique.0.push(col as f64);
                unique.1.push(row as f64);
            }
        }
    }

    let height = 85.0 + ROW_MM * subs.len() as f64;
    let mut fig = Figure::new()
        .size_mm(WIDTH_MM, height)
        .tiles(3, 1)
        .fit_row(1)
        .fit_row(2)
        .title("Branches each case of the minimal finite set takes")
        .parameter("study", "branch-coverage")
        .parameter("figure", "case-matrix")
        .parameter("cases", cases.len() as i64)
        .parameter("subroutines", subs.len() as i64)
        .tag("validation");
    let ramp: Vec<Color> = (0..RAMP)
        .map(|k| figures::ylgnbu(k as f64 / (RAMP - 1) as f64))
        .collect();
    let taken = fig
        .color_scale()
        .colormap(ironlab::ir::Colormap::Listed { colors: ramp })
        .limits(0.0, 1.0)
        .tick_step(0.25)
        .non_finite(NOT_ENTERED)
        .key_label(ColorScaleKey::NonFinite, "subroutine not entered by the case")
        .show_key(ColorScaleKey::NonFinite)
        .title("Fraction of the subroutine's branches the case took")
        .id();
    let key = fig.panel(2, 0).legend().boxed(false).id();
    let mut ax = fig.axes(0, 0);
    ax.cells(&fraction)
        .color_scale(taken)
        .edge_color(Color::WHITE)
        .edge_width(0.4);
    ax.scatter(&unique.0, &unique.1)
        .label("takes a branch no other finite candidate takes")
        .marker(Marker::Circle)
        .marker_size(4.0)
        .marker_face(Color::WHITE)
        .marker_edge(Color::BLACK);
    ax.xcategories(cases.iter().map(|c| tt(c)))
        .xtick_label_angle(90.0)
        .ycategories(subs.iter().map(|k| tt(&name(k))))
        .box_on(false)
        .legend_in(key);
    fig.panel(1, 0)
        .colorbar(taken)
        .orientation(ColorbarOrientation::Horizontal)
        .title_placement(ColorbarTitlePlacement::Above)
        .length(0.8);
    fig
}

/// The per-subroutine counts the matplotlib case matrix carried beside its rows: every branch,
/// the reachable ones (not annotated unreachable, not a never-taken DO-line edge) and those taken
/// on any measured run (the finite set, the non-finite set and the measured pathological cases).
pub fn subroutine_counts_md(s: &Value) -> String {
    let mut t =
        String::from("| file | subroutine | branches | reachable | taken, all cases |\n|---|---|---:|---:|---:|\n");
    for r in s["subroutines"].as_array().unwrap() {
        let n = |k: &str| r[k].as_u64().unwrap_or(0);
        if n("branches") == 0 {
            continue;
        }
        t += &format!(
            "| `{}` | `{}` | {} | {} | {} |\n",
            str_of(r, "file"),
            str_of(r, "name"),
            n("branches"),
            n("branches") - n("unreachable") - n("loop_entry"),
            n("cover") + n("nonfinite") + n("pathological"),
        );
    }
    t
}

/// Draw both figures into the run folder (`.fig` and `.pdf`). A figure that does not validate or
/// lay out is reported and not written, and the publish check then fails on the page that embeds
/// it.
pub fn render(run_dir: &Path, s: &Value) -> bool {
    let mut ok = true;
    for (name, fig) in [("branch-map", branch_map(s)), ("case-matrix", case_matrix(s))] {
        if let Err(e) = figures::save_with_pdf(&fig, &run_dir.join(name)) {
            eprintln!("figure not drawn: {e}");
            ok = false;
        }
    }
    ok
}
