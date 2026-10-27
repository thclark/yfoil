//! The run folder's tables (Markdown fragments in `publish/tables/`) and its `README.md`.

use crate::analysis::{
    decades, flagged, gap, groups, histogram, ill_conditioned, may_be_ill_conditioned, measured, plateau_thresholds,
    points, sci, thresholds, worst_spreads, worst_spreads_of, DECADES,
};
use crate::drive::val;
use crate::utilities::publish::Bundle;
use serde_json::Value;
use std::fmt::Write;
use std::path::Path;

/// The thresholds under study, (name, value), warning first.
pub fn named_thresholds(summary: &Value) -> [(&'static str, f64); 2] {
    let t = thresholds(summary);
    [("warning", t.warning), ("error", t.error)]
}

/// ` ← warning threshold` when `t` is one of the thresholds (to rounding), else empty.
fn threshold_mark(summary: &Value, t: f64) -> String {
    named_thresholds(summary)
        .iter()
        .filter(|(_, v)| (t - v).abs() <= 1e-9 * v)
        .map(|(n, _)| format!(" ← {n} threshold"))
        .collect()
}

fn decade_label(i: usize) -> String {
    let k = *DECADES.start() + i as i32;
    if i == 0 {
        format!("< 1e{}", k + 1)
    } else if k > *DECADES.end() {
        "≥ 1".into()
    } else {
        format!("[1e{k}, 1e{})", k + 1)
    }
}

pub fn histogram_md(s: &Value) -> String {
    let gs = groups(s);
    let all = worst_spreads(s);
    let (h, nan) = histogram(&all);
    let per: Vec<(Vec<usize>, usize, usize)> = gs
        .iter()
        .map(|g| {
            let w = worst_spreads_of(s, Some(g));
            let (h, nan) = histogram(&w);
            (h, nan, w.len())
        })
        .collect();
    let mut o = String::from("| Largest spread |");
    for (g, (_, _, n)) in gs.iter().zip(&per) {
        write!(o, " {g} ({n}) |").unwrap();
    }
    writeln!(o, " all ({}) |", all.len()).unwrap();
    o.push_str(&format!("|---|{}---:|\n", "---:|".repeat(gs.len())));
    // a threshold of 1e-k is the bottom edge of the decade [1e-k, 1e-k+1)
    for (i, n) in h.iter().enumerate() {
        let k = *DECADES.start() + i as i32;
        let mark: String = named_thresholds(s)
            .iter()
            .filter(|(_, v)| cr_libm::log10(*v).floor() as i32 == k)
            .map(|(name, v)| format!(" ← {name} threshold {}", sci(*v)))
            .collect();
        if *n == 0 && mark.is_empty() {
            continue;
        }
        write!(o, "| {}{mark} |", decade_label(i)).unwrap();
        for (hg, _, _) in &per {
            write!(o, " {} |", hg[i]).unwrap();
        }
        writeln!(o, " {n} |").unwrap();
    }
    write!(o, "| not a number |").unwrap();
    for (_, ng, _) in &per {
        write!(o, " {ng} |").unwrap();
    }
    writeln!(o, " {nan} |").unwrap();
    o
}

pub fn plateau_md(s: &Value) -> String {
    let w = worst_spreads(s);
    let mut o = String::from("| Threshold | Points whose largest spread exceeds it |\n|---:|---:|\n");
    for t in plateau_thresholds() {
        let f = w.iter().filter(|&&p| flagged(p, t)).count();
        writeln!(o, "| {}{} | {f} of {} |", sci(t), threshold_mark(s, t), w.len()).unwrap();
    }
    o
}

pub fn gap_md(s: &Value) -> String {
    let mut o = String::from(
        "| Threshold | Cases | Points measured | Beyond it | Largest spread at or below | Smallest spread above | Same points beyond it for any threshold in | Decades of clearance below, above |\n\
         |---|---|---:|---:|---:|---:|---|---|\n",
    );
    let gs = groups(s);
    let rows: Vec<(String, Vec<f64>)> = gs
        .iter()
        .map(|g| (g.clone(), worst_spreads_of(s, Some(g))))
        .chain([("**all**".to_string(), worst_spreads(s))])
        .collect();
    for (tname, t) in named_thresholds(s) {
        for (name, w) in &rows {
            let g = gap(w.iter().copied(), t);
            let (lo, hi) = g.interval();
            let (below, above) = g.clearance(t);
            writeln!(
                o,
                "| {tname} {} | {name} | {} | {} ({} not a number) | {} | {} | [{}, {}) | {}, {} |",
                sci(t),
                g.n,
                g.n_above,
                g.n_nan,
                g.below.map_or("—".into(), sci),
                g.above.map_or("—".into(), sci),
                sci(lo),
                sci(hi),
                decades(below),
                decades(above),
            )
            .unwrap();
        }
    }
    o
}

pub fn gap_by_quantity_md(s: &Value) -> String {
    let t = thresholds(s).error;
    let names = [
        "alpha",
        "cl",
        "cm",
        "cd_pressure",
        "cd",
        "cd_friction",
        "transition_upper",
        "transition_lower",
    ];
    let mut o = format!(
        "| Quantity | Points reporting it | Largest spread, points within the error threshold | Decades below the error threshold {} | Ill-conditioned points this quantity alone flags | Its spread over the ill-conditioned points |\n|---|---:|---:|---:|---:|---|\n",
        sci(t)
    );
    let measured: Vec<&Value> = points(s).filter(|p| measured(p)).collect();
    let n_ill = measured.iter().filter(|p| ill_conditioned(p, s)).count();
    for n in names {
        let (mut reporting, mut well, mut flags) = (0usize, 0.0_f64, 0usize);
        let (mut lo, mut hi) = (f64::INFINITY, 0.0_f64);
        let mut ill_nan = 0usize;
        for p in &measured {
            let Some(v) = p["spread"].get(n).map(val) else { continue };
            reporting += 1;
            if ill_conditioned(p, s) {
                if flagged(v, t) {
                    flags += 1;
                }
                if v.is_nan() {
                    ill_nan += 1;
                } else {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            } else {
                well = well.max(v);
            }
        }
        let range = if lo.is_finite() {
            format!(
                "{} … {}{}",
                sci(lo),
                sci(hi),
                if ill_nan > 0 {
                    format!(" ({ill_nan} not a number)")
                } else {
                    String::new()
                }
            )
        } else if ill_nan > 0 {
            format!("{ill_nan} not a number")
        } else {
            "—".into()
        };
        writeln!(
            o,
            "| `{n}` | {reporting} | {} | {} | {flags} of {n_ill} | {range} |",
            sci(well),
            if well > 0.0 {
                decades(cr_libm::log10(t / well))
            } else {
                "∞".into()
            }
        )
        .unwrap();
    }
    o
}

fn operating_point(p: &Value) -> String {
    match p["command"].as_str() {
        Some("cl") => format!("CL {} (α {:.2}°)", p["input"], val(&p["alpha_deg"])),
        _ => format!("α {}°", p["input"]),
    }
}

/// The points `pick` selects, one row each.
fn point_list(s: &Value, pick: fn(&Value, &Value) -> bool) -> String {
    let mut o = String::from(
        "| Group | Case | Leg | Operating point | Worst quantity | Spread | Twins not valid |\n|---|---|---:|---|---|---:|---:|\n",
    );
    let mut n = 0;
    for c in s["cases"].as_array().unwrap() {
        for p in c["points"].as_array().unwrap() {
            if pick(p, s) {
                n += 1;
                writeln!(
                    o,
                    "| {} | `{}` | {} | {} | `{}` | {} | {} of 5 |",
                    c["group"].as_str().unwrap_or(""),
                    c["name"].as_str().unwrap(),
                    p["leg"],
                    operating_point(p),
                    p["worst"].as_str().unwrap_or(""),
                    sci(val(&p["worst_spread"])),
                    p["twins_not_valid"]
                )
                .unwrap();
            }
        }
    }
    if n == 0 {
        o.push_str("| — | | | none | | | |\n");
    }
    o
}

pub fn ill_conditioned_md(s: &Value) -> String {
    point_list(s, ill_conditioned)
}

pub fn may_be_ill_conditioned_md(s: &Value) -> String {
    point_list(s, may_be_ill_conditioned)
}

pub fn cases_md(s: &Value) -> String {
    let mut o = String::from(
        "| Group | Case | Points solved | Invalid on a single-run reason (not measured) | Measured | May be ill-conditioned (warning) | Ill-conditioned (withheld) |\n|---|---|---:|---:|---:|---:|---:|\n",
    );
    let (mut tp, mut ti, mut tm, mut tw, mut tc) = (0, 0, 0, 0, 0);
    let mut per_group: Vec<(String, usize, [usize; 5])> = groups(s).into_iter().map(|g| (g, 0, [0; 5])).collect();
    for c in s["cases"].as_array().unwrap() {
        let pts = c["points"].as_array().unwrap();
        let invalid = pts
            .iter()
            .filter(|p| !p["single_run_reasons"].as_array().unwrap().is_empty())
            .count();
        let n_measured = pts.iter().filter(|p| measured(p)).count();
        let warned = pts.iter().filter(|p| may_be_ill_conditioned(p, s)).count();
        let ill = pts.iter().filter(|p| ill_conditioned(p, s)).count();
        (tp, ti, tm, tw, tc) = (tp + pts.len(), ti + invalid, tm + n_measured, tw + warned, tc + ill);
        let group = c["group"].as_str().unwrap_or("");
        if let Some((_, n, t)) = per_group.iter_mut().find(|(g, _, _)| g == group) {
            *n += 1;
            for (a, b) in t.iter_mut().zip([pts.len(), invalid, n_measured, warned, ill]) {
                *a += b;
            }
        }
        writeln!(
            o,
            "| {group} | `{}` | {} | {invalid} | {n_measured} | {warned} | {ill} |",
            c["name"].as_str().unwrap(),
            pts.len()
        )
        .unwrap();
    }
    for (g, n, t) in &per_group {
        writeln!(
            o,
            "| **{g}** | **{n} cases** | {} | {} | {} | {} | {} |",
            t[0], t[1], t[2], t[3], t[4]
        )
        .unwrap();
    }
    writeln!(
        o,
        "| **all** | **{} cases** | {tp} | {ti} | {tm} | {tw} | {tc} |",
        s["cases"].as_array().unwrap().len()
    )
    .unwrap();
    let not_driven = s["not_driven"].as_array().cloned().unwrap_or_default();
    if !not_driven.is_empty() {
        o.push_str("\nNot driven:\n\n");
        for c in not_driven {
            writeln!(
                o,
                "- {} `{}`: {}",
                c["group"].as_str().unwrap_or(""),
                c["name"].as_str().unwrap_or(""),
                c["why"].as_str().unwrap_or("")
            )
            .unwrap();
        }
    }
    o
}

/// A table of the bundle: its file name and how it is drawn from `summary.json`.
type Table = (&'static str, fn(&Value) -> String);

pub const TABLES: [Table; 7] = [
    ("histogram", histogram_md),
    ("plateau", plateau_md),
    ("gap", gap_md),
    ("gap-by-quantity", gap_by_quantity_md),
    ("ill-conditioned", ill_conditioned_md),
    ("may-be-ill-conditioned", may_be_ill_conditioned_md),
    ("cases", cases_md),
];

pub fn write_readme(run_dir: &Path, s: &Value) {
    let mut o = String::from(
        "# Conditioning\n\n\
         The thresholds for yFoil's own twins: how far five 1-ULP jitters of the panel coordinates move each \
         point yFoil solves; where between the noise floor and a twin that reached another solution the \
         error threshold has to sit for the set of withheld points not to depend on it; and how many points \
         the warning threshold marks inside the noise floor.\n\n\
         The points are those of two groups of cases of `xtask/fixtures-config/cases.toml`: the series cases \
         (`group = \"series\"`: the aerofoil series and flow regimes of interest under normal conditions) and \
         the branch-coverage cases (`group = \"branch-coverage\"`: single operating points and short sequences \
         that take XFOIL's rarely-taken branches, viscous and inviscid, fixed CL, TYPE 2 and 3, forced \
         transition, DAMP). Each case contributes its panels and its OPER procedure; nothing the reference \
         computed is read. A case whose procedure the driver cannot express is listed under `cases`, not \
         driven. A point invalid on a single-run reason (Kármán–Tsien out \
         of domain, not converged, …) is reported invalid without its twins being compared, so it is counted \
         in `cases` but not measured.\n\n",
    );
    for (name, f) in TABLES {
        writeln!(o, "## {name}\n\n{}", f(s)).unwrap();
    }
    std::fs::write(run_dir.join("README.md"), o).unwrap();
}

pub fn write_publish(run_dir: &Path, s: &Value) {
    let meta: Value = serde_json::from_str(&std::fs::read_to_string(run_dir.join("metadata.json")).unwrap()).unwrap();
    let bundle = Bundle::create(run_dir);
    bundle.provenance(
        "cargo run --release -p conditioning",
        meta["generated"].as_str().unwrap_or(""),
        meta["git_describe"].as_str().unwrap_or(""),
        meta["git_dirty"].as_bool().unwrap_or(false),
        meta["host"].as_str().unwrap_or(""),
    );
    for (name, f) in TABLES {
        bundle.table(name, &f(s));
    }
    for f in crate::figures::NAMES {
        bundle.figure(&run_dir.join(format!("{f}.fig")));
    }
}
