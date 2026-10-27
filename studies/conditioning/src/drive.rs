//! Each case's OPER script, driven through yFoil with its per-point twins.

use serde_json::{json, Value};
use series_cases::Case;
use std::path::{Path, PathBuf};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::output::{PointRecord, Reason, Thresholds, Warning};
use yfoil::solver::analysis::{PointResult, Session};
use yfoil::solver::Twins;

/// Settings of a `cases.toml` case that change its procedure or its inputs and that this driver
/// does not express: a case setting one is not driven. (Everything else a case may set either is
/// expressed — `alphas`, `alphas_after_reinit`, `polar`, `cls`, `inviscid`, `matyp`, `xtr`, `damp`
/// and the flow conditions — or only controls what the fixture pipeline dumps and keeps.)
pub const UNEXPRESSED: [&str; 10] = [
    "geometry_only",
    "ppar",
    "tgap",
    "stop_after_setbl",
    "stop_on_nonfinite",
    "cvpar",
    "cterat",
    "ctrrat",
    "xsref",
    "xpref",
];

/// NSEQEX: XFOIL's ASEQ halts once this many consecutive points fail to converge.
pub const NSEQEX: usize = 4;

/// One OPER command of a case's script.
#[derive(Clone, Copy)]
enum Op {
    /// `ALFA a`
    Alfa(f64),
    /// one point of `ASEQ`
    Seq(f64),
    /// `CL c`
    Cl(f64),
}

/// The case's script as `cargo xtask fixtures` writes it (`xtask/src/main.rs`): legs separated by
/// `INIT`, each a list of commands; the CL points follow in the last leg.
fn script(case: &Case) -> Vec<Vec<Op>> {
    let seq = |alphas: &[f64]| -> Vec<Op> {
        alphas
            .iter()
            .enumerate()
            .map(|(i, &a)| if i == 0 { Op::Alfa(a) } else { Op::Seq(a) })
            .collect()
    };
    let alfa = |alphas: &[f64]| -> Vec<Op> { alphas.iter().map(|&a| Op::Alfa(a)).collect() };
    let mut legs = if case.polar {
        vec![seq(&case.alphas)]
    } else {
        vec![alfa(&case.alphas)]
    };
    if !case.alphas_after_reinit.is_empty() {
        legs.push(if case.polar {
            seq(&case.alphas_after_reinit)
        } else {
            alfa(&case.alphas_after_reinit)
        });
    }
    legs.last_mut().unwrap().extend(case.cls.iter().map(|&c| Op::Cl(c)));
    legs
}

/// A JSON number, or `null` for a non-finite one.
pub fn num(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else {
        Value::Null
    }
}

/// A JSON number back, `null` (a non-finite spread) as NaN.
pub fn val(v: &Value) -> f64 {
    v.as_f64().unwrap_or(f64::NAN)
}

/// The case's panels: its tracked input (`tests/fixtures/xfoil/<case>/panels.json`) or, for a
/// case not tracked, the same file the fixture pipeline generated (`target/fixtures/<case>/`).
pub fn panels(root: &Path, case: &str) -> Option<PathBuf> {
    [
        root.join("tests/fixtures/xfoil").join(case).join("panels.json"),
        root.join("target/fixtures").join(case).join("panels.json"),
    ]
    .into_iter()
    .find(|p| p.exists())
}

/// yFoil's record of one point: its single-run reasons, and the twins' spreads when the point
/// is valid on those.
fn point(p: &PointResult, viscous: bool, leg: usize, command: &str, input: f64) -> Value {
    let r = PointRecord::from_point(p, viscous, true, Thresholds::default());
    let single_run: Vec<String> = r
        .reasons
        .iter()
        .filter(|x| **x != Reason::IllConditioned)
        .map(|x| x.to_string())
        .collect();
    let c = r.diagnostics.as_ref().and_then(|d| d.conditioning.as_ref());
    let spread = c.map(|c| {
        c.spread
            .entries()
            .into_iter()
            .map(|(n, s)| (n.to_string(), num(s)))
            .collect::<serde_json::Map<_, _>>()
    });
    json!({
        "leg": leg,
        "command": command,
        "input": input,
        "alpha_deg": num(p.alpha.to_degrees()),
        "converged": p.converged,
        "single_run_reasons": single_run,
        "measured": c.is_some(),
        "spread": spread,
        "worst": c.map(|c| c.worst.clone()),
        "worst_spread": c.map(|c| num(c.worst_spread())),
        "twins_not_valid": c.map(|c| c.twins_not_valid),
        "ill_conditioned": r.reasons.contains(&Reason::IllConditioned),
        "may_be_ill_conditioned": r.warnings.iter().any(|w| matches!(w, Warning::IllConditioned { .. })),
    })
}

/// Drive one case's script through yFoil with twins: the case's JSON record.
pub fn case(panels: &Path, case: &Case) -> Value {
    let geometry = read_geometry_from_file(panels.to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    let viscous = !case.inviscid;
    let mut session = Session::new(&airfoil, case.conditions()).with_twins(Some(Twins::new(&geometry)));
    let mut points: Vec<Value> = vec![];
    for (l, ops) in script(case).iter().enumerate() {
        if l > 0 {
            session.init();
        }
        let mut failures = 0usize;
        for op in ops {
            if failures >= NSEQEX && matches!(op, Op::Seq(_)) {
                continue; // the ASEQ halted: its remaining points are not attempted
            }
            let (p, command, input) = match *op {
                Op::Alfa(a) => (session.alpha(a.to_radians()), "alfa", a),
                Op::Seq(a) => (session.sequence_point(a.to_radians()), "aseq", a),
                Op::Cl(c) => (session.cl(c), "cl", c),
            };
            // ASEQ's halting rule on yFoil's own convergence; ALFA and CL commands do not halt
            if matches!(op, Op::Seq(_)) && viscous && !p.converged {
                failures += 1;
            } else {
                failures = 0;
            }
            points.push(point(&p, viscous, l + 1, command, input));
        }
    }
    json!({
        "name": case.name,
        "group": case.group,
        "section": case.section(),
        "script": case.script(),
        "inviscid": case.inviscid,
        "points": points,
    })
}
