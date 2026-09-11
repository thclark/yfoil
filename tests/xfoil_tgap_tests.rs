//! GDES `TGAP`, reproduced: XFOIL LOADs yFoil's panels, the instrumented TGAP dumps the buffer
//! airfoil with its spline slopes, the leading-edge point and the gap direction before the
//! change, and the moved coordinates after it (`xfoil_tgap.dat`, ES24.16). The "before" block
//! must be the LOADed panels bitwise (the Rule 4 handoff, seen from inside GDES); `set_te_gap`
//! must reproduce the "after" block within `TOL_PURE`, the spline-derived quantities (SBLE,
//! XBLE, YBLE, the sharp-TE direction) likewise. Two tracked cases: a closed 6-series trailing
//! edge (the direction comes from the end slopes) and the blunt 4-digit one (the direction is
//! the existing gap).

mod fixtures;
mod utilities;

use std::collections::HashMap;
use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{arc_coordinate, find_le, set_te_gap, spline_segmented, spline_value, Geometry};

struct TgapDump {
    header: HashMap<String, f64>,
    before: Vec<[f64; 5]>,
    after: Vec<[f64; 2]>,
}

fn load(case: &str) -> (Geometry, TgapDump, [f64; 2]) {
    let dir = format!("tests/fixtures/xfoil/{case}");
    let panels: Geometry = serde_json::from_str(
        &std::fs::read_to_string(fixtures::require_fixture(&format!("{dir}/panels.json"))).unwrap(),
    )
    .unwrap();
    let manifest: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(fixtures::require_fixture(&format!("{dir}/manifest.json"))).unwrap(),
    )
    .unwrap();
    let tgap = manifest["case"]["tgap"].as_array().unwrap();
    let args = [tgap[0].as_f64().unwrap(), tgap[1].as_f64().unwrap()];
    let mut d = TgapDump {
        header: HashMap::new(),
        before: Vec::new(),
        after: Vec::new(),
    };
    let text = std::fs::read_to_string(fixtures::require_fixture(&format!("{dir}/xfoil_tgap.dat"))).unwrap();
    for l in text.lines() {
        let nums = |r: &str| -> Vec<f64> { r.split_whitespace().skip(1).map(|t| t.parse().unwrap()).collect() };
        if let Some(r) = l.strip_prefix('B') {
            let v = nums(r);
            d.before.push([v[0], v[1], v[2], v[3], v[4]]);
        } else if let Some(r) = l.strip_prefix('A') {
            let v = nums(r);
            d.after.push([v[0], v[1]]);
        } else if let Some((k, v)) = l.split_once('=') {
            if let Ok(f) = v.trim().parse::<f64>() {
                d.header.insert(k.trim().to_string(), f);
            }
        }
    }
    (panels, d, args)
}

fn check(case: &str) {
    let (panels, d, [gap, blend]) = load(case);
    let nb = d.header["NB"] as usize;
    assert_eq!(nb, panels.x.len(), "{case}: NB");
    assert_eq!(d.before.len(), nb, "{case}: before rows");
    assert_eq!(d.after.len(), nb, "{case}: after rows");

    // the buffer airfoil inside GDES is the LOADed file, bitwise
    for (i, r) in d.before.iter().enumerate() {
        assert_eq!(r[0].to_bits(), panels.x[i].to_bits(), "{case}: XB({}) handoff", i + 1);
        assert_eq!(r[1].to_bits(), panels.y[i].to_bits(), "{case}: YB({}) handoff", i + 1);
    }

    // the spline-derived inputs of TGAP
    let sb = arc_coordinate(&panels.x, &panels.y);
    let xbp = spline_segmented(&panels.x, &sb);
    let ybp = spline_segmented(&panels.y, &sb);
    for (i, r) in d.before.iter().enumerate() {
        assert_within(sb[i], r[2], TOL_PURE, 1.0, &format!("{case}: SB({})", i + 1));
        assert_within(xbp[i], r[3], TOL_PURE, 1.0, &format!("{case}: XBP({})", i + 1));
        assert_within(ybp[i], r[4], TOL_PURE, 1.0, &format!("{case}: YBP({})", i + 1));
    }
    let sble = find_le(&panels.x, &xbp, &panels.y, &ybp, &sb);
    assert_within(sble, d.header["SBLE"], TOL_PURE, 1.0, &format!("{case}: SBLE"));
    assert_within(
        spline_value(sble, &panels.x, &xbp, &sb),
        d.header["XBLE"],
        TOL_PURE,
        1.0,
        &format!("{case}: XBLE"),
    );
    assert_within(
        spline_value(sble, &panels.y, &ybp, &sb),
        d.header["YBLE"],
        TOL_PURE,
        1.0,
        &format!("{case}: YBLE"),
    );
    assert_eq!(d.header["GAPNEW"], gap, "{case}: GAPNEW is the case's gap");
    assert_eq!(d.header["DOC"], blend, "{case}: DOC is the case's blend");

    // the moved buffer airfoil
    let moved = set_te_gap(&panels, gap, blend);
    for (i, r) in d.after.iter().enumerate() {
        assert_within(
            moved.x[i],
            r[0],
            TOL_PURE,
            1.0,
            &format!("{case}: XB({}) after TGAP", i + 1),
        );
        assert_within(
            moved.y[i],
            r[1],
            TOL_PURE,
            1.0,
            &format!("{case}: YB({}) after TGAP", i + 1),
        );
    }
    // The resulting gap is the requested one only when a gap existed: on a sharp trailing edge
    // XFOIL's direction is the mean of the two end tangents, whose length is cos(half the TE
    // angle) times the spline-parameter stretch, not 1 (docs/xfoil-known-issues.md, TGAP). The
    // ratio is reported; the nodes above are the gate.
    let new_gap = (moved.x[0] - moved.x[nb - 1]).hypot(moved.y[0] - moved.y[nb - 1]);
    if d.header["GAP"] > 0.0 {
        assert_within(new_gap, gap, TOL_PURE, 1.0, &format!("{case}: resulting gap"));
    } else {
        assert!(
            (new_gap / gap - 1.0).abs() < 1e-2,
            "{case}: resulting gap {new_gap:e} for requested {gap:e}"
        );
    }
    let bitwise = d
        .after
        .iter()
        .enumerate()
        .filter(|(i, r)| r[0].to_bits() == moved.x[*i].to_bits() && r[1].to_bits() == moved.y[*i].to_bits())
        .count();
    println!(
        "{case}: {bitwise} / {nb} moved nodes bitwise identical, gap {:.3e} → {:.3e}",
        d.header["GAP"], new_gap
    );
}

#[test]
fn tgap_on_a_closed_six_series_trailing_edge() {
    // GAP = 0 exactly: the direction is the mean of the end slopes (XFOIL's sharp branch)
    let (_, d, _) = load("naca63-415_n160_tgap");
    assert_eq!(d.header["GAP"], 0.0, "the 6-series trailing edge closes exactly");
    check("naca63-415_n160_tgap");
}

#[test]
fn tgap_on_the_blunt_four_digit_trailing_edge() {
    let (_, d, _) = load("naca0012_n60_tgap");
    assert!(d.header["GAP"] > 0.0, "the 4-digit trailing edge is open");
    check("naca0012_n60_tgap");
}
