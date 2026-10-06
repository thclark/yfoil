//! (a) Subroutine equivalence: GDES `TGAP`, reproduced: XFOIL LOADs yFoil's panels, the instrumented TGAP dumps the buffer
//! airfoil with its spline slopes, the leading-edge point and the gap direction before the
//! change, and the moved coordinates after it (`xfoil_tgap.dat`, ES24.16). `set_te_gap` must
//! reproduce the "after" block within `TOL_PURE`, the spline-derived quantities (SBLE, XBLE, YBLE)
//! likewise. Two cases: a closed 6-series trailing edge (the direction comes from the end slopes)
//! and the blunt 4-digit one (the direction is the existing gap). The handoff is checked in
//! `apparatus/tgap_handoff.rs`, the resulting gap in `known_issues/tgap_6_5.rs`.
//!
//! Fixtures: `tests/fixtures/xfoil/*_tgap/` — `cargo xtask fixtures --case <name>`.

use crate::fixtures::tgap_fixtures::load;
use crate::utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{arc_coordinate, find_le, set_te_gap, spline_segmented, spline_value};

fn check(case: &str) {
    let (panels, d, [gap, blend]) = load(case);
    let nb = d.header["NB"] as usize;
    assert_eq!(nb, panels.x.len(), "{case}: NB");
    assert_eq!(d.before.len(), nb, "{case}: before rows");
    assert_eq!(d.after.len(), nb, "{case}: after rows");

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
