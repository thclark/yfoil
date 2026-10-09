//! (d) Known XFOIL weakness: known issues §6.5: TGAP on a sharp trailing edge does not deliver the requested gap. XFOIL's
//! direction there is the mean of the two end tangents, whose length is cos(half the TE angle)
//! times the spline-parameter stretch, not 1. yFoil replicates it (the nodes are gated against
//! XFOIL in `subroutine/tgap.rs`); this defines the resulting gap on both branches.
//!
//! Fixtures: `tests/fixtures/xfoil/*_tgap/` — `cargo xtask fixtures --case <name>`.

use crate::fixtures::tgap_fixtures::load;
use crate::utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::set_te_gap;

fn resulting_gap(case: &str) -> (f64, f64, f64) {
    let (panels, d, [gap, blend]) = load(case);
    let moved = set_te_gap(&panels, gap, blend);
    let nb = moved.x.len();
    let new_gap = (moved.x[0] - moved.x[nb - 1]).hypot(moved.y[0] - moved.y[nb - 1]);
    (d.header["GAP"], gap, new_gap)
}

#[test]
fn sharp_trailing_edge_gap_misses_the_request_by_under_one_percent() {
    let (old, requested, got) = resulting_gap("naca63-415_n160_tgap");
    assert_eq!(old, 0.0);
    assert!(
        (got / requested - 1.0).abs() < 1e-2,
        "resulting gap {got:e} for requested {requested:e}"
    );
}

#[test]
fn blunt_trailing_edge_gap_is_the_request() {
    let (old, requested, got) = resulting_gap("naca0012_n60_tgap");
    assert!(old > 0.0);
    assert_within(got, requested, TOL_PURE, 1.0, "resulting gap");
}
