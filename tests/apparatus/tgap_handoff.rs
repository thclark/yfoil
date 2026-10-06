//! (f) Reference integrity: the Rule 4 handoff seen from inside GDES: the buffer airfoil XFOIL's TGAP starts from is the
//! LOADed panel file, bitwise, on both `*_tgap` cases.
//!
//! Fixtures: `tests/fixtures/xfoil/*_tgap/` — `cargo xtask fixtures --case <name>`.

use crate::fixtures::tgap_fixtures::load;

fn check(case: &str) {
    let (panels, d, _) = load(case);
    assert_eq!(d.before.len(), panels.x.len(), "{case}: before rows");
    for (i, r) in d.before.iter().enumerate() {
        assert_eq!(r[0].to_bits(), panels.x[i].to_bits(), "{case}: XB({}) handoff", i + 1);
        assert_eq!(r[1].to_bits(), panels.y[i].to_bits(), "{case}: YB({}) handoff", i + 1);
    }
}

#[test]
fn tgap_buffer_is_the_loaded_six_series_panels() {
    check("naca63-415_n160_tgap");
}

#[test]
fn tgap_buffer_is_the_loaded_four_digit_panels() {
    check("naca0012_n60_tgap");
}
