//! yFoil's own arc-length cosine repanelling (`repanel_cosine`), frozen. It has no XFOIL
//! equivalent and no external reference; what it must do is keep producing exactly what it
//! produced when the tracked fixtures were derived with it. The golden files under
//! `tests/fixtures/repanel_cosine/` were written once by this test (`UPDATE_GOLDEN=1`), and every
//! node is gated at `TOL_PURE`; the node count (n + 1, its historic behaviour) is asserted too.

mod fixtures;
mod utilities;

use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{naca_4digit, repanel_cosine, Geometry, Thickness};

fn golden(name: &str, produced: &Geometry) {
    let rel = format!("tests/fixtures/repanel_cosine/{name}.json");
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&rel);
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        let mut plain = produced.clone();
        plain.generator = None;
        std::fs::write(&path, serde_json::to_string_pretty(&plain).unwrap()).unwrap();
        eprintln!("wrote {}", path.display());
    }
    let expected: Geometry =
        serde_json::from_str(&std::fs::read_to_string(fixtures::require_fixture(&rel)).unwrap()).unwrap();
    assert_eq!(produced.x.len(), expected.x.len(), "{name}: node count");
    for i in 0..expected.x.len() {
        assert_within(produced.x[i], expected.x[i], TOL_PURE, 1.0, &format!("{name}: x[{i}]"));
        assert_within(produced.y[i], expected.y[i], TOL_PURE, 1.0, &format!("{name}: y[{i}]"));
    }
}

#[test]
fn cosine_repanelling_is_frozen() {
    let source = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();
    for (bias, tag) in [(0.15, "015"), (1.0, "100"), (1.8, "180")] {
        let g = repanel_cosine(&source, 100, bias);
        assert_eq!(g.x.len(), 101, "bias {bias}: the cosine method writes n + 1 nodes");
        golden(&format!("naca0012_n100_bias{tag}"), &g);
    }
}
