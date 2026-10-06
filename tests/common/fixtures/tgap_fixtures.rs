//! The TGAP dump (`xfoil_tgap.dat`) of the `*_tgap` cases: the buffer airfoil before GDES `TGAP`
//! with its spline slopes, the leading-edge point and gap direction, and the moved coordinates.

use std::collections::HashMap;
use yfoil::geometry::Geometry;

/// Header scalars, the "before" rows (X Y SB XBP YBP) and the "after" rows (X Y).
pub struct TgapDump {
    pub header: HashMap<String, f64>,
    pub before: Vec<[f64; 5]>,
    pub after: Vec<[f64; 2]>,
}

/// The case's panels, its dump, and its `[gap, blend]` arguments.
pub fn load(case: &str) -> (Geometry, TgapDump, [f64; 2]) {
    let dir = format!("tests/fixtures/xfoil/{case}");
    let panels: Geometry =
        serde_json::from_str(&std::fs::read_to_string(super::require_fixture(&format!("{dir}/panels.json"))).unwrap())
            .unwrap();
    let manifest: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(super::require_fixture(&format!("{dir}/manifest.json"))).unwrap(),
    )
    .unwrap();
    let tgap = manifest["case"]["tgap"].as_array().unwrap();
    let args = [tgap[0].as_f64().unwrap(), tgap[1].as_f64().unwrap()];
    let mut d = TgapDump {
        header: HashMap::new(),
        before: Vec::new(),
        after: Vec::new(),
    };
    let text = std::fs::read_to_string(super::require_fixture(&format!("{dir}/xfoil_tgap.dat"))).unwrap();
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
