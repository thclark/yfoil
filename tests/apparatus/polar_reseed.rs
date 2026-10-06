//! (f) Reference integrity: the reference agrees with itself: in the polar case (`ALFA 0 / ASEQ 1 5 1 / INIT / ALFA 0 /
//! ASEQ -1 -5 -1`) the re-solved 0° after `INIT` (VISCAL call 7) is the first 0° again — textually
//! identical ES24.16 records. That is what makes it a seed for the downward leg rather than a
//! second polar point, and what lets yFoil's content-addressed ids cite the first solve.

use crate::fixtures;
use std::collections::HashMap;

fn parse_points(path: &std::path::Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out: Vec<HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "CALL" {
            out.push(HashMap::new());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

#[test]
fn test_xfoil_resolved_zero_after_init_is_the_first_zero() {
    let xf = parse_points(&fixtures::require_fixture(
        "tests/fixtures/xfoil/naca0012_n60_polar_re1e6/viscal_points.dat",
    ));
    assert_eq!(xf.len(), 12, "expected 12 VISCAL calls (0..5, INIT, 0, -1..-5)");
    for key in [
        "NITDONE", "LVCONV", "RMSBL", "CL", "CM", "CD", "CDF", "CDP", "CL_ALF", "IST", "ITRAN1", "ITRAN2", "XOCTR1",
        "XOCTR2",
    ] {
        assert_eq!(xf[0][key], xf[6][key], "XFOIL re-solved 0° after INIT: {key}");
    }
}
