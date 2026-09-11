//! Stage G gate: XFOIL's own NACA generator and PANGEN, reproduced.
//!
//! Each case runs XFOIL's `NACA dddd` then `PPAR / N n`; the instrumented PANGEN dumps the
//! 245-point buffer airfoil (XB/YB/SB), the paneling parameters and the N paneled nodes
//! (X/Y/S) at ES24.16. yFoil's `naca_*digit_vertical` must reproduce the buffer and
//! `repanel_by_curvature` the nodes, within `TOL_PURE`; the bitwise-identical counts are reported.
//! (This is the optional track: the solver equivalence never depends on it, because yFoil
//! generates the panels and XFOIL LOADs them — CLAUDE.md Rule 4.)

mod fixtures;
mod utilities;

use std::collections::HashMap;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{
    arc_coordinate, find_le, naca_4digit_vertical, naca_5digit_vertical, panel_foil, repanel_by_curvature,
    spline_segmented, PangenConfig, XFOIL_NACA_NSIDE,
};

struct PangenDump {
    header: HashMap<String, String>,
    buffer: Vec<[f64; 3]>,
    panels: Vec<[f64; 3]>,
}

fn load(case: &str) -> PangenDump {
    let p: PathBuf = fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}/xfoil_pangen.dat"));
    let mut d = PangenDump {
        header: HashMap::new(),
        buffer: Vec::new(),
        panels: Vec::new(),
    };
    for l in std::fs::read_to_string(p).unwrap().lines() {
        let row = |r: &str| -> [f64; 3] {
            let v: Vec<f64> = r
                .split_once(")=")
                .unwrap()
                .1
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            [v[0], v[1], v[2]]
        };
        if let Some(r) = l.strip_prefix("B(") {
            d.buffer.push(row(r));
        } else if let Some(r) = l.strip_prefix("P(") {
            d.panels.push(row(r));
        } else if let Some((k, v)) = l.split_once('=') {
            d.header.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    d
}

fn check(case: &str, spec: &str, npan: usize) {
    let d = load(case);
    assert_eq!(
        d.header["NB"].parse::<usize>().unwrap(),
        2 * XFOIL_NACA_NSIDE - 1,
        "{case}: NB"
    );
    assert_eq!(d.header["NPAN"].parse::<usize>().unwrap(), npan, "{case}: NPAN");
    assert_eq!(
        d.header["N"].parse::<usize>().unwrap(),
        npan,
        "{case}: N (no corners inserted)"
    );
    // the PPAR parameters XFOIL ran with, from the dump header (XREF = XSREF1 XSREF2 XPREF1 XPREF2;
    // a window of 1.0 1.0 is XFOIL's "off")
    let xref: Vec<f64> = d.header["XREF"]
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect();
    assert_eq!(xref.len(), 4, "{case}: XREF");
    let window = |a: f64, b: f64| if a == 1.0 && b == 1.0 { None } else { Some([a, b]) };
    let config = PangenConfig {
        curvature_bunching: d.header["CVPAR"].parse().unwrap(),
        te_curvature_ratio: d.header["CTERAT"].parse().unwrap(),
        refined_curvature_ratio: d.header["CTRRAT"].parse().unwrap(),
        refine_upper: window(xref[0], xref[1]),
        refine_lower: window(xref[2], xref[3]),
    };
    println!("{case}: PPAR {config:?}");

    // the buffer airfoil
    let buffer = if spec.len() == 4 {
        naca_4digit_vertical(spec).unwrap()
    } else {
        naca_5digit_vertical(spec).unwrap()
    };
    assert_eq!(buffer.x.len(), d.buffer.len(), "{case}: buffer point count");
    let sb = arc_coordinate(&buffer.x, &buffer.y);
    let (mut bits, mut worst) = (0usize, 0.0_f64);
    for (i, b) in d.buffer.iter().enumerate() {
        for (name, ours, theirs) in [
            ("XB", buffer.x[i], b[0]),
            ("YB", buffer.y[i], b[1]),
            ("SB", sb[i], b[2]),
        ] {
            if ours.to_bits() == theirs.to_bits() {
                bits += 1;
            }
            worst = worst.max((ours - theirs).abs());
            assert_within(ours, theirs, TOL_PURE, 1.0, &format!("{case}: {name}({})", i + 1));
        }
    }
    println!(
        "{case}: buffer {} points — {bits}/{} values bitwise identical, worst |diff| {worst:.2e}",
        d.buffer.len(),
        3 * d.buffer.len()
    );

    // LEFIND on the buffer (SBLE), then PANGEN
    let xbp = spline_segmented(&buffer.x, &sb);
    let ybp = spline_segmented(&buffer.y, &sb);
    let sble = find_le(&buffer.x, &xbp, &buffer.y, &ybp, &sb);
    assert_within(
        sble,
        d.header["SBLE"].parse().unwrap(),
        TOL_PURE,
        1.0,
        &format!("{case}: SBLE"),
    );

    let paneled = repanel_by_curvature(&buffer, npan, &config);
    assert_eq!(paneled.x.len(), d.panels.len(), "{case}: node count");
    let af = panel_foil(&paneled);
    let (mut bits, mut worst) = (0usize, 0.0_f64);
    for (i, p) in d.panels.iter().enumerate() {
        for (name, ours, theirs) in [("X", af.x[i], p[0]), ("Y", af.y[i], p[1]), ("S", af.s[i], p[2])] {
            if ours.to_bits() == theirs.to_bits() {
                bits += 1;
            }
            worst = worst.max((ours - theirs).abs());
            assert_within(ours, theirs, TOL_PURE, 1.0, &format!("{case}: {name}({})", i + 1));
        }
    }
    assert_within(
        af.s_le,
        d.header["SLE"].parse().unwrap(),
        TOL_PURE,
        1.0,
        &format!("{case}: SLE"),
    );
    assert_eq!(af.sharp_te, d.header["SHARP"] == "T", "{case}: SHARP");
    println!(
        "{case}: PANGEN {} nodes — {bits}/{} values bitwise identical, worst |diff| {worst:.2e}",
        d.panels.len(),
        3 * d.panels.len()
    );
}

#[test]
fn test_naca0012_n160_matches_xfoil() {
    check("pangen_naca0012_n160", "0012", 160);
}

#[test]
fn test_naca0012_n81_matches_xfoil() {
    check("pangen_naca0012_n81", "0012", 81);
}

#[test]
fn test_naca4412_n160_matches_xfoil() {
    check("pangen_naca4412_n160", "4412", 160);
}

#[test]
fn test_naca4412_n81_matches_xfoil() {
    check("pangen_naca4412_n81", "4412", 81);
}

#[test]
fn test_naca23012_n160_matches_xfoil() {
    check("pangen_naca23012_n160", "23012", 160);
}

#[test]
fn test_naca4412_n160_refined_ppar_matches_xfoil() {
    // PPAR off default with both refinement windows: CVPAR 1.5, CTERAT 0.3, CTRRAT 0.5,
    // XT 0.2 0.4, XB 0.3 0.6 — the XSREF/XPREF branches of PANGEN
    check("pangen_naca4412_n160_refined", "4412", 160);
}

#[test]
fn test_naca0012_n120_bunching_ppar_matches_xfoil() {
    // CVPAR 0.5, CTERAT 0.05, no windows
    check("pangen_naca0012_n120_bunching", "0012", 120);
}
