//! The NACA series generators against the public-domain NASA/PDAS `naca456` reference
//! (`tests/fixtures/naca456/`, written by `scripts/naca456-fixtures.sh` from naca456's own
//! module procedures at 17 significant figures). Every case in `scripts/naca456-fixtures/cases.txt`
//! is gated: the thickness form and its slope, the mean line and its slope, and the combined
//! upper and lower surfaces at naca456's 98 "very fine" stations. A missing fixture fails.
//!
//! Tolerances are in `tests/utilities/tolerances.rs`; the 6-series thickness gate is derived
//! from the reference's own 1e-6 root tolerance (see there).

mod fixtures;
mod utilities;

use std::collections::HashMap;
use utilities::tolerances::{
    assert_within, within, NACA456_ROOT_TOL, TOL_NACA456_CLOSED_FORM, TOL_NACA456_SIX_SERIES_MEAN_LINE,
    TOL_NACA456_SIX_SERIES_THICKNESS,
};
use yfoil::geometry::{MeanLine, Section, ThicknessForm};

struct Fixture {
    header: HashMap<String, String>,
    rows: Vec<[f64; 9]>,
}

fn load(slug: &str) -> Fixture {
    let p = fixtures::require_fixture(&format!("tests/fixtures/naca456/{slug}.dat"));
    let text = std::fs::read_to_string(&p).unwrap();
    let mut header = HashMap::new();
    let mut rows = Vec::new();
    for l in text.lines() {
        if let Some(h) = l.strip_prefix("# ") {
            if let Some((k, v)) = h.split_once(' ') {
                header.insert(k.to_string(), v.trim().to_string());
            }
        } else if !l.trim().is_empty() {
            let v: Vec<f64> = l.split_whitespace().map(|t| t.parse().unwrap()).collect();
            assert_eq!(v.len(), 9, "{slug}: row width");
            rows.push([v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]]);
        }
    }
    assert_eq!(
        rows.len(),
        header["stations"].parse::<usize>().unwrap(),
        "{slug}: station count"
    );
    Fixture { header, rows }
}

/// The section the fixture describes, built from its `name` (and `a`), and checked against
/// the namelist parameters the driver echoed
fn section_of(f: &Fixture, slug: &str) -> Section {
    let name = &f.header["name"];
    let a: f64 = f.header["a"].parse().unwrap();
    let mut s = Section::from_designation(name).unwrap_or_else(|e| panic!("{slug}: {e}"));
    if let MeanLine::SixSeries { .. } = s.mean_line {
        if a != 1.0 {
            s = s.with_a(a).unwrap();
        }
    }
    let h = |k: &str| -> f64 { f.header[k].parse().unwrap() };
    assert_eq!(s.thickness_form.t(), h("toc"), "{slug}: t");
    match (&s.thickness_form, f.header["profile"].as_str()) {
        (ThicknessForm::FourDigit { .. }, "4") => {}
        (
            ThicknessForm::FourDigitModified {
                le_radius_index,
                x_max_thickness,
                ..
            },
            "4M",
        ) => {
            assert_eq!(*le_radius_index, h("leindex"), "{slug}: I");
            assert_eq!(*x_max_thickness, h("xmaxt"), "{slug}: M");
        }
        (ThicknessForm::SixSeries(form), p) => assert_eq!(form.family.label(), p, "{slug}: family"),
        (t, p) => panic!("{slug}: thickness form {t:?} vs profile {p}"),
    }
    match (&s.mean_line, f.header["camber"].as_str()) {
        // a symmetric section of any family: the driver was given no mean line
        (m, "0") if m.is_symmetric() => {}
        (MeanLine::TwoDigit { m, p }, "2") => {
            assert_eq!(*m, h("cmax"), "{slug}: m");
            assert_eq!(*p, h("xmaxc"), "{slug}: p");
        }
        (MeanLine::ThreeDigit { cl, p }, "3") | (MeanLine::ThreeDigitReflex { cl, p }, "3R") => {
            assert_eq!(*cl, h("cl"), "{slug}: cl");
            assert_eq!(*p, h("xmaxc"), "{slug}: p");
        }
        (MeanLine::SixSeries { a, cl }, "6") => {
            assert_eq!(*cl, h("cl"), "{slug}: cl");
            assert_eq!(*a, h("a"), "{slug}: a");
        }
        (MeanLine::SixSeriesModified { cl }, "6A") => assert_eq!(*cl, h("cl"), "{slug}: cl"),
        (m, c) => panic!("{slug}: mean line {m:?} vs camber {c}"),
    }
    s
}

/// Worst relative error over the case, per column, for the report in the assertion message
#[derive(Default)]
struct Worst {
    yt: f64,
    ytp: f64,
    ym: f64,
    ymp: f64,
    surface: f64,
}

fn check(slug: &str) -> Worst {
    let f = load(slug);
    let s = section_of(&f, slug);
    let six_series_thickness = matches!(s.thickness_form, ThicknessForm::SixSeries(_));
    let six_series_mean_line = matches!(
        s.mean_line,
        MeanLine::SixSeries { .. } | MeanLine::SixSeriesModified { .. }
    );
    let tol_t = if six_series_thickness {
        TOL_NACA456_SIX_SERIES_THICKNESS
    } else {
        TOL_NACA456_CLOSED_FORM
    };
    let tol_m = if six_series_mean_line {
        TOL_NACA456_SIX_SERIES_MEAN_LINE
    } else {
        TOL_NACA456_CLOSED_FORM
    };
    let aft_4m = |x: f64| -> bool {
        matches!(s.thickness_form, ThicknessForm::FourDigitModified { x_max_thickness, .. } if x >= x_max_thickness)
    };
    let mut worst = Worst::default();
    let rel = |a: f64, b: f64, scale: f64| (a - b).abs() / a.abs().max(b.abs()).max(scale);
    for r in &f.rows {
        let [x, yt, ytp, ym, ymp, xu, yu, xl, yl] = *r;
        let (yt_y, ytp_y) = s.thickness_form.at(x);
        let (ym_y, ymp_y) = s.mean_line.at(x);
        let at = |what: &str| format!("{slug} at x = {x}: {what}");

        // thickness: the 6-series ordinate carries the reference's root-finding error
        let tol_yt = if six_series_thickness {
            tol_t + NACA456_ROOT_TOL * ytp_y.abs()
        } else {
            tol_t
        };
        assert_within(yt_y, yt, tol_yt, 1.0, &at("y_t"));
        worst.yt = worst.yt.max(rel(yt_y, yt, 1.0));
        // slope: infinite at the leading edge (naca456 writes 1e20/1e22); naca456 reports the
        // aft slope of the 4-digit modified form as dy/d(1 − x), i.e. with the sign reversed
        if x > 0.0 && ytp.abs() < 1e10 {
            let ytp_ref = if aft_4m(x) { -ytp } else { ytp };
            if six_series_thickness {
                // the reference slope is y'(s)/x'(s) at the station reached, not requested; the
                // slope of the slope scales the error the same way as the ordinate
                let (_, ytp_h) = s.thickness_form.at(x + 1e-6);
                let (_, ytp_l) = s.thickness_form.at((x - 1e-6).max(0.0));
                let curvature = ((ytp_h - ytp_l) / 2e-6).abs();
                let tol = tol_t + NACA456_ROOT_TOL * curvature;
                assert!(
                    within(ytp_y, ytp_ref, tol, 1.0),
                    "{}: {ytp_y:.17e} vs {ytp_ref:.17e} (tol {tol:.2e})",
                    at("dy_t/dx")
                );
            } else {
                assert_within(ytp_y, ytp_ref, tol_t, 1.0, &at("dy_t/dx"));
            }
            worst.ytp = worst.ytp.max(rel(ytp_y, ytp_ref, 1.0));
        }
        assert_within(ym_y, ym, tol_m, 1.0, &at("y_c"));
        assert_within(ymp_y, ymp, tol_m, 1.0, &at("dy_c/dx"));
        worst.ym = worst.ym.max(rel(ym_y, ym, 1.0));
        worst.ymp = worst.ymp.max(rel(ymp_y, ymp, 1.0));

        // the combined surfaces, perpendicular application: the reference's thickness error
        // enters both coordinates through y_t
        let (xu_y, yu_y, xl_y, yl_y) = s.surface_at(x);
        let tol_s = tol_yt.max(tol_m);
        for (a, b, what) in [
            (xu_y, xu, "x_upper"),
            (yu_y, yu, "y_upper"),
            (xl_y, xl, "x_lower"),
            (yl_y, yl, "y_lower"),
        ] {
            assert_within(a, b, tol_s, 1.0, &at(what));
            worst.surface = worst.surface.max(rel(a, b, 1.0));
        }
    }
    worst
}

/// Every case of `scripts/naca456-fixtures/cases.txt`, as listed by its manifest
fn cases() -> Vec<String> {
    let p = fixtures::require_fixture("tests/fixtures/naca456/manifest.json");
    let m: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    m["cases"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn every_naca456_case_is_reproduced() {
    let cases = cases();
    assert!(cases.len() >= 30, "expected the full case list, got {}", cases.len());
    for slug in &cases {
        let w = check(slug);
        println!(
            "{slug:16} worst rel err: y_t {:.1e}  dy_t/dx {:.1e}  y_c {:.1e}  dy_c/dx {:.1e}  surface {:.1e}",
            w.yt, w.ytp, w.ym, w.ymp, w.surface
        );
    }
}

#[test]
fn closed_form_families_agree_to_round_off() {
    // the gate above already applies TOL_NACA456_CLOSED_FORM; this records that the achieved
    // agreement is round-off, an order of magnitude inside it
    for slug in [
        "naca0012",
        "naca4412",
        "naca0012-34",
        "naca23018",
        "naca23112",
        "naca16-009",
    ] {
        let w = check(slug);
        assert!(
            w.yt < 1e-13 && w.ym < 1e-13 && w.surface < 1e-13,
            "{slug}: {:.1e} {:.1e} {:.1e}",
            w.yt,
            w.ym,
            w.surface
        );
    }
}

#[test]
fn six_series_families_agree_within_the_reference_root_tolerance() {
    for slug in [
        "naca63-415",
        "naca64-010",
        "naca65-410",
        "naca66-021",
        "naca67-021",
        "naca64A010",
        "naca63A415",
        "naca65A210",
    ] {
        let w = check(slug);
        // the ordinates differ by the reference's own 1e-6 station error times the slope
        assert!(w.yt < 1e-4, "{slug}: y_t worst {:.1e}", w.yt);
    }
}

#[test]
fn six_series_mean_lines_agree_to_the_reference_pi() {
    // naca456's PI = 3.141592654 puts its 6-series lines 1.1e-10 (relative) from yFoil's
    for slug in [
        "naca16-212",
        "naca63-415",
        "naca63-415-a05",
        "naca63-415-a00",
        "naca63A415",
    ] {
        let w = check(slug);
        assert!(
            w.ym < 1e-9 && w.ymp < 1e-9,
            "{slug}: y_c {:.1e} dy_c/dx {:.1e}",
            w.ym,
            w.ymp
        );
    }
}
