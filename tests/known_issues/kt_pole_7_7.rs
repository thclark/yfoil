//! (d) Known XFOIL weakness: known issues §7.7 and §4 (CPCALC/CLCALC): TYPE 2 drives itself onto the Kármán–Tsien pole,
//! and XFOIL reports numbers past it without a warning. yFoil's defined behaviour is to replicate
//! the numbers and classify them as invalid (`docs/guide/validity.md`); these tests check that
//! classification against XFOIL's own record of the same occasions.
//!
//! The instrumented reference logs every occasion the validity record describes as an event in
//! `events.dat` (`xfoil/instrumentation/instrument/19-validity-events.patch`). Each SPECAL call opens
//! a span with a `SPECAL_ENTER` marker, so the events of one operating point are the lines between
//! its marker and the next. An event fires on *any* call within the span — including the
//! intermediate iterates of SPECAL's CL(M) Newton — while yFoil's record describes the *reported*
//! state, the last call. So where XFOIL's record is per point (the CL(M) Newton exhausted, logged
//! once after the loop) yFoil's flag must equal the event; where it is per call (MRCL's
//! substitutions, the Kármán–Tsien tests) a yFoil flag must have its event in the span.
//!
//! Fixtures (`cargo xtask fixtures --case <name>`): `naca4412_n60_inviscid_m07_a10`,
//! `naca64a010_n60_inviscid_sweep30_type2_m03`.

use crate::fixtures;
use std::path::Path;
use yfoil::bl::system::{MachClDependence, ReClDependence};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};
use yfoil::solver::validity::ValidityRecord;

fn case_dir(case: &str) -> std::path::PathBuf {
    fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}"))
}

/// The events of one SPECAL call: the sites logged between its `SPECAL_ENTER` and the next.
#[derive(Debug, Default)]
struct Span {
    alpha_deg: f64,
    sites: Vec<String>,
}

impl Span {
    fn has(&self, site: &str) -> bool {
        self.sites.iter().any(|s| s == site)
    }
}

fn spans(dir: &Path) -> Vec<Span> {
    let text = std::fs::read_to_string(dir.join("events.dat")).expect("events.dat");
    let mut out: Vec<Span> = vec![];
    for line in text.lines() {
        let t: Vec<&str> = line.split_whitespace().collect();
        assert!(t.len() == 6 && t[0] == "EV", "malformed event line: {line}");
        if t[1] == "SPECAL_ENTER" {
            let alfa: f64 = t[5].parse().unwrap();
            out.push(Span {
                alpha_deg: alfa.to_degrees(),
                sites: vec![],
            });
        } else {
            out.last_mut()
                .expect("an event before the first SPECAL_ENTER")
                .sites
                .push(t[1].to_string());
        }
    }
    out
}

/// Run the case's ALFA sequence inviscid through one session, as the reference did, and pair
/// each point's validity record with XFOIL's span of events.
fn run_inviscid(case: &str, mach: f64, type2: bool) -> Vec<(Span, ValidityRecord)> {
    let dir = case_dir(case);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let mut session = Session::new(
        &panel_foil(&geometry),
        FlowConditions {
            re: None,
            mach,
            mach_cl_dependence: if type2 {
                MachClDependence::InverseSqrtCl
            } else {
                MachClDependence::Fixed
            },
            re_cl_dependence: if type2 {
                ReClDependence::InverseSqrtCl
            } else {
                ReClDependence::Fixed
            },
            ..FlowConditions::default()
        },
    );
    spans(&dir)
        .into_iter()
        .map(|span| {
            let p = session.alpha(span.alpha_deg.to_radians());
            (span, p.validity)
        })
        .collect()
}

/// The per-call flags: a flag on the reported state needs its event in XFOIL's span.
fn assert_flags_have_events(points: &[(Span, ValidityRecord)]) {
    for (span, v) in points {
        let a = span.alpha_deg;
        let pairs = [
            (v.karman_tsien_margin_forces <= 0.0, "CLCALC_KT_DOMAIN"),
            (v.karman_tsien_margin_pressure <= 0.0, "CPCALC_KT_DOMAIN"),
            (v.cl_floored, "MRCL_CL_FLOOR"),
            (v.mach_limited, "MRCL_MACH_LIMIT"),
            (v.re_limited, "MRCL_RE_LIMIT"),
        ];
        for (flag, site) in pairs {
            assert!(
                !flag || span.has(site),
                "alpha {a}: yFoil's reported state took the {site} branch, XFOIL's span never did"
            );
        }
    }
}

#[test]
fn test_point_past_the_karman_tsien_pole_matches_xfoil() {
    let points = run_inviscid("naca4412_n60_inviscid_m07_a10", 0.7, false);
    assert_eq!(points.len(), 1);
    let (span, v) = &points[0];
    // at a fixed Mach every call of the span sees the same denominator, so here the event and
    // the flag must agree in both directions
    assert!(
        span.has("CLCALC_KT_DOMAIN") && span.has("CPCALC_KT_DOMAIN"),
        "the reference crossed the pole"
    );
    assert!(
        v.karman_tsien_margin_forces <= 0.0 && v.karman_tsien_margin_pressure <= 0.0,
        "yFoil must record the crossing XFOIL logged: {v:?}"
    );
    assert!(v.out_of_domain() && !v.iteration_exhausted() && !v.conditions_substituted());
}

#[test]
fn test_type2_sweep_mach_newton_exhaustion_matches_xfoil_exactly() {
    let points = run_inviscid("naca64a010_n60_inviscid_sweep30_type2_m03", 0.3, true);
    assert_eq!(points.len(), 61, "one SPECAL per degree from -30 to 30");
    let ours: Vec<i64> = points
        .iter()
        .filter(|(_, v)| v.mach_cl_newton_exhausted)
        .map(|(s, _)| s.alpha_deg.round() as i64)
        .collect();
    let xfoil: Vec<i64> = points
        .iter()
        .filter(|(s, _)| s.has("SPECAL_MINF_FAIL"))
        .map(|(s, _)| s.alpha_deg.round() as i64)
        .collect();
    assert!(!xfoil.is_empty(), "the sweep must exercise the exhaustion");
    assert_eq!(ours, xfoil, "SPECAL's CL(M) Newton exhausted at exactly XFOIL's points");
}

#[test]
fn test_type2_sweep_flags_are_taken_only_where_xfoil_took_them() {
    let points = run_inviscid("naca64a010_n60_inviscid_sweep30_type2_m03", 0.3, true);
    assert_flags_have_events(&points);
    // and the sweep is not vacuous: every kind of flag is raised somewhere
    assert!(points.iter().any(|(_, v)| v.out_of_domain()));
    assert!(points.iter().any(|(_, v)| v.cl_floored && v.mach_limited));
}
