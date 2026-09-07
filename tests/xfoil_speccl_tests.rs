//! S11 gates: fixed-CL mode (SPECCL + VISCAL with alpha as the unknown) and MATYP/RETYP ≠ 1
//! (MRCL: Re and Mach scale with CL), each against a minimal reference case that records the
//! VISCAL-level state: the inviscid arrays entering the loop (`viscal_inviscid.dat`, i.e.
//! SPECCL's converged alpha), every iteration (`viscal_iters_all.dat`, with ALFA/MINF/REINF)
//! and the converged point (`viscal_points.dat`).

mod fixtures;
mod utilities;

use std::collections::HashMap;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file, PaneledAirfoil};
use yfoil::solver::analysis::{FlowConditions, PointResult, Session};
use yfoil::solver::specal::cl_command;

fn case_path(case: &str, name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}/{name}"))
}

fn header(path: &PathBuf) -> HashMap<String, String> {
    let mut h = HashMap::new();
    for l in std::fs::read_to_string(path).unwrap().lines() {
        if let Some((k, v)) = l.split_once('=') {
            h.entry(k.trim().to_string()).or_insert(v.trim().to_string());
        }
    }
    h
}

/// `IT k i RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2 ALFA MINF REINF` lines of call 1.
fn iters_call_1(path: &PathBuf) -> Vec<Vec<f64>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| l.strip_prefix("IT "))
        .map(|r| {
            r.split_whitespace()
                .map(|t| t.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .filter(|v| v[0] as usize == 1)
        .collect()
}

fn airfoil(case: &str) -> PaneledAirfoil {
    let g = read_geometry_from_file(case_path(case, "panels.json").to_str().unwrap()).expect("panels.json");
    create_paneled_airfoil(&g)
}

/// Compare one converged point and its iteration history with the reference records.
fn check(case: &str, p: &PointResult, st: &yfoil::solver::blstate::SolverState) {
    let pt = header(&case_path(case, "viscal_points.dat"));
    let its = iters_call_1(&case_path(case, "viscal_iters_all.dat"));

    assert_eq!(
        p.iterations,
        pt["NITDONE"].parse::<usize>().unwrap(),
        "{case}: iteration count"
    );
    assert_eq!(p.converged, pt["LVCONV"] == "T", "{case}: LVCONV");
    assert_eq!(p.iterations, its.len(), "{case}: per-iteration record length");
    for (y, x) in p.iteration_records.iter().zip(&its) {
        let ctx = format!("{case} iteration {}", y.iteration);
        for (name, ours, theirs, scale) in [
            ("RMSBL", y.residual, x[2], 1.0),
            ("RLX", y.relaxation, x[3], 1.0),
            ("CL", y.cl, x[4], 1.0),
            ("CD", y.cd, x[5], 1.0),
            ("CM", y.cm, x[6], 1.0),
            ("ALFA", y.alpha, x[11], 1.0),
            ("MINF", y.mach, x[12], 1.0),
            ("REINF", y.re, x[13], x[13]),
        ] {
            assert_within(ours, theirs, TOL_SOLVER, scale, &format!("{ctx}: {name}"));
        }
        assert_eq!(y.i_stagnation_node, x[8] as usize, "{ctx}: IST");
        assert_eq!(
            y.i_transition_station[1..],
            [x[9] as usize, x[10] as usize],
            "{ctx}: ITRAN"
        );
        println!(
            "  {ctx}: rms {:.4e} rlx {:.4} alpha {:.8}° CL {:.8} CD {:.8} Re {:.1} — match",
            y.residual,
            y.relaxation,
            y.alpha.to_degrees(),
            y.cl,
            y.cd,
            y.re
        );
    }
    for (name, ours, scale) in [
        ("ALFA", p.alpha, 1.0),
        ("CL", p.cl, 1.0),
        ("CM", p.cm, 1.0),
        ("CD", p.cd, 1.0),
        ("CDF", p.cd_friction, 1.0),
        ("CDP", p.cd_pressure, 1.0),
        ("XOCTR1", p.transition_upper[0], 1.0),
        ("XOCTR2", p.transition_lower[0], 1.0),
        ("MINF", st.mach, 1.0),
        ("REINF", st.re, st.re),
    ] {
        assert_within(
            ours,
            pt[name].parse().unwrap(),
            TOL_SOLVER,
            scale,
            &format!("{case}: {name}"),
        );
    }
    assert_eq!(st.i_stagnation_node, pt["IST"].parse::<usize>().unwrap(), "{case}: IST");
    assert_eq!(
        p.i_transition_station[1..],
        [
            pt["ITRAN1"].parse::<usize>().unwrap(),
            pt["ITRAN2"].parse::<usize>().unwrap()
        ],
        "{case}: ITRAN"
    );
    println!(
        "{case}: {} iterations, alpha {:.6}° CL {:.8} CD {:.8} CM {:+.8} Re {:.1} XTR {:.5}/{:.5} — all within {TOL_SOLVER:.0e}",
        p.iterations,
        p.alpha.to_degrees(),
        p.cl,
        p.cd,
        p.cm,
        st.re,
        p.transition_upper[0],
        p.transition_lower[0]
    );
}

#[test]
fn test_fixed_cl_point_matches_xfoil() {
    let case = "naca0012_n60_cl03_re1e6";
    let af = airfoil(case);
    let mut session = Session::new(&af, FlowConditions::default());

    // SPECCL alone: the inviscid alpha for CL = 0.3 that VISCAL starts from
    cl_command(&mut session.state, &mut session.inviscid, 0.3);
    let inv = header(&case_path(case, "viscal_inviscid.dat"));
    assert_within(
        session.state.alpha,
        inv["ALFA"].parse().unwrap(),
        TOL_SOLVER,
        1.0,
        "SPECCL alpha",
    );
    assert_within(
        session.state.cl,
        inv["CL"].parse().unwrap(),
        TOL_SOLVER,
        1.0,
        "SPECCL CL",
    );
    assert!(!session.state.alpha_specified);

    // the full CL command: SPECCL again (identical), then VISCAL with alpha as the unknown
    let p = session.cl(0.3);
    check(case, &p, &session.state);
    assert_within(p.cl, 0.3, TOL_SOLVER, 1.0, "viscous CL meets CLSPEC");
}

#[test]
fn test_matyp_retyp_2_point_matches_xfoil() {
    let case = "naca0012_n60_a2_re1e6_type2";
    let af = airfoil(case);
    let spec = FlowConditions {
        mach_cl_dependence: 2,
        re_cl_dependence: 2,
        ..FlowConditions::default()
    };
    let mut session = Session::new(&af, spec);
    let p = session.alpha(2.0_f64.to_radians());
    // Re = Re1 / sqrt(CL) through MRCL, updated every iteration from the current CL
    assert_within(
        session.state.re,
        1.0e6 / p.cl.sqrt(),
        TOL_SOLVER,
        session.state.re,
        "REINF = REINF1/sqrt(CL)",
    );
    check(case, &p, &session.state);
}
