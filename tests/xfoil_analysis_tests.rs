//! S9 end-to-end gate: from YFoil's own panels (the fixture pipeline's `panels.json`) through
//! the inviscid solve, wake, pointers, DIJ and the VISCAL loop to the converged point, against
//! the instrumented reference's `viscal_final.dat` / `viscal_iter.dat` on the same panels.
//!
//! Everything upstream of the Newton loop is YFoil's here (GGCALC, XYWAKE, QDCALC, MRCHUE), so
//! this is the propagated-floor comparison: the S3/S4 gates put the wake position and DIJ at
//! ~5e-11, which is what bounds the final forces (TOL_SOLVER is at the measured noise floor).

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::parse_bl_dump;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

#[test]
fn test_end_to_end_alpha_2_matches_xfoil_converged_point() {
    let geometry = read_geometry_from_file(fixture_path("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = create_paneled_airfoil(&geometry);
    let spec = FlowConditions {
        re: 1.0e6,
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    };
    let mut session = Session::new(&airfoil, spec);
    let p = session.alpha(2.0_f64.to_radians());

    let fin = parse_bl_dump(&fixture_path("viscal_final.dat"));
    let iters = std::fs::read_to_string(fixture_path("viscal_iter.dat"))
        .unwrap()
        .lines()
        .filter(|l| l.starts_with("ITER="))
        .count();

    // identical iteration count and convergence
    assert_eq!(p.iterations, fin.int("NITDONE"), "VISCAL iteration count");
    assert_eq!(p.iterations, iters);
    assert_eq!(p.converged, fin.logical("LVCONV"), "LVCONV");

    // the converged point
    for (name, ours) in [
        ("CL", p.cl),
        ("CM", p.cm),
        ("CD", p.cd),
        ("CDF", p.cd_friction),
        ("CDP", p.cd_pressure),
        ("CL_ALF", p.cl_d_alpha),
        ("XOCTR1", p.transition_upper[0]),
        ("XOCTR2", p.transition_lower[0]),
    ] {
        assert_within(ours, fin.real(name), TOL_SOLVER, 1.0, &format!("end-to-end: {name}"));
    }
    assert_eq!(session.state.i_stagnation_node, fin.int("IST"), "IST");
    assert_eq!(
        p.i_transition_station[1..],
        [fin.int("ITRAN1"), fin.int("ITRAN2")],
        "ITRAN"
    );

    // per-iteration RMSBL/RLX/CL/CD
    let text = std::fs::read_to_string(fixture_path("viscal_iter.dat")).unwrap();
    let mut blocks: Vec<std::collections::HashMap<String, String>> = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "ITER" {
            blocks.push(Default::default());
        }
        if let Some(b) = blocks.last_mut() {
            b.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    for (y, x) in p.iteration_records.iter().zip(&blocks) {
        for (name, ours) in [
            ("RMSBL", y.residual),
            ("RLX", y.relaxation),
            ("CL", y.cl),
            ("CD", y.cd),
            ("CM", y.cm),
        ] {
            assert_within(
                ours,
                x[name].parse().unwrap(),
                TOL_SOLVER,
                1.0,
                &format!("end-to-end iteration {}: {name}", y.iteration),
            );
        }
        assert_eq!(
            y.i_stagnation_node,
            x["IST"].parse::<usize>().unwrap(),
            "iteration {}: IST",
            y.iteration
        );
    }

    // node arrays
    let mut nodes = 0;
    for l in std::fs::read_to_string(fixture_path("viscal_final.dat"))
        .unwrap()
        .lines()
    {
        let Some(r) = l.strip_prefix("NODE(") else { continue };
        let (idx, vals) = r.split_once(")=").unwrap();
        let i: usize = idx.trim().parse().unwrap();
        let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
        assert_within(
            session.state.cp_inviscid[i],
            v[0],
            TOL_SOLVER,
            1.0,
            &format!("end-to-end: CPI({i})"),
        );
        assert_within(
            session.state.cp_viscous[i],
            v[1],
            TOL_SOLVER,
            1.0,
            &format!("end-to-end: CPV({i})"),
        );
        nodes += 1;
    }
    assert_eq!(nodes, session.state.n_foil_nodes + session.state.n_wake_nodes);

    println!(
        "end-to-end alpha=2: {} iterations, CL {:.10} CD {:.10} CM {:.10} XTR {:.6}/{:.6} — all within {TOL_SOLVER:.0e} of XFOIL",
        p.iterations, p.cl, p.cd, p.cm, p.transition_upper[0], p.transition_lower[0]
    );
}
