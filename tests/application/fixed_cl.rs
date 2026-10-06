//! (c) yFoil functionality: the operating-point commands' own contracts: a `CL` point leaves alpha unspecified and meets
//! the requested lift; under TYPE 2 the Reynolds number follows the lift as Re₁/√CL.
//!
//! Fixtures: the panels of `tests/fixtures/xfoil/naca0012_n60_cl03_re1e6/` — `cargo xtask fixtures --case naca0012_n60_cl03_re1e6`.

use crate::fixtures;
use crate::utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::bl::system::{MachClDependence, ReClDependence};
use yfoil::geometry::{panel_foil, read_geometry_from_file, PanelledFoil};
use yfoil::solver::analysis::{FlowConditions, Session};

fn airfoil() -> PanelledFoil {
    let p = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_cl03_re1e6/panels.json");
    panel_foil(&read_geometry_from_file(p.to_str().unwrap()).expect("panels.json"))
}

#[test]
fn test_cl_point_meets_the_requested_lift() {
    let mut session = Session::new(&airfoil(), FlowConditions::default());
    let p = session.cl(0.3);
    assert!(!session.state().alpha_specified, "a CL point solves for alpha");
    assert_within(p.cl, 0.3, TOL_SOLVER, 1.0, "viscous CL meets CLSPEC");
}

#[test]
fn test_type2_reynolds_number_follows_the_lift() {
    let spec = FlowConditions {
        mach_cl_dependence: MachClDependence::InverseSqrtCl,
        re_cl_dependence: ReClDependence::InverseSqrtCl,
        ..FlowConditions::default()
    };
    let mut session = Session::new(&airfoil(), spec);
    let p = session.alpha(2.0_f64.to_radians());
    assert_within(
        session.state().re,
        1.0e6 / p.cl.sqrt(),
        TOL_SOLVER,
        session.state().re,
        "REINF = REINF1/sqrt(CL)",
    );
}
