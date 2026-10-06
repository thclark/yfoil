//! (a) Subroutine equivalence: SPECCL alone: the inviscid alpha for a specified CL, which VISCAL then starts from, against
//! the reference's `viscal_inviscid.dat` (ALFA, CL entering the loop) on `naca0012_n60_cl03_re1e6`.
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_cl03_re1e6/` — `cargo xtask fixtures --case naca0012_n60_cl03_re1e6`.

use crate::fixtures;
use crate::utilities::tolerances::{assert_within, TOL_SOLVER};
use std::collections::HashMap;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};
use yfoil::solver::specal::cl_command;

#[test]
fn test_speccl_alpha_for_cl_matches_xfoil() {
    let dir = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_cl03_re1e6");
    let g = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let mut session = Session::new(&panel_foil(&g), FlowConditions::default());
    let (st, sys) = session.parts_mut();
    cl_command(st, sys, 0.3);
    let mut inv = HashMap::new();
    for l in std::fs::read_to_string(dir.join("viscal_inviscid.dat"))
        .unwrap()
        .lines()
    {
        if let Some((k, v)) = l.split_once('=') {
            inv.entry(k.trim().to_string()).or_insert(v.trim().to_string());
        }
    }
    assert_within(
        session.state().alpha,
        inv["ALFA"].parse().unwrap(),
        TOL_SOLVER,
        1.0,
        "SPECCL alpha",
    );
    assert_within(
        session.state().cl,
        inv["CL"].parse().unwrap(),
        TOL_SOLVER,
        1.0,
        "SPECCL CL",
    );
}
