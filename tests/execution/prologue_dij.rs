//! (b) Execution equivalence: the VISCAL prologue from yFoil's own geometry — GGCALC → QISET → SPECAL's GAM → XYWAKE → QWCALC → QDCALC — against the reference run's QINVU and DIJ.
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_a2_re1e6/` — `cargo xtask fixtures --case naca0012_n60_a2_re1e6`.

use crate::fixtures;
use crate::utilities;

use fixtures::pointers_fixtures::{parse_dij, parse_pointers, parse_uinv};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_PURE, TOL_SOLVER};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::ggcalc::{build_inviscid_system, InviscidSystem};
use yfoil::solver::qdcalc::build_dij;
use yfoil::solver::velocity::set_q_inviscid;
use yfoil::solver::xywake::{build_wake, set_wake_q_basis};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// GGCALC on yFoil's geometry, then the alpha superposition of SPECAL.
fn prologue_to_specal() -> (SolverState, InviscidSystem, fixtures::pointers_fixtures::UinvFixture) {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    let mut st = SolverState::from_foil(&af, f.nw);
    let sys = build_inviscid_system(&mut st);
    st.alpha = u.alfa;
    set_q_inviscid(&mut st, u.alfa);
    // SPECAL: GAM(I) = COSA*GAMU(I,1) + SINA*GAMU(I,2)  (= QINV on the airfoil)
    for i in 1..=st.n_foil_nodes {
        st.gamma[i] = st.q_inviscid[i];
    }
    (st, sys, u)
}

/// Compare a DIJ against the reference with the row-scaled metric; returns the worst scaled error.
fn check_dij(st: &SolverState, tol: f64, what: &str) -> (f64, usize, usize) {
    let (n, nw, xd) = parse_dij(&fixture_path("xfoil_dij.dat"));
    assert_eq!((n, nw), (st.n_foil_nodes, st.n_wake_nodes));
    let np = n + nw;
    let mut worst = (0.0_f64, 0usize, 0usize);
    for i in 1..=np {
        let rowmax = (1..=np).map(|j| xd[i][j].abs()).fold(0.0_f64, f64::max);
        for j in 1..=np {
            let (a, b) = (st.dij[i][j], xd[i][j]);
            let scaled = (a - b).abs() / a.abs().max(b.abs()).max(rowmax);
            if scaled > worst.0 {
                worst = (scaled, i, j);
            }
            assert_within(
                a,
                b,
                tol,
                rowmax,
                &format!("{what}: DIJ({i},{j}) [rowmax {rowmax:.3e}]"),
            );
        }
    }
    // structural: first wake row equals the TE row, bitwise
    for j in 1..=np {
        assert_eq!(
            st.dij[n + 1][j].to_bits(),
            st.dij[n][j].to_bits(),
            "DIJ(N+1,{j}) != DIJ(N,{j})"
        );
    }
    println!(
        "{what}: {np}x{np} DIJ within {tol:.0e} row-scaled (worst {:.2e} at ({},{}))",
        worst.0, worst.1, worst.2
    );
    worst
}

/// The whole prologue from yFoil's geometry. DIJ's sensitivity to a wake node position is
/// O(1/distance) ≈ 3e2 here, and XYWAKE places the nodes within ~1e-12 of XFOIL's (its input
/// GAM is 4e-14 off after the LU solve), so the DIJ floor for this chained test is ~3e-10:
/// gated at TOL_SOLVER, not TOL_LINALG. The isolated test above is the translation gate.
#[test]
fn test_prologue_dij_from_yfoil_geometry() {
    let (mut st, mut sys, u) = prologue_to_specal();
    build_wake(&mut st, 1.0);
    set_wake_q_basis(&mut st);
    for i in (st.n_foil_nodes + 1)..=(st.n_foil_nodes + st.n_wake_nodes) {
        assert_within(
            st.q_inviscid_basis[1][i],
            u.qinvu1[i],
            TOL_PURE,
            1.0,
            &format!("wake QINVU({i},1)"),
        );
    }
    build_dij(&mut st, &mut sys);
    check_dij(&st, TOL_SOLVER, "qdcalc (yFoil wake)");
}
