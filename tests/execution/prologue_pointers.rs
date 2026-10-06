//! (b) Execution equivalence: the VISCAL prologue from yFoil's own geometry — TECALC → QISET → STFIND → IBLPAN → XICALC → IBLSYS → UICALC — against `xfoil_pointers.dat` / `xfoil_uinv.dat`.
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_a2_re1e6/` — `cargo xtask fixtures --case naca0012_n60_a2_re1e6`.

use crate::fixtures;
use crate::utilities;

use fixtures::pointers_fixtures::{parse_pointers, parse_uinv, PointersFixture};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::pointers::{find_stagnation, map_stations_to_nodes, map_stations_to_rows, set_station_xi};
use yfoil::solver::velocity::{set_q_inviscid, set_ue_inviscid};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn pointers() -> PointersFixture {
    parse_pointers(&fixture_path("xfoil_pointers.dat"), 1)
}

/// The whole VISCAL prologue from yFoil's own geometry: TECALC → (wake nodes from XFOIL, until
/// XYWAKE lands in S3) → QISET → GAM=QINV → STFIND → IBLPAN → XICALC → IBLSYS → UICALC.
#[test]
fn test_prologue_from_yfoil_geometry_matches_xfoil() {
    let f = pointers();
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    let mut st = SolverState::from_foil(&af, f.nw);
    st.set_wake_nodes(&f.x[f.n + 1..], &f.y[f.n + 1..], &f.s[f.n + 1..]);
    st.q_inviscid_basis[1] = u.qinvu1.clone();
    st.q_inviscid_basis[2] = u.qinvu2.clone();
    set_q_inviscid(&mut st, u.alfa);
    // before the first viscous iteration GAM is the inviscid solution at this alpha
    for i in 1..=f.n {
        st.gamma[i] = st.q_inviscid[i];
    }
    find_stagnation(&mut st);
    map_stations_to_nodes(&mut st);
    set_station_xi(&mut st);
    map_stations_to_rows(&mut st);
    set_ue_inviscid(&mut st);

    assert_eq!(st.i_stagnation_node, f.ist, "IST");
    assert_within(st.s_stagnation, f.sst, TOL_PURE, 1.0, "SST");
    assert_eq!((st.n_stations, st.i_te_station, st.n_rows), (f.nbl, f.iblte, f.nsys));
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            assert_eq!(st.i_node[is][ibl], f.ipan[is][ibl]);
            assert_eq!(st.i_row[is][ibl], f.isys[is][ibl]);
            assert_within(
                st.xi[is][ibl],
                f.xssi[is][ibl],
                TOL_PURE,
                1.0,
                &format!("XSSI({ibl},{is})"),
            );
            assert_within(
                st.ue_inviscid[is][ibl],
                u.uinv[is][ibl],
                TOL_PURE,
                1.0,
                &format!("UINV({ibl},{is})"),
            );
        }
    }
    println!(
        "prologue: IST={} NBL={:?} NSYS={} — pointers exact, XSSI/UINV within {TOL_PURE:.0e}",
        st.i_stagnation_node,
        &st.n_stations[1..],
        st.n_rows
    );
}
