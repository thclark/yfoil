//! S2 gate: the BL pointer layer (STFIND, IBLPAN, XICALC, IBLSYS) and the inviscid velocity
//! layer (QISET, UICALC) against `xfoil_pointers.dat` / `xfoil_uinv.dat` from the tracked
//! reference case. Integer pointers must be exactly equal; pure-arithmetic reals bit-identical;
//! quantities that pass through libm (`cos`, `sin`) or YFoil's own spline derivatives are
//! held to `TOL_PURE`.

mod fixtures;
mod utilities;

use fixtures::pointers_fixtures::{parse_inviscid_gam, parse_pointers, parse_uinv, PointersFixture};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::pointers::{
    find_stagnation, map_stations_to_nodes, map_stations_to_rows, set_station_xi, set_te_thickness,
};
use yfoil::solver::velocity::{set_q_inviscid, set_ue_inviscid};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn pointers() -> PointersFixture {
    parse_pointers(&fixture_path("xfoil_pointers.dat"), 1)
}

/// A state whose panel arrays are XFOIL's own (nodes incl. wake, GAM at iteration 0), so each
/// subroutine is tested on exactly the inputs XFOIL gave it.
fn state_from_fixture(f: &PointersFixture) -> SolverState {
    let mut st = SolverState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    st.te_thickness_normal = f.ante;
    st.sharp_te = f.sharp;
    let gam = parse_inviscid_gam(&fixture_path("xfoil_inviscid.dat"));
    assert_eq!(gam.len(), f.n + 1, "GAM block length");
    st.gamma = gam;
    // xp/yp only enter XICALC's wake-gap cubic; take YFoil's spline for those
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    for i in 1..=f.n {
        st.dxds[i] = af.xp[i - 1];
        st.dyds[i] = af.yp[i - 1];
    }
    st
}

fn assert_bits(a: f64, b: f64, what: &str) {
    assert_eq!(a.to_bits(), b.to_bits(), "{what}: yfoil={a:.17e} xfoil={b:.17e}");
}

#[test]
fn test_tecalc_matches_xfoil() {
    let f = pointers();
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    let st = SolverState::from_foil(&af, f.nw);
    // arc length is pure arithmetic (SCALC): bit-identical
    for i in 1..=f.n {
        assert_bits(st.s[i], f.s[i], &format!("S({i})"));
    }
    assert_eq!(st.sharp_te, f.sharp, "SHARP");
    // ANTE depends on spline derivatives XP/YP at the TE
    assert_within(st.te_thickness_normal, f.ante, TOL_PURE, 1.0, "ANTE");
    let mut st2 = st.clone();
    set_te_thickness(&mut st2);
    assert_bits(st2.te_thickness_normal, st.te_thickness_normal, "TECALC idempotent");
}

#[test]
fn test_stfind_matches_xfoil() {
    let f = pointers();
    let mut st = state_from_fixture(&f);
    find_stagnation(&mut st);
    assert_eq!(st.i_stagnation_node, f.ist, "IST");
    assert_bits(st.s_stagnation, f.sst, "SST");
    assert_bits(st.s_stagnation_d_gamma_node0, f.sst_go, "SST_GO");
    assert_bits(st.s_stagnation_d_gamma_node1, f.sst_gp, "SST_GP");
}

#[test]
fn test_iblpan_and_iblsys_match_xfoil() {
    let f = pointers();
    let mut st = state_from_fixture(&f);
    st.i_stagnation_node = f.ist;
    map_stations_to_nodes(&mut st);
    map_stations_to_rows(&mut st);
    assert_eq!(st.i_te_station[1..], f.iblte[1..], "IBLTE");
    assert_eq!(st.n_stations[1..], f.nbl[1..], "NBL");
    assert_eq!(st.n_rows, f.nsys, "NSYS");
    for is in 1..=2 {
        for ibl in 1..=(f.iblte[is] + f.nw) {
            assert_eq!(st.i_node[is][ibl], f.ipan[is][ibl], "IPAN({ibl},{is})");
            assert_eq!(
                st.velocity_sign[is][ibl].to_bits(),
                f.vti[is][ibl].to_bits(),
                "VTI({ibl},{is})"
            );
            assert_eq!(st.i_row[is][ibl], f.isys[is][ibl], "ISYS({ibl},{is})");
        }
    }
    // reference-case sanity (from IBLPAN's definition): NBL(2) = IBLTE(2) + NW, NSYS = ΣNBL-2
    assert_eq!(f.nbl[2], f.iblte[2] + f.nw);
    assert_eq!(f.nsys, f.nbl[1] + f.nbl[2] - 2);
}

#[test]
fn test_xicalc_matches_xfoil() {
    let f = pointers();
    let mut st = state_from_fixture(&f);
    st.i_stagnation_node = f.ist;
    st.s_stagnation = f.sst;
    map_stations_to_nodes(&mut st);
    set_station_xi(&mut st);
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            assert_bits(st.xi[is][ibl], f.xssi[is][ibl], &format!("XSSI({ibl},{is})"));
        }
        // upper-side wake mirror (plotting pointers) is filled by XICALC too
        if is == 1 {
            for iw in 1..=f.nw {
                let ibl = f.iblte[1] + iw;
                assert_bits(st.xi[1][ibl], f.xssi[1][ibl], &format!("XSSI({ibl},1) wake mirror"));
            }
        }
    }
    // WGAP goes through YFoil's XP/YP (spline derivatives): tolerance, scale = ANTE
    for iw in 1..=f.nw {
        assert_within(st.wake_gap[iw], f.wgap[iw], TOL_PURE, f.ante, &format!("WGAP({iw})"));
    }
}

#[test]
fn test_qiset_and_uicalc_match_xfoil() {
    let f = pointers();
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    assert_eq!((u.n, u.nw), (f.n, f.nw));
    let mut st = state_from_fixture(&f);
    st.i_stagnation_node = f.ist;
    map_stations_to_nodes(&mut st);
    st.q_inviscid_basis[1] = u.qinvu1.clone();
    st.q_inviscid_basis[2] = u.qinvu2.clone();
    set_q_inviscid(&mut st, u.alfa);
    // cos/sin come from libm: same host gives identical bits, other hosts an ULP — TOL_PURE
    for i in 1..=(f.n + f.nw) {
        assert_within(st.q_inviscid[i], u.qinv[i], TOL_PURE, 1.0, &format!("QINV({i})"));
        assert_within(
            st.q_inviscid_d_alpha[i],
            u.qinv_a[i],
            TOL_PURE,
            1.0,
            &format!("QINV_A({i})"),
        );
    }
    // UICALC is a sign flip: feed XFOIL's QINV and require bit-identity
    st.q_inviscid = u.qinv.clone();
    st.q_inviscid_d_alpha = u.qinv_a.clone();
    set_ue_inviscid(&mut st);
    for is in 1..=2 {
        assert_eq!(u.uinv[is].len(), f.nbl[is] + 1, "UINV rows side {is}");
        for ibl in 1..=f.nbl[is] {
            assert_bits(st.ue_inviscid[is][ibl], u.uinv[is][ibl], &format!("UINV({ibl},{is})"));
            assert_bits(
                st.ue_inviscid_d_alpha[is][ibl],
                u.uinv_a[is][ibl],
                &format!("UINV_A({ibl},{is})"),
            );
        }
    }
}

/// The whole VISCAL prologue from YFoil's own geometry: TECALC → (wake nodes from XFOIL, until
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

/// NCALC / APCALC: YFoil's node normals and panel angles for the airfoil against XFOIL's.
/// Normals come from spline derivatives, angles from atan2 — both held to TOL_PURE.
#[test]
fn test_airfoil_normals_and_panel_angles_match_xfoil() {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    for i in 1..=f.n {
        assert_within(af.nx[i - 1], f.nx[i], TOL_PURE, 1.0, &format!("NX({i})"));
        assert_within(af.ny[i - 1], f.ny[i], TOL_PURE, 1.0, &format!("NY({i})"));
        assert_within(af.apanel[i - 1], f.apanel[i], TOL_PURE, 1.0, &format!("APANEL({i})"));
    }
}
