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
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::solver::blstate::BlState;
use yfoil::solver::pointers::{iblpan, iblsys, stfind, tecalc, xicalc};
use yfoil::solver::velocity::{qiset, uicalc};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn pointers() -> PointersFixture {
    parse_pointers(&fixture_path("xfoil_pointers.dat"), 1)
}

/// A state whose panel arrays are XFOIL's own (nodes incl. wake, GAM at iteration 0), so each
/// subroutine is tested on exactly the inputs XFOIL gave it.
fn state_from_fixture(f: &PointersFixture) -> BlState {
    let mut st = BlState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    st.ante = f.ante;
    st.sharp = f.sharp;
    let gam = parse_inviscid_gam(&fixture_path("xfoil_inviscid.dat"));
    assert_eq!(gam.len(), f.n + 1, "GAM block length");
    st.gam = gam;
    // xp/yp only enter XICALC's wake-gap cubic; take YFoil's spline for those
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = create_paneled_airfoil(&geom);
    for i in 1..=f.n {
        st.xp[i] = af.xp[i - 1];
        st.yp[i] = af.yp[i - 1];
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
    let af = create_paneled_airfoil(&geom);
    let st = BlState::from_airfoil(&af, f.nw);
    // arc length is pure arithmetic (SCALC): bit-identical
    for i in 1..=f.n {
        assert_bits(st.s[i], f.s[i], &format!("S({i})"));
    }
    assert_eq!(st.sharp, f.sharp, "SHARP");
    // ANTE depends on spline derivatives XP/YP at the TE
    assert_within(st.ante, f.ante, TOL_PURE, 1.0, "ANTE");
    let mut st2 = st.clone();
    tecalc(&mut st2);
    assert_bits(st2.ante, st.ante, "TECALC idempotent");
}

#[test]
fn test_stfind_matches_xfoil() {
    let f = pointers();
    let mut st = state_from_fixture(&f);
    stfind(&mut st);
    assert_eq!(st.ist, f.ist, "IST");
    assert_bits(st.sst, f.sst, "SST");
    assert_bits(st.sst_go, f.sst_go, "SST_GO");
    assert_bits(st.sst_gp, f.sst_gp, "SST_GP");
}

#[test]
fn test_iblpan_and_iblsys_match_xfoil() {
    let f = pointers();
    let mut st = state_from_fixture(&f);
    st.ist = f.ist;
    iblpan(&mut st);
    iblsys(&mut st);
    assert_eq!(st.iblte[1..], f.iblte[1..], "IBLTE");
    assert_eq!(st.nbl[1..], f.nbl[1..], "NBL");
    assert_eq!(st.nsys, f.nsys, "NSYS");
    for is in 1..=2 {
        for ibl in 1..=(f.iblte[is] + f.nw) {
            assert_eq!(st.ipan[is][ibl], f.ipan[is][ibl], "IPAN({ibl},{is})");
            assert_eq!(st.vti[is][ibl].to_bits(), f.vti[is][ibl].to_bits(), "VTI({ibl},{is})");
            assert_eq!(st.isys[is][ibl], f.isys[is][ibl], "ISYS({ibl},{is})");
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
    st.ist = f.ist;
    st.sst = f.sst;
    iblpan(&mut st);
    xicalc(&mut st);
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            assert_bits(st.xssi[is][ibl], f.xssi[is][ibl], &format!("XSSI({ibl},{is})"));
        }
        // upper-side wake mirror (plotting pointers) is filled by XICALC too
        if is == 1 {
            for iw in 1..=f.nw {
                let ibl = f.iblte[1] + iw;
                assert_bits(st.xssi[1][ibl], f.xssi[1][ibl], &format!("XSSI({ibl},1) wake mirror"));
            }
        }
    }
    // WGAP goes through YFoil's XP/YP (spline derivatives): tolerance, scale = ANTE
    for iw in 1..=f.nw {
        assert_within(st.wgap[iw], f.wgap[iw], TOL_PURE, f.ante, &format!("WGAP({iw})"));
    }
}

#[test]
fn test_qiset_and_uicalc_match_xfoil() {
    let f = pointers();
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    assert_eq!((u.n, u.nw), (f.n, f.nw));
    let mut st = state_from_fixture(&f);
    st.ist = f.ist;
    iblpan(&mut st);
    st.qinvu[1] = u.qinvu1.clone();
    st.qinvu[2] = u.qinvu2.clone();
    qiset(&mut st, u.alfa);
    // cos/sin come from libm: same host gives identical bits, other hosts an ULP — TOL_PURE
    for i in 1..=(f.n + f.nw) {
        assert_within(st.qinv[i], u.qinv[i], TOL_PURE, 1.0, &format!("QINV({i})"));
        assert_within(st.qinv_a[i], u.qinv_a[i], TOL_PURE, 1.0, &format!("QINV_A({i})"));
    }
    // UICALC is a sign flip: feed XFOIL's QINV and require bit-identity
    st.qinv = u.qinv.clone();
    st.qinv_a = u.qinv_a.clone();
    uicalc(&mut st);
    for is in 1..=2 {
        assert_eq!(u.uinv[is].len(), f.nbl[is] + 1, "UINV rows side {is}");
        for ibl in 1..=f.nbl[is] {
            assert_bits(st.uinv[is][ibl], u.uinv[is][ibl], &format!("UINV({ibl},{is})"));
            assert_bits(st.uinv_a[is][ibl], u.uinv_a[is][ibl], &format!("UINV_A({ibl},{is})"));
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
    let af = create_paneled_airfoil(&geom);
    let mut st = BlState::from_airfoil(&af, f.nw);
    st.set_wake_nodes(&f.x[f.n + 1..], &f.y[f.n + 1..], &f.s[f.n + 1..]);
    st.qinvu[1] = u.qinvu1.clone();
    st.qinvu[2] = u.qinvu2.clone();
    qiset(&mut st, u.alfa);
    // before the first viscous iteration GAM is the inviscid solution at this alpha
    for i in 1..=f.n {
        st.gam[i] = st.qinv[i];
    }
    stfind(&mut st);
    iblpan(&mut st);
    xicalc(&mut st);
    iblsys(&mut st);
    uicalc(&mut st);

    assert_eq!(st.ist, f.ist, "IST");
    assert_within(st.sst, f.sst, TOL_PURE, 1.0, "SST");
    assert_eq!((st.nbl, st.iblte, st.nsys), (f.nbl, f.iblte, f.nsys));
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            assert_eq!(st.ipan[is][ibl], f.ipan[is][ibl]);
            assert_eq!(st.isys[is][ibl], f.isys[is][ibl]);
            assert_within(
                st.xssi[is][ibl],
                f.xssi[is][ibl],
                TOL_PURE,
                1.0,
                &format!("XSSI({ibl},{is})"),
            );
            assert_within(
                st.uinv[is][ibl],
                u.uinv[is][ibl],
                TOL_PURE,
                1.0,
                &format!("UINV({ibl},{is})"),
            );
        }
    }
    println!(
        "prologue: IST={} NBL={:?} NSYS={} — pointers exact, XSSI/UINV within {TOL_PURE:.0e}",
        st.ist,
        &st.nbl[1..],
        st.nsys
    );
}

/// NCALC / APCALC: YFoil's node normals and panel angles for the airfoil against XFOIL's.
/// Normals come from spline derivatives, angles from atan2 — both held to TOL_PURE.
#[test]
fn test_airfoil_normals_and_panel_angles_match_xfoil() {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = create_paneled_airfoil(&geom);
    for i in 1..=f.n {
        assert_within(af.nx[i - 1], f.nx[i], TOL_PURE, 1.0, &format!("NX({i})"));
        assert_within(af.ny[i - 1], f.ny[i], TOL_PURE, 1.0, &format!("NY({i})"));
        assert_within(af.apanel[i - 1], f.apanel[i], TOL_PURE, 1.0, &format!("APANEL({i})"));
    }
}
