//! S3 gate: XYWAKE (wake trajectory by streamline tracing with PSILIN), SETEXP and QWCALC
//! against the wake nodes and wake QINVU recorded by the instrumented reference.
//!
//! Inputs are XFOIL's own airfoil arrays (nodes, normals, panel angles, GAM/GAMU at the
//! first VISCAL call) so that the wake code is tested in isolation; YFoil's spline provides
//! only XP/YP at the two TE nodes (the wake's initial direction). Everything goes through
//! libm (log, atan2, sqrt) and a 100-iteration Newton in SETEXP, so the gate is TOL_PURE.

mod fixtures;
mod utilities;

use fixtures::pointers_fixtures::{parse_inviscid_gam, parse_pointers, parse_uinv};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_PURE};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::xywake::{build_wake, exponential_spacing, set_wake_q_basis};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn state_at_first_viscal_call() -> (
    SolverState,
    fixtures::pointers_fixtures::PointersFixture,
    fixtures::pointers_fixtures::UinvFixture,
) {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    let mut st = SolverState::empty(f.n, f.nw);
    for i in 1..=f.n {
        st.x[i] = f.x[i];
        st.y[i] = f.y[i];
        st.s[i] = f.s[i];
        st.normal_x[i] = f.nx[i];
        st.normal_y[i] = f.ny[i];
        st.panel_angle[i] = f.apanel[i];
    }
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = panel_foil(&geom);
    for i in 1..=f.n {
        st.dxds[i] = af.xp[i - 1];
        st.dyds[i] = af.yp[i - 1];
    }
    st.chord = f.chord;
    st.te_thickness_normal = f.ante;
    st.te_thickness_parallel = f.aste;
    st.te_gap = f.dste;
    st.sharp_te = f.sharp;
    st.gamma = parse_inviscid_gam(&fixture_path("xfoil_inviscid.dat"));
    st.q_inviscid_basis[1] = u.qinvu1.clone();
    st.q_inviscid_basis[2] = u.qinvu2.clone();
    st.alpha = u.alfa;
    st.qinf = 1.0;
    (st, f, u)
}

#[test]
fn test_setexp_reproduces_wake_spacing() {
    // The dumped wake S values are S(N) + SNEW; back out SNEW and compare with the port.
    let (st, f, _) = state_at_first_viscal_call();
    let n = f.n;
    let ds1 = 0.5 * (st.s[2] - st.s[1] + st.s[n] - st.s[n - 1]);
    let snew = exponential_spacing(ds1, 1.0 * st.chord, f.nw);
    for iw in 1..=f.nw {
        let xfoil = f.s[n + iw] - f.s[n];
        assert_within(snew[iw], xfoil, TOL_PURE, 1.0, &format!("SNEW({})", n + iw));
    }
}

#[test]
fn test_xywake_matches_xfoil_wake_nodes() {
    let (mut st, f, _) = state_at_first_viscal_call();
    build_wake(&mut st, 1.0);
    let n = f.n;
    for iw in 1..=f.nw {
        let i = n + iw;
        assert_within(st.x[i], f.x[i], TOL_PURE, 1.0, &format!("X({i})"));
        assert_within(st.y[i], f.y[i], TOL_PURE, 1.0, &format!("Y({i})"));
        assert_within(st.s[i], f.s[i], TOL_PURE, 1.0, &format!("S({i})"));
        assert_within(st.normal_x[i], f.nx[i], TOL_PURE, 1.0, &format!("NX({i})"));
        assert_within(st.normal_y[i], f.ny[i], TOL_PURE, 1.0, &format!("NY({i})"));
        // APANEL is not set at the last wake node (XYWAKE exits the loop before it)
        if iw < f.nw {
            assert_within(st.panel_angle[i], f.apanel[i], TOL_PURE, 1.0, &format!("APANEL({i})"));
        }
    }
    println!("xywake: {} wake nodes within {TOL_PURE:.0e} of XFOIL", f.nw);
}

#[test]
fn test_qwcalc_matches_xfoil_wake_qinvu() {
    let (mut st, f, u) = state_at_first_viscal_call();
    build_wake(&mut st, 1.0);
    set_wake_q_basis(&mut st);
    for i in (f.n + 1)..=(f.n + f.nw) {
        assert_within(
            st.q_inviscid_basis[1][i],
            u.qinvu1[i],
            TOL_PURE,
            1.0,
            &format!("QINVU({i},1)"),
        );
        assert_within(
            st.q_inviscid_basis[2][i],
            u.qinvu2[i],
            TOL_PURE,
            1.0,
            &format!("QINVU({i},2)"),
        );
    }
}
