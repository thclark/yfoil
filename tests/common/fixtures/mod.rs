#![allow(dead_code)] // shared test-support module; each test crate uses a subset
//! Fixture loading utilities for validation tests
//!
//! This module provides structures and functions for loading XFOIL-generated
//! test fixtures for numerical validation of yFoil.

use mrchdu_fixtures::{parse_bl_dump, BlDump};
use mrchue_fixtures::parse_bl_state;
use pointers_fixtures::{parse_dij, parse_pointers, parse_uinv};
use yfoil::bl::mrchue::march_direct;
use yfoil::bl::system::FlowParameters;
use yfoil::solver::blstate::SolverState;

pub mod blsolv_fixtures;
pub mod cases;
pub mod mrchdu_fixtures;
pub mod mrchue_fixtures;
pub mod pointers_fixtures;
pub mod tgap_fixtures;

/// Assert that two values match within tolerance
pub fn assert_xfoil_match(name: &str, yfoil_val: f64, xfoil_val: f64, rel_tol: f64, abs_tol: f64) {
    let diff = (yfoil_val - xfoil_val).abs();
    let max_val = xfoil_val.abs().max(abs_tol);
    let rel_err = diff / max_val;

    if rel_err > rel_tol && diff > abs_tol {
        panic!(
            "{}: yFoil={:.10e}, XFOIL={:.10e}, rel_err={:.2e} (tol={:.2e})",
            name, yfoil_val, xfoil_val, rel_err, rel_tol
        );
    }
}

/// Resolve a fixture path under the tracked reference case and **fail** if it is missing.
/// Never skip silently (CLAUDE.md Rule 7). Regenerate with `cargo xtask fixtures`.
pub fn require_fixture(rel: &str) -> std::path::PathBuf {
    let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel);
    let producer = if rel.contains("fixtures/naca456") {
        "`scripts/naca456-fixtures.sh`"
    } else if rel.contains("fixtures/repanel_cosine") {
        "`UPDATE_SNAPSHOT=1 cargo test --test application repanel_cosine`"
    } else if rel.contains("fixtures/xfoil/") || rel.contains("fixtures/subroutines") {
        "`cargo xtask fixtures` (xtask/fixtures-config/cases.toml)"
    } else {
        "a hand-written file — see tests/fixtures/README.md"
    };
    assert!(
        p.exists(),
        "fixture missing: {} — its producer is {producer}",
        p.display()
    );
    p
}

/// The tracked CI reference case (NACA 0012, N=60, alpha=2, Re=1e6, M=0, Ncrit=9, ITER 20).
pub const REF_CASE: &str = "tests/fixtures/xfoil/naca0012_n60_a2_re1e6";

/// XFOIL's complete state at the start of MRCHUE on the reference case: the pointer layer,
/// UEDG = UINV, and the BL parameters SETBL derives (asserted bitwise against the dump).
#[allow(dead_code)]
pub fn state_before_mrchue() -> (SolverState, FlowParameters, [f64; 3]) {
    let f = parse_pointers(&require_fixture(&format!("{}/{}", REF_CASE, "xfoil_pointers.dat")), 1);
    let u = parse_uinv(&require_fixture(&format!("{}/{}", REF_CASE, "xfoil_uinv.dat")), 1);
    let d = parse_bl_state(&require_fixture(&format!("{}/{}", REF_CASE, "mrchdu_input_1.dat")));
    let mut st = SolverState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    st.i_stagnation_node = f.ist;
    st.s_stagnation = f.sst;
    st.n_stations = f.nbl;
    st.i_te_station = f.iblte;
    st.i_node = f.ipan.clone();
    st.velocity_sign = f.vti.clone();
    st.i_row = f.isys.clone();
    st.xi = f.xssi.clone();
    st.wake_gap = f.wgap.clone();
    st.te_thickness_normal = f.ante;
    st.te_thickness_parallel = f.aste;
    st.te_gap = f.dste;
    st.sharp_te = f.sharp;
    st.chord = f.chord;
    st.s_le = f.sle;
    st.x_le = f.xle;
    st.y_le = f.yle;
    st.x_te = f.xte;
    st.y_te = f.yte;
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            st.ue_inviscid[is][ibl] = u.uinv[is][ibl];
            st.ue[is][ibl] = u.uinv[is][ibl];
        }
    }
    let params = FlowParameters::new(d.minf, d.reinf, 1.4);
    // the BL parameters SETBL derives must match the reference bitwise before marching
    assert_eq!(params.re.to_bits(), d.reybl.to_bits(), "REYBL");
    assert_eq!(params.h_stagnation_inv.to_bits(), d.hstinv.to_bits(), "HSTINV");
    assert_eq!(params.gamma_gas_m1.to_bits(), d.gm1bl.to_bits(), "GM1BL");
    (st, params, [0.0, d.acrit[1], d.acrit[2]])
}

/// XFOIL's complete state entering the BL march inside SETBL call `k` (after the parameter
/// setup and, on call 1, after MRCHUE): the pointer layer (geometry, IPAN, WGAP), DIJ,
/// UINV/UINV_A, the BL arrays exactly as SETBL had them, that iteration's SST/SST_GO/SST_GP and
/// the SETBL control flags. On call 1 the COMMON state (COM1/COM2/XT) is what MRCHUE left
/// behind, obtained by replaying MRCHUE from its own gated input.
#[allow(dead_code)]
pub fn state_before_setbl_march(k: usize) -> (SolverState, FlowParameters, BlDump) {
    let fx = |name: &str| require_fixture(&format!("{}/{}", REF_CASE, name));
    let f = parse_pointers(&fx("xfoil_pointers.dat"), 1);
    let u = parse_uinv(&fx("xfoil_uinv.dat"), 1);
    let (n, nw, dij) = parse_dij(&fx("xfoil_dij.dat"));
    assert_eq!((n, nw), (f.n, f.nw), "DIJ dump vs pointer dump size");
    let d = parse_bl_dump(&fx(&format!("mrchdu_input_{k}.dat")));
    assert_eq!(
        [d.int("NBL1"), d.int("NBL2")],
        [f.nbl[1], f.nbl[2]],
        "call {k}: NBL vs pointer dump"
    );
    assert_eq!(
        [d.int("IBLTE1"), d.int("IBLTE2")],
        [f.iblte[1], f.iblte[2]],
        "call {k}: IBLTE"
    );
    let mut st = SolverState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    // XP/YP are not dumped; they are SEGSPL of the dumped X/Y/S (identical construction to
    // panel_foil / SolverState::from_foil)
    {
        let n = f.n;
        let xp = yfoil::geometry::spline_derivatives(&f.x[1..=n], &f.s[1..=n]);
        let yp = yfoil::geometry::spline_derivatives(&f.y[1..=n], &f.s[1..=n]);
        st.dxds[1..=n].copy_from_slice(&xp);
        st.dyds[1..=n].copy_from_slice(&yp);
    }
    st.normal_x = f.nx.clone();
    st.normal_y = f.ny.clone();
    st.panel_angle = f.apanel.clone();
    st.n_stations = f.nbl;
    st.i_te_station = f.iblte;
    st.n_rows = f.nsys;
    st.i_node = f.ipan.clone();
    st.velocity_sign = f.vti.clone();
    st.i_row = f.isys.clone();
    st.wake_gap = f.wgap.clone();
    st.te_thickness_normal = f.ante;
    st.te_thickness_parallel = f.aste;
    st.te_gap = f.dste;
    st.sharp_te = f.sharp;
    st.chord = f.chord;
    st.s_le = f.sle;
    st.x_le = f.xle;
    st.y_le = f.yle;
    st.x_te = f.xte;
    st.y_te = f.yte;
    st.dij = dij;
    // this iteration's stagnation point
    st.i_stagnation_node = d.int("IST");
    st.s_stagnation = d.real("SST");
    st.s_stagnation_d_gamma_node0 = d.real("SST_GO");
    st.s_stagnation_d_gamma_node1 = d.real("SST_GP");
    // SETBL controls
    st.alpha_specified = d.logical("LALFA");
    st.mach_cl_dependence = yfoil::bl::system::MachClDependence::from_xfoil(d.int("MATYP"));
    st.re_cl_dependence = yfoil::bl::system::ReClDependence::from_xfoil(d.int("RETYP"));
    st.mach_cl1 = d.real("MINF");
    st.re_cl1 = d.real("REINF");
    st.cl = d.real("CLMR");
    st.cl_specified = d.real("CLMR");
    st.qinf = d.real("QINF");
    st.alpha = u.alfa;
    st.elimination_threshold = d.real("VACCEL");
    st.ncrit = [0.0, d.real("ACRIT1"), d.real("ACRIT2")];
    st.x_trip = [0.0, d.real("XSTRIP1"), d.real("XSTRIP2")];
    st.i_transition_station = [0, d.int("ITRAN1"), d.int("ITRAN2")];
    st.bl_initialised = true;
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            let r = d.bl[is][ibl];
            st.xi[is][ibl] = r[0];
            st.ue[is][ibl] = r[1];
            st.theta[is][ibl] = r[2];
            st.dstar[is][ibl] = r[3];
            st.sqrtctau[is][ibl] = r[4];
            st.mass_defect[is][ibl] = r[5];
            st.ue_inviscid[is][ibl] = u.uinv[is][ibl];
            st.ue_inviscid_d_alpha[is][ibl] = u.uinv_a[is][ibl];
            // STMOVE recomputes XSSI every VISCAL iteration (SST moves); the pointer dump is
            // from the first IBLSYS call, so only call 1 can be checked against it.
            if k == 1 {
                assert_eq!(f.xssi[is][ibl].to_bits(), r[0].to_bits(), "call {k}: XSSI({ibl},{is})");
            }
        }
    }
    let params = FlowParameters::new(d.real("MINF"), d.real("REINF"), 1.4);
    assert_eq!(params.re.to_bits(), d.real("REYBL").to_bits(), "REYBL");
    assert_eq!(params.h_stagnation_inv.to_bits(), d.real("HSTINV").to_bits(), "HSTINV");
    assert_eq!(params.gamma_gas_m1.to_bits(), d.real("GM1BL").to_bits(), "GM1BL");
    if k == 1 {
        let (mut pre, p2, a2) = state_before_mrchue();
        march_direct(&mut pre, &p2, a2, None);
        st.station1 = pre.station1;
        st.station2 = pre.station2;
        st.transition = pre.transition.clone();
    }
    (st, params, d)
}
