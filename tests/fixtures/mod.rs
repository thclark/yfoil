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
pub mod mrchdu_fixtures;
pub mod mrchue_fixtures;
pub mod pointers_fixtures;

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

/// Test case input parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestInput {
    pub airfoil: String,
    pub alpha_deg: f64,
    pub reynolds: f64,
    pub mach: f64,
    pub n_panels: usize,
    pub ncrit: f64,
}

/// Panel geometry from XFOIL
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryFixture {
    pub n: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub s: Vec<f64>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
}

/// Inviscid solution from XFOIL
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InviscidFixture {
    pub alpha_rad: f64,
    pub gamma: Vec<f64>,
    pub qinv: Vec<f64>,
    pub cpi: Vec<f64>,
    pub cl_inv: f64,
    pub cm_inv: f64,
}

/// BL station data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FixtureStation {
    pub ibl: usize,
    pub x: f64,
    pub s: f64,
    pub ue: f64,
    pub dstar: f64,
    pub theta: f64,
    pub hk: f64,
    pub cf: f64,
    pub ctau: f64,
    pub mass: f64,
    pub regime: String,
}

/// BL data for one surface
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FixtureSide {
    pub n_stations: usize,
    pub stations: Vec<FixtureStation>,
    pub itran: usize,
}

/// Full BL fixture
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLFixture {
    pub upper: FixtureSide,
    pub lower: FixtureSide,
}

/// VISCAL iteration data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViscalIteration {
    pub iter: usize,
    pub alpha_deg: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cdp: f64,
    pub cm: f64,
    pub rmsbl: f64,
    pub rmxbl: f64,
    pub rlx: f64,
}

/// VISCAL iteration history
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViscalFixture {
    pub n_iterations: usize,
    pub converged: bool,
    pub iterations: Vec<ViscalIteration>,
}

/// Final converged results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FinalFixture {
    pub alpha_deg: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cdp: f64,
    pub cm: f64,
    pub xtr_upper: f64,
    pub xtr_lower: f64,
}

/// Complete test fixture
#[derive(Debug, Clone)]
pub struct TestFixture {
    pub input: TestInput,
    pub geometry: Option<GeometryFixture>,
    pub inviscid: Option<InviscidFixture>,
    pub bl_stations: Option<BLFixture>,
    pub viscal_iters: Option<ViscalFixture>,
    pub final_result: Option<FinalFixture>,
}

impl TestFixture {
    /// Load a test fixture from a directory
    pub fn load(dir: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let input: TestInput = load_json(&dir.join("input.json"))?;

        let geometry = load_json_opt(&dir.join("geometry.json"))?;
        let inviscid = load_json_opt(&dir.join("inviscid.json"))?;
        let bl_stations = load_json_opt(&dir.join("bl_stations.json"))?;
        let viscal_iters = load_json_opt(&dir.join("viscal_iters.json"))?;
        let final_result = load_json_opt(&dir.join("final.json"))?;

        Ok(TestFixture {
            input,
            geometry,
            inviscid,
            bl_stations,
            viscal_iters,
            final_result,
        })
    }
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(path)?;
    let data: T = serde_json::from_str(&content)?;
    Ok(data)
}

fn load_json_opt<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<Option<T>, Box<dyn std::error::Error>> {
    if path.exists() {
        Ok(Some(load_json(path)?))
    } else {
        Ok(None)
    }
}

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

/// Default relative tolerance for XFOIL comparisons
pub const DEFAULT_REL_TOL: f64 = 1e-10;

/// Default absolute tolerance for XFOIL comparisons
pub const DEFAULT_ABS_TOL: f64 = 1e-14;

/// Wake initial conditions from XFOIL
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeInitFixture {
    /// Theta from upper TE
    pub theta_te_station1: f64,
    /// Theta from lower TE
    pub theta_te_station2: f64,
    /// Delta* from upper TE
    pub dstar_te_station1: f64,
    /// Delta* from lower TE
    pub dstar_te_station2: f64,
    /// Ctau from upper TE
    pub sqrtctau_te_station1: f64,
    /// Ctau from lower TE
    pub sqrtctau_te_station2: f64,
    /// Ue at upper TE
    pub ue_te_station1: f64,
    /// Ue at lower TE
    pub ue_te_station2: f64,
    /// TE base thickness (gap)
    pub te_thickness_normal: f64,
    /// Combined theta at wake start
    pub theta_te: f64,
    /// Combined delta* at wake start (includes ANTE)
    pub dstar_te: f64,
    /// Weighted Ctau at wake start
    pub sqrtctau_te: f64,
}

/// Wake station data from XFOIL
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeStationFixture {
    /// Wake station index (1-based)
    pub i_wake: i32,
    /// Arc length from LE
    pub xi: f64,
    /// Edge velocity
    pub ue: f64,
    /// Momentum thickness
    pub theta: f64,
    /// Displacement thickness
    pub dstar: f64,
    /// Mass defect
    pub mass_defect: f64,
    /// Ctau (squared shear stress coeff)
    pub sqrtctau: f64,
    /// Theta (same as THET)
    pub theta_station2: f64,
    /// Delta* (may differ from DSTR by wake thickness)
    pub dstar_station2: f64,
    /// Compressibility-corrected Ue
    pub ue_station2: f64,
    /// Kinematic shape factor
    pub hk_station2: f64,
    /// Energy shape factor
    pub hstar_station2: f64,
    /// Skin friction (0 in wake)
    pub cf_station2: f64,
    /// Dissipation coefficient
    pub cdiss_station2: f64,
}

/// Complete wake fixture
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeFixture {
    pub description: String,
    pub foil: String,
    pub alpha_deg: f64,
    pub re: f64,
    pub mach: f64,
    pub wake_init: WakeInitFixture,
    pub wake_stations: Vec<WakeStationFixture>,
}

impl WakeFixture {
    /// Load a wake fixture from a JSON file
    pub fn load(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        load_json(path)
    }
}

/// Resolve a fixture path under the tracked reference case and **fail** if it is missing.
/// Never skip silently (CLAUDE.md Rule 7). Regenerate with `cargo xtask fixtures`.
pub fn require_fixture(rel: &str) -> std::path::PathBuf {
    let p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel);
    assert!(
        p.exists(),
        "fixture missing: {} — regenerate with `cargo xtask fixtures` (see xtask/fixtures-config/cases.toml)",
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
