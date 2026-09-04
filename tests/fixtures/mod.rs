#![allow(dead_code)] // shared test-support module; each test crate uses a subset
//! Fixture loading utilities for validation tests
//!
//! This module provides structures and functions for loading XFOIL-generated
//! test fixtures for numerical validation of YFoil.

use mrchdu_fixtures::{parse_bl_dump, BlDump};
use mrchue_fixtures::parse_bl_state;
use pointers_fixtures::{parse_dij, parse_pointers, parse_uinv};
use yfoil::bl::mrchue::mrchue;
use yfoil::bl::system::BLGlobalParams;
use yfoil::solver::blstate::BlState;

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
    pub n_crit: f64,
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
pub struct BLStation {
    pub ibl: usize,
    pub x: f64,
    pub s: f64,
    pub ue: f64,
    pub delta_star: f64,
    pub theta: f64,
    pub hk: f64,
    pub cf: f64,
    pub ctau: f64,
    pub mass: f64,
    pub regime: String,
}

/// BL data for one surface
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLSide {
    pub n_stations: usize,
    pub stations: Vec<BLStation>,
    pub itran: usize,
}

/// Full BL fixture
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BLFixture {
    pub upper: BLSide,
    pub lower: BLSide,
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
            "{}: YFoil={:.10e}, XFOIL={:.10e}, rel_err={:.2e} (tol={:.2e})",
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
#[allow(non_snake_case)]
pub struct WakeInitFixture {
    /// Theta from upper TE
    pub THET_TE1: f64,
    /// Theta from lower TE
    pub THET_TE2: f64,
    /// Delta* from upper TE
    pub DSTR_TE1: f64,
    /// Delta* from lower TE
    pub DSTR_TE2: f64,
    /// Ctau from upper TE
    pub CTAU_TE1: f64,
    /// Ctau from lower TE
    pub CTAU_TE2: f64,
    /// Ue at upper TE
    pub UEDG_TE1: f64,
    /// Ue at lower TE
    pub UEDG_TE2: f64,
    /// TE base thickness (gap)
    pub ANTE: f64,
    /// Combined theta at wake start
    pub TTE: f64,
    /// Combined delta* at wake start (includes ANTE)
    pub DTE: f64,
    /// Weighted Ctau at wake start
    pub CTE: f64,
}

/// Wake station data from XFOIL
#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct WakeStationFixture {
    /// Wake station index (1-based)
    pub IW: i32,
    /// Arc length from LE
    pub XSSI: f64,
    /// Edge velocity
    pub UEDG: f64,
    /// Momentum thickness
    pub THET: f64,
    /// Displacement thickness
    pub DSTR: f64,
    /// Mass defect
    pub MASS: f64,
    /// Ctau (squared shear stress coeff)
    pub CTAU: f64,
    /// Theta (same as THET)
    pub T2: f64,
    /// Delta* (may differ from DSTR by wake thickness)
    pub D2: f64,
    /// Compressibility-corrected Ue
    pub U2: f64,
    /// Kinematic shape factor
    pub HK2: f64,
    /// Energy shape factor
    pub HS2: f64,
    /// Skin friction (0 in wake)
    pub CF2: f64,
    /// Dissipation coefficient
    pub DI2: f64,
}

/// Complete wake fixture
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WakeFixture {
    pub description: String,
    pub airfoil: String,
    pub alpha_deg: f64,
    pub reynolds: f64,
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
pub fn state_before_mrchue() -> (BlState, BLGlobalParams, [f64; 3]) {
    let f = parse_pointers(&require_fixture(&format!("{}/{}", REF_CASE, "xfoil_pointers.dat")), 1);
    let u = parse_uinv(&require_fixture(&format!("{}/{}", REF_CASE, "xfoil_uinv.dat")), 1);
    let d = parse_bl_state(&require_fixture(&format!("{}/{}", REF_CASE, "mrchdu_input_1.dat")));
    let mut st = BlState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    st.ist = f.ist;
    st.sst = f.sst;
    st.nbl = f.nbl;
    st.iblte = f.iblte;
    st.ipan = f.ipan.clone();
    st.vti = f.vti.clone();
    st.isys = f.isys.clone();
    st.xssi = f.xssi.clone();
    st.wgap = f.wgap.clone();
    st.ante = f.ante;
    st.aste = f.aste;
    st.dste = f.dste;
    st.sharp = f.sharp;
    st.chord = f.chord;
    st.sle = f.sle;
    st.xle = f.xle;
    st.yle = f.yle;
    st.xte = f.xte;
    st.yte = f.yte;
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            st.uinv[is][ibl] = u.uinv[is][ibl];
            st.uedg[is][ibl] = u.uinv[is][ibl];
        }
    }
    let params = BLGlobalParams::new(d.minf, d.reinf, 1.4);
    // the BL parameters SETBL derives must match the reference bitwise before marching
    assert_eq!(params.reybl.to_bits(), d.reybl.to_bits(), "REYBL");
    assert_eq!(params.hstinv.to_bits(), d.hstinv.to_bits(), "HSTINV");
    assert_eq!(params.gm1.to_bits(), d.gm1bl.to_bits(), "GM1BL");
    (st, params, [0.0, d.acrit[1], d.acrit[2]])
}

/// XFOIL's complete state entering the BL march inside SETBL call `k` (after the parameter
/// setup and, on call 1, after MRCHUE): the pointer layer (geometry, IPAN, WGAP), DIJ,
/// UINV/UINV_A, the BL arrays exactly as SETBL had them, that iteration's SST/SST_GO/SST_GP and
/// the SETBL control flags. On call 1 the COMMON state (COM1/COM2/XT) is what MRCHUE left
/// behind, obtained by replaying MRCHUE from its own gated input.
#[allow(dead_code)]
pub fn state_before_setbl_march(k: usize) -> (BlState, BLGlobalParams, BlDump) {
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
    let mut st = BlState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    // XP/YP are not dumped; they are SEGSPL of the dumped X/Y/S (identical construction to
    // create_paneled_airfoil / BlState::from_airfoil)
    {
        let n = f.n;
        let xp = yfoil::geometry::spline(&f.x[1..=n], &f.s[1..=n]);
        let yp = yfoil::geometry::spline(&f.y[1..=n], &f.s[1..=n]);
        st.xp[1..=n].copy_from_slice(&xp);
        st.yp[1..=n].copy_from_slice(&yp);
    }
    st.nx = f.nx.clone();
    st.ny = f.ny.clone();
    st.apanel = f.apanel.clone();
    st.nbl = f.nbl;
    st.iblte = f.iblte;
    st.nsys = f.nsys;
    st.ipan = f.ipan.clone();
    st.vti = f.vti.clone();
    st.isys = f.isys.clone();
    st.wgap = f.wgap.clone();
    st.ante = f.ante;
    st.aste = f.aste;
    st.dste = f.dste;
    st.sharp = f.sharp;
    st.chord = f.chord;
    st.sle = f.sle;
    st.xle = f.xle;
    st.yle = f.yle;
    st.xte = f.xte;
    st.yte = f.yte;
    st.dij = dij;
    // this iteration's stagnation point
    st.ist = d.int("IST");
    st.sst = d.real("SST");
    st.sst_go = d.real("SST_GO");
    st.sst_gp = d.real("SST_GP");
    // SETBL controls
    st.lalfa = d.logical("LALFA");
    st.matyp = d.int("MATYP");
    st.retyp = d.int("RETYP");
    st.minf1 = d.real("MINF");
    st.reinf1 = d.real("REINF");
    st.cl = d.real("CLMR");
    st.clspec = d.real("CLMR");
    st.qinf = d.real("QINF");
    st.alfa = u.alfa;
    st.vaccel = d.real("VACCEL");
    st.acrit = [0.0, d.real("ACRIT1"), d.real("ACRIT2")];
    st.xstrip = [0.0, d.real("XSTRIP1"), d.real("XSTRIP2")];
    st.itran = [0, d.int("ITRAN1"), d.int("ITRAN2")];
    st.lblini = true;
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            let r = d.bl[is][ibl];
            st.xssi[is][ibl] = r[0];
            st.uedg[is][ibl] = r[1];
            st.thet[is][ibl] = r[2];
            st.dstr[is][ibl] = r[3];
            st.ctau[is][ibl] = r[4];
            st.mass[is][ibl] = r[5];
            st.uinv[is][ibl] = u.uinv[is][ibl];
            st.uinv_a[is][ibl] = u.uinv_a[is][ibl];
            // STMOVE recomputes XSSI every VISCAL iteration (SST moves); the pointer dump is
            // from the first IBLSYS call, so only call 1 can be checked against it.
            if k == 1 {
                assert_eq!(f.xssi[is][ibl].to_bits(), r[0].to_bits(), "call {k}: XSSI({ibl},{is})");
            }
        }
    }
    let params = BLGlobalParams::new(d.real("MINF"), d.real("REINF"), 1.4);
    assert_eq!(params.reybl.to_bits(), d.real("REYBL").to_bits(), "REYBL");
    assert_eq!(params.hstinv.to_bits(), d.real("HSTINV").to_bits(), "HSTINV");
    assert_eq!(params.gm1.to_bits(), d.real("GM1BL").to_bits(), "GM1BL");
    if k == 1 {
        let (mut pre, p2, a2) = state_before_mrchue();
        mrchue(&mut pre, &p2, a2, None);
        st.com1 = pre.com1;
        st.com2 = pre.com2;
        st.trloc = pre.trloc.clone();
    }
    (st, params, d)
}
