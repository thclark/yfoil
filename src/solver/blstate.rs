//! The boundary-layer state, laid out exactly as XFOIL's COMMON blocks (XFOIL.INC / XBL.INC).
//!
//! Conventions in this module and in `pointers.rs` / `velocity.rs` (CLAUDE.md, Architecture):
//! - **1-based indexing with a dummy slot 0**, so `state.ipan[is][ibl]` reads as `IPAN(IBL,IS)`.
//! - Two sides, `is = 1` (upper) and `is = 2` (lower); the wake is appended to side 2
//!   (`NBL(2) = IBLTE(2) + NW`). Side 1 carries duplicate wake pointers for plotting only.
//! - Panel arrays run `1..=n+nw`: airfoil nodes `1..=n`, wake nodes `n+1..=n+nw`.
//!
//! Nothing here computes anything; the translated subroutines live beside it.

use crate::bl::system::{StationState, Transition};
use crate::geometry::PanelledFoil;

/// BL and panel state (see module docs for indexing).
#[derive(Debug, Clone)]
pub struct SolverState {
    /// Number of airfoil panel nodes (N)
    pub n_foil_nodes: usize,
    /// Number of wake nodes (NW)
    pub n_wake_nodes: usize,

    // ---- panel-level arrays, index 1..=n+nw (airfoil then wake) ------------------------
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub s: Vec<f64>,
    /// Spline derivatives dX/dS, dY/dS at airfoil nodes (1..=n)
    pub dxds: Vec<f64>,
    pub dyds: Vec<f64>,
    /// Node unit normals NX/NY and panel angles APANEL (1..=n+nw; wake part set by XYWAKE)
    pub normal_x: Vec<f64>,
    pub normal_y: Vec<f64>,
    pub panel_angle: Vec<f64>,
    /// Viscous source strengths SIG(I) (1..=n+nw)
    pub sigma: Vec<f64>,
    /// Freestream speed QINF (1.0) and angle of attack ALFA (radians)
    pub qinf: f64,
    pub alpha: f64,
    /// Surface vorticity / tangential velocity GAM(I), GAM_A(I) (1..=n)
    pub gamma: Vec<f64>,
    pub gamma_d_alpha: Vec<f64>,
    /// Inviscid tangential velocity for alpha = 0, 90 deg: `qinvu[1][i]`, `qinvu[2][i]`
    pub q_inviscid_basis: [Vec<f64>; 3],
    /// Inviscid tangential velocity at current alpha and its alpha-derivative (1..=n+nw)
    pub q_inviscid: Vec<f64>,
    pub q_inviscid_d_alpha: Vec<f64>,
    /// Viscous tangential velocity (1..=n+nw)
    pub q_viscous: Vec<f64>,

    // ---- airfoil scalars (GEOPAR / TECALC) ---------------------------------------------
    pub chord: f64,
    pub s_le: f64,
    pub x_le: f64,
    pub y_le: f64,
    pub x_te: f64,
    pub y_te: f64,
    /// TE gap projected normal / parallel to the bisector, total gap, sharp flag (TECALC)
    pub te_thickness_normal: f64,
    pub te_thickness_parallel: f64,
    pub te_gap: f64,
    pub sharp_te: bool,

    // ---- stagnation point (STFIND) ----------------------------------------------------
    pub i_stagnation_node: usize,
    pub s_stagnation: f64,
    pub s_stagnation_d_gamma_node0: f64,
    pub s_stagnation_d_gamma_node1: f64,

    // ---- BL pointer layer (IBLPAN / IBLSYS), sides indexed 1..=2 -----------------------
    pub n_stations: [usize; 3],
    pub i_te_station: [usize; 3],
    pub i_transition_station: [usize; 3],
    pub n_rows: usize,
    pub i_node: [Vec<usize>; 3],
    pub velocity_sign: [Vec<f64>; 3],
    pub i_row: [Vec<usize>; 3],

    // ---- BL station arrays, `[is][ibl]` -----------------------------------------------
    pub xi: [Vec<f64>; 3],
    pub ue: [Vec<f64>; 3],
    pub ue_inviscid: [Vec<f64>; 3],
    pub ue_inviscid_d_alpha: [Vec<f64>; 3],
    pub mass_defect: [Vec<f64>; 3],
    pub theta: [Vec<f64>; 3],
    pub dstar: [Vec<f64>; 3],
    pub sqrtctau: [Vec<f64>; 3],
    pub delta: [Vec<f64>; 3],
    pub thetastar: [Vec<f64>; 3],
    pub us_plot_scale: [Vec<f64>; 3],
    pub tau: [Vec<f64>; 3],
    pub dissipation: [Vec<f64>; 3],
    pub sqrtctaueq: [Vec<f64>; 3],
    /// Wake gap ("dead air" thickness) WGAP(IW), 1..=nw
    pub wake_gap: Vec<f64>,
    /// Source influence matrix DIJ(I,J) = dQtan(I)/dSig(J), 1-based (N+NW)×(N+NW) (QDCALC)
    pub dij: Vec<Vec<f64>>,
    /// Forced-transition x/c per side (XSTRIP); >= 1.0 means free transition
    pub x_trip: [f64; 3],
    /// XSSITR(IS): arc length of transition, TFORCE(IS): transition was forced (set by MRCHDU/MRCHUE)
    pub xi_transition: [f64; 3],
    pub transition_forced: [bool; 3],
    /// COM1/COM2 (the "1" and "2" station COMMON blocks) and XT persist across MRCHUE/MRCHDU/SETBL
    /// calls in XFOIL; the marches take them from here and put them back.
    pub station1: StationState,
    pub station2: StationState,
    /// XT and its XT_* sensitivities (TRCHEK2's COMMON outputs)
    pub transition: Transition,
    /// MINF1/REINF1 (unit-CL values), MATYP/RETYP, and the current MINF/REINF set by MRCL
    pub mach_cl1: f64,
    pub re_cl1: f64,
    pub mach_cl_dependence: usize,
    pub re_cl_dependence: usize,
    /// IDAMP: 0 = envelope e^n (DAMPL), 1 = modified envelope method (DAMPL2), OPER `DAMP`
    pub amplification_model: usize,
    pub mach: f64,
    pub re: f64,
    /// LALFA (fixed alpha; else fixed CL = CLSPEC), CL, CLSPEC
    pub alpha_specified: bool,
    pub cl: f64,
    pub cl_specified: f64,
    /// LBLINI: BL arrays initialised (MRCHUE done)
    pub bl_initialised: bool,
    /// ACRIT(IS), VACCEL, GAMMA
    pub ncrit: [f64; 3],
    pub elimination_threshold: f64,
    pub gamma_gas: f64,
    /// XOCTR/YOCTR/TINDEX (transition x/c, y/c and station index, per side)
    pub x_transition: [f64; 3],
    pub y_transition: [f64; 3],
    pub transition_node_fraction: [f64; 3],
    /// VISCAL-level flags (XFOIL.INC): wake/pointers/DIJ built, viscous mode, converged
    pub wake_built: bool,
    pub pointers_built: bool,
    pub dij_wake_built: bool,
    pub viscous: bool,
    pub converged: bool,
    /// AWAKE/AVISC/MVISC: alpha the wake was built for, alpha and Mach of the converged point
    pub alpha_wake: f64,
    pub alpha_converged: f64,
    pub mach_converged: f64,
    /// TKLAM/TKL_MSQ (COMSET), MINF_CL/REINF_CL (VISCAL's MRCL)
    pub karman_tsien: f64,
    pub karman_tsien_d_machsqd: f64,
    pub mach_d_cl: f64,
    pub re_d_cl: f64,
    /// Force coefficients (CLCALC/CDCALC) and the moment reference point
    pub cm: f64,
    pub cd_pressure: f64,
    pub cd: f64,
    pub cd_friction: f64,
    pub cl_d_alpha: f64,
    pub cl_d_machsqd: f64,
    pub cm_ref_x: f64,
    pub cm_ref_y: f64,
    /// CPI/CPV: inviscid and viscous Cp at the nodes (1..=N+NW)
    pub cp_inviscid: Vec<f64>,
    pub cp_viscous: Vec<f64>,
}

impl SolverState {
    /// Upper bound on BL stations per side: IBLTE(IS) <= N, plus NW wake stations.
    /// (XFOIL sizes these arrays IVX = IQX/2 + IWX + 50; we size to what this case needs.)
    fn ivx(n: usize, nw: usize) -> usize {
        n + nw + 2
    }

    /// An all-zero state with arrays sized for `n` airfoil nodes and `nw` wake nodes.
    pub fn empty(n: usize, nw: usize) -> Self {
        let np = n + nw + 1;
        let ivx = Self::ivx(n, nw);
        let side = || [Vec::new(), vec![0.0; ivx], vec![0.0; ivx]];
        let side_u = || [Vec::new(), vec![0usize; ivx], vec![0usize; ivx]];
        Self {
            n_foil_nodes: n,
            n_wake_nodes: nw,
            x: vec![0.0; np],
            y: vec![0.0; np],
            s: vec![0.0; np],
            dxds: vec![0.0; n + 1],
            dyds: vec![0.0; n + 1],
            normal_x: vec![0.0; np],
            normal_y: vec![0.0; np],
            panel_angle: vec![0.0; np],
            sigma: vec![0.0; np],
            qinf: 1.0,
            alpha: 0.0,
            gamma: vec![0.0; n + 1],
            gamma_d_alpha: vec![0.0; n + 1],
            q_inviscid_basis: [Vec::new(), vec![0.0; np], vec![0.0; np]],
            q_inviscid: vec![0.0; np],
            q_inviscid_d_alpha: vec![0.0; np],
            q_viscous: vec![0.0; np],
            chord: 0.0,
            s_le: 0.0,
            x_le: 0.0,
            y_le: 0.0,
            x_te: 0.0,
            y_te: 0.0,
            te_thickness_normal: 0.0,
            te_thickness_parallel: 0.0,
            te_gap: 0.0,
            sharp_te: false,
            i_stagnation_node: 0,
            s_stagnation: 0.0,
            s_stagnation_d_gamma_node0: 0.0,
            s_stagnation_d_gamma_node1: 0.0,
            n_stations: [0; 3],
            i_te_station: [0; 3],
            i_transition_station: [0; 3],
            n_rows: 0,
            i_node: side_u(),
            velocity_sign: side(),
            i_row: side_u(),
            xi: side(),
            ue: side(),
            ue_inviscid: side(),
            ue_inviscid_d_alpha: side(),
            mass_defect: side(),
            theta: side(),
            dstar: side(),
            sqrtctau: side(),
            delta: side(),
            thetastar: side(),
            us_plot_scale: side(),
            tau: side(),
            dissipation: side(),
            sqrtctaueq: side(),
            wake_gap: vec![0.0; nw + 1],
            dij: Vec::new(),
            x_trip: [0.0, 1.0, 1.0],
            xi_transition: [0.0; 3],
            transition_forced: [false; 3],
            station1: StationState::default(),
            station2: StationState::default(),
            transition: Transition::default(),
            mach_cl1: 0.0,
            re_cl1: 0.0,
            mach_cl_dependence: 1,
            re_cl_dependence: 1,
            amplification_model: 0,
            mach: 0.0,
            re: 0.0,
            alpha_specified: true,
            cl: 0.0,
            cl_specified: 0.0,
            bl_initialised: false,
            ncrit: [0.0, 9.0, 9.0],
            elimination_threshold: 0.01,
            gamma_gas: 1.4,
            x_transition: [0.0; 3],
            y_transition: [0.0; 3],
            transition_node_fraction: [0.0; 3],
            wake_built: false,
            pointers_built: false,
            dij_wake_built: false,
            viscous: false,
            converged: false,
            alpha_wake: 0.0,
            alpha_converged: 0.0,
            mach_converged: 0.0,
            karman_tsien: 0.0,
            karman_tsien_d_machsqd: 0.0,
            mach_d_cl: 0.0,
            re_d_cl: 0.0,
            cm: 0.0,
            cd_pressure: 0.0,
            cd: 0.0,
            cd_friction: 0.0,
            cl_d_alpha: 0.0,
            cl_d_machsqd: 0.0,
            cm_ref_x: 0.25,
            cm_ref_y: 0.0,
            cp_inviscid: vec![0.0; n + nw + 1],
            cp_viscous: vec![0.0; n + nw + 1],
        }
    }

    /// State for a paneled airfoil with `nw` wake nodes to be set later (XYWAKE, stage S3).
    /// Copies the airfoil geometry into the 1-based panel arrays and evaluates TECALC.
    pub fn from_foil(airfoil: &PanelledFoil, nw: usize) -> Self {
        let n = airfoil.n_foil_nodes;
        let mut st = Self::empty(n, nw);
        for i in 1..=n {
            st.x[i] = airfoil.x[i - 1];
            st.y[i] = airfoil.y[i - 1];
            st.s[i] = airfoil.s[i - 1];
            st.dxds[i] = airfoil.dxds[i - 1];
            st.dyds[i] = airfoil.dyds[i - 1];
            st.normal_x[i] = airfoil.normal_x[i - 1];
            st.normal_y[i] = airfoil.normal_y[i - 1];
            st.panel_angle[i] = airfoil.panel_angle[i - 1];
        }
        st.chord = airfoil.chord;
        st.s_le = airfoil.s_le;
        st.x_le = crate::geometry::spline_value(airfoil.s_le, &airfoil.x, &airfoil.dxds, &airfoil.s);
        st.y_le = crate::geometry::spline_value(airfoil.s_le, &airfoil.y, &airfoil.dyds, &airfoil.s);
        st.x_te = 0.5 * (st.x[1] + st.x[n]);
        st.y_te = 0.5 * (st.y[1] + st.y[n]);
        crate::solver::pointers::set_te_thickness(&mut st);
        st
    }

    /// Install wake node coordinates `n+1..=n+nw` (from XYWAKE, or from a fixture).
    pub fn set_wake_nodes(&mut self, x: &[f64], y: &[f64], s: &[f64]) {
        assert_eq!(x.len(), self.n_wake_nodes);
        for iw in 1..=self.n_wake_nodes {
            self.x[self.n_foil_nodes + iw] = x[iw - 1];
            self.y[self.n_foil_nodes + iw] = y[iw - 1];
            self.s[self.n_foil_nodes + iw] = s[iw - 1];
        }
    }
}
