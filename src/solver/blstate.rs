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
use crate::geometry::PaneledAirfoil;

/// BL and panel state (see module docs for indexing).
#[derive(Debug, Clone)]
pub struct BlState {
    /// Number of airfoil panel nodes (N)
    pub n: usize,
    /// Number of wake nodes (NW)
    pub nw: usize,

    // ---- panel-level arrays, index 1..=n+nw (airfoil then wake) ------------------------
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub s: Vec<f64>,
    /// Spline derivatives dX/dS, dY/dS at airfoil nodes (1..=n)
    pub xp: Vec<f64>,
    pub yp: Vec<f64>,
    /// Node unit normals NX/NY and panel angles APANEL (1..=n+nw; wake part set by XYWAKE)
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub apanel: Vec<f64>,
    /// Viscous source strengths SIG(I) (1..=n+nw)
    pub sig: Vec<f64>,
    /// Freestream speed QINF (1.0) and angle of attack ALFA (radians)
    pub qinf: f64,
    pub alfa: f64,
    /// Surface vorticity / tangential velocity GAM(I), GAM_A(I) (1..=n)
    pub gam: Vec<f64>,
    pub gam_a: Vec<f64>,
    /// Inviscid tangential velocity for alpha = 0, 90 deg: `qinvu[1][i]`, `qinvu[2][i]`
    pub qinvu: [Vec<f64>; 3],
    /// Inviscid tangential velocity at current alpha and its alpha-derivative (1..=n+nw)
    pub qinv: Vec<f64>,
    pub qinv_a: Vec<f64>,
    /// Viscous tangential velocity (1..=n+nw)
    pub qvis: Vec<f64>,

    // ---- airfoil scalars (GEOPAR / TECALC) ---------------------------------------------
    pub chord: f64,
    pub sle: f64,
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    /// TE gap projected normal / parallel to the bisector, total gap, sharp flag (TECALC)
    pub ante: f64,
    pub aste: f64,
    pub dste: f64,
    pub sharp: bool,

    // ---- stagnation point (STFIND) ----------------------------------------------------
    pub ist: usize,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,

    // ---- BL pointer layer (IBLPAN / IBLSYS), sides indexed 1..=2 -----------------------
    pub nbl: [usize; 3],
    pub iblte: [usize; 3],
    pub itran: [usize; 3],
    pub nsys: usize,
    pub ipan: [Vec<usize>; 3],
    pub vti: [Vec<f64>; 3],
    pub isys: [Vec<usize>; 3],

    // ---- BL station arrays, `[is][ibl]` -----------------------------------------------
    pub xssi: [Vec<f64>; 3],
    pub uedg: [Vec<f64>; 3],
    pub uinv: [Vec<f64>; 3],
    pub uinv_a: [Vec<f64>; 3],
    pub mass: [Vec<f64>; 3],
    pub thet: [Vec<f64>; 3],
    pub dstr: [Vec<f64>; 3],
    pub ctau: [Vec<f64>; 3],
    pub delt: [Vec<f64>; 3],
    pub tstr: [Vec<f64>; 3],
    pub uslp: [Vec<f64>; 3],
    pub guxq: [Vec<f64>; 3],
    pub guxd: [Vec<f64>; 3],
    pub tau: [Vec<f64>; 3],
    pub dis: [Vec<f64>; 3],
    pub ctq: [Vec<f64>; 3],
    /// Wake gap ("dead air" thickness) WGAP(IW), 1..=nw
    pub wgap: Vec<f64>,
    /// Source influence matrix DIJ(I,J) = dQtan(I)/dSig(J), 1-based (N+NW)×(N+NW) (QDCALC)
    pub dij: Vec<Vec<f64>>,
    /// Forced-transition x/c per side (XSTRIP); >= 1.0 means free transition
    pub xstrip: [f64; 3],
    /// XSSITR(IS): arc length of transition, TFORCE(IS): transition was forced (set by MRCHDU/MRCHUE)
    pub xssitr: [f64; 3],
    pub tforce: [bool; 3],
    /// COM1/COM2 (the "1" and "2" station COMMON blocks) and XT persist across MRCHUE/MRCHDU/SETBL
    /// calls in XFOIL; the marches take them from here and put them back.
    pub com1: StationState,
    pub com2: StationState,
    /// XT and its XT_* sensitivities (TRCHEK2's COMMON outputs)
    pub trloc: Transition,
    /// MINF1/REINF1 (unit-CL values), MATYP/RETYP, and the current MINF/REINF set by MRCL
    pub minf1: f64,
    pub reinf1: f64,
    pub matyp: usize,
    pub retyp: usize,
    /// IDAMP: 0 = envelope e^n (DAMPL), 1 = modified envelope method (DAMPL2), OPER `DAMP`
    pub idamp: usize,
    pub minf: f64,
    pub reinf: f64,
    /// LALFA (fixed alpha; else fixed CL = CLSPEC), CL, CLSPEC
    pub lalfa: bool,
    pub cl: f64,
    pub clspec: f64,
    /// LBLINI: BL arrays initialised (MRCHUE done)
    pub lblini: bool,
    /// ACRIT(IS), VACCEL, GAMMA
    pub acrit: [f64; 3],
    pub vaccel: f64,
    pub gamma: f64,
    /// XOCTR/YOCTR/TINDEX (transition x/c, y/c and station index, per side)
    pub xoctr: [f64; 3],
    pub yoctr: [f64; 3],
    pub tindex: [f64; 3],
    /// VISCAL-level flags (XFOIL.INC): wake/pointers/DIJ built, viscous mode, converged
    pub lwake: bool,
    pub lipan: bool,
    pub lwdij: bool,
    pub lvisc: bool,
    pub lvconv: bool,
    /// AWAKE/AVISC/MVISC: alpha the wake was built for, alpha and Mach of the converged point
    pub awake: f64,
    pub avisc: f64,
    pub mvisc: f64,
    /// TKLAM/TKL_MSQ (COMSET), MINF_CL/REINF_CL (VISCAL's MRCL)
    pub tklam: f64,
    pub tkl_msq: f64,
    pub minf_cl: f64,
    pub reinf_cl: f64,
    /// Force coefficients (CLCALC/CDCALC) and the moment reference point
    pub cm: f64,
    pub cdp: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cl_alf: f64,
    pub cl_msq: f64,
    pub xcmref: f64,
    pub ycmref: f64,
    /// CPI/CPV: inviscid and viscous Cp at the nodes (1..=N+NW)
    pub cpi: Vec<f64>,
    pub cpv: Vec<f64>,
}

impl BlState {
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
            n,
            nw,
            x: vec![0.0; np],
            y: vec![0.0; np],
            s: vec![0.0; np],
            xp: vec![0.0; n + 1],
            yp: vec![0.0; n + 1],
            nx: vec![0.0; np],
            ny: vec![0.0; np],
            apanel: vec![0.0; np],
            sig: vec![0.0; np],
            qinf: 1.0,
            alfa: 0.0,
            gam: vec![0.0; n + 1],
            gam_a: vec![0.0; n + 1],
            qinvu: [Vec::new(), vec![0.0; np], vec![0.0; np]],
            qinv: vec![0.0; np],
            qinv_a: vec![0.0; np],
            qvis: vec![0.0; np],
            chord: 0.0,
            sle: 0.0,
            xle: 0.0,
            yle: 0.0,
            xte: 0.0,
            yte: 0.0,
            ante: 0.0,
            aste: 0.0,
            dste: 0.0,
            sharp: false,
            ist: 0,
            sst: 0.0,
            sst_go: 0.0,
            sst_gp: 0.0,
            nbl: [0; 3],
            iblte: [0; 3],
            itran: [0; 3],
            nsys: 0,
            ipan: side_u(),
            vti: side(),
            isys: side_u(),
            xssi: side(),
            uedg: side(),
            uinv: side(),
            uinv_a: side(),
            mass: side(),
            thet: side(),
            dstr: side(),
            ctau: side(),
            delt: side(),
            tstr: side(),
            uslp: side(),
            guxq: side(),
            guxd: side(),
            tau: side(),
            dis: side(),
            ctq: side(),
            wgap: vec![0.0; nw + 1],
            dij: Vec::new(),
            xstrip: [0.0, 1.0, 1.0],
            xssitr: [0.0; 3],
            tforce: [false; 3],
            com1: StationState::default(),
            com2: StationState::default(),
            trloc: Transition::default(),
            minf1: 0.0,
            reinf1: 0.0,
            matyp: 1,
            retyp: 1,
            idamp: 0,
            minf: 0.0,
            reinf: 0.0,
            lalfa: true,
            cl: 0.0,
            clspec: 0.0,
            lblini: false,
            acrit: [0.0, 9.0, 9.0],
            vaccel: 0.01,
            gamma: 1.4,
            xoctr: [0.0; 3],
            yoctr: [0.0; 3],
            tindex: [0.0; 3],
            lwake: false,
            lipan: false,
            lwdij: false,
            lvisc: false,
            lvconv: false,
            awake: 0.0,
            avisc: 0.0,
            mvisc: 0.0,
            tklam: 0.0,
            tkl_msq: 0.0,
            minf_cl: 0.0,
            reinf_cl: 0.0,
            cm: 0.0,
            cdp: 0.0,
            cd: 0.0,
            cdf: 0.0,
            cl_alf: 0.0,
            cl_msq: 0.0,
            xcmref: 0.25,
            ycmref: 0.0,
            cpi: vec![0.0; n + nw + 1],
            cpv: vec![0.0; n + nw + 1],
        }
    }

    /// State for a paneled airfoil with `nw` wake nodes to be set later (XYWAKE, stage S3).
    /// Copies the airfoil geometry into the 1-based panel arrays and evaluates TECALC.
    pub fn from_airfoil(airfoil: &PaneledAirfoil, nw: usize) -> Self {
        let n = airfoil.n;
        let mut st = Self::empty(n, nw);
        for i in 1..=n {
            st.x[i] = airfoil.x[i - 1];
            st.y[i] = airfoil.y[i - 1];
            st.s[i] = airfoil.s[i - 1];
            st.xp[i] = airfoil.xp[i - 1];
            st.yp[i] = airfoil.yp[i - 1];
            st.nx[i] = airfoil.nx[i - 1];
            st.ny[i] = airfoil.ny[i - 1];
            st.apanel[i] = airfoil.apanel[i - 1];
        }
        st.chord = airfoil.chord;
        st.sle = airfoil.sle;
        st.xle = crate::geometry::seval(airfoil.sle, &airfoil.x, &airfoil.xp, &airfoil.s);
        st.yle = crate::geometry::seval(airfoil.sle, &airfoil.y, &airfoil.yp, &airfoil.s);
        st.xte = 0.5 * (st.x[1] + st.x[n]);
        st.yte = 0.5 * (st.y[1] + st.y[n]);
        crate::solver::pointers::tecalc(&mut st);
        st
    }

    /// Install wake node coordinates `n+1..=n+nw` (from XYWAKE, or from a fixture).
    pub fn set_wake_nodes(&mut self, x: &[f64], y: &[f64], s: &[f64]) {
        assert_eq!(x.len(), self.nw);
        for iw in 1..=self.nw {
            self.x[self.n + iw] = x[iw - 1];
            self.y[self.n + iw] = y[iw - 1];
            self.s[self.n + iw] = s[iw - 1];
        }
    }
}
