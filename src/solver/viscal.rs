//! VISCAL (xoper.f): converges the viscous operating point. Line-for-line on `BlState`:
//! the prologue (XYWAKE, QWCALC, QISET, STFIND/IBLPAN/XICALC/IBLSYS, UICALC, QDCALC) and the
//! Newton loop SETBL → BLSOLV → UPDATE → (MRCL+COMSET | QISET+UICALC) → QVFUE → GAMQV →
//! STMOVE → CLCALC → CDCALC, converging on RMSBL < EPS1.

use crate::bl::blsolv::blsolv;
use crate::solver::blstate::BlState;
use crate::solver::clcalc::{cdcalc, clcalc, comset, cpcalc};
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::pointers::{iblpan, iblsys, stfind, stmove, xicalc};
use crate::solver::qdcalc::qdcalc;
use crate::solver::setbl::{mrcl, setbl};
use crate::solver::update::update;
use crate::solver::velocity::{gamqv, qiset, qvfue, uicalc};
use crate::solver::xywake::{qwcalc, xywake};

/// Convergence tolerance EPS1
pub const EPS1: f64 = 1.0e-4;

/// One VISCAL iteration as XFOIL reports it (after CLCALC/CDCALC, before the convergence test).
#[derive(Debug, Clone, Default)]
pub struct ViscalIter {
    pub iter: usize,
    pub rmsbl: f64,
    pub rmxbl: f64,
    pub vmxbl: char,
    pub imxbl: usize,
    pub ismxbl: usize,
    pub rlx: f64,
    pub alfa: f64,
    pub minf: f64,
    pub reinf: f64,
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    pub cdf: f64,
    pub cdp: f64,
    pub cl_alf: f64,
    pub cl_msq: f64,
    pub ist: usize,
    pub sst: f64,
    pub itran: [usize; 3],
    pub xoctr: [f64; 3],
    /// RMSBL < EPS1 at this iteration
    pub converged: bool,
}

/// VISCAL. `sys` is the inviscid system (AIJ factors and BIJ) QDCALC needs when the source
/// influence matrix does not exist yet; `waklen` is WAKLEN (1.0 in XFOIL). Returns whether the
/// point converged; `st.lvconv/avisc/mvisc` are set as XFOIL sets them.
pub fn viscal(
    st: &mut BlState,
    mut sys: Option<&mut InviscidSystem>,
    niter: usize,
    waklen: f64,
    mut trace: Option<&mut Vec<ViscalIter>>,
) -> bool {
    // calculate wake trajectory from current inviscid solution if necessary
    if !st.lwake {
        xywake(st, waklen);
    }

    // set velocities on wake from airfoil vorticity for alpha=0, 90
    qwcalc(st);

    // set velocities on airfoil and wake for initial alpha
    qiset(st, st.alfa);

    if !st.lipan {
        if st.lblini {
            gamqv(st);
        }
        // locate stagnation point arc length position and panel index
        stfind(st);
        // set  BL position -> panel position  pointers
        iblpan(st);
        // calculate surface arc length array for current stagnation point location
        xicalc(st);
        // set  BL position -> system line  pointers
        iblsys(st);
    }

    // set inviscid BL edge velocity UINV from QINV
    uicalc(st);

    if !st.lblini {
        // set initial Ue from inviscid Ue
        for is in 1..=2 {
            for ibl in 1..=st.nbl[is] {
                st.uedg[is][ibl] = st.uinv[is][ibl];
            }
        }
    }

    if st.lvconv {
        // set correct CL if converged point exists
        qvfue(st);
        let nt = st.n + st.nw;
        if st.lvisc {
            st.cpv = cpcalc(nt, &st.qvis, st.qinf, st.minf);
            st.cpi = cpcalc(nt, &st.qinv, st.qinf, st.minf);
        } else {
            st.cpi = cpcalc(st.n, &st.qinv, st.qinf, st.minf);
        }
        gamqv(st);
        clcalc(st);
        cdcalc(st);
    }

    // set up source influence matrix if it doesn't exist
    let ladij = sys.as_ref().map(|s| s.ladij).unwrap_or(true);
    if !st.lwdij || !ladij {
        let s = sys
            .as_mut()
            .expect("VISCAL: the source influence matrix does not exist and no inviscid system was given");
        qdcalc(st, s);
    }

    // Newton iteration for entire BL solution
    let mut converged = false;
    for iter in 1..=niter {
        // fill Newton system for BL variables
        let r = setbl(st);
        // solve Newton system with custom solver
        let sol = blsolv(r.sys);
        // update BL variables
        let u = update(st, &sol.vdel, st.minf_cl);

        if st.lalfa {
            // set new freestream Mach, Re from new CL
            let (m_cl, re_cl) = mrcl(st, st.cl);
            st.minf_cl = m_cl;
            st.reinf_cl = re_cl;
            comset(st);
        } else {
            // set new inviscid speeds QINV and UINV for new alpha
            qiset(st, st.alfa);
            uicalc(st);
        }

        // calculate edge velocities QVIS(.) from UEDG(..)
        qvfue(st);
        // set GAM distribution from QVIS
        gamqv(st);
        // relocate stagnation point
        stmove(st);
        // set updated CL,CD
        clcalc(st);
        cdcalc(st);

        // display changes and test for convergence
        let conv = u.rmsbl < EPS1;
        if let Some(t) = trace.as_mut() {
            t.push(ViscalIter {
                iter,
                rmsbl: u.rmsbl,
                rmxbl: u.rmxbl,
                vmxbl: u.vmxbl,
                imxbl: u.imxbl,
                ismxbl: u.ismxbl,
                rlx: u.rlx,
                alfa: st.alfa,
                minf: st.minf,
                reinf: st.reinf,
                cl: st.cl,
                cm: st.cm,
                cd: st.cd,
                cdf: st.cdf,
                cdp: st.cdp,
                cl_alf: st.cl_alf,
                cl_msq: st.cl_msq,
                ist: st.ist,
                sst: st.sst,
                itran: st.itran,
                xoctr: st.xoctr,
                converged: conv,
            });
        }
        if conv {
            st.lvconv = true;
            st.avisc = st.alfa;
            st.mvisc = st.minf;
            converged = true;
            break;
        }
    }
    // 'VISCAL:  Convergence failed' if the loop ran out

    let nt = st.n + st.nw;
    st.cpi = cpcalc(nt, &st.qinv, st.qinf, st.minf);
    st.cpv = cpcalc(nt, &st.qvis, st.qinf, st.minf);
    converged
}
