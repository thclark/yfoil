//! VISCAL (xoper.f): converges the viscous operating point. Line-for-line on `BlState`:
//! the prologue (XYWAKE, QWCALC, QISET, STFIND/IBLPAN/XICALC/IBLSYS, UICALC, QDCALC) and the
//! Newton loop SETBL → BLSOLV → UPDATE → (MRCL+COMSET | QISET+UICALC) → QVFUE → GAMQV →
//! STMOVE → CLCALC → CDCALC, converging on RMSBL < EPS1.

use crate::bl::blsolv::solve_newton_system;
use crate::solver::blstate::BlState;
use crate::solver::clcalc::{cdcalc, clcalc, comset, cpcalc};
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::pointers::{iblpan, iblsys, stfind, stmove, xicalc};
use crate::solver::qdcalc::qdcalc;
use crate::solver::setbl::{assemble_newton_system, set_mach_re_from_cl};
use crate::solver::update::apply_newton_update;
use crate::solver::velocity::{gamqv, qiset, qvfue, uicalc};
use crate::solver::xywake::{qwcalc, xywake};

/// Convergence tolerance EPS1
pub const CONVERGENCE_TOLERANCE: f64 = 1.0e-4;

/// One VISCAL iteration as XFOIL reports it (after CLCALC/CDCALC, before the convergence test).
#[derive(Debug, Clone, Default)]
pub struct IterationRecord {
    pub iteration: usize,
    pub residual: f64,
    pub residual_max: f64,
    pub residual_max_variable: char,
    pub i_residual_max_station: usize,
    pub residual_max_side: usize,
    pub relaxation: f64,
    pub alpha: f64,
    pub mach: f64,
    pub re: f64,
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    pub cd_friction: f64,
    pub cd_pressure: f64,
    pub cl_d_alpha: f64,
    pub cl_d_machsqd: f64,
    pub i_stagnation_node: usize,
    pub s_stagnation: f64,
    pub i_transition_station: [usize; 3],
    pub x_transition: [f64; 3],
    /// RMSBL < EPS1 at this iteration
    pub converged: bool,
}

/// VISCAL. `sys` is the inviscid system (AIJ factors and BIJ) QDCALC needs when the source
/// influence matrix does not exist yet; `waklen` is WAKLEN (1.0 in XFOIL). Returns whether the
/// point converged; `st.lvconv/avisc/mvisc` are set as XFOIL sets them.
pub fn solve_viscous(
    st: &mut BlState,
    mut sys: Option<&mut InviscidSystem>,
    niter: usize,
    waklen: f64,
    mut trace: Option<&mut Vec<IterationRecord>>,
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
        let r = assemble_newton_system(st);
        // solve Newton system with custom solver
        let sol = solve_newton_system(r.newton);
        // update BL variables
        let u = apply_newton_update(st, &sol.deltas, st.minf_cl);

        if st.lalfa {
            // set new freestream Mach, Re from new CL
            let (m_cl, re_cl) = set_mach_re_from_cl(st, st.cl);
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
        let conv = u.residual < CONVERGENCE_TOLERANCE;
        if let Some(t) = trace.as_mut() {
            t.push(IterationRecord {
                iteration: iter,
                residual: u.residual,
                residual_max: u.residual_max,
                residual_max_variable: u.residual_max_variable,
                i_residual_max_station: u.i_residual_max_station,
                residual_max_side: u.residual_max_side,
                relaxation: u.relaxation,
                alpha: st.alfa,
                mach: st.minf,
                re: st.reinf,
                cl: st.cl,
                cm: st.cm,
                cd: st.cd,
                cd_friction: st.cdf,
                cd_pressure: st.cdp,
                cl_d_alpha: st.cl_alf,
                cl_d_machsqd: st.cl_msq,
                i_stagnation_node: st.ist,
                s_stagnation: st.sst,
                i_transition_station: st.itran,
                x_transition: st.xoctr,
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
