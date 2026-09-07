//! VISCAL (xoper.f): converges the viscous operating point. Line-for-line on `SolverState`:
//! the prologue (XYWAKE, QWCALC, QISET, STFIND/IBLPAN/XICALC/IBLSYS, UICALC, QDCALC) and the
//! Newton loop SETBL → BLSOLV → UPDATE → (MRCL+COMSET | QISET+UICALC) → QVFUE → GAMQV →
//! STMOVE → CLCALC → CDCALC, converging on RMSBL < EPS1.

use crate::bl::blsolv::solve_newton_system;
use crate::solver::blstate::SolverState;
use crate::solver::clcalc::{compute_cd, compute_cl_cm, compute_cp, set_compressibility};
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::pointers::{
    find_stagnation, map_stations_to_nodes, map_stations_to_rows, move_stagnation, set_station_xi,
};
use crate::solver::qdcalc::build_dij;
use crate::solver::setbl::{assemble_newton_system, set_mach_re_from_cl};
use crate::solver::update::{apply_newton_update, ResidualMaxVariable};
use crate::solver::velocity::{set_gamma_from_q_viscous, set_q_inviscid, set_q_viscous_from_ue, set_ue_inviscid};
use crate::solver::xywake::{build_wake, set_wake_q_basis};

/// Convergence tolerance EPS1
pub const CONVERGENCE_TOLERANCE: f64 = 1.0e-4;

/// One VISCAL iteration as XFOIL reports it (after CLCALC/CDCALC, before the convergence test).
#[derive(Debug, Clone, Default)]
pub struct IterationRecord {
    pub iteration: usize,
    pub residual: f64,
    pub residual_max: f64,
    pub residual_max_variable: ResidualMaxVariable,
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
#[doc(alias = "VISCAL")]
pub fn solve_viscous(
    st: &mut SolverState,
    mut sys: Option<&mut InviscidSystem>,
    niter: usize,
    waklen: f64,
    mut trace: Option<&mut Vec<IterationRecord>>,
) -> bool {
    // calculate wake trajectory from current inviscid solution if necessary
    if !st.wake_built {
        build_wake(st, waklen);
    }

    // set velocities on wake from airfoil vorticity for alpha=0, 90
    set_wake_q_basis(st);

    // set velocities on airfoil and wake for initial alpha
    set_q_inviscid(st, st.alpha);

    if !st.pointers_built {
        if st.bl_initialised {
            set_gamma_from_q_viscous(st);
        }
        // locate stagnation point arc length position and panel index
        find_stagnation(st);
        // set  BL position -> panel position  pointers
        map_stations_to_nodes(st);
        // calculate surface arc length array for current stagnation point location
        set_station_xi(st);
        // set  BL position -> system line  pointers
        map_stations_to_rows(st);
    }

    // set inviscid BL edge velocity UINV from QINV
    set_ue_inviscid(st);

    if !st.bl_initialised {
        // set initial Ue from inviscid Ue
        for is in 1..=2 {
            for ibl in 1..=st.n_stations[is] {
                st.ue[is][ibl] = st.ue_inviscid[is][ibl];
            }
        }
    }

    if st.converged {
        // set correct CL if converged point exists
        set_q_viscous_from_ue(st);
        let nt = st.n_foil_nodes + st.n_wake_nodes;
        if st.viscous {
            st.cp_viscous = compute_cp(nt, &st.q_viscous, st.qinf, st.mach);
            st.cp_inviscid = compute_cp(nt, &st.q_inviscid, st.qinf, st.mach);
        } else {
            st.cp_inviscid = compute_cp(st.n_foil_nodes, &st.q_inviscid, st.qinf, st.mach);
        }
        set_gamma_from_q_viscous(st);
        compute_cl_cm(st);
        compute_cd(st);
    }

    // set up source influence matrix if it doesn't exist
    let ladij = sys.as_ref().map(|s| s.dij_foil_built).unwrap_or(true);
    if !st.dij_wake_built || !ladij {
        let s = sys
            .as_mut()
            .expect("VISCAL: the source influence matrix does not exist and no inviscid system was given");
        build_dij(st, s);
    }

    // Newton iteration for entire BL solution
    let mut converged = false;
    for iter in 1..=niter {
        // fill Newton system for BL variables
        let r = assemble_newton_system(st);
        // solve Newton system with custom solver
        let sol = solve_newton_system(r.newton);
        // update BL variables
        let u = apply_newton_update(st, &sol.deltas, st.mach_d_cl);

        if st.alpha_specified {
            // set new freestream Mach, Re from new CL
            let (m_cl, re_cl) = set_mach_re_from_cl(st, st.cl);
            st.mach_d_cl = m_cl;
            st.re_d_cl = re_cl;
            set_compressibility(st);
        } else {
            // set new inviscid speeds QINV and UINV for new alpha
            set_q_inviscid(st, st.alpha);
            set_ue_inviscid(st);
        }

        // calculate edge velocities QVIS(.) from UEDG(..)
        set_q_viscous_from_ue(st);
        // set GAM distribution from QVIS
        set_gamma_from_q_viscous(st);
        // relocate stagnation point
        move_stagnation(st);
        // set updated CL,CD
        compute_cl_cm(st);
        compute_cd(st);

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
                alpha: st.alpha,
                mach: st.mach,
                re: st.re,
                cl: st.cl,
                cm: st.cm,
                cd: st.cd,
                cd_friction: st.cd_friction,
                cd_pressure: st.cd_pressure,
                cl_d_alpha: st.cl_d_alpha,
                cl_d_machsqd: st.cl_d_machsqd,
                i_stagnation_node: st.i_stagnation_node,
                s_stagnation: st.s_stagnation,
                i_transition_station: st.i_transition_station,
                x_transition: st.x_transition,
                converged: conv,
            });
        }
        if conv {
            st.converged = true;
            st.alpha_converged = st.alpha;
            st.mach_converged = st.mach;
            converged = true;
            break;
        }
    }
    // 'VISCAL:  Convergence failed' if the loop ran out

    let nt = st.n_foil_nodes + st.n_wake_nodes;
    st.cp_inviscid = compute_cp(nt, &st.q_inviscid, st.qinf, st.mach);
    st.cp_viscous = compute_cp(nt, &st.q_viscous, st.qinf, st.mach);
    converged
}
