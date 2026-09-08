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
    state: &mut SolverState,
    mut sys: Option<&mut InviscidSystem>,
    niter: usize,
    waklen: f64,
    mut trace: Option<&mut Vec<IterationRecord>>,
) -> bool {
    // calculate wake trajectory from current inviscid solution if necessary
    if !state.wake_built {
        build_wake(state, waklen);
    }

    // set velocities on wake from airfoil vorticity for alpha=0, 90
    set_wake_q_basis(state);

    // set velocities on airfoil and wake for initial alpha
    set_q_inviscid(state, state.alpha);

    if !state.pointers_built {
        if state.bl_initialised {
            set_gamma_from_q_viscous(state);
        }
        // locate stagnation point arc length position and panel index
        find_stagnation(state);
        // set  BL position -> panel position  pointers
        map_stations_to_nodes(state);
        // calculate surface arc length array for current stagnation point location
        set_station_xi(state);
        // set  BL position -> system line  pointers
        map_stations_to_rows(state);
    }

    // set inviscid BL edge velocity UINV from QINV
    set_ue_inviscid(state);

    if !state.bl_initialised {
        // set initial Ue from inviscid Ue
        for side in 1..=2 {
            for i_station in 1..=state.n_stations[side] {
                state.ue[side][i_station] = state.ue_inviscid[side][i_station];
            }
        }
    }

    if state.converged {
        // set correct CL if converged point exists
        set_q_viscous_from_ue(state);
        let nt = state.n_foil_nodes + state.n_wake_nodes;
        if state.viscous {
            state.cp_viscous = compute_cp(nt, &state.q_viscous, state.qinf, state.mach);
            state.cp_inviscid = compute_cp(nt, &state.q_inviscid, state.qinf, state.mach);
        } else {
            state.cp_inviscid = compute_cp(state.n_foil_nodes, &state.q_inviscid, state.qinf, state.mach);
        }
        set_gamma_from_q_viscous(state);
        compute_cl_cm(state);
        compute_cd(state);
    }

    // set up source influence matrix if it doesn't exist
    let ladij = sys.as_ref().map(|s| s.dij_foil_built).unwrap_or(true);
    if !state.dij_wake_built || !ladij {
        let s = sys
            .as_mut()
            .expect("VISCAL: the source influence matrix does not exist and no inviscid system was given");
        build_dij(state, s);
    }

    // Newton iteration for entire BL solution
    let mut converged = false;
    for iter in 1..=niter {
        // fill Newton system for BL variables
        let r = assemble_newton_system(state);
        // solve Newton system with custom solver
        let sol = solve_newton_system(r.newton);
        // update BL variables
        let u = apply_newton_update(state, &sol.deltas, state.mach_d_cl);

        if state.alpha_specified {
            // set new freestream Mach, Re from new CL
            let (m_cl, re_cl) = set_mach_re_from_cl(state, state.cl);
            state.mach_d_cl = m_cl;
            state.re_d_cl = re_cl;
            set_compressibility(state);
        } else {
            // set new inviscid speeds QINV and UINV for new alpha
            set_q_inviscid(state, state.alpha);
            set_ue_inviscid(state);
        }

        // calculate edge velocities QVIS(.) from UEDG(..)
        set_q_viscous_from_ue(state);
        // set GAM distribution from QVIS
        set_gamma_from_q_viscous(state);
        // relocate stagnation point
        move_stagnation(state);
        // set updated CL,CD
        compute_cl_cm(state);
        compute_cd(state);

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
                alpha: state.alpha,
                mach: state.mach,
                re: state.re,
                cl: state.cl,
                cm: state.cm,
                cd: state.cd,
                cd_friction: state.cd_friction,
                cd_pressure: state.cd_pressure,
                cl_d_alpha: state.cl_d_alpha,
                cl_d_machsqd: state.cl_d_machsqd,
                i_stagnation_node: state.i_stagnation_node,
                s_stagnation: state.s_stagnation,
                i_transition_station: state.i_transition_station,
                x_transition: state.x_transition,
                converged: conv,
            });
        }
        if conv {
            state.converged = true;
            state.alpha_converged = state.alpha;
            state.mach_converged = state.mach;
            converged = true;
            break;
        }
    }
    // 'VISCAL:  Convergence failed' if the loop ran out

    let nt = state.n_foil_nodes + state.n_wake_nodes;
    state.cp_inviscid = compute_cp(nt, &state.q_inviscid, state.qinf, state.mach);
    state.cp_viscous = compute_cp(nt, &state.q_viscous, state.qinf, state.mach);
    converged
}
