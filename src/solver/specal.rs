//! SPECAL (xoper.f): converges to the specified alpha — superimposes the alpha = 0, 90
//! vorticity solutions, sets QINV, iterates M(CL) for MATYP ≠ 1, and sets CL/CM/CDP and the Cp
//! distributions. Also the OPER-level bookkeeping that precedes it when a new alpha is specified.

use crate::bl::system::MachClDependence;
use crate::solver::blstate::SolverState;
use crate::solver::clcalc::{compute_cl_cm, compute_cp, set_compressibility};
use crate::solver::ggcalc::{build_inviscid_system, InviscidSystem};
use crate::solver::pointers::set_te_thickness;
use crate::solver::setbl::set_mach_re_from_cl;
use crate::solver::velocity::set_q_inviscid;

/// OPER's `ALFA` command: sets LALFA, ALFA (radians) and QINF = 1, runs SPECAL, then
/// invalidates the wake and the converged flag when alpha or Mach moved by more than 1e-5
/// (in exactly XFOIL's order: SPECAL first, then the tests).
#[doc(alias = "ALFA")]
pub fn alpha_command(state: &mut SolverState, sys: &mut Option<InviscidSystem>, alfa: f64) {
    state.alpha_specified = true;
    state.alpha = alfa;
    state.qinf = 1.0;
    solve_inviscid_at_alpha(state, sys);
    if (state.alpha - state.alpha_wake).abs() > 1.0e-5 {
        state.wake_built = false;
    }
    if (state.alpha - state.alpha_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
    if (state.mach - state.mach_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
}

/// One point of OPER's `ASEQ`: sets ALFA, invalidates the wake/converged flags (ASEQ tests
/// them *before* SPECAL, the ALFA command after — the numbers do not depend on the order), then
/// SPECAL. The caller runs VISCAL with ITMAX + 5, as ASEQ does.
#[doc(alias = "ASEQ")]
pub fn sequence_command(state: &mut SolverState, sys: &mut Option<InviscidSystem>, alfa: f64) {
    state.alpha = alfa;
    if (state.alpha - state.alpha_wake).abs() > 1.0e-5 {
        state.wake_built = false;
    }
    if (state.alpha - state.alpha_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
    if (state.mach - state.mach_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
    solve_inviscid_at_alpha(state, sys);
}

/// SPECAL. `sys` holds the inviscid system (GGCALC's AIJ factors, BIJ, GAMU in `st.qinvu`);
/// it is built here when absent (LGAMU/LQAIJ false).
#[doc(alias = "SPECAL")]
pub fn solve_inviscid_at_alpha(state: &mut SolverState, sys: &mut Option<InviscidSystem>) {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if sys.is_none() {
        *sys = Some(build_inviscid_system(state));
    }

    let cosa = state.alpha.cos();
    let sina = state.alpha.sin();

    // superimpose suitably weighted  alpha = 0, 90  distributions (GAMU(I,·) = QINVU(I,·) on
    // the airfoil nodes; PSIO is not kept)
    for i in 1..=state.n_foil_nodes {
        state.gamma[i] = cosa * state.q_inviscid_basis[1][i] + sina * state.q_inviscid_basis[2][i];
        state.gamma_d_alpha[i] = -sina * state.q_inviscid_basis[1][i] + cosa * state.q_inviscid_basis[2][i];
    }

    set_te_thickness(state);
    set_q_inviscid(state, state.alpha);

    // set initial guess for the Newton variable CLM
    let mut clm = 1.0;

    // set corresponding  M(CLM), Re(CLM)
    let (mut minf_clm, _reinf_clm) = set_mach_re_from_cl(state, clm);
    set_compressibility(state);

    // set corresponding CL(M)
    compute_cl_cm(state);

    // iterate on CLM
    let mut exhausted = true;
    for _itcl in 1..=20 {
        let msq_clm = 2.0 * state.mach * minf_clm;
        let dclm = (state.cl - clm) / (1.0 - state.cl_d_machsqd * msq_clm);

        let clm1 = clm;
        let mut rlx = 1.0;

        // under-relaxation loop to avoid driving M(CL) above 1
        for _irlx in 1..=12 {
            clm = clm1 + rlx * dclm;

            // set new freestream Mach M(CLM)
            let (m, _) = set_mach_re_from_cl(state, clm);
            minf_clm = m;

            // if Mach is OK, go do next Newton iteration
            if state.mach_cl_dependence == MachClDependence::Fixed || state.mach == 0.0 || minf_clm != 0.0 {
                break;
            }
            rlx = 0.5 * rlx;
        }

        // set new CL(M)
        set_compressibility(state);
        compute_cl_cm(state);

        if dclm.abs() <= 1.0e-6 {
            exhausted = false;
            break;
        }
    }
    // 'SPECAL:  Minf convergence failed' if the loop ran out
    // observation only (never read by the solver)
    state.validity.mach_cl_newton_exhausted = exhausted;

    // set final Mach, CL, Cp distributions, and hinge moment
    let (m_cl, re_cl) = set_mach_re_from_cl(state, state.cl);
    state.mach_d_cl = m_cl;
    state.re_d_cl = re_cl;
    set_compressibility(state);
    compute_cl_cm(state);
    (state.cp_inviscid, _) = compute_cp(state.n_foil_nodes, &state.q_inviscid, state.qinf, state.mach);
    let margin = if state.viscous {
        let nt = state.n_foil_nodes + state.n_wake_nodes;
        let (cpv, mv) = compute_cp(nt, &state.q_viscous, state.qinf, state.mach);
        state.cp_viscous = cpv;
        let (cpi, mi) = compute_cp(nt, &state.q_inviscid, state.qinf, state.mach);
        state.cp_inviscid = cpi;
        mv.min(mi)
    } else {
        let (cpi, mi) = compute_cp(state.n_foil_nodes, &state.q_inviscid, state.qinf, state.mach);
        state.cp_inviscid = cpi;
        mi
    };
    state.validity.karman_tsien_margin_pressure = margin;
}

/// SPECCL: converges to the specified inviscid CL by Newton iteration on alpha (20 iterations,
/// |DALFA| ≤ 1e-6), with MINF/REINF set from CLSPEC by MRCL and held fixed. Returns XFOIL's
/// ITAL as the loop leaves it: the iteration that met the tolerance, or 21 when the cap was
/// exhausted ('SPECCL:  CL convergence failed').
#[doc(alias = "SPECCL")]
pub fn solve_inviscid_at_cl(state: &mut SolverState, sys: &mut Option<InviscidSystem>) -> usize {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if sys.is_none() {
        *sys = Some(build_inviscid_system(state));
    }

    // set freestream Mach from specified CL -- Mach will be held fixed
    let (m_cl, re_cl) = set_mach_re_from_cl(state, state.cl_specified);
    state.mach_d_cl = m_cl;
    state.re_d_cl = re_cl;
    set_compressibility(state);

    // current alpha is the initial guess for Newton variable ALFA
    let set_gam = |state: &mut SolverState| {
        let cosa = state.alpha.cos();
        let sina = state.alpha.sin();
        for i in 1..=state.n_foil_nodes {
            state.gamma[i] = cosa * state.q_inviscid_basis[1][i] + sina * state.q_inviscid_basis[2][i];
            state.gamma_d_alpha[i] = -sina * state.q_inviscid_basis[1][i] + cosa * state.q_inviscid_basis[2][i];
        }
    };
    set_gam(state);

    // get corresponding CL, CL_alpha, CL_Mach
    compute_cl_cm(state);

    // Newton loop for alpha to get specified inviscid CL
    let mut ital = 21;
    for it in 1..=20 {
        let dalfa = (state.cl_specified - state.cl) / state.cl_d_alpha;
        let rlx = 1.0;
        state.alpha += rlx * dalfa;

        // set new surface speed distribution
        set_gam(state);

        // set new CL(alpha)
        compute_cl_cm(state);

        if dalfa.abs() <= 1.0e-6 {
            ital = it;
            break;
        }
    }
    // 'SPECCL:  CL convergence failed' if the loop ran out (ital stays 21)

    // set final surface speed and Cp distributions
    set_te_thickness(state);
    set_q_inviscid(state, state.alpha);
    let margin = if state.viscous {
        let nt = state.n_foil_nodes + state.n_wake_nodes;
        let (cpv, mv) = compute_cp(nt, &state.q_viscous, state.qinf, state.mach);
        state.cp_viscous = cpv;
        let (cpi, mi) = compute_cp(nt, &state.q_inviscid, state.qinf, state.mach);
        state.cp_inviscid = cpi;
        mv.min(mi)
    } else {
        let (cpi, mi) = compute_cp(state.n_foil_nodes, &state.q_inviscid, state.qinf, state.mach);
        state.cp_inviscid = cpi;
        mi
    };
    state.validity.karman_tsien_margin_pressure = margin;
    ital
}

/// OPER's `CL` command: LALFA false, ALFA reset to 0 as the Newton initial guess, QINF = 1,
/// SPECCL, then the wake/converged invalidations.
/// Returns SPECCL's exit iteration (see `solve_inviscid_at_cl`).
#[doc(alias = "CL")]
pub fn cl_command(state: &mut SolverState, sys: &mut Option<InviscidSystem>, clspec: f64) -> usize {
    state.cl_specified = clspec;
    state.alpha_specified = false;
    state.alpha = 0.0;
    state.qinf = 1.0;
    let ital = solve_inviscid_at_cl(state, sys);
    if (state.alpha - state.alpha_wake).abs() > 1.0e-5 {
        state.wake_built = false;
    }
    if (state.alpha - state.alpha_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
    if (state.mach - state.mach_converged).abs() > 1.0e-5 {
        state.converged = false;
    }
    ital
}
