//! SPECAL (xoper.f): converges to the specified alpha — superimposes the alpha = 0, 90
//! vorticity solutions, sets QINV, iterates M(CL) for MATYP ≠ 1, and sets CL/CM/CDP and the Cp
//! distributions. Also the OPER-level bookkeeping that precedes it when a new alpha is specified.

use crate::solver::blstate::BlState;
use crate::solver::clcalc::{clcalc, comset, cpcalc};
use crate::solver::ggcalc::{ggcalc, InviscidSystem};
use crate::solver::pointers::tecalc;
use crate::solver::setbl::set_mach_re_from_cl;
use crate::solver::velocity::qiset;

/// OPER's `ALFA` command: sets LALFA, ALFA (radians) and QINF = 1, runs SPECAL, then
/// invalidates the wake and the converged flag when alpha or Mach moved by more than 1e-5
/// (in exactly XFOIL's order: SPECAL first, then the tests).
pub fn alfa_command(st: &mut BlState, sys: &mut Option<InviscidSystem>, alfa: f64) {
    st.lalfa = true;
    st.alfa = alfa;
    st.qinf = 1.0;
    specal(st, sys);
    if (st.alfa - st.awake).abs() > 1.0e-5 {
        st.lwake = false;
    }
    if (st.alfa - st.avisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
    if (st.minf - st.mvisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
}

/// One point of OPER's `ASEQ`: sets ALFA, invalidates the wake/converged flags (ASEQ tests
/// them *before* SPECAL, the ALFA command after — the numbers do not depend on the order), then
/// SPECAL. The caller runs VISCAL with ITMAX + 5, as ASEQ does.
pub fn aseq_point(st: &mut BlState, sys: &mut Option<InviscidSystem>, alfa: f64) {
    st.alfa = alfa;
    if (st.alfa - st.awake).abs() > 1.0e-5 {
        st.lwake = false;
    }
    if (st.alfa - st.avisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
    if (st.minf - st.mvisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
    specal(st, sys);
}

/// SPECAL. `sys` holds the inviscid system (GGCALC's AIJ factors, BIJ, GAMU in `st.qinvu`);
/// it is built here when absent (LGAMU/LQAIJ false).
pub fn specal(st: &mut BlState, sys: &mut Option<InviscidSystem>) {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if sys.is_none() {
        *sys = Some(ggcalc(st));
    }

    let cosa = st.alfa.cos();
    let sina = st.alfa.sin();

    // superimpose suitably weighted  alpha = 0, 90  distributions (GAMU(I,·) = QINVU(I,·) on
    // the airfoil nodes; PSIO is not kept)
    for i in 1..=st.n {
        st.gam[i] = cosa * st.qinvu[1][i] + sina * st.qinvu[2][i];
        st.gam_a[i] = -sina * st.qinvu[1][i] + cosa * st.qinvu[2][i];
    }

    tecalc(st);
    qiset(st, st.alfa);

    // set initial guess for the Newton variable CLM
    let mut clm = 1.0;

    // set corresponding  M(CLM), Re(CLM)
    let (mut minf_clm, _reinf_clm) = set_mach_re_from_cl(st, clm);
    comset(st);

    // set corresponding CL(M)
    clcalc(st);

    // iterate on CLM
    for _itcl in 1..=20 {
        let msq_clm = 2.0 * st.minf * minf_clm;
        let dclm = (st.cl - clm) / (1.0 - st.cl_msq * msq_clm);

        let clm1 = clm;
        let mut rlx = 1.0;

        // under-relaxation loop to avoid driving M(CL) above 1
        for _irlx in 1..=12 {
            clm = clm1 + rlx * dclm;

            // set new freestream Mach M(CLM)
            let (m, _) = set_mach_re_from_cl(st, clm);
            minf_clm = m;

            // if Mach is OK, go do next Newton iteration
            if st.matyp == 1 || st.minf == 0.0 || minf_clm != 0.0 {
                break;
            }
            rlx = 0.5 * rlx;
        }

        // set new CL(M)
        comset(st);
        clcalc(st);

        if dclm.abs() <= 1.0e-6 {
            break;
        }
    }
    // 'SPECAL:  Minf convergence failed' if the loop ran out

    // set final Mach, CL, Cp distributions, and hinge moment
    let (m_cl, re_cl) = set_mach_re_from_cl(st, st.cl);
    st.minf_cl = m_cl;
    st.reinf_cl = re_cl;
    comset(st);
    clcalc(st);
    st.cpi = cpcalc(st.n, &st.qinv, st.qinf, st.minf);
    if st.lvisc {
        let nt = st.n + st.nw;
        st.cpv = cpcalc(nt, &st.qvis, st.qinf, st.minf);
        st.cpi = cpcalc(nt, &st.qinv, st.qinf, st.minf);
    } else {
        st.cpi = cpcalc(st.n, &st.qinv, st.qinf, st.minf);
    }
}

/// SPECCL: converges to the specified inviscid CL by Newton iteration on alpha (20 iterations,
/// |DALFA| ≤ 1e-6), with MINF/REINF set from CLSPEC by MRCL and held fixed.
pub fn speccl(st: &mut BlState, sys: &mut Option<InviscidSystem>) {
    // calculate surface vorticity distributions for alpha = 0, 90 degrees
    if sys.is_none() {
        *sys = Some(ggcalc(st));
    }

    // set freestream Mach from specified CL -- Mach will be held fixed
    let (m_cl, re_cl) = set_mach_re_from_cl(st, st.clspec);
    st.minf_cl = m_cl;
    st.reinf_cl = re_cl;
    comset(st);

    // current alpha is the initial guess for Newton variable ALFA
    let set_gam = |st: &mut BlState| {
        let cosa = st.alfa.cos();
        let sina = st.alfa.sin();
        for i in 1..=st.n {
            st.gam[i] = cosa * st.qinvu[1][i] + sina * st.qinvu[2][i];
            st.gam_a[i] = -sina * st.qinvu[1][i] + cosa * st.qinvu[2][i];
        }
    };
    set_gam(st);

    // get corresponding CL, CL_alpha, CL_Mach
    clcalc(st);

    // Newton loop for alpha to get specified inviscid CL
    for _ital in 1..=20 {
        let dalfa = (st.clspec - st.cl) / st.cl_alf;
        let rlx = 1.0;
        st.alfa += rlx * dalfa;

        // set new surface speed distribution
        set_gam(st);

        // set new CL(alpha)
        clcalc(st);

        if dalfa.abs() <= 1.0e-6 {
            break;
        }
    }
    // 'SPECCL:  CL convergence failed' if the loop ran out

    // set final surface speed and Cp distributions
    tecalc(st);
    qiset(st, st.alfa);
    if st.lvisc {
        let nt = st.n + st.nw;
        st.cpv = cpcalc(nt, &st.qvis, st.qinf, st.minf);
        st.cpi = cpcalc(nt, &st.qinv, st.qinf, st.minf);
    } else {
        st.cpi = cpcalc(st.n, &st.qinv, st.qinf, st.minf);
    }
}

/// OPER's `CL` command: LALFA false, ALFA reset to 0 as the Newton initial guess, QINF = 1,
/// SPECCL, then the wake/converged invalidations.
pub fn cl_command(st: &mut BlState, sys: &mut Option<InviscidSystem>, clspec: f64) {
    st.clspec = clspec;
    st.lalfa = false;
    st.alfa = 0.0;
    st.qinf = 1.0;
    speccl(st, sys);
    if (st.alfa - st.awake).abs() > 1.0e-5 {
        st.lwake = false;
    }
    if (st.alfa - st.avisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
    if (st.minf - st.mvisc).abs() > 1.0e-5 {
        st.lvconv = false;
    }
}
