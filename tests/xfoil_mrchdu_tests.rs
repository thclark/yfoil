//! S6 gate: MRCHDU (mixed-mode march along the Ue-Hk characteristic) against the reference.
//!
//! Two views of the same run: for SETBL calls 1..3, the BL state before MRCHDU
//! (`mrchdu_input_<k>.dat`) is loaded, marched, and compared with the state after
//! (`mrchdu_output_<k>.dat`); and on call 1 the per-iteration Newton trace
//! (`xfoil_mrchdu_trace.dat`) is compared station by station, iteration by iteration.

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::{parse_bl_dump, parse_mrchdu_trace};

use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::bl::mrchdu::{march_prescribed_dstar, MrchduTrace};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn check_state_after(k: usize) {
    let (mut st, params, _) = fixtures::state_before_setbl_march(k);
    let acrit = st.ncrit;
    let o = parse_bl_dump(&fixture_path(&format!("mrchdu_output_{k}.dat")));
    march_prescribed_dstar(&mut st, &params, acrit, None);
    assert_eq!(
        st.i_transition_station[1..],
        [o.int("ITRAN1"), o.int("ITRAN2")],
        "call {k}: ITRAN"
    );
    assert_eq!(
        st.transition_forced[1..],
        [o.logical("TFORCE1"), o.logical("TFORCE2")],
        "call {k}: TFORCE"
    );
    for is in 1..=2 {
        assert_within(
            st.xi_transition[is],
            o.real(&format!("XSSITR{is}")),
            TOL_SOLVER,
            1.0,
            &format!("call {k}: XSSITR({is})"),
        );
    }
    let names = [
        "XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS", "TAU", "DIS", "CTQ", "DELT", "TSTR",
    ];
    let mut worst = (0.0_f64, "", 0, 0);
    for is in 1..=2 {
        let nbl = o.nbl(is);
        assert_eq!(nbl, st.n_stations[is], "call {k}: NBL({is})");
        for ibl in 2..=nbl {
            let ours = [
                st.xi[is][ibl],
                st.ue[is][ibl],
                st.theta[is][ibl],
                st.dstar[is][ibl],
                st.sqrtctau[is][ibl],
                st.mass_defect[is][ibl],
                st.tau[is][ibl],
                st.dissipation[is][ibl],
                st.sqrtctaueq[is][ibl],
                st.delta[is][ibl],
                st.thetastar[is][ibl],
            ];
            let theirs = |j: usize, m: usize| if m < 6 { o.bl[is][j][m] } else { o.blx[is][j][m - 6] };
            for (m, name) in names.iter().enumerate() {
                let (a, b) = (ours[m], theirs(ibl, m));
                // scale: the largest value of that variable on the side
                let scale = (2..=nbl).map(|j| theirs(j, m).abs()).fold(0.0_f64, f64::max);
                let e = (a - b).abs() / a.abs().max(b.abs()).max(scale);
                if e > worst.0 {
                    worst = (e, name, is, ibl);
                }
                assert_within(a, b, TOL_SOLVER, scale, &format!("call {k}: {name}({ibl},{is})"));
            }
        }
    }
    println!(
        "mrchdu call {k}: state after march within {TOL_SOLVER:.0e} (worst {:.2e} {} at ({},{})); ITRAN={:?}",
        worst.0,
        worst.1,
        worst.3,
        worst.2,
        &st.i_transition_station[1..]
    );
}

#[test]
fn test_mrchdu_reproduces_xfoil_state_call_1() {
    check_state_after(1);
}

#[test]
fn test_mrchdu_reproduces_xfoil_state_call_2() {
    check_state_after(2);
}

#[test]
fn test_mrchdu_reproduces_xfoil_state_call_3() {
    check_state_after(3);
}

/// Compare two vectors with the row-scaled metric, appending mismatches.
fn chk(mism: &mut Vec<String>, ctx: &str, a: &[f64], b: &[f64], what: &str, floor: f64) {
    let scale = b
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(a.iter().fold(0.0_f64, |m, v| m.max(v.abs())))
        .max(floor);
    for i in 0..a.len().min(b.len()) {
        let e = (a[i] - b[i]).abs() / scale;
        if e > TOL_SOLVER {
            mism.push(format!(
                "{ctx} {what}[{i}] yfoil={:+.10e} xfoil={:+.10e} scaled-err={e:.2e}",
                a[i], b[i]
            ));
        }
    }
}

#[test]
fn test_mrchdu_newton_trace_matches_xfoil_iteration_by_iteration() {
    let (mut st, params, _) = fixtures::state_before_setbl_march(1);
    let acrit = st.ncrit;
    let xf = parse_mrchdu_trace(&fixture_path("xfoil_mrchdu_trace.dat"));
    let mut tr = MrchduTrace::default();
    march_prescribed_dstar(&mut st, &params, acrit, Some(&mut tr));

    let mut mism = Vec::new();
    let n = xf.iters.len().min(tr.iters.len());
    if xf.iters.len() != tr.iters.len() {
        mism.push(format!(
            "iteration count: yfoil {} vs xfoil {} station-iterations",
            tr.iters.len(),
            xf.iters.len()
        ));
    }
    for i in 0..n {
        let (y, x) = (&tr.iters[i], &xf.iters[i]);
        let ctx = format!("IBL={:3} IS={} ITBL={:2}", x.ibl, x.is, x.itbl);
        if (y.ibl, y.is, y.itbl) != (x.ibl, x.is, x.itbl) {
            mism.push(format!(
                "{ctx}: sequence diverged (yfoil at IBL={} IS={} ITBL={})",
                y.ibl, y.is, y.itbl
            ));
            break;
        }
        chk(&mut mism, &ctx, &y.ampl, &x.ampl, "AMPL(A1,A2,XT,ACRIT)", 1.0);
        if y.tran != x.tran || y.itran != x.itran {
            mism.push(format!(
                "{ctx} TRAN/ITRAN yfoil=({},{}) xfoil=({},{})",
                y.tran, y.itran, x.tran, x.itran
            ));
        }
        for (name, a, b) in [
            ("PRIMARY", &y.primary, &x.primary),
            ("KINEMATIC", &y.kinematic, &x.kinematic),
            ("CLOSURE", &y.closure, &x.closure),
        ] {
            for j in 0..5 {
                chk(&mut mism, &ctx, &[a[j]], &[b[j]], &format!("{name}[{j}]"), 0.0);
            }
        }
        chk(
            &mut mism,
            &ctx,
            &[y.ueref, y.hkref],
            &[x.ueref, x.hkref],
            "REF(UEREF,HKREF)",
            0.0,
        );
        match (y.sens, x.sens) {
            (Some(a), Some(b)) => chk(&mut mism, &ctx, &a, &b, "SENS(SENNEW,SENS)", 0.0),
            (None, None) => {}
            _ => mism.push(format!("{ctx} SENS present in one trace only")),
        }
        // residuals are O(1) nondimensional equations: compare on the unit scale
        chk(&mut mism, &ctx, &y.residual, &x.residual, "RESIDUAL", 1.0);
        for k in 0..4 {
            chk(&mut mism, &ctx, &y.vs2[k], &x.vs2[k], &format!("VS2_{}", k + 1), 0.0);
        }
        // solution scaled by the variable it updates
        let sc = [x.primary[4].abs().max(1e-3), x.primary[2], x.primary[3], x.primary[1]];
        for j in 0..4 {
            chk(
                &mut mism,
                &ctx,
                &[y.solution[j] / sc[j]],
                &[x.solution[j] / sc[j]],
                &format!("SOLUTION[{j}]"),
                1.0,
            );
        }
        chk(&mut mism, &ctx, &[y.dmax, y.rlx], &[x.dmax, x.rlx], "RELAX", 1.0);
        for j in 0..5 {
            chk(
                &mut mism,
                &ctx,
                &[y.updated[j]],
                &[x.updated[j]],
                &format!("UPDATED[{j}]"),
                0.0,
            );
        }
        if y.converged != x.converged {
            mism.push(format!("{ctx} CONVERGED yfoil={} xfoil={}", y.converged, x.converged));
        }
        if mism.len() > 40 {
            break;
        }
    }
    if tr.failed != xf.failed {
        mism.push(format!("FAILED stations: yfoil {:?} xfoil {:?}", tr.failed, xf.failed));
    }
    println!(
        "mrchdu trace: {} station-iterations compared, {} mismatching values",
        n,
        mism.len()
    );
    for m in &mism {
        println!("  {m}");
    }
    assert!(
        mism.is_empty(),
        "{} mismatches against the reference Newton trace (first: {})",
        mism.len(),
        mism[0]
    );
}
