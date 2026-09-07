//! S5 gate: MRCHUE (first-pass direct march with Ue prescribed) against the reference.
//!
//! Two views of the same run: the final BL state after MRCHUE (`mrchdu_input_1.dat`, which the
//! instrumentation writes just before the first MRCHDU) and the per-iteration Newton trace
//! (`xfoil_newton_trace.dat`), which lets the march be forward-stepped station by station so
//! the *first* divergent quantity is reported, per CLAUDE.md Rule 3.

mod fixtures;
mod utilities;

use fixtures::mrchue_fixtures::{parse_bl_state, parse_newton_trace};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::bl::mrchue::{march_direct, MrchueTrace};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// State exactly as VISCAL hands it to SETBL/MRCHUE on the first call: pointer layer and wake
/// from the reference, UEDG = UINV, transition free.
use fixtures::state_before_mrchue;

#[test]
fn test_mrchue_reproduces_xfoil_state_after_first_march() {
    let (mut st, params, acrit) = state_before_mrchue();
    let d = parse_bl_state(&fixture_path("mrchdu_input_1.dat"));
    march_direct(&mut st, &params, acrit, None);
    assert_eq!(st.i_transition_station[1..], d.itran[1..], "ITRAN");
    let names = ["XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS"];
    let mut worst = (0.0_f64, "", 0, 0);
    for is in 1..=2 {
        for ibl in 2..=d.nbl[is] {
            let ours = [
                st.xi[is][ibl],
                st.ue[is][ibl],
                st.theta[is][ibl],
                st.dstar[is][ibl],
                st.sqrtctau[is][ibl],
                st.mass_defect[is][ibl],
            ];
            for (k, name) in names.iter().enumerate() {
                let (a, b) = (ours[k], d.bl[is][ibl][k]);
                // scale: the largest value of that variable on the side
                let scale = (2..=d.nbl[is]).map(|j| d.bl[is][j][k].abs()).fold(0.0_f64, f64::max);
                let e = (a - b).abs() / a.abs().max(b.abs()).max(scale);
                if e > worst.0 {
                    worst = (e, name, is, ibl);
                }
                assert_within(a, b, TOL_SOLVER, scale, &format!("{name}({ibl},{is})"));
            }
        }
    }
    println!(
        "mrchue: state after first march within {TOL_SOLVER:.0e} (worst {:.2e} {} at ({},{})); ITRAN={:?}",
        worst.0,
        worst.1,
        worst.3,
        worst.2,
        &st.i_transition_station[1..]
    );
}

/// Compare two vectors with the row-scaled metric, appending mismatches.
fn chk(mism: &mut Vec<String>, ctx: &str, a: &[f64], b: &[f64], what: &str, floor: f64) {
    let scale = b
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(floor)
        .max(f64::MIN_POSITIVE);
    for (i, (p, q)) in a.iter().zip(b).enumerate() {
        let e = (p - q).abs() / p.abs().max(q.abs()).max(scale);
        if e > TOL_SOLVER {
            mism.push(format!(
                "{ctx} {what}[{i}] yfoil={p:+.10e} xfoil={q:+.10e} scaled-err={e:.2e}"
            ));
        }
    }
}

/// Forward-step the march against XFOIL's per-iteration trace. Every station-iteration is
/// compared field by field; all mismatches (not just the first) are listed so the pattern
/// is visible, and the test fails if there are any.
#[test]
fn test_mrchue_newton_trace_matches_xfoil_iteration_by_iteration() {
    let (mut st, params, acrit) = state_before_mrchue();
    let mut trace = MrchueTrace::default();
    march_direct(&mut st, &params, acrit, Some(&mut trace));
    let xf = parse_newton_trace(&fixture_path("xfoil_newton_trace.dat"));
    assert!(!xf.is_empty());
    let mut mism: Vec<String> = Vec::new();
    let mut compared = 0;
    for (k, x) in xf.iter().enumerate() {
        let Some(y) = trace.iters.get(k) else {
            mism.push(format!(
                "YFoil trace ended after {k} iterations; XFOIL has {} (next: IBL={} IS={} ITBL={})",
                xf.len(),
                x.ibl,
                x.is,
                x.itbl
            ));
            break;
        };
        if (y.is, y.ibl, y.itbl) != (x.is, x.ibl, x.itbl) {
            mism.push(format!(
                "sequence diverged at record {k}: yfoil=(IS={},IBL={},ITBL={}) xfoil=(IS={},IBL={},ITBL={})",
                y.is, y.ibl, y.itbl, x.is, x.ibl, x.itbl
            ));
            break;
        }
        let ctx = format!("IBL={:3} IS={} ITBL={:2}", x.ibl, x.is, x.itbl);
        if (y.tran, y.itran) != (x.tran, x.itran) {
            mism.push(format!(
                "{ctx} TRAN/ITRAN yfoil=({},{}) xfoil=({},{})",
                y.tran, y.itran, x.tran, x.itran
            ));
        }
        chk(&mut mism, &ctx, &y.ampl, &x.ampl, "AMPL(A1,A2,XT,ACRIT)", 0.0);
        chk(&mut mism, &ctx, &y.primary, &x.primary, "PRIMARY", 0.0);
        chk(&mut mism, &ctx, &y.kinematic, &x.kinematic, "KINEMATIC", 0.0);
        chk(&mut mism, &ctx, &y.closure, &x.closure, "CLOSURE", 0.0);
        // residuals are O(1) nondimensional equations: compare on the unit scale
        // VSREZ(4) is stale at the dump point (set only in the direct/inverse branch after it)
        chk(&mut mism, &ctx, &y.residual[..3], &x.residual[..3], "RESIDUAL", 1.0);
        for r in 0..3 {
            chk(&mut mism, &ctx, &y.vs2[r], &x.vs2[r], &format!("VS2_{}", r + 1), 0.0);
        }
        // Newton corrections: scale each by the magnitude of the variable it corrects
        let scales = [
            x.primary[4].abs().max(0.03),
            x.primary[2].abs(),
            x.primary[3].abs(),
            x.primary[1].abs(),
        ];
        for i in 0..4 {
            chk(
                &mut mism,
                &ctx,
                &[y.solution[i]],
                &[x.solution[i]],
                &format!("SOLUTION[{i}]"),
                scales[i],
            );
        }
        chk(&mut mism, &ctx, &[y.dmax, y.rlx], &[x.dmax, x.rlx], "RELAX", 0.0);
        if y.has_update != x.has_update {
            mism.push(format!(
                "{ctx} update path: yfoil has_update={} xfoil has_update={}",
                y.has_update, x.has_update
            ));
        } else if x.has_update {
            chk(&mut mism, &ctx, &y.updated, &x.updated, "UPDATED", 0.0);
        }
        if y.converged != x.converged {
            mism.push(format!(
                "{ctx} CONVERGED flag yfoil={} xfoil={}",
                y.converged, x.converged
            ));
        }
        compared += 1;
        if mism.len() > 40 {
            break;
        }
    }
    if trace.iters.len() != xf.len() {
        mism.push(format!(
            "iteration count: yfoil={} xfoil={}",
            trace.iters.len(),
            xf.len()
        ));
    }
    println!(
        "mrchue trace: {compared} station-iterations compared, {} mismatching values",
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
