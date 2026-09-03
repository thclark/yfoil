//! S6 gate: MRCHDU (mixed-mode march along the Ue-Hk characteristic) against the reference.
//!
//! Two views of the same run: for SETBL calls 1..3, the BL state before MRCHDU
//! (`mrchdu_input_<k>.dat`) is loaded, marched, and compared with the state after
//! (`mrchdu_output_<k>.dat`); and on call 1 the per-iteration Newton trace
//! (`xfoil_mrchdu_trace.dat`) is compared station by station, iteration by iteration.

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::{parse_bl_dump, parse_mrchdu_trace, BlDump};
use fixtures::pointers_fixtures::parse_pointers;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::bl::mrchdu::{mrchdu, MrchduTrace};
use yfoil::bl::mrchue::mrchue;
use yfoil::bl::system::BLGlobalParams;
use yfoil::solver::blstate::BlState;

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// XFOIL's complete MRCHDU input on SETBL call `k`: the pointer layer (geometry, IPAN, WGAP —
/// unchanged while IST is unchanged, which the NBL check below asserts) and the BL arrays
/// exactly as SETBL had them.
fn state_before_mrchdu(k: usize) -> (BlState, BLGlobalParams, [f64; 3], BlDump) {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let d = parse_bl_dump(&fixture_path(&format!("mrchdu_input_{k}.dat")));
    assert_eq!(
        [d.int("NBL1"), d.int("NBL2")],
        [f.nbl[1], f.nbl[2]],
        "call {k}: NBL vs pointer dump"
    );
    assert_eq!(
        [d.int("IBLTE1"), d.int("IBLTE2")],
        [f.iblte[1], f.iblte[2]],
        "call {k}: IBLTE"
    );
    let mut st = BlState::empty(f.n, f.nw);
    st.x = f.x.clone();
    st.y = f.y.clone();
    st.s = f.s.clone();
    st.ist = f.ist;
    st.sst = f.sst;
    st.nbl = f.nbl;
    st.iblte = f.iblte;
    st.ipan = f.ipan.clone();
    st.vti = f.vti.clone();
    st.isys = f.isys.clone();
    st.wgap = f.wgap.clone();
    st.ante = f.ante;
    st.aste = f.aste;
    st.dste = f.dste;
    st.sharp = f.sharp;
    st.chord = f.chord;
    st.sle = f.sle;
    st.xle = f.xle;
    st.yle = f.yle;
    st.xte = f.xte;
    st.yte = f.yte;
    st.xstrip = [0.0, d.real("XSTRIP1"), d.real("XSTRIP2")];
    st.itran = [0, d.int("ITRAN1"), d.int("ITRAN2")];
    for is in 1..=2 {
        for ibl in 1..=f.nbl[is] {
            let r = d.bl[is][ibl];
            st.xssi[is][ibl] = r[0];
            st.uedg[is][ibl] = r[1];
            st.thet[is][ibl] = r[2];
            st.dstr[is][ibl] = r[3];
            st.ctau[is][ibl] = r[4];
            st.mass[is][ibl] = r[5];
            // STMOVE recomputes XSSI every VISCAL iteration (SST moves); the pointer dump is
            // from the first IBLSYS call, so only call 1 can be checked against it.
            if k == 1 {
                assert_eq!(
                    f.xssi[is][ibl].to_bits(),
                    r[0].to_bits(),
                    "call {k}: XSSI({ibl},{is}) differs between the pointer dump and the BL dump"
                );
            }
        }
    }
    let params = BLGlobalParams::new(d.real("MINF"), d.real("REINF"), 1.4);
    assert_eq!(params.reybl.to_bits(), d.real("REYBL").to_bits(), "REYBL");
    assert_eq!(params.hstinv.to_bits(), d.real("HSTINV").to_bits(), "HSTINV");
    assert_eq!(params.gm1.to_bits(), d.real("GM1BL").to_bits(), "GM1BL");
    let acrit = [0.0, d.real("ACRIT1"), d.real("ACRIT2")];
    if k == 1 {
        // Chained replay: XFOIL's COM1/COM2/XT COMMON state entering the first MRCHDU is what
        // MRCHUE left behind, so run MRCHUE from its own gated input and take that state
        // (the BL arrays it produces equal the dump to 1e-10 — the S5 gate — but the dump is
        // authoritative for the arrays themselves).
        let (mut pre, p2, a2) = fixtures::state_before_mrchue();
        mrchue(&mut pre, &p2, a2, None);
        st.com1 = pre.com1;
        st.com2 = pre.com2;
        st.xt = pre.xt;
    }
    (st, params, acrit, d)
}

fn check_state_after(k: usize) {
    let (mut st, params, acrit, _) = state_before_mrchdu(k);
    let o = parse_bl_dump(&fixture_path(&format!("mrchdu_output_{k}.dat")));
    mrchdu(&mut st, &params, acrit, None);
    assert_eq!(st.itran[1..], [o.int("ITRAN1"), o.int("ITRAN2")], "call {k}: ITRAN");
    assert_eq!(
        st.tforce[1..],
        [o.logical("TFORCE1"), o.logical("TFORCE2")],
        "call {k}: TFORCE"
    );
    for is in 1..=2 {
        assert_within(
            st.xssitr[is],
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
        assert_eq!(nbl, st.nbl[is], "call {k}: NBL({is})");
        for ibl in 2..=nbl {
            let ours = [
                st.xssi[is][ibl],
                st.uedg[is][ibl],
                st.thet[is][ibl],
                st.dstr[is][ibl],
                st.ctau[is][ibl],
                st.mass[is][ibl],
                st.tau[is][ibl],
                st.dis[is][ibl],
                st.ctq[is][ibl],
                st.delt[is][ibl],
                st.tstr[is][ibl],
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
        &st.itran[1..]
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
    let (mut st, params, acrit, _) = state_before_mrchdu(1);
    let xf = parse_mrchdu_trace(&fixture_path("xfoil_mrchdu_trace.dat"));
    let mut tr = MrchduTrace::default();
    mrchdu(&mut st, &params, acrit, Some(&mut tr));

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
