//! S7 gate: SETBL's global Newton system (VA/VB/VDEL/VM/VZ) against the reference.
//!
//! `blsolv_input.dat` is SETBL's output as BLSOLV receives it, dumped on VISCAL iterations 1..3.
//! For each, XFOIL's complete state entering SETBL's march (`mrchdu_input_<k>.dat` plus the
//! pointer layer, DIJ and UINV) is loaded, SETBL is run, and every coefficient is compared.
//! `setbl_output_<k>.dat` covers what SETBL writes besides the system: ITRAN/XSSITR/TFORCE,
//! XOCTR/YOCTR/TINDEX, DULE, USAV and the per-station TAU/DIS/CTQ/DELT/USLP.

mod fixtures;
mod utilities;

use fixtures::blsolv_fixtures::parse_blsolv_input;
use fixtures::mrchdu_fixtures::parse_bl_dump;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_LINALG, TOL_SOLVER};
use yfoil::solver::setbl::setbl;

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// Row-scaled comparison of one system block: `|a-b| <= tol * max(|a|, |b|, scale)` where the
/// scale is the largest magnitude in the same (equation, column) slot over all stations.
fn check_block(name: &str, k: usize, ours: &[[[f64; 2]; 3]], theirs: &[[[f64; 2]; 3]], tol: f64) -> f64 {
    assert_eq!(ours.len(), theirs.len(), "call {k}: {name} length");
    let mut sc = [[0.0_f64; 2]; 3];
    for row in theirs {
        for eq in 0..3 {
            for c in 0..2 {
                sc[eq][c] = sc[eq][c].max(row[eq][c].abs());
            }
        }
    }
    let mut worst = 0.0_f64;
    for (iv, (a, b)) in ours.iter().zip(theirs).enumerate() {
        for eq in 0..3 {
            for c in 0..2 {
                if a[eq][c] == 0.0 && b[eq][c] == 0.0 {
                    continue; // structurally empty slot (e.g. the momentum row's Ctau column)
                }
                let e = (a[eq][c] - b[eq][c]).abs() / a[eq][c].abs().max(b[eq][c].abs()).max(sc[eq][c]);
                worst = worst.max(e);
                assert_within(
                    a[eq][c],
                    b[eq][c],
                    tol,
                    sc[eq][c],
                    &format!("call {k}: {name}(eq {}, col {}, iv {})", eq + 1, c + 1, iv + 1),
                );
            }
        }
    }
    worst
}

fn check_call(k: usize) {
    let (mut st, _params, d) = fixtures::state_before_setbl_march(k);
    let xf = parse_blsolv_input(&fixture_path("blsolv_input.dat"), k).expect("BLSOLV input fixture");
    let out = parse_bl_dump(&fixture_path(&format!("setbl_output_{k}.dat")));

    let r = setbl(&mut st);

    // the CL-dependence sensitivities SETBL forms first
    assert_eq!(r.re_clmr.to_bits(), d.real("RE_CLMR").to_bits(), "call {k}: RE_CLMR");
    assert_eq!(r.msq_clmr.to_bits(), d.real("MSQ_CLMR").to_bits(), "call {k}: MSQ_CLMR");

    // system shape
    assert_eq!(r.sys.nsys, xf.nsys, "call {k}: NSYS");
    assert_eq!(r.sys.ivte1, Some(xf.ivte1_0based()), "call {k}: IVTE1");
    assert_eq!(r.sys.ivz, Some(xf.ivz_0based()), "call {k}: IVZ");
    assert_eq!(r.sys.vaccel.to_bits(), xf.vaccel.to_bits(), "call {k}: VACCEL");
    assert_eq!(
        r.sys.arc_length.unwrap().to_bits(),
        xf.arc_length.unwrap().to_bits(),
        "call {k}: ARC_LENGTH"
    );

    // the blocks
    let w_va = check_block("VA", k, &r.sys.va, &xf.va, TOL_LINALG);
    let w_vb = check_block("VB", k, &r.sys.vb, &xf.vb, TOL_LINALG);
    let w_vdel = check_block("VDEL", k, &r.sys.vdel, &xf.vdel_in, TOL_SOLVER);

    // VZ (3x2)
    let mut w_vz = 0.0_f64;
    let sc_vz = xf.vz.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    for eq in 0..3 {
        for c in 0..2 {
            let (a, b) = (r.sys.vz[eq][c], xf.vz[eq][c]);
            if a == 0.0 && b == 0.0 {
                continue;
            }
            w_vz = w_vz.max((a - b).abs() / a.abs().max(b.abs()).max(sc_vz));
            assert_within(
                a,
                b,
                TOL_LINALG,
                sc_vz,
                &format!("call {k}: VZ(eq {}, col {})", eq + 1, c + 1),
            );
        }
    }

    // VM: row-scaled — each (iv, equation) row is scaled by its own largest |VM| over jv
    let mut w_vm = 0.0_f64;
    for iv in 0..xf.nsys {
        for eq in 0..3 {
            let sc = (0..xf.nsys).map(|jv| xf.vm[iv][jv][eq].abs()).fold(0.0_f64, f64::max);
            for jv in 0..xf.nsys {
                let (a, b) = (r.sys.vm[iv][jv][eq], xf.vm[iv][jv][eq]);
                if a == 0.0 && b == 0.0 {
                    continue;
                }
                w_vm = w_vm.max((a - b).abs() / a.abs().max(b.abs()).max(sc));
                assert_within(
                    a,
                    b,
                    TOL_SOLVER,
                    sc,
                    &format!("call {k}: VM(eq {}, jv {}, iv {})", eq + 1, jv + 1, iv + 1),
                );
            }
        }
    }

    // what SETBL writes besides the system
    assert_eq!(st.itran[1..], [out.int("ITRAN1"), out.int("ITRAN2")], "call {k}: ITRAN");
    assert_eq!(
        st.tforce[1..],
        [out.logical("TFORCE1"), out.logical("TFORCE2")],
        "call {k}: TFORCE"
    );
    for is in 1..=2 {
        for (name, ours) in [
            ("XSSITR", st.xssitr[is]),
            ("XOCTR", st.xoctr[is]),
            ("YOCTR", st.yoctr[is]),
            ("TINDEX", st.tindex[is]),
            ("DULE", r.dule[is]),
        ] {
            assert_within(
                ours,
                out.real(&format!("{name}{is}")),
                TOL_SOLVER,
                1.0,
                &format!("call {k}: {name}({is})"),
            );
        }
    }
    let names = ["UEDG", "TAU", "DIS", "CTQ", "DELT", "USLP"];
    let mut w_st = 0.0_f64;
    for is in 1..=2 {
        let nbl = out.nbl(is);
        assert_eq!(nbl, st.nbl[is], "call {k}: NBL({is})");
        let theirs = |j: usize, m: usize| {
            if m == 0 {
                out.bl[is][j][1]
            } else {
                out.blx[is][j][m - 1]
            }
        };
        for ibl in 2..=nbl {
            let ours = [
                st.uedg[is][ibl],
                st.tau[is][ibl],
                st.dis[is][ibl],
                st.ctq[is][ibl],
                st.delt[is][ibl],
                st.uslp[is][ibl],
            ];
            for (m, name) in names.iter().enumerate() {
                let (a, b) = (ours[m], theirs(ibl, m));
                let scale = (2..=nbl).map(|j| theirs(j, m).abs()).fold(0.0_f64, f64::max);
                w_st = w_st.max((a - b).abs() / a.abs().max(b.abs()).max(scale));
                assert_within(a, b, TOL_SOLVER, scale, &format!("call {k}: {name}({ibl},{is})"));
            }
        }
    }
    println!(
        "setbl call {k}: NSYS={} worst scaled errors VA {w_va:.2e} VB {w_vb:.2e} VDEL {w_vdel:.2e} VZ {w_vz:.2e} VM {w_vm:.2e} state {w_st:.2e}; ITRAN={:?}",
        xf.nsys,
        &st.itran[1..]
    );
}

#[test]
fn test_setbl_system_matches_xfoil_call_1() {
    check_call(1);
}

#[test]
fn test_setbl_system_matches_xfoil_call_2() {
    check_call(2);
}

#[test]
fn test_setbl_system_matches_xfoil_call_3() {
    check_call(3);
}
