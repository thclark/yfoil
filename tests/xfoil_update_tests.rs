//! S8 gate: UPDATE against the reference, as the last step of a replayed VISCAL iteration.
//!
//! For VISCAL iterations 1..3, XFOIL's state entering SETBL's march is loaded, then
//! SETBL → BLSOLV → UPDATE is run (the first two already gated at 1e-10 / bit-identity), and
//! the post-UPDATE state and diagnostics (RLX, RMSBL, RMXBL/VMXBL/IMXBL/ISMXBL, DAC, CL) are
//! compared with `update_output_<k>.dat`.

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::parse_bl_dump;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::bl::blsolv::solve_newton_system;
use yfoil::solver::setbl::assemble_newton_system;
use yfoil::solver::update::apply_newton_update;

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn check_call(k: usize) {
    let (mut st, _params, _d) = fixtures::state_before_setbl_march(k);
    let o = parse_bl_dump(&fixture_path(&format!("update_output_{k}.dat")));

    let r = assemble_newton_system(&mut st);
    let sol = solve_newton_system(r.newton);
    let u = apply_newton_update(&mut st, &sol.deltas, r.mach_d_cl);

    // diagnostics
    for (name, ours) in [
        ("RLX", u.relaxation),
        ("RMSBL", u.residual),
        ("RMXBL", u.residual_max),
        ("DAC", u.free_variable_change),
        ("CLNEW", u.cl_new),
        ("CL_A", u.cl_d_alpha),
        ("CL_MS", u.cl_d_machsqd),
        ("CL_AC", u.cl_d_free),
        ("CL", st.cl),
        ("ALFA", st.alpha),
    ] {
        assert_within(ours, o.real(name), TOL_SOLVER, 1.0, &format!("call {k}: {name}"));
    }
    assert_eq!(
        u.residual_max_variable.to_string(),
        o.header["VMXBL"],
        "call {k}: VMXBL"
    );
    assert_eq!(u.i_residual_max_station, o.int("IMXBL"), "call {k}: IMXBL");
    assert_eq!(u.residual_max_side, o.int("ISMXBL"), "call {k}: ISMXBL");
    assert_eq!(r.mach_d_cl.to_bits(), o.real("MINF_CL").to_bits(), "call {k}: MINF_CL");

    // state: side 1 includes the wake stations equated from side 2
    let names = [
        "UEDG", "THET", "DSTR", "CTAU", "MASS", "TAU", "DIS", "CTQ", "DELT", "TSTR",
    ];
    let mut worst = (0.0_f64, "", 0, 0);
    for is in 1..=2 {
        let nrows = o.nbl(is);
        let theirs = |j: usize, m: usize| {
            if m < 5 {
                o.bl[is][j][m + 1]
            } else {
                o.blx[is][j][m - 5]
            }
        };
        for ibl in 2..=nrows {
            let ours = [
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
            for (m, name) in names.iter().enumerate() {
                let (a, b) = (ours[m], theirs(ibl, m));
                let scale = (2..=nrows).map(|j| theirs(j, m).abs()).fold(0.0_f64, f64::max);
                let e = (a - b).abs() / a.abs().max(b.abs()).max(scale);
                if e > worst.0 {
                    worst = (e, name, is, ibl);
                }
                assert_within(a, b, TOL_SOLVER, scale, &format!("call {k}: {name}({ibl},{is})"));
            }
        }
    }
    println!(
        "update call {k}: RLX={:.6} RMSBL={:.6e} RMXBL={:+.4e} {}@({},{}) DAC={:+.4e}; state within {TOL_SOLVER:.0e} (worst {:.2e} {} at ({},{}))",
        u.relaxation, u.residual, u.residual_max, u.residual_max_variable, u.i_residual_max_station, u.residual_max_side, u.free_variable_change, worst.0, worst.1, worst.3, worst.2
    );
}

#[test]
fn test_update_matches_xfoil_iteration_1() {
    check_call(1);
}

#[test]
fn test_update_matches_xfoil_iteration_2() {
    check_call(2);
}

#[test]
fn test_update_matches_xfoil_iteration_3() {
    check_call(3);
}
