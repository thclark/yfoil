//! S9 gate: the VISCAL Newton loop against the reference, as a replay.
//!
//! Seeded with XFOIL's complete state entering the first SETBL (BL arrays after MRCHUE, the
//! pointer layer, DIJ, the inviscid arrays of `viscal_inviscid.dat`), the loop
//! SETBL → BLSOLV → UPDATE → MRCL/COMSET → QVFUE → GAMQV → STMOVE → CLCALC → CDCALC is run
//! to convergence and every iteration's RMSBL/RLX/RMXBL/CL/CM/CD/CDF/CDP/IST/ITRAN/XOCTR is
//! compared with `viscal_iter.dat`; the iteration count must be identical. The final Cp,
//! velocity and BL arrays are compared with `viscal_final.dat`.

mod fixtures;
mod utilities;

use fixtures::mrchdu_fixtures::parse_bl_dump;
use fixtures::pointers_fixtures::parse_viscal_inviscid;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::clcalc::set_compressibility;
use yfoil::solver::viscal::{solve_viscous, IterationRecord};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// One `ITER=` block of viscal_iter.dat as a key/value map.
fn parse_iter_blocks(path: &PathBuf) -> Vec<std::collections::HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = Vec::new();
    for l in text.lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "ITER" {
            out.push(std::collections::HashMap::new());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

fn real(m: &std::collections::HashMap<String, String>, k: &str) -> f64 {
    m[k].parse().unwrap()
}
fn int(m: &std::collections::HashMap<String, String>, k: &str) -> usize {
    m[k].parse().unwrap()
}

/// XFOIL's state entering VISCAL's Newton loop on the reference run.
fn state_entering_newton_loop() -> SolverState {
    let (mut st, _params, _d) = fixtures::state_before_setbl_march(1);
    let inv = parse_viscal_inviscid(&fixture_path("viscal_inviscid.dat"));
    let nt = st.n_foil_nodes + st.n_wake_nodes;
    assert_eq!(inv.qinv.len(), nt + 1, "viscal_inviscid.dat node count");
    for i in 1..=nt {
        st.q_inviscid_basis[1][i] = inv.qinvu1[i];
        st.q_inviscid_basis[2][i] = inv.qinvu2[i];
        st.q_inviscid[i] = inv.qinv[i];
        st.q_inviscid_d_alpha[i] = inv.qinv_a[i];
        if i <= st.n_foil_nodes {
            st.gamma[i] = inv.gam[i];
            st.gamma_d_alpha[i] = inv.gam_a[i];
        }
    }
    st.alpha = inv.header["ALFA"].parse().unwrap();
    st.qinf = inv.header["QINF"].parse().unwrap();
    st.mach = inv.header["MINF"].parse().unwrap();
    st.mach_cl1 = st.mach;
    st.cl = inv.header["CL"].parse().unwrap();
    st.mach_d_cl = inv.header["MINF_CL"].parse().unwrap();
    // the dump is taken after UICALC, before QDCALC and SETBL's MRCHUE; the replay starts
    // with both done (DIJ loaded, BL arrays from the post-MRCHUE dump)
    st.wake_built = inv.header["LWAKE"] == "T";
    st.pointers_built = inv.header["LIPAN"] == "T";
    st.converged = inv.header["LVCONV"] == "T";
    st.bl_initialised = true;
    st.dij_wake_built = true;
    st.viscous = true;
    set_compressibility(&mut st);
    assert_eq!(
        st.karman_tsien.to_bits(),
        inv.header["TKLAM"].parse::<f64>().unwrap().to_bits(),
        "TKLAM"
    );
    assert!(
        st.wake_built && st.pointers_built && st.bl_initialised && st.dij_wake_built,
        "replay expects the prologue done"
    );
    st
}

#[test]
fn test_viscal_replay_matches_xfoil_iteration_by_iteration() {
    let mut st = state_entering_newton_loop();
    let xf = parse_iter_blocks(&fixture_path("viscal_iter.dat"));
    let fin = parse_bl_dump(&fixture_path("viscal_final.dat"));

    let mut tr: Vec<IterationRecord> = Vec::new();
    let converged = solve_viscous(&mut st, None, 20, 1.0, Some(&mut tr));

    // identical iteration count and convergence outcome
    assert_eq!(tr.len(), xf.len(), "VISCAL iteration count");
    assert_eq!(converged, fin.logical("LVCONV"), "LVCONV");
    assert_eq!(tr.len(), fin.int("NITDONE"), "NITDONE");

    for (y, x) in tr.iter().zip(&xf) {
        let it = int(x, "ITER");
        assert_eq!(y.iteration, it);
        let ctx = format!("iteration {it}");
        for (name, ours) in [
            ("RMSBL", y.residual),
            ("RMXBL", y.residual_max),
            ("RLX", y.relaxation),
            ("ALFA", y.alpha),
            ("MINF", y.mach),
            ("REINF", y.re),
            ("CL", y.cl),
            ("CM", y.cm),
            ("CD", y.cd),
            ("CDF", y.cd_friction),
            ("CDP", y.cd_pressure),
            ("CL_ALF", y.cl_d_alpha),
            ("CL_MSQ", y.cl_d_machsqd),
            ("SST", y.s_stagnation),
            ("XOCTR1", y.x_transition[1]),
            ("XOCTR2", y.x_transition[2]),
        ] {
            let scale = if name == "REINF" { real(x, name) } else { 1.0 };
            assert_within(ours, real(x, name), TOL_SOLVER, scale, &format!("{ctx}: {name}"));
        }
        assert_eq!(y.residual_max_variable.to_string(), x["VMXBL"], "{ctx}: VMXBL");
        assert_eq!(y.i_residual_max_station, int(x, "IMXBL"), "{ctx}: IMXBL");
        assert_eq!(y.residual_max_side, int(x, "ISMXBL"), "{ctx}: ISMXBL");
        assert_eq!(y.i_stagnation_node, int(x, "IST"), "{ctx}: IST");
        assert_eq!(
            y.i_transition_station[1..],
            [int(x, "ITRAN1"), int(x, "ITRAN2")],
            "{ctx}: ITRAN"
        );
        assert_eq!(y.converged, x["LVCONV"] == "T", "{ctx}: LVCONV");
        println!(
            "viscal iter {it}: rms {:.4e} max {:+.4e} {} at {:3}{:2}  CL {:.6} CD {:.6}  (RLX {:.3}, IST {})",
            y.residual,
            y.residual_max,
            y.residual_max_variable,
            y.i_residual_max_station,
            y.residual_max_side,
            y.cl,
            y.cd,
            y.relaxation,
            y.i_stagnation_node
        );
    }

    // final state
    for (name, ours) in [
        ("ALFA", st.alpha),
        ("CL", st.cl),
        ("CM", st.cm),
        ("CD", st.cd),
        ("CDF", st.cd_friction),
        ("CDP", st.cd_pressure),
        ("CL_ALF", st.cl_d_alpha),
        ("CL_MSQ", st.cl_d_machsqd),
        ("AVISC", st.alpha_converged),
        ("MVISC", st.mach_converged),
        ("XOCTR1", st.x_transition[1]),
        ("XOCTR2", st.x_transition[2]),
    ] {
        assert_within(ours, fin.real(name), TOL_SOLVER, 1.0, &format!("final: {name}"));
    }
    assert_eq!(st.i_stagnation_node, fin.int("IST"), "final: IST");
    assert_eq!(
        st.i_transition_station[1..],
        [fin.int("ITRAN1"), fin.int("ITRAN2")],
        "final: ITRAN"
    );

    // per node: CPI CPV QINV QVIS GAM
    let text = std::fs::read_to_string(fixture_path("viscal_final.dat")).unwrap();
    let mut nodes = 0;
    for l in text.lines() {
        let Some(r) = l.strip_prefix("NODE(") else { continue };
        let (idx, vals) = r.split_once(")=").unwrap();
        let i: usize = idx.trim().parse().unwrap();
        let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
        for (name, ours, theirs) in [
            ("CPI", st.cp_inviscid[i], v[0]),
            ("CPV", st.cp_viscous[i], v[1]),
            ("QINV", st.q_inviscid[i], v[2]),
            ("QVIS", st.q_viscous[i], v[3]),
        ] {
            assert_within(ours, theirs, TOL_SOLVER, 1.0, &format!("final: {name}({i})"));
        }
        if i <= st.n_foil_nodes {
            // GAM is only defined on the airfoil nodes
            assert_within(st.gamma[i], v[4], TOL_SOLVER, 1.0, &format!("final: GAM({i})"));
        }
        nodes += 1;
    }
    assert_eq!(nodes, st.n_foil_nodes + st.n_wake_nodes, "final: node count");

    // per station: XSSI UEDG THET DSTR CTAU MASS
    let names = ["XSSI", "UEDG", "THET", "DSTR", "CTAU", "MASS"];
    for is in 1..=2 {
        let nbl = fin.nbl(is);
        assert_eq!(nbl, st.n_stations[is], "final: NBL({is})");
        for ibl in 2..=nbl {
            let ours = [
                st.xi[is][ibl],
                st.ue[is][ibl],
                st.theta[is][ibl],
                st.dstar[is][ibl],
                st.sqrtctau[is][ibl],
                st.mass_defect[is][ibl],
            ];
            for (m, name) in names.iter().enumerate() {
                let scale = (2..=nbl).map(|j| fin.bl[is][j][m].abs()).fold(0.0_f64, f64::max);
                assert_within(
                    ours[m],
                    fin.bl[is][ibl][m],
                    TOL_SOLVER,
                    scale,
                    &format!("final: {name}({ibl},{is})"),
                );
            }
        }
    }
    println!(
        "viscal: {} iterations, converged={converged}; CL {:.8} CD {:.8} CDF {:.8} CDP {:.8} CM {:.8}; XTR {:.6} {:.6}",
        tr.len(),
        st.cl,
        st.cd,
        st.cd_friction,
        st.cd_pressure,
        st.cm,
        st.x_transition[1],
        st.x_transition[2]
    );
}
