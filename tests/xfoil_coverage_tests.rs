//! Rule 6 coverage cases: each exercises branches the reference case leaves dead, and each is
//! gated the same way — the `Session` run must reproduce every VISCAL call's iteration count,
//! converged flag, IST/ITRAN and per-iteration RMSBL/RLX/CL/CD/CM/ALFA, and the converged point,
//! against that case's `viscal_points.dat` / `viscal_iters_all.dat`.
//!
//! | case | dead branch it wakes up |
//! |---|---|
//! | `naca0012_n60_sharp_a2_re1e6` | SHARP TE: GGCALC/QDCALC bisector rows, TECALC with zero gap |
//! | `naca0012_n60_a2_re1e6_m03` | M = 0.3: Kármán–Tsien, every `*_MS` sensitivity, HSTINV/REYBL_MS |
//! | `naca0012_n60_a12_re1e6` | high alpha: RLX limiting, Hk clamps, negative-Ue fix-up |
//! | `naca0012_n60_a4_re1e5` | low Re: laminar separation, MRCHUE DIRECT switch (HLMAX/HTMAX) |
//! | `naca0012_n60_a2_re1e6_xtr03` | XTR 0.3 0.3: XIFSET with XSTRIP < 1, TRCHEK2 forced transition |
//! | `naca0012_n60_a2_repeat_re1e6` | the same alpha twice: VISCAL entered with LWAKE/LWDIJ/LVCONV already set |
//! | `naca0012_n60_a2_cl03_re1e6` | CL after a converged alpha: SPECCL with LGAMU/LQAIJ valid |
//! | `naca0012_n60_a2_re1e6_type3` | TYPE 3 (MATYP 1, RETYP 3): MRCL's 1/CL Re scaling |
//! | `naca0012_n60_a2_re1e6_m03_type2` | MATYP = 2 with M > 0: SPECAL's Mach–CL Newton loop |
//! | `naca0012_n60_a2_re1e6_xtr_coinc` | trip in the natural-transition interval: TRCHEK2 TRFREE .AND. TRFORC |
//! | `naca0012_n60_a2_re1e6_damp` | OPER DAMP: IDAMPV = 1, the modified-envelope amplification (DAMPL2) |
//!
//! The last five were added from the gcov measurement (`cargo xtask coverage`,
//! `docs/validation/coverage.md`): each wakes branches the earlier set left open.

use yfoil::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
mod fixtures;
mod utilities;

use std::path::PathBuf;
use utilities::records::{check_call, load, transient_tol, Outcome};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};

fn case_dir(case: &str) -> PathBuf {
    fixtures::require_fixture(&format!("tests/fixtures/xfoil/{case}"))
}

/// Run the case's ALFA sequence through one session and check every VISCAL call.
fn run_case(case: &str, spec: FlowConditions, alphas_deg: &[f64]) -> Vec<Outcome> {
    let dir = case_dir(case);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    let rec = load(&dir);
    assert_eq!(rec.points.len(), alphas_deg.len(), "{case}: VISCAL call count");
    let mut session = Session::new(&airfoil, spec);
    let mut outcomes = Vec::new();
    for (i, a) in alphas_deg.iter().enumerate() {
        let p = session.alpha(a.to_radians());
        outcomes.push(check_call(&rec, i + 1, &p, session.state(), transient_tol(i + 1)));
    }
    outcomes
}

#[test]
fn test_sharp_trailing_edge_matches_xfoil() {
    let dir = case_dir("naca0012_n60_sharp_a2_re1e6");
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    assert!(airfoil.sharp_te, "the case must actually take the SHARP path");
    let outcomes = run_case("naca0012_n60_sharp_a2_re1e6", FlowConditions::default(), &[2.0]);
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_mach_0_3_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_m03",
        FlowConditions {
            mach: 0.3,
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_high_alpha_separated_matches_xfoil() {
    // The reference's own 1-ULP trajectory becomes hypersensitive from iteration 18 (floor 5e-5
    // → 1.4e-3, see noise_floor.json); yFoil tracks it like a ~2-ULP perturbation through
    // iteration 18 and parts at 19. The two replay tests below show the steps themselves are
    // faithful, so the documented outcome is threshold-straddling at iteration 19.
    let outcomes = run_case("naca0012_n60_a12_re1e6", FlowConditions::default(), &[12.0]);
    match &outcomes[0] {
        Outcome::Straddling { at_iteration, .. } => assert_eq!(*at_iteration, 19, "{outcomes:?}"),
        Outcome::Match => println!("12° case now matches to the end (the straddle closed)"),
    }
}

#[test]
fn test_low_re_laminar_separation_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a4_re1e5",
        FlowConditions {
            re: Some(1.0e5),
            ..FlowConditions::default()
        },
        &[4.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_forced_transition_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_xtr03",
        FlowConditions {
            x_trip: [0.3, 0.3],
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

/// Replay one VISCAL iteration of the 12° case from XFOIL's exact state entering SETBL call
/// `k` (`mrchdu_input_<k>.dat`) and compare with XFOIL's state after UPDATE `k`
/// (`update_output_<k>.dat`). This isolates the iteration where the unconverged trajectories
/// part (18 → 19) from everything accumulated before it: if the one-step map agrees, the
/// divergence is accumulated differences crossing a threshold (threshold-straddling); if not,
/// the step itself is wrong.
fn replay_iteration(case: &str, alpha_deg: f64, k: usize) {
    use fixtures::mrchdu_fixtures::parse_bl_dump;
    use utilities::tolerances::{assert_within, TOL_SOLVER};
    use yfoil::bl::blsolv::solve_newton_system;
    use yfoil::solver::clcalc::set_compressibility;
    use yfoil::solver::pointers::{map_stations_to_nodes, map_stations_to_rows, set_station_xi};
    use yfoil::solver::setbl::assemble_newton_system;
    use yfoil::solver::update::apply_newton_update;
    use yfoil::solver::velocity::set_ue_inviscid;
    use yfoil::solver::viscal::solve_viscous;

    let dir = case_dir(case);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let d = parse_bl_dump(&dir.join(format!("mrchdu_input_{k}.dat")));
    let o = parse_bl_dump(&dir.join(format!("update_output_{k}.dat")));

    // prologue only (wake, QINVU/QINV, pointers, UINV, DIJ) — then XFOIL's state at call k
    let mut session = Session::new(&airfoil, FlowConditions::default());
    {
        let (st, sys) = session.parts_mut();
        yfoil::solver::specal::alpha_command(st, sys, alpha_deg.to_radians());
        solve_viscous(st, sys.as_mut(), 0, 1.0, None);
    }
    let st = session.state_mut();
    st.bl_initialised = true;
    st.i_stagnation_node = d.int("IST");
    st.s_stagnation = d.real("SST");
    st.s_stagnation_d_gamma_node0 = d.real("SST_GO");
    st.s_stagnation_d_gamma_node1 = d.real("SST_GP");
    map_stations_to_nodes(st);
    set_station_xi(st);
    map_stations_to_rows(st);
    set_ue_inviscid(st);
    assert_eq!(
        [d.int("NBL1"), d.int("NBL2")],
        [st.n_stations[1], st.n_stations[2]],
        "NBL from the dumped IST"
    );
    st.i_transition_station = [0, d.int("ITRAN1"), d.int("ITRAN2")];
    st.cl = d.real("CLMR");
    st.mach_d_cl = 0.0;
    set_compressibility(st);
    for is in 1..=2 {
        for ibl in 1..=st.n_stations[is] {
            let r = d.bl[is][ibl];
            assert_within(
                st.xi[is][ibl],
                r[0],
                TOL_SOLVER,
                1.0,
                &format!("XSSI({ibl},{is}) rebuilt from IST/SST"),
            );
            st.xi[is][ibl] = r[0];
            st.ue[is][ibl] = r[1];
            st.theta[is][ibl] = r[2];
            st.dstar[is][ibl] = r[3];
            st.sqrtctau[is][ibl] = r[4];
            st.mass_defect[is][ibl] = r[5];
        }
    }

    // one iteration
    let r = assemble_newton_system(st);
    let sol = solve_newton_system(r.newton);
    let u = apply_newton_update(st, &sol.deltas, st.mach_d_cl);

    // the reference's own 1-ULP spread of this iteration's values bounds what one step can be
    // expected to reproduce (the replay's only foreign input is yFoil's DIJ, ~5e-11)
    let rec = load(&dir);
    let fl = |name: &str| {
        rec.floor
            .as_ref()
            .and_then(|f| f.calls.first())
            .and_then(|c| c.iterations.get(k - 1))
            .and_then(|m| m.get(name))
            .copied()
            .unwrap_or(0.0)
    };
    for (name, ours) in [
        ("RLX", u.relaxation),
        ("RMSBL", u.residual),
        ("CL", st.cl),
        ("RMXBL", u.residual_max),
    ] {
        let (a, b) = (ours, o.real(name));
        let lim = utilities::records::allowed(a, b, TOL_SOLVER, 1.0, fl(name));
        assert!(
            (a - b).abs() <= lim,
            "{case} replay call {k}: {name}: yfoil={a:.17e} xfoil={b:.17e} |diff|={:.3e} > allowed {lim:.3e}",
            (a - b).abs()
        );
    }
    // DAC has no recorded floor; its effect is gated through CL (CL += RLX*DAC) and RLX above
    let _ = u.free_variable_change;
    assert_eq!(
        u.residual_max_variable.to_string(),
        o.header["VMXBL"],
        "{case} replay call {k}: VMXBL"
    );
    assert_eq!(
        (u.i_residual_max_station, u.residual_max_side),
        (o.int("IMXBL"), o.int("ISMXBL")),
        "{case} replay call {k}: IMXBL/ISMXBL"
    );
    let names = ["UEDG", "THET", "DSTR", "CTAU", "MASS"];
    let arr_floor = |name: &str| {
        rec.floor
            .as_ref()
            .and_then(|f| f.update_output.get(&k))
            .and_then(|m| m.get(name))
            .copied()
            .unwrap_or(0.0)
    };
    let mut worst = (0.0_f64, "", 0, 0);
    for is in 1..=2 {
        let nrows = o.nbl(is);
        for ibl in 2..=nrows {
            let ours = [
                st.ue[is][ibl],
                st.theta[is][ibl],
                st.dstar[is][ibl],
                st.sqrtctau[is][ibl],
                st.mass_defect[is][ibl],
            ];
            for (m, name) in names.iter().enumerate() {
                let b = o.bl[is][ibl][m + 1];
                let scale = (2..=nrows).map(|j| o.bl[is][j][m + 1].abs()).fold(0.0_f64, f64::max);
                let e = (ours[m] - b).abs() / ours[m].abs().max(b.abs()).max(scale);
                if e > worst.0 {
                    worst = (e, name, is, ibl);
                }
                let lim = utilities::records::allowed(ours[m], b, TOL_SOLVER, scale, arr_floor(name));
                assert!(
                    (ours[m] - b).abs() <= lim,
                    "{case} replay call {k}: {name}({ibl},{is}): yfoil={:.17e} xfoil={b:.17e} |diff|={:.3e} > allowed {lim:.3e} (array floor {:.2e})",
                    ours[m],
                    (ours[m] - b).abs(),
                    arr_floor(name)
                );
            }
        }
    }
    println!(
        "{case} replay of iteration {k} from XFOIL's state: RLX {:.6} RMSBL {:.6e} {}@({},{}) CL {:.8} — every array within max({TOL_SOLVER:.0e}·scale, 4× its own 1-ULP floor) (worst rel {:.2e} {} at ({},{}))",
        u.relaxation, u.residual, u.residual_max_variable, u.i_residual_max_station, u.residual_max_side, st.cl, worst.0, worst.1, worst.3, worst.2
    );
}

#[test]
fn test_high_alpha_iteration_18_replay_from_xfoil_state() {
    replay_iteration("naca0012_n60_a12_re1e6", 12.0, 18);
}

#[test]
fn test_high_alpha_iteration_19_replay_from_xfoil_state() {
    replay_iteration("naca0012_n60_a12_re1e6", 12.0, 19);
}

#[test]
fn test_repeated_alpha_matches_xfoil() {
    let outcomes = run_case("naca0012_n60_a2_repeat_re1e6", FlowConditions::default(), &[2.0, 2.0]);
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_cl_after_alpha_matches_xfoil() {
    let case = "naca0012_n60_a2_cl03_re1e6";
    let dir = case_dir(case);
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let rec = load(&dir);
    assert_eq!(rec.points.len(), 2, "{case}: VISCAL call count");
    let mut session = Session::new(&airfoil, FlowConditions::default());
    let p = session.alpha(2.0_f64.to_radians());
    let o1 = check_call(&rec, 1, &p, session.state(), transient_tol(1));
    let p = session.cl(0.3);
    let o2 = check_call(&rec, 2, &p, session.state(), transient_tol(2));
    assert_eq!((o1, o2), (Outcome::Match, Outcome::Match));
}

#[test]
fn test_type_3_matches_xfoil() {
    // XFOIL's TYPE 3 is MATYP = 1, RETYP = 3 (xoper.f:362); MATYP = 3 is never set by TYPE
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_type3",
        FlowConditions {
            mach_cl_dependence: MachClDependence::Fixed,
            re_cl_dependence: ReClDependence::InverseCl,
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_matyp_2_with_mach_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_m03_type2",
        FlowConditions {
            mach: 0.3,
            mach_cl_dependence: MachClDependence::InverseSqrtCl,
            re_cl_dependence: ReClDependence::InverseSqrtCl,
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_trip_in_transition_interval_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_xtr_coinc",
        FlowConditions {
            x_trip: [0.48, 0.87],
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}

#[test]
fn test_modified_amplification_damp_matches_xfoil() {
    let outcomes = run_case(
        "naca0012_n60_a2_re1e6_damp",
        FlowConditions {
            amplification_model: AmplificationModel::ModifiedEnvelope,
            ..FlowConditions::default()
        },
        &[2.0],
    );
    assert!(outcomes.iter().all(|o| *o == Outcome::Match), "{outcomes:?}");
}
