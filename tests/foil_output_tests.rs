//! The foil output of an operating point (`AnalysisOutput::from_session`): geometry with the
//! wake, every per-station BL quantity on both sides and the wake, the markers, and the polar
//! sweep's embedded distributions. The station layout is checked against the pointer layer, the
//! live closure columns against XFOIL's identities (M = 0: Ue = UEDG/QINF and Hk = H bit-exact),
//! the CPDISP wake split against its own formula, and the live-versus-stored lag against the final
//! Newton residual (see the module docs of `yfoil::output::foil`).

mod fixtures;
mod utilities;

use std::path::PathBuf;
use utilities::tolerances::TOL_SOLVER;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::output::{AnalysisOutput, BlQuantity, PolarOutput, SeparationKind};
use yfoil::solver::analysis::{compute_polar_with, FlowConditions, PolarConfig, Session};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

fn spec() -> FlowConditions {
    FlowConditions {
        re: 1.0e6,
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    }
}

fn ref_case_session() -> Session {
    let geometry = read_geometry_from_file(fixture_path("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);
    Session::new(&airfoil, spec())
}

#[test]
fn analysis_output_carries_geometry_wake_and_every_station() {
    let mut session = ref_case_session();
    let p = session.alpha(2.0_f64.to_radians());
    assert!(p.converged, "the reference case converges");
    let out = AnalysisOutput::from_session(&session, &p, "naca0012", &spec(), false);
    let st = &session.state;

    // geometry: airfoil nodes and the wake, as the solver holds them
    assert_eq!(out.geometry.n(), st.n_foil_nodes);
    assert_eq!(out.geometry.x, st.x[1..=st.n_foil_nodes]);
    assert_eq!(out.geometry.nx, st.normal_x[1..=st.n_foil_nodes]);
    let wake = out.geometry.wake.as_ref().expect("viscous point has a wake");
    assert_eq!(wake.x.len(), st.n_wake_nodes);
    assert_eq!(wake.x, st.x[st.n_foil_nodes + 1..=st.n_foil_nodes + st.n_wake_nodes]);
    assert_eq!(
        wake.nx,
        st.normal_x[st.n_foil_nodes + 1..=st.n_foil_nodes + st.n_wake_nodes]
    );

    // stations: side 1 and 2 up to IBLTE, the wake after it, one entry per column
    let bl = out.boundary_layer.as_ref().expect("viscous point has a boundary layer");
    for side in [&bl.upper, &bl.lower, &bl.wake] {
        side.check_lengths().unwrap();
    }
    assert_eq!(bl.iblte, [st.i_te_station[1], st.i_te_station[2]]);
    assert_eq!(bl.upper.len(), st.i_te_station[1] - 1);
    assert_eq!(bl.lower.len(), st.i_te_station[2] - 1);
    assert_eq!(bl.wake.len(), st.n_wake_nodes);
    assert_eq!(bl.nw, st.n_wake_nodes);
    assert_eq!(bl.itran, [st.i_transition_station[1], st.i_transition_station[2]]);
    assert_eq!(bl.upper.ibl[0], 2);
    assert_eq!(*bl.wake.ibl.first().unwrap(), st.i_te_station[2] + 1);

    // node pointers map onto the geometry bit-exactly (wake nodes onto the wake)
    for side in [&bl.upper, &bl.lower] {
        for (k, &node) in side.node.iter().enumerate() {
            assert!(node >= 1 && node <= st.n_foil_nodes);
            assert_eq!(side.x[k].to_bits(), out.geometry.x[node - 1].to_bits());
            assert_eq!(side.y[k].to_bits(), out.geometry.y[node - 1].to_bits());
        }
    }
    for (k, &node) in bl.wake.node.iter().enumerate() {
        assert!(node > st.n_foil_nodes);
        assert_eq!(bl.wake.x[k].to_bits(), wake.x[node - st.n_foil_nodes - 1].to_bits());
    }

    // primaries are the state's arrays verbatim
    for (k, &ibl) in bl.upper.ibl.iter().enumerate() {
        assert_eq!(bl.upper.thet[k].to_bits(), st.theta[1][ibl].to_bits());
        assert_eq!(bl.upper.dstr[k].to_bits(), st.dstar[1][ibl].to_bits());
        assert_eq!(bl.upper.uedg[k].to_bits(), st.ue[1][ibl].to_bits());
        assert_eq!(bl.upper.stored.tau[k].to_bits(), st.tau[1][ibl].to_bits());
        assert_eq!(bl.upper.stored.tstr[k].to_bits(), st.thetastar[1][ibl].to_bits());
    }

    // M = 0: the Karman–Tsien transformation is the identity and Hk = H, bit for bit
    for side in [&bl.upper, &bl.lower, &bl.wake] {
        for k in 0..side.len() {
            assert_eq!(side.ue[k].to_bits(), (side.uedg[k] / bl.qinf).to_bits(), "ue at M=0");
            assert_eq!(side.hk[k].to_bits(), side.h[k].to_bits(), "hk = h at M=0");
            assert_eq!(side.msq[k], 0.0);
        }
    }
    // the wake has no wall shear
    assert!(bl.wake.cf.iter().all(|&c| c == 0.0));
    assert!(bl.wake.stored.cf_dump.iter().all(|&c| c == 0.0));

    // every plottable column is finite and sized
    for q in BlQuantity::ALL {
        for side in [&bl.upper, &bl.lower, &bl.wake] {
            let col = side.column(q);
            assert_eq!(col.len(), side.len());
            assert!(col.iter().all(|v| v.is_finite()), "{q} is finite");
        }
        assert!(bl.max_abs(q) > 0.0, "{q} has a scale");
    }
}

#[test]
fn live_closures_lag_the_stored_arrays_by_the_final_newton_correction() {
    let mut session = ref_case_session();
    let p = session.alpha(2.0_f64.to_radians());
    let out = AnalysisOutput::from_session(&session, &p, "naca0012", &spec(), false);
    let bl = out.boundary_layer.as_ref().unwrap();

    // The stored arrays were evaluated on the state entering the final Newton iteration; the live
    // columns on the state after its correction. RMSBL is the rms of the normalised corrections,
    // so the two agree to a modest multiple of it, and a wrong flow type (laminar closures on a
    // turbulent station) would show as an O(1) difference.
    let bound = 50.0 * p.residual.max(1e-6);
    let mut largest: f64 = 0.0;
    for side in [&bl.upper, &bl.lower, &bl.wake] {
        for k in 0..side.len() {
            let rel = |a: f64, b: f64| (a - b).abs() / a.abs().max(b.abs()).max(1e-3);
            largest = largest.max(rel(side.hs[k], side.stored.hs_dump[k]));
            largest = largest.max(rel(side.cf[k], side.stored.cf_dump[k]));
            largest = largest.max(rel(side.delta[k], side.stored.delt[k]));
            largest = largest.max(rel(side.ctq[k], side.stored.ctq[k]));
        }
    }
    assert!(
        largest > 0.0,
        "live and stored closures are not identical after a Newton correction"
    );
    assert!(
        largest < bound,
        "live vs stored closure mismatch {largest:.3e} exceeds {bound:.3e}"
    );
}

#[test]
fn markers_and_wake_split_follow_xfoil() {
    let mut session = ref_case_session();
    let p = session.alpha(2.0_f64.to_radians());
    let out = AnalysisOutput::from_session(&session, &p, "naca0012", &spec(), false);
    let bl = out.boundary_layer.as_ref().unwrap();
    let st = &session.state;

    // transition: XOCTR is what the operating point reports; the point lies on the surface
    assert_eq!(bl.transition[0].x_c.to_bits(), out.result.xtr_upper.to_bits());
    assert_eq!(bl.transition[1].x_c.to_bits(), out.result.xtr_lower.to_bits());
    assert_eq!(bl.transition[0].station, st.i_transition_station[1]);
    assert!(!bl.transition[0].forced);
    assert_eq!(bl.transition[0].s, st.s_stagnation - st.xi_transition[1]);
    assert_eq!(bl.transition[1].s, st.s_stagnation + st.xi_transition[2]);
    for t in &bl.transition {
        assert!(
            (t.x - t.x_c).abs() < 1e-3,
            "transition x on the surface ≈ x/c for a unit chord"
        );
    }

    // stagnation near the LE, on the spline
    assert_eq!(bl.stagnation.ist, st.i_stagnation_node);
    assert_eq!(bl.stagnation.sst, st.s_stagnation);
    assert!(bl.stagnation.x.abs() < 0.01);

    // CPDISP's split recomputes from the emitted columns; TESYS closes it at convergence
    let dstrte = bl.wake.dstr[0];
    let f1 = (bl.upper.dstr[bl.upper.len() - 1] + 0.5 * bl.ante) / dstrte;
    let f2 = (bl.lower.dstr[bl.lower.len() - 1] + 0.5 * bl.ante) / dstrte;
    assert_eq!(bl.wake_split[0].to_bits(), f1.to_bits());
    assert_eq!(bl.wake_split[1].to_bits(), f2.to_bits());
    // TESYS's DSTR(wake 1) = DSTR1 + DSTR2 + ANTE is a Newton equation, satisfied to the final
    // residual, not to the noise floor
    let closure = (bl.wake_split[0] + bl.wake_split[1] - 1.0).abs();
    assert!(
        closure <= p.residual.max(TOL_SOLVER),
        "DSF1 + DSF2 - 1 = {closure:.3e} exceeds the final RMSBL {:.3e}",
        p.residual
    );

    // derived separation: consistent with the cf sign pattern it was read from
    for m in &bl.derived_separation {
        let side = if m.side == 1 { &bl.upper } else { &bl.lower };
        let k = side.ibl.iter().position(|&i| i == m.station_before).unwrap();
        match m.kind {
            SeparationKind::Separation => assert!(side.cf[k] >= 0.0 && side.cf[k + 1] < 0.0),
            SeparationKind::Reattachment => assert!(side.cf[k] < 0.0 && side.cf[k + 1] >= 0.0),
        }
        assert!(m.s > st.s[side.node[k]].min(st.s[side.node[k + 1]]));
    }
}

#[test]
fn json_round_trips_and_inviscid_has_no_boundary_layer() {
    let mut session = ref_case_session();
    let p = session.alpha(2.0_f64.to_radians());
    let out = AnalysisOutput::from_session(&session, &p, "naca0012", &spec(), false);
    let json = out.to_json().unwrap();
    let back: AnalysisOutput = serde_json::from_str(&json).unwrap();
    assert_eq!(back.geometry, out.geometry);
    assert_eq!(back.boundary_layer, out.boundary_layer);

    let inviscid = FlowConditions { re: 0.0, ..spec() };
    let geometry = read_geometry_from_file(fixture_path("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let mut session = Session::new(&airfoil, inviscid.clone());
    let p = session.alpha(2.0_f64.to_radians());
    let out = AnalysisOutput::from_session(&session, &p, "naca0012", &inviscid, true);
    assert!(out.boundary_layer.is_none());
    assert!(out.geometry.wake.is_none());
    assert_eq!(out.geometry.n(), airfoil.n);
}

#[test]
fn polar_observer_sees_every_visited_point_in_its_own_state() {
    let geometry = read_geometry_from_file(fixture_path("panels.json").to_str().unwrap()).unwrap();
    let airfoil = panel_foil(&geometry);
    let config = PolarConfig {
        alpha_max: 2.0,
        alpha_min: -2.0,
        alpha_step: 1.0,
        conditions: spec(),
        ..Default::default()
    };
    let mut records: Vec<AnalysisOutput> = Vec::new();
    let mut visited: Vec<f64> = Vec::new();
    let result = compute_polar_with(&airfoil, &config, &mut |session, p| {
        visited.push(p.alpha.to_degrees());
        assert_eq!(session.state.alpha, p.alpha, "the session is in the point's state");
        records.push(AnalysisOutput::from_session(
            session,
            p,
            "naca0012",
            &config.conditions,
            false,
        ));
    });
    // 0, 1, 2 then (after INIT) -1, -2: the sweep order, every point once
    assert_eq!(visited, vec![0.0, 1.0, 2.0, -1.0, -2.0]);
    assert_eq!(records.len(), result.results.len() + result.failed_alphas.len());

    records.sort_by(|a, b| a.result.alpha_deg.partial_cmp(&b.result.alpha_deg).unwrap());
    let mut polar = PolarOutput::from_polar(&result, "naca0012");
    polar.distributions = records;
    let json = polar.to_json().unwrap();
    let back: PolarOutput = serde_json::from_str(&json).unwrap();
    assert_eq!(back.distributions.len(), 5);
    for d in &back.distributions {
        let pt = back
            .points
            .iter()
            .find(|q| q.alpha_deg == d.result.alpha_deg)
            .expect("every converged point has its distribution");
        assert_eq!(pt.cl.to_bits(), d.result.cl.to_bits());
        assert!(d.boundary_layer.is_some());
    }
    // a polar without distributions serialises without the field
    let plain = PolarOutput::from_polar(&result, "naca0012");
    assert!(!plain.to_json().unwrap().contains("\"distributions\""));
}
