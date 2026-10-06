//! (a) Subroutine equivalence: the BL closure functions — HKIN, CFL, HSL, DIL, HST, CFT, DAMPL,
//! AXSET, BLKIN, BLVAR — against every distinct call the reference logged, inputs to outputs.
//!
//! Each output is compared within `TOL_PURE` by the Rule 1 metric, on the scale of the largest
//! magnitude that output takes in the fixture set (several are near zero at some inputs: a
//! derivative with respect to MSQ at M = 0, the amplification rate below critical). BLVAR's log
//! records its values and not their derivatives, so the values are what is gated there.
//!
//! Fixtures: `tests/fixtures/subroutines/<closure>/<case>_NNN.json` and `manifest.json` — `cargo
//! xtask fixtures --case closures_naca0012_n60_re1e6 --case closures_naca0012_n60_re1e6_m03`
//! (`xtask/src/closures.rs`).

use crate::utilities::tolerances::{assert_within, TOL_PURE};
use serde_json::Value;
use std::path::Path;
use yfoil::bl::system::{
    amplification_rate, interval_amplification_rate, AmplificationModel, FlowParameters, FlowRegime, StationState,
};
use yfoil::bl::{cdiss_laminar, cf_laminar, cf_turbulent, hk_from_h, hstar_laminar, hstar_turbulent};

/// Every fixture of one closure: (file name, input, output). Fails if there are none, or if any
/// file does not parse.
fn cases(closure: &str) -> Vec<(String, Value, Value)> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/subroutines")
        .join(closure);
    let mut out = vec![];
    for e in std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .flatten()
    {
        let name = e.file_name().to_string_lossy().to_string();
        let text = std::fs::read_to_string(e.path()).unwrap();
        let v: Value = serde_json::from_str(&text).unwrap_or_else(|err| panic!("{name}: {err}"));
        out.push((name, v["input"].clone(), v["output"].clone()));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    assert!(!out.is_empty(), "no {closure} fixtures in {}", dir.display());
    out
}

fn f(v: &Value, key: &str) -> f64 {
    v[key].as_f64().unwrap_or_else(|| panic!("missing `{key}`"))
}

/// Run `eval` on every fixture of `closure` and compare each output it returns with XFOIL's.
fn check(closure: &str, eval: impl Fn(&Value) -> Vec<(&'static str, f64)>) {
    let all = cases(closure);
    let results: Vec<Vec<(&str, f64)>> = all.iter().map(|(_, i, _)| eval(i)).collect();
    for (k, (name, _)) in results[0].iter().enumerate() {
        // an output identically zero over the set (CF_MSQ of CFL, HS_RT of HSL) is compared on 1
        let scale = all.iter().map(|(_, _, o)| f(o, name).abs()).fold(0.0_f64, f64::max);
        let scale = if scale > 0.0 { scale } else { 1.0 };
        for ((file, _, o), ours) in all.iter().zip(&results) {
            assert_within(
                ours[k].1,
                f(o, name),
                TOL_PURE,
                scale,
                &format!("{closure} {file}: {name}"),
            );
        }
    }
}

#[test]
fn test_hkin_matches_xfoil() {
    check("hkin", |i| {
        let (hk, hk_h, hk_msq) = hk_from_h(f(i, "h"), f(i, "msq"));
        vec![("hk", hk), ("hk_h", hk_h), ("hk_msq", hk_msq)]
    });
}

#[test]
fn test_cfl_matches_xfoil() {
    check("cfl", |i| {
        let r = cf_laminar(f(i, "hk"), f(i, "rt"), f(i, "msq"));
        vec![
            ("cf", r.value),
            ("cf_hk", r.value_d_hk),
            ("cf_rt", r.value_d_retheta),
            ("cf_msq", r.value_d_machsqd),
        ]
    });
}

#[test]
fn test_hsl_matches_xfoil() {
    check("hsl", |i| {
        let r = hstar_laminar(f(i, "hk"), f(i, "rt"), f(i, "msq"));
        vec![
            ("hs", r.value),
            ("hs_hk", r.value_d_hk),
            ("hs_rt", r.value_d_retheta),
            ("hs_msq", r.value_d_machsqd),
        ]
    });
}

#[test]
fn test_dil_matches_xfoil() {
    check("dil", |i| {
        let r = cdiss_laminar(f(i, "hk"), f(i, "rt"));
        vec![("di", r.value), ("di_hk", r.value_d_hk), ("di_rt", r.value_d_retheta)]
    });
}

#[test]
fn test_hst_matches_xfoil() {
    check("hst", |i| {
        let r = hstar_turbulent(f(i, "hk"), f(i, "rt"), f(i, "msq"));
        vec![
            ("hs", r.value),
            ("hs_hk", r.value_d_hk),
            ("hs_rt", r.value_d_retheta),
            ("hs_msq", r.value_d_machsqd),
        ]
    });
}

#[test]
fn test_cft_matches_xfoil() {
    check("cft", |i| {
        let r = cf_turbulent(f(i, "hk"), f(i, "rt"), f(i, "msq"), 1.0);
        vec![
            ("cf", r.value),
            ("cf_hk", r.value_d_hk),
            ("cf_rt", r.value_d_retheta),
            ("cf_msq", r.value_d_machsqd),
        ]
    });
}

#[test]
fn test_dampl_matches_xfoil() {
    check("dampl", |i| {
        let r = amplification_rate(f(i, "hk"), f(i, "th"), f(i, "rt"));
        vec![
            ("ax", r.rate),
            ("ax_hk", r.rate_d_hk),
            ("ax_th", r.rate_d_theta),
            ("ax_rt", r.rate_d_retheta),
        ]
    });
}

#[test]
fn test_axset_matches_xfoil() {
    check("axset", |i| {
        let r = interval_amplification_rate(
            f(i, "hk1"),
            f(i, "t1"),
            f(i, "rt1"),
            f(i, "a1"),
            f(i, "hk2"),
            f(i, "t2"),
            f(i, "rt2"),
            f(i, "a2"),
            f(i, "acrit"),
            AmplificationModel::Envelope,
        );
        vec![
            ("ax", r.rate),
            ("ax_hk1", r.rate_d_hk_station1),
            ("ax_t1", r.rate_d_theta_station1),
            ("ax_rt1", r.rate_d_retheta_station1),
            ("ax_a1", r.rate_d_ampl_station1),
            ("ax_hk2", r.rate_d_hk_station2),
            ("ax_t2", r.rate_d_theta_station2),
            ("ax_rt2", r.rate_d_retheta_station2),
            ("ax_a2", r.rate_d_ampl_station2),
        ]
    });
}

#[test]
fn test_blkin_matches_xfoil() {
    check("blkin", |i| {
        let mut st = StationState {
            theta: f(i, "t2"),
            dstar: f(i, "d2"),
            ue: f(i, "u2"),
            ..StationState::default()
        };
        let mut params = FlowParameters::new(0.0, f(i, "reybl"), 1.0 + f(i, "gm1bl"));
        params.h_stagnation_inv = f(i, "hstinv");
        params.gamma_gas_m1 = f(i, "gm1bl");
        params.rho_stagnation = f(i, "rstbl");
        params.sutherland_ratio = f(i, "hvrat");
        params.re = f(i, "reybl");
        st.set_kinematic_variables(&params);
        vec![
            ("m2", st.machsqd_edge),
            ("h2", st.h),
            ("hk2", st.hk),
            ("rt2", st.retheta),
            ("hk2_t2", st.hk_d_theta),
            ("hk2_d2", st.hk_d_dstar),
            ("hk2_u2", st.hk_d_ue),
            ("rt2_t2", st.retheta_d_theta),
            ("rt2_u2", st.retheta_d_ue),
        ]
    });
}

#[test]
fn test_blvar_matches_xfoil() {
    check("blvar", |i| {
        let mut st = StationState {
            hk: f(i, "hk2"),
            retheta: f(i, "rt2"),
            machsqd_edge: f(i, "m2"),
            theta: f(i, "t2"),
            dstar: f(i, "d2"),
            sqrtctau: f(i, "s2"),
            ..StationState::default()
        };
        st.h = st.dstar / st.theta;
        let regime = match i["ityp"].as_i64().unwrap() {
            1 => FlowRegime::Laminar,
            2 => FlowRegime::Turbulent,
            3 => FlowRegime::Wake,
            t => panic!("ITYP {t}"),
        };
        st.set_closure_variables(regime, &FlowParameters::new(0.0, 1e6, 1.4));
        vec![
            ("hs2", st.hstar),
            ("cf2", st.cf),
            ("di2", st.cdiss),
            ("us2", st.us),
            ("cq2", st.sqrtctaueq),
            ("de2", st.delta),
        ]
    });
}

/// The fixtures' provenance: every file belongs to a case the manifest records, with the reference
/// build that produced it.
#[test]
fn test_closure_fixtures_have_a_manifest() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/subroutines");
    let m: Value = serde_json::from_str(&std::fs::read_to_string(root.join("manifest.json")).unwrap()).unwrap();
    let cases = m["cases"].as_object().expect("manifest cases");
    for closure in [
        "hkin", "cfl", "hsl", "dil", "hst", "cft", "dampl", "axset", "blkin", "blvar",
    ] {
        for (file, _, _) in self::cases(closure) {
            assert!(
                cases.keys().any(|c| file.starts_with(&format!("{c}_"))),
                "{closure}/{file}: no case in manifest.json"
            );
        }
    }
    assert!(cases
        .values()
        .all(|c| c["xfoil_ref"].as_array().is_some_and(|a| !a.is_empty())));
}
