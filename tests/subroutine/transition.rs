//! (a) Subroutine equivalence: TRCHEK2's N2 Newton, iterate by iterate, from the inputs XFOIL
//! handed it. For every dumped SETBL call the instrumented reference writes each TRCHEK2 call's
//! interval inputs and every iterate of its Newton just before the convergence test
//! (`trchek2_<k>.dat`: CALL X1 X2 T1 T2 D1 D2 U1 U2 AMPL1 AMPL2 AMCRIT XIFORC; GUESS; ITAM n AMPL2 XT
//! AX AX_A2 RES RES_A2 DA2 RLX; DONE n AMPL2 XT). yFoil's `check_transition_traced` runs on the
//! same inputs; the iterate count must be equal, the initial guess within `TOL_SOLVER` and the
//! final XT within `TOL_PURE`. The calls are the upper-side transition intervals of the eight
//! iterations of the 7° point of the NACA 0012 polar (known issues §7.9).
//!
//! Fixtures: `tests/fixtures/xfoil/naca0012_n60_polar30_re1e6/` (`trchek2_<k>.dat`,
//! `mrchdu_output_<k>.dat`, k = 40–47) — `cargo xtask fixtures --case naca0012_n60_polar30_re1e6`.

use crate::fixtures::cases::load;
use crate::fixtures::mrchdu_fixtures::parse_bl_dump;
use crate::utilities::tolerances::{assert_within, TOL_PURE, TOL_SOLVER};
use yfoil::bl::params::FlowParameters;
use yfoil::bl::station::StationState;
use yfoil::bl::transition::{check_transition_traced, TransitionIterate};

const CASE: &str = "naca0012_n60_polar30_re1e6";
/// the 7° point
const CALL: usize = 8;

/// One TRCHEK2 call of the reference: the interval's inputs, the iterates, and where it stopped.
struct TrchekCall {
    /// X1 X2 T1 T2 D1 D2 U1 U2 AMPL1 AMPL2 (as the march handed it over) AMCRIT XIFORC
    inputs: [f64; 12],
    /// the initial guess AMPL1 + AX·(X2 − X1)
    guess: f64,
    iterates: Vec<TransitionIterate>,
    done: (usize, f64, f64),
}

fn parse_trchek2(path: &std::path::Path) -> Vec<TrchekCall> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| panic!("{}", path.display()));
    let mut calls: Vec<TrchekCall> = vec![];
    for l in text.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        let f = |i: usize| -> f64 { t[i].parse().unwrap() };
        match t[0] {
            "CALL" => calls.push(TrchekCall {
                inputs: std::array::from_fn(|i| f(i + 1)),
                guess: f64::NAN,
                iterates: vec![],
                done: (0, 0.0, 0.0),
            }),
            "GUESS" => calls.last_mut().unwrap().guess = f(1),
            "ITAM" => calls.last_mut().unwrap().iterates.push(TransitionIterate {
                itam: t[1].parse().unwrap(),
                ampl2: f(2),
                xt: f(3),
                ax: f(4),
                ax_a2: f(5),
                res: f(6),
                res_a2: f(7),
                da2: f(8),
                rlx: f(9),
            }),
            "DONE" => calls.last_mut().unwrap().done = (t[1].parse().unwrap(), f(2), f(3)),
            _ => {}
        }
    }
    calls
}

#[test]
fn test_trchek2_newton_path_matches_xfoil_at_the_7deg_point() {
    let c = load(CASE);
    let nit: Vec<usize> = c.points.iter().map(|p| p["NITDONE"].parse().unwrap()).collect();
    let before: usize = nit[..CALL - 1].iter().sum();
    let mut params = FlowParameters::new(c.spec.mach, c.spec.re.unwrap(), 1.4);
    params.amplification_model = c.spec.amplification_model;
    for i in 0..nit[CALL - 1] {
        let k = before + i + 1;
        let calls = parse_trchek2(&c.dir.join(format!("trchek2_{k}.dat")));
        let m = parse_bl_dump(&c.dir.join(format!("mrchdu_output_{k}.dat")));
        let x_xt: f64 = m.header["XSSITR1"].trim().parse().unwrap();
        // the call whose result is XSSITR(1): the transition interval on the upper side
        let call = calls
            .iter()
            .find(|t| t.done.2 == x_xt)
            .unwrap_or_else(|| panic!("{CASE} SETBL {k}: no TRCHEK2 call ends with XSSITR(1) = {x_xt:.17e}"));
        let [x1, x2, t1, t2, d1, d2, u1, u2, ampl1, ampl2, amcrit, xiforc] = call.inputs;
        let mut s1 = StationState::default();
        let mut s2 = StationState::default();
        // laminar stations: CTAU carries the amplification; no wake gap on the airfoil
        s1.set_primary_variables(x1, ampl1, ampl1, t1, d1, 0.0, u1, &params);
        s1.set_kinematic_variables(&params);
        s2.set_primary_variables(x2, ampl2, ampl2, t2, d2, 0.0, u2, &params);
        s2.set_kinematic_variables(&params);
        let mut trace = vec![];
        let _ = check_transition_traced(&s1, &s2, ampl1, amcrit, xiforc, &params, Some(&mut trace));
        let what = format!("{CASE} SETBL {k}");
        assert_within(
            trace[0].ampl2,
            call.guess,
            TOL_SOLVER,
            1.0,
            &format!("{what}: initial AMPL2"),
        );
        let ours: Vec<&TransitionIterate> = trace.iter().filter(|t| t.itam > 0).collect();
        assert_eq!(ours.len(), call.iterates.len(), "{what}: N2 Newton iterate count");
        assert_within(
            ours.last().unwrap().xt,
            call.done.2,
            TOL_PURE,
            1.0,
            &format!("{what}: final XT"),
        );
    }
}
