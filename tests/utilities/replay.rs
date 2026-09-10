//! The one-step replay harness (CLAUDE.md Rule 3): seed yFoil with XFOIL's exact state entering
//! SETBL call `k` (`mrchdu_input_<k>.dat`), run one VISCAL iteration, and compare with XFOIL's
//! state after UPDATE `k` (`update_output_<k>.dat`). This isolates one iteration from everything
//! accumulated before it: if the one-step map agrees, a divergence of the full runs is
//! accumulated differences crossing a threshold (threshold-straddling); if not, the step itself
//! is wrong.

#![allow(dead_code)]

use super::records::{allowed, load};
use super::tolerances::{assert_within, TOL_SOLVER};
use std::collections::HashMap;
use std::path::Path;
use yfoil::bl::blsolv::solve_newton_system;
use yfoil::solver::analysis::Session;
use yfoil::solver::clcalc::set_compressibility;
use yfoil::solver::pointers::{map_stations_to_nodes, map_stations_to_rows, set_station_xi};
use yfoil::solver::setbl::assemble_newton_system;
use yfoil::solver::update::apply_newton_update;
use yfoil::solver::velocity::set_ue_inviscid;

/// A parsed BL state dump (`fixtures::mrchdu_fixtures::BlDump`), borrowed so this shared module
/// needs no `fixtures` dependency: the header key/values and the per-side station rows
/// `XSSI UEDG THET DSTR CTAU MASS`.
pub struct DumpView<'a> {
    pub header: &'a HashMap<String, String>,
    pub bl: &'a [Vec<[f64; 6]>; 3],
}

impl DumpView<'_> {
    pub fn int(&self, key: &str) -> usize {
        self.header[key].parse().unwrap_or_else(|_| panic!("header `{key}`"))
    }
    pub fn real(&self, key: &str) -> f64 {
        self.header[key].parse().unwrap_or_else(|_| panic!("header `{key}`"))
    }
    /// NBL(is): the number of station rows (rows are 1-based with a dummy slot 0)
    pub fn nbl(&self, is: usize) -> usize {
        self.bl[is].len() - 1
    }
}

/// Replay SETBL/UPDATE call `k`, which is iteration `iteration` of VISCAL call `call`, from
/// XFOIL's dumped state `d` (`mrchdu_input_<k>.dat`) against `o` (`update_output_<k>.dat`),
/// with the case's records and noise floor read from `dir`. `prepare` must bring `session` to the state in which XFOIL
/// entered that VISCAL call with the prologue done (SPECAL, then `solve_viscous` with zero
/// iterations: wake, QINV, pointers, UINV, DIJ) — everything the dump does not carry.
pub fn replay_iteration(
    case: &str,
    dir: &Path,
    d: &DumpView,
    o: &DumpView,
    mut session: Session,
    prepare: impl FnOnce(&mut Session),
    call: usize,
    iteration: usize,
    k: usize,
) {
    prepare(&mut session);
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
    let rec = load(dir);
    let fl = |name: &str| {
        rec.floor
            .as_ref()
            .and_then(|f| f.calls.get(call - 1))
            .and_then(|c| c.iterations.get(iteration - 1))
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
        let lim = allowed(a, b, TOL_SOLVER, 1.0, fl(name));
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
    // Side 1's rows beyond NBL(1) are the "upper wake arrays" UPDATE equates to side 2's wake
    // (xbl.f:1546-1557) — CTAU, THET, DSTR, UEDG and the closures, never MASS. MASS(IBLTE(1)+k, 1)
    // is written by nothing and read by nothing (SETBL's mass-defect sums run to NBL(1)), so it
    // holds whatever each code last left there and is not compared.
    let nbl1 = st.n_stations[1];
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
                if is == 1 && ibl > nbl1 && *name == "MASS" {
                    continue;
                }
                let b = o.bl[is][ibl][m + 1];
                let scale = (2..=nrows).map(|j| o.bl[is][j][m + 1].abs()).fold(0.0_f64, f64::max);
                let e = (ours[m] - b).abs() / ours[m].abs().max(b.abs()).max(scale);
                if e > worst.0 {
                    worst = (e, name, is, ibl);
                }
                let lim = allowed(ours[m], b, TOL_SOLVER, scale, arr_floor(name));
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
        "{case} replay of SETBL call {k} (VISCAL call {call} iteration {iteration}) from XFOIL's state: RLX {:.6} RMSBL {:.6e} {}@({},{}) CL {:.8} — every array within max({TOL_SOLVER:.0e}·scale, 4× its own 1-ULP floor) (worst rel {:.2e} {} at ({},{}))",
        u.relaxation, u.residual, u.residual_max_variable, u.i_residual_max_station, u.residual_max_side, st.cl, worst.0, worst.1, worst.3, worst.2
    );
}
