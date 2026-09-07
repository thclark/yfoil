//! Velocity layer: QISET, UICALC, UECALC, QVFUE, GAMQV, UESET, DSSET (xpanel.f).
//! Line-for-line translations on the 1-based `SolverState`.

use crate::solver::blstate::SolverState;

/// QISET: inviscid panel tangential velocity for the current alpha from the alpha=0,90 solutions.
#[doc(alias = "QISET")]
pub fn set_q_inviscid(state: &mut SolverState, alfa: f64) {
    let cosa = alfa.cos();
    let sina = alfa.sin();
    for i in 1..=(state.n_foil_nodes + state.n_wake_nodes) {
        state.q_inviscid[i] = cosa * state.q_inviscid_basis[1][i] + sina * state.q_inviscid_basis[2][i];
        state.q_inviscid_d_alpha[i] = -sina * state.q_inviscid_basis[1][i] + cosa * state.q_inviscid_basis[2][i];
    }
}

/// UICALC: inviscid Ue from panel inviscid tangential velocity.
#[doc(alias = "UICALC")]
pub fn set_ue_inviscid(state: &mut SolverState) {
    for side in 1..=2 {
        state.ue_inviscid[side][1] = 0.0;
        state.ue_inviscid_d_alpha[side][1] = 0.0;
        for i_station in 2..=state.n_stations[side] {
            let i = state.i_node[side][i_station];
            state.ue_inviscid[side][i_station] = state.velocity_sign[side][i_station] * state.q_inviscid[i];
            state.ue_inviscid_d_alpha[side][i_station] =
                state.velocity_sign[side][i_station] * state.q_inviscid_d_alpha[i];
        }
    }
}

/// UECALC: viscous Ue from panel viscous tangential velocity.
#[doc(alias = "UECALC")]
pub fn set_ue_from_q_viscous(state: &mut SolverState) {
    for side in 1..=2 {
        state.ue[side][1] = 0.0;
        for i_station in 2..=state.n_stations[side] {
            let i = state.i_node[side][i_station];
            state.ue[side][i_station] = state.velocity_sign[side][i_station] * state.q_viscous[i];
        }
    }
}

/// QVFUE: panel viscous tangential velocity from viscous Ue.
pub fn set_q_viscous_from_ue(state: &mut SolverState) {
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            let i = state.i_node[side][i_station];
            state.q_viscous[i] = state.velocity_sign[side][i_station] * state.ue[side][i_station];
        }
    }
}

/// GAMQV: GAM from QVIS (airfoil nodes only), GAM_A from QINV_A.
#[doc(alias = "GAMQV")]
pub fn set_gamma_from_q_viscous(state: &mut SolverState) {
    for i in 1..=state.n_foil_nodes {
        state.gamma[i] = state.q_viscous[i];
        state.gamma_d_alpha[i] = state.q_inviscid_d_alpha[i];
    }
}

/// UESET: Ue from inviscid Ue plus all source (mass defect) influence through `st.dij`.
pub fn set_ue_with_sources(state: &mut SolverState) {
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            let i = state.i_node[side][i_station];
            let mut dui = 0.0;
            for j_side in 1..=2 {
                for j_station in 2..=state.n_stations[j_side] {
                    let j = state.i_node[j_side][j_station];
                    let ue_m = -state.velocity_sign[side][i_station]
                        * state.velocity_sign[j_side][j_station]
                        * state.dij[i][j];
                    dui += ue_m * state.mass_defect[j_side][j_station];
                }
            }
            state.ue[side][i_station] = state.ue_inviscid[side][i_station] + dui;
        }
    }
}

/// DSSET: displacement thickness from mass defect and Ue.
#[doc(alias = "DSSET")]
pub fn set_dstar_from_mass(state: &mut SolverState) {
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            state.dstar[side][i_station] = state.mass_defect[side][i_station] / state.ue[side][i_station];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::pointers::{map_stations_to_nodes, map_stations_to_rows};

    fn small_state() -> SolverState {
        let (n, nw) = (12, 4);
        let mut state = SolverState::empty(n, nw);
        state.i_stagnation_node = 6;
        map_stations_to_nodes(&mut state);
        map_stations_to_rows(&mut state);
        state
    }

    /// QVFUE then UECALC is the identity on UEDG (VTI*VTI = 1), and GAMQV copies QVIS into GAM.
    #[test]
    fn qvfue_uecalc_round_trip_is_identity() {
        let mut state = small_state();
        for side in 1..=2 {
            for i_station in 2..=state.n_stations[side] {
                state.ue[side][i_station] = 0.1 * i_station as f64 + side as f64;
            }
        }
        let before = state.ue.clone();
        set_q_viscous_from_ue(&mut state);
        set_gamma_from_q_viscous(&mut state);
        set_ue_from_q_viscous(&mut state);
        for side in 1..=2 {
            for i_station in 2..=state.n_stations[side] {
                assert_eq!(state.ue[side][i_station].to_bits(), before[side][i_station].to_bits());
                let i = state.i_node[side][i_station];
                if i <= state.n_foil_nodes {
                    assert_eq!(state.gamma[i].to_bits(), state.q_viscous[i].to_bits());
                }
            }
        }
    }

    /// With DIJ = 0, UESET reduces to UEDG = UINV; with MASS = 0 likewise; DSSET inverts MASS.
    #[test]
    fn ueset_and_dsset_degenerate_cases() {
        let mut state = small_state();
        let np = state.n_foil_nodes + state.n_wake_nodes;
        state.dij = vec![vec![0.0; np + 1]; np + 1];
        for side in 1..=2 {
            for i_station in 2..=state.n_stations[side] {
                state.ue_inviscid[side][i_station] = 1.0 + 0.01 * i_station as f64;
                state.mass_defect[side][i_station] = 0.5 * i_station as f64;
            }
        }
        set_ue_with_sources(&mut state);
        set_dstar_from_mass(&mut state);
        for side in 1..=2 {
            for i_station in 2..=state.n_stations[side] {
                assert_eq!(
                    state.ue[side][i_station].to_bits(),
                    state.ue_inviscid[side][i_station].to_bits()
                );
                assert_eq!(
                    state.dstar[side][i_station].to_bits(),
                    (state.mass_defect[side][i_station] / state.ue[side][i_station]).to_bits()
                );
            }
        }
    }
}
