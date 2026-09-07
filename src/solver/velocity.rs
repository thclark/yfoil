//! Velocity layer: QISET, UICALC, UECALC, QVFUE, GAMQV, UESET, DSSET (xpanel.f).
//! Line-for-line translations on the 1-based `SolverState`.

use crate::solver::blstate::SolverState;

/// QISET: inviscid panel tangential velocity for the current alpha from the alpha=0,90 solutions.
#[doc(alias = "QISET")]
pub fn set_q_inviscid(st: &mut SolverState, alfa: f64) {
    let cosa = alfa.cos();
    let sina = alfa.sin();
    for i in 1..=(st.n_foil_nodes + st.n_wake_nodes) {
        st.q_inviscid[i] = cosa * st.q_inviscid_basis[1][i] + sina * st.q_inviscid_basis[2][i];
        st.q_inviscid_d_alpha[i] = -sina * st.q_inviscid_basis[1][i] + cosa * st.q_inviscid_basis[2][i];
    }
}

/// UICALC: inviscid Ue from panel inviscid tangential velocity.
#[doc(alias = "UICALC")]
pub fn set_ue_inviscid(st: &mut SolverState) {
    for is in 1..=2 {
        st.ue_inviscid[is][1] = 0.0;
        st.ue_inviscid_d_alpha[is][1] = 0.0;
        for ibl in 2..=st.n_stations[is] {
            let i = st.i_node[is][ibl];
            st.ue_inviscid[is][ibl] = st.velocity_sign[is][ibl] * st.q_inviscid[i];
            st.ue_inviscid_d_alpha[is][ibl] = st.velocity_sign[is][ibl] * st.q_inviscid_d_alpha[i];
        }
    }
}

/// UECALC: viscous Ue from panel viscous tangential velocity.
#[doc(alias = "UECALC")]
pub fn set_ue_from_q_viscous(st: &mut SolverState) {
    for is in 1..=2 {
        st.ue[is][1] = 0.0;
        for ibl in 2..=st.n_stations[is] {
            let i = st.i_node[is][ibl];
            st.ue[is][ibl] = st.velocity_sign[is][ibl] * st.q_viscous[i];
        }
    }
}

/// QVFUE: panel viscous tangential velocity from viscous Ue.
pub fn set_q_viscous_from_ue(st: &mut SolverState) {
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            let i = st.i_node[is][ibl];
            st.q_viscous[i] = st.velocity_sign[is][ibl] * st.ue[is][ibl];
        }
    }
}

/// GAMQV: GAM from QVIS (airfoil nodes only), GAM_A from QINV_A.
#[doc(alias = "GAMQV")]
pub fn set_gamma_from_q_viscous(st: &mut SolverState) {
    for i in 1..=st.n_foil_nodes {
        st.gamma[i] = st.q_viscous[i];
        st.gamma_d_alpha[i] = st.q_inviscid_d_alpha[i];
    }
}

/// UESET: Ue from inviscid Ue plus all source (mass defect) influence through `st.dij`.
pub fn set_ue_with_sources(st: &mut SolverState) {
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            let i = st.i_node[is][ibl];
            let mut dui = 0.0;
            for js in 1..=2 {
                for jbl in 2..=st.n_stations[js] {
                    let j = st.i_node[js][jbl];
                    let ue_m = -st.velocity_sign[is][ibl] * st.velocity_sign[js][jbl] * st.dij[i][j];
                    dui += ue_m * st.mass_defect[js][jbl];
                }
            }
            st.ue[is][ibl] = st.ue_inviscid[is][ibl] + dui;
        }
    }
}

/// DSSET: displacement thickness from mass defect and Ue.
#[doc(alias = "DSSET")]
pub fn set_dstar_from_mass(st: &mut SolverState) {
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            st.dstar[is][ibl] = st.mass_defect[is][ibl] / st.ue[is][ibl];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::pointers::{map_stations_to_nodes, map_stations_to_rows};

    fn small_state() -> SolverState {
        let (n, nw) = (12, 4);
        let mut st = SolverState::empty(n, nw);
        st.i_stagnation_node = 6;
        map_stations_to_nodes(&mut st);
        map_stations_to_rows(&mut st);
        st
    }

    /// QVFUE then UECALC is the identity on UEDG (VTI*VTI = 1), and GAMQV copies QVIS into GAM.
    #[test]
    fn qvfue_uecalc_round_trip_is_identity() {
        let mut st = small_state();
        for is in 1..=2 {
            for ibl in 2..=st.n_stations[is] {
                st.ue[is][ibl] = 0.1 * ibl as f64 + is as f64;
            }
        }
        let before = st.ue.clone();
        set_q_viscous_from_ue(&mut st);
        set_gamma_from_q_viscous(&mut st);
        set_ue_from_q_viscous(&mut st);
        for is in 1..=2 {
            for ibl in 2..=st.n_stations[is] {
                assert_eq!(st.ue[is][ibl].to_bits(), before[is][ibl].to_bits());
                let i = st.i_node[is][ibl];
                if i <= st.n_foil_nodes {
                    assert_eq!(st.gamma[i].to_bits(), st.q_viscous[i].to_bits());
                }
            }
        }
    }

    /// With DIJ = 0, UESET reduces to UEDG = UINV; with MASS = 0 likewise; DSSET inverts MASS.
    #[test]
    fn ueset_and_dsset_degenerate_cases() {
        let mut st = small_state();
        let np = st.n_foil_nodes + st.n_wake_nodes;
        st.dij = vec![vec![0.0; np + 1]; np + 1];
        for is in 1..=2 {
            for ibl in 2..=st.n_stations[is] {
                st.ue_inviscid[is][ibl] = 1.0 + 0.01 * ibl as f64;
                st.mass_defect[is][ibl] = 0.5 * ibl as f64;
            }
        }
        set_ue_with_sources(&mut st);
        set_dstar_from_mass(&mut st);
        for is in 1..=2 {
            for ibl in 2..=st.n_stations[is] {
                assert_eq!(st.ue[is][ibl].to_bits(), st.ue_inviscid[is][ibl].to_bits());
                assert_eq!(
                    st.dstar[is][ibl].to_bits(),
                    (st.mass_defect[is][ibl] / st.ue[is][ibl]).to_bits()
                );
            }
        }
    }
}
