//! Velocity layer: QISET, UICALC, UECALC, QVFUE, GAMQV, UESET, DSSET (xpanel.f).
//! Line-for-line translations on the 1-based `BlState`.

use crate::solver::blstate::BlState;

/// QISET: inviscid panel tangential velocity for the current alpha from the alpha=0,90 solutions.
pub fn qiset(st: &mut BlState, alfa: f64) {
    let cosa = alfa.cos();
    let sina = alfa.sin();
    for i in 1..=(st.n + st.nw) {
        st.qinv[i] = cosa * st.qinvu[1][i] + sina * st.qinvu[2][i];
        st.qinv_a[i] = -sina * st.qinvu[1][i] + cosa * st.qinvu[2][i];
    }
}

/// UICALC: inviscid Ue from panel inviscid tangential velocity.
pub fn uicalc(st: &mut BlState) {
    for is in 1..=2 {
        st.uinv[is][1] = 0.0;
        st.uinv_a[is][1] = 0.0;
        for ibl in 2..=st.nbl[is] {
            let i = st.ipan[is][ibl];
            st.uinv[is][ibl] = st.vti[is][ibl] * st.qinv[i];
            st.uinv_a[is][ibl] = st.vti[is][ibl] * st.qinv_a[i];
        }
    }
}

/// UECALC: viscous Ue from panel viscous tangential velocity.
pub fn uecalc(st: &mut BlState) {
    for is in 1..=2 {
        st.uedg[is][1] = 0.0;
        for ibl in 2..=st.nbl[is] {
            let i = st.ipan[is][ibl];
            st.uedg[is][ibl] = st.vti[is][ibl] * st.qvis[i];
        }
    }
}

/// QVFUE: panel viscous tangential velocity from viscous Ue.
pub fn qvfue(st: &mut BlState) {
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            let i = st.ipan[is][ibl];
            st.qvis[i] = st.vti[is][ibl] * st.uedg[is][ibl];
        }
    }
}

/// GAMQV: GAM from QVIS (airfoil nodes only), GAM_A from QINV_A.
pub fn gamqv(st: &mut BlState) {
    for i in 1..=st.n {
        st.gam[i] = st.qvis[i];
        st.gam_a[i] = st.qinv_a[i];
    }
}

/// UESET: Ue from inviscid Ue plus all source (mass defect) influence through `st.dij`.
pub fn ueset(st: &mut BlState) {
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            let i = st.ipan[is][ibl];
            let mut dui = 0.0;
            for js in 1..=2 {
                for jbl in 2..=st.nbl[js] {
                    let j = st.ipan[js][jbl];
                    let ue_m = -st.vti[is][ibl] * st.vti[js][jbl] * st.dij[i][j];
                    dui += ue_m * st.mass[js][jbl];
                }
            }
            st.uedg[is][ibl] = st.uinv[is][ibl] + dui;
        }
    }
}

/// DSSET: displacement thickness from mass defect and Ue.
pub fn dsset(st: &mut BlState) {
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            st.dstr[is][ibl] = st.mass[is][ibl] / st.uedg[is][ibl];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::pointers::{iblpan, iblsys};

    fn small_state() -> BlState {
        let (n, nw) = (12, 4);
        let mut st = BlState::empty(n, nw);
        st.ist = 6;
        iblpan(&mut st);
        iblsys(&mut st);
        st
    }

    /// QVFUE then UECALC is the identity on UEDG (VTI*VTI = 1), and GAMQV copies QVIS into GAM.
    #[test]
    fn qvfue_uecalc_round_trip_is_identity() {
        let mut st = small_state();
        for is in 1..=2 {
            for ibl in 2..=st.nbl[is] {
                st.uedg[is][ibl] = 0.1 * ibl as f64 + is as f64;
            }
        }
        let before = st.uedg.clone();
        qvfue(&mut st);
        gamqv(&mut st);
        uecalc(&mut st);
        for is in 1..=2 {
            for ibl in 2..=st.nbl[is] {
                assert_eq!(st.uedg[is][ibl].to_bits(), before[is][ibl].to_bits());
                let i = st.ipan[is][ibl];
                if i <= st.n {
                    assert_eq!(st.gam[i].to_bits(), st.qvis[i].to_bits());
                }
            }
        }
    }

    /// With DIJ = 0, UESET reduces to UEDG = UINV; with MASS = 0 likewise; DSSET inverts MASS.
    #[test]
    fn ueset_and_dsset_degenerate_cases() {
        let mut st = small_state();
        let np = st.n + st.nw;
        st.dij = vec![vec![0.0; np + 1]; np + 1];
        for is in 1..=2 {
            for ibl in 2..=st.nbl[is] {
                st.uinv[is][ibl] = 1.0 + 0.01 * ibl as f64;
                st.mass[is][ibl] = 0.5 * ibl as f64;
            }
        }
        ueset(&mut st);
        dsset(&mut st);
        for is in 1..=2 {
            for ibl in 2..=st.nbl[is] {
                assert_eq!(st.uedg[is][ibl].to_bits(), st.uinv[is][ibl].to_bits());
                assert_eq!(
                    st.dstr[is][ibl].to_bits(),
                    (st.mass[is][ibl] / st.uedg[is][ibl]).to_bits()
                );
            }
        }
    }
}
