//! BL pointer layer: STFIND, IBLPAN, XICALC, IBLSYS, XIFSET (xpanel.f / xbl.f), plus TECALC.
//!
//! Each function is a line-for-line translation; the Fortran statement is quoted where the
//! translation is not obvious. Indexing is 1-based (see `blstate.rs`).

use crate::geometry::{deval, seval, spline};
use crate::solver::blstate::SolverState;

/// TECALC (xfoil.f): TE gap areas and the sharp-TE flag.
pub fn tecalc(st: &mut SolverState) {
    let n = st.n_foil_nodes;
    // set TE base vector and TE bisector components
    let dxte = st.x[1] - st.x[n];
    let dyte = st.y[1] - st.y[n];
    let dxs = 0.5 * (-st.dxds[1] + st.dxds[n]);
    let dys = 0.5 * (-st.dyds[1] + st.dyds[n]);
    // normal and streamwise projected TE gap areas
    st.te_thickness_normal = dxs * dyte - dys * dxte;
    st.te_thickness_parallel = dxs * dxte + dys * dyte;
    // total TE gap area
    st.te_gap = (dxte * dxte + dyte * dyte).sqrt();
    st.sharp_te = st.te_gap < 0.0001 * st.chord;
}

/// STFIND (xpanel.f): stagnation point arc length SST, panel index IST, and the
/// sensitivities SST_GO = dSST/dGAM(IST), SST_GP = dSST/dGAM(IST+1).
pub fn stfind(st: &mut SolverState) {
    let n = st.n_foil_nodes;
    let mut i = n / 2; // fallback if no sign change is found ("Stagnation point not found")
    for ii in 1..n {
        if st.gamma[ii] >= 0.0 && st.gamma[ii + 1] < 0.0 {
            i = ii;
            break;
        }
    }
    st.i_stagnation_node = i;
    let dgam = st.gamma[i + 1] - st.gamma[i];
    let ds = st.s[i + 1] - st.s[i];

    // evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
    let mut sst = if st.gamma[i] < -st.gamma[i + 1] {
        st.s[i] - ds * (st.gamma[i] / dgam)
    } else {
        st.s[i + 1] - ds * (st.gamma[i + 1] / dgam)
    };

    // tweak stagnation point if it falls right on a node (very unlikely)
    if sst <= st.s[i] {
        sst = st.s[i] + 1.0e-7;
    }
    if sst >= st.s[i + 1] {
        sst = st.s[i + 1] - 1.0e-7;
    }
    st.s_stagnation = sst;
    st.s_stagnation_d_gamma_node0 = (sst - st.s[i + 1]) / dgam;
    st.s_stagnation_d_gamma_node1 = (st.s[i] - sst) / dgam;
}

/// IBLPAN (xpanel.f): BL station -> panel node pointers IPAN, the VTI sign, IBLTE and NBL.
pub fn iblpan(st: &mut SolverState) {
    let (n, nw, ist) = (st.n_foil_nodes, st.n_wake_nodes, st.i_stagnation_node);

    // top surface first
    let is = 1;
    let mut ibl = 1;
    for i in (1..=ist).rev() {
        ibl += 1;
        st.i_node[is][ibl] = i;
        st.velocity_sign[is][ibl] = 1.0;
    }
    st.i_te_station[is] = ibl;
    st.n_stations[is] = ibl;

    // bottom surface next
    let is = 2;
    let mut ibl = 1;
    for i in (ist + 1)..=n {
        ibl += 1;
        st.i_node[is][ibl] = i;
        st.velocity_sign[is][ibl] = -1.0;
    }

    // wake
    st.i_te_station[is] = ibl;
    for iw in 1..=nw {
        let i = n + iw;
        let ibl = st.i_te_station[is] + iw;
        st.i_node[is][ibl] = i;
        st.velocity_sign[is][ibl] = -1.0;
    }
    st.n_stations[is] = st.i_te_station[is] + nw;

    // upper wake pointers (for plotting only)
    for iw in 1..=nw {
        st.i_node[1][st.i_te_station[1] + iw] = st.i_node[2][st.i_te_station[2] + iw];
        st.velocity_sign[1][st.i_te_station[1] + iw] = 1.0;
    }
    // LIPAN = .TRUE.
    st.pointers_built = true;
}

/// XICALC (xpanel.f): BL arc length XSSI on each side and along the wake, and the TE
/// "dead air" gap array WGAP.
pub fn xicalc(st: &mut SolverState) {
    let n = st.n_foil_nodes;
    let xfeps = 1.0e-7;
    // minimum xi node arc length near stagnation point
    let xeps = xfeps * (st.s[n] - st.s[1]);

    let is = 1;
    st.xi[is][1] = 0.0;
    for ibl in 2..=st.i_te_station[is] {
        let i = st.i_node[is][ibl];
        st.xi[is][ibl] = (st.s_stagnation - st.s[i]).max(xeps);
    }

    let is = 2;
    st.xi[is][1] = 0.0;
    for ibl in 2..=st.i_te_station[is] {
        let i = st.i_node[is][ibl];
        st.xi[is][ibl] = (st.s[i] - st.s_stagnation).max(xeps);
    }

    let (is1, is2) = (1, 2);
    let ibl1 = st.i_te_station[is1] + 1;
    st.xi[is1][ibl1] = st.xi[is1][ibl1 - 1];
    let ibl2 = st.i_te_station[is2] + 1;
    st.xi[is2][ibl2] = st.xi[is2][ibl2 - 1];

    // DO 25 IBL=IBLTE(IS)+2, NBL(IS)   with IS = 2 from the preceding loop
    for ibl in (st.i_te_station[is] + 2)..=st.n_stations[is] {
        let i = st.i_node[is][ibl];
        let dxssi = ((st.x[i] - st.x[i - 1]).powi(2) + (st.y[i] - st.y[i - 1]).powi(2)).sqrt();
        let ibl1 = st.i_te_station[is1] + ibl - st.i_te_station[is];
        let ibl2 = st.i_te_station[is2] + ibl - st.i_te_station[is];
        st.xi[is1][ibl1] = st.xi[is1][ibl1 - 1] + dxssi;
        st.xi[is2][ibl2] = st.xi[is2][ibl2 - 1] + dxssi;
    }

    // trailing edge flap length to TE gap ratio
    let telrat = 2.50;

    // set up parameters for TE flap cubics
    let crosp = (st.dxds[1] * st.dyds[n] - st.dyds[1] * st.dxds[n])
        / ((st.dxds[1].powi(2) + st.dyds[1].powi(2)) * (st.dxds[n].powi(2) + st.dyds[n].powi(2))).sqrt();
    let mut dwdxte = crosp / (1.0 - crosp.powi(2)).sqrt();

    // limit cubic to avoid absurd TE gap widths
    dwdxte = dwdxte.max(-3.0 / telrat);
    dwdxte = dwdxte.min(3.0 / telrat);

    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;

    if st.sharp_te {
        for iw in 1..=st.n_wake_nodes {
            st.wake_gap[iw] = 0.0;
        }
    } else {
        // set TE flap (wake gap) array
        let is = 2;
        for iw in 1..=st.n_wake_nodes {
            let ibl = st.i_te_station[is] + iw;
            let zn = 1.0 - (st.xi[is][ibl] - st.xi[is][st.i_te_station[is]]) / (telrat * st.te_thickness_normal);
            st.wake_gap[iw] = 0.0;
            if zn >= 0.0 {
                st.wake_gap[iw] = st.te_thickness_normal * (aa + bb * zn) * zn.powi(2);
            }
        }
    }
}

/// IBLSYS (xbl.f): Newton-system row number ISYS for each BL station, and NSYS.
pub fn iblsys(st: &mut SolverState) {
    let mut iv = 0;
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            iv += 1;
            st.i_row[is][ibl] = iv;
        }
    }
    st.n_rows = iv;
}

/// SINVRT (spline.f): inverse spline S(X) by Newton iteration from the initial guess `si`.
/// Returns the input value if 10 iterations do not converge (XFOIL prints a warning).
pub fn sinvrt(mut si: f64, xi: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let n = s.len();
    let sisav = si;
    for _ in 0..10 {
        let res = seval(si, x, xs, s) - xi;
        let resp = deval(si, x, xs, s);
        let ds = -res / resp;
        si += ds;
        if (ds / (s[n - 1] - s[0])).abs() < 1.0e-5 {
            return si;
        }
    }
    sisav
}

/// XIFSET (xbl.f): forced-transition BL coordinate XIFORC for side `is`.
pub fn xifset(st: &SolverState, is: usize) -> f64 {
    if st.x_trip[is] >= 1.0 {
        return st.xi[is][st.i_te_station[is]];
    }
    let n = st.n_foil_nodes;
    let chx = st.x_te - st.x_le;
    let chy = st.y_te - st.y_le;
    let chsq = chx * chx + chy * chy;

    // calculate chord-based x/c, y/c (0-based slices for the spline helpers)
    let s: Vec<f64> = st.s[1..=n].to_vec();
    let mut w1 = vec![0.0; n];
    let mut w2 = vec![0.0; n];
    for i in 1..=n {
        w1[i - 1] = ((st.x[i] - st.x_le) * chx + (st.y[i] - st.y_le) * chy) / chsq;
        w2[i - 1] = ((st.y[i] - st.y_le) * chx - (st.x[i] - st.x_le) * chy) / chsq;
    }
    let w3 = spline(&w1, &s); // SPLIND(W1,W3,S,N,-999.0,-999.0)
    let _w4 = spline(&w2, &s); // SPLIND(W2,W4,S,N,-999.0,-999.0) — computed, unused by XFOIL too

    let mut xiforc = if is == 1 {
        // set approximate arc length of forced transition point for SINVRT
        let str0 = st.s_le + (st.s[1] - st.s_le) * st.x_trip[is];
        let str_ = sinvrt(str0, st.x_trip[is], &w1, &w3, &s);
        (st.s_stagnation - str_).min(st.xi[is][st.i_te_station[is]])
    } else {
        // same for bottom side
        let str0 = st.s_le + (st.s[n] - st.s_le) * st.x_trip[is];
        let str_ = sinvrt(str0, st.x_trip[is], &w1, &w3, &s);
        (str_ - st.s_stagnation).min(st.xi[is][st.i_te_station[is]])
    };

    if xiforc < 0.0 {
        // "Stagnation point is past trip on side IS"
        xiforc = st.xi[is][st.i_te_station[is]];
    }
    xiforc
}

/// STMOVE: moves the stagnation point location to a new panel. Re-runs STFIND on the current
/// GAM; if IST is unchanged only XICALC is redone, otherwise the pointer layer is rebuilt and
/// the BL arrays and ITRAN are shifted by IDIF. Always refreshes MASS = DSTR*UEDG.
pub fn stmove(st: &mut SolverState) {
    // locate new stagnation point arc length SST from GAM distribution
    let istold = st.i_stagnation_node;
    stfind(st);

    if istold == st.i_stagnation_node {
        // recalculate new arc length array
        xicalc(st);
    } else {
        // set new BL position -> panel position pointers
        iblpan(st);
        // set new inviscid BL edge velocity UINV from QINV
        crate::solver::velocity::uicalc(st);
        // recalculate new arc length array
        xicalc(st);
        // set BL position -> system line pointers
        iblsys(st);

        if st.i_stagnation_node > istold {
            // increase in number of points on top side (IS=1)
            let idif = st.i_stagnation_node - istold;
            st.i_transition_station[1] += idif;
            st.i_transition_station[2] -= idif;

            // move top side BL variables downstream
            for ibl in (idif + 2..=st.n_stations[1]).rev() {
                st.sqrtctau[1][ibl] = st.sqrtctau[1][ibl - idif];
                st.theta[1][ibl] = st.theta[1][ibl - idif];
                st.dstar[1][ibl] = st.dstar[1][ibl - idif];
                st.ue[1][ibl] = st.ue[1][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = st.ue[1][idif + 2] / st.xi[1][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                st.sqrtctau[1][ibl] = st.sqrtctau[1][idif + 2];
                st.theta[1][ibl] = st.theta[1][idif + 2];
                st.dstar[1][ibl] = st.dstar[1][idif + 2];
                st.ue[1][ibl] = dudx * st.xi[1][ibl];
            }

            // move bottom side BL variables upstream
            for ibl in 2..=st.n_stations[2] {
                st.sqrtctau[2][ibl] = st.sqrtctau[2][ibl + idif];
                st.theta[2][ibl] = st.theta[2][ibl + idif];
                st.dstar[2][ibl] = st.dstar[2][ibl + idif];
                st.ue[2][ibl] = st.ue[2][ibl + idif];
            }
        } else {
            // increase in number of points on bottom side (IS=2)
            let idif = istold - st.i_stagnation_node;
            st.i_transition_station[1] -= idif;
            st.i_transition_station[2] += idif;

            // move bottom side BL variables downstream
            for ibl in (idif + 2..=st.n_stations[2]).rev() {
                st.sqrtctau[2][ibl] = st.sqrtctau[2][ibl - idif];
                st.theta[2][ibl] = st.theta[2][ibl - idif];
                st.dstar[2][ibl] = st.dstar[2][ibl - idif];
                st.ue[2][ibl] = st.ue[2][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = st.ue[2][idif + 2] / st.xi[2][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                st.sqrtctau[2][ibl] = st.sqrtctau[2][idif + 2];
                st.theta[2][ibl] = st.theta[2][idif + 2];
                st.dstar[2][ibl] = st.dstar[2][idif + 2];
                st.ue[2][ibl] = dudx * st.xi[2][ibl];
            }

            // move top side BL variables upstream
            for ibl in 2..=st.n_stations[1] {
                st.sqrtctau[1][ibl] = st.sqrtctau[1][ibl + idif];
                st.theta[1][ibl] = st.theta[1][ibl + idif];
                st.dstar[1][ibl] = st.dstar[1][ibl + idif];
                st.ue[1][ibl] = st.ue[1][ibl + idif];
            }
        }

        // tweak Ue so it's not zero, in case stag. point is right on node
        let ueps = 1.0e-7;
        for is in 1..=2 {
            for ibl in 2..=st.n_stations[is] {
                let i = st.i_node[is][ibl];
                if st.ue[is][ibl] <= ueps {
                    st.ue[is][ibl] = ueps;
                    st.q_viscous[i] = st.velocity_sign[is][ibl] * ueps;
                    st.gamma[i] = st.velocity_sign[is][ibl] * ueps;
                }
            }
        }
    }

    // set new mass array since Ue has been tweaked
    for is in 1..=2 {
        for ibl in 2..=st.n_stations[is] {
            st.mass_defect[is][ibl] = st.dstar[is][ibl] * st.ue[is][ibl];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// XSTRIP >= 1 (free transition) returns the TE arc length without touching the spline.
    #[test]
    fn xifset_free_transition_returns_te_arc_length() {
        let mut st = SolverState::empty(10, 3);
        st.i_stagnation_node = 5;
        iblpan(&mut st);
        st.xi[1][st.i_te_station[1]] = 1.25;
        st.xi[2][st.i_te_station[2]] = 1.5;
        st.x_trip = [0.0, 1.0, 1.0];
        assert_eq!(xifset(&st, 1), 1.25);
        assert_eq!(xifset(&st, 2), 1.5);
    }

    /// SINVRT on the identity spline x(s) = s returns s = xi.
    #[test]
    fn sinvrt_identity_spline() {
        let s: Vec<f64> = (0..11).map(|i| i as f64 * 0.1).collect();
        let x = s.clone();
        let xs = spline(&x, &s);
        let si = sinvrt(0.3, 0.42, &x, &xs, &s);
        assert!((si - 0.42).abs() < 1e-5 * (s[10] - s[0]) * 10.0, "si = {si}");
    }
}
