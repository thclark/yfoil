//! BL pointer layer: STFIND, IBLPAN, XICALC, IBLSYS, XIFSET (xpanel.f / xbl.f), plus TECALC.
//!
//! Each function is a line-for-line translation; the Fortran statement is quoted where the
//! translation is not obvious. Indexing is 1-based (see `blstate.rs`).

use crate::geometry::{spline_derivatives, spline_slope, spline_value};
use crate::solver::blstate::SolverState;

/// TECALC (xfoil.f): TE gap areas and the sharp-TE flag.
#[doc(alias = "TECALC")]
pub fn set_te_thickness(state: &mut SolverState) {
    let n = state.n_foil_nodes;
    // set TE base vector and TE bisector components
    let dxte = state.x[1] - state.x[n];
    let dyte = state.y[1] - state.y[n];
    let dxs = 0.5 * (-state.dxds[1] + state.dxds[n]);
    let dys = 0.5 * (-state.dyds[1] + state.dyds[n]);
    // normal and streamwise projected TE gap areas
    state.te_thickness_normal = dxs * dyte - dys * dxte;
    state.te_thickness_parallel = dxs * dxte + dys * dyte;
    // total TE gap area
    state.te_gap = (dxte * dxte + dyte * dyte).sqrt();
    state.sharp_te = state.te_gap < 0.0001 * state.chord;
}

/// STFIND (xpanel.f): stagnation point arc length SST, panel index IST, and the
/// sensitivities SST_GO = dSST/dGAM(IST), SST_GP = dSST/dGAM(IST+1).
#[doc(alias = "STFIND")]
pub fn find_stagnation(state: &mut SolverState) {
    let n = state.n_foil_nodes;
    let mut i = n / 2; // fallback if no sign change side found ("Stagnation point not found")
    for ii in 1..n {
        if state.gamma[ii] >= 0.0 && state.gamma[ii + 1] < 0.0 {
            i = ii;
            break;
        }
    }
    state.i_stagnation_node = i;
    let dgam = state.gamma[i + 1] - state.gamma[i];
    let ds = state.s[i + 1] - state.s[i];

    // evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
    let mut sst = if state.gamma[i] < -state.gamma[i + 1] {
        state.s[i] - ds * (state.gamma[i] / dgam)
    } else {
        state.s[i + 1] - ds * (state.gamma[i + 1] / dgam)
    };

    // tweak stagnation point if it falls right on a node (very unlikely)
    if sst <= state.s[i] {
        sst = state.s[i] + 1.0e-7;
    }
    if sst >= state.s[i + 1] {
        sst = state.s[i + 1] - 1.0e-7;
    }
    state.s_stagnation = sst;
    state.s_stagnation_d_gamma_node0 = (sst - state.s[i + 1]) / dgam;
    state.s_stagnation_d_gamma_node1 = (state.s[i] - sst) / dgam;
}

/// IBLPAN (xpanel.f): BL station -> panel node pointers IPAN, the VTI sign, IBLTE and NBL.
#[doc(alias = "IBLPAN")]
pub fn map_stations_to_nodes(state: &mut SolverState) {
    let (n, nw, ist) = (state.n_foil_nodes, state.n_wake_nodes, state.i_stagnation_node);

    // top surface first
    let side = 1;
    let mut i_station = 1;
    for i in (1..=ist).rev() {
        i_station += 1;
        state.i_node[side][i_station] = i;
        state.velocity_sign[side][i_station] = 1.0;
    }
    state.i_te_station[side] = i_station;
    state.n_stations[side] = i_station;

    // bottom surface next
    let side = 2;
    let mut i_station = 1;
    for i in (ist + 1)..=n {
        i_station += 1;
        state.i_node[side][i_station] = i;
        state.velocity_sign[side][i_station] = -1.0;
    }

    // wake
    state.i_te_station[side] = i_station;
    for i_wake in 1..=nw {
        let i = n + i_wake;
        let i_station = state.i_te_station[side] + i_wake;
        state.i_node[side][i_station] = i;
        state.velocity_sign[side][i_station] = -1.0;
    }
    state.n_stations[side] = state.i_te_station[side] + nw;

    // upper wake pointers (for plotting only)
    for i_wake in 1..=nw {
        state.i_node[1][state.i_te_station[1] + i_wake] = state.i_node[2][state.i_te_station[2] + i_wake];
        state.velocity_sign[1][state.i_te_station[1] + i_wake] = 1.0;
    }
    // LIPAN = .TRUE.
    state.pointers_built = true;
}

/// XICALC (xpanel.f): BL arc length XSSI on each side and along the wake, and the TE
/// "dead air" gap array WGAP.
#[doc(alias = "XICALC")]
pub fn set_station_xi(state: &mut SolverState) {
    let n = state.n_foil_nodes;
    let xfeps = 1.0e-7;
    // minimum xi node arc length near stagnation point
    let xeps = xfeps * (state.s[n] - state.s[1]);

    let side = 1;
    state.xi[side][1] = 0.0;
    for i_station in 2..=state.i_te_station[side] {
        let i = state.i_node[side][i_station];
        state.xi[side][i_station] = (state.s_stagnation - state.s[i]).max(xeps);
    }

    let side = 2;
    state.xi[side][1] = 0.0;
    for i_station in 2..=state.i_te_station[side] {
        let i = state.i_node[side][i_station];
        state.xi[side][i_station] = (state.s[i] - state.s_stagnation).max(xeps);
    }

    let (is1, is2) = (1, 2);
    let ibl1 = state.i_te_station[is1] + 1;
    state.xi[is1][ibl1] = state.xi[is1][ibl1 - 1];
    let ibl2 = state.i_te_station[is2] + 1;
    state.xi[is2][ibl2] = state.xi[is2][ibl2 - 1];

    // DO 25 IBL=IBLTE(IS)+2, NBL(IS)   with IS = 2 from the preceding loop
    for i_station in (state.i_te_station[side] + 2)..=state.n_stations[side] {
        let i = state.i_node[side][i_station];
        let dxssi = ((state.x[i] - state.x[i - 1]).powi(2) + (state.y[i] - state.y[i - 1]).powi(2)).sqrt();
        let ibl1 = state.i_te_station[is1] + i_station - state.i_te_station[side];
        let ibl2 = state.i_te_station[is2] + i_station - state.i_te_station[side];
        state.xi[is1][ibl1] = state.xi[is1][ibl1 - 1] + dxssi;
        state.xi[is2][ibl2] = state.xi[is2][ibl2 - 1] + dxssi;
    }

    // trailing edge flap length to TE gap ratio
    let telrat = 2.50;

    // set up parameters for TE flap cubics
    let crosp = (state.dxds[1] * state.dyds[n] - state.dyds[1] * state.dxds[n])
        / ((state.dxds[1].powi(2) + state.dyds[1].powi(2)) * (state.dxds[n].powi(2) + state.dyds[n].powi(2))).sqrt();
    let mut dwdxte = crosp / (1.0 - crosp.powi(2)).sqrt();

    // limit cubic to avoid absurd TE gap widths
    dwdxte = dwdxte.max(-3.0 / telrat);
    dwdxte = dwdxte.min(3.0 / telrat);

    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;

    if state.sharp_te {
        for i_wake in 1..=state.n_wake_nodes {
            state.wake_gap[i_wake] = 0.0;
        }
    } else {
        // set TE flap (wake gap) array
        let side = 2;
        for i_wake in 1..=state.n_wake_nodes {
            let i_station = state.i_te_station[side] + i_wake;
            let zn = 1.0
                - (state.xi[side][i_station] - state.xi[side][state.i_te_station[side]])
                    / (telrat * state.te_thickness_normal);
            state.wake_gap[i_wake] = 0.0;
            if zn >= 0.0 {
                state.wake_gap[i_wake] = state.te_thickness_normal * (aa + bb * zn) * zn.powi(2);
            }
        }
    }
}

/// IBLSYS (xbl.f): Newton-system row number ISYS for each BL station, and NSYS.
#[doc(alias = "IBLSYS")]
pub fn map_stations_to_rows(state: &mut SolverState) {
    let mut i_row = 0;
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            i_row += 1;
            state.i_row[side][i_station] = i_row;
        }
    }
    state.n_rows = i_row;
}

/// SINVRT (spline.f): inverse spline S(X) by Newton iteration from the initial guess `si`.
/// Returns the input value if 10 iterations do not converge (XFOIL prints a warning).
#[doc(alias = "SINVRT")]
pub fn s_at_x(mut si: f64, xi: f64, x: &[f64], xs: &[f64], s: &[f64]) -> f64 {
    let n = s.len();
    let sisav = si;
    for _ in 0..10 {
        let res = spline_value(si, x, xs, s) - xi;
        let resp = spline_slope(si, x, xs, s);
        let ds = -res / resp;
        si += ds;
        if (ds / (s[n - 1] - s[0])).abs() < 1.0e-5 {
            return si;
        }
    }
    sisav
}

/// XIFSET (xbl.f): forced-transition BL coordinate XIFORC for side `is`.
#[doc(alias = "XIFSET")]
pub fn xi_trip(state: &SolverState, side: usize) -> f64 {
    if state.x_trip[side] >= 1.0 {
        return state.xi[side][state.i_te_station[side]];
    }
    let n = state.n_foil_nodes;
    let chx = state.x_te - state.x_le;
    let chy = state.y_te - state.y_le;
    let chsq = chx * chx + chy * chy;

    // calculate chord-based x/c, y/c (0-based slices for the spline helpers)
    let s: Vec<f64> = state.s[1..=n].to_vec();
    let mut w1 = vec![0.0; n];
    let mut w2 = vec![0.0; n];
    for i in 1..=n {
        w1[i - 1] = ((state.x[i] - state.x_le) * chx + (state.y[i] - state.y_le) * chy) / chsq;
        w2[i - 1] = ((state.y[i] - state.y_le) * chx - (state.x[i] - state.x_le) * chy) / chsq;
    }
    let w3 = spline_derivatives(&w1, &s); // SPLIND(W1,W3,S,N,-999.0,-999.0)
    let _w4 = spline_derivatives(&w2, &s); // SPLIND(W2,W4,S,N,-999.0,-999.0) — computed, unused by XFOIL too

    let mut xiforc = if side == 1 {
        // set approximate arc length of forced transition point for SINVRT
        let str0 = state.s_le + (state.s[1] - state.s_le) * state.x_trip[side];
        let str_ = s_at_x(str0, state.x_trip[side], &w1, &w3, &s);
        (state.s_stagnation - str_).min(state.xi[side][state.i_te_station[side]])
    } else {
        // same for bottom side
        let str0 = state.s_le + (state.s[n] - state.s_le) * state.x_trip[side];
        let str_ = s_at_x(str0, state.x_trip[side], &w1, &w3, &s);
        (str_ - state.s_stagnation).min(state.xi[side][state.i_te_station[side]])
    };

    if xiforc < 0.0 {
        // "Stagnation point is past trip on side IS"
        xiforc = state.xi[side][state.i_te_station[side]];
    }
    xiforc
}

/// STMOVE: moves the stagnation point location to a new panel. Re-runs STFIND on the current
/// GAM; if IST is unchanged only XICALC is redone, otherwise the pointer layer is rebuilt and
/// the BL arrays and ITRAN are shifted by IDIF. Always refreshes MASS = DSTR*UEDG.
#[doc(alias = "STMOVE")]
pub fn move_stagnation(state: &mut SolverState) {
    // locate new stagnation point arc length SST from GAM distribution
    let istold = state.i_stagnation_node;
    find_stagnation(state);

    if istold == state.i_stagnation_node {
        // recalculate new arc length array
        set_station_xi(state);
    } else {
        // set new BL position -> panel position pointers
        map_stations_to_nodes(state);
        // set new inviscid BL edge velocity UINV from QINV
        crate::solver::velocity::set_ue_inviscid(state);
        // recalculate new arc length array
        set_station_xi(state);
        // set BL position -> system line pointers
        map_stations_to_rows(state);

        if state.i_stagnation_node > istold {
            // increase in number of points on top side (IS=1)
            let idif = state.i_stagnation_node - istold;
            state.i_transition_station[1] += idif;
            state.i_transition_station[2] -= idif;

            // move top side BL variables downstream
            for i_station in (idif + 2..=state.n_stations[1]).rev() {
                state.sqrtctau[1][i_station] = state.sqrtctau[1][i_station - idif];
                state.theta[1][i_station] = state.theta[1][i_station - idif];
                state.dstar[1][i_station] = state.dstar[1][i_station - idif];
                state.ue[1][i_station] = state.ue[1][i_station - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = state.ue[1][idif + 2] / state.xi[1][idif + 2];
            for i_station in (2..=idif + 1).rev() {
                state.sqrtctau[1][i_station] = state.sqrtctau[1][idif + 2];
                state.theta[1][i_station] = state.theta[1][idif + 2];
                state.dstar[1][i_station] = state.dstar[1][idif + 2];
                state.ue[1][i_station] = dudx * state.xi[1][i_station];
            }

            // move bottom side BL variables upstream
            for i_station in 2..=state.n_stations[2] {
                state.sqrtctau[2][i_station] = state.sqrtctau[2][i_station + idif];
                state.theta[2][i_station] = state.theta[2][i_station + idif];
                state.dstar[2][i_station] = state.dstar[2][i_station + idif];
                state.ue[2][i_station] = state.ue[2][i_station + idif];
            }
        } else {
            // increase in number of points on bottom side (IS=2)
            let idif = istold - state.i_stagnation_node;
            state.i_transition_station[1] -= idif;
            state.i_transition_station[2] += idif;

            // move bottom side BL variables downstream
            for i_station in (idif + 2..=state.n_stations[2]).rev() {
                state.sqrtctau[2][i_station] = state.sqrtctau[2][i_station - idif];
                state.theta[2][i_station] = state.theta[2][i_station - idif];
                state.dstar[2][i_station] = state.dstar[2][i_station - idif];
                state.ue[2][i_station] = state.ue[2][i_station - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = state.ue[2][idif + 2] / state.xi[2][idif + 2];
            for i_station in (2..=idif + 1).rev() {
                state.sqrtctau[2][i_station] = state.sqrtctau[2][idif + 2];
                state.theta[2][i_station] = state.theta[2][idif + 2];
                state.dstar[2][i_station] = state.dstar[2][idif + 2];
                state.ue[2][i_station] = dudx * state.xi[2][i_station];
            }

            // move top side BL variables upstream
            for i_station in 2..=state.n_stations[1] {
                state.sqrtctau[1][i_station] = state.sqrtctau[1][i_station + idif];
                state.theta[1][i_station] = state.theta[1][i_station + idif];
                state.dstar[1][i_station] = state.dstar[1][i_station + idif];
                state.ue[1][i_station] = state.ue[1][i_station + idif];
            }
        }

        // tweak Ue so it's not zero, in case stag. point is right on node
        let ueps = 1.0e-7;
        for side in 1..=2 {
            for i_station in 2..=state.n_stations[side] {
                let i = state.i_node[side][i_station];
                if state.ue[side][i_station] <= ueps {
                    state.ue[side][i_station] = ueps;
                    state.q_viscous[i] = state.velocity_sign[side][i_station] * ueps;
                    state.gamma[i] = state.velocity_sign[side][i_station] * ueps;
                }
            }
        }
    }

    // set new mass array since Ue has been tweaked
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            state.mass_defect[side][i_station] = state.dstar[side][i_station] * state.ue[side][i_station];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// XSTRIP >= 1 (free transition) returns the TE arc length without touching the spline.
    #[test]
    fn xifset_free_transition_returns_te_arc_length() {
        let mut state = SolverState::empty(10, 3);
        state.i_stagnation_node = 5;
        map_stations_to_nodes(&mut state);
        state.xi[1][state.i_te_station[1]] = 1.25;
        state.xi[2][state.i_te_station[2]] = 1.5;
        state.x_trip = [0.0, 1.0, 1.0];
        assert_eq!(xi_trip(&state, 1), 1.25);
        assert_eq!(xi_trip(&state, 2), 1.5);
    }

    /// SINVRT on the identity spline x(s) = s returns s = xi.
    #[test]
    fn sinvrt_identity_spline() {
        let s: Vec<f64> = (0..11).map(|i| i as f64 * 0.1).collect();
        let x = s.clone();
        let xs = spline_derivatives(&x, &s);
        let si = s_at_x(0.3, 0.42, &x, &xs, &s);
        assert!((si - 0.42).abs() < 1e-5 * (s[10] - s[0]) * 10.0, "si = {si}");
    }
}
