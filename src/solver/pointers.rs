//! BL pointer layer: STFIND, IBLPAN, XICALC, IBLSYS, XIFSET (xpanel.f / xbl.f), plus TECALC.
//!
//! Each function is a line-for-line translation; the Fortran statement is quoted where the
//! translation is not obvious. Indexing is 1-based (see `blstate.rs`).

use crate::geometry::{deval, seval, spline};
use crate::solver::blstate::BlState;

/// TECALC (xfoil.f): TE gap areas and the sharp-TE flag.
pub fn tecalc(st: &mut BlState) {
    let n = st.n;
    // set TE base vector and TE bisector components
    let dxte = st.x[1] - st.x[n];
    let dyte = st.y[1] - st.y[n];
    let dxs = 0.5 * (-st.xp[1] + st.xp[n]);
    let dys = 0.5 * (-st.yp[1] + st.yp[n]);
    // normal and streamwise projected TE gap areas
    st.ante = dxs * dyte - dys * dxte;
    st.aste = dxs * dxte + dys * dyte;
    // total TE gap area
    st.dste = (dxte * dxte + dyte * dyte).sqrt();
    st.sharp = st.dste < 0.0001 * st.chord;
}

/// STFIND (xpanel.f): stagnation point arc length SST, panel index IST, and the
/// sensitivities SST_GO = dSST/dGAM(IST), SST_GP = dSST/dGAM(IST+1).
pub fn stfind(st: &mut BlState) {
    let n = st.n;
    let mut i = n / 2; // fallback if no sign change is found ("Stagnation point not found")
    for ii in 1..n {
        if st.gam[ii] >= 0.0 && st.gam[ii + 1] < 0.0 {
            i = ii;
            break;
        }
    }
    st.ist = i;
    let dgam = st.gam[i + 1] - st.gam[i];
    let ds = st.s[i + 1] - st.s[i];

    // evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
    let mut sst = if st.gam[i] < -st.gam[i + 1] {
        st.s[i] - ds * (st.gam[i] / dgam)
    } else {
        st.s[i + 1] - ds * (st.gam[i + 1] / dgam)
    };

    // tweak stagnation point if it falls right on a node (very unlikely)
    if sst <= st.s[i] {
        sst = st.s[i] + 1.0e-7;
    }
    if sst >= st.s[i + 1] {
        sst = st.s[i + 1] - 1.0e-7;
    }
    st.sst = sst;
    st.sst_go = (sst - st.s[i + 1]) / dgam;
    st.sst_gp = (st.s[i] - sst) / dgam;
}

/// IBLPAN (xpanel.f): BL station -> panel node pointers IPAN, the VTI sign, IBLTE and NBL.
pub fn iblpan(st: &mut BlState) {
    let (n, nw, ist) = (st.n, st.nw, st.ist);

    // top surface first
    let is = 1;
    let mut ibl = 1;
    for i in (1..=ist).rev() {
        ibl += 1;
        st.ipan[is][ibl] = i;
        st.vti[is][ibl] = 1.0;
    }
    st.iblte[is] = ibl;
    st.nbl[is] = ibl;

    // bottom surface next
    let is = 2;
    let mut ibl = 1;
    for i in (ist + 1)..=n {
        ibl += 1;
        st.ipan[is][ibl] = i;
        st.vti[is][ibl] = -1.0;
    }

    // wake
    st.iblte[is] = ibl;
    for iw in 1..=nw {
        let i = n + iw;
        let ibl = st.iblte[is] + iw;
        st.ipan[is][ibl] = i;
        st.vti[is][ibl] = -1.0;
    }
    st.nbl[is] = st.iblte[is] + nw;

    // upper wake pointers (for plotting only)
    for iw in 1..=nw {
        st.ipan[1][st.iblte[1] + iw] = st.ipan[2][st.iblte[2] + iw];
        st.vti[1][st.iblte[1] + iw] = 1.0;
    }
    // LIPAN = .TRUE.
    st.lipan = true;
}

/// XICALC (xpanel.f): BL arc length XSSI on each side and along the wake, and the TE
/// "dead air" gap array WGAP.
pub fn xicalc(st: &mut BlState) {
    let n = st.n;
    let xfeps = 1.0e-7;
    // minimum xi node arc length near stagnation point
    let xeps = xfeps * (st.s[n] - st.s[1]);

    let is = 1;
    st.xssi[is][1] = 0.0;
    for ibl in 2..=st.iblte[is] {
        let i = st.ipan[is][ibl];
        st.xssi[is][ibl] = (st.sst - st.s[i]).max(xeps);
    }

    let is = 2;
    st.xssi[is][1] = 0.0;
    for ibl in 2..=st.iblte[is] {
        let i = st.ipan[is][ibl];
        st.xssi[is][ibl] = (st.s[i] - st.sst).max(xeps);
    }

    let (is1, is2) = (1, 2);
    let ibl1 = st.iblte[is1] + 1;
    st.xssi[is1][ibl1] = st.xssi[is1][ibl1 - 1];
    let ibl2 = st.iblte[is2] + 1;
    st.xssi[is2][ibl2] = st.xssi[is2][ibl2 - 1];

    // DO 25 IBL=IBLTE(IS)+2, NBL(IS)   with IS = 2 from the preceding loop
    for ibl in (st.iblte[is] + 2)..=st.nbl[is] {
        let i = st.ipan[is][ibl];
        let dxssi = ((st.x[i] - st.x[i - 1]).powi(2) + (st.y[i] - st.y[i - 1]).powi(2)).sqrt();
        let ibl1 = st.iblte[is1] + ibl - st.iblte[is];
        let ibl2 = st.iblte[is2] + ibl - st.iblte[is];
        st.xssi[is1][ibl1] = st.xssi[is1][ibl1 - 1] + dxssi;
        st.xssi[is2][ibl2] = st.xssi[is2][ibl2 - 1] + dxssi;
    }

    // trailing edge flap length to TE gap ratio
    let telrat = 2.50;

    // set up parameters for TE flap cubics
    let crosp = (st.xp[1] * st.yp[n] - st.yp[1] * st.xp[n])
        / ((st.xp[1].powi(2) + st.yp[1].powi(2)) * (st.xp[n].powi(2) + st.yp[n].powi(2))).sqrt();
    let mut dwdxte = crosp / (1.0 - crosp.powi(2)).sqrt();

    // limit cubic to avoid absurd TE gap widths
    dwdxte = dwdxte.max(-3.0 / telrat);
    dwdxte = dwdxte.min(3.0 / telrat);

    let aa = 3.0 + telrat * dwdxte;
    let bb = -2.0 - telrat * dwdxte;

    if st.sharp {
        for iw in 1..=st.nw {
            st.wgap[iw] = 0.0;
        }
    } else {
        // set TE flap (wake gap) array
        let is = 2;
        for iw in 1..=st.nw {
            let ibl = st.iblte[is] + iw;
            let zn = 1.0 - (st.xssi[is][ibl] - st.xssi[is][st.iblte[is]]) / (telrat * st.ante);
            st.wgap[iw] = 0.0;
            if zn >= 0.0 {
                st.wgap[iw] = st.ante * (aa + bb * zn) * zn.powi(2);
            }
        }
    }
}

/// IBLSYS (xbl.f): Newton-system row number ISYS for each BL station, and NSYS.
pub fn iblsys(st: &mut BlState) {
    let mut iv = 0;
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            iv += 1;
            st.isys[is][ibl] = iv;
        }
    }
    st.nsys = iv;
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
pub fn xifset(st: &BlState, is: usize) -> f64 {
    if st.xstrip[is] >= 1.0 {
        return st.xssi[is][st.iblte[is]];
    }
    let n = st.n;
    let chx = st.xte - st.xle;
    let chy = st.yte - st.yle;
    let chsq = chx * chx + chy * chy;

    // calculate chord-based x/c, y/c (0-based slices for the spline helpers)
    let s: Vec<f64> = st.s[1..=n].to_vec();
    let mut w1 = vec![0.0; n];
    let mut w2 = vec![0.0; n];
    for i in 1..=n {
        w1[i - 1] = ((st.x[i] - st.xle) * chx + (st.y[i] - st.yle) * chy) / chsq;
        w2[i - 1] = ((st.y[i] - st.yle) * chx - (st.x[i] - st.xle) * chy) / chsq;
    }
    let w3 = spline(&w1, &s); // SPLIND(W1,W3,S,N,-999.0,-999.0)
    let _w4 = spline(&w2, &s); // SPLIND(W2,W4,S,N,-999.0,-999.0) — computed, unused by XFOIL too

    let mut xiforc = if is == 1 {
        // set approximate arc length of forced transition point for SINVRT
        let str0 = st.sle + (st.s[1] - st.sle) * st.xstrip[is];
        let str_ = sinvrt(str0, st.xstrip[is], &w1, &w3, &s);
        (st.sst - str_).min(st.xssi[is][st.iblte[is]])
    } else {
        // same for bottom side
        let str0 = st.sle + (st.s[n] - st.sle) * st.xstrip[is];
        let str_ = sinvrt(str0, st.xstrip[is], &w1, &w3, &s);
        (str_ - st.sst).min(st.xssi[is][st.iblte[is]])
    };

    if xiforc < 0.0 {
        // "Stagnation point is past trip on side IS"
        xiforc = st.xssi[is][st.iblte[is]];
    }
    xiforc
}

/// STMOVE: moves the stagnation point location to a new panel. Re-runs STFIND on the current
/// GAM; if IST is unchanged only XICALC is redone, otherwise the pointer layer is rebuilt and
/// the BL arrays and ITRAN are shifted by IDIF. Always refreshes MASS = DSTR*UEDG.
pub fn stmove(st: &mut BlState) {
    // locate new stagnation point arc length SST from GAM distribution
    let istold = st.ist;
    stfind(st);

    if istold == st.ist {
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

        if st.ist > istold {
            // increase in number of points on top side (IS=1)
            let idif = st.ist - istold;
            st.itran[1] += idif;
            st.itran[2] -= idif;

            // move top side BL variables downstream
            for ibl in (idif + 2..=st.nbl[1]).rev() {
                st.ctau[1][ibl] = st.ctau[1][ibl - idif];
                st.thet[1][ibl] = st.thet[1][ibl - idif];
                st.dstr[1][ibl] = st.dstr[1][ibl - idif];
                st.uedg[1][ibl] = st.uedg[1][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = st.uedg[1][idif + 2] / st.xssi[1][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                st.ctau[1][ibl] = st.ctau[1][idif + 2];
                st.thet[1][ibl] = st.thet[1][idif + 2];
                st.dstr[1][ibl] = st.dstr[1][idif + 2];
                st.uedg[1][ibl] = dudx * st.xssi[1][ibl];
            }

            // move bottom side BL variables upstream
            for ibl in 2..=st.nbl[2] {
                st.ctau[2][ibl] = st.ctau[2][ibl + idif];
                st.thet[2][ibl] = st.thet[2][ibl + idif];
                st.dstr[2][ibl] = st.dstr[2][ibl + idif];
                st.uedg[2][ibl] = st.uedg[2][ibl + idif];
            }
        } else {
            // increase in number of points on bottom side (IS=2)
            let idif = istold - st.ist;
            st.itran[1] -= idif;
            st.itran[2] += idif;

            // move bottom side BL variables downstream
            for ibl in (idif + 2..=st.nbl[2]).rev() {
                st.ctau[2][ibl] = st.ctau[2][ibl - idif];
                st.thet[2][ibl] = st.thet[2][ibl - idif];
                st.dstr[2][ibl] = st.dstr[2][ibl - idif];
                st.uedg[2][ibl] = st.uedg[2][ibl - idif];
            }

            // set BL variables between old and new stagnation point
            let dudx = st.uedg[2][idif + 2] / st.xssi[2][idif + 2];
            for ibl in (2..=idif + 1).rev() {
                st.ctau[2][ibl] = st.ctau[2][idif + 2];
                st.thet[2][ibl] = st.thet[2][idif + 2];
                st.dstr[2][ibl] = st.dstr[2][idif + 2];
                st.uedg[2][ibl] = dudx * st.xssi[2][ibl];
            }

            // move top side BL variables upstream
            for ibl in 2..=st.nbl[1] {
                st.ctau[1][ibl] = st.ctau[1][ibl + idif];
                st.thet[1][ibl] = st.thet[1][ibl + idif];
                st.dstr[1][ibl] = st.dstr[1][ibl + idif];
                st.uedg[1][ibl] = st.uedg[1][ibl + idif];
            }
        }

        // tweak Ue so it's not zero, in case stag. point is right on node
        let ueps = 1.0e-7;
        for is in 1..=2 {
            for ibl in 2..=st.nbl[is] {
                let i = st.ipan[is][ibl];
                if st.uedg[is][ibl] <= ueps {
                    st.uedg[is][ibl] = ueps;
                    st.qvis[i] = st.vti[is][ibl] * ueps;
                    st.gam[i] = st.vti[is][ibl] * ueps;
                }
            }
        }
    }

    // set new mass array since Ue has been tweaked
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            st.mass[is][ibl] = st.dstr[is][ibl] * st.uedg[is][ibl];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// XSTRIP >= 1 (free transition) returns the TE arc length without touching the spline.
    #[test]
    fn xifset_free_transition_returns_te_arc_length() {
        let mut st = BlState::empty(10, 3);
        st.ist = 5;
        iblpan(&mut st);
        st.xssi[1][st.iblte[1]] = 1.25;
        st.xssi[2][st.iblte[2]] = 1.5;
        st.xstrip = [0.0, 1.0, 1.0];
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
