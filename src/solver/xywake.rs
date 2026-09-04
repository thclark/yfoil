//! XYWAKE and QWCALC (xpanel.f) and SETEXP (xutils.f): wake node coordinates traced along
//! streamlines of the current solution, and the alpha = 0, 90 tangential velocities on it.

use crate::solver::blstate::BlState;
use crate::solver::psilin::psilin;

/// SETEXP: geometrically stretched array S(1..=nn) with S(1) = 0, S(2)-S(1) = ds1, S(nn) = smax.
/// Ported character for character: the ratio is found by Newton iteration to |dRatio| < 1e-5
/// — that tolerance is load-bearing for wake node positions (plan stage S3).
pub fn setexp(ds1: f64, smax: f64, nn: usize) -> Vec<f64> {
    let sigma = smax / ds1;
    let nex = nn - 1;
    let rnex = nex as f64;
    let rni = 1.0 / rnex;

    // solve quadratic for initial geometric ratio guess
    let aaa = rnex * (rnex - 1.0) * (rnex - 2.0) / 6.0;
    let bbb = rnex * (rnex - 1.0) / 2.0;
    let ccc = rnex - sigma;

    let mut disc = bbb * bbb - 4.0 * aaa * ccc;
    disc = disc.max(0.0);

    let mut ratio = if nex <= 1 {
        panic!("SETEXP: Cannot fill array.  N too small.");
    } else if nex == 2 {
        -ccc / bbb + 1.0
    } else {
        (-bbb + disc.sqrt()) / (2.0 * aaa) + 1.0
    };

    if ratio != 1.0 {
        // Newton iteration for actual geometric ratio
        let mut converged = false;
        for _ in 0..100 {
            let sigman = (ratio.powi(nex as i32) - 1.0) / (ratio - 1.0);
            let res = sigman.powf(rni) - sigma.powf(rni);
            let dresdr =
                rni * sigman.powf(rni) * (rnex * ratio.powi(nex as i32 - 1) - sigman) / (ratio.powi(nex as i32) - 1.0);
            let dratio = -res / dresdr;
            ratio += dratio;
            if dratio.abs() < 1.0e-5 {
                converged = true;
                break;
            }
        }
        if !converged {
            // XFOIL: 'SETEXP: Convergence failed.  Continuing anyway ...'
        }
    }

    // set up stretched array using converged geometric ratio
    let mut s = vec![0.0; nn + 1];
    s[1] = 0.0;
    let mut ds = ds1;
    for k in 2..=nn {
        s[k] = s[k - 1] + ds;
        ds *= ratio;
    }
    s
}

/// XYWAKE: sets wake node coordinates X/Y/S, normals NX/NY and panel angles APANEL for
/// nodes n+1..=n+nw from the current GAM (and SIG, though SIGLIN is off here) distribution.
/// `waklen` is XFOIL's WAKLEN (chords). Requires `st.nw == n/12 + 10*INT(WAKLEN)`.
pub fn xywake(st: &mut BlState, waklen: f64) {
    let n = st.n;
    let nw = st.nw;
    debug_assert_eq!(nw, n / 12 + 10 * (waklen as usize), "NW must follow XYWAKE's formula");

    let ds1 = 0.5 * (st.s[2] - st.s[1] + st.s[n] - st.s[n - 1]);
    let snew = setexp(ds1, waklen * st.chord, nw); // SNEW(N+1..N+NW) as snew[1..=nw]

    let xte = 0.5 * (st.x[1] + st.x[n]);
    let yte = 0.5 * (st.y[1] + st.y[n]);

    // set first wake point a tiny distance behind TE
    let i = n + 1;
    let sx = 0.5 * (st.yp[n] - st.yp[1]);
    let sy = 0.5 * (st.xp[1] - st.xp[n]);
    let smod = (sx * sx + sy * sy).sqrt();
    st.nx[i] = sx / smod;
    st.ny[i] = sy / smod;
    st.x[i] = xte - 0.0001 * st.ny[i];
    st.y[i] = yte + 0.0001 * st.nx[i];
    st.s[i] = st.s[n];

    // calculate streamfunction gradient components at first point
    let psi_x = psilin(st, i, st.x[i], st.y[i], 1.0, 0.0, false).psi_ni;
    let psi_y = psilin(st, i, st.x[i], st.y[i], 0.0, 1.0, false).psi_ni;

    // set unit vector normal to wake at first point
    st.nx[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
    st.ny[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

    // set angle of wake panel normal
    st.apanel[i] = psi_y.atan2(psi_x);

    // set rest of wake points
    for i in (n + 2)..=(n + nw) {
        let ds = snew[i - n] - snew[i - n - 1];

        // set new point DS downstream of last point
        st.x[i] = st.x[i - 1] - ds * st.ny[i];
        st.y[i] = st.y[i - 1] + ds * st.nx[i];
        st.s[i] = st.s[i - 1] + ds;

        if i == n + nw {
            break;
        }

        // calculate normal vector for next point
        let psi_x = psilin(st, i, st.x[i], st.y[i], 1.0, 0.0, false).psi_ni;
        let psi_y = psilin(st, i, st.x[i], st.y[i], 0.0, 1.0, false).psi_ni;

        st.nx[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
        st.ny[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

        // set angle of wake panel normal
        st.apanel[i] = psi_y.atan2(psi_x);
    }
    // LWAKE = .TRUE., AWAKE = ALFA, LWDIJ = .FALSE. (new wake geometry invalidates the wake DIJ)
    st.lwake = true;
    st.awake = st.alfa;
    st.lwdij = false;
}

/// QWCALC: inviscid tangential velocity for alpha = 0, 90 on the wake due to freestream and
/// airfoil surface vorticity.
pub fn qwcalc(st: &mut BlState) {
    let n = st.n;
    // first wake point (same as TE)
    st.qinvu[1][n + 1] = st.qinvu[1][n];
    st.qinvu[2][n + 1] = st.qinvu[2][n];
    // rest of wake
    for i in (n + 2)..=(n + st.nw) {
        let p = psilin(st, i, st.x[i], st.y[i], st.nx[i], st.ny[i], false);
        st.qinvu[1][i] = p.qtan1;
        st.qinvu[2][i] = p.qtan2;
    }
}
