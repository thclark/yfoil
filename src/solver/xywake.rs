//! XYWAKE and QWCALC (xpanel.f) and SETEXP (xutils.f): wake node coordinates traced along
//! streamlines of the current solution, and the alpha = 0, 90 tangential velocities on it.

use crate::solver::blstate::SolverState;
use crate::solver::psilin::panel_influence;

/// SETEXP: geometrically stretched array S(1..=nn) with S(1) = 0, S(2)-S(1) = ds1, S(nn) = smax.
/// Ported character for character: the ratio is found by Newton iteration to |dRatio| < 1e-5
/// — that tolerance is load-bearing for wake node positions (plan stage S3).
#[doc(alias = "SETEXP")]
pub fn exponential_spacing(ds1: f64, smax: f64, nn: usize) -> Vec<f64> {
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
#[doc(alias = "XYWAKE")]
pub fn build_wake(st: &mut SolverState, waklen: f64) {
    let n = st.n_foil_nodes;
    let nw = st.n_wake_nodes;
    debug_assert_eq!(nw, n / 12 + 10 * (waklen as usize), "NW must follow XYWAKE's formula");

    let ds1 = 0.5 * (st.s[2] - st.s[1] + st.s[n] - st.s[n - 1]);
    let snew = exponential_spacing(ds1, waklen * st.chord, nw); // SNEW(N+1..N+NW) as snew[1..=nw]

    let xte = 0.5 * (st.x[1] + st.x[n]);
    let yte = 0.5 * (st.y[1] + st.y[n]);

    // set first wake point a tiny distance behind TE
    let i = n + 1;
    let sx = 0.5 * (st.dyds[n] - st.dyds[1]);
    let sy = 0.5 * (st.dxds[1] - st.dxds[n]);
    let smod = (sx * sx + sy * sy).sqrt();
    st.normal_x[i] = sx / smod;
    st.normal_y[i] = sy / smod;
    st.x[i] = xte - 0.0001 * st.normal_y[i];
    st.y[i] = yte + 0.0001 * st.normal_x[i];
    st.s[i] = st.s[n];

    // calculate streamfunction gradient components at first point
    let psi_x = panel_influence(st, i, st.x[i], st.y[i], 1.0, 0.0, false).psi_d_n;
    let psi_y = panel_influence(st, i, st.x[i], st.y[i], 0.0, 1.0, false).psi_d_n;

    // set unit vector normal to wake at first point
    st.normal_x[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
    st.normal_y[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

    // set angle of wake panel normal
    st.panel_angle[i] = psi_y.atan2(psi_x);

    // set rest of wake points
    for i in (n + 2)..=(n + nw) {
        let ds = snew[i - n] - snew[i - n - 1];

        // set new point DS downstream of last point
        st.x[i] = st.x[i - 1] - ds * st.normal_y[i];
        st.y[i] = st.y[i - 1] + ds * st.normal_x[i];
        st.s[i] = st.s[i - 1] + ds;

        if i == n + nw {
            break;
        }

        // calculate normal vector for next point
        let psi_x = panel_influence(st, i, st.x[i], st.y[i], 1.0, 0.0, false).psi_d_n;
        let psi_y = panel_influence(st, i, st.x[i], st.y[i], 0.0, 1.0, false).psi_d_n;

        st.normal_x[i + 1] = -psi_x / (psi_x * psi_x + psi_y * psi_y).sqrt();
        st.normal_y[i + 1] = -psi_y / (psi_x * psi_x + psi_y * psi_y).sqrt();

        // set angle of wake panel normal
        st.panel_angle[i] = psi_y.atan2(psi_x);
    }
    // LWAKE = .TRUE., AWAKE = ALFA, LWDIJ = .FALSE. (new wake geometry invalidates the wake DIJ)
    st.wake_built = true;
    st.alpha_wake = st.alpha;
    st.dij_wake_built = false;
}

/// QWCALC: inviscid tangential velocity for alpha = 0, 90 on the wake due to freestream and
/// airfoil surface vorticity.
#[doc(alias = "QWCALC")]
pub fn set_wake_q_basis(st: &mut SolverState) {
    let n = st.n_foil_nodes;
    // first wake point (same as TE)
    st.q_inviscid_basis[1][n + 1] = st.q_inviscid_basis[1][n];
    st.q_inviscid_basis[2][n + 1] = st.q_inviscid_basis[2][n];
    // rest of wake
    for i in (n + 2)..=(n + st.n_wake_nodes) {
        let p = panel_influence(st, i, st.x[i], st.y[i], st.normal_x[i], st.normal_y[i], false);
        st.q_inviscid_basis[1][i] = p.qtan_alpha0;
        st.q_inviscid_basis[2][i] = p.qtan_alpha90;
    }
}
