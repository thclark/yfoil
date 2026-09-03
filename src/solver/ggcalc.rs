//! GGCALC (xpanel.f): the two unit vorticity distributions (alpha = 0, 90 deg) from the
//! Psi = Psio system, with the Kutta condition and the sharp-TE bisector condition; keeps the
//! factored AIJ and the source-influence matrix BIJ that QDCALC needs.

use crate::solver::blstate::BlState;
use crate::solver::ludcmp::{baksub, ludcmp, LuFactors};
use crate::solver::psilin::psilin;

/// The linear system GGCALC leaves behind: factored dPsi/dGam (AIJ, (N+1)×(N+1)) and
/// dPsi/dSig (BIJ, (N+1)×(N+NW); wake columns filled by QDCALC). 1-based.
#[derive(Debug, Clone)]
pub struct InviscidSystem {
    pub aij: LuFactors,
    pub bij: Vec<Vec<f64>>,
    /// LADIJ: airfoil block of DIJ already computed
    pub ladij: bool,
}

/// ATANC (xutils.f): ATAN2 with branch-cut checking against a previous angle.
pub fn atanc(y: f64, x: f64, thold: f64) -> f64 {
    let pi: f64 = "3.1415926535897932384".parse().unwrap();
    let tpi: f64 = "6.2831853071795864769".parse().unwrap();
    let thnew = y.atan2(x);
    let dthet = thnew - thold;
    let dtcorr = dthet - tpi * ((dthet + pi.copysign(dthet)) / tpi).trunc();
    thold + dtcorr
}

/// GGCALC. Sets `st.gam = 0`, `st.qinvu[1..=2][1..=n]`, and returns the factored system.
pub fn ggcalc(st: &mut BlState) -> InviscidSystem {
    let n = st.n;
    let np = n + st.nw;
    // distance of internal control point ahead of sharp TE (fraction of smaller panel length)
    let bwt = 0.1;

    for i in 1..=n {
        st.gam[i] = 0.0;
        st.qinvu[1][i] = 0.0;
        st.qinvu[2][i] = 0.0;
    }

    let mut aij = vec![vec![0.0; n + 2]; n + 2];
    let mut bij = vec![vec![0.0; np + 1]; n + 2];
    let mut gamu1 = vec![0.0; n + 2];
    let mut gamu2 = vec![0.0; n + 2];

    // Set up matrix system for Psi = Psio on airfoil surface; unknowns (dGamma)i and dPsio.
    for i in 1..=n {
        let p = psilin(st, i, st.x[i], st.y[i], st.nx[i], st.ny[i], true);
        // RES1 = PSI( 0) - PSIO,  RES2 = PSI(90) - PSIO
        let res1 = st.qinf * st.y[i];
        let res2 = -st.qinf * st.x[i];
        for j in 1..=n {
            aij[i][j] = p.dzdg[j];
        }
        for j in 1..=n {
            bij[i][j] = -p.dzdm[j];
        }
        aij[i][n + 1] = -1.0;
        gamu1[i] = -res1;
        gamu2[i] = -res2;
    }

    // Kutta condition: RES = GAM(1) + GAM(N)
    let res = 0.0;
    for j in 1..=(n + 1) {
        aij[n + 1][j] = 0.0;
    }
    aij[n + 1][1] = 1.0;
    aij[n + 1][n] = 1.0;
    gamu1[n + 1] = -res;
    gamu2[n + 1] = -res;
    // no direct source influence on the Kutta condition
    for j in 1..=n {
        bij[n + 1][j] = 0.0;
    }

    if st.sharp {
        // zero internal velocity in TE corner: TE bisector angle
        let ag1 = (-st.yp[1]).atan2(-st.xp[1]);
        let ag2 = atanc(st.yp[n], st.xp[n], ag1);
        let abis = 0.5 * (ag1 + ag2);
        let cbis = abis.cos();
        let sbis = abis.sin();
        // minimum panel length adjacent to TE
        let ds1 = ((st.x[1] - st.x[2]).powi(2) + (st.y[1] - st.y[2]).powi(2)).sqrt();
        let ds2 = ((st.x[n] - st.x[n - 1]).powi(2) + (st.y[n] - st.y[n - 1]).powi(2)).sqrt();
        let dsmin = ds1.min(ds2);
        // control point on bisector just ahead of TE point
        let xbis = st.xte - bwt * dsmin * cbis;
        let ybis = st.yte - bwt * dsmin * sbis;
        // velocity component along bisector line (I = 0: off-surface point)
        let p = psilin(st, 0, xbis, ybis, -sbis, cbis, true);
        for j in 1..=n {
            aij[n][j] = p.dqdg[j];
        }
        for j in 1..=n {
            bij[n][j] = -p.dqdm[j];
        }
        aij[n][n + 1] = 0.0;
        gamu1[n] = -cbis;
        gamu2[n] = -sbis;
    }

    // LU-factor coefficient matrix AIJ and solve for the two vorticity distributions
    let lu = ludcmp(n + 1, aij);
    baksub(&lu, &mut gamu1);
    baksub(&lu, &mut gamu2);

    // inviscid alpha=0,90 surface speeds for this geometry
    st.qinvu[1][1..=n].copy_from_slice(&gamu1[1..=n]);
    st.qinvu[2][1..=n].copy_from_slice(&gamu2[1..=n]);
    // GAMU(N+1) is PSIO for each solution; keep it with the airfoil arrays via qinvu's slot n+1
    // only if there is no wake node there — the wake overwrites it in QWCALC as XFOIL does too.
    InviscidSystem {
        aij: lu,
        bij,
        ladij: false,
    }
}
