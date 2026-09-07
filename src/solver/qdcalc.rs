//! PSWLIN and QDCALC (xpanel.f): wake-source streamfunction sensitivities and the full
//! (N+NW)×(N+NW) source influence matrix DIJ = dQtan/dSig.

use crate::solver::blstate::SolverState;
use crate::solver::ggcalc::InviscidSystem;
use crate::solver::ludcmp::lu_back_substitute;
use crate::solver::psilin::{panel_influence, pi_consts};

/// What PSWLIN returns: Psi, dPsi/dn, and the wake-source sensitivities (indexed n+1..=n+nw).
#[derive(Debug, Clone)]
pub struct WakeSourceInfluence {
    pub psi: f64,
    pub psi_d_n: f64,
    pub psi_d_sigma: Vec<f64>,
    pub qtan_d_sigma: Vec<f64>,
}

/// PSWLIN(I, XI, YI, NXI, NYI, PSI, PSI_NI): streamfunction at node I due to the wake
/// sources. Note the branch-cut correction is `- (0.5-0.5*SGN)*PI` here (PSILIN has `+`).
#[doc(alias = "PSWLIN")]
pub fn wake_source_influence(
    state: &SolverState,
    i: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
) -> WakeSourceInfluence {
    let n = state.n_foil_nodes;
    let nw = state.n_wake_nodes;
    let np = n + nw;
    let (pi, _hopi, qopi) = pi_consts();
    let (x, y) = (&state.x, &state.y);
    let io = i;

    let mut out = WakeSourceInfluence {
        psi: 0.0,
        psi_d_n: 0.0,
        psi_d_sigma: vec![0.0; np + 1],
        qtan_d_sigma: vec![0.0; np + 1],
    };

    for jo in (n + 1)..=(n + nw - 1) {
        let jp = jo + 1;
        let mut jm = jo - 1;
        let mut jq = jp + 1;
        if jo == n + 1 {
            jm = jo;
        } else if jo == n + nw - 1 {
            jq = jp;
        }

        let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
        let dsio = 1.0 / dso;
        let apan = state.panel_angle[jo];

        let rx1 = xi - x[jo];
        let ry1 = yi - y[jo];
        let rx2 = xi - x[jp];
        let ry2 = yi - y[jp];

        let sx = (x[jp] - x[jo]) * dsio;
        let sy = (y[jp] - y[jo]) * dsio;

        let x1 = sx * rx1 + sy * ry1;
        let x2 = sx * rx2 + sy * ry2;
        let yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        let sgn = if io > n && io <= n + nw {
            1.0
        } else {
            1.0_f64.copysign(yy)
        };

        let (g1, t1) = if io != jo && rs1 > 0.0 {
            (rs1.ln(), (sgn * x1).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * pi)
        } else {
            (0.0, 0.0)
        };
        let (g2, t2) = if io != jp && rs2 > 0.0 {
            (rs2.ln(), (sgn * x2).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * pi)
        } else {
            (0.0, 0.0)
        };

        let x1i = sx * nxi + sy * nyi;
        let x2i = sx * nxi + sy * nyi;
        let yyi = sx * nyi - sy * nxi;

        // set up midpoint quantities
        let x0 = 0.5 * (x1 + x2);
        let rs0 = x0 * x0 + yy * yy;
        let g0 = rs0.ln();
        let t0 = (sgn * x0).atan2(sgn * yy) - (0.5 - 0.5 * sgn) * pi;

        // calculate source contribution to Psi for 1-0 half-panel
        let dxinv = 1.0 / (x1 - x0);
        let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
        let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan) + (x0 - x1) * yy) * dxinv;

        let psx1 = -(t1 - apan);
        let psx0 = t0 - apan;
        let psyy = 0.5 * (g1 - g0);

        let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
        let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
        let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;

        let dsm = ((x[jp] - x[jm]).powi(2) + (y[jp] - y[jm]).powi(2)).sqrt();
        let dsim = 1.0 / dsm;

        let ssum = (state.sigma[jp] - state.sigma[jo]) * dsio + (state.sigma[jp] - state.sigma[jm]) * dsim;
        let sdif = (state.sigma[jp] - state.sigma[jo]) * dsio - (state.sigma[jp] - state.sigma[jm]) * dsim;

        out.psi += qopi * (psum * ssum + pdif * sdif);

        // dPsi/dm
        out.psi_d_sigma[jm] += qopi * (-psum * dsim + pdif * dsim);
        out.psi_d_sigma[jo] += qopi * (-psum * dsio - pdif * dsio);
        out.psi_d_sigma[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

        // dPsi/dni
        let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
        let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
        out.psi_d_n += qopi * (psni * ssum + pdni * sdif);
        out.qtan_d_sigma[jm] += qopi * (-psni * dsim + pdni * dsim);
        out.qtan_d_sigma[jo] += qopi * (-psni * dsio - pdni * dsio);
        out.qtan_d_sigma[jp] += qopi * (psni * (dsio + dsim) + pdni * (dsio - dsim));

        // calculate source contribution to Psi for 0-2 half-panel
        let dxinv = 1.0 / (x0 - x2);
        let psum = x2 * (t2 - apan) - x0 * (t0 - apan) + 0.5 * yy * (g0 - g2);
        let pdif = ((x0 + x2) * psum + rs0 * (t0 - apan) - rs2 * (t2 - apan) + (x2 - x0) * yy) * dxinv;

        let psx0 = -(t0 - apan);
        let psx2 = t2 - apan;
        let psyy = 0.5 * (g0 - g2);

        let pdx0 = ((x0 + x2) * psx0 + psum + 2.0 * x0 * (t0 - apan) - pdif) * dxinv;
        let pdx2 = ((x0 + x2) * psx2 + psum - 2.0 * x2 * (t2 - apan) + pdif) * dxinv;
        let pdyy = ((x0 + x2) * psyy + 2.0 * (x2 - x0 + yy * (t0 - t2))) * dxinv;

        let dsp = ((x[jq] - x[jo]).powi(2) + (y[jq] - y[jo]).powi(2)).sqrt();
        let dsip = 1.0 / dsp;

        let ssum = (state.sigma[jq] - state.sigma[jo]) * dsip + (state.sigma[jp] - state.sigma[jo]) * dsio;
        let sdif = (state.sigma[jq] - state.sigma[jo]) * dsip - (state.sigma[jp] - state.sigma[jo]) * dsio;

        out.psi += qopi * (psum * ssum + pdif * sdif);

        // dPsi/dm
        out.psi_d_sigma[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
        out.psi_d_sigma[jp] += qopi * (psum * dsio - pdif * dsio);
        out.psi_d_sigma[jq] += qopi * (psum * dsip + pdif * dsip);

        // dPsi/dni
        let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
        let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
        out.psi_d_n += qopi * (psni * ssum + pdni * sdif);
        out.qtan_d_sigma[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
        out.qtan_d_sigma[jp] += qopi * (psni * dsio - pdni * dsio);
        out.qtan_d_sigma[jq] += qopi * (psni * dsip + pdni * dsip);
    }
    out
}

/// QDCALC: source panel influence coefficient matrix for the current airfoil and wake
/// geometry, stored 1-based in `st.dij[i][j]` for i, j in 1..=N+NW.
#[doc(alias = "QDCALC")]
pub fn build_dij(state: &mut SolverState, sys: &mut InviscidSystem) {
    let n = state.n_foil_nodes;
    let nw = state.n_wake_nodes;
    let np = n + nw;
    // DIJ persists in COMMON: the airfoil block (1..N, 1..N) is computed once (LADIJ) and only the
    // wake rows/columns are refreshed here when the wake moves (LWDIJ)
    if state.dij.len() != np + 1 || state.dij.iter().any(|r| r.len() != np + 1) {
        state.dij = vec![vec![0.0; np + 1]; np + 1];
    }

    if !sys.dij_foil_built {
        // source influence matrix for airfoil surface: multiply each dPsi/dSig vector by the
        // inverse of the factored dPsi/dGam matrix
        for j in 1..=n {
            let mut col: Vec<f64> = (0..=(n + 1)).map(|i| sys.bij[i][j]).collect();
            lu_back_substitute(&sys.aij_lu, &mut col);
            for i in 1..=(n + 1) {
                sys.bij[i][j] = col[i];
            }
            // store resulting dGam/dSig = dQtan/dSig vector
            for i in 1..=n {
                state.dij[i][j] = sys.bij[i][j];
            }
        }
        sys.dij_foil_built = true;
    }

    // set up coefficient matrix of dPsi/dm on airfoil surface
    for i in 1..=n {
        let p = wake_source_influence(state, i, state.x[i], state.y[i], state.normal_x[i], state.normal_y[i]);
        for j in (n + 1)..=np {
            sys.bij[i][j] = -p.psi_d_sigma[j];
        }
    }

    // Kutta condition (no direct source influence)
    for j in (n + 1)..=np {
        sys.bij[n + 1][j] = 0.0;
    }
    // sharp TE gamma extrapolation also has no source influence
    if state.sharp_te {
        for j in (n + 1)..=np {
            sys.bij[n][j] = 0.0;
        }
    }

    // multiply by inverse of factored dPsi/dGam matrix
    for j in (n + 1)..=np {
        let mut col: Vec<f64> = (0..=(n + 1)).map(|i| sys.bij[i][j]).collect();
        lu_back_substitute(&sys.aij_lu, &mut col);
        for i in 1..=(n + 1) {
            sys.bij[i][j] = col[i];
        }
    }

    // set the source influence matrix for the wake sources
    for i in 1..=n {
        for j in (n + 1)..=np {
            state.dij[i][j] = sys.bij[i][j];
        }
    }

    // influence of sources on the wake velocities: dQtan/dGam and dQtan/dSig at wake points
    let mut cij = vec![vec![0.0; n + 1]; nw + 1];
    for i in (n + 1)..=np {
        let i_wake = i - n;
        // airfoil contribution at wake panel node
        let p = panel_influence(
            state,
            i,
            state.x[i],
            state.y[i],
            state.normal_x[i],
            state.normal_y[i],
            true,
        );
        for j in 1..=n {
            cij[i_wake][j] = p.qtan_d_gamma[j];
        }
        for j in 1..=n {
            state.dij[i][j] = p.qtan_d_sigma[j];
        }
        // wake contribution
        let w = wake_source_influence(state, i, state.x[i], state.y[i], state.normal_x[i], state.normal_y[i]);
        for j in (n + 1)..=np {
            state.dij[i][j] = w.qtan_d_sigma[j];
        }
    }

    // add on effect of all sources on airfoil vorticity which affects wake Qtan
    for i in (n + 1)..=np {
        let i_wake = i - n;
        // airfoil surface source contribution first
        for j in 1..=n {
            let mut sum = 0.0;
            for k in 1..=n {
                sum += cij[i_wake][k] * state.dij[k][j];
            }
            state.dij[i][j] += sum;
        }
        // wake source contribution next
        for j in (n + 1)..=np {
            let mut sum = 0.0;
            for k in 1..=n {
                sum += cij[i_wake][k] * sys.bij[k][j];
            }
            state.dij[i][j] += sum;
        }
    }

    // make sure first wake point has same velocity as trailing edge
    for j in 1..=np {
        state.dij[n + 1][j] = state.dij[n][j];
    }
    // LWDIJ = .TRUE.
    state.dij_wake_built = true;
}
