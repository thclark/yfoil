//! PSILIN (xpanel.f): streamfunction at a panel or wake node due to freestream, all bound
//! vorticity and (optionally) the viscous source distribution, with the sensitivity vectors
//! XFOIL returns through COMMON: dPsi/dGam (DZDG), dPsi/dSig (DZDM), and the tangential
//! velocity sensitivities DQDG/DQDM, plus QTAN1/QTAN2 (alpha = 0, 90 tangential velocities).
//!
//! Not translated: the GEOLIN branch (geometric sensitivities for inverse design: DZDN,
//! Z_QDOF*) and the LIMAGE ground-image branch — neither is reachable in analysis mode.

use crate::solver::blstate::SolverState;

/// Everything PSILIN leaves in COMMON that analysis mode reads.
#[derive(Debug, Clone)]
pub struct PanelInfluence {
    pub psi: f64,
    pub psi_d_n: f64,
    pub qtan_alpha0: f64,
    pub qtan_alpha90: f64,
    pub qtan_sigma: f64,
    pub psi_d_qinf: f64,
    pub psi_d_alpha: f64,
    /// 1-based, 1..=n
    pub psi_d_gamma: Vec<f64>,
    pub qtan_d_gamma: Vec<f64>,
    pub psi_d_sigma: Vec<f64>,
    pub qtan_d_sigma: Vec<f64>,
}

/// XFOIL's INIT: PI = 4.0*ATAN(1.0), HOPI = 0.5/PI, QOPI = 0.25/PI — computed the same way.
pub fn pi_consts() -> (f64, f64, f64) {
    let pi = 4.0 * (1.0_f64).atan();
    (pi, 0.50 / pi, 0.25 / pi)
}

/// PSILIN(I, XI, YI, NXI, NYI, PSI, PSI_NI, GEOLIN=.FALSE., SIGLIN).
/// `i` is the 1-based node index the point belongs to (airfoil 1..=n, wake n+1..=n+nw);
/// it only affects the self-influence skips and the arctan reflection flag.
#[doc(alias = "PSILIN")]
pub fn panel_influence(
    st: &SolverState,
    i: usize,
    xi: f64,
    yi: f64,
    nxi: f64,
    nyi: f64,
    siglin: bool,
) -> PanelInfluence {
    let n = st.n_foil_nodes;
    let (pi, hopi, qopi) = pi_consts();
    let (x, y, s) = (&st.x, &st.y, &st.s);
    let gamu1 = &st.q_inviscid_basis[1];
    let gamu2 = &st.q_inviscid_basis[2];

    // distance tolerance for determining if two points are the same
    let seps = (s[n] - s[1]) * 1.0e-5;
    let io = i;
    let cosa = st.alpha.cos();
    let sina = st.alpha.sin();

    let mut out = PanelInfluence {
        psi: 0.0,
        psi_d_n: 0.0,
        qtan_alpha0: 0.0,
        qtan_alpha90: 0.0,
        qtan_sigma: 0.0,
        psi_d_qinf: 0.0,
        psi_d_alpha: 0.0,
        psi_d_gamma: vec![0.0; n + 1],
        qtan_d_gamma: vec![0.0; n + 1],
        psi_d_sigma: vec![0.0; n + 1],
        qtan_d_sigma: vec![0.0; n + 1],
    };

    let (scs, sds) = if st.sharp_te {
        (1.0, 0.0)
    } else {
        (st.te_thickness_normal / st.te_gap, st.te_thickness_parallel / st.te_gap)
    };

    // carried out of the loop for the TE panel (labels 11/12)
    let (mut x1, mut x2, mut yy) = (0.0, 0.0, 0.0);
    let (mut g1, mut g2, mut t1, mut t2) = (0.0, 0.0, 0.0, 0.0);
    let mut apan = 0.0;
    let (mut x1i, mut x2i, mut yyi) = (0.0, 0.0, 0.0);
    let (mut jo, mut jp) = (n, 1usize);
    let mut te_closed = false; // GO TO 12: sharp/closed TE, skip the TE panel entirely

    for jo_ in 1..=n {
        jo = jo_;
        jp = jo + 1;
        let mut jm = jo.wrapping_sub(1);
        let mut jq = jp + 1;
        if jo == 1 {
            jm = jo;
        } else if jo == n - 1 {
            jq = jp;
        } else if jo == n {
            jp = 1;
            if (x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2) < seps * seps {
                te_closed = true;
                break;
            }
        }

        let dso = ((x[jo] - x[jp]).powi(2) + (y[jo] - y[jp]).powi(2)).sqrt();
        // skip null panel
        if dso == 0.0 {
            continue;
        }
        let dsio = 1.0 / dso;
        apan = st.panel_angle[jo];

        let rx1 = xi - x[jo];
        let ry1 = yi - y[jo];
        let rx2 = xi - x[jp];
        let ry2 = yi - y[jp];

        let sx = (x[jp] - x[jo]) * dsio;
        let sy = (y[jp] - y[jo]) * dsio;

        x1 = sx * rx1 + sy * ry1;
        x2 = sx * rx2 + sy * ry2;
        yy = sx * ry1 - sy * rx1;

        let rs1 = rx1 * rx1 + ry1 * ry1;
        let rs2 = rx2 * rx2 + ry2 * ry2;

        // set reflection flag SGN to avoid branch problems with arctan
        let sgn = if io >= 1 && io <= n {
            1.0 // no problem on airfoil surface
        } else {
            1.0_f64.copysign(yy) // make sure arctan falls between -/+ Pi/2
        };

        // set log(r^2) and arctan(x/y), correcting for reflection if any
        if io != jo && rs1 > 0.0 {
            g1 = rs1.ln();
            t1 = (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * pi;
        } else {
            g1 = 0.0;
            t1 = 0.0;
        }
        if io != jp && rs2 > 0.0 {
            g2 = rs2.ln();
            t2 = (sgn * x2).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * pi;
        } else {
            g2 = 0.0;
            t2 = 0.0;
        }

        x1i = sx * nxi + sy * nyi;
        x2i = sx * nxi + sy * nyi;
        yyi = sx * nyi - sy * nxi;

        if jo == n {
            break; // GO TO 11: TE panel handled below
        }

        if siglin {
            // set up midpoint quantities
            let x0 = 0.5 * (x1 + x2);
            let rs0 = x0 * x0 + yy * yy;
            let g0 = rs0.ln();
            let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * pi;

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

            let ssum = (st.sigma[jp] - st.sigma[jo]) * dsio + (st.sigma[jp] - st.sigma[jm]) * dsim;
            let sdif = (st.sigma[jp] - st.sigma[jo]) * dsio - (st.sigma[jp] - st.sigma[jm]) * dsim;

            out.psi += qopi * (psum * ssum + pdif * sdif);

            // dPsi/dm
            out.psi_d_sigma[jm] += qopi * (-psum * dsim + pdif * dsim);
            out.psi_d_sigma[jo] += qopi * (-psum * dsio - pdif * dsio);
            out.psi_d_sigma[jp] += qopi * (psum * (dsio + dsim) + pdif * (dsio - dsim));

            // dPsi/dni
            let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
            let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
            out.psi_d_n += qopi * (psni * ssum + pdni * sdif);
            out.qtan_sigma += qopi * (psni * ssum + pdni * sdif);

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

            let ssum = (st.sigma[jq] - st.sigma[jo]) * dsip + (st.sigma[jp] - st.sigma[jo]) * dsio;
            let sdif = (st.sigma[jq] - st.sigma[jo]) * dsip - (st.sigma[jp] - st.sigma[jo]) * dsio;

            out.psi += qopi * (psum * ssum + pdif * sdif);

            // dPsi/dm
            out.psi_d_sigma[jo] += qopi * (-psum * (dsip + dsio) - pdif * (dsip - dsio));
            out.psi_d_sigma[jp] += qopi * (psum * dsio - pdif * dsio);
            out.psi_d_sigma[jq] += qopi * (psum * dsip + pdif * dsip);

            // dPsi/dni
            let psni = psx0 * (x1i + x2i) * 0.5 + psx2 * x2i + psyy * yyi;
            let pdni = pdx0 * (x1i + x2i) * 0.5 + pdx2 * x2i + pdyy * yyi;
            out.psi_d_n += qopi * (psni * ssum + pdni * sdif);
            out.qtan_sigma += qopi * (psni * ssum + pdni * sdif);

            out.qtan_d_sigma[jo] += qopi * (-psni * (dsip + dsio) - pdni * (dsip - dsio));
            out.qtan_d_sigma[jp] += qopi * (psni * dsio - pdni * dsio);
            out.qtan_d_sigma[jq] += qopi * (psni * dsip + pdni * dsip);
        }

        // calculate vortex panel contribution to Psi
        let dxinv = 1.0 / (x1 - x2);
        let psis = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);
        let psid = ((x1 + x2) * psis + 0.5 * (rs2 * g2 - rs1 * g1 + x1 * x1 - x2 * x2)) * dxinv;

        let psx1 = 0.5 * g1;
        let psx2 = -0.5 * g2;
        let psyy = t1 - t2;

        let pdx1 = ((x1 + x2) * psx1 + psis - x1 * g1 - psid) * dxinv;
        let pdx2 = ((x1 + x2) * psx2 + psis + x2 * g2 + psid) * dxinv;
        let pdyy = ((x1 + x2) * psyy - yy * (g1 - g2)) * dxinv;

        let gsum1 = gamu1[jp] + gamu1[jo];
        let gsum2 = gamu2[jp] + gamu2[jo];
        let gdif1 = gamu1[jp] - gamu1[jo];
        let gdif2 = gamu2[jp] - gamu2[jo];

        let gsum = st.gamma[jp] + st.gamma[jo];
        let gdif = st.gamma[jp] - st.gamma[jo];

        out.psi += qopi * (psis * gsum + psid * gdif);

        // dPsi/dGam
        out.psi_d_gamma[jo] += qopi * (psis - psid);
        out.psi_d_gamma[jp] += qopi * (psis + psid);

        // dPsi/dni
        let psni = psx1 * x1i + psx2 * x2i + psyy * yyi;
        let pdni = pdx1 * x1i + pdx2 * x2i + pdyy * yyi;
        out.psi_d_n += qopi * (gsum * psni + gdif * pdni);
        out.qtan_alpha0 += qopi * (gsum1 * psni + gdif1 * pdni);
        out.qtan_alpha90 += qopi * (gsum2 * psni + gdif2 * pdni);

        out.qtan_d_gamma[jo] += qopi * (psni - pdni);
        out.qtan_d_gamma[jp] += qopi * (psni + pdni);
    }

    if !te_closed {
        // label 11: TE panel (JO = N, JP = 1) using the last loop's X1/X2/YY/G/T/APAN
        let psig = 0.5 * yy * (g1 - g2) + x2 * (t2 - apan) - x1 * (t1 - apan);
        let pgam = 0.5 * x1 * g1 - 0.5 * x2 * g2 + x2 - x1 + yy * (t1 - t2);

        let psigx1 = -(t1 - apan);
        let psigx2 = t2 - apan;
        let psigyy = 0.5 * (g1 - g2);
        let pgamx1 = 0.5 * g1;
        let pgamx2 = -0.5 * g2;
        let pgamyy = t1 - t2;

        let psigni = psigx1 * x1i + psigx2 * x2i + psigyy * yyi;
        let pgamni = pgamx1 * x1i + pgamx2 * x2i + pgamyy * yyi;

        // TE panel source and vortex strengths
        let sigte1 = 0.5 * scs * (gamu1[jp] - gamu1[jo]);
        let sigte2 = 0.5 * scs * (gamu2[jp] - gamu2[jo]);
        let gamte1 = -0.5 * sds * (gamu1[jp] - gamu1[jo]);
        let gamte2 = -0.5 * sds * (gamu2[jp] - gamu2[jo]);

        let sigte = 0.5 * scs * (st.gamma[jp] - st.gamma[jo]);
        let gamte = -0.5 * sds * (st.gamma[jp] - st.gamma[jo]);

        // TE panel contribution to Psi
        out.psi += hopi * (psig * sigte + pgam * gamte);

        // dPsi/dGam
        out.psi_d_gamma[jo] -= hopi * psig * scs * 0.5;
        out.psi_d_gamma[jp] += hopi * psig * scs * 0.5;
        out.psi_d_gamma[jo] += hopi * pgam * sds * 0.5;
        out.psi_d_gamma[jp] -= hopi * pgam * sds * 0.5;

        // dPsi/dni
        out.psi_d_n += hopi * (psigni * sigte + pgamni * gamte);
        out.qtan_alpha0 += hopi * (psigni * sigte1 + pgamni * gamte1);
        out.qtan_alpha90 += hopi * (psigni * sigte2 + pgamni * gamte2);

        out.qtan_d_gamma[jo] -= hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
        out.qtan_d_gamma[jp] += hopi * (psigni * 0.5 * scs - pgamni * 0.5 * sds);
    }

    // label 12: freestream terms
    out.psi += st.qinf * (cosa * yi - sina * xi);
    out.psi_d_n += st.qinf * (cosa * nyi - sina * nxi);
    out.qtan_alpha0 += st.qinf * nyi;
    out.qtan_alpha90 -= st.qinf * nxi;
    out.psi_d_qinf += cosa * yi - sina * xi;
    out.psi_d_alpha -= st.qinf * (sina * yi + cosa * xi);

    out
}
