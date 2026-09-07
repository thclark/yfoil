//! UPDATE (xbl.f): adds the Newton deltas to the boundary layer variables, checks for
//! excessive changes and under-relaxes if necessary, calculates the max and rms changes, and
//! the change in the global variable "AC" (CL if LALFA, alpha otherwise). Line-for-line on
//! `BlState`.
//!
//! XFOIL aliases `UNEW/U_AC` onto `VA` and `QNEW/Q_AC` onto `VB` with EQUIVALENCE; here they
//! are plain locals — BLSOLV's input was consumed by `blsolv`, so nothing can read VA/VB after
//! the solve by construction.

use crate::bl::system::limit_dstar;
use crate::solver::blstate::BlState;

/// What UPDATE reports besides the arrays it writes into `BlState`.
#[derive(Debug, Clone, Default)]
pub struct UpdateSummary {
    pub relaxation: f64,
    pub residual: f64,
    /// Largest normalised change (signed; for Ue XFOIL stores the raw DUEDG)
    pub residual_max: f64,
    /// 'n' (amplification), 'C' (Ctau), 'T' (theta), 'D' (delta*) or 'U' (Ue)
    pub residual_max_variable: char,
    pub i_residual_max_station: usize,
    pub residual_max_side: usize,
    /// Change in the global variable AC (CL or alpha) before under-relaxation
    pub free_variable_change: f64,
    pub cl_new: f64,
    pub cl_d_alpha: f64,
    pub cl_d_machsqd: f64,
    pub cl_d_free: f64,
}

/// UPDATE. `vdel[iv-1][k][0..2]` is BLSOLV's solution (residual and AC-sensitivity columns);
/// `minf_cl` is MINF_CL, d(MINF)/d(CL) from the last MRCL (0 for MATYP = 1).
pub fn apply_newton_update(st: &mut BlState, vdel: &[[[f64; 2]; 3]], minf_cl: f64) -> UpdateSummary {
    let pi = 4.0 * (1.0_f64).atan();
    let dtor = pi / 180.0;
    let gamm1 = st.gamma - 1.0;

    // max allowable alpha changes per iteration
    let dalmax = 0.5 * dtor;
    let dalmin = -0.5 * dtor;

    // max allowable CL change per iteration
    let dclmax = 0.5;
    let mut dclmin = -0.5;
    if st.matyp != 1 {
        dclmin = (-0.5_f64).max(-0.9 * st.cl);
    }

    let hstinv = gamm1 * ((st.minf / st.qinf) * (st.minf / st.qinf)) / (1.0 + 0.5 * gamm1 * (st.minf * st.minf));

    let nmax = st.uedg[1].len().max(st.uedg[2].len());
    let mut unew: [Vec<f64>; 3] = [Vec::new(), vec![0.0; nmax], vec![0.0; nmax]];
    let mut u_ac: [Vec<f64>; 3] = [Vec::new(), vec![0.0; nmax], vec![0.0; nmax]];

    // calculate new Ue distribution assuming no under-relaxation
    // also set the sensitivity of Ue wrt to alpha or Re
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            let i = st.ipan[is][ibl];
            let mut dui = 0.0;
            let mut dui_ac = 0.0;
            for js in 1..=2 {
                for jbl in 2..=st.nbl[js] {
                    let j = st.ipan[js][jbl];
                    let jv = st.isys[js][jbl];
                    let ue_m = -st.vti[is][ibl] * st.vti[js][jbl] * st.dij[i][j];
                    dui += ue_m * (st.mass[js][jbl] + vdel[jv - 1][2][0]);
                    dui_ac += ue_m * (-vdel[jv - 1][2][1]);
                }
            }
            // UINV depends on "AC" only if "AC" is alpha
            let uinv_ac = if st.lalfa { 0.0 } else { st.uinv_a[is][ibl] };
            unew[is][ibl] = st.uinv[is][ibl] + dui;
            u_ac[is][ibl] = uinv_ac + dui_ac;
        }
    }

    // set new Qtan from new Ue with appropriate sign change
    let n = st.n;
    let mut qnew = vec![0.0; n + 2];
    let mut q_ac = vec![0.0; n + 2];
    for is in 1..=2 {
        for ibl in 2..=st.iblte[is] {
            let i = st.ipan[is][ibl];
            qnew[i] = st.vti[is][ibl] * unew[is][ibl];
            q_ac[i] = st.vti[is][ibl] * u_ac[is][ibl];
        }
    }

    // calculate new CL from this new Qtan
    let sa = st.alfa.sin();
    let ca = st.alfa.cos();
    let beta = (1.0 - st.minf * st.minf).sqrt();
    let beta_msq = -0.5 / beta;
    let bfac = 0.5 * (st.minf * st.minf) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    let mut clnew = 0.0;
    let mut cl_a = 0.0;
    let mut cl_ms = 0.0;
    let mut cl_ac = 0.0;

    let qinf = st.qinf;
    let mut i = 1;
    let mut cginc = 1.0 - (qnew[i] / qinf) * (qnew[i] / qinf);
    let mut cpg1 = cginc / (beta + bfac * cginc);
    let mut cpg1_ms = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
    let mut cpi_q = -2.0 * qnew[i] / (qinf * qinf);
    let mut cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    let mut cpg1_ac = cpc_cpi * cpi_q * q_ac[i];

    while i <= n {
        let ip = if i == n { 1 } else { i + 1 };
        cginc = 1.0 - (qnew[ip] / qinf) * (qnew[ip] / qinf);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_ms = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        cpi_q = -2.0 * qnew[ip] / (qinf * qinf);
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_ac = cpc_cpi * cpi_q * q_ac[ip];

        let dx = (st.x[ip] - st.x[i]) * ca + (st.y[ip] - st.y[i]) * sa;
        let dx_a = -(st.x[ip] - st.x[i]) * sa + (st.y[ip] - st.y[i]) * ca;

        let ag = 0.5 * (cpg2 + cpg1);
        let ag_ms = 0.5 * (cpg2_ms + cpg1_ms);
        let ag_ac = 0.5 * (cpg2_ac + cpg1_ac);

        clnew += dx * ag;
        cl_a += dx_a * ag;
        cl_ms += dx * ag_ms;
        cl_ac += dx * ag_ac;

        cpg1 = cpg2;
        cpg1_ms = cpg2_ms;
        cpg1_ac = cpg2_ac;
        i += 1;
    }

    // initialize under-relaxation factor
    let mut rlx = 1.0;

    let dac;
    if st.lalfa {
        // alpha is prescribed: AC is CL
        // set change in Re to account for CL changing, since Re = Re(CL)
        dac = (clnew - st.cl) / (1.0 - cl_ac - cl_ms * 2.0 * st.minf * minf_cl);
        // set under-relaxation factor if Re change is too large
        if rlx * dac > dclmax {
            rlx = dclmax / dac;
        }
        if rlx * dac < dclmin {
            rlx = dclmin / dac;
        }
    } else {
        // CL is prescribed: AC is alpha
        // set change in alpha to drive CL to prescribed value
        dac = (clnew - st.clspec) / (0.0 - cl_ac - cl_a);
        // set under-relaxation factor if alpha change is too large
        if rlx * dac > dalmax {
            rlx = dalmax / dac;
        }
        if rlx * dac < dalmin {
            rlx = dalmin / dac;
        }
    }

    let mut rmsbl = 0.0_f64;
    let mut rmxbl = 0.0_f64;
    let mut vmxbl = ' ';
    let mut imxbl = 0;
    let mut ismxbl = 0;

    let dhi = 1.5;
    let dlo = -0.5;

    // calculate changes in BL variables and under-relaxation if needed
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            let iv = st.isys[is][ibl];

            // set changes without underrelaxation
            let dctau = vdel[iv - 1][0][0] - dac * vdel[iv - 1][0][1];
            let dthet = vdel[iv - 1][1][0] - dac * vdel[iv - 1][1][1];
            let dmass = vdel[iv - 1][2][0] - dac * vdel[iv - 1][2][1];
            let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - st.uedg[is][ibl];
            let ddstr = (dmass - st.dstr[is][ibl] * duedg) / st.uedg[is][ibl];

            // normalize changes
            let dn1 = if ibl < st.itran[is] {
                dctau / 10.0
            } else {
                dctau / st.ctau[is][ibl]
            };
            let dn2 = dthet / st.thet[is][ibl];
            let dn3 = ddstr / st.dstr[is][ibl];
            let dn4 = duedg.abs() / 0.25;

            // accumulate for rms change
            rmsbl = rmsbl + dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

            // see if Ctau needs underrelaxation
            let rdn1 = rlx * dn1;
            if dn1.abs() > rmxbl.abs() {
                rmxbl = dn1;
                vmxbl = if ibl < st.itran[is] { 'n' } else { 'C' };
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn1 > dhi {
                rlx = dhi / dn1;
            }
            if rdn1 < dlo {
                rlx = dlo / dn1;
            }

            // see if Theta needs underrelaxation
            let rdn2 = rlx * dn2;
            if dn2.abs() > rmxbl.abs() {
                rmxbl = dn2;
                vmxbl = 'T';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn2 > dhi {
                rlx = dhi / dn2;
            }
            if rdn2 < dlo {
                rlx = dlo / dn2;
            }

            // see if Dstar needs underrelaxation
            let rdn3 = rlx * dn3;
            if dn3.abs() > rmxbl.abs() {
                rmxbl = dn3;
                vmxbl = 'D';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn3 > dhi {
                rlx = dhi / dn3;
            }
            if rdn3 < dlo {
                rlx = dlo / dn3;
            }

            // see if Ue needs underrelaxation
            let rdn4 = rlx * dn4;
            if dn4.abs() > rmxbl.abs() {
                rmxbl = duedg;
                vmxbl = 'U';
                imxbl = ibl;
                ismxbl = is;
            }
            if rdn4 > dhi {
                rlx = dhi / dn4;
            }
            if rdn4 < dlo {
                rlx = dlo / dn4;
            }
        }
    }

    // set true rms change
    rmsbl = (rmsbl / (4.0 * ((st.nbl[1] + st.nbl[2]) as f64))).sqrt();

    if st.lalfa {
        // set underrelaxed change in Reynolds number from change in lift
        st.cl += rlx * dac;
    } else {
        // set underrelaxed change in alpha
        st.alfa += rlx * dac;
    }

    // update BL variables with underrelaxed changes
    for is in 1..=2 {
        for ibl in 2..=st.nbl[is] {
            let iv = st.isys[is][ibl];
            let dctau = vdel[iv - 1][0][0] - dac * vdel[iv - 1][0][1];
            let dthet = vdel[iv - 1][1][0] - dac * vdel[iv - 1][1][1];
            let dmass = vdel[iv - 1][2][0] - dac * vdel[iv - 1][2][1];
            let duedg = unew[is][ibl] + dac * u_ac[is][ibl] - st.uedg[is][ibl];
            let ddstr = (dmass - st.dstr[is][ibl] * duedg) / st.uedg[is][ibl];

            st.ctau[is][ibl] += rlx * dctau;
            st.thet[is][ibl] += rlx * dthet;
            st.dstr[is][ibl] += rlx * ddstr;
            st.uedg[is][ibl] += rlx * duedg;

            let dswaki = if ibl > st.iblte[is] {
                st.wgap[ibl - st.iblte[is]]
            } else {
                0.0
            };

            // eliminate absurd transients
            if ibl >= st.itran[is] {
                st.ctau[is][ibl] = st.ctau[is][ibl].min(0.25);
            }
            let hklim = if ibl <= st.iblte[is] { 1.02 } else { 1.00005 };
            let ue = st.uedg[is][ibl];
            let msq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
            let mut dsw = st.dstr[is][ibl] - dswaki;
            limit_dstar(&mut dsw, st.thet[is][ibl], ue, msq, hklim);
            st.dstr[is][ibl] = dsw + dswaki;

            // set new mass defect (nonlinear update)
            st.mass[is][ibl] = st.dstr[is][ibl] * st.uedg[is][ibl];
        }

        // make sure there are no "islands" of negative Ue
        for ibl in 3..=st.iblte[is] {
            if st.uedg[is][ibl - 1] > 0.0 && st.uedg[is][ibl] <= 0.0 {
                st.uedg[is][ibl] = st.uedg[is][ibl - 1];
                st.mass[is][ibl] = st.dstr[is][ibl] * st.uedg[is][ibl];
            }
        }
    }

    // equate upper wake arrays to lower wake arrays
    for kbl in 1..=(st.nbl[2] - st.iblte[2]) {
        let (i1, i2) = (st.iblte[1] + kbl, st.iblte[2] + kbl);
        st.ctau[1][i1] = st.ctau[2][i2];
        st.thet[1][i1] = st.thet[2][i2];
        st.dstr[1][i1] = st.dstr[2][i2];
        st.uedg[1][i1] = st.uedg[2][i2];
        st.tau[1][i1] = st.tau[2][i2];
        st.dis[1][i1] = st.dis[2][i2];
        st.ctq[1][i1] = st.ctq[2][i2];
        st.delt[1][i1] = st.delt[2][i2];
        st.tstr[1][i1] = st.tstr[2][i2];
    }

    UpdateSummary {
        relaxation: rlx,
        residual: rmsbl,
        residual_max: rmxbl,
        residual_max_variable: vmxbl,
        i_residual_max_station: imxbl,
        residual_max_side: ismxbl,
        free_variable_change: dac,
        cl_new: clnew,
        cl_d_alpha: cl_a,
        cl_d_machsqd: cl_ms,
        cl_d_free: cl_ac,
    }
}
