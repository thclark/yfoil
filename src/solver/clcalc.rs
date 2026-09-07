//! CPCALC, CLCALC, CDCALC (xfoil.f) and COMSET on `SolverState`, line-for-line.

use crate::solver::blstate::SolverState;

/// COMSET: Kármán–Tsien parameter TKLAM and its M² derivative for the current MINF.
/// (CPSTAR/QSTAR, the sonic Cp and speed, are plotting quantities and are not kept.)
pub fn set_compressibility(st: &mut SolverState) {
    let beta = (1.0 - st.mach * st.mach).sqrt();
    let beta_msq = -0.5 / beta;
    st.karman_tsien = (st.mach * st.mach) / ((1.0 + beta) * (1.0 + beta));
    st.karman_tsien_d_machsqd = 1.0 / ((1.0 + beta) * (1.0 + beta)) - 2.0 * st.karman_tsien / (1.0 + beta) * beta_msq;
}

/// CPCALC: compressible Cp from speed, for `q[1..=n]` (1-based, slot 0 unused). Returns the
/// Cp array in the same layout. (XFOIL only warns when the Kármán–Tsien denominator goes
/// non-positive; the values are still those it computes.)
pub fn compute_cp(n: usize, q: &[f64], qinf: f64, minf: f64) -> Vec<f64> {
    let beta = (1.0 - minf * minf).sqrt();
    let bfac = 0.5 * (minf * minf) / (1.0 + beta);
    let mut cp = vec![0.0; n + 1];
    for i in 1..=n {
        let cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
        let den = beta + bfac * cpinc;
        cp[i] = cpinc / den;
    }
    cp
}

/// CLCALC: integrates surface pressures from GAM to get CL, CM and CDP, and dCL/dalpha,
/// dCL/dM² for the prescribed-CL routines. Uses the moment reference `st.xcmref/ycmref`.
pub fn compute_cl_cm(st: &mut SolverState) {
    let n = st.n_foil_nodes;
    let (x, y, gam, gam_a) = (&st.x, &st.y, &st.gamma, &st.gamma_d_alpha);
    let (xref, yref) = (st.cm_ref_x, st.cm_ref_y);
    let (minf, qinf) = (st.mach, st.qinf);

    let sa = st.alpha.sin();
    let ca = st.alpha.cos();

    let beta = (1.0 - minf * minf).sqrt();
    let beta_msq = -0.5 / beta;
    let bfac = 0.5 * (minf * minf) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    let mut cl = 0.0;
    let mut cm = 0.0;
    let mut cdp = 0.0;
    let mut cl_alf = 0.0;
    let mut cl_msq = 0.0;

    let mut i = 1;
    let mut cginc = 1.0 - (gam[i] / qinf) * (gam[i] / qinf);
    let mut cpg1 = cginc / (beta + bfac * cginc);
    let mut cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
    let mut cpi_gam = -2.0 * gam[i] / (qinf * qinf);
    let mut cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    let mut cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];

    while i <= n {
        let ip = if i == n { 1 } else { i + 1 };
        cginc = 1.0 - (gam[ip] / qinf) * (gam[ip] / qinf);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        cpi_gam = -2.0 * gam[ip] / (qinf * qinf);
        cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        let dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        let dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        let dg = cpg2 - cpg1;

        let ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        let ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        let ag = 0.5 * (cpg2 + cpg1);

        let dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
        let ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        let ag_msq = 0.5 * (cpg2_msq + cpg1_msq);

        cl += dx * ag;
        cdp -= dy * ag;
        cm = cm - dx * (ag * ax + dg * dx / 12.0) - dy * (ag * ay + dg * dy / 12.0);

        cl_alf = cl_alf + dx * ag_alf + ag * dx_alf;
        cl_msq += dx * ag_msq;

        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
        i += 1;
    }

    st.cl = cl;
    st.cm = cm;
    st.cd_pressure = cdp;
    st.cl_d_alpha = cl_alf;
    st.cl_d_machsqd = cl_msq;
}

/// CDCALC: total CD from the wake end by the Squire–Young extrapolation (with the
/// Kármán–Tsien correction) and the friction drag CDF from the surface TAU integral.
pub fn compute_cd(st: &mut SolverState) {
    let sa = st.alpha.sin();
    let ca = st.alpha.cos();

    if st.viscous && st.bl_initialised {
        // set variables at the end of the wake
        let nbl2 = st.n_stations[2];
        let thwake = st.theta[2][nbl2];
        let urat = st.ue[2][nbl2] / st.qinf;
        let uewake = st.ue[2][nbl2] * (1.0 - st.karman_tsien) / (1.0 - st.karman_tsien * (urat * urat));
        let shwake = st.dstar[2][nbl2] / st.theta[2][nbl2];

        // extrapolate wake to downstream infinity using Squire-Young relation
        // (reduces errors of the wake not being long enough)
        st.cd = 2.0 * thwake * (uewake / st.qinf).powf(0.5 * (5.0 + shwake));
    } else {
        st.cd = 0.0;
    }

    // calculate friction drag coefficient
    let mut cdf = 0.0;
    for is in 1..=2 {
        for ibl in 3..=st.i_te_station[is] {
            let i = st.i_node[is][ibl];
            let im = st.i_node[is][ibl - 1];
            let dx = (st.x[i] - st.x[im]) * ca + (st.y[i] - st.y[im]) * sa;
            cdf += 0.5 * (st.tau[is][ibl] + st.tau[is][ibl - 1]) * dx * 2.0 / (st.qinf * st.qinf);
        }
    }
    st.cd_friction = cdf;
}
