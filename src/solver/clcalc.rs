//! CPCALC, CLCALC, CDCALC (xfoil.f) and COMSET on `SolverState`, line-for-line.

use crate::solver::blstate::SolverState;

/// COMSET: Kármán–Tsien parameter TKLAM and its M² derivative for the current MINF.
/// (CPSTAR/QSTAR, the sonic Cp and speed, are plotting quantities and are not kept.)
#[doc(alias = "COMSET")]
pub fn set_compressibility(state: &mut SolverState) {
    let beta = (1.0 - state.mach * state.mach).sqrt();
    let beta_msq = -0.5 / beta;
    state.karman_tsien = (state.mach * state.mach) / ((1.0 + beta) * (1.0 + beta));
    state.karman_tsien_d_machsqd =
        1.0 / ((1.0 + beta) * (1.0 + beta)) - 2.0 * state.karman_tsien / (1.0 + beta) * beta_msq;
}

/// CPCALC: compressible Cp from speed, for `q[1..=n]` (1-based, slot 0 unused). Returns the
/// Cp array in the same layout. (XFOIL only warns when the Kármán–Tsien denominator goes
/// non-positive; the values are still those it computes.)
#[doc(alias = "CPCALC")]
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
#[doc(alias = "CLCALC")]
pub fn compute_cl_cm(state: &mut SolverState) {
    let n = state.n_foil_nodes;
    let (x, y, gam, gam_a) = (&state.x, &state.y, &state.gamma, &state.gamma_d_alpha);
    let (xref, yref) = (state.cm_ref_x, state.cm_ref_y);
    let (minf, qinf) = (state.mach, state.qinf);

    let sa = state.alpha.sin();
    let ca = state.alpha.cos();

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

    state.cl = cl;
    state.cm = cm;
    state.cd_pressure = cdp;
    state.cl_d_alpha = cl_alf;
    state.cl_d_machsqd = cl_msq;
}

/// CDCALC: total CD from the wake end by the Squire–Young extrapolation (with the
/// Kármán–Tsien correction) and the friction drag CDF from the surface TAU integral.
#[doc(alias = "CDCALC")]
pub fn compute_cd(state: &mut SolverState) {
    let sa = state.alpha.sin();
    let ca = state.alpha.cos();

    if state.viscous && state.bl_initialised {
        // set variables at the end of the wake
        let nbl2 = state.n_stations[2];
        let thwake = state.theta[2][nbl2];
        let urat = state.ue[2][nbl2] / state.qinf;
        let uewake = state.ue[2][nbl2] * (1.0 - state.karman_tsien) / (1.0 - state.karman_tsien * (urat * urat));
        let shwake = state.dstar[2][nbl2] / state.theta[2][nbl2];

        // extrapolate wake to downstream infinity using Squire-Young relation
        // (reduces errors of the wake not being long enough)
        state.cd = 2.0 * thwake * (uewake / state.qinf).powf(0.5 * (5.0 + shwake));
    } else {
        state.cd = 0.0;
    }

    // calculate friction drag coefficient
    let mut cdf = 0.0;
    for side in 1..=2 {
        for i_station in 3..=state.i_te_station[side] {
            let i = state.i_node[side][i_station];
            let im = state.i_node[side][i_station - 1];
            let dx = (state.x[i] - state.x[im]) * ca + (state.y[i] - state.y[im]) * sa;
            cdf += 0.5 * (state.tau[side][i_station] + state.tau[side][i_station - 1]) * dx * 2.0
                / (state.qinf * state.qinf);
        }
    }
    state.cd_friction = cdf;
}
