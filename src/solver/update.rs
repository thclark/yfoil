//! UPDATE (xbl.f): adds the Newton deltas to the boundary layer variables, checks for
//! excessive changes and under-relaxes if necessary, calculates the max and rms changes, and
//! the change in the global variable "AC" (CL if LALFA, alpha otherwise). Line-for-line on
//! `SolverState`.
//!
//! XFOIL aliases `UNEW/U_AC` onto `VA` and `QNEW/Q_AC` onto `VB` with EQUIVALENCE; here they
//! are plain locals — BLSOLV's input was consumed by `blsolv`, so nothing can read VA/VB after
//! the solve by construction.

use crate::bl::system::limit_dstar;
use crate::bl::system::MachClDependence;
use crate::solver::blstate::SolverState;

/// VMXBL: which primary variable had the largest normalised Newton change. Displays as
/// XFOIL's character (`n`, `C`, `T`, `D`, `U`; blank before any change).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ResidualMaxVariable {
    #[default]
    Unset,
    /// Amplification factor N (laminar)
    Ampl,
    /// Cτ^½ (turbulent)
    Sqrtctau,
    Theta,
    Dstar,
    Ue,
}

impl std::fmt::Display for ResidualMaxVariable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Unset => " ",
            Self::Ampl => "n",
            Self::Sqrtctau => "C",
            Self::Theta => "T",
            Self::Dstar => "D",
            Self::Ue => "U",
        })
    }
}

/// What UPDATE reports besides the arrays it writes into `SolverState`.
#[derive(Debug, Clone, Default)]
pub struct UpdateSummary {
    pub relaxation: f64,
    pub residual: f64,
    /// Largest normalised change (signed; for Ue XFOIL stores the raw DUEDG)
    pub residual_max: f64,
    /// 'n' (amplification), 'C' (Ctau), 'T' (theta), 'D' (delta*) or 'U' (Ue)
    pub residual_max_variable: ResidualMaxVariable,
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
#[doc(alias = "UPDATE")]
pub fn apply_newton_update(state: &mut SolverState, vdel: &[[[f64; 2]; 3]], minf_cl: f64) -> UpdateSummary {
    let pi = 4.0 * (1.0_f64).atan();
    let dtor = pi / 180.0;
    let gamm1 = state.gamma_gas - 1.0;

    // max allowable alpha changes per iteration
    let dalmax = 0.5 * dtor;
    let dalmin = -0.5 * dtor;

    // max allowable CL change per iteration
    let dclmax = 0.5;
    let mut dclmin = -0.5;
    if state.mach_cl_dependence != MachClDependence::Fixed {
        dclmin = (-0.5_f64).max(-0.9 * state.cl);
    }

    let hstinv = gamm1 * ((state.mach / state.qinf) * (state.mach / state.qinf))
        / (1.0 + 0.5 * gamm1 * (state.mach * state.mach));

    let nmax = state.ue[1].len().max(state.ue[2].len());
    let mut unew: [Vec<f64>; 3] = [Vec::new(), vec![0.0; nmax], vec![0.0; nmax]];
    let mut u_ac: [Vec<f64>; 3] = [Vec::new(), vec![0.0; nmax], vec![0.0; nmax]];

    // calculate new Ue distribution assuming no under-relaxation
    // also set the sensitivity of Ue wrt to alpha or Re
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            let i = state.i_node[side][i_station];
            let mut dui = 0.0;
            let mut dui_ac = 0.0;
            for j_side in 1..=2 {
                for j_station in 2..=state.n_stations[j_side] {
                    let j = state.i_node[j_side][j_station];
                    let j_row = state.i_row[j_side][j_station];
                    let ue_m = -state.velocity_sign[side][i_station]
                        * state.velocity_sign[j_side][j_station]
                        * state.dij[i][j];
                    dui += ue_m * (state.mass_defect[j_side][j_station] + vdel[j_row - 1][2][0]);
                    dui_ac += ue_m * (-vdel[j_row - 1][2][1]);
                }
            }
            // UINV depends on "AC" only if "AC" is alpha
            let uinv_ac = if state.alpha_specified {
                0.0
            } else {
                state.ue_inviscid_d_alpha[side][i_station]
            };
            unew[side][i_station] = state.ue_inviscid[side][i_station] + dui;
            u_ac[side][i_station] = uinv_ac + dui_ac;
        }
    }

    // set new Qtan from new Ue with appropriate sign change
    let n = state.n_foil_nodes;
    let mut qnew = vec![0.0; n + 2];
    let mut q_ac = vec![0.0; n + 2];
    for side in 1..=2 {
        for i_station in 2..=state.i_te_station[side] {
            let i = state.i_node[side][i_station];
            qnew[i] = state.velocity_sign[side][i_station] * unew[side][i_station];
            q_ac[i] = state.velocity_sign[side][i_station] * u_ac[side][i_station];
        }
    }

    // calculate new CL from this new Qtan
    let sa = state.alpha.sin();
    let ca = state.alpha.cos();
    let beta = (1.0 - state.mach * state.mach).sqrt();
    let beta_msq = -0.5 / beta;
    let bfac = 0.5 * (state.mach * state.mach) / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;

    let mut clnew = 0.0;
    let mut cl_a = 0.0;
    let mut cl_ms = 0.0;
    let mut cl_ac = 0.0;

    let qinf = state.qinf;
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

        let dx = (state.x[ip] - state.x[i]) * ca + (state.y[ip] - state.y[i]) * sa;
        let dx_a = -(state.x[ip] - state.x[i]) * sa + (state.y[ip] - state.y[i]) * ca;

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
    if state.alpha_specified {
        // alpha is prescribed: AC is CL
        // set change in Re to account for CL changing, since Re = Re(CL)
        dac = (clnew - state.cl) / (1.0 - cl_ac - cl_ms * 2.0 * state.mach * minf_cl);
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
        dac = (clnew - state.cl_specified) / (0.0 - cl_ac - cl_a);
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
    let mut vmxbl = ResidualMaxVariable::Unset;
    let mut imxbl = 0;
    let mut ismxbl = 0;

    let dhi = 1.5;
    let dlo = -0.5;

    // calculate changes in BL variables and under-relaxation if needed
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            let i_row = state.i_row[side][i_station];

            // set changes without underrelaxation
            let dctau = vdel[i_row - 1][0][0] - dac * vdel[i_row - 1][0][1];
            let dthet = vdel[i_row - 1][1][0] - dac * vdel[i_row - 1][1][1];
            let dmass = vdel[i_row - 1][2][0] - dac * vdel[i_row - 1][2][1];
            let duedg = unew[side][i_station] + dac * u_ac[side][i_station] - state.ue[side][i_station];
            let ddstr = (dmass - state.dstar[side][i_station] * duedg) / state.ue[side][i_station];

            // normalize changes
            let dn1 = if i_station < state.i_transition_station[side] {
                dctau / 10.0
            } else {
                dctau / state.sqrtctau[side][i_station]
            };
            let dn2 = dthet / state.theta[side][i_station];
            let dn3 = ddstr / state.dstar[side][i_station];
            let dn4 = duedg.abs() / 0.25;

            // accumulate for rms change
            rmsbl = rmsbl + dn1 * dn1 + dn2 * dn2 + dn3 * dn3 + dn4 * dn4;

            // see if Ctau needs underrelaxation
            let rdn1 = rlx * dn1;
            if dn1.abs() > rmxbl.abs() {
                rmxbl = dn1;
                vmxbl = if i_station < state.i_transition_station[side] {
                    ResidualMaxVariable::Ampl
                } else {
                    ResidualMaxVariable::Sqrtctau
                };
                imxbl = i_station;
                ismxbl = side;
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
                vmxbl = ResidualMaxVariable::Theta;
                imxbl = i_station;
                ismxbl = side;
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
                vmxbl = ResidualMaxVariable::Dstar;
                imxbl = i_station;
                ismxbl = side;
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
                vmxbl = ResidualMaxVariable::Ue;
                imxbl = i_station;
                ismxbl = side;
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
    rmsbl = (rmsbl / (4.0 * ((state.n_stations[1] + state.n_stations[2]) as f64))).sqrt();

    if state.alpha_specified {
        // set underrelaxed change in Reynolds number from change in lift
        state.cl += rlx * dac;
    } else {
        // set underrelaxed change in alpha
        state.alpha += rlx * dac;
    }

    // update BL variables with underrelaxed changes
    for side in 1..=2 {
        for i_station in 2..=state.n_stations[side] {
            let i_row = state.i_row[side][i_station];
            let dctau = vdel[i_row - 1][0][0] - dac * vdel[i_row - 1][0][1];
            let dthet = vdel[i_row - 1][1][0] - dac * vdel[i_row - 1][1][1];
            let dmass = vdel[i_row - 1][2][0] - dac * vdel[i_row - 1][2][1];
            let duedg = unew[side][i_station] + dac * u_ac[side][i_station] - state.ue[side][i_station];
            let ddstr = (dmass - state.dstar[side][i_station] * duedg) / state.ue[side][i_station];

            state.sqrtctau[side][i_station] += rlx * dctau;
            state.theta[side][i_station] += rlx * dthet;
            state.dstar[side][i_station] += rlx * ddstr;
            state.ue[side][i_station] += rlx * duedg;

            let dswaki = if i_station > state.i_te_station[side] {
                state.wake_gap[i_station - state.i_te_station[side]]
            } else {
                0.0
            };

            // eliminate absurd transients
            if i_station >= state.i_transition_station[side] {
                state.sqrtctau[side][i_station] = state.sqrtctau[side][i_station].min(0.25);
            }
            let hklim = if i_station <= state.i_te_station[side] {
                1.02
            } else {
                1.00005
            };
            let ue = state.ue[side][i_station];
            let msq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
            let mut dsw = state.dstar[side][i_station] - dswaki;
            limit_dstar(&mut dsw, state.theta[side][i_station], ue, msq, hklim);
            state.dstar[side][i_station] = dsw + dswaki;

            // set new mass defect (nonlinear update)
            state.mass_defect[side][i_station] = state.dstar[side][i_station] * state.ue[side][i_station];
        }

        // make sure there are no "islands" of negative Ue
        for i_station in 3..=state.i_te_station[side] {
            if state.ue[side][i_station - 1] > 0.0 && state.ue[side][i_station] <= 0.0 {
                state.ue[side][i_station] = state.ue[side][i_station - 1];
                state.mass_defect[side][i_station] = state.dstar[side][i_station] * state.ue[side][i_station];
            }
        }
    }

    // equate upper wake arrays to lower wake arrays
    for k_station in 1..=(state.n_stations[2] - state.i_te_station[2]) {
        let (i1, i2) = (state.i_te_station[1] + k_station, state.i_te_station[2] + k_station);
        state.sqrtctau[1][i1] = state.sqrtctau[2][i2];
        state.theta[1][i1] = state.theta[2][i2];
        state.dstar[1][i1] = state.dstar[2][i2];
        state.ue[1][i1] = state.ue[2][i2];
        state.tau[1][i1] = state.tau[2][i2];
        state.dissipation[1][i1] = state.dissipation[2][i2];
        state.sqrtctaueq[1][i1] = state.sqrtctaueq[2][i2];
        state.delta[1][i1] = state.delta[2][i2];
        state.thetastar[1][i1] = state.thetastar[2][i2];
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
