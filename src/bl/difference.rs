//! The local 4×5 Newton system of XBL.INC's /V_SYS/ block: BLDIF (with its upwinding) and
//! TRDIF (xblsys.f).

use super::params::*;
use super::station::{MidpointCf, StationState};
use super::transition::{interval_amplification_rate, Transition};

// ============================================================================
// Local BL Equation Coefficients
// ============================================================================

/// Translates XFOIL's `BLDIF`.
///
/// Upwinding parameter calculation
///
/// Returns (upw, upw_u1, upw_t1, upw_d1, upw_u2, upw_t2, upw_d2, upw_ms)
#[doc(alias = "BLDIF")]
fn upwinding(s1: &StationState, s2: &StationState, is_wake: bool) -> Upwinding {
    let hk1 = s1.hk;
    let hk2 = s2.hk;

    // Upwinding constant (less in wake)
    let hupwt = 1.0;
    let hdcon = if is_wake {
        hupwt / (hk2 * hk2)
    } else {
        5.0 * hupwt / (hk2 * hk2)
    };
    let hd_hk1 = 0.0;
    let hd_hk2 = -hdcon * 2.0 / hk2;

    // Local upwinding based on log(Hk-1) change
    let arg = ((hk2 - 1.0) / (hk1 - 1.0)).abs();
    let hl = arg.ln();
    let hl_hk1 = -1.0 / (hk1 - 1.0);
    let hl_hk2 = 1.0 / (hk2 - 1.0);

    // Upwinding parameter: 0.5 = trapezoidal, 1.0 = backward Euler
    let hlsq = (hl * hl).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;
    let upw_hl = ehh * hl * hdcon;
    let upw_hd = 0.5 * ehh * hlsq;

    let upw_hk1 = upw_hl * hl_hk1 + upw_hd * hd_hk1;
    let upw_hk2 = upw_hl * hl_hk2 + upw_hd * hd_hk2;

    Upwinding {
        weight: upw,
        weight_d_ue_station1: upw_hk1 * s1.hk_d_ue,
        weight_d_theta_station1: upw_hk1 * s1.hk_d_theta,
        weight_d_dstar_station1: upw_hk1 * s1.hk_d_dstar,
        weight_d_ue_station2: upw_hk2 * s2.hk_d_ue,
        weight_d_theta_station2: upw_hk2 * s2.hk_d_theta,
        weight_d_dstar_station2: upw_hk2 * s2.hk_d_dstar,
        weight_d_machsqd: upw_hk1 * s1.hk_d_machsqd + upw_hk2 * s2.hk_d_machsqd,
    }
}

/// Upwinding parameters
#[derive(Debug, Clone, Default)]
struct Upwinding {
    weight: f64,
    weight_d_ue_station1: f64,
    weight_d_theta_station1: f64,
    weight_d_dstar_station1: f64,
    weight_d_ue_station2: f64,
    weight_d_theta_station2: f64,
    weight_d_dstar_station2: f64,
    weight_d_machsqd: f64,
}

/// Mirrors XFOIL's `V_SYS`.
///
/// Local BL equation coefficients (from XFOIL's V_SYS)
///
/// These are the Jacobian entries for a single station pair (1→2).
#[derive(Debug, Clone, Default)]
#[doc(alias = "V_SYS")]
pub struct IntervalSystem {
    /// Jacobian w.r.t. previous station: VS1(4,5)
    /// Rows: 4 equations (momentum, shape, lag, auxiliary)
    /// Cols: 5 variables (Ctau, Theta, Dstar, Ue, X)
    pub jacobian_station1: [[f64; 5]; 4],

    /// Jacobian w.r.t. current station: VS2(4,5)
    pub jacobian_station2: [[f64; 5]; 4],

    /// Residual vector: VSREZ(4)
    pub residual: [f64; 4],

    /// Sensitivity to Reynolds number: VSR(4)
    pub residual_d_re: [f64; 4],

    /// Sensitivity to Mach squared: VSM(4)
    pub residual_d_machsqd: [f64; 4],

    /// Sensitivity to arc length: VSX(4)
    pub residual_d_xi: [f64; 4],
}

impl IntervalSystem {
    /// Set up the Newton system for a BL interval (BLDIF equivalent)
    ///
    /// This sets up the Jacobian and residual for the BL equations between
    /// two stations. The equations are:
    /// - Row 1: Amplification (laminar) or Shear lag (turbulent/wake)
    /// - Row 2: Momentum integral
    /// - Row 3: Shape parameter (energy)
    ///
    /// # Arguments
    /// * `s1` - Station 1 state (upstream)
    /// * `s2` - Station 2 state (downstream)
    /// * `cfm` - Midpoint skin friction
    /// * `flow_type` - Type of BL flow
    /// * `is_similarity` - True if station 2 is a similarity station (LE)
    #[doc(alias = "BLDIF")]
    pub fn assemble_interval_equations(
        &mut self,
        s1: &StationState,
        s2: &StationState,
        cfm: &MidpointCf,
        flow_type: FlowRegime,
        is_similarity: bool,
        acrit: f64,
        amplification_model: AmplificationModel,
    ) {
        // Initialize to zero
        for k in 0..4 {
            self.residual[k] = 0.0;
            self.residual_d_machsqd[k] = 0.0;
            self.residual_d_re[k] = 0.0;
            self.residual_d_xi[k] = 0.0;
            for l in 0..5 {
                self.jacobian_station1[k][l] = 0.0;
                self.jacobian_station2[k][l] = 0.0;
            }
        }

        // Logarithmic differences
        let (xlog, ulog, tlog, hlog, ddlog) = if is_similarity {
            // Similarity station: prescribed differences
            (1.0, BULE, 0.5 * (1.0 - BULE), 0.0, 0.0)
        } else {
            // Normal station: compute from values
            let xlog = (s2.xi / s1.xi).ln();
            let ulog = (s2.ue / s1.ue).ln();
            let tlog = (s2.theta / s1.theta).ln();
            let hlog = (s2.hstar / s1.hstar).ln();
            (xlog, ulog, tlog, hlog, 1.0)
        };

        // Compute upwinding parameters
        let is_wake = flow_type == FlowRegime::Wake;
        let upw = upwinding(s1, s2, is_wake);

        // Equation 1: Amplification (laminar) or Shear lag (turbulent/wake)
        match flow_type {
            FlowRegime::Laminar if is_similarity => {
                // LE point: set zero amplification factor (XFOIL: VS2(1,1) = 1.0)
                // This ensures the pivot in BLSOLV is well-conditioned
                self.jacobian_station2[0][0] = 1.0;
                self.residual[0] = -s2.ampl;
            }
            FlowRegime::Laminar => {
                // laminar part --> set amplification equation (BLDIF ITYP=1), verbatim:
                // set average amplification AX over interval X1..X2
                let r = interval_amplification_rate(
                    s1.hk,
                    s1.theta,
                    s1.retheta,
                    s1.ampl,
                    s2.hk,
                    s2.theta,
                    s2.retheta,
                    s2.ampl,
                    acrit,
                    amplification_model,
                );
                let ax = r.rate;
                let rezc = s2.ampl - s1.ampl - ax * (s2.xi - s1.xi);
                let z_ax = -(s2.xi - s1.xi);

                self.jacobian_station1[0][0] = z_ax * r.rate_d_ampl_station1 - 1.0;
                self.jacobian_station1[0][1] = z_ax
                    * (r.rate_d_hk_station1 * s1.hk_d_theta
                        + r.rate_d_theta_station1
                        + r.rate_d_retheta_station1 * s1.retheta_d_theta);
                self.jacobian_station1[0][2] = z_ax * (r.rate_d_hk_station1 * s1.hk_d_dstar);
                self.jacobian_station1[0][3] =
                    z_ax * (r.rate_d_hk_station1 * s1.hk_d_ue + r.rate_d_retheta_station1 * s1.retheta_d_ue);
                self.jacobian_station1[0][4] = ax;
                self.jacobian_station2[0][0] = z_ax * r.rate_d_ampl_station2 + 1.0;
                self.jacobian_station2[0][1] = z_ax
                    * (r.rate_d_hk_station2 * s2.hk_d_theta
                        + r.rate_d_theta_station2
                        + r.rate_d_retheta_station2 * s2.retheta_d_theta);
                self.jacobian_station2[0][2] = z_ax * (r.rate_d_hk_station2 * s2.hk_d_dstar);
                self.jacobian_station2[0][3] =
                    z_ax * (r.rate_d_hk_station2 * s2.hk_d_ue + r.rate_d_retheta_station2 * s2.retheta_d_ue);
                self.jacobian_station2[0][4] = -ax;
                self.residual_d_machsqd[0] = z_ax
                    * (r.rate_d_hk_station1 * s1.hk_d_machsqd
                        + r.rate_d_retheta_station1 * s1.retheta_d_machsqd
                        + r.rate_d_hk_station2 * s2.hk_d_machsqd
                        + r.rate_d_retheta_station2 * s2.retheta_d_machsqd);
                self.residual_d_re[0] =
                    z_ax * (r.rate_d_retheta_station1 * s1.retheta_d_re + r.rate_d_retheta_station2 * s2.retheta_d_re);
                self.residual_d_xi[0] = 0.0;
                self.residual[0] = -rezc;
            }
            FlowRegime::Turbulent | FlowRegime::Wake => {
                // Shear lag equation
                self.shear_lag_equation(s1, s2, &upw, flow_type);
            }
        }

        // Equation 2: Momentum integral equation
        self.momentum_equation(s1, s2, cfm, xlog, ulog, tlog, ddlog);

        // Equation 3: Shape parameter equation
        self.shape_equation(s1, s2, &upw, xlog, ulog, hlog, ddlog);
    }

    /// BLDIF (xblsys.f) "turbulent part --> set shear lag equation" (row 1), line for line.
    #[doc(alias = "BLDIF")]
    fn shear_lag_equation(&mut self, s1: &StationState, s2: &StationState, upw: &Upwinding, flow_type: FlowRegime) {
        let u = upw.weight;
        let sa = (1.0 - u) * s1.sqrtctau + u * s2.sqrtctau;
        let cqa = (1.0 - u) * s1.sqrtctaueq + u * s2.sqrtctaueq;
        let cfa = (1.0 - u) * s1.cf + u * s2.cf;
        let hka = (1.0 - u) * s1.hk + u * s2.hk;
        let usa = 0.5 * (s1.us + s2.us);
        let rta = 0.5 * (s1.retheta + s2.retheta);
        let dea = 0.5 * (s1.delta + s2.delta);
        let da = 0.5 * (s1.dstar + s2.dstar);
        // increased dissipation length in wake (decrease its reciprocal)
        let ald = if flow_type == FlowRegime::Wake {
            WAKE_DISSIPATION_LENGTH_RATIO
        } else {
            1.0
        };

        // set and linearize  equilibrium 1/Ue dUe/dx   ...  NEW  12 Oct 94
        let (hkc, hkc_hka, hkc_rta) = if flow_type == FlowRegime::Turbulent {
            let gcc = GBETA_LOCUS_WALL;
            let mut hkc = hka - 1.0 - gcc / rta;
            let mut hkc_hka = 1.0;
            let mut hkc_rta = gcc / (rta * rta);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hka = 0.0;
                hkc_rta = 0.0;
            }
            (hkc, hkc_hka, hkc_rta)
        } else {
            (hka - 1.0, 1.0, 0.0)
        };
        let hr = hkc / (GBETA_LOCUS_A * ald * hka);
        let hr_hka = hkc_hka / (GBETA_LOCUS_A * ald * hka) - hr / hka;
        let _hr_rta = hkc_rta / (GBETA_LOCUS_A * ald * hka);
        let uq = (0.5 * cfa - hr * hr) / (GBETA_LOCUS_B * da);
        let uq_hka = -2.0 * hr * hr_hka / (GBETA_LOCUS_B * da);
        let uq_cfa = 0.5 / (GBETA_LOCUS_B * da);
        let uq_da = -uq / da;
        // (XFOIL also forms UQ_RTA and UQ_T1..UQ_RE here; none of them enter the Jacobian below)

        let scc = LAG_CONSTANT * 1.333 / (1.0 + usa);
        let scc_usa = -scc / (1.0 + usa);
        let scc_us1 = scc_usa * 0.5;
        let scc_us2 = scc_usa * 0.5;
        let _ = (scc_us1, scc_us2);

        let slog = (s2.sqrtctau / s1.sqrtctau).ln();
        let dxi = s2.xi - s1.xi;
        let ulog = (s2.ue / s1.ue).ln();

        let rezc = scc * (cqa - sa * ald) * dxi - dea * 2.0 * slog
            + dea * 2.0 * (uq * dxi - ulog) * LAG_PRESSURE_GRADIENT_WEIGHT;

        let z_cfa = dea * 2.0 * uq_cfa * dxi * LAG_PRESSURE_GRADIENT_WEIGHT;
        let z_hka = dea * 2.0 * uq_hka * dxi * LAG_PRESSURE_GRADIENT_WEIGHT;
        let z_da = dea * 2.0 * uq_da * dxi * LAG_PRESSURE_GRADIENT_WEIGHT;
        let z_sl = -dea * 2.0;
        let z_ul = -dea * 2.0 * LAG_PRESSURE_GRADIENT_WEIGHT;
        let z_dxi = scc * (cqa - sa * ald) + dea * 2.0 * uq * LAG_PRESSURE_GRADIENT_WEIGHT;
        let z_usa = scc_usa * (cqa - sa * ald) * dxi;
        let z_cqa = scc * dxi;
        let z_sa = -scc * dxi * ald;
        let z_dea = 2.0 * ((uq * dxi - ulog) * LAG_PRESSURE_GRADIENT_WEIGHT - slog);
        let z_upw = z_cqa * (s2.sqrtctaueq - s1.sqrtctaueq)
            + z_sa * (s2.sqrtctau - s1.sqrtctau)
            + z_cfa * (s2.cf - s1.cf)
            + z_hka * (s2.hk - s1.hk);

        let z_de1 = 0.5 * z_dea;
        let z_de2 = 0.5 * z_dea;
        let z_us1 = 0.5 * z_usa;
        let z_us2 = 0.5 * z_usa;
        let z_d1 = 0.5 * z_da;
        let z_d2 = 0.5 * z_da;
        let z_u1 = -z_ul / s1.ue;
        let z_u2 = z_ul / s2.ue;
        let z_x1 = -z_dxi;
        let z_x2 = z_dxi;
        let z_s1 = (1.0 - u) * z_sa - z_sl / s1.sqrtctau;
        let z_s2 = u * z_sa + z_sl / s2.sqrtctau;
        let z_cq1 = (1.0 - u) * z_cqa;
        let z_cq2 = u * z_cqa;
        let z_cf1 = (1.0 - u) * z_cfa;
        let z_cf2 = u * z_cfa;
        let z_hk1 = (1.0 - u) * z_hka;
        let z_hk2 = u * z_hka;

        self.jacobian_station1[0][0] = z_s1;
        self.jacobian_station1[0][1] =
            z_upw * upw.weight_d_theta_station1 + z_de1 * s1.delta_d_theta + z_us1 * s1.us_d_theta;
        self.jacobian_station1[0][2] =
            z_d1 + z_upw * upw.weight_d_dstar_station1 + z_de1 * s1.delta_d_dstar + z_us1 * s1.us_d_dstar;
        self.jacobian_station1[0][3] =
            z_u1 + z_upw * upw.weight_d_ue_station1 + z_de1 * s1.delta_d_ue + z_us1 * s1.us_d_ue;
        self.jacobian_station1[0][4] = z_x1;
        self.jacobian_station2[0][0] = z_s2;
        self.jacobian_station2[0][1] =
            z_upw * upw.weight_d_theta_station2 + z_de2 * s2.delta_d_theta + z_us2 * s2.us_d_theta;
        self.jacobian_station2[0][2] =
            z_d2 + z_upw * upw.weight_d_dstar_station2 + z_de2 * s2.delta_d_dstar + z_us2 * s2.us_d_dstar;
        self.jacobian_station2[0][3] =
            z_u2 + z_upw * upw.weight_d_ue_station2 + z_de2 * s2.delta_d_ue + z_us2 * s2.us_d_ue;
        self.jacobian_station2[0][4] = z_x2;
        self.residual_d_machsqd[0] = z_upw * upw.weight_d_machsqd
            + z_de1 * s1.delta_d_machsqd
            + z_us1 * s1.us_d_machsqd
            + z_de2 * s2.delta_d_machsqd
            + z_us2 * s2.us_d_machsqd;

        self.jacobian_station1[0][1] = self.jacobian_station1[0][1]
            + z_cq1 * s1.sqrtctaueq_d_theta
            + z_cf1 * s1.cf_d_theta
            + z_hk1 * s1.hk_d_theta;
        self.jacobian_station1[0][2] = self.jacobian_station1[0][2]
            + z_cq1 * s1.sqrtctaueq_d_dstar
            + z_cf1 * s1.cf_d_dstar
            + z_hk1 * s1.hk_d_dstar;
        self.jacobian_station1[0][3] =
            self.jacobian_station1[0][3] + z_cq1 * s1.sqrtctaueq_d_ue + z_cf1 * s1.cf_d_ue + z_hk1 * s1.hk_d_ue;
        self.jacobian_station2[0][1] = self.jacobian_station2[0][1]
            + z_cq2 * s2.sqrtctaueq_d_theta
            + z_cf2 * s2.cf_d_theta
            + z_hk2 * s2.hk_d_theta;
        self.jacobian_station2[0][2] = self.jacobian_station2[0][2]
            + z_cq2 * s2.sqrtctaueq_d_dstar
            + z_cf2 * s2.cf_d_dstar
            + z_hk2 * s2.hk_d_dstar;
        self.jacobian_station2[0][3] =
            self.jacobian_station2[0][3] + z_cq2 * s2.sqrtctaueq_d_ue + z_cf2 * s2.cf_d_ue + z_hk2 * s2.hk_d_ue;
        self.residual_d_machsqd[0] = self.residual_d_machsqd[0]
            + z_cq1 * s1.sqrtctaueq_d_machsqd
            + z_cf1 * s1.cf_d_machsqd
            + z_hk1 * s1.hk_d_machsqd
            + z_cq2 * s2.sqrtctaueq_d_machsqd
            + z_cf2 * s2.cf_d_machsqd
            + z_hk2 * s2.hk_d_machsqd;
        self.residual_d_re[0] =
            z_cq1 * s1.sqrtctaueq_d_re + z_cf1 * s1.cf_d_re + z_cq2 * s2.sqrtctaueq_d_re + z_cf2 * s2.cf_d_re;
        self.residual_d_xi[0] = 0.0;
        self.residual[0] = -rezc;
    }

    /// Translates XFOIL's `BLDIF`.
    ///
    /// Set up the momentum integral equation (row 2)
    #[doc(alias = "BLDIF")]
    fn momentum_equation(
        &mut self,
        s1: &StationState,
        s2: &StationState,
        cfm: &MidpointCf,
        xlog: f64,
        ulog: f64,
        tlog: f64,
        ddlog: f64,
    ) {
        // Averaged values
        let ha = 0.5 * (s1.h + s2.h);
        let ma = 0.5 * (s1.machsqd_edge + s2.machsqd_edge);
        let xa = 0.5 * (s1.xi + s2.xi);
        let ta = 0.5 * (s1.theta + s2.theta);
        let hwa = 0.5 * (s1.wake_gap / s1.theta + s2.wake_gap / s2.theta);

        // Cf term using central CFM for accuracy
        let cfx = 0.5 * cfm.cf * xa / ta + 0.25 * (s1.cf * s1.xi / s1.theta + s2.cf * s2.xi / s2.theta);
        let cfx_xa = 0.5 * cfm.cf / ta;
        let cfx_ta = -0.5 * cfm.cf * xa / (ta * ta);
        let cfx_x1 = 0.25 * s1.cf / s1.theta + cfx_xa * 0.5;
        let cfx_x2 = 0.25 * s2.cf / s2.theta + cfx_xa * 0.5;
        let cfx_t1 = -0.25 * s1.cf * s1.xi / (s1.theta * s1.theta) + cfx_ta * 0.5;
        let cfx_t2 = -0.25 * s2.cf * s2.xi / (s2.theta * s2.theta) + cfx_ta * 0.5;
        let cfx_cf1 = 0.25 * s1.xi / s1.theta;
        let cfx_cf2 = 0.25 * s2.xi / s2.theta;
        let cfx_cfm = 0.5 * xa / ta;

        let btmp = ha + 2.0 - ma + hwa;

        // Momentum equation residual
        let rezt = tlog + btmp * ulog - xlog * 0.5 * cfx;

        // Z coefficients
        let z_cfx = -xlog * 0.5;
        let z_ha = ulog;
        let z_hwa = ulog;
        let z_ma = -ulog;
        let z_xl = -ddlog * 0.5 * cfx;
        let z_ul = ddlog * btmp;
        let z_tl = ddlog;

        let z_cfm = z_cfx * cfx_cfm;
        let z_cf1 = z_cfx * cfx_cf1;
        let z_cf2 = z_cfx * cfx_cf2;

        let z_t1 = -z_tl / s1.theta + z_cfx * cfx_t1 + z_hwa * 0.5 * (-s1.wake_gap / (s1.theta * s1.theta));
        let z_t2 = z_tl / s2.theta + z_cfx * cfx_t2 + z_hwa * 0.5 * (-s2.wake_gap / (s2.theta * s2.theta));
        let z_x1 = -z_xl / s1.xi + z_cfx * cfx_x1;
        let z_x2 = z_xl / s2.xi + z_cfx * cfx_x2;
        let z_u1 = -z_ul / s1.ue;
        let z_u2 = z_ul / s2.ue;

        // Jacobian entries for row 2
        self.jacobian_station1[1][1] =
            0.5 * z_ha * s1.h_d_theta + z_cfm * cfm.cf_d_theta_station1 + z_cf1 * s1.cf_d_theta + z_t1;
        self.jacobian_station1[1][2] =
            0.5 * z_ha * s1.h_d_dstar + z_cfm * cfm.cf_d_dstar_station1 + z_cf1 * s1.cf_d_dstar;
        self.jacobian_station1[1][3] =
            0.5 * z_ma * s1.machsqd_edge_d_ue + z_cfm * cfm.cf_d_ue_station1 + z_cf1 * s1.cf_d_ue + z_u1;
        self.jacobian_station1[1][4] = z_x1;

        self.jacobian_station2[1][1] =
            0.5 * z_ha * s2.h_d_theta + z_cfm * cfm.cf_d_theta_station2 + z_cf2 * s2.cf_d_theta + z_t2;
        self.jacobian_station2[1][2] =
            0.5 * z_ha * s2.h_d_dstar + z_cfm * cfm.cf_d_dstar_station2 + z_cf2 * s2.cf_d_dstar;
        self.jacobian_station2[1][3] =
            0.5 * z_ma * s2.machsqd_edge_d_ue + z_cfm * cfm.cf_d_ue_station2 + z_cf2 * s2.cf_d_ue + z_u2;
        self.jacobian_station2[1][4] = z_x2;

        self.residual_d_machsqd[1] = 0.5 * z_ma * s1.machsqd_edge_d_machsqd
            + z_cfm * cfm.cf_d_machsqd
            + z_cf1 * s1.cf_d_machsqd
            + 0.5 * z_ma * s2.machsqd_edge_d_machsqd
            + z_cf2 * s2.cf_d_machsqd;
        self.residual_d_re[1] = z_cfm * cfm.cf_d_re + z_cf1 * s1.cf_d_re + z_cf2 * s2.cf_d_re;

        self.residual[1] = -rezt;
    }

    /// Translates XFOIL's `BLDIF`.
    ///
    /// Set up the shape parameter equation (row 3)
    #[doc(alias = "BLDIF")]
    fn shape_equation(
        &mut self,
        s1: &StationState,
        s2: &StationState,
        upw: &Upwinding,
        xlog: f64,
        ulog: f64,
        hlog: f64,
        ddlog: f64,
    ) {
        let u = upw.weight;
        let xot1 = s1.xi / s1.theta;
        let xot2 = s2.xi / s2.theta;

        // Averaged values
        let ha = 0.5 * (s1.h + s2.h);
        let hsa = 0.5 * (s1.hstar + s2.hstar);
        let hca = 0.5 * (s1.hstarstar + s2.hstarstar);
        let hwa = 0.5 * (s1.wake_gap / s1.theta + s2.wake_gap / s2.theta);

        // Upwind-weighted DI and CF
        let dix = (1.0 - u) * s1.cdiss * xot1 + u * s2.cdiss * xot2;
        let cfx = (1.0 - u) * s1.cf * xot1 + u * s2.cf * xot2;
        let dix_upw = s2.cdiss * xot2 - s1.cdiss * xot1;
        let cfx_upw = s2.cf * xot2 - s1.cf * xot1;

        let btmp = 2.0 * hca / hsa + 1.0 - ha - hwa;

        // Shape equation residual
        let rezh = hlog + btmp * ulog + xlog * (0.5 * cfx - dix);

        // Z coefficients
        let z_cfx = xlog * 0.5;
        let z_dix = -xlog;
        let z_hca = 2.0 * ulog / hsa;
        let z_ha = -ulog;
        let z_hwa = -ulog;
        let z_xl = ddlog * (0.5 * cfx - dix);
        let z_ul = ddlog * btmp;
        let z_hl = ddlog;

        let z_upw = z_cfx * cfx_upw + z_dix * dix_upw;

        let z_hs1 = -hca * ulog / (hsa * hsa) - z_hl / s1.hstar;
        let z_hs2 = -hca * ulog / (hsa * hsa) + z_hl / s2.hstar;

        let z_cf1 = (1.0 - u) * z_cfx * xot1;
        let z_cf2 = u * z_cfx * xot2;
        let z_di1 = (1.0 - u) * z_dix * xot1;
        let z_di2 = u * z_dix * xot2;

        let z_t1 = (1.0 - u) * (z_cfx * s1.cf + z_dix * s1.cdiss) * (-xot1 / s1.theta)
            + z_hwa * 0.5 * (-s1.wake_gap / (s1.theta * s1.theta));
        let z_t2 = u * (z_cfx * s2.cf + z_dix * s2.cdiss) * (-xot2 / s2.theta)
            + z_hwa * 0.5 * (-s2.wake_gap / (s2.theta * s2.theta));
        let z_x1 = (1.0 - u) * (z_cfx * s1.cf + z_dix * s1.cdiss) / s1.theta - z_xl / s1.xi;
        let z_x2 = u * (z_cfx * s2.cf + z_dix * s2.cdiss) / s2.theta + z_xl / s2.xi;
        let z_u1 = -z_ul / s1.ue;
        let z_u2 = z_ul / s2.ue;

        // Jacobian entries for row 3
        self.jacobian_station1[2][0] = z_di1 * s1.cdiss_d_sqrtctau;
        self.jacobian_station1[2][1] = z_hs1 * s1.hstar_d_theta
            + z_cf1 * s1.cf_d_theta
            + z_di1 * s1.cdiss_d_theta
            + z_t1
            + 0.5 * (z_hca * s1.hstarstar_d_theta + z_ha * s1.h_d_theta)
            + z_upw * upw.weight_d_theta_station1;
        self.jacobian_station1[2][2] = z_hs1 * s1.hstar_d_dstar
            + z_cf1 * s1.cf_d_dstar
            + z_di1 * s1.cdiss_d_dstar
            + 0.5 * (z_hca * s1.hstarstar_d_dstar + z_ha * s1.h_d_dstar)
            + z_upw * upw.weight_d_dstar_station1;
        self.jacobian_station1[2][3] = z_hs1 * s1.hstar_d_ue
            + z_cf1 * s1.cf_d_ue
            + z_di1 * s1.cdiss_d_ue
            + z_u1
            + 0.5 * z_hca * s1.hstarstar_d_ue
            + z_upw * upw.weight_d_ue_station1;
        self.jacobian_station1[2][4] = z_x1;

        self.jacobian_station2[2][0] = z_di2 * s2.cdiss_d_sqrtctau;
        self.jacobian_station2[2][1] = z_hs2 * s2.hstar_d_theta
            + z_cf2 * s2.cf_d_theta
            + z_di2 * s2.cdiss_d_theta
            + z_t2
            + 0.5 * (z_hca * s2.hstarstar_d_theta + z_ha * s2.h_d_theta)
            + z_upw * upw.weight_d_theta_station2;
        self.jacobian_station2[2][2] = z_hs2 * s2.hstar_d_dstar
            + z_cf2 * s2.cf_d_dstar
            + z_di2 * s2.cdiss_d_dstar
            + 0.5 * (z_hca * s2.hstarstar_d_dstar + z_ha * s2.h_d_dstar)
            + z_upw * upw.weight_d_dstar_station2;
        self.jacobian_station2[2][3] = z_hs2 * s2.hstar_d_ue
            + z_cf2 * s2.cf_d_ue
            + z_di2 * s2.cdiss_d_ue
            + z_u2
            + 0.5 * z_hca * s2.hstarstar_d_ue
            + z_upw * upw.weight_d_ue_station2;
        self.jacobian_station2[2][4] = z_x2;

        self.residual_d_machsqd[2] = z_hs1 * s1.hstar_d_machsqd
            + z_cf1 * s1.cf_d_machsqd
            + z_di1 * s1.cdiss_d_machsqd
            + z_hs2 * s2.hstar_d_machsqd
            + z_cf2 * s2.cf_d_machsqd
            + z_di2 * s2.cdiss_d_machsqd
            + 0.5 * (z_hca * s1.hstarstar_d_machsqd + z_hca * s2.hstarstar_d_machsqd)
            + z_upw * upw.weight_d_machsqd;
        self.residual_d_re[2] = z_hs1 * s1.hstar_d_re
            + z_cf1 * s1.cf_d_re
            + z_di1 * s1.cdiss_d_re
            + z_hs2 * s2.hstar_d_re
            + z_cf2 * s2.cf_d_re
            + z_di2 * s2.cdiss_d_re;

        self.residual[2] = -rezh;
    }

    /// Translates XFOIL's `TRDIF`.
    ///
    /// Set up Newton system for transition interval (TRDIF equivalent)
    ///
    /// Handles intervals that span laminar-turbulent transition by:
    /// 1. Setting up laminar equations from X1 to XT (transition)
    /// 2. Setting up turbulent equations from XT to X2
    /// 3. Summing the contributions
    ///
    /// # Arguments
    /// * `s1` - Station 1 state (upstream, laminar)
    /// * `s2` - Station 2 state (downstream, turbulent)
    /// * `trans` - Transition location and derivatives
    /// * `acrit` - Critical amplification factor
    /// * `params` - Global BL parameters
    #[allow(clippy::too_many_lines)]
    #[doc(alias = "TRDIF")]
    pub fn assemble_transition_equations(
        &mut self,
        s1: &StationState,
        s2: &StationState,
        trans: &Transition,
        acrit: f64,
        params: &FlowParameters,
    ) {
        // Weighting factors for linear interpolation to transition point
        let wf2 = (trans.xi_transition - s1.xi) / (s2.xi - s1.xi);
        let wf2_xt = 1.0 / (s2.xi - s1.xi);
        let wf1 = 1.0 - wf2;

        // Derivatives of weighting factors w.r.t. station variables
        let wf2_a1 = wf2_xt * trans.xi_transition_d_ampl_station1;
        let wf2_x1 = wf2_xt * trans.xi_transition_d_xi_station1 + (wf2 - 1.0) / (s2.xi - s1.xi);
        let wf2_x2 = wf2_xt * trans.xi_transition_d_xi_station2 - wf2 / (s2.xi - s1.xi);
        let wf2_t1 = wf2_xt * trans.xi_transition_d_theta_station1;
        let wf2_t2 = wf2_xt * trans.xi_transition_d_theta_station2;
        let wf2_d1 = wf2_xt * trans.xi_transition_d_dstar_station1;
        let wf2_d2 = wf2_xt * trans.xi_transition_d_dstar_station2;
        let wf2_u1 = wf2_xt * trans.xi_transition_d_ue_station1;
        let wf2_u2 = wf2_xt * trans.xi_transition_d_ue_station2;
        let wf2_ms = wf2_xt * trans.xi_transition_d_machsqd;
        let wf2_re = wf2_xt * trans.xi_transition_d_re;
        let wf2_xf = wf2_xt * trans.xi_transition_d_x_trip;

        let wf1_a1 = -wf2_a1;
        let wf1_x1 = -wf2_x1;
        let wf1_x2 = -wf2_x2;
        let wf1_t1 = -wf2_t1;
        let wf1_t2 = -wf2_t2;
        let wf1_d1 = -wf2_d1;
        let wf1_d2 = -wf2_d2;
        let wf1_u1 = -wf2_u1;
        let wf1_u2 = -wf2_u2;
        let wf1_ms = -wf2_ms;
        let wf1_re = -wf2_re;
        let wf1_xf = -wf2_xf;

        // *** PART 1: Laminar from X1 to XT ***

        // Interpolate primary variables to transition point
        let tt = s1.theta * wf1 + s2.theta * wf2;
        let tt_a1 = s1.theta * wf1_a1 + s2.theta * wf2_a1;
        let tt_x1 = s1.theta * wf1_x1 + s2.theta * wf2_x1;
        let tt_x2 = s1.theta * wf1_x2 + s2.theta * wf2_x2;
        let tt_t1 = s1.theta * wf1_t1 + s2.theta * wf2_t1 + wf1;
        let tt_t2 = s1.theta * wf1_t2 + s2.theta * wf2_t2 + wf2;
        let tt_d1 = s1.theta * wf1_d1 + s2.theta * wf2_d1;
        let tt_d2 = s1.theta * wf1_d2 + s2.theta * wf2_d2;
        let tt_u1 = s1.theta * wf1_u1 + s2.theta * wf2_u1;
        let tt_u2 = s1.theta * wf1_u2 + s2.theta * wf2_u2;
        let tt_ms = s1.theta * wf1_ms + s2.theta * wf2_ms;
        let tt_re = s1.theta * wf1_re + s2.theta * wf2_re;
        let tt_xf = s1.theta * wf1_xf + s2.theta * wf2_xf;

        let dt = s1.dstar * wf1 + s2.dstar * wf2;
        let dt_a1 = s1.dstar * wf1_a1 + s2.dstar * wf2_a1;
        let dt_x1 = s1.dstar * wf1_x1 + s2.dstar * wf2_x1;
        let dt_x2 = s1.dstar * wf1_x2 + s2.dstar * wf2_x2;
        let dt_t1 = s1.dstar * wf1_t1 + s2.dstar * wf2_t1;
        let dt_t2 = s1.dstar * wf1_t2 + s2.dstar * wf2_t2;
        let dt_d1 = s1.dstar * wf1_d1 + s2.dstar * wf2_d1 + wf1;
        let dt_d2 = s1.dstar * wf1_d2 + s2.dstar * wf2_d2 + wf2;
        let dt_u1 = s1.dstar * wf1_u1 + s2.dstar * wf2_u1;
        let dt_u2 = s1.dstar * wf1_u2 + s2.dstar * wf2_u2;
        let dt_ms = s1.dstar * wf1_ms + s2.dstar * wf2_ms;
        let dt_re = s1.dstar * wf1_re + s2.dstar * wf2_re;
        let dt_xf = s1.dstar * wf1_xf + s2.dstar * wf2_xf;

        let ut = s1.ue * wf1 + s2.ue * wf2;
        let ut_a1 = s1.ue * wf1_a1 + s2.ue * wf2_a1;
        let ut_x1 = s1.ue * wf1_x1 + s2.ue * wf2_x1;
        let ut_x2 = s1.ue * wf1_x2 + s2.ue * wf2_x2;
        let ut_t1 = s1.ue * wf1_t1 + s2.ue * wf2_t1;
        let ut_t2 = s1.ue * wf1_t2 + s2.ue * wf2_t2;
        let ut_d1 = s1.ue * wf1_d1 + s2.ue * wf2_d1;
        let ut_d2 = s1.ue * wf1_d2 + s2.ue * wf2_d2;
        let ut_u1 = s1.ue * wf1_u1 + s2.ue * wf2_u1 + wf1;
        let ut_u2 = s1.ue * wf1_u2 + s2.ue * wf2_u2 + wf2;
        let ut_ms = s1.ue * wf1_ms + s2.ue * wf2_ms;
        let ut_re = s1.ue * wf1_re + s2.ue * wf2_re;
        let ut_xf = s1.ue * wf1_xf + s2.ue * wf2_xf;

        // Create transition-point state for laminar part
        // set primary "T" variables at XT (really placed into "2" variables): XFOIL overwrites
        // X2/T2/D2/U2/AMPL2/S2 on the saved station-2 COMMON, so U2_UEI, U2_MS and DW2 are
        // those of station 2 — no BLPRV here.
        let mut state = s2.clone();
        state.xi = trans.xi_transition;
        state.theta = tt;
        state.dstar = dt;
        state.ue = ut;
        state.ampl = acrit;
        state.sqrtctau = 0.0;
        state.set_kinematic_variables(params);
        state.set_closure_variables(FlowRegime::Laminar, params);

        // Calculate midpoint Cf for X1-XT
        let cfm_lam = MidpointCf::compute(s1, &state, FlowRegime::Laminar, false);

        // Call BLDIF for laminar part (X1 to XT)
        let mut lam_sys = IntervalSystem::default();
        lam_sys.assemble_interval_equations(
            s1,
            &state,
            &cfm_lam,
            FlowRegime::Laminar,
            false,
            acrit,
            params.amplification_model,
        );

        // Convert laminar system sensitivities from "T" variables to "1" and "2" variables
        // Using chain rule for derivatives
        let mut bl1: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut bl2: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut blrez: [f64; 4] = [0.0; 4];
        let mut blm: [f64; 4] = [0.0; 4];
        let mut blr: [f64; 4] = [0.0; 4];
        let mut blx: [f64; 4] = [0.0; 4];

        for k in 1..3 {
            // Row 2 and 3 (momentum and shape)
            blrez[k] = lam_sys.residual[k];
            blm[k] = lam_sys.residual_d_machsqd[k]
                + lam_sys.jacobian_station2[k][1] * tt_ms
                + lam_sys.jacobian_station2[k][2] * dt_ms
                + lam_sys.jacobian_station2[k][3] * ut_ms
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_machsqd;
            blr[k] = lam_sys.residual_d_re[k]
                + lam_sys.jacobian_station2[k][1] * tt_re
                + lam_sys.jacobian_station2[k][2] * dt_re
                + lam_sys.jacobian_station2[k][3] * ut_re
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_re;
            blx[k] = lam_sys.residual_d_xi[k]
                + lam_sys.jacobian_station2[k][1] * tt_xf
                + lam_sys.jacobian_station2[k][2] * dt_xf
                + lam_sys.jacobian_station2[k][3] * ut_xf
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_x_trip;

            bl1[k][0] = lam_sys.jacobian_station1[k][0]
                + lam_sys.jacobian_station2[k][1] * tt_a1
                + lam_sys.jacobian_station2[k][2] * dt_a1
                + lam_sys.jacobian_station2[k][3] * ut_a1
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_ampl_station1;
            bl1[k][1] = lam_sys.jacobian_station1[k][1]
                + lam_sys.jacobian_station2[k][1] * tt_t1
                + lam_sys.jacobian_station2[k][2] * dt_t1
                + lam_sys.jacobian_station2[k][3] * ut_t1
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_theta_station1;
            bl1[k][2] = lam_sys.jacobian_station1[k][2]
                + lam_sys.jacobian_station2[k][1] * tt_d1
                + lam_sys.jacobian_station2[k][2] * dt_d1
                + lam_sys.jacobian_station2[k][3] * ut_d1
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_dstar_station1;
            bl1[k][3] = lam_sys.jacobian_station1[k][3]
                + lam_sys.jacobian_station2[k][1] * tt_u1
                + lam_sys.jacobian_station2[k][2] * dt_u1
                + lam_sys.jacobian_station2[k][3] * ut_u1
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_ue_station1;
            bl1[k][4] = lam_sys.jacobian_station1[k][4]
                + lam_sys.jacobian_station2[k][1] * tt_x1
                + lam_sys.jacobian_station2[k][2] * dt_x1
                + lam_sys.jacobian_station2[k][3] * ut_x1
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_xi_station1;

            bl2[k][0] = 0.0; // No dA2 dependence (A2 side turbulent Ctau)
            bl2[k][1] = lam_sys.jacobian_station2[k][1] * tt_t2
                + lam_sys.jacobian_station2[k][2] * dt_t2
                + lam_sys.jacobian_station2[k][3] * ut_t2
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_theta_station2;
            bl2[k][2] = lam_sys.jacobian_station2[k][1] * tt_d2
                + lam_sys.jacobian_station2[k][2] * dt_d2
                + lam_sys.jacobian_station2[k][3] * ut_d2
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_dstar_station2;
            bl2[k][3] = lam_sys.jacobian_station2[k][1] * tt_u2
                + lam_sys.jacobian_station2[k][2] * dt_u2
                + lam_sys.jacobian_station2[k][3] * ut_u2
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_ue_station2;
            bl2[k][4] = lam_sys.jacobian_station2[k][1] * tt_x2
                + lam_sys.jacobian_station2[k][2] * dt_x2
                + lam_sys.jacobian_station2[k][3] * ut_x2
                + lam_sys.jacobian_station2[k][4] * trans.xi_transition_d_xi_station2;
        }

        // *** PART 2: Turbulent from XT to X2 ***

        // Calculate equilibrium shear coefficient CQT at transition
        state.set_closure_variables(FlowRegime::Turbulent, params);

        // Set initial shear stress: ST = CTR * CQ
        // where CTR = CTRCON * exp(-CTRCEX/(HK-1))
        let hk_minus_one = state.hk - 1.0;
        let ctr = TRANSITION_SQRTCTAU_FACTOR * (-TRANSITION_SQRTCTAU_EXPONENT / hk_minus_one).exp();
        let ctr_hk = ctr * TRANSITION_SQRTCTAU_EXPONENT / (hk_minus_one * hk_minus_one);

        let s_t = ctr * state.sqrtctaueq;
        let st_tt = ctr * state.sqrtctaueq_d_theta + state.sqrtctaueq * ctr_hk * state.hk_d_theta;
        let st_dt = ctr * state.sqrtctaueq_d_dstar + state.sqrtctaueq * ctr_hk * state.hk_d_dstar;
        let st_ut = ctr * state.sqrtctaueq_d_ue + state.sqrtctaueq * ctr_hk * state.hk_d_ue;
        let st_ms = ctr * state.sqrtctaueq_d_machsqd + state.sqrtctaueq * ctr_hk * state.hk_d_machsqd;
        let st_re = ctr * state.sqrtctaueq_d_re;

        // ST sensitivities w.r.t. actual "1" and "2" variables
        let st_a1 = st_tt * tt_a1 + st_dt * dt_a1 + st_ut * ut_a1;
        let st_x1 = st_tt * tt_x1 + st_dt * dt_x1 + st_ut * ut_x1;
        let st_x2 = st_tt * tt_x2 + st_dt * dt_x2 + st_ut * ut_x2;
        let st_t1 = st_tt * tt_t1 + st_dt * dt_t1 + st_ut * ut_t1;
        let st_t2 = st_tt * tt_t2 + st_dt * dt_t2 + st_ut * ut_t2;
        let st_d1 = st_tt * tt_d1 + st_dt * dt_d1 + st_ut * ut_d1;
        let st_d2 = st_tt * tt_d2 + st_dt * dt_d2 + st_ut * ut_d2;
        let st_u1 = st_tt * tt_u1 + st_dt * dt_u1 + st_ut * ut_u1;
        let st_u2 = st_tt * tt_u2 + st_dt * dt_u2 + st_ut * ut_u2;
        let st_ms_total = st_tt * tt_ms + st_dt * dt_ms + st_ut * ut_ms + st_ms;
        let st_re_total = st_tt * tt_re + st_dt * dt_re + st_ut * ut_re + st_re;
        let st_xf = st_tt * tt_xf + st_dt * dt_xf + st_ut * ut_xf;

        // Update transition station with turbulent initial condition
        state.sqrtctau = s_t;

        // Recalculate turbulent secondary variables with proper CTI
        state.set_closure_variables(FlowRegime::Turbulent, params);

        // Calculate midpoint Cf for XT-X2
        let cfm_turb = MidpointCf::compute(&state, s2, FlowRegime::Turbulent, false);

        // Call BLDIF for turbulent part (XT to X2)
        let mut turb_sys = IntervalSystem::default();
        turb_sys.assemble_interval_equations(
            &state,
            s2,
            &cfm_turb,
            FlowRegime::Turbulent,
            false,
            acrit,
            params.amplification_model,
        );

        // Convert turbulent system sensitivities from "T" variables to "1" and "2" variables
        let mut bt1: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut bt2: [[f64; 5]; 4] = [[0.0; 5]; 4];
        let mut btrez: [f64; 4] = [0.0; 4];
        let mut btm: [f64; 4] = [0.0; 4];
        let mut btr: [f64; 4] = [0.0; 4];
        let mut btx: [f64; 4] = [0.0; 4];

        for k in 0..3 {
            btrez[k] = turb_sys.residual[k];
            btm[k] = turb_sys.residual_d_machsqd[k]
                + turb_sys.jacobian_station1[k][0] * st_ms_total
                + turb_sys.jacobian_station1[k][1] * tt_ms
                + turb_sys.jacobian_station1[k][2] * dt_ms
                + turb_sys.jacobian_station1[k][3] * ut_ms
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_machsqd;
            btr[k] = turb_sys.residual_d_re[k]
                + turb_sys.jacobian_station1[k][0] * st_re_total
                + turb_sys.jacobian_station1[k][1] * tt_re
                + turb_sys.jacobian_station1[k][2] * dt_re
                + turb_sys.jacobian_station1[k][3] * ut_re
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_re;
            btx[k] = turb_sys.residual_d_xi[k]
                + turb_sys.jacobian_station1[k][0] * st_xf
                + turb_sys.jacobian_station1[k][1] * tt_xf
                + turb_sys.jacobian_station1[k][2] * dt_xf
                + turb_sys.jacobian_station1[k][3] * ut_xf
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_x_trip;

            bt1[k][0] = turb_sys.jacobian_station1[k][0] * st_a1
                + turb_sys.jacobian_station1[k][1] * tt_a1
                + turb_sys.jacobian_station1[k][2] * dt_a1
                + turb_sys.jacobian_station1[k][3] * ut_a1
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_ampl_station1;
            bt1[k][1] = turb_sys.jacobian_station1[k][0] * st_t1
                + turb_sys.jacobian_station1[k][1] * tt_t1
                + turb_sys.jacobian_station1[k][2] * dt_t1
                + turb_sys.jacobian_station1[k][3] * ut_t1
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_theta_station1;
            bt1[k][2] = turb_sys.jacobian_station1[k][0] * st_d1
                + turb_sys.jacobian_station1[k][1] * tt_d1
                + turb_sys.jacobian_station1[k][2] * dt_d1
                + turb_sys.jacobian_station1[k][3] * ut_d1
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_dstar_station1;
            bt1[k][3] = turb_sys.jacobian_station1[k][0] * st_u1
                + turb_sys.jacobian_station1[k][1] * tt_u1
                + turb_sys.jacobian_station1[k][2] * dt_u1
                + turb_sys.jacobian_station1[k][3] * ut_u1
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_ue_station1;
            bt1[k][4] = turb_sys.jacobian_station1[k][0] * st_x1
                + turb_sys.jacobian_station1[k][1] * tt_x1
                + turb_sys.jacobian_station1[k][2] * dt_x1
                + turb_sys.jacobian_station1[k][3] * ut_x1
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_xi_station1;

            bt2[k][0] = turb_sys.jacobian_station2[k][0];
            bt2[k][1] = turb_sys.jacobian_station2[k][1]
                + turb_sys.jacobian_station1[k][0] * st_t2
                + turb_sys.jacobian_station1[k][1] * tt_t2
                + turb_sys.jacobian_station1[k][2] * dt_t2
                + turb_sys.jacobian_station1[k][3] * ut_t2
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_theta_station2;
            bt2[k][2] = turb_sys.jacobian_station2[k][2]
                + turb_sys.jacobian_station1[k][0] * st_d2
                + turb_sys.jacobian_station1[k][1] * tt_d2
                + turb_sys.jacobian_station1[k][2] * dt_d2
                + turb_sys.jacobian_station1[k][3] * ut_d2
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_dstar_station2;
            bt2[k][3] = turb_sys.jacobian_station2[k][3]
                + turb_sys.jacobian_station1[k][0] * st_u2
                + turb_sys.jacobian_station1[k][1] * tt_u2
                + turb_sys.jacobian_station1[k][2] * dt_u2
                + turb_sys.jacobian_station1[k][3] * ut_u2
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_ue_station2;
            bt2[k][4] = turb_sys.jacobian_station2[k][4]
                + turb_sys.jacobian_station1[k][0] * st_x2
                + turb_sys.jacobian_station1[k][1] * tt_x2
                + turb_sys.jacobian_station1[k][2] * dt_x2
                + turb_sys.jacobian_station1[k][3] * ut_x2
                + turb_sys.jacobian_station1[k][4] * trans.xi_transition_d_xi_station2;
        }

        // *** COMBINE: Add laminar and turbulent parts ***

        // Row 1: Shear stress (from turbulent part only - laminar row is amplification)
        self.residual[0] = btrez[0];
        self.residual_d_machsqd[0] = btm[0];
        self.residual_d_re[0] = btr[0];
        self.residual_d_xi[0] = btx[0];
        for l in 0..5 {
            self.jacobian_station1[0][l] = bt1[0][l];
            self.jacobian_station2[0][l] = bt2[0][l];
        }

        // Rows 2 and 3: Sum laminar and turbulent contributions
        for k in 1..3 {
            self.residual[k] = blrez[k] + btrez[k];
            self.residual_d_machsqd[k] = blm[k] + btm[k];
            self.residual_d_re[k] = blr[k] + btr[k];
            self.residual_d_xi[k] = blx[k] + btx[k];
            for l in 0..5 {
                self.jacobian_station1[k][l] = bl1[k][l] + bt1[k][l];
                self.jacobian_station2[k][l] = bl2[k][l] + bt2[k][l];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    // ========================================================================
    // BLDIF Tests - validate momentum equation against XFOIL Fortran
    // ========================================================================

    #[test]
    fn test_bldif_momentum_equation() {
        // Test case from Fortran test_bldif.f90
        // Typical turbulent BL stations

        // Station 1
        let mut s1 = StationState::default();
        s1.xi = 0.10;
        s1.ue = 1.15;
        s1.theta = 0.0018;
        s1.dstar = 0.0045;
        s1.wake_gap = 0.0;
        s1.h = s1.dstar / s1.theta;
        s1.machsqd_edge = 0.0;
        s1.hstar = 1.75;
        s1.cf = 0.0025;
        s1.hk = 1.40;
        s1.retheta = 2070.0;
        s1.delta = 0.012;
        s1.sqrtctau = 0.015;
        s1.sqrtctaueq = 0.04;
        s1.us = 0.5;
        s1.cdiss = 0.0008;

        // Station 1 derivatives
        s1.h_d_theta = -s1.h / s1.theta;
        s1.h_d_dstar = 1.0 / s1.theta;
        s1.machsqd_edge_d_ue = 0.0;
        s1.machsqd_edge_d_machsqd = 0.0;
        s1.cf_d_theta = 0.0;
        s1.cf_d_dstar = 0.0;
        s1.cf_d_ue = 0.0;
        s1.cf_d_machsqd = 0.0;
        s1.cf_d_re = -s1.cf / 1e6;
        s1.hk_d_theta = 0.0;
        s1.hk_d_dstar = 0.0;
        s1.hk_d_ue = 0.0;
        s1.hk_d_machsqd = 0.0;
        s1.hstar_d_theta = 0.0;
        s1.hstar_d_dstar = 0.0;
        s1.hstar_d_ue = 0.0;
        s1.hstar_d_machsqd = 0.0;
        s1.hstar_d_re = 0.0;
        s1.delta_d_theta = 0.0;
        s1.delta_d_dstar = 0.0;
        s1.delta_d_ue = 0.0;
        s1.delta_d_machsqd = 0.0;
        s1.us_d_theta = 0.0;
        s1.us_d_dstar = 0.0;
        s1.us_d_ue = 0.0;
        s1.us_d_machsqd = 0.0;
        s1.us_d_re = 0.0;
        s1.sqrtctaueq_d_theta = 0.0;
        s1.sqrtctaueq_d_dstar = 0.0;
        s1.sqrtctaueq_d_ue = 0.0;
        s1.sqrtctaueq_d_machsqd = 0.0;
        s1.sqrtctaueq_d_re = 0.0;
        s1.cdiss_d_theta = 0.0;
        s1.cdiss_d_dstar = 0.0;
        s1.cdiss_d_ue = 0.0;
        s1.cdiss_d_sqrtctau = 0.0;
        s1.cdiss_d_machsqd = 0.0;
        s1.cdiss_d_re = 0.0;
        s1.hstarstar = 0.0;
        s1.hstarstar_d_theta = 0.0;
        s1.hstarstar_d_dstar = 0.0;
        s1.hstarstar_d_ue = 0.0;
        s1.hstarstar_d_machsqd = 0.0;
        s1.retheta_d_theta = 0.0;
        s1.retheta_d_ue = 0.0;
        s1.retheta_d_machsqd = 0.0;
        s1.retheta_d_re = 0.0;

        // Station 2
        let mut s2 = StationState::default();
        s2.xi = 0.12;
        s2.ue = 1.12;
        s2.theta = 0.0022;
        s2.dstar = 0.0052;
        s2.wake_gap = 0.0;
        s2.h = s2.dstar / s2.theta;
        s2.machsqd_edge = 0.0;
        s2.hstar = 1.76;
        s2.cf = 0.0024;
        s2.hk = 1.38;
        s2.retheta = 2460.0;
        s2.delta = 0.015;
        s2.sqrtctau = 0.014;
        s2.sqrtctaueq = 0.038;
        s2.us = 0.52;
        s2.cdiss = 0.00075;

        // Station 2 derivatives
        s2.h_d_theta = -s2.h / s2.theta;
        s2.h_d_dstar = 1.0 / s2.theta;
        s2.machsqd_edge_d_ue = 0.0;
        s2.machsqd_edge_d_machsqd = 0.0;
        s2.cf_d_theta = 0.0;
        s2.cf_d_dstar = 0.0;
        s2.cf_d_ue = 0.0;
        s2.cf_d_machsqd = 0.0;
        s2.cf_d_re = -s2.cf / 1e6;
        s2.hk_d_theta = 0.0;
        s2.hk_d_dstar = 0.0;
        s2.hk_d_ue = 0.0;
        s2.hk_d_machsqd = 0.0;
        s2.hstar_d_theta = 0.0;
        s2.hstar_d_dstar = 0.0;
        s2.hstar_d_ue = 0.0;
        s2.hstar_d_machsqd = 0.0;
        s2.hstar_d_re = 0.0;
        s2.delta_d_theta = 0.0;
        s2.delta_d_dstar = 0.0;
        s2.delta_d_ue = 0.0;
        s2.delta_d_machsqd = 0.0;
        s2.us_d_theta = 0.0;
        s2.us_d_dstar = 0.0;
        s2.us_d_ue = 0.0;
        s2.us_d_machsqd = 0.0;
        s2.us_d_re = 0.0;
        s2.sqrtctaueq_d_theta = 0.0;
        s2.sqrtctaueq_d_dstar = 0.0;
        s2.sqrtctaueq_d_ue = 0.0;
        s2.sqrtctaueq_d_machsqd = 0.0;
        s2.sqrtctaueq_d_re = 0.0;
        s2.cdiss_d_theta = 0.0;
        s2.cdiss_d_dstar = 0.0;
        s2.cdiss_d_ue = 0.0;
        s2.cdiss_d_sqrtctau = 0.0;
        s2.cdiss_d_machsqd = 0.0;
        s2.cdiss_d_re = 0.0;
        s2.hstarstar = 0.0;
        s2.hstarstar_d_theta = 0.0;
        s2.hstarstar_d_dstar = 0.0;
        s2.hstarstar_d_ue = 0.0;
        s2.hstarstar_d_machsqd = 0.0;
        s2.retheta_d_theta = 0.0;
        s2.retheta_d_ue = 0.0;
        s2.retheta_d_machsqd = 0.0;
        s2.retheta_d_re = 0.0;

        // Create CFM (simplified average)
        let mut cfm = MidpointCf::default();
        cfm.cf = 0.5 * (s1.cf + s2.cf);
        cfm.cf_d_re = 0.5 * (s1.cf_d_re + s2.cf_d_re);

        // Run BLDIF
        let mut sys = IntervalSystem::default();
        sys.assemble_interval_equations(
            &s1,
            &s2,
            &cfm,
            FlowRegime::Turbulent,
            false,
            9.0,
            AmplificationModel::Envelope,
        );

        // Check momentum equation residual (row 2)
        // VSREZ[2] = -0.7123274356e-1 from Fortran
        assert_relative_eq!(sys.residual[1], -0.07123274356, epsilon = 1e-4);

        // Check Jacobian entries for row 2 (momentum)
        // Note: indexing is [row][col] where col is 0=S, 1=T, 2=D, 3=U, 4=X
        // VS1[2,1] (dT1) = -0.5339051514e+3
        assert_relative_eq!(sys.jacobian_station1[1][1], -533.9, epsilon = 1.0);

        // VS1[2,2] (dD1) = -0.7342562675e+1
        assert_relative_eq!(sys.jacobian_station1[1][2], -7.34, epsilon = 0.1);

        // VS1[2,3] (dU1) = -0.3853754759e+1
        assert_relative_eq!(sys.jacobian_station1[1][3], -3.85, epsilon = 0.1);

        // VS2[2,1] (dT2) = 0.4716367493e+3
        assert_relative_eq!(sys.jacobian_station2[1][1], 471.6, epsilon = 1.0);

        // VS2[2,2] (dD2) = -0.6007551670e+1
        assert_relative_eq!(sys.jacobian_station2[1][2], -6.0, epsilon = 0.1);

        // VS2[2,3] (dU2) = 0.3956980467e+1
        assert_relative_eq!(sys.jacobian_station2[1][3], 3.96, epsilon = 0.1);
    }
}
