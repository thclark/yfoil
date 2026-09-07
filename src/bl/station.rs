//! The BL station state of XBL.INC's /V_VAR1/ and /V_VAR2/ blocks (COM1/COM2) with BLPRV,
//! BLKIN and BLVAR (xblsys.f), the midpoint skin friction BLMID, and DSLIM (xbl.f).

use super::closure::{
    cdiss_laminar, cdiss_wake, cf_laminar, cf_turbulent, hk_from_h, hstar_laminar, hstar_turbulent, hstarstar,
};
use super::params::*;

/// Mirrors XFOIL's `COM1`, `COM2`, `V_VAR1`, `V_VAR2`.
///
/// BL station state variables (from XFOIL's V_VAR1/V_VAR2)
///
/// This holds all primary and derived variables at a single BL station,
/// along with their derivatives with respect to the primary variables.
///
/// The naming follows XFOIL conventions where derivatives are denoted by
/// suffixes like _u (w.r.t. U), _t (w.r.t. θ), _d (w.r.t. δ*), etc.
#[derive(Debug, Clone, Default)]
#[doc(alias = "COM1")]
#[doc(alias = "COM2")]
#[doc(alias = "V_VAR1")]
#[doc(alias = "V_VAR2")]
pub struct StationState {
    // ========================================================================
    // Primary variables (set by blprv)
    // ========================================================================
    /// Arc length position (X2 in XFOIL)
    pub xi: f64,
    /// Edge velocity (compressible) (U2 in XFOIL)
    pub ue: f64,
    /// Momentum thickness θ (T2 in XFOIL)
    pub theta: f64,
    /// Displacement thickness δ* without wake gap (D2 in XFOIL)
    pub dstar: f64,
    /// Shear stress coefficient Ctau (turbulent) (S2 in XFOIL)
    pub sqrtctau: f64,
    /// Amplification factor (laminar) (AMPL2 in XFOIL)
    pub ampl: f64,
    /// Wake gap contribution to δ* (DW2 in XFOIL)
    pub wake_gap: f64,

    // Velocity transformation (compressible to incompressible)
    /// ∂U/∂Uei (compressible w.r.t. incompressible)
    pub ue_d_uei: f64,
    /// ∂U/∂M²
    pub ue_d_machsqd: f64,

    // ========================================================================
    // Kinematic secondary variables (set by blkin)
    // ========================================================================
    /// Shape factor H = δ*/θ
    pub h: f64,
    pub h_d_theta: f64, // ∂H/∂θ
    pub h_d_dstar: f64, // ∂H/∂δ*

    /// Edge Mach number squared M²
    pub machsqd_edge: f64,
    pub machsqd_edge_d_ue: f64,      // ∂M²/∂U
    pub machsqd_edge_d_machsqd: f64, // ∂M²/∂(M∞²)

    /// Density ratio ρ/ρ∞ (static to freestream)
    pub rho: f64,
    pub rho_d_ue: f64,
    pub rho_d_machsqd: f64,

    /// Kinematic viscosity ν ratio
    pub nu: f64,
    pub nu_d_ue: f64,
    pub nu_d_machsqd: f64,
    pub nu_d_re: f64,

    /// Kinematic shape factor Hk
    pub hk: f64,
    pub hk_d_ue: f64,      // ∂Hk/∂U
    pub hk_d_theta: f64,   // ∂Hk/∂θ
    pub hk_d_dstar: f64,   // ∂Hk/∂δ*
    pub hk_d_machsqd: f64, // ∂Hk/∂M²

    /// Momentum Reynolds number Rθ = ρ·Ue·θ/μ
    pub retheta: f64,
    pub retheta_d_ue: f64,
    pub retheta_d_theta: f64,
    pub retheta_d_machsqd: f64,
    pub retheta_d_re: f64,

    // ========================================================================
    // Turbulence-dependent secondary variables (set by blvar)
    // ========================================================================
    /// Density thickness shape factor H**
    pub hstarstar: f64,
    pub hstarstar_d_ue: f64,
    pub hstarstar_d_theta: f64,
    pub hstarstar_d_dstar: f64,
    pub hstarstar_d_machsqd: f64,

    /// Energy shape factor H*
    pub hstar: f64,
    pub hstar_d_ue: f64,
    pub hstar_d_theta: f64,
    pub hstar_d_dstar: f64,
    pub hstar_d_machsqd: f64,
    pub hstar_d_re: f64,

    /// Normalized slip velocity Us
    pub us: f64,
    pub us_d_ue: f64,
    pub us_d_theta: f64,
    pub us_d_dstar: f64,
    pub us_d_machsqd: f64,
    pub us_d_re: f64,

    /// Equilibrium shear stress coefficient CQ (Ctau^(1/2)_eq)
    pub sqrtctaueq: f64,
    pub sqrtctaueq_d_ue: f64,
    pub sqrtctaueq_d_theta: f64,
    pub sqrtctaueq_d_dstar: f64,
    pub sqrtctaueq_d_machsqd: f64,
    pub sqrtctaueq_d_re: f64,

    /// Skin friction coefficient Cf
    pub cf: f64,
    pub cf_d_ue: f64,
    pub cf_d_theta: f64,
    pub cf_d_dstar: f64,
    pub cf_d_machsqd: f64,
    pub cf_d_re: f64,

    /// Dissipation coefficient 2*CD/H*
    pub cdiss: f64,
    pub cdiss_d_ue: f64,
    pub cdiss_d_theta: f64,
    pub cdiss_d_dstar: f64,
    pub cdiss_d_sqrtctau: f64, // ∂Di/∂Ctau
    pub cdiss_d_machsqd: f64,
    pub cdiss_d_re: f64,

    /// BL thickness δ (from Green's correlation)
    pub delta: f64,
    pub delta_d_ue: f64,
    pub delta_d_theta: f64,
    pub delta_d_dstar: f64,
    pub delta_d_machsqd: f64,

    // ========================================================================
    // For convenience
    // ========================================================================
    /// Mass defect M = Ue·δ* (stored, not computed with derivatives here)
    pub mass_defect: f64,
}

impl StationState {
    /// Set primary BL variables from input parameters (BLPRV equivalent)
    ///
    /// This converts incompressible edge velocity Uei to compressible U
    /// using the Karman-Tsien transformation.
    ///
    /// # Arguments
    /// * `xsi` - Arc length position
    /// * `ami` - Amplification factor (laminar)
    /// * `cti` - Shear stress coefficient (turbulent)
    /// * `thi` - Momentum thickness θ
    /// * `dsi` - Total displacement thickness δ* (including wake gap)
    /// * `dswaki` - Wake gap contribution
    /// * `uei` - Edge velocity (incompressible)
    /// * `params` - Global BL parameters
    #[doc(alias = "BLPRV")]
    pub fn set_primary_variables(
        &mut self,
        xsi: f64,
        ami: f64,
        cti: f64,
        thi: f64,
        dsi: f64,
        dswaki: f64,
        uei: f64,
        params: &FlowParameters,
    ) {
        self.xi = xsi;
        self.ampl = ami;
        self.sqrtctau = cti;
        self.theta = thi;
        self.dstar = dsi - dswaki; // D2 is delta* without wake gap
        self.wake_gap = dswaki;

        // Karman-Tsien velocity transformation: Ue_compressible from Ue_incompressible
        // U2 = UEI*(1-TKBL) / (1 - TKBL*(UEI/QINFBL)^2)
        let tk = params.karman_tsien;
        let qinf = params.qinf;
        let uei_q = uei / qinf;
        let uei_q2 = uei_q * uei_q;
        let denom = 1.0 - tk * uei_q2;

        self.ue = uei * (1.0 - tk) / denom;

        // Derivative: ∂U/∂Uei
        // U2_UEI = (1 + TKBL*(2*U2*UEI/QINFBL^2 - 1)) / (1 - TKBL*(UEI/QINFBL)^2)
        self.ue_d_uei = (1.0 + tk * (2.0 * self.ue * uei / (qinf * qinf) - 1.0)) / denom;

        // Derivative: ∂U/∂M² (through TK)
        // U2_MS = (U2*(UEI/QINFBL)^2 - UEI) * TKBL_MS / denom
        self.ue_d_machsqd = (self.ue * uei_q2 - uei) * params.karman_tsien_d_machsqd / denom;
    }

    /// Calculate turbulence-independent secondary variables (BLKIN equivalent)
    ///
    /// This computes the kinematic BL quantities that don't depend on whether
    /// the flow is laminar or turbulent: M², density ratio, shape factor H,
    /// kinematic shape factor Hk, and momentum Reynolds number Rθ.
    ///
    /// # Arguments
    /// * `params` - Global BL parameters
    #[doc(alias = "BLKIN")]
    pub fn set_kinematic_variables(&mut self, params: &FlowParameters) {
        let u = self.ue;
        let t = self.theta;
        let d = self.dstar;
        let gm1 = params.gamma_gas_m1;
        let hstinv = params.h_stagnation_inv;
        let hstinv_ms = params.h_stagnation_inv_d_machsqd;
        let hvrat = params.sutherland_ratio;
        let rst = params.rho_stagnation;
        let rst_ms = params.rho_stagnation_d_machsqd;
        let reybl = params.re;
        let reybl_ms = params.re_d_machsqd;
        let reybl_re = params.re_d_re;

        // Edge Mach number squared
        // M² = U²·HSTINV / (γ-1)·(1 - 0.5·U²·HSTINV))
        let u2_hstinv = u * u * hstinv;
        let denom = gm1 * (1.0 - 0.5 * u2_hstinv);
        self.machsqd_edge = u2_hstinv / denom;

        let tr = 1.0 + 0.5 * gm1 * self.machsqd_edge;
        self.machsqd_edge_d_ue = 2.0 * self.machsqd_edge * tr / u;
        self.machsqd_edge_d_machsqd = u * u * tr / denom * hstinv_ms;

        // Edge static density (isentropic)
        // R = RST * TR^(-1/(γ-1))
        self.rho = rst * tr.powf(-1.0 / gm1);
        self.rho_d_ue = -self.rho / tr * 0.5 * self.machsqd_edge_d_ue;
        self.rho_d_machsqd = -self.rho / tr * 0.5 * self.machsqd_edge_d_machsqd + rst_ms * tr.powf(-1.0 / gm1);

        // Shape factor H = δ*/θ
        self.h = d / t;
        self.h_d_dstar = 1.0 / t;
        self.h_d_theta = -self.h / t;

        // Edge static/stagnation enthalpy ratio
        let herat = 1.0 - 0.5 * u * u * hstinv;
        let he_u = -u * hstinv;
        let he_ms = -0.5 * u * u * hstinv_ms;

        // Molecular viscosity ratio (Sutherland-type)
        // V = sqrt(HERAT^3) * (1+HVRAT)/(HERAT+HVRAT) / REYBL
        let herat_32 = (herat * herat * herat).sqrt(); // SQRT(HERAT**3): cube first, as the Fortran does
        self.nu = herat_32 * (1.0 + hvrat) / (herat + hvrat) / reybl;
        let v_he = self.nu * (1.5 / herat - 1.0 / (herat + hvrat));

        self.nu_d_ue = v_he * he_u;
        self.nu_d_machsqd = -self.nu / reybl * reybl_ms + v_he * he_ms;
        self.nu_d_re = -self.nu / reybl * reybl_re;

        // Kinematic shape factor Hk (compressibility correction)
        let (hk, hk_h, hk_msq) = hk_from_h(self.h, self.machsqd_edge);
        self.hk = hk;

        self.hk_d_ue = hk_msq * self.machsqd_edge_d_ue;
        self.hk_d_theta = hk_h * self.h_d_theta;
        self.hk_d_dstar = hk_h * self.h_d_dstar;
        self.hk_d_machsqd = hk_msq * self.machsqd_edge_d_machsqd;

        // Momentum thickness Reynolds number Rθ = ρ·U·θ/μ = R·U·T/V
        self.retheta = self.rho * u * t / self.nu;
        self.retheta_d_ue = self.retheta * (1.0 / u + self.rho_d_ue / self.rho - self.nu_d_ue / self.nu);
        self.retheta_d_theta = self.retheta / t;
        self.retheta_d_machsqd = self.retheta * (self.rho_d_machsqd / self.rho - self.nu_d_machsqd / self.nu);
        self.retheta_d_re = self.retheta * (-self.nu_d_re / self.nu);
    }

    /// Copy all values from another station state
    ///
    /// Used for the COM1 = COM2 operation in XFOIL
    pub fn copy_from(&mut self, other: &StationState) {
        *self = other.clone();
    }

    /// Calculate all secondary BL variables (BLVAR equivalent)
    ///
    /// This calculates the turbulence-dependent secondary variables
    /// based on the flow type (laminar, turbulent, wake).
    ///
    /// # Arguments
    /// * `flow_type` - Type of BL flow (Laminar, Turbulent, or Wake)
    /// * `params` - Global BL parameters
    #[doc(alias = "BLVAR")]
    pub fn set_closure_variables(&mut self, flow_type: FlowRegime, _params: &FlowParameters) {
        let hk = self.hk;
        let rt = self.retheta;
        let msq = self.machsqd_edge;
        let h = self.h;
        let t = self.theta;
        let d = self.dstar;
        let s = self.sqrtctau; // Ctau (or amplification for laminar)

        // XFOIL clamps HK2 itself in COMMON (BLVAR: `HK2 = MAX(HK2, ...)`) and leaves the Hk
        // derivatives as BLKIN set them. The clamped value persists: BLMID, the next station's
        // HK1 (after COM1 = COM2) and MRCHUE's HTARG all read it, so it is written back here.
        let hk = match flow_type {
            FlowRegime::Wake => hk.max(1.00005),
            _ => hk.max(1.05),
        };
        self.hk = hk;

        // ====================================================================
        // Density thickness shape parameter H** (from HCT)
        // ====================================================================
        let (hc, hc_hk, hc_msq) = hstarstar(hk, msq);

        self.hstarstar = hc;
        self.hstarstar_d_ue = hc_hk * self.hk_d_ue + hc_msq * self.machsqd_edge_d_ue;
        self.hstarstar_d_theta = hc_hk * self.hk_d_theta;
        self.hstarstar_d_dstar = hc_hk * self.hk_d_dstar;
        self.hstarstar_d_machsqd = hc_hk * self.hk_d_machsqd + hc_msq * self.machsqd_edge_d_machsqd;

        // ====================================================================
        // Energy shape factor H* (from HSL or HST)
        // ====================================================================
        let hs_result = match flow_type {
            FlowRegime::Laminar => hstar_laminar(hk, rt, msq),
            FlowRegime::Turbulent | FlowRegime::Wake => hstar_turbulent(hk, rt, msq),
        };

        self.hstar = hs_result.value;
        self.hstar_d_ue = hs_result.value_d_hk * self.hk_d_ue
            + hs_result.value_d_retheta * self.retheta_d_ue
            + hs_result.value_d_machsqd * self.machsqd_edge_d_ue;
        self.hstar_d_theta = hs_result.value_d_hk * self.hk_d_theta + hs_result.value_d_retheta * self.retheta_d_theta;
        self.hstar_d_dstar = hs_result.value_d_hk * self.hk_d_dstar;
        self.hstar_d_machsqd = hs_result.value_d_hk * self.hk_d_machsqd
            + hs_result.value_d_retheta * self.retheta_d_machsqd
            + hs_result.value_d_machsqd * self.machsqd_edge_d_machsqd;
        self.hstar_d_re = hs_result.value_d_retheta * self.retheta_d_re;

        // ---- normalized slip velocity  Us
        let us = 0.5 * self.hstar * (1.0 - (hk - 1.0) / (GBETA_LOCUS_B * h));
        let us_hs = 0.5 * (1.0 - (hk - 1.0) / (GBETA_LOCUS_B * h));
        let us_hk = 0.5 * self.hstar * (-1.0 / (GBETA_LOCUS_B * h));
        let us_h = 0.5 * self.hstar * (hk - 1.0) / (GBETA_LOCUS_B * (h * h));
        self.us = us;
        self.us_d_ue = us_hs * self.hstar_d_ue + us_hk * self.hk_d_ue;
        self.us_d_theta = us_hs * self.hstar_d_theta + us_hk * self.hk_d_theta + us_h * self.h_d_theta;
        self.us_d_dstar = us_hs * self.hstar_d_dstar + us_hk * self.hk_d_dstar + us_h * self.h_d_dstar;
        self.us_d_machsqd = us_hs * self.hstar_d_machsqd + us_hk * self.hk_d_machsqd;
        self.us_d_re = us_hs * self.hstar_d_re;
        if flow_type != FlowRegime::Wake && self.us > 0.95 {
            self.us = 0.98;
            self.us_d_ue = 0.0;
            self.us_d_theta = 0.0;
            self.us_d_dstar = 0.0;
            self.us_d_machsqd = 0.0;
            self.us_d_re = 0.0;
        }
        if flow_type == FlowRegime::Wake && self.us > 0.99995 {
            self.us = 0.99995;
            self.us_d_ue = 0.0;
            self.us_d_theta = 0.0;
            self.us_d_dstar = 0.0;
            self.us_d_machsqd = 0.0;
            self.us_d_re = 0.0;
        }
        let us = self.us;

        // ---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
        let mut hkc = hk - 1.0;
        let mut hkc_hk = 1.0;
        let mut hkc_rt = 0.0;
        if flow_type == FlowRegime::Turbulent {
            let gcc = GBETA_LOCUS_WALL;
            hkc = hk - 1.0 - gcc / rt;
            hkc_hk = 1.0;
            hkc_rt = gcc / (rt * rt);
            if hkc < 0.01 {
                hkc = 0.01;
                hkc_hk = 0.0;
                hkc_rt = 0.0;
            }
        }
        let hkb = hk - 1.0;
        let usb = 1.0 - us;
        let cq = (SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * (hkc * hkc) / (usb * h * (hk * hk))).sqrt();
        let cq_hs = SQRTCTAUEQ_COEFFICIENT * hkb * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq;
        let cq_us = SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / usb * 0.5 / cq;
        let cq_hk = SQRTCTAUEQ_COEFFICIENT * self.hstar * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq
            - SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * (hkc * hkc) / (usb * h * ((hk * hk) * hk)) * 2.0 * 0.5 / cq
            + SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_hk;
        let cq_rt = SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_rt;
        let cq_h = -(SQRTCTAUEQ_COEFFICIENT * self.hstar * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / h * 0.5 / cq);
        self.sqrtctaueq = cq;
        self.sqrtctaueq_d_ue = cq_hs * self.hstar_d_ue + cq_us * self.us_d_ue + cq_hk * self.hk_d_ue;
        self.sqrtctaueq_d_theta = cq_hs * self.hstar_d_theta + cq_us * self.us_d_theta + cq_hk * self.hk_d_theta;
        self.sqrtctaueq_d_dstar = cq_hs * self.hstar_d_dstar + cq_us * self.us_d_dstar + cq_hk * self.hk_d_dstar;
        self.sqrtctaueq_d_machsqd =
            cq_hs * self.hstar_d_machsqd + cq_us * self.us_d_machsqd + cq_hk * self.hk_d_machsqd;
        self.sqrtctaueq_d_re = cq_hs * self.hstar_d_re + cq_us * self.us_d_re;
        self.sqrtctaueq_d_ue = self.sqrtctaueq_d_ue + cq_rt * self.retheta_d_ue;
        self.sqrtctaueq_d_theta = self.sqrtctaueq_d_theta + cq_h * self.h_d_theta + cq_rt * self.retheta_d_theta;
        self.sqrtctaueq_d_dstar = self.sqrtctaueq_d_dstar + cq_h * self.h_d_dstar;
        self.sqrtctaueq_d_machsqd = self.sqrtctaueq_d_machsqd + cq_rt * self.retheta_d_machsqd;
        self.sqrtctaueq_d_re = self.sqrtctaueq_d_re + cq_rt * self.retheta_d_re;

        // ---- set skin friction coefficient
        let (cf, cf_hk, cf_rt, cf_m) = match flow_type {
            // wake
            FlowRegime::Wake => (0.0, 0.0, 0.0, 0.0),
            // laminar
            FlowRegime::Laminar => {
                let r = cf_laminar(hk, rt, msq);
                (r.value, r.value_d_hk, r.value_d_retheta, r.value_d_machsqd)
            }
            // turbulent
            FlowRegime::Turbulent => {
                let r = cf_turbulent(hk, rt, msq, CF_TURBULENT_FACTOR);
                let l = cf_laminar(hk, rt, msq);
                if l.value > r.value {
                    // laminar Cf is greater than turbulent Cf -- use laminar
                    // (this will only occur for unreasonably small Rtheta)
                    (l.value, l.value_d_hk, l.value_d_retheta, l.value_d_machsqd)
                } else {
                    (r.value, r.value_d_hk, r.value_d_retheta, r.value_d_machsqd)
                }
            }
        };
        self.cf = cf;
        self.cf_d_ue = cf_hk * self.hk_d_ue + cf_rt * self.retheta_d_ue + cf_m * self.machsqd_edge_d_ue;
        self.cf_d_theta = cf_hk * self.hk_d_theta + cf_rt * self.retheta_d_theta;
        self.cf_d_dstar = cf_hk * self.hk_d_dstar;
        self.cf_d_machsqd =
            cf_hk * self.hk_d_machsqd + cf_rt * self.retheta_d_machsqd + cf_m * self.machsqd_edge_d_machsqd;
        self.cf_d_re = cf_rt * self.retheta_d_re;

        // ---- dissipation function    2 CD / H*
        match flow_type {
            FlowRegime::Laminar => {
                // laminar
                let r = cdiss_laminar(hk, rt);
                self.cdiss = r.value;
                self.cdiss_d_ue = r.value_d_hk * self.hk_d_ue + r.value_d_retheta * self.retheta_d_ue;
                self.cdiss_d_theta = r.value_d_hk * self.hk_d_theta + r.value_d_retheta * self.retheta_d_theta;
                self.cdiss_d_dstar = r.value_d_hk * self.hk_d_dstar;
                self.cdiss_d_sqrtctau = 0.0;
                self.cdiss_d_machsqd = r.value_d_hk * self.hk_d_machsqd + r.value_d_retheta * self.retheta_d_machsqd;
                self.cdiss_d_re = r.value_d_retheta * self.retheta_d_re;
            }
            FlowRegime::Turbulent => {
                // turbulent wall contribution
                let c = cf_turbulent(hk, rt, msq, CF_TURBULENT_FACTOR);
                let cf2t = c.value;
                let cf2t_u = c.value_d_hk * self.hk_d_ue
                    + c.value_d_retheta * self.retheta_d_ue
                    + c.value_d_machsqd * self.machsqd_edge_d_ue;
                let cf2t_t = c.value_d_hk * self.hk_d_theta + c.value_d_retheta * self.retheta_d_theta;
                let cf2t_d = c.value_d_hk * self.hk_d_dstar;
                let cf2t_ms = c.value_d_hk * self.hk_d_machsqd
                    + c.value_d_retheta * self.retheta_d_machsqd
                    + c.value_d_machsqd * self.machsqd_edge_d_machsqd;
                let cf2t_re = c.value_d_retheta * self.retheta_d_re;
                let mut di = (0.5 * cf2t * us) * 2.0 / self.hstar;
                let di_hs = -((0.5 * cf2t * us) * 2.0 / (self.hstar * self.hstar));
                let di_us = (0.5 * cf2t) * 2.0 / self.hstar;
                let di_cf2t = (0.5 * us) * 2.0 / self.hstar;
                let mut di_s = 0.0;
                let mut di_u = di_hs * self.hstar_d_ue + di_us * self.us_d_ue + di_cf2t * cf2t_u;
                let mut di_t = di_hs * self.hstar_d_theta + di_us * self.us_d_theta + di_cf2t * cf2t_t;
                let mut di_d = di_hs * self.hstar_d_dstar + di_us * self.us_d_dstar + di_cf2t * cf2t_d;
                let mut di_ms = di_hs * self.hstar_d_machsqd + di_us * self.us_d_machsqd + di_cf2t * cf2t_ms;
                let mut di_re = di_hs * self.hstar_d_re + di_us * self.us_d_re + di_cf2t * cf2t_re;

                // set minimum Hk for wake layer to still exist
                let grt = rt.ln();
                let hmin = 1.0 + 2.1 / grt;
                let hm_rt = -(2.1 / (grt * grt)) / rt;

                // set factor DFAC for correcting wall dissipation for very low Hk
                let fl = (hk - 1.0) / (hmin - 1.0);
                let fl_hk = 1.0 / (hmin - 1.0);
                let fl_rt = (-fl / (hmin - 1.0)) * hm_rt;
                let tfl = fl.tanh();
                let dfac = 0.5 + 0.5 * tfl;
                let df_fl = 0.5 * (1.0 - tfl * tfl);
                let df_hk = df_fl * fl_hk;
                let df_rt = df_fl * fl_rt;

                di_s *= dfac;
                di_u = di_u * dfac + di * (df_hk * self.hk_d_ue + df_rt * self.retheta_d_ue);
                di_t = di_t * dfac + di * (df_hk * self.hk_d_theta + df_rt * self.retheta_d_theta);
                di_d = di_d * dfac + di * (df_hk * self.hk_d_dstar);
                di_ms = di_ms * dfac + di * (df_hk * self.hk_d_machsqd + df_rt * self.retheta_d_machsqd);
                di_re = di_re * dfac + di * (df_rt * self.retheta_d_re);
                di *= dfac;

                self.cdiss = di;
                self.cdiss_d_sqrtctau = di_s;
                self.cdiss_d_ue = di_u;
                self.cdiss_d_theta = di_t;
                self.cdiss_d_dstar = di_d;
                self.cdiss_d_machsqd = di_ms;
                self.cdiss_d_re = di_re;
            }
            FlowRegime::Wake => {
                // zero wall contribution for wake
                self.cdiss = 0.0;
                self.cdiss_d_sqrtctau = 0.0;
                self.cdiss_d_ue = 0.0;
                self.cdiss_d_theta = 0.0;
                self.cdiss_d_dstar = 0.0;
                self.cdiss_d_machsqd = 0.0;
                self.cdiss_d_re = 0.0;
            }
        }

        // ---- Add on turbulent outer layer contribution
        if flow_type != FlowRegime::Laminar {
            let dd = (s * s) * (0.995 - us) * 2.0 / self.hstar;
            let dd_hs = -((s * s) * (0.995 - us) * 2.0 / (self.hstar * self.hstar));
            let dd_us = -((s * s) * 2.0 / self.hstar);
            let dd_s = s * 2.0 * (0.995 - us) * 2.0 / self.hstar;
            self.cdiss = self.cdiss + dd;
            self.cdiss_d_sqrtctau = dd_s;
            self.cdiss_d_ue = self.cdiss_d_ue + dd_hs * self.hstar_d_ue + dd_us * self.us_d_ue;
            self.cdiss_d_theta = self.cdiss_d_theta + dd_hs * self.hstar_d_theta + dd_us * self.us_d_theta;
            self.cdiss_d_dstar = self.cdiss_d_dstar + dd_hs * self.hstar_d_dstar + dd_us * self.us_d_dstar;
            self.cdiss_d_machsqd = self.cdiss_d_machsqd + dd_hs * self.hstar_d_machsqd + dd_us * self.us_d_machsqd;
            self.cdiss_d_re = self.cdiss_d_re + dd_hs * self.hstar_d_re + dd_us * self.us_d_re;

            // add laminar stress contribution to outer layer CD
            let dd = 0.15 * ((0.995 - us) * (0.995 - us)) / rt * 2.0 / self.hstar;
            let dd_us = -0.15 * (0.995 - us) * 2.0 / rt * 2.0 / self.hstar;
            let dd_hs = -dd / self.hstar;
            let dd_rt = -dd / rt;
            self.cdiss = self.cdiss + dd;
            self.cdiss_d_ue =
                self.cdiss_d_ue + dd_hs * self.hstar_d_ue + dd_us * self.us_d_ue + dd_rt * self.retheta_d_ue;
            self.cdiss_d_theta = self.cdiss_d_theta
                + dd_hs * self.hstar_d_theta
                + dd_us * self.us_d_theta
                + dd_rt * self.retheta_d_theta;
            self.cdiss_d_dstar = self.cdiss_d_dstar + dd_hs * self.hstar_d_dstar + dd_us * self.us_d_dstar;
            self.cdiss_d_machsqd = self.cdiss_d_machsqd
                + dd_hs * self.hstar_d_machsqd
                + dd_us * self.us_d_machsqd
                + dd_rt * self.retheta_d_machsqd;
            self.cdiss_d_re =
                self.cdiss_d_re + dd_hs * self.hstar_d_re + dd_us * self.us_d_re + dd_rt * self.retheta_d_re;
        }

        if flow_type == FlowRegime::Turbulent {
            let l = cdiss_laminar(hk, rt);
            if l.value > self.cdiss {
                // laminar CD is greater than turbulent CD -- use laminar
                // (this will only occur for unreasonably small Rtheta)
                self.cdiss = l.value;
                self.cdiss_d_sqrtctau = 0.0;
                self.cdiss_d_ue = l.value_d_hk * self.hk_d_ue + l.value_d_retheta * self.retheta_d_ue;
                self.cdiss_d_theta = l.value_d_hk * self.hk_d_theta + l.value_d_retheta * self.retheta_d_theta;
                self.cdiss_d_dstar = l.value_d_hk * self.hk_d_dstar;
                self.cdiss_d_machsqd = l.value_d_hk * self.hk_d_machsqd + l.value_d_retheta * self.retheta_d_machsqd;
                self.cdiss_d_re = l.value_d_retheta * self.retheta_d_re;
            }
        }

        if flow_type == FlowRegime::Wake {
            // laminar wake CD
            let l = cdiss_wake(hk, rt);
            if l.value > self.cdiss {
                // laminar wake CD is greater than turbulent CD -- use laminar
                self.cdiss = l.value;
                self.cdiss_d_sqrtctau = 0.0;
                self.cdiss_d_ue = l.value_d_hk * self.hk_d_ue + l.value_d_retheta * self.retheta_d_ue;
                self.cdiss_d_theta = l.value_d_hk * self.hk_d_theta + l.value_d_retheta * self.retheta_d_theta;
                self.cdiss_d_dstar = l.value_d_hk * self.hk_d_dstar;
                self.cdiss_d_machsqd = l.value_d_hk * self.hk_d_machsqd + l.value_d_retheta * self.retheta_d_machsqd;
                self.cdiss_d_re = l.value_d_retheta * self.retheta_d_re;
            }
        }

        if flow_type == FlowRegime::Wake {
            // double dissipation for the wake (two wake halves)
            self.cdiss *= 2.0;
            self.cdiss_d_sqrtctau *= 2.0;
            self.cdiss_d_ue *= 2.0;
            self.cdiss_d_theta *= 2.0;
            self.cdiss_d_dstar *= 2.0;
            self.cdiss_d_machsqd *= 2.0;
            self.cdiss_d_re *= 2.0;
        }

        // ====================================================================
        // BL thickness Delta (Green's correlation)
        // DE = (3.15 + 1.72/(HK-1)) * T + D
        // ====================================================================
        let de = (3.15 + 1.72 / (hk - 1.0)) * t + d;
        let de_hk = -1.72 / ((hk - 1.0) * (hk - 1.0)) * t;

        self.delta = de;
        self.delta_d_ue = de_hk * self.hk_d_ue;
        self.delta_d_theta = de_hk * self.hk_d_theta + 3.15 + 1.72 / (hk - 1.0);
        self.delta_d_dstar = de_hk * self.hk_d_dstar + 1.0;
        self.delta_d_machsqd = de_hk * self.hk_d_machsqd;

        // Clamp DE to reasonable values
        let hdmax = 12.0;
        if self.delta > hdmax * t {
            self.delta = hdmax * t;
            self.delta_d_ue = 0.0;
            self.delta_d_theta = hdmax;
            self.delta_d_dstar = 0.0;
            self.delta_d_machsqd = 0.0;
        }
    }
}

// ============================================================================
// Midpoint Skin Friction (BLMID)
// ============================================================================

/// Mirrors XFOIL's `CFM`.
///
/// Midpoint skin friction result (from XFOIL's BLMID)
///
/// This holds the skin friction coefficient at the midpoint between two
/// stations, along with its derivatives.
#[derive(Debug, Clone, Default)]
#[doc(alias = "CFM")]
pub struct MidpointCf {
    /// Midpoint skin friction coefficient
    pub cf: f64,

    /// Derivatives w.r.t. station 1
    pub cf_d_ue_station1: f64,
    pub cf_d_theta_station1: f64,
    pub cf_d_dstar_station1: f64,

    /// Derivatives w.r.t. station 2
    pub cf_d_ue_station2: f64,
    pub cf_d_theta_station2: f64,
    pub cf_d_dstar_station2: f64,

    /// Derivatives w.r.t. global parameters
    pub cf_d_machsqd: f64,
    pub cf_d_re: f64,
}

impl MidpointCf {
    /// Calculate midpoint skin friction (BLMID equivalent)
    ///
    /// Calculates the skin friction coefficient at the midpoint between
    /// two BL stations. For turbulent flow, uses the maximum of turbulent
    /// and laminar Cf.
    ///
    /// # Arguments
    /// * `s1` - Station 1 state
    /// * `s2` - Station 2 state
    /// * `flow_type` - Type of BL flow
    /// * `is_similarity` - True if this is a similarity station (copy s2→s1)
    #[doc(alias = "BLMID")]
    pub fn compute(s1: &StationState, s2: &StationState, flow_type: FlowRegime, is_similarity: bool) -> Self {
        let mut result = Self::default();

        // For similarity station, station 1 equals station 2
        let (hk1, rt1, m1) = if is_similarity {
            (s2.hk, s2.retheta, s2.machsqd_edge)
        } else {
            (s1.hk, s1.retheta, s1.machsqd_edge)
        };

        let (hk1_u1, hk1_t1, hk1_d1, hk1_ms) = if is_similarity {
            (s2.hk_d_ue, s2.hk_d_theta, s2.hk_d_dstar, s2.hk_d_machsqd)
        } else {
            (s1.hk_d_ue, s1.hk_d_theta, s1.hk_d_dstar, s1.hk_d_machsqd)
        };

        let (rt1_u1, rt1_t1, rt1_ms, rt1_re) = if is_similarity {
            (
                s2.retheta_d_ue,
                s2.retheta_d_theta,
                s2.retheta_d_machsqd,
                s2.retheta_d_re,
            )
        } else {
            (
                s1.retheta_d_ue,
                s1.retheta_d_theta,
                s1.retheta_d_machsqd,
                s1.retheta_d_re,
            )
        };

        let (m1_u1, m1_ms) = if is_similarity {
            (s2.machsqd_edge_d_ue, s2.machsqd_edge_d_machsqd)
        } else {
            (s1.machsqd_edge_d_ue, s1.machsqd_edge_d_machsqd)
        };

        // Midpoint averages
        let hka = 0.5 * (hk1 + s2.hk);
        let rta = 0.5 * (rt1 + s2.retheta);
        let ma = 0.5 * (m1 + s2.machsqd_edge);

        // Midpoint skin friction coefficient
        let (cfm, cfm_hka, cfm_rta, cfm_ma) = match flow_type {
            FlowRegime::Wake => {
                // Zero skin friction in wake
                (0.0, 0.0, 0.0, 0.0)
            }
            FlowRegime::Laminar => {
                // Laminar Cf
                let cf_result = cf_laminar(hka, rta, ma);
                (cf_result.value, cf_result.value_d_hk, cf_result.value_d_retheta, 0.0)
            }
            FlowRegime::Turbulent => {
                // Turbulent Cf
                let cf_turb_result = cf_turbulent(hka, rta, ma, CF_TURBULENT_FACTOR);
                // Check if laminar is higher
                let cf_lam_result = cf_laminar(hka, rta, ma);

                if cf_lam_result.value > cf_turb_result.value {
                    (
                        cf_lam_result.value,
                        cf_lam_result.value_d_hk,
                        cf_lam_result.value_d_retheta,
                        0.0,
                    )
                } else {
                    (
                        cf_turb_result.value,
                        cf_turb_result.value_d_hk,
                        cf_turb_result.value_d_retheta,
                        cf_turb_result.value_d_machsqd,
                    )
                }
            }
        };

        result.cf = cfm;

        // Derivatives w.r.t. station 1 variables (factor of 0.5 from averaging)
        result.cf_d_ue_station1 = 0.5 * (cfm_hka * hk1_u1 + cfm_ma * m1_u1 + cfm_rta * rt1_u1);
        result.cf_d_theta_station1 = 0.5 * (cfm_hka * hk1_t1 + cfm_rta * rt1_t1);
        result.cf_d_dstar_station1 = 0.5 * cfm_hka * hk1_d1;

        // Derivatives w.r.t. station 2 variables
        result.cf_d_ue_station2 =
            0.5 * (cfm_hka * s2.hk_d_ue + cfm_ma * s2.machsqd_edge_d_ue + cfm_rta * s2.retheta_d_ue);
        result.cf_d_theta_station2 = 0.5 * (cfm_hka * s2.hk_d_theta + cfm_rta * s2.retheta_d_theta);
        result.cf_d_dstar_station2 = 0.5 * cfm_hka * s2.hk_d_dstar;

        // Derivatives w.r.t. global parameters
        result.cf_d_machsqd = 0.5
            * (cfm_hka * hk1_ms
                + cfm_ma * m1_ms
                + cfm_rta * rt1_ms
                + cfm_hka * s2.hk_d_machsqd
                + cfm_ma * s2.machsqd_edge_d_machsqd
                + cfm_rta * s2.retheta_d_machsqd);
        result.cf_d_re = 0.5 * (cfm_rta * rt1_re + cfm_rta * s2.retheta_d_re);

        result
    }
}

/// Limit displacement thickness to keep Hk above minimum (DSLIM equivalent)
///
/// Adjusts δ* to prevent kinematic shape factor from dropping below HKLIM.
///
/// # Arguments
/// * `dstr` - Displacement thickness (modified in place)
/// * `thet` - Momentum thickness
/// * `uedg` - Edge velocity
/// * `msq` - Edge Mach number squared
/// * `hklim` - Minimum kinematic shape factor
#[doc(alias = "DSLIM")]
pub fn limit_dstar(dstr: &mut f64, thet: f64, _uedg: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet;
    let (hk, hk_h, _hk_m) = hk_from_h(h, msq);

    // If Hk is below limit, adjust δ*
    let dh = (hklim - hk).max(0.0) / hk_h;
    *dstr += dh * thet;
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_station_state_default() {
        let state = StationState::default();
        assert_eq!(state.theta, 0.0);
        assert_eq!(state.cf, 0.0);
    }
    // ========================================================================
    // BLPRV Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blprv_incompressible() {
        // Test case: M=0, Re=1e6
        // XFOIL reference values from Fortran test
        let params = FlowParameters::incompressible(1e6);
        let mut state = StationState::default();

        // Inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // Check primary variables
        assert_eq!(state.xi, 0.1);
        assert_eq!(state.ampl, 3.0);
        assert_eq!(state.sqrtctau, 0.015);
        assert_eq!(state.theta, 0.002);
        assert_eq!(state.dstar, 0.005);
        assert_eq!(state.wake_gap, 0.0);

        // At M=0, U2 = Uei (no transformation)
        assert_relative_eq!(state.ue, 1.2, epsilon = 1e-10);
        assert_relative_eq!(state.ue_d_uei, 1.0, epsilon = 1e-10);
        // Note: u_ms is the sensitivity d(U2)/d(M²). Even at M=0 this is non-zero because
        // COMSET's TKL_MSQ = 1/(1+BETA)² = 0.25 at M=0 (BETA = 1).
        // U2_MS = (U2*UEI² - UEI) * TKBL_MS = (1.2*1.44 - 1.2) * 0.25 = 0.132
        assert_relative_eq!(state.ue_d_machsqd, 0.132, epsilon = 1e-6);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (derived with TKBL = 1/beta - 1, not COMSET's TKLAM) — regenerate from the M=0.3 coverage case"]
    fn test_blprv_compressible() {
        // Test case: M=0.5, Re=1e6
        // XFOIL reference values from Fortran test
        let params = FlowParameters::new(0.5, 1e6, 1.4);
        let mut state = StationState::default();

        // Same inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.set_primary_variables(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // XFOIL reference: U2 = 1.305093527
        assert_relative_eq!(state.ue, 1.305093527, epsilon = 1e-5);
        // XFOIL reference: U2_UEI = 1.711017609
        assert_relative_eq!(state.ue_d_uei, 1.711017609, epsilon = 1e-5);
        // XFOIL reference: U2_MS = 0.6728398800
        assert_relative_eq!(state.ue_d_machsqd, 0.6728398800, epsilon = 1e-5);
    }

    // ========================================================================
    // BLKIN Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blkin_incompressible() {
        // Test case: M=0, Re=1e6
        let params = FlowParameters::incompressible(1e6);
        let mut state = StationState::default();

        // Run BLPRV first
        state.set_primary_variables(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.set_kinematic_variables(&params);

        // XFOIL reference values (M=0 case)
        // M2 = 0
        assert_relative_eq!(state.machsqd_edge, 0.0, epsilon = 1e-10);

        // R2 = 1.0
        assert_relative_eq!(state.rho, 1.0, epsilon = 1e-10);

        // H2 = D2/T2 = 0.005/0.002 = 2.5
        assert_relative_eq!(state.h, 2.5, epsilon = 1e-10);
        assert_relative_eq!(state.h_d_theta, -1250.0, epsilon = 1e-6); // -H/T
        assert_relative_eq!(state.h_d_dstar, 500.0, epsilon = 1e-6); // 1/T

        // HK2 = 2.5 (same as H at M=0)
        assert_relative_eq!(state.hk, 2.5, epsilon = 1e-6);

        // RT2 = 2400 (R*U*T/V = 1*1.2*0.002/1e-6)
        assert_relative_eq!(state.retheta, 2400.0, epsilon = 1.0);

        // RT2_U2 = 2000 (RT/U = 2400/1.2)
        assert_relative_eq!(state.retheta_d_ue, 2000.0, epsilon = 1.0);

        // RT2_T2 = 1200000 (RT/T = 2400/0.002)
        assert_relative_eq!(state.retheta_d_theta, 1200000.0, epsilon = 100.0);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (assumed HVRAT=0.35; XFOIL's analysis path leaves HVRAT=0) — regenerate from the M=0.3 coverage case"]
    fn test_blkin_compressible() {
        // Test case: M=0.5, Re=1e6
        let params = FlowParameters::new(0.5, 1e6, 1.4);
        let mut state = StationState::default();

        // Run BLPRV first
        state.set_primary_variables(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.set_kinematic_variables(&params);

        // XFOIL reference values (M=0.5 case)
        // M2 = 0.4413362443
        assert_relative_eq!(state.machsqd_edge, 0.4413362443, epsilon = 1e-5);

        // R2 = 0.9143959880
        assert_relative_eq!(state.rho, 0.9143959880, epsilon = 1e-5);

        // HK2 = 2.259336710
        assert_relative_eq!(state.hk, 2.259336710, epsilon = 1e-5);

        // RT2 = 2453.646240
        assert_relative_eq!(state.retheta, 2453.646240, epsilon = 1.0);
    }

    #[test]
    fn test_blkin_shape_factor_derivatives() {
        // Verify shape factor derivatives are computed correctly
        let params = FlowParameters::incompressible(1e6);
        let mut state = StationState::default();

        state.set_primary_variables(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);
        state.set_kinematic_variables(&params);

        // H = D/T = 2.5
        // dH/dD = 1/T = 1/0.002 = 500
        // dH/dT = -D/T² = -0.005/0.002² = -1250
        assert_relative_eq!(state.h_d_dstar, 1.0 / 0.002, epsilon = 1e-10);
        assert_relative_eq!(state.h_d_theta, -0.005 / 0.002_f64.powi(2), epsilon = 1e-10);
    }

    // ========================================================================
    // BLVAR Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blvar_laminar() {
        // Test case: Laminar BL
        // XFOIL reference values from Fortran test
        let params = FlowParameters::incompressible(1e6);
        let mut state = StationState::default();

        // Set up: HK=2.5, RT=2400, M=0, H=2.5, T=0.002, D=0.005, S=3.0
        // We need to set the kinematic variables directly since blprv/blkin
        // uses its own computations
        state.hk = 2.5;
        state.retheta = 2400.0;
        state.machsqd_edge = 0.0;
        state.h = 2.5;
        state.theta = 0.002;
        state.dstar = 0.005;
        state.sqrtctau = 3.0; // Amplification for laminar

        // Set the derivatives (simplified - zero for test)
        state.hk_d_ue = 0.0;
        state.hk_d_theta = -1250.0;
        state.hk_d_dstar = 500.0;
        state.hk_d_machsqd = 0.0;
        state.retheta_d_ue = 2000.0;
        state.retheta_d_theta = 1200000.0;
        state.retheta_d_machsqd = 0.0;
        state.retheta_d_re = 0.0024;
        state.h_d_theta = -1250.0;
        state.h_d_dstar = 500.0;
        state.machsqd_edge_d_ue = 0.0;
        state.machsqd_edge_d_machsqd = 0.0;

        state.set_closure_variables(FlowRegime::Laminar, &params);

        // XFOIL reference values:
        // HC2 = 0.0 (M=0)
        assert_relative_eq!(state.hstarstar, 0.0, epsilon = 1e-10);

        // HS2 = 1.584867239
        assert_relative_eq!(state.hstar, 1.584867239, epsilon = 1e-5);

        // US2 = 0.1584867090
        assert_relative_eq!(state.us, 0.1584867090, epsilon = 1e-5);

        // CQ2 = 0.07772709429
        assert_relative_eq!(state.sqrtctaueq, 0.07772709429, epsilon = 1e-4);

        // CF2 = 0.2045119036e-3
        assert_relative_eq!(state.cf, 0.2045119036e-3, epsilon = 1e-7);

        // DI2 = 0.9419409616e-4
        assert_relative_eq!(state.cdiss, 0.9419409616e-4, epsilon = 1e-7);

        // DE2 = 0.01359333377
        assert_relative_eq!(state.delta, 0.01359333377, epsilon = 1e-6);
    }

    #[test]
    fn test_blvar_turbulent() {
        // Test case: Turbulent BL
        // XFOIL reference values from Fortran test
        let params = FlowParameters::incompressible(1e6);
        let mut state = StationState::default();

        // Set up: HK=1.4, RT=10000, M=0, H=1.4, T=0.005, D=0.007, S=0.015
        state.hk = 1.4;
        state.retheta = 10000.0;
        state.machsqd_edge = 0.0;
        state.h = 1.4;
        state.theta = 0.005;
        state.dstar = 0.007;
        state.sqrtctau = 0.015;

        // Set derivatives (simplified)
        state.hk_d_ue = 0.0;
        state.hk_d_theta = -280.0;
        state.hk_d_dstar = 200.0;
        state.hk_d_machsqd = 0.0;
        state.retheta_d_ue = 0.0;
        state.retheta_d_theta = 2000000.0;
        state.retheta_d_machsqd = 0.0;
        state.retheta_d_re = 0.01;
        state.h_d_theta = -280.0;
        state.h_d_dstar = 200.0;
        state.machsqd_edge_d_ue = 0.0;
        state.machsqd_edge_d_machsqd = 0.0;
        state.hstar_d_ue = 0.0;
        state.hstar_d_theta = 0.0;
        state.hstar_d_dstar = 0.0;
        state.hstar_d_machsqd = 0.0;
        state.hstar_d_re = 0.0;
        state.us_d_ue = 0.0;
        state.us_d_theta = 0.0;
        state.us_d_dstar = 0.0;
        state.us_d_machsqd = 0.0;
        state.us_d_re = 0.0;

        state.set_closure_variables(FlowRegime::Turbulent, &params);

        // XFOIL reference values:
        // HS2 = 1.755310297
        assert_relative_eq!(state.hstar, 1.755310297, epsilon = 1e-4);

        // US2 = 0.5433103442
        assert_relative_eq!(state.us, 0.5433103442, epsilon = 1e-4);

        // CQ2 = 0.03632329032
        assert_relative_eq!(state.sqrtctaueq, 0.03632329032, epsilon = 1e-4);

        // CF2 = 0.2286923816e-2
        assert_relative_eq!(state.cf, 0.2286923816e-2, epsilon = 1e-5);

        // DI2 = 0.8065673755e-3 (turbulent DI is more complex, allow larger tolerance)
        assert_relative_eq!(state.cdiss, 0.8065673755e-3, epsilon = 1e-4);

        // DE2 = 0.04425000027
        assert_relative_eq!(state.delta, 0.04425000027, epsilon = 1e-5);
    }

    // ========================================================================
    // BLMID Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blmid_laminar() {
        // Test case: Laminar (ITYP=1)
        // XFOIL reference values from Fortran test
        // HKA = 2.4, RTA = 2200, MA = 0

        // Station 1
        let mut s1 = StationState::default();
        s1.hk = 2.3;
        s1.retheta = 2000.0;
        s1.machsqd_edge = 0.0;
        s1.hk_d_ue = 0.0;
        s1.hk_d_theta = -1150.0;
        s1.hk_d_dstar = 500.0;
        s1.hk_d_machsqd = -0.29;
        s1.retheta_d_ue = 1666.67;
        s1.retheta_d_theta = 1000000.0;
        s1.retheta_d_machsqd = 0.0;
        s1.retheta_d_re = 0.002;
        s1.machsqd_edge_d_ue = 0.0;
        s1.machsqd_edge_d_machsqd = 0.0;

        // Station 2
        let mut s2 = StationState::default();
        s2.hk = 2.5;
        s2.retheta = 2400.0;
        s2.machsqd_edge = 0.0;
        s2.hk_d_ue = 0.0;
        s2.hk_d_theta = -1250.0;
        s2.hk_d_dstar = 500.0;
        s2.hk_d_machsqd = -0.29;
        s2.retheta_d_ue = 2000.0;
        s2.retheta_d_theta = 1200000.0;
        s2.retheta_d_machsqd = 0.0;
        s2.retheta_d_re = 0.0024;
        s2.machsqd_edge_d_ue = 0.0;
        s2.machsqd_edge_d_machsqd = 0.0;

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Laminar, false);

        // XFOIL reference values
        assert_relative_eq!(result.cf, 0.2577280102e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cf_d_ue_station1, -0.9762444097e-4, epsilon = 1e-7);
        assert_relative_eq!(result.cf_d_theta_station1, 0.1515112668, epsilon = 1e-3);
        assert_relative_eq!(result.cf_d_dstar_station1, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cf_d_ue_station2, -0.1171490949e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cf_d_theta_station2, 0.1580646932, epsilon = 1e-3);
        assert_relative_eq!(result.cf_d_dstar_station2, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cf_d_machsqd, 0.1059563219e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cf_d_re, -0.2577280056e-9, epsilon = 1e-12);
    }

    #[test]
    fn test_blmid_turbulent() {
        // Test case: Turbulent (ITYP=2)
        // XFOIL reference values from Fortran test
        // HKA = 1.375, RTA = 9000, MA = 0

        // Station 1
        let mut s1 = StationState::default();
        s1.hk = 1.35;
        s1.retheta = 8000.0;
        s1.machsqd_edge = 0.0;
        s1.hk_d_ue = 0.0;
        s1.hk_d_theta = -270.0;
        s1.hk_d_dstar = 200.0;
        s1.hk_d_machsqd = -0.29;
        s1.retheta_d_ue = 0.0;
        s1.retheta_d_theta = 1600000.0;
        s1.retheta_d_machsqd = 0.0;
        s1.retheta_d_re = 0.008;
        s1.machsqd_edge_d_ue = 0.0;
        s1.machsqd_edge_d_machsqd = 0.0;

        // Station 2
        let mut s2 = StationState::default();
        s2.hk = 1.40;
        s2.retheta = 10000.0;
        s2.machsqd_edge = 0.0;
        s2.hk_d_ue = 0.0;
        s2.hk_d_theta = -280.0;
        s2.hk_d_dstar = 200.0;
        s2.hk_d_machsqd = -0.29;
        s2.retheta_d_ue = 0.0;
        s2.retheta_d_theta = 2000000.0;
        s2.retheta_d_machsqd = 0.0;
        s2.retheta_d_re = 0.01;
        s2.machsqd_edge_d_ue = 0.0;
        s2.machsqd_edge_d_machsqd = 0.0;

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Turbulent, false);

        // XFOIL reference values (turbulent Cf used)
        assert_relative_eq!(result.cf, 0.2450317144e-2, epsilon = 1e-5);
        assert_relative_eq!(result.cf_d_theta_station1, 0.5299984217, epsilon = 1e-2);
        assert_relative_eq!(result.cf_d_dstar_station1, -0.4310033321, epsilon = 1e-3);
        assert_relative_eq!(result.cf_d_theta_station2, 0.5385844707, epsilon = 1e-2);
        assert_relative_eq!(result.cf_d_dstar_station2, -0.4310033321, epsilon = 1e-3);
    }

    #[test]
    fn test_blmid_wake() {
        // Test case: Wake (ITYP=3)
        // Wake should have zero skin friction
        let s1 = StationState::default();
        let s2 = StationState::default();

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Wake, false);

        assert_eq!(result.cf, 0.0);
        assert_eq!(result.cf_d_ue_station1, 0.0);
        assert_eq!(result.cf_d_theta_station1, 0.0);
        assert_eq!(result.cf_d_dstar_station1, 0.0);
    }
    // ========================================================================
    // DSLIM Tests
    // ========================================================================

    #[test]
    fn test_dslim_no_change_needed() {
        // Test DSLIM when Hk is already above limit
        let mut dstr = 0.005;
        let thet = 0.002;
        let uedg = 1.0;
        let msq = 0.0;
        let hklim = 1.02;

        // H = 2.5, Hk = 2.5 at M=0 (well above 1.02)
        limit_dstar(&mut dstr, thet, uedg, msq, hklim);

        // Should not change
        assert_relative_eq!(dstr, 0.005, epsilon = 1e-10);
    }

    #[test]
    fn test_dslim_adjustment_needed() {
        // Test DSLIM when Hk is below limit
        let mut dstr = 0.00102; // H = 1.02/2 = 0.51 initially
        let thet = 0.001;
        let uedg = 1.0;
        let msq = 0.0;
        let hklim = 1.5; // Hk = H at M=0, so we need H >= 1.5

        // H = 1.02, Hk = 1.02 at M=0 (below 1.5)
        limit_dstar(&mut dstr, thet, uedg, msq, hklim);

        // Should increase dstr to raise Hk to at least hklim
        let new_h = dstr / thet;
        assert!(new_h >= 1.5, "H should be raised to at least hklim");
    }
}
