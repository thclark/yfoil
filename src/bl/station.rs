//! The BL station state of XBL.INC's /V_VAR1/ and /V_VAR2/ blocks (COM1/COM2) with BLPRV,
//! BLKIN and BLVAR (xblsys.f), the midpoint skin friction BLMID, and DSLIM (xbl.f).

use super::closure::{cf_lam, cf_turb, di_lam, dilw, hc_turb, hkin, hs_lam, hs_turb};
use super::params::*;

/// BL station state variables (from XFOIL's V_VAR1/V_VAR2)
///
/// This holds all primary and derived variables at a single BL station,
/// along with their derivatives with respect to the primary variables.
///
/// The naming follows XFOIL conventions where derivatives are denoted by
/// suffixes like _u (w.r.t. U), _t (w.r.t. θ), _d (w.r.t. δ*), etc.
#[derive(Debug, Clone, Default)]
pub struct BLStationState {
    // ========================================================================
    // Primary variables (set by blprv)
    // ========================================================================
    /// Arc length position (X2 in XFOIL)
    pub x: f64,
    /// Edge velocity (compressible) (U2 in XFOIL)
    pub u: f64,
    /// Momentum thickness θ (T2 in XFOIL)
    pub theta: f64,
    /// Displacement thickness δ* without wake gap (D2 in XFOIL)
    pub dstar: f64,
    /// Shear stress coefficient Ctau (turbulent) (S2 in XFOIL)
    pub ctau: f64,
    /// Amplification factor (laminar) (AMPL2 in XFOIL)
    pub ampl: f64,
    /// Wake gap contribution to δ* (DW2 in XFOIL)
    pub dw: f64,

    // Velocity transformation (compressible to incompressible)
    /// ∂U/∂Uei (compressible w.r.t. incompressible)
    pub u_uei: f64,
    /// ∂U/∂M²
    pub u_ms: f64,

    // ========================================================================
    // Kinematic secondary variables (set by blkin)
    // ========================================================================
    /// Shape factor H = δ*/θ
    pub h: f64,
    pub h_t: f64, // ∂H/∂θ
    pub h_d: f64, // ∂H/∂δ*

    /// Edge Mach number squared M²
    pub msq: f64,
    pub msq_u: f64,  // ∂M²/∂U
    pub msq_ms: f64, // ∂M²/∂(M∞²)

    /// Density ratio ρ/ρ∞ (static to freestream)
    pub r: f64,
    pub r_u: f64,
    pub r_ms: f64,

    /// Kinematic viscosity ν ratio
    pub v: f64,
    pub v_u: f64,
    pub v_ms: f64,
    pub v_re: f64,

    /// Kinematic shape factor Hk
    pub hk: f64,
    pub hk_u: f64,  // ∂Hk/∂U
    pub hk_t: f64,  // ∂Hk/∂θ
    pub hk_d: f64,  // ∂Hk/∂δ*
    pub hk_ms: f64, // ∂Hk/∂M²

    /// Momentum Reynolds number Rθ = ρ·Ue·θ/μ
    pub rt: f64,
    pub rt_u: f64,
    pub rt_t: f64,
    pub rt_ms: f64,
    pub rt_re: f64,

    // ========================================================================
    // Turbulence-dependent secondary variables (set by blvar)
    // ========================================================================
    /// Density thickness shape factor H**
    pub hc: f64,
    pub hc_u: f64,
    pub hc_t: f64,
    pub hc_d: f64,
    pub hc_ms: f64,

    /// Energy shape factor H*
    pub hs: f64,
    pub hs_u: f64,
    pub hs_t: f64,
    pub hs_d: f64,
    pub hs_ms: f64,
    pub hs_re: f64,

    /// Normalized slip velocity Us
    pub us: f64,
    pub us_u: f64,
    pub us_t: f64,
    pub us_d: f64,
    pub us_ms: f64,
    pub us_re: f64,

    /// Equilibrium shear stress coefficient CQ (Ctau^(1/2)_eq)
    pub cq: f64,
    pub cq_u: f64,
    pub cq_t: f64,
    pub cq_d: f64,
    pub cq_ms: f64,
    pub cq_re: f64,

    /// Skin friction coefficient Cf
    pub cf: f64,
    pub cf_u: f64,
    pub cf_t: f64,
    pub cf_d: f64,
    pub cf_ms: f64,
    pub cf_re: f64,

    /// Dissipation coefficient 2*CD/H*
    pub di: f64,
    pub di_u: f64,
    pub di_t: f64,
    pub di_d: f64,
    pub di_s: f64, // ∂Di/∂Ctau
    pub di_ms: f64,
    pub di_re: f64,

    /// BL thickness δ (from Green's correlation)
    pub de: f64,
    pub de_u: f64,
    pub de_t: f64,
    pub de_d: f64,
    pub de_ms: f64,

    // ========================================================================
    // For convenience
    // ========================================================================
    /// Mass defect M = Ue·δ* (stored, not computed with derivatives here)
    pub mass: f64,
}

impl BLStationState {
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
    pub fn blprv(
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
        self.x = xsi;
        self.ampl = ami;
        self.ctau = cti;
        self.theta = thi;
        self.dstar = dsi - dswaki; // D2 is delta* without wake gap
        self.dw = dswaki;

        // Karman-Tsien velocity transformation: Ue_compressible from Ue_incompressible
        // U2 = UEI*(1-TKBL) / (1 - TKBL*(UEI/QINFBL)^2)
        let tk = params.karman_tsien;
        let qinf = params.qinf;
        let uei_q = uei / qinf;
        let uei_q2 = uei_q * uei_q;
        let denom = 1.0 - tk * uei_q2;

        self.u = uei * (1.0 - tk) / denom;

        // Derivative: ∂U/∂Uei
        // U2_UEI = (1 + TKBL*(2*U2*UEI/QINFBL^2 - 1)) / (1 - TKBL*(UEI/QINFBL)^2)
        self.u_uei = (1.0 + tk * (2.0 * self.u * uei / (qinf * qinf) - 1.0)) / denom;

        // Derivative: ∂U/∂M² (through TK)
        // U2_MS = (U2*(UEI/QINFBL)^2 - UEI) * TKBL_MS / denom
        self.u_ms = (self.u * uei_q2 - uei) * params.karman_tsien_d_machsqd / denom;
    }

    /// Calculate turbulence-independent secondary variables (BLKIN equivalent)
    ///
    /// This computes the kinematic BL quantities that don't depend on whether
    /// the flow is laminar or turbulent: M², density ratio, shape factor H,
    /// kinematic shape factor Hk, and momentum Reynolds number Rθ.
    ///
    /// # Arguments
    /// * `params` - Global BL parameters
    pub fn blkin(&mut self, params: &FlowParameters) {
        let u = self.u;
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
        self.msq = u2_hstinv / denom;

        let tr = 1.0 + 0.5 * gm1 * self.msq;
        self.msq_u = 2.0 * self.msq * tr / u;
        self.msq_ms = u * u * tr / denom * hstinv_ms;

        // Edge static density (isentropic)
        // R = RST * TR^(-1/(γ-1))
        self.r = rst * tr.powf(-1.0 / gm1);
        self.r_u = -self.r / tr * 0.5 * self.msq_u;
        self.r_ms = -self.r / tr * 0.5 * self.msq_ms + rst_ms * tr.powf(-1.0 / gm1);

        // Shape factor H = δ*/θ
        self.h = d / t;
        self.h_d = 1.0 / t;
        self.h_t = -self.h / t;

        // Edge static/stagnation enthalpy ratio
        let herat = 1.0 - 0.5 * u * u * hstinv;
        let he_u = -u * hstinv;
        let he_ms = -0.5 * u * u * hstinv_ms;

        // Molecular viscosity ratio (Sutherland-type)
        // V = sqrt(HERAT^3) * (1+HVRAT)/(HERAT+HVRAT) / REYBL
        let herat_32 = (herat * herat * herat).sqrt(); // SQRT(HERAT**3): cube first, as the Fortran does
        self.v = herat_32 * (1.0 + hvrat) / (herat + hvrat) / reybl;
        let v_he = self.v * (1.5 / herat - 1.0 / (herat + hvrat));

        self.v_u = v_he * he_u;
        self.v_ms = -self.v / reybl * reybl_ms + v_he * he_ms;
        self.v_re = -self.v / reybl * reybl_re;

        // Kinematic shape factor Hk (compressibility correction)
        let (hk, hk_h, hk_msq) = hkin(self.h, self.msq);
        self.hk = hk;

        self.hk_u = hk_msq * self.msq_u;
        self.hk_t = hk_h * self.h_t;
        self.hk_d = hk_h * self.h_d;
        self.hk_ms = hk_msq * self.msq_ms;

        // Momentum thickness Reynolds number Rθ = ρ·U·θ/μ = R·U·T/V
        self.rt = self.r * u * t / self.v;
        self.rt_u = self.rt * (1.0 / u + self.r_u / self.r - self.v_u / self.v);
        self.rt_t = self.rt / t;
        self.rt_ms = self.rt * (self.r_ms / self.r - self.v_ms / self.v);
        self.rt_re = self.rt * (-self.v_re / self.v);
    }

    /// Copy all values from another station state
    ///
    /// Used for the COM1 = COM2 operation in XFOIL
    pub fn copy_from(&mut self, other: &BLStationState) {
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
    pub fn blvar(&mut self, flow_type: FlowRegime, _params: &FlowParameters) {
        let hk = self.hk;
        let rt = self.rt;
        let msq = self.msq;
        let h = self.h;
        let t = self.theta;
        let d = self.dstar;
        let s = self.ctau; // Ctau (or amplification for laminar)

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
        let (hc, hc_hk, hc_msq) = hc_turb(hk, msq);

        self.hc = hc;
        self.hc_u = hc_hk * self.hk_u + hc_msq * self.msq_u;
        self.hc_t = hc_hk * self.hk_t;
        self.hc_d = hc_hk * self.hk_d;
        self.hc_ms = hc_hk * self.hk_ms + hc_msq * self.msq_ms;

        // ====================================================================
        // Energy shape factor H* (from HSL or HST)
        // ====================================================================
        let hs_result = match flow_type {
            FlowRegime::Laminar => hs_lam(hk, rt, msq),
            FlowRegime::Turbulent | FlowRegime::Wake => hs_turb(hk, rt, msq),
        };

        self.hs = hs_result.val;
        self.hs_u = hs_result.val_hk * self.hk_u + hs_result.val_rt * self.rt_u + hs_result.val_msq * self.msq_u;
        self.hs_t = hs_result.val_hk * self.hk_t + hs_result.val_rt * self.rt_t;
        self.hs_d = hs_result.val_hk * self.hk_d;
        self.hs_ms = hs_result.val_hk * self.hk_ms + hs_result.val_rt * self.rt_ms + hs_result.val_msq * self.msq_ms;
        self.hs_re = hs_result.val_rt * self.rt_re;

        // ---- normalized slip velocity  Us
        let us = 0.5 * self.hs * (1.0 - (hk - 1.0) / (GBETA_LOCUS_B * h));
        let us_hs = 0.5 * (1.0 - (hk - 1.0) / (GBETA_LOCUS_B * h));
        let us_hk = 0.5 * self.hs * (-1.0 / (GBETA_LOCUS_B * h));
        let us_h = 0.5 * self.hs * (hk - 1.0) / (GBETA_LOCUS_B * (h * h));
        self.us = us;
        self.us_u = us_hs * self.hs_u + us_hk * self.hk_u;
        self.us_t = us_hs * self.hs_t + us_hk * self.hk_t + us_h * self.h_t;
        self.us_d = us_hs * self.hs_d + us_hk * self.hk_d + us_h * self.h_d;
        self.us_ms = us_hs * self.hs_ms + us_hk * self.hk_ms;
        self.us_re = us_hs * self.hs_re;
        if flow_type != FlowRegime::Wake && self.us > 0.95 {
            self.us = 0.98;
            self.us_u = 0.0;
            self.us_t = 0.0;
            self.us_d = 0.0;
            self.us_ms = 0.0;
            self.us_re = 0.0;
        }
        if flow_type == FlowRegime::Wake && self.us > 0.99995 {
            self.us = 0.99995;
            self.us_u = 0.0;
            self.us_t = 0.0;
            self.us_d = 0.0;
            self.us_ms = 0.0;
            self.us_re = 0.0;
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
        let cq = (SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk))).sqrt();
        let cq_hs = SQRTCTAUEQ_COEFFICIENT * hkb * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq;
        let cq_us = SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / usb * 0.5 / cq;
        let cq_hk = SQRTCTAUEQ_COEFFICIENT * self.hs * (hkc * hkc) / (usb * h * (hk * hk)) * 0.5 / cq
            - SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * (hkc * hkc) / (usb * h * ((hk * hk) * hk)) * 2.0 * 0.5 / cq
            + SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_hk;
        let cq_rt = SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * hkc / (usb * h * (hk * hk)) * 2.0 * 0.5 / cq * hkc_rt;
        let cq_h = -(SQRTCTAUEQ_COEFFICIENT * self.hs * hkb * (hkc * hkc) / (usb * h * (hk * hk)) / h * 0.5 / cq);
        self.cq = cq;
        self.cq_u = cq_hs * self.hs_u + cq_us * self.us_u + cq_hk * self.hk_u;
        self.cq_t = cq_hs * self.hs_t + cq_us * self.us_t + cq_hk * self.hk_t;
        self.cq_d = cq_hs * self.hs_d + cq_us * self.us_d + cq_hk * self.hk_d;
        self.cq_ms = cq_hs * self.hs_ms + cq_us * self.us_ms + cq_hk * self.hk_ms;
        self.cq_re = cq_hs * self.hs_re + cq_us * self.us_re;
        self.cq_u = self.cq_u + cq_rt * self.rt_u;
        self.cq_t = self.cq_t + cq_h * self.h_t + cq_rt * self.rt_t;
        self.cq_d = self.cq_d + cq_h * self.h_d;
        self.cq_ms = self.cq_ms + cq_rt * self.rt_ms;
        self.cq_re = self.cq_re + cq_rt * self.rt_re;

        // ---- set skin friction coefficient
        let (cf, cf_hk, cf_rt, cf_m) = match flow_type {
            // wake
            FlowRegime::Wake => (0.0, 0.0, 0.0, 0.0),
            // laminar
            FlowRegime::Laminar => {
                let r = cf_lam(hk, rt, msq);
                (r.val, r.val_hk, r.val_rt, r.val_msq)
            }
            // turbulent
            FlowRegime::Turbulent => {
                let r = cf_turb(hk, rt, msq, CF_TURBULENT_FACTOR);
                let l = cf_lam(hk, rt, msq);
                if l.val > r.val {
                    // laminar Cf is greater than turbulent Cf -- use laminar
                    // (this will only occur for unreasonably small Rtheta)
                    (l.val, l.val_hk, l.val_rt, l.val_msq)
                } else {
                    (r.val, r.val_hk, r.val_rt, r.val_msq)
                }
            }
        };
        self.cf = cf;
        self.cf_u = cf_hk * self.hk_u + cf_rt * self.rt_u + cf_m * self.msq_u;
        self.cf_t = cf_hk * self.hk_t + cf_rt * self.rt_t;
        self.cf_d = cf_hk * self.hk_d;
        self.cf_ms = cf_hk * self.hk_ms + cf_rt * self.rt_ms + cf_m * self.msq_ms;
        self.cf_re = cf_rt * self.rt_re;

        // ---- dissipation function    2 CD / H*
        match flow_type {
            FlowRegime::Laminar => {
                // laminar
                let r = di_lam(hk, rt);
                self.di = r.val;
                self.di_u = r.val_hk * self.hk_u + r.val_rt * self.rt_u;
                self.di_t = r.val_hk * self.hk_t + r.val_rt * self.rt_t;
                self.di_d = r.val_hk * self.hk_d;
                self.di_s = 0.0;
                self.di_ms = r.val_hk * self.hk_ms + r.val_rt * self.rt_ms;
                self.di_re = r.val_rt * self.rt_re;
            }
            FlowRegime::Turbulent => {
                // turbulent wall contribution
                let c = cf_turb(hk, rt, msq, CF_TURBULENT_FACTOR);
                let cf2t = c.val;
                let cf2t_u = c.val_hk * self.hk_u + c.val_rt * self.rt_u + c.val_msq * self.msq_u;
                let cf2t_t = c.val_hk * self.hk_t + c.val_rt * self.rt_t;
                let cf2t_d = c.val_hk * self.hk_d;
                let cf2t_ms = c.val_hk * self.hk_ms + c.val_rt * self.rt_ms + c.val_msq * self.msq_ms;
                let cf2t_re = c.val_rt * self.rt_re;
                let mut di = (0.5 * cf2t * us) * 2.0 / self.hs;
                let di_hs = -((0.5 * cf2t * us) * 2.0 / (self.hs * self.hs));
                let di_us = (0.5 * cf2t) * 2.0 / self.hs;
                let di_cf2t = (0.5 * us) * 2.0 / self.hs;
                let mut di_s = 0.0;
                let mut di_u = di_hs * self.hs_u + di_us * self.us_u + di_cf2t * cf2t_u;
                let mut di_t = di_hs * self.hs_t + di_us * self.us_t + di_cf2t * cf2t_t;
                let mut di_d = di_hs * self.hs_d + di_us * self.us_d + di_cf2t * cf2t_d;
                let mut di_ms = di_hs * self.hs_ms + di_us * self.us_ms + di_cf2t * cf2t_ms;
                let mut di_re = di_hs * self.hs_re + di_us * self.us_re + di_cf2t * cf2t_re;

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
                di_u = di_u * dfac + di * (df_hk * self.hk_u + df_rt * self.rt_u);
                di_t = di_t * dfac + di * (df_hk * self.hk_t + df_rt * self.rt_t);
                di_d = di_d * dfac + di * (df_hk * self.hk_d);
                di_ms = di_ms * dfac + di * (df_hk * self.hk_ms + df_rt * self.rt_ms);
                di_re = di_re * dfac + di * (df_rt * self.rt_re);
                di *= dfac;

                self.di = di;
                self.di_s = di_s;
                self.di_u = di_u;
                self.di_t = di_t;
                self.di_d = di_d;
                self.di_ms = di_ms;
                self.di_re = di_re;
            }
            FlowRegime::Wake => {
                // zero wall contribution for wake
                self.di = 0.0;
                self.di_s = 0.0;
                self.di_u = 0.0;
                self.di_t = 0.0;
                self.di_d = 0.0;
                self.di_ms = 0.0;
                self.di_re = 0.0;
            }
        }

        // ---- Add on turbulent outer layer contribution
        if flow_type != FlowRegime::Laminar {
            let dd = (s * s) * (0.995 - us) * 2.0 / self.hs;
            let dd_hs = -((s * s) * (0.995 - us) * 2.0 / (self.hs * self.hs));
            let dd_us = -((s * s) * 2.0 / self.hs);
            let dd_s = s * 2.0 * (0.995 - us) * 2.0 / self.hs;
            self.di = self.di + dd;
            self.di_s = dd_s;
            self.di_u = self.di_u + dd_hs * self.hs_u + dd_us * self.us_u;
            self.di_t = self.di_t + dd_hs * self.hs_t + dd_us * self.us_t;
            self.di_d = self.di_d + dd_hs * self.hs_d + dd_us * self.us_d;
            self.di_ms = self.di_ms + dd_hs * self.hs_ms + dd_us * self.us_ms;
            self.di_re = self.di_re + dd_hs * self.hs_re + dd_us * self.us_re;

            // add laminar stress contribution to outer layer CD
            let dd = 0.15 * ((0.995 - us) * (0.995 - us)) / rt * 2.0 / self.hs;
            let dd_us = -0.15 * (0.995 - us) * 2.0 / rt * 2.0 / self.hs;
            let dd_hs = -dd / self.hs;
            let dd_rt = -dd / rt;
            self.di = self.di + dd;
            self.di_u = self.di_u + dd_hs * self.hs_u + dd_us * self.us_u + dd_rt * self.rt_u;
            self.di_t = self.di_t + dd_hs * self.hs_t + dd_us * self.us_t + dd_rt * self.rt_t;
            self.di_d = self.di_d + dd_hs * self.hs_d + dd_us * self.us_d;
            self.di_ms = self.di_ms + dd_hs * self.hs_ms + dd_us * self.us_ms + dd_rt * self.rt_ms;
            self.di_re = self.di_re + dd_hs * self.hs_re + dd_us * self.us_re + dd_rt * self.rt_re;
        }

        if flow_type == FlowRegime::Turbulent {
            let l = di_lam(hk, rt);
            if l.val > self.di {
                // laminar CD is greater than turbulent CD -- use laminar
                // (this will only occur for unreasonably small Rtheta)
                self.di = l.val;
                self.di_s = 0.0;
                self.di_u = l.val_hk * self.hk_u + l.val_rt * self.rt_u;
                self.di_t = l.val_hk * self.hk_t + l.val_rt * self.rt_t;
                self.di_d = l.val_hk * self.hk_d;
                self.di_ms = l.val_hk * self.hk_ms + l.val_rt * self.rt_ms;
                self.di_re = l.val_rt * self.rt_re;
            }
        }

        if flow_type == FlowRegime::Wake {
            // laminar wake CD
            let l = dilw(hk, rt);
            if l.val > self.di {
                // laminar wake CD is greater than turbulent CD -- use laminar
                self.di = l.val;
                self.di_s = 0.0;
                self.di_u = l.val_hk * self.hk_u + l.val_rt * self.rt_u;
                self.di_t = l.val_hk * self.hk_t + l.val_rt * self.rt_t;
                self.di_d = l.val_hk * self.hk_d;
                self.di_ms = l.val_hk * self.hk_ms + l.val_rt * self.rt_ms;
                self.di_re = l.val_rt * self.rt_re;
            }
        }

        if flow_type == FlowRegime::Wake {
            // double dissipation for the wake (two wake halves)
            self.di *= 2.0;
            self.di_s *= 2.0;
            self.di_u *= 2.0;
            self.di_t *= 2.0;
            self.di_d *= 2.0;
            self.di_ms *= 2.0;
            self.di_re *= 2.0;
        }

        // ====================================================================
        // BL thickness Delta (Green's correlation)
        // DE = (3.15 + 1.72/(HK-1)) * T + D
        // ====================================================================
        let de = (3.15 + 1.72 / (hk - 1.0)) * t + d;
        let de_hk = -1.72 / ((hk - 1.0) * (hk - 1.0)) * t;

        self.de = de;
        self.de_u = de_hk * self.hk_u;
        self.de_t = de_hk * self.hk_t + 3.15 + 1.72 / (hk - 1.0);
        self.de_d = de_hk * self.hk_d + 1.0;
        self.de_ms = de_hk * self.hk_ms;

        // Clamp DE to reasonable values
        let hdmax = 12.0;
        if self.de > hdmax * t {
            self.de = hdmax * t;
            self.de_u = 0.0;
            self.de_t = hdmax;
            self.de_d = 0.0;
            self.de_ms = 0.0;
        }
    }
}

// ============================================================================
// Midpoint Skin Friction (BLMID)
// ============================================================================

/// Midpoint skin friction result (from XFOIL's BLMID)
///
/// This holds the skin friction coefficient at the midpoint between two
/// stations, along with its derivatives.
#[derive(Debug, Clone, Default)]
pub struct MidpointCf {
    /// Midpoint skin friction coefficient
    pub cfm: f64,

    /// Derivatives w.r.t. station 1
    pub cfm_u1: f64,
    pub cfm_t1: f64,
    pub cfm_d1: f64,

    /// Derivatives w.r.t. station 2
    pub cfm_u2: f64,
    pub cfm_t2: f64,
    pub cfm_d2: f64,

    /// Derivatives w.r.t. global parameters
    pub cfm_ms: f64,
    pub cfm_re: f64,
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
    pub fn compute(s1: &BLStationState, s2: &BLStationState, flow_type: FlowRegime, is_similarity: bool) -> Self {
        let mut result = Self::default();

        // For similarity station, station 1 equals station 2
        let (hk1, rt1, m1) = if is_similarity {
            (s2.hk, s2.rt, s2.msq)
        } else {
            (s1.hk, s1.rt, s1.msq)
        };

        let (hk1_u1, hk1_t1, hk1_d1, hk1_ms) = if is_similarity {
            (s2.hk_u, s2.hk_t, s2.hk_d, s2.hk_ms)
        } else {
            (s1.hk_u, s1.hk_t, s1.hk_d, s1.hk_ms)
        };

        let (rt1_u1, rt1_t1, rt1_ms, rt1_re) = if is_similarity {
            (s2.rt_u, s2.rt_t, s2.rt_ms, s2.rt_re)
        } else {
            (s1.rt_u, s1.rt_t, s1.rt_ms, s1.rt_re)
        };

        let (m1_u1, m1_ms) = if is_similarity {
            (s2.msq_u, s2.msq_ms)
        } else {
            (s1.msq_u, s1.msq_ms)
        };

        // Midpoint averages
        let hka = 0.5 * (hk1 + s2.hk);
        let rta = 0.5 * (rt1 + s2.rt);
        let ma = 0.5 * (m1 + s2.msq);

        // Midpoint skin friction coefficient
        let (cfm, cfm_hka, cfm_rta, cfm_ma) = match flow_type {
            FlowRegime::Wake => {
                // Zero skin friction in wake
                (0.0, 0.0, 0.0, 0.0)
            }
            FlowRegime::Laminar => {
                // Laminar Cf
                let cf_result = cf_lam(hka, rta, ma);
                (cf_result.val, cf_result.val_hk, cf_result.val_rt, 0.0)
            }
            FlowRegime::Turbulent => {
                // Turbulent Cf
                let cf_turb_result = cf_turb(hka, rta, ma, CF_TURBULENT_FACTOR);
                // Check if laminar is higher
                let cf_lam_result = cf_lam(hka, rta, ma);

                if cf_lam_result.val > cf_turb_result.val {
                    (cf_lam_result.val, cf_lam_result.val_hk, cf_lam_result.val_rt, 0.0)
                } else {
                    (
                        cf_turb_result.val,
                        cf_turb_result.val_hk,
                        cf_turb_result.val_rt,
                        cf_turb_result.val_msq,
                    )
                }
            }
        };

        result.cfm = cfm;

        // Derivatives w.r.t. station 1 variables (factor of 0.5 from averaging)
        result.cfm_u1 = 0.5 * (cfm_hka * hk1_u1 + cfm_ma * m1_u1 + cfm_rta * rt1_u1);
        result.cfm_t1 = 0.5 * (cfm_hka * hk1_t1 + cfm_rta * rt1_t1);
        result.cfm_d1 = 0.5 * cfm_hka * hk1_d1;

        // Derivatives w.r.t. station 2 variables
        result.cfm_u2 = 0.5 * (cfm_hka * s2.hk_u + cfm_ma * s2.msq_u + cfm_rta * s2.rt_u);
        result.cfm_t2 = 0.5 * (cfm_hka * s2.hk_t + cfm_rta * s2.rt_t);
        result.cfm_d2 = 0.5 * cfm_hka * s2.hk_d;

        // Derivatives w.r.t. global parameters
        result.cfm_ms = 0.5
            * (cfm_hka * hk1_ms
                + cfm_ma * m1_ms
                + cfm_rta * rt1_ms
                + cfm_hka * s2.hk_ms
                + cfm_ma * s2.msq_ms
                + cfm_rta * s2.rt_ms);
        result.cfm_re = 0.5 * (cfm_rta * rt1_re + cfm_rta * s2.rt_re);

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
pub fn dslim(dstr: &mut f64, thet: f64, _uedg: f64, msq: f64, hklim: f64) {
    let h = *dstr / thet;
    let (hk, hk_h, _hk_m) = hkin(h, msq);

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
        let state = BLStationState::default();
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
        let mut state = BLStationState::default();

        // Inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // Check primary variables
        assert_eq!(state.x, 0.1);
        assert_eq!(state.ampl, 3.0);
        assert_eq!(state.ctau, 0.015);
        assert_eq!(state.theta, 0.002);
        assert_eq!(state.dstar, 0.005);
        assert_eq!(state.dw, 0.0);

        // At M=0, U2 = Uei (no transformation)
        assert_relative_eq!(state.u, 1.2, epsilon = 1e-10);
        assert_relative_eq!(state.u_uei, 1.0, epsilon = 1e-10);
        // Note: u_ms is the sensitivity d(U2)/d(M²). Even at M=0 this is non-zero because
        // COMSET's TKL_MSQ = 1/(1+BETA)² = 0.25 at M=0 (BETA = 1).
        // U2_MS = (U2*UEI² - UEI) * TKBL_MS = (1.2*1.44 - 1.2) * 0.25 = 0.132
        assert_relative_eq!(state.u_ms, 0.132, epsilon = 1e-6);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (derived with TKBL = 1/beta - 1, not COMSET's TKLAM) — regenerate from the M=0.3 coverage case"]
    fn test_blprv_compressible() {
        // Test case: M=0.5, Re=1e6
        // XFOIL reference values from Fortran test
        let params = FlowParameters::new(0.5, 1e6, 1.4);
        let mut state = BLStationState::default();

        // Same inputs
        let xsi = 0.1;
        let ami = 3.0;
        let cti = 0.015;
        let thi = 0.002;
        let dsi = 0.005;
        let dswaki = 0.0;
        let uei = 1.2;

        state.blprv(xsi, ami, cti, thi, dsi, dswaki, uei, &params);

        // XFOIL reference: U2 = 1.305093527
        assert_relative_eq!(state.u, 1.305093527, epsilon = 1e-5);
        // XFOIL reference: U2_UEI = 1.711017609
        assert_relative_eq!(state.u_uei, 1.711017609, epsilon = 1e-5);
        // XFOIL reference: U2_MS = 0.6728398800
        assert_relative_eq!(state.u_ms, 0.6728398800, epsilon = 1e-5);
    }

    // ========================================================================
    // BLKIN Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blkin_incompressible() {
        // Test case: M=0, Re=1e6
        let params = FlowParameters::incompressible(1e6);
        let mut state = BLStationState::default();

        // Run BLPRV first
        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.blkin(&params);

        // XFOIL reference values (M=0 case)
        // M2 = 0
        assert_relative_eq!(state.msq, 0.0, epsilon = 1e-10);

        // R2 = 1.0
        assert_relative_eq!(state.r, 1.0, epsilon = 1e-10);

        // H2 = D2/T2 = 0.005/0.002 = 2.5
        assert_relative_eq!(state.h, 2.5, epsilon = 1e-10);
        assert_relative_eq!(state.h_t, -1250.0, epsilon = 1e-6); // -H/T
        assert_relative_eq!(state.h_d, 500.0, epsilon = 1e-6); // 1/T

        // HK2 = 2.5 (same as H at M=0)
        assert_relative_eq!(state.hk, 2.5, epsilon = 1e-6);

        // RT2 = 2400 (R*U*T/V = 1*1.2*0.002/1e-6)
        assert_relative_eq!(state.rt, 2400.0, epsilon = 1.0);

        // RT2_U2 = 2000 (RT/U = 2400/1.2)
        assert_relative_eq!(state.rt_u, 2000.0, epsilon = 1.0);

        // RT2_T2 = 1200000 (RT/T = 2400/0.002)
        assert_relative_eq!(state.rt_t, 1200000.0, epsilon = 100.0);
    }

    #[test]
    #[ignore = "S10: expected values were unsourced (assumed HVRAT=0.35; XFOIL's analysis path leaves HVRAT=0) — regenerate from the M=0.3 coverage case"]
    fn test_blkin_compressible() {
        // Test case: M=0.5, Re=1e6
        let params = FlowParameters::new(0.5, 1e6, 1.4);
        let mut state = BLStationState::default();

        // Run BLPRV first
        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);

        // Run BLKIN
        state.blkin(&params);

        // XFOIL reference values (M=0.5 case)
        // M2 = 0.4413362443
        assert_relative_eq!(state.msq, 0.4413362443, epsilon = 1e-5);

        // R2 = 0.9143959880
        assert_relative_eq!(state.r, 0.9143959880, epsilon = 1e-5);

        // HK2 = 2.259336710
        assert_relative_eq!(state.hk, 2.259336710, epsilon = 1e-5);

        // RT2 = 2453.646240
        assert_relative_eq!(state.rt, 2453.646240, epsilon = 1.0);
    }

    #[test]
    fn test_blkin_shape_factor_derivatives() {
        // Verify shape factor derivatives are computed correctly
        let params = FlowParameters::incompressible(1e6);
        let mut state = BLStationState::default();

        state.blprv(0.1, 3.0, 0.015, 0.002, 0.005, 0.0, 1.2, &params);
        state.blkin(&params);

        // H = D/T = 2.5
        // dH/dD = 1/T = 1/0.002 = 500
        // dH/dT = -D/T² = -0.005/0.002² = -1250
        assert_relative_eq!(state.h_d, 1.0 / 0.002, epsilon = 1e-10);
        assert_relative_eq!(state.h_t, -0.005 / 0.002_f64.powi(2), epsilon = 1e-10);
    }

    // ========================================================================
    // BLVAR Tests - validate against XFOIL Fortran output
    // ========================================================================

    #[test]
    fn test_blvar_laminar() {
        // Test case: Laminar BL
        // XFOIL reference values from Fortran test
        let params = FlowParameters::incompressible(1e6);
        let mut state = BLStationState::default();

        // Set up: HK=2.5, RT=2400, M=0, H=2.5, T=0.002, D=0.005, S=3.0
        // We need to set the kinematic variables directly since blprv/blkin
        // uses its own computations
        state.hk = 2.5;
        state.rt = 2400.0;
        state.msq = 0.0;
        state.h = 2.5;
        state.theta = 0.002;
        state.dstar = 0.005;
        state.ctau = 3.0; // Amplification for laminar

        // Set the derivatives (simplified - zero for test)
        state.hk_u = 0.0;
        state.hk_t = -1250.0;
        state.hk_d = 500.0;
        state.hk_ms = 0.0;
        state.rt_u = 2000.0;
        state.rt_t = 1200000.0;
        state.rt_ms = 0.0;
        state.rt_re = 0.0024;
        state.h_t = -1250.0;
        state.h_d = 500.0;
        state.msq_u = 0.0;
        state.msq_ms = 0.0;

        state.blvar(FlowRegime::Laminar, &params);

        // XFOIL reference values:
        // HC2 = 0.0 (M=0)
        assert_relative_eq!(state.hc, 0.0, epsilon = 1e-10);

        // HS2 = 1.584867239
        assert_relative_eq!(state.hs, 1.584867239, epsilon = 1e-5);

        // US2 = 0.1584867090
        assert_relative_eq!(state.us, 0.1584867090, epsilon = 1e-5);

        // CQ2 = 0.07772709429
        assert_relative_eq!(state.cq, 0.07772709429, epsilon = 1e-4);

        // CF2 = 0.2045119036e-3
        assert_relative_eq!(state.cf, 0.2045119036e-3, epsilon = 1e-7);

        // DI2 = 0.9419409616e-4
        assert_relative_eq!(state.di, 0.9419409616e-4, epsilon = 1e-7);

        // DE2 = 0.01359333377
        assert_relative_eq!(state.de, 0.01359333377, epsilon = 1e-6);
    }

    #[test]
    fn test_blvar_turbulent() {
        // Test case: Turbulent BL
        // XFOIL reference values from Fortran test
        let params = FlowParameters::incompressible(1e6);
        let mut state = BLStationState::default();

        // Set up: HK=1.4, RT=10000, M=0, H=1.4, T=0.005, D=0.007, S=0.015
        state.hk = 1.4;
        state.rt = 10000.0;
        state.msq = 0.0;
        state.h = 1.4;
        state.theta = 0.005;
        state.dstar = 0.007;
        state.ctau = 0.015;

        // Set derivatives (simplified)
        state.hk_u = 0.0;
        state.hk_t = -280.0;
        state.hk_d = 200.0;
        state.hk_ms = 0.0;
        state.rt_u = 0.0;
        state.rt_t = 2000000.0;
        state.rt_ms = 0.0;
        state.rt_re = 0.01;
        state.h_t = -280.0;
        state.h_d = 200.0;
        state.msq_u = 0.0;
        state.msq_ms = 0.0;
        state.hs_u = 0.0;
        state.hs_t = 0.0;
        state.hs_d = 0.0;
        state.hs_ms = 0.0;
        state.hs_re = 0.0;
        state.us_u = 0.0;
        state.us_t = 0.0;
        state.us_d = 0.0;
        state.us_ms = 0.0;
        state.us_re = 0.0;

        state.blvar(FlowRegime::Turbulent, &params);

        // XFOIL reference values:
        // HS2 = 1.755310297
        assert_relative_eq!(state.hs, 1.755310297, epsilon = 1e-4);

        // US2 = 0.5433103442
        assert_relative_eq!(state.us, 0.5433103442, epsilon = 1e-4);

        // CQ2 = 0.03632329032
        assert_relative_eq!(state.cq, 0.03632329032, epsilon = 1e-4);

        // CF2 = 0.2286923816e-2
        assert_relative_eq!(state.cf, 0.2286923816e-2, epsilon = 1e-5);

        // DI2 = 0.8065673755e-3 (turbulent DI is more complex, allow larger tolerance)
        assert_relative_eq!(state.di, 0.8065673755e-3, epsilon = 1e-4);

        // DE2 = 0.04425000027
        assert_relative_eq!(state.de, 0.04425000027, epsilon = 1e-5);
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
        let mut s1 = BLStationState::default();
        s1.hk = 2.3;
        s1.rt = 2000.0;
        s1.msq = 0.0;
        s1.hk_u = 0.0;
        s1.hk_t = -1150.0;
        s1.hk_d = 500.0;
        s1.hk_ms = -0.29;
        s1.rt_u = 1666.67;
        s1.rt_t = 1000000.0;
        s1.rt_ms = 0.0;
        s1.rt_re = 0.002;
        s1.msq_u = 0.0;
        s1.msq_ms = 0.0;

        // Station 2
        let mut s2 = BLStationState::default();
        s2.hk = 2.5;
        s2.rt = 2400.0;
        s2.msq = 0.0;
        s2.hk_u = 0.0;
        s2.hk_t = -1250.0;
        s2.hk_d = 500.0;
        s2.hk_ms = -0.29;
        s2.rt_u = 2000.0;
        s2.rt_t = 1200000.0;
        s2.rt_ms = 0.0;
        s2.rt_re = 0.0024;
        s2.msq_u = 0.0;
        s2.msq_ms = 0.0;

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Laminar, false);

        // XFOIL reference values
        assert_relative_eq!(result.cfm, 0.2577280102e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_u1, -0.9762444097e-4, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_t1, 0.1515112668, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_d1, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cfm_u2, -0.1171490949e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_t2, 0.1580646932, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_d2, -0.9134165943e-1, epsilon = 1e-4);
        assert_relative_eq!(result.cfm_ms, 0.1059563219e-3, epsilon = 1e-7);
        assert_relative_eq!(result.cfm_re, -0.2577280056e-9, epsilon = 1e-12);
    }

    #[test]
    fn test_blmid_turbulent() {
        // Test case: Turbulent (ITYP=2)
        // XFOIL reference values from Fortran test
        // HKA = 1.375, RTA = 9000, MA = 0

        // Station 1
        let mut s1 = BLStationState::default();
        s1.hk = 1.35;
        s1.rt = 8000.0;
        s1.msq = 0.0;
        s1.hk_u = 0.0;
        s1.hk_t = -270.0;
        s1.hk_d = 200.0;
        s1.hk_ms = -0.29;
        s1.rt_u = 0.0;
        s1.rt_t = 1600000.0;
        s1.rt_ms = 0.0;
        s1.rt_re = 0.008;
        s1.msq_u = 0.0;
        s1.msq_ms = 0.0;

        // Station 2
        let mut s2 = BLStationState::default();
        s2.hk = 1.40;
        s2.rt = 10000.0;
        s2.msq = 0.0;
        s2.hk_u = 0.0;
        s2.hk_t = -280.0;
        s2.hk_d = 200.0;
        s2.hk_ms = -0.29;
        s2.rt_u = 0.0;
        s2.rt_t = 2000000.0;
        s2.rt_ms = 0.0;
        s2.rt_re = 0.01;
        s2.msq_u = 0.0;
        s2.msq_ms = 0.0;

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Turbulent, false);

        // XFOIL reference values (turbulent Cf used)
        assert_relative_eq!(result.cfm, 0.2450317144e-2, epsilon = 1e-5);
        assert_relative_eq!(result.cfm_t1, 0.5299984217, epsilon = 1e-2);
        assert_relative_eq!(result.cfm_d1, -0.4310033321, epsilon = 1e-3);
        assert_relative_eq!(result.cfm_t2, 0.5385844707, epsilon = 1e-2);
        assert_relative_eq!(result.cfm_d2, -0.4310033321, epsilon = 1e-3);
    }

    #[test]
    fn test_blmid_wake() {
        // Test case: Wake (ITYP=3)
        // Wake should have zero skin friction
        let s1 = BLStationState::default();
        let s2 = BLStationState::default();

        let result = MidpointCf::compute(&s1, &s2, FlowRegime::Wake, false);

        assert_eq!(result.cfm, 0.0);
        assert_eq!(result.cfm_u1, 0.0);
        assert_eq!(result.cfm_t1, 0.0);
        assert_eq!(result.cfm_d1, 0.0);
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
        dslim(&mut dstr, thet, uedg, msq, hklim);

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
        dslim(&mut dstr, thet, uedg, msq, hklim);

        // Should increase dstr to raise Hk to at least hklim
        let new_h = dstr / thet;
        assert!(new_h >= 1.5, "H should be raised to at least hklim");
    }
}
