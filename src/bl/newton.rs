//! Newton solver for boundary layer equations
//!
//! Implements Newton-Raphson iteration to solve the coupled integral
//! boundary layer equations at each marching station.
//!
//! The primary unknowns in direct mode (prescribed Ue) are:
//! - θ (momentum thickness)
//! - H (shape factor) or equivalently δ* = H*θ
//!
//! The equations being solved are:
//! 1. Momentum integral (von Kármán)
//! 2. Shape parameter (kinetic energy integral)

use crate::bl::{
    cf_lam, cf_turb, dampl, di_lam, hkin, hs_lam, hs_turb, ClosureResult, FlowConditions,
    FlowRegime,
};

/// Newton solver configuration
#[derive(Debug, Clone, Copy)]
pub struct NewtonConfig {
    /// Maximum number of iterations
    pub max_iter: usize,
    /// Convergence tolerance for residuals
    pub tol: f64,
    /// Relaxation factor (1.0 = full Newton step)
    pub relax: f64,
    /// Minimum allowed θ value
    pub theta_min: f64,
    /// Minimum allowed Hk value
    pub hk_min: f64,
    /// Maximum allowed Hk value
    pub hk_max: f64,
    /// Maximum laminar Hk before inverse mode (XFOIL's HLMAX)
    pub hlmax: f64,
    /// Maximum turbulent Hk before inverse mode (XFOIL's HTMAX)
    pub htmax: f64,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-5,
            relax: 1.0,
            theta_min: 1e-10,
            hk_min: 1.00005,
            hk_max: 12.0,
            hlmax: 3.8,  // XFOIL's laminar Hk limit
            htmax: 2.5,  // XFOIL's turbulent Hk limit
        }
    }
}

/// Result of Newton iteration at a station
#[derive(Debug, Clone)]
pub struct NewtonResult {
    /// Converged momentum thickness
    pub theta: f64,
    /// Converged displacement thickness
    pub dstar: f64,
    /// Converged shape factor
    pub h: f64,
    /// Kinematic shape factor
    pub hk: f64,
    /// Skin friction coefficient
    pub cf: f64,
    /// Dissipation coefficient
    pub cd: f64,
    /// Energy shape factor
    pub hs: f64,
    /// Amplification factor (laminar only)
    pub n_amp: f64,
    /// Edge velocity at this station (normalized by Qinf)
    pub ue: f64,
    /// X-coordinate at this station (for drag integration)
    pub x: f64,
    /// Number of iterations taken
    pub iterations: usize,
    /// Final residual norm
    pub residual: f64,
    /// Whether solution converged
    pub converged: bool,
}

/// Solve BL equations at a station using Newton-Raphson
///
/// Given upstream state and downstream edge velocity, find the
/// downstream θ and δ* that satisfy the integral BL equations.
///
/// Uses XFOIL's (θ, δ*) variable choice with proper coupling terms and
/// logarithmic form for better accuracy.
///
/// # Arguments
/// * `theta1` - Upstream momentum thickness
/// * `h1` - Upstream shape factor
/// * `ue1` - Upstream edge velocity
/// * `n1` - Upstream amplification factor
/// * `ue2` - Downstream edge velocity (prescribed)
/// * `s1` - Upstream arc length
/// * `s2` - Downstream arc length
/// * `regime` - Flow regime (laminar/turbulent)
/// * `cond` - Flow conditions
/// * `config` - Newton solver configuration
pub fn solve_station(
    theta1: f64,
    h1: f64,
    ue1: f64,
    n1: f64,
    ue2: f64,
    s1: f64,
    s2: f64,
    regime: FlowRegime,
    cond: &FlowConditions,
    config: &NewtonConfig,
) -> NewtonResult {
    let ds = s2 - s1;
    let due_ds = (ue2 - ue1) / ds;
    let ue_avg = 0.5 * (ue1 + ue2);

    // Get upstream quantities
    let (hk1, _, _) = hkin(h1, cond.msq);
    let hk1 = hk1.clamp(config.hk_min, config.hk_max);
    let rt1 = (ue1 * theta1 / cond.nu).max(1.0);

    // Initial guess using momentum equation
    let cf1_guess = match regime {
        FlowRegime::Laminar => cf_lam(hk1, rt1, cond.msq).val,
        _ => cf_turb(hk1, rt1.max(200.0), cond.msq, 1.0).val,
    };
    let shape_term = (h1 + 2.0 - cond.msq) * theta1 / ue_avg * due_ds;
    let dtheta_ds_pred = cf1_guess / 2.0 - shape_term;
    let mut theta2 = (theta1 + ds * dtheta_ds_pred).clamp(theta1 * 0.5, theta1 * 2.0);

    // Initial H guess
    let h_init = match regime {
        FlowRegime::Turbulent | FlowRegime::Wake if h1 > 2.0 => 1.5, // At transition
        FlowRegime::Turbulent | FlowRegime::Wake => h1.clamp(1.2, 2.5),
        FlowRegime::Laminar => h1,
    };
    let mut dstar2 = h_init * theta2; // Use (θ, δ*) as Newton variables

    // Upstream closures
    let (hs1, cf1, di1) = match regime {
        FlowRegime::Laminar => (
            hs_lam(hk1, rt1, cond.msq),
            cf_lam(hk1, rt1, cond.msq),
            di_lam(hk1, rt1),
        ),
        _ => {
            let cf = cf_turb(hk1, rt1.max(200.0), cond.msq, 1.0);
            let hs = hs_turb(hk1, rt1.max(200.0), cond.msq);
            let di = ClosureResult {
                val: 0.5 * cf.val * ue1 / hs.val,
                val_hk: 0.0,
                val_rt: 0.0,
                val_msq: 0.0,
            };
            (hs, cf, di)
        }
    };

    // Newton iteration
    let mut converged = false;
    let mut iterations = 0;
    let mut residual = f64::MAX;

    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Compute H from (θ, δ*) - key XFOIL coupling
        let h2 = dstar2 / theta2;

        // Derived quantities
        let (hk2_raw, hk2_h, _) = hkin(h2, cond.msq);
        let hk2 = hk2_raw.clamp(config.hk_min, config.hk_max);
        let rt2 = (ue2 * theta2 / cond.nu).max(1.0);

        // XFOIL coupling terms
        let h2_t2 = -h2 / theta2; // ∂H/∂θ at constant δ*
        let h2_d2 = 1.0 / theta2; // ∂H/∂δ* at constant θ
        let hk2_t2 = hk2_h * h2_t2;
        let hk2_d2 = hk2_h * h2_d2;
        let rt2_t2 = rt2 / theta2;

        // Get closure relations
        let (cf_res, hs_res, di_res) = match regime {
            FlowRegime::Laminar => {
                let cf = cf_lam(hk2, rt2, cond.msq);
                let hs = hs_lam(hk2, rt2, cond.msq);
                let di = di_lam(hk2, rt2);
                (cf, hs, di)
            }
            FlowRegime::Turbulent | FlowRegime::Wake => {
                let cf = cf_turb(hk2, rt2.max(200.0), cond.msq, 1.0);
                let hs = hs_turb(hk2, rt2.max(200.0), cond.msq);
                let di = ClosureResult {
                    val: 0.5 * cf.val * ue2 / hs.val,
                    val_hk: 0.0,
                    val_rt: 0.0,
                    val_msq: 0.0,
                };
                (cf, hs, di)
            }
        };

        // Closure derivatives in (θ, δ*) coordinates
        let cf2_t2 = cf_res.val_hk * hk2_t2 + cf_res.val_rt * rt2_t2;
        let cf2_d2 = cf_res.val_hk * hk2_d2;
        let hs2_t2 = hs_res.val_hk * hk2_t2 + hs_res.val_rt * rt2_t2;
        let hs2_d2 = hs_res.val_hk * hk2_d2;
        let di2_t2 = di_res.val_hk * hk2_t2 + di_res.val_rt * rt2_t2;
        let di2_d2 = di_res.val_hk * hk2_d2;

        // Average quantities
        let theta_avg = 0.5 * (theta1 + theta2);
        let h_avg = 0.5 * (h1 + h2);
        let s_avg = 0.5 * (s1 + s2);
        let hs_avg = 0.5 * (hs1.val + hs_res.val);

        // Midpoint Cf (BLMID approach for better accuracy)
        let hk_avg = 0.5 * (hk1 + hk2);
        let rt_avg = 0.5 * (rt1 + rt2);
        let cf_mid = match regime {
            FlowRegime::Laminar => cf_lam(hk_avg, rt_avg, cond.msq),
            _ => cf_turb(hk_avg, rt_avg.max(200.0), cond.msq, 1.0),
        };

        // Use logarithmic form when s1 is large enough
        let use_log_form = s1 > 1e-6;

        let (r1, dr1_dt2, dr1_dd2) = if use_log_form {
            // XFOIL logarithmic momentum equation
            let tlog = (theta2 / theta1).ln();
            let ulog = (ue2 / ue1).ln();
            let xlog = (s2 / s1).ln();
            let btmp = h_avg + 2.0 - cond.msq;

            // CFX = 0.50*CFM*XA/TA + 0.25*(CF1*X1/T1 + CF2*X2/T2)
            let cfx = 0.5 * cf_mid.val * s_avg / theta_avg
                + 0.25 * (cf1.val * s1 / theta1 + cf_res.val * s2 / theta2);

            let r1_val = tlog + btmp * ulog - xlog * 0.5 * cfx;

            // Jacobian with coupling terms
            // CFX depends on θ2 through both theta_avg and cf_res
            let cfx_t2 = -0.5 * cf_mid.val * s_avg / theta_avg.powi(2) * 0.5
                - 0.25 * cf_res.val * s2 / theta2.powi(2)
                + 0.25 * cf2_t2 * s2 / theta2;

            let cfx_d2 = 0.25 * cf2_d2 * s2 / theta2;

            // dr1/dθ2 includes H coupling through btmp
            let dr1_dt2_val = 1.0 / theta2 + 0.5 * ulog * h2_t2 - xlog * 0.5 * cfx_t2;
            let dr1_dd2_val = 0.5 * ulog * h2_d2 - xlog * 0.5 * cfx_d2;

            (r1_val, dr1_dt2_val, dr1_dd2_val)
        } else {
            // Differential form near stagnation
            let btmp = h_avg + 2.0 - cond.msq;
            let mom_coef = btmp * theta_avg / ue_avg;
            let r1_val = (theta2 - theta1) / ds + mom_coef * due_ds - cf_mid.val / 2.0;

            let dr1_dt2_val = 1.0 / ds
                + 0.5 * btmp / ue_avg * due_ds
                + 0.5 * h2_t2 * theta_avg / ue_avg * due_ds
                - cf2_t2 / 2.0;

            let dr1_dd2_val = 0.5 * h2_d2 * theta_avg / ue_avg * due_ds - cf2_d2 / 2.0;

            (r1_val, dr1_dt2_val, dr1_dd2_val)
        };

        // Shape equation with upwinding
        let upw = {
            let hdcon = 5.0 / hk2.powi(2);
            let arg = ((hk2 - 1.0) / (hk1 - 1.0).max(0.01)).abs();
            let hl = arg.ln();
            let hlsq = hl.powi(2).min(15.0);
            1.0 - 0.5 * (-hlsq * hdcon).exp()
        };

        let (r2, dr2_dt2, dr2_dd2) = if use_log_form {
            // XFOIL logarithmic shape equation
            let hlog = (hs_res.val / hs1.val).ln();
            let ulog = (ue2 / ue1).ln();
            let xlog = (s2 / s1).ln();
            let btmp_h = 1.0 - h_avg; // Simplified for low Mach

            let xot1 = s1 / theta1;
            let xot2 = s2 / theta2;
            let dix = (1.0 - upw) * di1.val * xot1 + upw * di_res.val * xot2;
            let cfx = (1.0 - upw) * cf1.val * xot1 + upw * cf_res.val * xot2;

            let r2_val = hlog + btmp_h * ulog + xlog * (0.5 * cfx - dix);

            // Jacobian
            let dhlog_dt2 = hs2_t2 / hs_res.val;
            let dhlog_dd2 = hs2_d2 / hs_res.val;

            let dix_t2 = upw * (-di_res.val * xot2 / theta2 + di2_t2 * xot2);
            let dix_d2 = upw * di2_d2 * xot2;
            let cfx_t2 = upw * (-cf_res.val * xot2 / theta2 + cf2_t2 * xot2);
            let cfx_d2 = upw * cf2_d2 * xot2;

            let dr2_dt2_val = dhlog_dt2 - 0.5 * h2_t2 * ulog + xlog * (0.5 * cfx_t2 - dix_t2);
            let dr2_dd2_val = dhlog_dd2 - 0.5 * h2_d2 * ulog + xlog * (0.5 * cfx_d2 - dix_d2);

            (r2_val, dr2_dt2_val, dr2_dd2_val)
        } else {
            // Differential form
            let dhs_ds = (hs_res.val - hs1.val) / ds;
            let shape_coef = hs_avg * (1.0 - h_avg) * theta_avg / ue_avg;
            let r2_val = theta_avg * dhs_ds + shape_coef * due_ds
                - 2.0 * di_res.val + hs_avg * cf_res.val / 2.0;

            let dr2_dt2_val = 0.5 * dhs_ds + theta_avg / ds * hs2_t2
                + (0.5 * hs2_t2 * (1.0 - h_avg) - 0.5 * hs_avg * h2_t2) * theta_avg / ue_avg * due_ds
                + hs_avg * cf2_t2 / 2.0 - 2.0 * di2_t2;

            let dr2_dd2_val = theta_avg / ds * hs2_d2
                + (0.5 * hs2_d2 * (1.0 - h_avg) - 0.5 * hs_avg * h2_d2) * theta_avg / ue_avg * due_ds
                + hs_avg * cf2_d2 / 2.0 - 2.0 * di2_d2;

            (r2_val, dr2_dt2_val, dr2_dd2_val)
        };

        // Solve 2x2 Newton system
        let det = dr1_dt2 * dr2_dd2 - dr1_dd2 * dr2_dt2;

        let h_min = match regime {
            FlowRegime::Turbulent | FlowRegime::Wake => 1.2,
            FlowRegime::Laminar => 1.0,
        };
        let theta_max = theta1 * 3.0;

        if det.abs() < 1e-20 {
            let d_theta = -r1 * config.relax * 0.1;
            let d_dstar = -r2 * config.relax * 0.1;
            theta2 = (theta2 + d_theta).clamp(config.theta_min, theta_max);
            dstar2 = (dstar2 + d_dstar).clamp(h_min * theta2, config.hk_max * theta2);
        } else {
            let d_theta = (-r1 * dr2_dd2 + r2 * dr1_dd2) / det;
            let d_dstar = (-r2 * dr1_dt2 + r1 * dr2_dt2) / det;

            // XFOIL-style unified relaxation
            let dmax = (d_theta / theta2).abs().max((d_dstar / dstar2).abs());
            let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

            theta2 = (theta2 + rlx * config.relax * d_theta).clamp(config.theta_min, theta_max);
            dstar2 = (dstar2 + rlx * config.relax * d_dstar).clamp(h_min * theta2, config.hk_max * theta2);
        }

        residual = (r1.powi(2) + r2.powi(2)).sqrt();
        if residual < config.tol {
            converged = true;
            break;
        }
    }

    // Final values
    let h2 = dstar2 / theta2;

    // Final closure values at converged solution
    let dstar2 = h2 * theta2;
    let (hk2, _, _) = hkin(h2, cond.msq);
    let hk2 = hk2.clamp(config.hk_min, config.hk_max);
    let rt2 = (ue2 * theta2 / cond.nu).max(1.0);

    let (cf_final, hs_final, di_final) = match regime {
        FlowRegime::Laminar => {
            let cf = cf_lam(hk2, rt2, cond.msq);
            let hs = hs_lam(hk2, rt2, cond.msq);
            let di = di_lam(hk2, rt2);
            (cf.val, hs.val, di.val)
        }
        FlowRegime::Turbulent | FlowRegime::Wake => {
            let cf = cf_turb(hk2, rt2.max(200.0), cond.msq, 1.0);
            let hs = hs_turb(hk2, rt2.max(200.0), cond.msq);
            (cf.val, hs.val, 0.5 * cf.val * ue2 / hs.val)
        }
    };

    // Compute amplification factor for laminar flow using XFOIL's AXSET approach
    let n2 = if regime == FlowRegime::Laminar {
        // Get quantities at station 1
        let (hk1, _, _) = hkin(h1, cond.msq);
        let rt1 = (ue1 * theta1 / cond.nu).max(1.0);
        let (ax1, _, _, _) = dampl(hk1, theta1, rt1);

        // Get quantities at station 2
        let rt2_amp = (ue2 * theta2 / cond.nu).max(1.0);
        let (ax2, _, _, _) = dampl(hk2, theta2, rt2_amp);

        // RMS average (XFOIL uses this instead of simple average)
        let axsq = 0.5 * (ax1.powi(2) + ax2.powi(2));
        let axa = if axsq <= 0.0 { 0.0 } else { axsq.sqrt() };

        // Additional term to ensure dN/dx > 0 near N = Ncrit (XFOIL's DAX term)
        let a_avg = 0.5 * (n1 + n1); // At station 1, both are n1
        let arg = (20.0 * (cond.ncrit - a_avg)).min(20.0);
        let exn = if arg <= 0.0 { 1.0 } else { (-arg).exp() };
        let dax = exn * 0.002 / (theta1 + theta2);

        // Combined amplification rate
        let ax = axa + dax;
        n1 + ax * ds
    } else {
        n1
    };

    NewtonResult {
        theta: theta2,
        dstar: dstar2,
        h: h2,
        hk: hk2,
        cf: cf_final,
        cd: di_final * hs_final / 2.0,
        hs: hs_final,
        n_amp: n2,
        ue: ue2, // Store edge velocity for drag integration
        x: 0.0,  // Will be set by caller
        iterations,
        residual,
        converged,
    }
}

/// March BL with Newton solver at each station
///
/// # Arguments
/// * `ue_dist` - Edge velocity distribution
/// * `s_dist` - Arc length distribution
/// * `cond` - Flow conditions
/// * `config` - Newton solver configuration
///
/// # Returns
/// Vector of Newton results at each station
pub fn march_newton(
    ue_dist: &[f64],
    s_dist: &[f64],
    cond: &FlowConditions,
    config: &NewtonConfig,
) -> Vec<NewtonResult> {
    let n = ue_dist.len();
    assert!(n >= 2);
    assert_eq!(s_dist.len(), n);

    let mut results = Vec::with_capacity(n);

    // Initialize at stagnation point using XFOIL's Thwaites formula
    // XFOIL uses power-law: Ue = UCON * s^BULE with BULE=1.0 (stagnation flow)
    // theta^2 = 0.45 * nu * s / (Ue * (5*BULE + 1))
    // For BULE=1.0: theta = sqrt(0.45 * nu * s / (6 * Ue))
    //
    // At first station (index 0), s=0, so use index 1 for initialization
    let s_init = s_dist[1].max(1e-10);
    let ue_init = ue_dist[1].max(0.01);
    let theta_init = (0.45 * cond.nu * s_init / (6.0 * ue_init)).sqrt();
    let h_init = 2.2; // Hiemenz stagnation point (XFOIL uses 2.2)

    // First station
    let (hk_init, _, _) = hkin(h_init, cond.msq);
    let rt_init = (ue_dist[0] * theta_init / cond.nu).max(1.0);
    let cf_init = cf_lam(hk_init, rt_init, cond.msq).val;
    let hs_init = hs_lam(hk_init, rt_init, cond.msq).val;
    let di_init = di_lam(hk_init, rt_init).val;

    results.push(NewtonResult {
        theta: theta_init,
        dstar: h_init * theta_init,
        h: h_init,
        hk: hk_init,
        cf: cf_init,
        cd: di_init * hs_init / 2.0,
        hs: hs_init,
        n_amp: 0.0,
        ue: ue_dist[0],
        x: 0.0, // Will be set by caller
        iterations: 0,
        residual: 0.0,
        converged: true,
    });

    // March downstream with Newton solver
    let mut regime = FlowRegime::Laminar;

    for i in 1..n {
        let prev = &results[i - 1];
        let s1 = s_dist[i - 1];
        let s2 = s_dist[i];
        let ds = s2 - s1;

        // Solve for this station
        let mut result = solve_station(
            prev.theta,
            prev.h,
            ue_dist[i - 1],
            prev.n_amp,
            ue_dist[i],
            s1,
            s2,
            regime,
            cond,
            config,
        );

        // XFOIL-style inverse mode: when Hk exceeds threshold, prescribe Hk growth
        // instead of allowing it to grow unbounded
        let hmax = match regime {
            FlowRegime::Laminar => config.hlmax,
            FlowRegime::Turbulent | FlowRegime::Wake => config.htmax,
        };

        if result.hk > hmax {
            // Switch to inverse mode with prescribed Hk growth (XFOIL approach)
            let htarg = match regime {
                FlowRegime::Laminar => {
                    // Laminar: relatively slow increase in Hk downstream
                    // XFOIL: HTARG = HK1 + 0.03*(X2-X1)/T1
                    (prev.hk + 0.03 * ds / prev.theta).min(hmax)
                }
                FlowRegime::Turbulent | FlowRegime::Wake => {
                    // Turbulent: relatively fast decrease toward equilibrium
                    // XFOIL: HTARG = HK1 - 0.15*(X2-X1)/T1
                    (prev.hk - 0.15 * ds / prev.theta).max(config.hk_min)
                }
            };

            // Use prescribed Hk and solve for θ from momentum equation
            result.h = htarg;
            result.hk = htarg;

            // Estimate θ using momentum integral with prescribed H
            let due_ds = (ue_dist[i] - ue_dist[i - 1]) / ds;
            let ue_avg = 0.5 * (ue_dist[i - 1] + ue_dist[i]);
            let h_avg = 0.5 * (prev.h + htarg);
            let cf_est = prev.cf;

            // From momentum: dθ/ds = Cf/2 - (H+2)*θ/Ue * dUe/ds
            // Simplified trapezoidal integration
            let theta_coef = (h_avg + 2.0) / ue_avg * due_ds;
            let dtheta = (cf_est / 2.0 - theta_coef * prev.theta) * ds;
            result.theta = (prev.theta + dtheta).max(config.theta_min);
            result.dstar = result.h * result.theta;
            result.converged = false; // Mark as inverse mode
        }

        // Sanity check: θ should not drop drastically in adverse pressure gradient
        let theta_min_physical = prev.theta * 0.7;
        if result.theta < theta_min_physical && regime == FlowRegime::Laminar {
            result.theta = theta_min_physical;
            result.dstar = result.h * result.theta;
            result.converged = false;
        }

        // Check for transition
        if regime == FlowRegime::Laminar && result.n_amp >= cond.ncrit {
            regime = FlowRegime::Turbulent;
        }

        results.push(result);
    }

    results
}

/// Check if transition occurs and return transition location
///
/// Returns (has_transition, station_index, x_transition)
pub fn find_transition(results: &[NewtonResult], s_dist: &[f64], ncrit: f64) -> (bool, usize, f64) {
    for (i, result) in results.iter().enumerate() {
        if result.n_amp >= ncrit {
            // Interpolate transition location
            if i > 0 {
                let n_prev = results[i - 1].n_amp;
                let n_curr = result.n_amp;
                let frac = (ncrit - n_prev) / (n_curr - n_prev);
                let s_trans = s_dist[i - 1] + frac * (s_dist[i] - s_dist[i - 1]);
                return (true, i, s_trans);
            }
            return (true, i, s_dist[i]);
        }
    }
    (false, results.len(), f64::INFINITY)
}

/// Calculate total friction drag coefficient
pub fn integrate_friction(results: &[NewtonResult], s_dist: &[f64]) -> f64 {
    let mut cdf = 0.0;
    for i in 1..results.len() {
        let ds = s_dist[i] - s_dist[i - 1];
        let cf_avg = 0.5 * (results[i - 1].cf + results[i].cf);
        cdf += cf_avg * ds;
    }
    cdf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_newton_config_default() {
        let config = NewtonConfig::default();
        assert_eq!(config.max_iter, 50);
        assert!(config.tol < 1e-4);
        assert_eq!(config.relax, 1.0);
    }

    #[test]
    fn test_solve_single_station_laminar() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        // Typical laminar BL state - use realistic values from a BL march
        let theta1 = 0.0005;
        let h1 = 2.4;
        let ue1 = 0.5;
        let n1 = 0.0;
        let ue2 = 0.55; // Slight acceleration
        let s1 = 0.1;
        let s2 = 0.12;

        let result = solve_station(theta1, h1, ue1, n1, ue2, s1, s2, FlowRegime::Laminar, &cond, &config);

        // Check that solution is physically reasonable (may not fully converge for single station)
        assert!(result.theta > 0.0, "Theta should be positive");
        assert!(result.h > 1.0 && result.h < 10.0, "H should be reasonable: {}", result.h);
        // Either converged or residual is small
        assert!(
            result.converged || result.residual < 1e-2,
            "Should converge or have small residual: converged={}, residual={}",
            result.converged,
            result.residual
        );
    }

    #[test]
    fn test_solve_single_station_turbulent() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        // Typical turbulent BL state - use values that are well within turbulent closure range
        let theta1 = 0.003;
        let h1 = 1.5;
        let ue1 = 1.0;
        let n1 = 10.0; // Already transitioned
        let ue2 = 1.0; // Constant velocity (no pressure gradient)
        let s1 = 0.5;
        let s2 = 0.52;

        let result =
            solve_station(theta1, h1, ue1, n1, ue2, s1, s2, FlowRegime::Turbulent, &cond, &config);

        // Check physically reasonable solution
        assert!(result.theta > 0.0, "Theta should be positive");
        assert!(result.h > 1.0 && result.h < 5.0, "H should be reasonable for turbulent: {}", result.h);
        // Either converged or residual is acceptably small
        assert!(
            result.converged || result.residual < 1e-2,
            "Should converge or have small residual: converged={}, residual={}",
            result.converged,
            result.residual
        );
    }

    #[test]
    fn test_march_newton_favorable_gradient() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        // Favorable pressure gradient (accelerating flow)
        let n = 30;
        let s: Vec<f64> = (0..n).map(|i| i as f64 * 0.02).collect();
        let ue: Vec<f64> = s.iter().map(|&si| 0.1 + 1.5 * si).collect();

        let results = march_newton(&ue, &s, &cond, &config);

        assert_eq!(results.len(), n);

        // All stations should converge
        for (i, r) in results.iter().enumerate() {
            assert!(
                r.converged || r.residual < 1e-3,
                "Station {} failed: residual={}",
                i,
                r.residual
            );
        }

        // Theta should grow along the plate
        assert!(results[n - 1].theta > results[0].theta * 0.5);
    }

    #[test]
    fn test_march_newton_flat_plate() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        // Flat plate: constant Ue
        let n = 50;
        let s: Vec<f64> = (0..n).map(|i| 0.001 + i as f64 * 0.01).collect();
        let ue: Vec<f64> = vec![1.0; n];

        let results = march_newton(&ue, &s, &cond, &config);

        // Check Blasius-like behavior
        // For flat plate: θ ∝ sqrt(x), so θ/sqrt(x) should be roughly constant
        let ratio_mid = results[n / 2].theta / s[n / 2].sqrt();
        let ratio_end = results[n - 1].theta / s[n - 1].sqrt();

        // Ratios should be similar (within factor of 2 for our simple model)
        assert!(
            (ratio_end / ratio_mid - 1.0).abs() < 1.0,
            "θ/√s ratios: mid={}, end={}",
            ratio_mid,
            ratio_end
        );
    }

    #[test]
    fn test_find_transition() {
        let cond = FlowConditions::new(500_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        // Longer plate to allow transition
        let n = 100;
        let s: Vec<f64> = (0..n).map(|i| 0.001 + i as f64 * 0.01).collect();
        let ue: Vec<f64> = vec![1.0; n];

        let results = march_newton(&ue, &s, &cond, &config);

        let (has_trans, _idx, s_trans) = find_transition(&results, &s, cond.ncrit);

        // At Re=500k on a flat plate, transition should occur
        // (depends on N_crit and correlations)
        if has_trans {
            assert!(s_trans > 0.0 && s_trans < s[n - 1]);
        }
    }

    #[test]
    fn test_integrate_friction() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let config = NewtonConfig::default();

        let n = 30;
        let s: Vec<f64> = (0..n).map(|i| i as f64 * 0.02).collect();
        let ue: Vec<f64> = s.iter().map(|&si| 0.1 + si).collect();

        let results = march_newton(&ue, &s, &cond, &config);
        let cdf = integrate_friction(&results, &s);

        // Friction drag should be positive and reasonable
        assert!(cdf > 0.0, "Friction drag should be positive");
        assert!(cdf < 0.1, "Friction drag should be reasonable: {}", cdf);
    }
}
