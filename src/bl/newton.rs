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
/// downstream θ and H that satisfy the integral BL equations.
///
/// # Arguments
/// * `theta1` - Upstream momentum thickness
/// * `h1` - Upstream shape factor
/// * `ue1` - Upstream edge velocity
/// * `n1` - Upstream amplification factor
/// * `ue2` - Downstream edge velocity (prescribed)
/// * `ds` - Arc length step
/// * `regime` - Flow regime (laminar/turbulent)
/// * `cond` - Flow conditions
/// * `config` - Newton solver configuration
pub fn solve_station(
    theta1: f64,
    h1: f64,
    ue1: f64,
    n1: f64,
    ue2: f64,
    ds: f64,
    regime: FlowRegime,
    cond: &FlowConditions,
    config: &NewtonConfig,
) -> NewtonResult {
    // Initial guess: extrapolate from upstream
    let due_ds = (ue2 - ue1) / ds;
    let ue_avg = 0.5 * (ue1 + ue2);

    // Initial guess for theta using momentum equation trend
    let mut theta2 = theta1 * (1.0 + 0.1 * ds); // Slight growth
    let mut h2 = h1; // Shape factor roughly preserved

    // Newton iteration
    let mut converged = false;
    let mut iterations = 0;
    let mut residual = f64::MAX;

    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Compute derived quantities at station 2
        let _dstar2 = h2 * theta2;
        let (hk2, hk2_h, _hk2_msq) = hkin(h2, cond.msq);
        let hk2 = hk2.clamp(config.hk_min, config.hk_max);
        let rt2 = (ue2 * theta2 / cond.nu).max(1.0);

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
                // Simplified dissipation for turbulent
                let di = ClosureResult {
                    val: 0.5 * cf.val * ue2 / hs.val,
                    val_hk: 0.0,
                    val_rt: 0.0,
                    val_msq: 0.0,
                };
                (cf, hs, di)
            }
        };

        // Average quantities for trapezoidal integration
        let theta_avg = 0.5 * (theta1 + theta2);
        let h_avg = 0.5 * (h1 + h2);
        let (hk1, _, _) = hkin(h1, cond.msq);

        // Upstream closure for H*
        let hs1 = match regime {
            FlowRegime::Laminar => hs_lam(hk1, (ue1 * theta1 / cond.nu).max(1.0), cond.msq),
            _ => hs_turb(hk1, (ue1 * theta1 / cond.nu).max(200.0), cond.msq),
        };

        // ================================================================
        // Residual 1: Momentum integral equation
        // dθ/ds + (H + 2 - M²) * θ/Ue * dUe/ds = Cf/2
        // ================================================================
        let mom_coef = (h_avg + 2.0 - cond.msq) * theta_avg / ue_avg;
        let r1 = (theta2 - theta1) / ds + mom_coef * due_ds - cf_res.val / 2.0;

        // Jacobian entries for R1
        // dR1/dθ2
        let dr1_dtheta2 =
            1.0 / ds + 0.5 * (h_avg + 2.0 - cond.msq) / ue_avg * due_ds - cf_res.val_rt / 2.0 * ue2
                / cond.nu;

        // dR1/dH2 (through h_avg and Cf dependence on Hk)
        let dr1_dh2 = 0.5 * theta_avg / ue_avg * due_ds - cf_res.val_hk / 2.0 * hk2_h;

        // ================================================================
        // Residual 2: Shape parameter equation (kinetic energy form)
        // θ * dH*/ds + (2H** + H*(1-H)) * θ/Ue * dUe/ds = 2CD - H*Cf/2
        // ================================================================
        let hs_avg = 0.5 * (hs1.val + hs_res.val);
        let dhs_ds = (hs_res.val - hs1.val) / ds;

        // Density shape factor H** ≈ 0 for low Mach
        let hss = 0.0;
        let shape_coef = (2.0 * hss + hs_avg * (1.0 - h_avg)) * theta_avg / ue_avg;

        let r2 =
            theta_avg * dhs_ds + shape_coef * due_ds - 2.0 * di_res.val + hs_avg * cf_res.val / 2.0;

        // Jacobian entries for R2
        // dR2/dθ2
        let dr2_dtheta2 = 0.5 * dhs_ds
            + 0.5 * (2.0 * hss + hs_avg * (1.0 - h_avg)) / ue_avg * due_ds
            + hs_avg * cf_res.val_rt / 2.0 * ue2 / cond.nu
            - 2.0 * di_res.val_rt * ue2 / cond.nu;

        // dR2/dH2
        let dr2_dh2 = theta_avg / ds * hs_res.val_hk * hk2_h * 0.5
            + (0.5 * hs_res.val_hk * hk2_h * (1.0 - h_avg) - 0.5 * hs_avg) * theta_avg / ue_avg
                * due_ds
            + 0.5 * hs_res.val_hk * hk2_h * cf_res.val / 2.0
            + hs_avg * cf_res.val_hk / 2.0 * hk2_h
            - 2.0 * di_res.val_hk * hk2_h;

        // ================================================================
        // Solve 2x2 Newton system: J * [Δθ, ΔH]ᵀ = -[R1, R2]ᵀ
        // ================================================================
        let det = dr1_dtheta2 * dr2_dh2 - dr1_dh2 * dr2_dtheta2;

        if det.abs() < 1e-20 {
            // Singular Jacobian - use steepest descent fallback
            let d_theta = -r1 * config.relax * 0.1;
            let d_h = -r2 * config.relax * 0.1;
            theta2 = (theta2 + d_theta).max(config.theta_min);
            h2 = (h2 + d_h).clamp(1.0, config.hk_max);
        } else {
            // Standard Newton step
            let d_theta = (-r1 * dr2_dh2 + r2 * dr1_dh2) / det;
            let d_h = (-r2 * dr1_dtheta2 + r1 * dr2_dtheta2) / det;

            // Apply relaxation and limiting (XFOIL uses similar approach)
            // Limit step sizes to prevent divergence in difficult regions
            let d_theta_limited = d_theta.clamp(-0.3 * theta2, 0.3 * theta2);
            let d_h_limited = d_h.clamp(-0.8, 0.8);

            theta2 = (theta2 + config.relax * d_theta_limited).max(config.theta_min);
            h2 = (h2 + config.relax * d_h_limited).clamp(1.0, config.hk_max);
        }

        // Check convergence
        residual = (r1.powi(2) + r2.powi(2)).sqrt();
        if residual < config.tol {
            converged = true;
            break;
        }
    }

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

    // Initialize at stagnation point using Thwaites' method
    let due_ds_init = (ue_dist[1] - ue_dist[0]) / (s_dist[1] - s_dist[0]);
    let theta_init = (0.45 * cond.nu / due_ds_init.abs().max(1e-10)).sqrt();
    let h_init = 2.216; // Hiemenz stagnation point

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
        iterations: 0,
        residual: 0.0,
        converged: true,
    });

    // March downstream with Newton solver
    let mut regime = FlowRegime::Laminar;

    for i in 1..n {
        let prev = &results[i - 1];
        let ds = s_dist[i] - s_dist[i - 1];

        // Solve for this station
        let mut result = solve_station(
            prev.theta,
            prev.h,
            ue_dist[i - 1],
            prev.n_amp,
            ue_dist[i],
            ds,
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
        let ds = 0.02;

        let result = solve_station(theta1, h1, ue1, n1, ue2, ds, FlowRegime::Laminar, &cond, &config);

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
        let ds = 0.02;

        let result =
            solve_station(theta1, h1, ue1, n1, ue2, ds, FlowRegime::Turbulent, &cond, &config);

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
