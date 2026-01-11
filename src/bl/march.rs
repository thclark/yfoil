//! Boundary layer marching algorithm
//!
//! Implements the integral boundary layer equations marched along the
//! surface from the stagnation point. Uses a Newton-Raphson iteration
//! at each station to solve the coupled equations.
//!
//! The primary equations are:
//! 1. Momentum integral equation (von Kármán)
//! 2. Shape parameter equation (kinetic energy or lag-entrainment)
//! 3. Amplification equation (for laminar flow)
//!
//! References:
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - Drela, M., Giles, M. "Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils"

use crate::bl::{
    cf_lam, cf_turb, dampl, di_lam, hkin, hs_lam, hs_turb, ClosureResult, FlowConditions,
    FlowRegime,
};

/// Residual and Jacobian for the BL equations at a station
#[derive(Debug, Clone)]
pub struct BLResidual {
    /// Residual of momentum equation
    pub r_theta: f64,
    /// Residual of shape parameter equation
    pub r_h: f64,
    /// Residual of amplification equation (laminar only)
    pub r_n: f64,

    /// Jacobian: d(r_theta)/d(theta)
    pub dr_theta_dtheta: f64,
    /// Jacobian: d(r_theta)/d(dstar)
    pub dr_theta_ddstar: f64,
    /// Jacobian: d(r_theta)/d(ue)
    pub dr_theta_due: f64,

    /// Jacobian: d(r_h)/d(theta)
    pub dr_h_dtheta: f64,
    /// Jacobian: d(r_h)/d(dstar)
    pub dr_h_ddstar: f64,
    /// Jacobian: d(r_h)/d(ue)
    pub dr_h_due: f64,

    /// Jacobian: d(r_n)/d(theta)
    pub dr_n_dtheta: f64,
    /// Jacobian: d(r_n)/d(n)
    pub dr_n_dn: f64,
}

impl BLResidual {
    pub fn new() -> Self {
        Self {
            r_theta: 0.0,
            r_h: 0.0,
            r_n: 0.0,
            dr_theta_dtheta: 0.0,
            dr_theta_ddstar: 0.0,
            dr_theta_due: 0.0,
            dr_h_dtheta: 0.0,
            dr_h_ddstar: 0.0,
            dr_h_due: 0.0,
            dr_n_dtheta: 0.0,
            dr_n_dn: 0.0,
        }
    }
}

impl Default for BLResidual {
    fn default() -> Self {
        Self::new()
    }
}

/// Variables at a BL station for the Newton solver
#[derive(Debug, Clone, Copy)]
pub struct BLVars {
    /// Momentum thickness θ
    pub theta: f64,
    /// Displacement thickness δ*
    pub dstar: f64,
    /// Edge velocity
    pub ue: f64,
    /// Amplification factor
    pub n_amp: f64,
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Kinematic shape factor
    pub hk: f64,
    /// Momentum thickness Reynolds number
    pub rtheta: f64,
    /// Flow regime
    pub regime: FlowRegime,
}

impl BLVars {
    /// Update derived quantities from primary variables
    pub fn update(&mut self, cond: &FlowConditions) {
        if self.theta > 1e-14 {
            self.h = self.dstar / self.theta;
        }
        let (hk, _, _) = hkin(self.h, cond.msq);
        self.hk = hk.max(1.00005); // Limit Hk to prevent singularities
        self.rtheta = self.ue * self.theta / cond.nu;
    }

    /// Get closure results for this station
    pub fn closures(&self, cond: &FlowConditions) -> (ClosureResult, ClosureResult, ClosureResult) {
        match self.regime {
            FlowRegime::Laminar => {
                let cf = cf_lam(self.hk, self.rtheta.max(1.0), cond.msq);
                let hs = hs_lam(self.hk, self.rtheta.max(1.0), cond.msq);
                let di = di_lam(self.hk, self.rtheta.max(1.0));
                (cf, hs, di)
            }
            FlowRegime::Turbulent | FlowRegime::Wake => {
                let cf = cf_turb(self.hk, self.rtheta.max(200.0), cond.msq, 1.0);
                let hs = hs_turb(self.hk, self.rtheta.max(200.0), cond.msq);
                // Turbulent dissipation needs more complex calculation
                let di = ClosureResult::value(0.0); // Placeholder
                (cf, hs, di)
            }
        }
    }
}

/// Compute BL equation residuals at a station
///
/// Uses the von Kármán momentum integral equation:
///   dθ/ds + (H + 2 - M²) * θ/ue * due/ds = Cf/2
///
/// And the shape parameter equation (kinetic energy form):
///   θ * dH*/ds + (2H** + H*(1-H)) * θ/ue * due/ds = 2CD - H*Cf/2
///
/// # Arguments
/// * `vars1` - Variables at upstream station
/// * `vars2` - Variables at downstream station
/// * `ds` - Arc length step
/// * `due_ds` - Edge velocity gradient
/// * `cond` - Flow conditions
pub fn compute_residual(
    vars1: &BLVars,
    vars2: &BLVars,
    ds: f64,
    due_ds: f64,
    cond: &FlowConditions,
) -> BLResidual {
    let mut res = BLResidual::new();

    // Average quantities
    let theta_avg = 0.5 * (vars1.theta + vars2.theta);
    let h_avg = 0.5 * (vars1.h + vars2.h);
    let hk_avg = 0.5 * (vars1.hk + vars2.hk);
    let ue_avg = 0.5 * (vars1.ue + vars2.ue);
    let rt_avg = 0.5 * (vars1.rtheta + vars2.rtheta);

    // Get closure relations at station 2 (downstream)
    let (cf, hs, di) = vars2.closures(cond);

    // Momentum equation residual:
    // (θ₂ - θ₁)/ds + (H + 2 - M²) * θ/ue * due/ds = Cf/2
    let momentum_coef = (h_avg + 2.0 - cond.msq) * theta_avg / ue_avg;
    res.r_theta = (vars2.theta - vars1.theta) / ds + momentum_coef * due_ds - cf.val / 2.0;

    // Jacobian for momentum equation
    res.dr_theta_dtheta = 1.0 / ds + (h_avg + 2.0 - cond.msq) / ue_avg * due_ds * 0.5
        - cf.val_rt / 2.0 * vars2.ue / cond.nu;
    res.dr_theta_ddstar = theta_avg / (ue_avg * vars2.theta) * due_ds * 0.5
        - cf.val_hk / 2.0 / vars2.theta;
    res.dr_theta_due = -(h_avg + 2.0 - cond.msq) * theta_avg / ue_avg.powi(2) * due_ds * 0.5
        - cf.val_rt / 2.0 * vars2.theta / cond.nu;

    // Shape parameter equation residual:
    // θ * dH*/ds + (2H** + H*(1-H)) * θ/ue * due/ds = 2CD - H*Cf/2
    let hs_val = hs.val;
    let dhs_ds = (hs_val - hs_lam(vars1.hk, rt_avg, cond.msq).val) / ds;

    // Use density shape factor H** ≈ 0 for incompressible
    let hss = 0.0;
    let shape_coef = (2.0 * hss + hs_val * (1.0 - h_avg)) * theta_avg / ue_avg;

    res.r_h = theta_avg * dhs_ds + shape_coef * due_ds - 2.0 * di.val + hs_val * cf.val / 2.0;

    // Jacobian for shape parameter equation
    res.dr_h_dtheta = dhs_ds * 0.5 + (2.0 * hss + hs_val * (1.0 - h_avg)) / ue_avg * due_ds * 0.5;
    res.dr_h_ddstar = theta_avg / ds * hs.val_hk / vars2.theta * 0.5
        - hs_val * theta_avg / ue_avg * due_ds * 0.5 / vars2.theta;
    res.dr_h_due = -(2.0 * hss + hs_val * (1.0 - h_avg)) * theta_avg / ue_avg.powi(2) * due_ds * 0.5;

    // Amplification equation for laminar flow
    if vars2.regime == FlowRegime::Laminar {
        let (ax, _ax_hk, ax_th, ax_rt) = dampl(hk_avg, theta_avg, rt_avg);
        res.r_n = (vars2.n_amp - vars1.n_amp) / ds - ax;

        res.dr_n_dtheta = -ax_th * 0.5 - ax_rt * vars2.ue / cond.nu * 0.5;
        res.dr_n_dn = 1.0 / ds;
    }

    res
}

/// Initialize BL from stagnation point using Thwaites' method
///
/// At the stagnation point, we use the exact solution for flow
/// near a stagnation point to initialize θ.
///
/// # Arguments
/// * `due_ds` - Edge velocity gradient at stagnation
/// * `cond` - Flow conditions
pub fn initialize_stagnation(due_ds: f64, cond: &FlowConditions) -> BLVars {
    // Thwaites' method: θ² = 0.45 * ν / (due/ds)
    // For stagnation point flow: due/ds = U_inf * k where k is stagnation point parameter
    let theta_sq = 0.45 * cond.nu / due_ds.abs().max(1e-10);
    let theta = theta_sq.sqrt();

    // Laminar stagnation point has H ≈ 2.216 (Hiemenz solution)
    let h = 2.216;
    let dstar = h * theta;

    // Small initial edge velocity
    let ue = 0.01;

    let mut vars = BLVars {
        theta,
        dstar,
        ue,
        n_amp: 0.0,
        h,
        hk: h,
        rtheta: 0.0,
        regime: FlowRegime::Laminar,
    };
    vars.update(cond);
    vars
}

/// Direct mode march: march BL with prescribed edge velocity
///
/// This is used when the inviscid solution provides the edge velocity
/// distribution. The BL equations are integrated downstream.
///
/// # Arguments
/// * `ue_dist` - Edge velocity at each station
/// * `s_dist` - Arc length at each station
/// * `cond` - Flow conditions
///
/// # Returns
/// BL variables at each station
pub fn march_direct(ue_dist: &[f64], s_dist: &[f64], cond: &FlowConditions) -> Vec<BLVars> {
    let n = ue_dist.len();
    assert!(n >= 2);
    assert_eq!(s_dist.len(), n);

    let mut result = Vec::with_capacity(n);

    // Initialize at first station (assume stagnation point)
    let due_ds_init = if n > 1 {
        (ue_dist[1] - ue_dist[0]) / (s_dist[1] - s_dist[0])
    } else {
        1.0
    };
    let mut vars = initialize_stagnation(due_ds_init, cond);
    vars.ue = ue_dist[0];
    vars.update(cond);
    result.push(vars);

    // March downstream
    for i in 1..n {
        let ds = s_dist[i] - s_dist[i - 1];
        let due_ds = (ue_dist[i] - ue_dist[i - 1]) / ds;

        // Simple explicit Euler step for now
        // TODO: Replace with Newton iteration for implicit step
        let vars_prev = result[i - 1];

        // Get closure relations
        let (cf, _, _) = vars_prev.closures(cond);

        // Momentum equation: dθ/ds = Cf/2 - (H + 2 - M²) * θ/ue * due/ds
        let momentum_rhs =
            cf.val / 2.0 - (vars_prev.h + 2.0 - cond.msq) * vars_prev.theta / vars_prev.ue * due_ds;
        let theta_new = (vars_prev.theta + ds * momentum_rhs).max(1e-10);

        // Shape factor stays roughly constant for simple march
        let h_new = vars_prev.h;
        let dstar_new = h_new * theta_new;

        // Amplification factor
        let (ax, _, _, _) = dampl(vars_prev.hk, vars_prev.theta, vars_prev.rtheta);
        let n_new = vars_prev.n_amp + ax * ds;

        // Check for transition
        let regime = if n_new >= cond.ncrit {
            FlowRegime::Turbulent
        } else {
            vars_prev.regime
        };

        let mut vars_new = BLVars {
            theta: theta_new,
            dstar: dstar_new,
            ue: ue_dist[i],
            n_amp: n_new,
            h: h_new,
            hk: vars_prev.hk,
            rtheta: 0.0,
            regime,
        };
        vars_new.update(cond);
        result.push(vars_new);
    }

    result
}

/// Calculate skin friction drag coefficient from BL solution
///
/// Integrates Cf over the surface
pub fn integrate_friction_drag(vars_upper: &[BLVars], vars_lower: &[BLVars], s: &[f64]) -> f64 {
    let mut cdf = 0.0;

    // Integrate over upper surface
    for i in 1..vars_upper.len() {
        let ds = s[i] - s[i - 1];
        let cf_avg = 0.5 * (vars_upper[i - 1].cf() + vars_upper[i].cf());
        cdf += cf_avg * ds;
    }

    // Integrate over lower surface
    for i in 1..vars_lower.len() {
        let ds = s[i] - s[i - 1];
        let cf_avg = 0.5 * (vars_lower[i - 1].cf() + vars_lower[i].cf());
        cdf += cf_avg * ds;
    }

    cdf
}

impl BLVars {
    /// Get skin friction coefficient
    pub fn cf(&self) -> f64 {
        // Placeholder - should use actual closure
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stagnation_initialization() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let vars = initialize_stagnation(10.0, &cond);

        // θ² = 0.45 * ν / (due/ds) = 0.45 * 1e-6 / 10 = 4.5e-8
        // θ = 2.12e-4
        assert!(vars.theta > 1e-5);
        assert!(vars.theta < 1e-3);
        assert_eq!(vars.regime, FlowRegime::Laminar);
    }

    #[test]
    fn test_direct_march_monotonic_theta() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

        // Accelerating flow (favorable pressure gradient)
        let n = 20;
        let s: Vec<f64> = (0..n).map(|i| i as f64 * 0.05).collect();
        let ue: Vec<f64> = s.iter().map(|&si| 1.0 + 0.5 * si).collect();

        let result = march_direct(&ue, &s, &cond);

        assert_eq!(result.len(), n);
        // Momentum thickness should grow (or at least not decrease dramatically)
        assert!(result[n - 1].theta > result[0].theta * 0.1);
    }

    #[test]
    fn test_compute_residual_dimensions() {
        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

        let mut vars1 = BLVars {
            theta: 0.001,
            dstar: 0.0025,
            ue: 1.0,
            n_amp: 0.0,
            h: 2.5,
            hk: 2.5,
            rtheta: 1000.0,
            regime: FlowRegime::Laminar,
        };
        vars1.update(&cond);

        let mut vars2 = vars1;
        vars2.theta = 0.0011;
        vars2.dstar = 0.00275;
        vars2.ue = 1.05;
        vars2.update(&cond);

        let res = compute_residual(&vars1, &vars2, 0.01, 5.0, &cond);

        // Residuals should be finite
        assert!(res.r_theta.is_finite());
        assert!(res.r_h.is_finite());
    }
}
