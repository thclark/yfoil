//! Wake boundary layer model
//!
//! The wake is the continuation of the boundary layer behind the trailing edge.
//! Unlike the airfoil surface, the wake has no solid boundary, so:
//! - Skin friction Cf = 0
//! - The wake carries the combined momentum defect from upper and lower surfaces
//! - Shape factor H approaches 1.0 far downstream (velocity defect profile)
//!
//! XFOIL models the wake using the same integral equations with Cf = 0,
//! marching the wake downstream for several chord lengths.
//!
//! References:
//! - Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
//! - Cebeci & Bradshaw, "Momentum Transfer in Boundary Layers"

use crate::bl::{hkin, hs_turb, FlowConditions, NewtonConfig, NewtonResult};

/// Wake configuration
#[derive(Debug, Clone, Copy)]
pub struct WakeConfig {
    /// Number of wake stations
    pub n_wake: usize,
    /// Wake length in chord lengths
    pub wake_length: f64,
    /// Initial wake spacing (fraction of last surface panel)
    pub initial_spacing: f64,
    /// Growth ratio for wake panel spacing
    pub growth_ratio: f64,
}

impl Default for WakeConfig {
    fn default() -> Self {
        Self {
            n_wake: 20,
            wake_length: 1.0, // 1 chord length downstream
            initial_spacing: 0.02,
            growth_ratio: 1.15,
        }
    }
}

/// Wake state at trailing edge (initial conditions)
#[derive(Debug, Clone, Copy)]
pub struct WakeInitialState {
    /// Combined momentum thickness θ = θ_upper + θ_lower
    pub theta: f64,
    /// Combined displacement thickness δ* = δ*_upper + δ*_lower
    pub dstar: f64,
    /// Shape factor H = δ*/θ
    pub h: f64,
    /// Edge velocity at TE (average of upper/lower)
    pub ue: f64,
    /// x-coordinate of TE
    pub x_te: f64,
    /// y-coordinate of TE
    pub y_te: f64,
}

impl WakeInitialState {
    /// Create wake initial state from upper and lower surface TE values
    pub fn from_surfaces(
        upper_te: &NewtonResult,
        lower_te: &NewtonResult,
        ue_upper: f64,
        ue_lower: f64,
        x_te: f64,
        y_te: f64,
    ) -> Self {
        // Combine momentum thickness (additive)
        let theta = upper_te.theta + lower_te.theta;

        // Combine displacement thickness (additive)
        let dstar = upper_te.dstar + lower_te.dstar;

        // Shape factor from combined values
        let h = if theta > 1e-10 { dstar / theta } else { 1.0 };

        // Average edge velocity
        let ue = 0.5 * (ue_upper.abs() + ue_lower.abs());

        Self {
            theta,
            dstar,
            h,
            ue,
            x_te,
            y_te,
        }
    }
}

/// Generate wake station coordinates
///
/// The wake extends downstream from the TE, initially following the
/// bisector of the TE angle, then curving towards the freestream direction.
///
/// # Arguments
/// * `x_te` - Trailing edge x-coordinate
/// * `y_te` - Trailing edge y-coordinate
/// * `alpha` - Angle of attack (radians)
/// * `chord` - Airfoil chord length
/// * `config` - Wake configuration
///
/// # Returns
/// (x, y, s) coordinates and arc lengths for wake stations
pub fn generate_wake_coordinates(
    x_te: f64,
    y_te: f64,
    alpha: f64,
    chord: f64,
    config: &WakeConfig,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut x = Vec::with_capacity(config.n_wake);
    let mut y = Vec::with_capacity(config.n_wake);
    let mut s = Vec::with_capacity(config.n_wake);

    // Wake direction: follows freestream (cos α, sin α)
    // For small α, wake goes roughly in +x direction
    let wake_angle = -alpha; // Wake deflects opposite to lift
    let cos_wake = wake_angle.cos();
    let sin_wake = wake_angle.sin();

    // Generate stations with geometric spacing
    let mut arc_len = 0.0;
    let mut ds = config.initial_spacing * chord;

    for i in 0..config.n_wake {
        if i == 0 {
            x.push(x_te);
            y.push(y_te);
            s.push(0.0);
        } else {
            let x_new = x[i - 1] + ds * cos_wake;
            let y_new = y[i - 1] + ds * sin_wake;
            arc_len += ds;

            x.push(x_new);
            y.push(y_new);
            s.push(arc_len);

            // Grow spacing
            ds *= config.growth_ratio;
        }

        // Limit wake length
        if arc_len > config.wake_length * chord {
            break;
        }
    }

    (x, y, s)
}

/// Solve wake boundary layer equations
///
/// The wake uses modified BL equations with Cf = 0:
/// - Momentum: dθ/ds + (H + 1) * θ/Ue * dUe/ds = 0
/// - Shape: θ * dH*/ds + ... = 2*CD (dissipation only)
///
/// # Arguments
/// * `init` - Wake initial state at TE
/// * `s_wake` - Arc length coordinates of wake stations
/// * `ue_wake` - Edge velocity at wake stations (from potential flow)
/// * `cond` - Flow conditions
/// * `config` - Newton solver configuration
///
/// # Returns
/// Wake BL solution at each station
pub fn solve_wake(
    init: &WakeInitialState,
    s_wake: &[f64],
    ue_wake: &[f64],
    cond: &FlowConditions,
    _config: &NewtonConfig,
) -> Vec<NewtonResult> {
    let n = s_wake.len();
    if n == 0 {
        return vec![];
    }

    let mut results = Vec::with_capacity(n);

    // First station: initial conditions from TE
    let (hk_init, _, _) = hkin(init.h, cond.msq);
    let rt_init = (init.ue * init.theta / cond.nu).max(200.0);
    let hs_init = hs_turb(hk_init.max(1.0001), rt_init, cond.msq).val;

    results.push(NewtonResult {
        theta: init.theta,
        dstar: init.dstar,
        h: init.h,
        hk: hk_init,
        cf: 0.0, // No skin friction in wake
        cd: 0.0,
        hs: hs_init,
        n_amp: 0.0, // Wake is always turbulent
        iterations: 0,
        residual: 0.0,
        converged: true,
    });

    // March wake downstream
    for i in 1..n {
        let prev = &results[i - 1];
        let ds = s_wake[i] - s_wake[i - 1];
        let ue_avg = 0.5 * (ue_wake[i - 1] + ue_wake[i]);
        let due_ds = (ue_wake[i] - ue_wake[i - 1]) / ds;

        // Wake momentum equation: dθ/ds = -(H + 1) * θ/Ue * dUe/ds
        // (Cf = 0 in wake)
        let theta_rhs = -(prev.h + 1.0 - cond.msq) * prev.theta / ue_avg.max(0.01) * due_ds;
        let theta_new = (prev.theta + ds * theta_rhs).max(1e-10);

        // Shape factor evolution in wake
        // Wake H tends toward 1.0 (far-field velocity defect wake)
        // Use simple relaxation model: H -> 1 + (H_te - 1) * exp(-s/s_scale)
        let h_scale = 0.5; // Characteristic length for H decay
        let h_asymptote = 1.0; // Far-field wake shape factor
        let h_excess = (prev.h - h_asymptote).max(0.0);
        let h_new = h_asymptote + h_excess * (-ds / h_scale).exp();
        let h_new = h_new.max(1.0001); // Keep H > 1

        let dstar_new = h_new * theta_new;

        // Update derived quantities
        let (hk_new, _, _) = hkin(h_new, cond.msq);
        let hk_new = hk_new.max(1.0001);
        let rt_new = (ue_wake[i] * theta_new / cond.nu).max(200.0);
        let hs_new = hs_turb(hk_new, rt_new, cond.msq).val;

        results.push(NewtonResult {
            theta: theta_new,
            dstar: dstar_new,
            h: h_new,
            hk: hk_new,
            cf: 0.0,
            cd: 0.0,
            hs: hs_new,
            n_amp: 0.0,
            iterations: 0,
            residual: 0.0,
            converged: true,
        });
    }

    results
}

/// Calculate wake edge velocity from potential flow
///
/// The wake sees a velocity that approaches freestream far downstream.
/// Near the TE, the velocity is affected by the airfoil circulation.
///
/// # Arguments
/// * `x_wake` - Wake x-coordinates
/// * `y_wake` - Wake y-coordinates
/// * `gamma_total` - Total circulation (from panel method)
/// * `alpha` - Angle of attack
/// * `chord` - Airfoil chord
///
/// # Returns
/// Edge velocity at each wake station
pub fn wake_edge_velocity(
    x_wake: &[f64],
    _y_wake: &[f64],
    gamma_total: f64,
    alpha: f64,
    chord: f64,
) -> Vec<f64> {
    let n = x_wake.len();
    let mut ue = Vec::with_capacity(n);

    // Freestream velocity (unit magnitude)
    let u_inf = alpha.cos();
    let v_inf = alpha.sin();

    for i in 0..n {
        // Distance from TE (approximate wake centerline)
        let x_rel = x_wake[i] - x_wake[0];

        // Induced velocity from bound vortex (simplified)
        // For a point vortex at origin, v_induced ~ Γ/(2π*r) perpendicular to r
        // Wake is roughly along x-axis, so induced v is primarily in y
        let r = (x_rel.powi(2) + 0.01 * chord.powi(2)).sqrt();
        let v_induced = gamma_total / (2.0 * std::f64::consts::PI * r) * 0.5;

        // Total velocity magnitude
        let u_total = u_inf;
        let v_total = v_inf + v_induced;
        let ue_mag = (u_total.powi(2) + v_total.powi(2)).sqrt();

        ue.push(ue_mag.max(0.1)); // Ensure positive
    }

    ue
}

/// Calculate momentum thickness at infinity (Squire-Young formula)
///
/// The drag coefficient is related to the far-wake momentum thickness:
/// CD = 2 * θ_inf / chord
///
/// Using Squire-Young: θ_inf = θ_te * (Ue_te / U_inf)^((H_te + 5)/2)
pub fn squire_young_drag(theta_te: f64, h_te: f64, ue_te: f64, chord: f64) -> f64 {
    let exponent = (h_te + 5.0) / 2.0;
    let theta_inf = theta_te * ue_te.powf(exponent);
    2.0 * theta_inf / chord
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wake_initial_state() {
        let upper = NewtonResult {
            theta: 0.003,
            dstar: 0.006,
            h: 2.0,
            hk: 2.0,
            cf: 0.002,
            cd: 0.001,
            hs: 1.5,
            n_amp: 10.0,
            iterations: 5,
            residual: 1e-6,
            converged: true,
        };
        let lower = NewtonResult {
            theta: 0.002,
            dstar: 0.004,
            h: 2.0,
            hk: 2.0,
            cf: 0.002,
            cd: 0.001,
            hs: 1.5,
            n_amp: 10.0,
            iterations: 5,
            residual: 1e-6,
            converged: true,
        };

        let init = WakeInitialState::from_surfaces(&upper, &lower, 1.0, 1.0, 1.0, 0.0);

        assert_eq!(init.theta, 0.005); // Combined
        assert_eq!(init.dstar, 0.010); // Combined
        assert_eq!(init.h, 2.0); // δ*/θ = 0.01/0.005
        assert_eq!(init.ue, 1.0);
    }

    #[test]
    fn test_generate_wake_coordinates() {
        let config = WakeConfig::default();
        let (x, y, s) = generate_wake_coordinates(1.0, 0.0, 0.0, 1.0, &config);

        assert!(!x.is_empty());
        assert_eq!(x.len(), y.len());
        assert_eq!(x.len(), s.len());

        // First point should be at TE
        assert_eq!(x[0], 1.0);
        assert_eq!(y[0], 0.0);
        assert_eq!(s[0], 0.0);

        // Wake should extend downstream (x increasing)
        for i in 1..x.len() {
            assert!(x[i] > x[i - 1], "Wake should extend downstream");
            assert!(s[i] > s[i - 1], "Arc length should increase");
        }
    }

    #[test]
    fn test_solve_wake() {
        let init = WakeInitialState {
            theta: 0.005,
            dstar: 0.010,
            h: 2.0,
            ue: 1.0,
            x_te: 1.0,
            y_te: 0.0,
        };

        let config = WakeConfig::default();
        let (_, _, s) = generate_wake_coordinates(1.0, 0.0, 0.0, 1.0, &config);
        let ue: Vec<f64> = vec![1.0; s.len()]; // Constant velocity

        let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
        let newton_config = NewtonConfig::default();

        let results = solve_wake(&init, &s, &ue, &cond, &newton_config);

        assert_eq!(results.len(), s.len());

        // All wake stations should have Cf = 0
        for r in &results {
            assert_eq!(r.cf, 0.0, "Wake should have zero skin friction");
        }

        // Shape factor should decrease towards 1.0
        let h_init = results[0].h;
        let h_final = results.last().unwrap().h;
        assert!(h_final < h_init, "Wake H should decrease towards 1.0");
        assert!(h_final >= 1.0, "Wake H should stay >= 1.0");
    }

    #[test]
    fn test_squire_young_drag() {
        // Typical values
        let theta_te = 0.005;
        let h_te = 2.0;
        let ue_te = 1.0;
        let chord = 1.0;

        let cd = squire_young_drag(theta_te, h_te, ue_te, chord);

        // CD should be positive and reasonable
        assert!(cd > 0.0);
        assert!(cd < 0.1);

        // At ue_te = 1, CD = 2 * θ_te / chord = 0.01
        assert!((cd - 0.01).abs() < 0.001);
    }

    #[test]
    fn test_wake_edge_velocity() {
        let x = vec![1.0, 1.5, 2.0, 3.0];
        let y = vec![0.0, 0.0, 0.0, 0.0];
        let gamma = 0.5; // Typical circulation
        let alpha = 0.0;
        let chord = 1.0;

        let ue = wake_edge_velocity(&x, &y, gamma, alpha, chord);

        assert_eq!(ue.len(), x.len());

        // All velocities should be positive
        for &u in &ue {
            assert!(u > 0.0);
        }

        // Velocity should approach freestream far downstream
        assert!((ue.last().unwrap() - 1.0).abs() < 0.5);
    }
}
