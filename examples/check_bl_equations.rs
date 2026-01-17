// Diagnostic to compare equation terms at a specific station
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::{FlowConditions, NewtonConfig, march_newton, hkin, cf_lam, hs_lam, di_lam};

fn main() {
    // Same setup as XFOIL comparison
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    // Find station closest to x/c = 0.50
    let target_x = 0.50;
    let idx = x_upper.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| ((*a - target_x).abs()).partial_cmp(&((*b - target_x).abs())).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    println!("=== Equation Terms at x/c = {:.4} (station {}) ===", x_upper[idx], idx);
    println!();

    // Previous station values
    let prev = &bl_upper[idx - 1];
    let curr = &bl_upper[idx];

    let s1 = s_upper[idx - 1];
    let s2 = s_upper[idx];
    let ue1 = ue_upper[idx - 1];
    let ue2 = ue_upper[idx];

    let theta1 = prev.theta;
    let theta2 = curr.theta;
    let h1 = prev.h;
    let h2 = curr.h;

    let (hk1, _, _) = hkin(h1, cond.msq);
    let (hk2, _, _) = hkin(h2, cond.msq);

    let rt1 = ue1 * theta1 / cond.nu;
    let rt2 = ue2 * theta2 / cond.nu;

    let cf1 = cf_lam(hk1, rt1, cond.msq);
    let cf2 = cf_lam(hk2, rt2, cond.msq);
    let hs1 = hs_lam(hk1, rt1, cond.msq);
    let hs2 = hs_lam(hk2, rt2, cond.msq);
    let di1 = di_lam(hk1, rt1);
    let di2 = di_lam(hk2, rt2);

    println!("Station values:");
    println!("  s1 = {:.6e}, s2 = {:.6e}, ds = {:.6e}", s1, s2, s2 - s1);
    println!("  Ue1 = {:.6}, Ue2 = {:.6}", ue1, ue2);
    println!("  θ1 = {:.6e}, θ2 = {:.6e}", theta1, theta2);
    println!("  H1 = {:.6}, H2 = {:.6}", h1, h2);
    println!("  Hk1 = {:.6}, Hk2 = {:.6}", hk1, hk2);
    println!("  Rt1 = {:.2}, Rt2 = {:.2}", rt1, rt2);
    println!();

    println!("Closure values:");
    println!("  Cf1 = {:.6e}, Cf2 = {:.6e}", cf1.val, cf2.val);
    println!("  H*1 = {:.6}, H*2 = {:.6}", hs1.val, hs2.val);
    println!("  DI1 = {:.6e}, DI2 = {:.6e}", di1.val, di2.val);
    println!();

    // Logarithmic terms
    let xlog = (s2 / s1).ln();
    let ulog = (ue2 / ue1).ln();
    let tlog = (theta2 / theta1).ln();
    let hlog = (hs2.val / hs1.val).ln();

    println!("Logarithmic terms:");
    println!("  XLOG = ln(s2/s1) = {:.6}", xlog);
    println!("  ULOG = ln(Ue2/Ue1) = {:.6}", ulog);
    println!("  TLOG = ln(θ2/θ1) = {:.6}", tlog);
    println!("  HLOG = ln(H*2/H*1) = {:.6}", hlog);
    println!();

    // Momentum equation terms
    let h_avg = 0.5 * (h1 + h2);
    let theta_avg = 0.5 * (theta1 + theta2);
    let s_avg = 0.5 * (s1 + s2);
    let btmp = h_avg + 2.0 - cond.msq; // For momentum equation

    // Midpoint Cf
    let hk_avg = 0.5 * (hk1 + hk2);
    let rt_avg = 0.5 * (rt1 + rt2);
    let cf_mid = cf_lam(hk_avg, rt_avg, cond.msq);

    // CFX for momentum equation
    let cfx_mom = 0.5 * cf_mid.val * s_avg / theta_avg
        + 0.25 * (cf1.val * s1 / theta1 + cf2.val * s2 / theta2);

    let r1_val = tlog + btmp * ulog - xlog * 0.5 * cfx_mom;

    println!("Momentum equation:");
    println!("  H_avg = {:.6}, θ_avg = {:.6e}", h_avg, theta_avg);
    println!("  BTMP = H + 2 - M² = {:.6}", btmp);
    println!("  CFM (midpoint) = {:.6e}", cf_mid.val);
    println!("  CFX = {:.6e}", cfx_mom);
    println!("  R1 = TLOG + BTMP*ULOG - XLOG*0.5*CFX = {:.6e}", r1_val);
    println!();

    // Shape equation terms (with upwinding)
    let hdcon = 5.0 / hk2.powi(2);
    let arg = ((hk2 - 1.0) / (hk1 - 1.0).max(0.01)).abs();
    let hl = arg.ln();
    let hlsq = hl.powi(2).min(15.0);
    let ehh = (-hlsq * hdcon).exp();
    let upw = 1.0 - 0.5 * ehh;

    let hs_avg = 0.5 * (hs1.val + hs2.val);
    let hss = 0.0; // density shape factor (incompressible)
    let btmp_h = 2.0 * hss / hs_avg + 1.0 - h_avg;

    let xot1 = s1 / theta1;
    let xot2 = s2 / theta2;
    let dix = (1.0 - upw) * di1.val * xot1 + upw * di2.val * xot2;
    let cfx_shape = (1.0 - upw) * cf1.val * xot1 + upw * cf2.val * xot2;

    let r2_val = hlog + btmp_h * ulog + xlog * (0.5 * cfx_shape - dix);

    println!("Shape parameter equation:");
    println!("  UPW = {:.6}", upw);
    println!("  H*_avg = {:.6}", hs_avg);
    println!("  BTMP_H = 2*H**/H* + 1 - H = {:.6}", btmp_h);
    println!("  X1/θ1 = {:.2}, X2/θ2 = {:.2}", xot1, xot2);
    println!("  DIX = {:.6e}", dix);
    println!("  CFX (shape) = {:.6e}", cfx_shape);
    println!("  R2 = HLOG + BTMP_H*ULOG + XLOG*(0.5*CFX - DIX) = {:.6e}", r2_val);
    println!();

    println!("Solution quality:");
    println!("  Iterations = {}", curr.iterations);
    println!("  Residual = {:.2e}", curr.residual);
    println!("  Converged = {}", curr.converged);
}
