use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{
    find_stagnation_point, solve_boundary_layer, calculate_friction_drag,
    ViscalConfig,
};

fn main() {
    // Generate geometry
    let geom = naca_4digit("0012", 160).expect("Failed");
    let airfoil = create_paneled_airfoil(&geom);

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);

    // Solve at alpha=2, Re=1e6
    let alpha_deg = 2.0f64;
    let alpha_rad = alpha_deg.to_radians();
    let cond = FlowConditions::new(1e6, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    // Get inviscid velocity at this alpha
    let velocity = inviscid.velocity_at_alpha(alpha_rad);

    // Find stagnation point
    let stag_idx = find_stagnation_point(&airfoil, &velocity);

    // Get total gamma for wake
    let gamma_total = inviscid.gam_0.iter().sum::<f64>() * alpha_rad.cos()
        + inviscid.gam_90.iter().sum::<f64>() * alpha_rad.sin();

    // Solve BL
    let bl = solve_boundary_layer(
        &airfoil, &velocity, stag_idx, alpha_rad, gamma_total,
        &cond, &config.newton, &config.wake
    );

    println!("=== YFoil BL Debug at alpha={} deg, Re={:.0e} ===", alpha_deg, cond.reynolds);
    println!("Stagnation index: {}", stag_idx);
    println!("Upper stations: {}", bl.upper.len());
    println!("Lower stations: {}", bl.lower.len());

    // Dump upper surface BL
    println!("\n--- Upper Surface BL (s, x, ue, theta, cf, h, tau=0.5*ue²*cf) ---");
    for (i, r) in bl.upper.iter().enumerate() {
        let x = bl.x_upper[i];
        let s = bl.s_upper[i];
        let tau = 0.5 * r.ue.powi(2) * r.cf;
        println!("{:3} s={:.5} x={:.5} ue={:.5} theta={:.2e} cf={:.4e} h={:.4} tau={:.4e}",
            i, s, x, r.ue, r.theta, r.cf, r.h, tau);
    }

    // Calculate friction drag both ways
    let cdf = calculate_friction_drag(&bl, alpha_rad);

    // Manual integration using arc length (old way) for comparison
    let mut cdf_arc = 0.0;
    for i in 1..bl.upper.len() {
        let ds = bl.s_upper[i] - bl.s_upper[i - 1];
        let cf_avg = 0.5 * (bl.upper[i - 1].cf + bl.upper[i].cf);
        cdf_arc += cf_avg * ds;
    }
    for i in 1..bl.lower.len() {
        let ds = bl.s_lower[i] - bl.s_lower[i - 1];
        let cf_avg = 0.5 * (bl.lower[i - 1].cf + bl.lower[i].cf);
        cdf_arc += cf_avg * ds;
    }

    // Manual integration using x-coord and Ue (new way)
    let ca = alpha_rad.cos();
    let mut cdf_manual = 0.0;
    for i in 1..bl.upper.len() {
        let dx = (bl.x_upper[i] - bl.x_upper[i - 1]) * ca;
        let tau_prev = 0.5 * bl.upper[i-1].ue.powi(2) * bl.upper[i-1].cf;
        let tau_curr = 0.5 * bl.upper[i].ue.powi(2) * bl.upper[i].cf;
        cdf_manual += (tau_prev + tau_curr) * dx;
    }
    for i in 1..bl.lower.len() {
        let dx = (bl.x_lower[i] - bl.x_lower[i - 1]) * ca;
        let tau_prev = 0.5 * bl.lower[i-1].ue.powi(2) * bl.lower[i-1].cf;
        let tau_curr = 0.5 * bl.lower[i].ue.powi(2) * bl.lower[i].cf;
        cdf_manual += (tau_prev + tau_curr) * dx;
    }

    println!("\n=== CDf comparisons ===");
    println!("CDf (api):        {:.6}", cdf);
    println!("CDf (manual new): {:.6}", cdf_manual);
    println!("CDf (arc length): {:.6}", cdf_arc);
    println!("XFOIL CDf:        0.00433");
}
