use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    
    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl = march_newton(&ue_upper, &s_upper, &cond, &config);
    
    // Print every station
    println!("{:>5} {:>10} {:>10} {:>12} {:>12} {:>8} {:>8} {:>8} {:>6}",
             "idx", "x", "s", "theta", "dstar", "H", "Cf", "Ue", "conv");
    
    for (i, r) in bl.iter().enumerate() {
        if i < x_upper.len() {
            let conv = if r.converged { "Y" } else { "N" };
            println!("{:>5} {:>10.5} {:>10.5} {:>12.5e} {:>12.5e} {:>8.4} {:>8.5} {:>8.4} {:>6}",
                     i, x_upper[i], s_upper[i], r.theta, r.dstar, r.h, r.cf, r.ue, conv);
        }
    }
    
    // Compare with XFOIL at key stations
    println!("\n=== XFOIL Reference (upper surface) ===");
    println!("Near LE (s~0.02): theta=3.3e-5, H=2.32, Cf=0.0056");
    println!("Mid (s~0.5): theta~5e-4, H~2.6");
    println!("TE (s~1.02): theta=0.002, dstar=0.003, H=1.59, Cf=0.002");
}
