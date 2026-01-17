use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::{FlowConditions, NewtonConfig, march_newton, rtheta_crit};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    println!("=== Early θ Growth Analysis ===");
    println!("Re = 1e6, ν = {:e}", cond.nu);
    println!();
    println!("{:>4} {:>8} {:>10} {:>10} {:>8} {:>8} {:>10} {:>10} {:>8}",
        "Stn", "x/c", "θ", "Ue", "Hk", "Rθ", "log(Rθ)", "log(Rcrit)", "Margin");
    println!("{:-<100}", "");
    
    for i in 0..25.min(bl_upper.len()) {
        let theta = bl_upper[i].theta;
        let ue = ue_upper[i];
        let hk = bl_upper[i].hk;
        let rt = ue * theta / cond.nu;
        let log_rt = rt.log10();
        let log_rcrit = rtheta_crit(hk);
        let margin = log_rt - log_rcrit;
        let marker = if margin > 0.0 { "*" } else { " " };
        
        println!("{:>4} {:>8.4} {:>10.2e} {:>10.4} {:>8.4} {:>8.1} {:>10.4} {:>10.4} {:>7.4}{}",
            i, x_upper[i], theta, ue, hk, rt, log_rt, log_rcrit, margin, marker);
    }
    
    // Compare with what XFOIL might have at x/c = 0.05 and 0.10
    println!();
    println!("=== Expected XFOIL values at key stations ===");
    println!("At x/c ≈ 0.05: θ ~ 1.2e-4, Rθ ~ 150, for Re=1e6");
    println!("At x/c ≈ 0.10: θ ~ 2.5e-4, Rθ ~ 300, for Re=1e6");
    println!();
    println!("The critical Rθ for Hk ≈ 2.5 is about 10^{:.2} ≈ {:.0}", 
        rtheta_crit(2.5), 10.0_f64.powf(rtheta_crit(2.5)));
}
