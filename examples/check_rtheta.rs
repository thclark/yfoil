use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::bl::{FlowConditions, NewtonConfig, rtheta_crit};
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::march_newton;

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

    println!("=== Rθ vs Critical Rθ Analysis ===\n");
    println!("{:>5} {:>8} {:>8} {:>8} {:>12} {:>12} {:>12} {:>10}",
        "Stn", "x/c", "Hk", "θ", "Rθ", "log10(Rθ)", "log10(Rcrit)", "Margin");
    println!("{:-<95}", "");

    for i in 20..65 {
        if i >= bl_upper.len() { break; }
        let hk = bl_upper[i].hk;
        let theta = bl_upper[i].theta;
        let ue = ue_upper[i];
        let rt = ue * theta / cond.nu;
        let log_rt = rt.log10();
        let log_rcrit = rtheta_crit(hk);
        let margin = log_rt - log_rcrit;

        let marker = if margin > 0.0 { "*" } else { " " };

        println!("{:>5} {:>8.4} {:>8.4} {:>8.2e} {:>12.1} {:>12.4} {:>12.4} {:>9.4}{}",
            i, x_upper[i], hk, theta, rt, log_rt, log_rcrit, margin, marker);
    }

    println!("\n* indicates Rθ > Rcrit (amplification active)");

    // Show what happens if we had Rθ growth matching XFOIL
    println!("\n=== Comparison with XFOIL-like Rθ values ===");
    println!("If we had larger Rθ values (as XFOIL typically computes),");
    println!("amplification would start earlier:\n");

    // At x/c = 0.1, XFOIL typically has θ ≈ 0.00025, Rθ ≈ 300 for Re=1e6
    // At x/c = 0.3, XFOIL typically has θ ≈ 0.00055, Rθ ≈ 650 for Re=1e6
    // At x/c = 0.5, XFOIL typically has θ ≈ 0.00085, Rθ ≈ 900 for Re=1e6

    let xfoil_estimates: [(f64, f64, f64, f64); 5] = [
        (0.10, 2.35, 0.00025, 300.0),
        (0.20, 2.45, 0.00040, 500.0),
        (0.30, 2.55, 0.00055, 650.0),
        (0.40, 2.70, 0.00070, 800.0),
        (0.50, 2.90, 0.00085, 920.0),
    ];

    println!("{:>8} {:>8} {:>10} {:>10} {:>12} {:>12}",
        "x/c", "Hk_exp", "θ_exp", "Rθ_exp", "log10(Rcrit)", "Margin");
    for (xc, hk, _th, rt) in xfoil_estimates {
        let log_rt = rt.log10();
        let log_rcrit = rtheta_crit(hk);
        let margin = log_rt - log_rcrit;
        let marker = if margin > 0.0 { "*" } else { " " };
        println!("{:>8.2} {:>8.2} {:>10.5} {:>10.1} {:>12.4} {:>11.4}{}",
            xc, hk, _th, rt, log_rcrit, margin, marker);
    }

    // Show our actual values vs expected
    println!("\n=== Our θ vs Expected XFOIL θ ===");
    println!("{:>8} {:>12} {:>12} {:>10}", "x/c", "θ_ours", "θ_expected", "Ratio");
    let expected = [
        (0.10, 2.5e-4),
        (0.20, 4.0e-4),
        (0.30, 5.5e-4),
        (0.40, 7.0e-4),
        (0.50, 8.5e-4),
    ];
    for (xc_exp, th_exp) in expected {
        // Find closest station
        let mut closest_idx = 0;
        let mut min_diff = f64::INFINITY;
        for (i, &x) in x_upper.iter().enumerate() {
            let diff = (x - xc_exp).abs();
            if diff < min_diff {
                min_diff = diff;
                closest_idx = i;
            }
        }
        if closest_idx < bl_upper.len() {
            let th_our = bl_upper[closest_idx].theta;
            let ratio = th_our / th_exp;
            println!("{:>8.2} {:>12.2e} {:>12.2e} {:>10.3}",
                x_upper[closest_idx], th_our, th_exp, ratio);
        }
    }
}
