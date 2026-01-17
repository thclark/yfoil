//! Compare yfoil polar with XFOIL for NACA 0012
//!
//! Run with: cargo run --example polar_comparison --release

use std::fs::File;
use std::io::Write;
use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit, write_dat_file};
use yfoil::solver::{solve_viscous_with_init, ViscalConfig, ViscousResult};

fn main() {
    println!("=== NACA 0012 Polar Comparison: yfoil vs XFOIL ===\n");

    // Generate NACA 0012 with 160 panels (same as XFOIL default)
    let n_panels = 160;
    let geom = naca_4digit("0012", n_panels).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    println!("Airfoil: NACA 0012");
    println!("Panels: {}", airfoil.n);
    println!("Sharp TE: {}", airfoil.sharp_te);

    // Save the geometry for XFOIL
    write_dat_file(&geom, "NACA 0012", "/tmp/naca0012_yfoil.dat").expect("Failed to write dat file");
    println!("Saved geometry to /tmp/naca0012_yfoil.dat\n");

    // Flow conditions
    let re = 1_000_000.0;
    let mach = 0.0;
    let ncrit = 9.0;
    let cond = FlowConditions::new(re, mach, ncrit, 1.0);

    println!("Flow conditions:");
    println!("  Re = {:.0}", re);
    println!("  Mach = {:.2}", mach);
    println!("  Ncrit = {:.1}", ncrit);
    println!();

    // Configure solver
    let mut config = ViscalConfig::default();
    config.max_iter = 100;

    // Polar sweep: 0 to 15, then -1 to -15
    let alphas_up: Vec<f64> = (0..=15).map(|i| i as f64).collect();
    let alphas_down: Vec<f64> = (-15..0).map(|i| i as f64).rev().collect();

    let mut results: Vec<ViscousResult> = Vec::new();

    println!("Running yfoil polar sweep...");
    println!("  Sweeping 0° to 15°...");

    // Sweep up from 0
    let mut prev_dq: Option<Vec<f64>> = None;
    for &alpha_deg in &alphas_up {
        let alpha_rad = alpha_deg.to_radians();
        let result = solve_viscous_with_init(&airfoil, alpha_rad, &cond, &config, prev_dq.as_deref());

        if result.converged {
            prev_dq = Some(result.dq_source.clone());
        }
        results.push(result);
    }

    println!("  Sweeping -1° to -15°...");

    // Sweep down from -1 (use solution from alpha=0 as starting point)
    prev_dq = results.first().map(|r| r.dq_source.clone());
    for &alpha_deg in &alphas_down {
        let alpha_rad = alpha_deg.to_radians();
        let result = solve_viscous_with_init(&airfoil, alpha_rad, &cond, &config, prev_dq.as_deref());

        if result.converged {
            prev_dq = Some(result.dq_source.clone());
        }
        results.push(result);
    }

    // Sort results by alpha for output
    results.sort_by(|a, b| a.alpha.partial_cmp(&b.alpha).unwrap());

    // Write yfoil results to file
    let mut file = File::create("/tmp/yfoil_polar.dat").expect("Failed to create output file");
    writeln!(file, "# yfoil NACA 0012 polar - Re={:.0}, M={:.2}, Ncrit={:.1}", re, mach, ncrit).unwrap();
    writeln!(file, "# alpha      CL         CD         CDf        CDp        CM       conv  iter").unwrap();
    for r in &results {
        writeln!(file, "{:7.2} {:10.6} {:10.6} {:10.6} {:10.6} {:10.6} {:5} {:4}",
            r.alpha.to_degrees(), r.cl, r.cd, r.cdf, r.cdp, r.cm,
            if r.converged { "Y" } else { "N" }, r.iterations
        ).unwrap();
    }
    println!("Saved yfoil polar to /tmp/yfoil_polar.dat\n");

    // Print yfoil results
    println!("yfoil Results:");
    println!("{:>7} {:>10} {:>10} {:>10} {:>10} {:>6} {:>4}",
        "alpha", "CL", "CD", "CDf", "CDp", "conv", "iter");
    println!("{}", "-".repeat(65));
    for r in &results {
        println!("{:7.2} {:10.6} {:10.6} {:10.6} {:10.6} {:>6} {:4}",
            r.alpha.to_degrees(), r.cl, r.cd, r.cdf, r.cdp,
            if r.converged { "Y" } else { "N" }, r.iterations);
    }
    println!();

    // Generate XFOIL script
    let xfoil_script = format!(r#"PLOP
G F

NACA 0012
PANE
PPAR
N {}


OPER
VISC {}
ITER 100
PACC
/tmp/xfoil_polar.dat

ASEQ 0 15 1
ASEQ -1 -15 -1

QUIT
"#, n_panels, re as i64);

    std::fs::write("/tmp/xfoil_polar.inp", &xfoil_script).expect("Failed to write XFOIL script");
    println!("Generated XFOIL script at /tmp/xfoil_polar.inp");
    println!("Run XFOIL with: cd /tmp && xfoil < xfoil_polar.inp");
}
