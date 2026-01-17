use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit, write_dat_file};
use yfoil::solver::{solve_viscous, calculate_friction_drag, ViscalConfig};

fn main() {
    // Generate geometry
    let geom = naca_4digit("0012", 160).expect("Failed");
    let airfoil = create_paneled_airfoil(&geom);

    // Save geometry
    write_dat_file(&geom, "NACA 0012 from yfoil", "/tmp/yfoil_test_airfoil.dat").expect("write");

    // Solve at alpha=2, Re=1e6
    let alpha_deg: f64 = 2.0;
    let alpha_rad = alpha_deg.to_radians();
    let cond = FlowConditions::new(1e6, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();
    let result = solve_viscous(&airfoil, alpha_rad, &cond, &config);

    println!("=== YFoil Viscous Solution at alpha={} deg, Re={:.0e} ===", alpha_deg, cond.reynolds);
    println!("CL = {:.6}", result.cl);
    println!("CD = {:.6}", result.cd);
    println!("CDf = {:.6}", result.cdf);
    println!("CDp = {:.6}", result.cdp);
    println!("CM = {:.6}", result.cm);
    println!("xtr_upper = {:.4}", result.xtr_upper);
    println!("xtr_lower = {:.4}", result.xtr_lower);
    println!("Converged: {}", result.converged);
    println!("Iterations: {}", result.iterations);
    println!("Residual: {:.6e}", result.residual);

    // BL details from converged solution
    let bl = &result.bl;
    println!("\n=== BL Details (converged) ===");
    println!("Upper surface stations: {}", bl.upper.len());
    println!("Lower surface stations: {}", bl.lower.len());

    // Dump upper surface BL (first 10 and last 10)
    println!("\n--- Upper Surface BL (first 10) ---");
    println!("{:>4} {:>8} {:>8} {:>12} {:>12} {:>8}", "i", "x", "ue", "theta", "cf", "h");
    for (i, r) in bl.upper.iter().take(10).enumerate() {
        let x = bl.x_upper[i];
        println!("{:>4} {:>8.5} {:>8.5} {:>12.4e} {:>12.4e} {:>8.4}", i, x, r.ue, r.theta, r.cf, r.h);
    }
    println!("... {} more stations ...", bl.upper.len().saturating_sub(20));
    println!("--- Upper Surface BL (last 10) ---");
    let start = bl.upper.len().saturating_sub(10);
    for (i, r) in bl.upper.iter().enumerate().skip(start) {
        let x = bl.x_upper[i];
        println!("{:>4} {:>8.5} {:>8.5} {:>12.4e} {:>12.4e} {:>8.4}", i, x, r.ue, r.theta, r.cf, r.h);
    }

    // Recalculate CDf to verify
    let cdf_check = calculate_friction_drag(bl, alpha_rad);
    println!("\nCDf recalculated: {:.6}", cdf_check);
    println!("CDf from result:  {:.6}", result.cdf);

    // Compare with XFOIL expected values (converged)
    println!("\n=== Comparison with XFOIL ===");
    println!("                  yfoil        XFOIL       diff (%)");
    let xfoil_cl = 0.2009;
    let xfoil_cd = 0.00572;
    let xfoil_cdf = 0.00433;
    let xfoil_cdp = 0.00139;

    println!("CL:          {:>10.6}  {:>10.6}  {:>+10.1}%", result.cl, xfoil_cl, 100.0 * (result.cl - xfoil_cl) / xfoil_cl);
    println!("CD:          {:>10.6}  {:>10.6}  {:>+10.1}%", result.cd, xfoil_cd, 100.0 * (result.cd - xfoil_cd) / xfoil_cd);
    println!("CDf:         {:>10.6}  {:>10.6}  {:>+10.1}%", result.cdf, xfoil_cdf, 100.0 * (result.cdf - xfoil_cdf) / xfoil_cdf);
    println!("CDp:         {:>10.6}  {:>10.6}  {:>+10.1}%", result.cdp, xfoil_cdp, 100.0 * (result.cdp - xfoil_cdp) / xfoil_cdp);
}
