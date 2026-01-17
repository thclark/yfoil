//! Investigate why BL theta is 2x too high
//!
//! Compare BL development along the surface with expected values

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Investigating BL theta discrepancy ===\n");

    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let mut config = ViscalConfig::default();
    config.tol_cl = 5e-4;

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    println!("Converged: {} in {} iterations\n", result.converged, result.iterations);

    // Print full upper surface BL development
    println!("=== Upper Surface BL Development ===");
    println!("{:>4} {:>8} {:>10} {:>10} {:>8} {:>8} {:>10} {:>8}",
             "i", "x/c", "s", "theta", "dstar", "H", "Cf", "Ue");
    println!("{}", "-".repeat(78));

    let s_total_upper = result.bl.s_upper.last().copied().unwrap_or(1.0);
    for (i, st) in result.bl.upper.iter().enumerate() {
        let s = result.bl.s_upper.get(i).copied().unwrap_or(0.0);
        let x = result.bl.x_upper.get(i).copied().unwrap_or(0.0);
        // Approximate x/c (assuming chord = 1)
        let xc = x;
        println!("{:4} {:8.4} {:10.6} {:10.6} {:10.6} {:8.3} {:10.6} {:8.4}",
                 i, xc, s, st.theta, st.dstar, st.h, st.cf, st.ue);
    }

    println!("\n=== Lower Surface BL Development ===");
    println!("{:>4} {:>8} {:>10} {:>10} {:>10} {:>8} {:>10} {:>8}",
             "i", "x/c", "s", "theta", "dstar", "H", "Cf", "Ue");
    println!("{}", "-".repeat(78));

    for (i, st) in result.bl.lower.iter().enumerate() {
        let s = result.bl.s_lower.get(i).copied().unwrap_or(0.0);
        let x = result.bl.x_lower.get(i).copied().unwrap_or(0.0);
        let xc = x;
        println!("{:4} {:8.4} {:10.6} {:10.6} {:10.6} {:8.3} {:10.6} {:8.4}",
                 i, xc, s, st.theta, st.dstar, st.h, st.cf, st.ue);
    }

    // Analysis
    println!("\n=== Analysis ===");

    // Check Blasius solution for reference
    // For flat plate at Re_x = 1e6:
    // theta/x = 0.664 / sqrt(Re_x) = 0.664 / 1000 = 0.000664
    // So at x = 0.5 (mid-chord), theta_Blasius = 0.000332
    println!("\nBlasius reference (flat plate, Re=1e6):");
    println!("  theta/x = 0.664/sqrt(Re_x)");
    println!("  At x=0.5: theta = 0.000332");
    println!("  At x=1.0: theta = 0.000664");

    // Check actual values at mid-chord
    let mid_idx_upper = result.bl.upper.len() / 2;
    if mid_idx_upper < result.bl.upper.len() {
        let st = &result.bl.upper[mid_idx_upper];
        let x = result.bl.x_upper.get(mid_idx_upper).copied().unwrap_or(0.5);
        println!("\nyfoil upper at x≈{:.2}: theta = {:.6}", x, st.theta);

        // Expected theta for laminar BL
        let re_x = 1e6 * x;
        let theta_blasius = 0.664 * x / re_x.sqrt();
        println!("  Blasius at x={:.2}: theta = {:.6}", x, theta_blasius);
        println!("  Ratio yfoil/Blasius = {:.2}", st.theta / theta_blasius);
    }

    // Check TE values
    println!("\n=== TE Values ===");
    if let Some(u) = result.bl.upper.last() {
        println!("Upper TE: theta={:.6}, H={:.3}", u.theta, u.h);
    }
    if let Some(l) = result.bl.lower.last() {
        println!("Lower TE: theta={:.6}, H={:.3}", l.theta, l.h);
    }

    // Expected TE theta for NACA 0012 at Re=1e6
    // XFOIL typically gives theta_upper ≈ 0.0015, theta_lower ≈ 0.0015
    // Combined ≈ 0.003
    println!("\nExpected (XFOIL reference):");
    println!("  Upper TE theta ≈ 0.0015");
    println!("  Lower TE theta ≈ 0.0015");
    println!("  Combined ≈ 0.003");

    // Check transition location
    println!("\n=== Transition ===");
    println!("xtr_upper = {:.3}", result.xtr_upper);
    println!("xtr_lower = {:.3}", result.xtr_lower);

    // Find where H drops below 1.6 (turbulent)
    for (i, st) in result.bl.upper.iter().enumerate() {
        if st.h < 1.6 && i > 0 && result.bl.upper[i-1].h >= 1.6 {
            let x = result.bl.x_upper.get(i).copied().unwrap_or(0.0);
            println!("Upper: H drops below 1.6 at station {}, x={:.4}", i, x);
            break;
        }
    }
    for (i, st) in result.bl.lower.iter().enumerate() {
        if st.h < 1.6 && i > 0 && result.bl.lower[i-1].h >= 1.6 {
            let x = result.bl.x_lower.get(i).copied().unwrap_or(0.0);
            println!("Lower: H drops below 1.6 at station {}, x={:.4}", i, x);
            break;
        }
    }
}
