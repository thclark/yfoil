//! Validate yfoil viscous solution against XFOIL reference values
//!
//! XFOIL reference (NACA 0012, Re=1e6, M=0, Ncrit=9):
//!   α = 0°:  CL ≈ 0,     CD ≈ 0.0064
//!   α = 5°:  CL ≈ 0.55,  CD ≈ 0.0076
//!   α = 10°: CL ≈ 1.08,  CD ≈ 0.0130

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Validating yfoil viscous solver against XFOIL ===\n");
    println!("NACA 0012, Re=1e6, M=0, Ncrit=9\n");

    let geom = naca_4digit("0012", 80).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    // XFOIL reference values (approximate)
    let cases: [(f64, f64, f64); 3] = [
        (0.0, 0.0, 0.0064),    // α=0°
        (5.0, 0.55, 0.0076),   // α=5°
        (10.0, 1.08, 0.0130),  // α=10°
    ];

    println!("{:>6} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
             "alpha", "CL", "CL_ref", "CL_err%", "CD", "CD_ref", "CD_err%", "Conv");
    println!("{:-<6} {:-<10} {:-<10} {:-<10} {:-<10} {:-<10} {:-<10} {:-<10}",
             "", "", "", "", "", "", "", "");

    for (alpha_deg, cl_ref, cd_ref) in cases {
        let alpha = alpha_deg.to_radians();
        let result = solve_viscous(&airfoil, alpha, &cond, &config);

        let cl_err = if cl_ref.abs() > 0.01 {
            100.0 * (result.cl - cl_ref).abs() / cl_ref.abs()
        } else {
            0.0 // Skip percentage for ~0 values
        };
        let cd_err = 100.0 * (result.cd - cd_ref).abs() / cd_ref;

        println!("{:6.1} {:10.4} {:10.4} {:10.1} {:10.5} {:10.5} {:10.1} {:>10}",
                 alpha_deg, result.cl, cl_ref, cl_err,
                 result.cd, cd_ref, cd_err,
                 if result.converged { "yes" } else { "NO" });

        // Also print component breakdown
        println!("        CDf={:.5} CDp={:.5} xtr_u={:.3} xtr_l={:.3} iter={}",
                 result.cdf, result.cdp, result.xtr_upper, result.xtr_lower, result.iterations);
    }

    println!("\n=== Summary ===");
    println!("Note: XFOIL reference values are approximate.");
    println!("Expected tolerances: CL within 10%, CD within 30%");
    println!("(Drag is more sensitive to coupling convergence and wake modeling)");
}
