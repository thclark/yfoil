use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();
    
    let result = solve_viscous(&airfoil, 0.0, &cond, &config);
    
    println!("Converged: {}, Iterations: {}", result.converged, result.iterations);
    println!("CL = {:.6}, CD = {:.6}", result.cl, result.cd);
    
    // Check dq_source at TE
    println!("\n=== Source velocity correction (dq_source) ===");
    println!("TE upper (idx 0): {:.6}", result.dq_source[0]);
    println!("TE lower (idx {}): {:.6}", result.dq_source.len()-1, result.dq_source[result.dq_source.len()-1]);
    
    // Near LE
    let le = 80;
    println!("Near LE (idx {}): {:.6}", le, result.dq_source[le.min(result.dq_source.len()-1)]);
    
    // Summary
    let max_dq = result.dq_source.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_dq = result.dq_source.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("\nMax dq_source: {:.6}", max_dq);
    println!("Min dq_source: {:.6}", min_dq);
    
    // Expected: to get from Qinv=0.755 to Ue=0.886, need dq = +0.13
    println!("\nExpected dq at TE to match XFOIL: +0.13");
    
    // Check the upper surface final Ue
    if !result.bl.upper.is_empty() {
        let upper_te = result.bl.upper.last().unwrap();
        println!("\nUpper TE Ue from BL: {:.6}", upper_te.ue);
    }
}
