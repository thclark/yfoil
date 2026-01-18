//! Test YFoil using XFOIL's exact panel coordinates
//!
//! This isolates whether the BL solver works correctly when given
//! identical geometry to XFOIL.

use std::fs::File;
use std::io::{BufRead, BufReader};

use yfoil::bl::{FlowConditions, NewtonConfig, march_newton};
use yfoil::geometry::{create_paneled_airfoil, Geometry};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, extract_lower_surface, find_stagnation_point};

fn main() {
    // Read XFOIL panel coordinates
    let file = File::open("/tmp/xfoil_panels.dat").expect("Failed to open XFOIL panels file");
    let reader = BufReader::new(file);
    
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    
    // Parse header to get N
    let n: usize = lines.iter()
        .find(|l| l.contains("N ="))
        .and_then(|l| l.split_whitespace().last())
        .and_then(|s| s.parse().ok())
        .expect("Failed to parse N");
    
    println!("=== Reading XFOIL Panel Geometry ===");
    println!("N panels: {}", n);
    
    // Find start of panel data
    let data_start = lines.iter().position(|l| l.contains("I, X, Y")).unwrap() + 1;
    
    // Parse panel coordinates
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    
    for line in &lines[data_start..data_start + n] {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let xi: f64 = parts[1].parse().expect("Failed to parse x");
            let yi: f64 = parts[2].parse().expect("Failed to parse y");
            x.push(xi);
            y.push(yi);
        }
    }
    
    println!("Read {} panel coordinates", x.len());
    println!("TE upper: x={:.10}, y={:.10}", x[0], y[0]);
    println!("TE lower: x={:.10}, y={:.10}", x[n-1], y[n-1]);
    
    // Create Geometry from XFOIL coordinates
    let geom = Geometry {
        reference: [0.25, 0.0],
        x_c: x,
        y_c: y,
    };
    
    // Create PaneledAirfoil
    let airfoil = create_paneled_airfoil(&geom);
    
    println!("\n=== Inviscid Solution ===");
    let inviscid = solve_inviscid(&airfoil);
    let alpha = 0.0_f64;
    let qinv = inviscid.velocity_at_nodes(alpha);
    
    println!("qinv[0] (upper TE): {:.10}", qinv[0]);
    println!("qinv[{}] (lower TE): {:.10}", n-1, qinv[n-1]);
    println!("XFOIL expected: ±0.7670558229");
    
    // Find stagnation point
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    println!("\nStagnation point: index {}, x={:.6}, y={:.6}",
             stag_idx, airfoil.x[stag_idx], airfoil.y[stag_idx]);
    
    // Extract surfaces
    let (x_upper, _y_upper, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    let (_x_lower, _y_lower, s_lower, ue_lower) = extract_lower_surface(&airfoil, &qinv, stag_idx);
    
    println!("\nUpper surface: {} stations", ue_upper.len());
    println!("Lower surface: {} stations", ue_lower.len());
    
    // Print edge velocities at TE
    println!("\nEdge velocity at TE:");
    println!("Upper: x={:.6}, ue={:.10}", x_upper.last().unwrap(), ue_upper.last().unwrap());
    println!("Lower: ue={:.10}", ue_lower.last().unwrap());
    
    // Flow conditions
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    
    // March BL
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &config);
    
    // Final TE values
    let upper_te = bl_upper.last().unwrap();
    let lower_te = bl_lower.last().unwrap();
    
    println!("\n=== BL Results ===");
    println!("Upper TE: theta={:.8e}, dstar={:.8e}, Ue={:.8e}",
             upper_te.theta, upper_te.dstar, upper_te.ue);
    println!("Lower TE: theta={:.8e}, dstar={:.8e}, Ue={:.8e}",
             lower_te.theta, lower_te.dstar, lower_te.ue);
    
    let tte = upper_te.theta + lower_te.theta;
    let dte = upper_te.dstar + lower_te.dstar;
    println!("\nCombined TTE: {:.8e}", tte);
    println!("Combined DTE: {:.8e}", dte);
    
    println!("\n=== XFOIL Expected ===");
    println!("XFOIL TTE: 3.9961520077e-3");
    println!("XFOIL DTE: 8.8492753897e-3 (includes ANTE)");
    
    let xfoil_tte = 3.9961520077e-3;
    let error = (tte - xfoil_tte).abs() / xfoil_tte * 100.0;
    println!("\nTTE error: {:.1}%", error);
}
