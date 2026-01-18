//! Debug BL initialization and first few stations

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
    
    let n: usize = lines.iter()
        .find(|l| l.contains("N ="))
        .and_then(|l| l.split_whitespace().last())
        .and_then(|s| s.parse().ok())
        .expect("Failed to parse N");
    
    let data_start = lines.iter().position(|l| l.contains("I, X, Y")).unwrap() + 1;
    
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
    
    let geom = Geometry {
        reference: [0.25, 0.0],
        x_c: x,
        y_c: y,
    };
    
    let airfoil = create_paneled_airfoil(&geom);
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    
    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    
    println!("=== Upper Surface First 10 Stations ===");
    println!("{:>4} {:>12} {:>14} {:>14}", "Idx", "x/c", "s", "Ue");
    for i in 0..10.min(ue_upper.len()) {
        println!("{:>4} {:>12.8} {:>14.8e} {:>14.8e}", i, x_upper[i], s_upper[i], ue_upper[i]);
    }
    
    // Show what march_newton will compute for initialization
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    println!("\n=== BL Initialization Debug ===");
    println!("cond.nu = {:e}", cond.nu);
    println!("s_dist[0] = {:e}", s_upper[0]);
    println!("s_dist[1] = {:e}", s_upper[1]);
    println!("ue_dist[0] = {:e}", ue_upper[0]);
    println!("ue_dist[1] = {:e}", ue_upper[1]);
    
    // Thwaites formula at station 1
    let s_init = s_upper[1].max(1e-10);
    let ue_init = ue_upper[1].max(0.01);
    let theta_init = (0.45 * cond.nu * s_init / (6.0 * ue_init)).sqrt();
    println!("\nThwaites initialization:");
    println!("theta_init = sqrt(0.45 * {} * {} / (6 * {})) = {:e}",
             cond.nu, s_init, ue_init, theta_init);
    
    // March BL
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);
    
    println!("\n=== BL March First 10 Stations ===");
    println!("{:>4} {:>14} {:>14} {:>10} {:>10}", "Idx", "theta", "dstar", "Hk", "ue");
    for (i, r) in bl_upper.iter().take(10).enumerate() {
        println!("{:>4} {:>14.8e} {:>14.8e} {:>10.4} {:>10.6}", 
                 i, r.theta, r.dstar, r.hk, r.ue);
    }
    
    // Last 5 stations
    println!("\n=== BL March Last 5 Stations ===");
    let start = bl_upper.len().saturating_sub(5);
    for (i, r) in bl_upper.iter().enumerate().skip(start) {
        println!("{:>4} {:>14.8e} {:>14.8e} {:>10.4} {:>10.6}", 
                 i, r.theta, r.dstar, r.hk, r.ue);
    }
    
    let te = bl_upper.last().unwrap();
    println!("\n=== Final TE Values ===");
    println!("theta = {:e}", te.theta);
    println!("dstar = {:e}", te.dstar);
    println!("Hk = {:.4}", te.hk);
    println!("XFOIL expected theta ≈ 2.0e-3 per surface");
}
