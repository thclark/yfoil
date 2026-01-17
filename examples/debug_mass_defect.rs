//! Debug mass defect magnitudes

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, read_dat_file};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("=== Debug Mass Defect Values ===\n");

    let (_, geom) =
        read_dat_file("/tmp/xfoil_naca0012_paneled.dat").expect("Run XFOIL first");

    let airfoil = create_paneled_airfoil(&geom);
    let n = airfoil.n;

    println!("Airfoil: {} panels", n);

    // Flow conditions
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    // Run solver with minimal iterations to get intermediate values
    let config = ViscalConfig {
        max_iter: 5,
        tol_cl: 1e-5,
        relax: 0.5,
        ..Default::default()
    };

    let result = solve_viscous(&airfoil, 0.0, &cond, &config);

    // Check dq_source values
    println!("\ndq_source statistics:");
    let mut sum = 0.0;
    let mut max_abs: f64 = 0.0;
    let mut max_idx = 0;
    for (i, &dq) in result.dq_source.iter().enumerate() {
        sum += dq.abs();
        if dq.abs() > max_abs {
            max_abs = dq.abs();
            max_idx = i;
        }
    }
    println!("Mean |dq|: {:.6}", sum / n as f64);
    println!("Max |dq|: {:.6} at index {}", max_abs, max_idx);

    println!("\nFirst 10 dq_source values:");
    for i in 0..10.min(n) {
        println!("dq[{}] = {:.6}", i, result.dq_source[i]);
    }

    println!("\nLE region dq_source values:");
    let le = airfoil.le_index;
    for i in (le - 3)..(le + 4).min(n) {
        println!("dq[{}] = {:.6}", i, result.dq_source[i]);
    }
}
