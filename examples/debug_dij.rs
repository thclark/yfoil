use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let inviscid = solve_inviscid(&airfoil);
    
    println!("Has DIJ: {}", inviscid.has_dij());
    
    if let Some(dij) = inviscid.get_dij() {
        let (rows, cols) = dij.shape();
        println!("DIJ matrix size: {}x{}", rows, cols);
        
        // Check some DIJ values
        println!("\nDIJ[0, 0..5]: {:?}", (0..5).map(|j| dij[(0, j)]).collect::<Vec<_>>());
        println!("DIJ[0, {}..{}]: {:?}", cols-5, cols, (cols-5..cols).map(|j| dij[(0, j)]).collect::<Vec<_>>());
        
        // Sum of each row (for debugging)
        println!("\nRow sums (first 5):");
        for i in 0..5 {
            let sum: f64 = (0..cols).map(|j| dij[(i, j)]).sum();
            println!("  Row {}: sum = {:.6}", i, sum);
        }
        
        // Diagonal values
        println!("\nDiagonal values:");
        for i in 0..5 {
            println!("  DIJ[{},{}] = {:.6}", i, i, dij[(i, i)]);
        }
        
        // Check magnitude
        let max_dij = dij.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_dij = dij.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("\nMax DIJ: {:.6}", max_dij);
        println!("Min DIJ: {:.6}", min_dij);
        
        // Test: compute dq for unit mass defect at panel 10
        let mut mass = vec![0.0; cols];
        mass[10] = 0.001; // Small mass defect at panel 10
        
        let mut dq_test = vec![0.0; rows];
        for i in 0..rows {
            for j in 0..cols {
                dq_test[i] += -dij[(i, j)] * mass[j]; // Note: XFOIL uses negative
            }
        }
        
        println!("\nTest dq for mass[10] = 0.001:");
        println!("  dq[8] = {:.6e}", dq_test[8]);
        println!("  dq[9] = {:.6e}", dq_test[9]);
        println!("  dq[10] = {:.6e}", dq_test[10]);
        println!("  dq[11] = {:.6e}", dq_test[11]);
        println!("  dq[12] = {:.6e}", dq_test[12]);
    } else {
        println!("DIJ matrix not computed!");
    }
}
