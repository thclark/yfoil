use yfoil::geometry::{panel_foil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::find_stagnation_point;

fn main() {
    // Use 80 panels to match debug_viscal_loop
    let geom = naca_4digit("0012", 80).unwrap();
    let airfoil = panel_foil(&geom);

    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    println!("=== DIJ Diagnostic for 80-panel NACA 0012 ===\n");
    println!("Stagnation index: {}", stag_idx);
    println!("Upper panels: 0 to {} ({} panels)", stag_idx - 1, stag_idx);
    println!("Lower panels: {} to {} ({} panels)", stag_idx, airfoil.n - 1, airfoil.n - stag_idx);

    if let Some(dij) = inviscid.get_dij() {
        let n = dij.shape().0;
        println!("\nDIJ matrix size: {}x{}", n, n);

        // Check magnitude
        let max_dij = dij.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_dij = dij.iter().cloned().fold(f64::INFINITY, f64::min);
        println!("Max DIJ: {:.6}", max_dij);
        println!("Min DIJ: {:.6}", min_dij);

        // Row 78 (lower TE area) - check DIJ values
        let row = 78;
        println!("\n=== Row {} (lower surface near TE) ===", row);
        println!("DIJ[{}, 0..5] (upper TE area): {:?}", row,
                 (0..5).map(|j| format!("{:.4}", dij[(row, j)])).collect::<Vec<_>>());
        println!("DIJ[{}, {}..{}] (around stag): {:?}", row, stag_idx-2, stag_idx+3,
                 (stag_idx-2..stag_idx+3).map(|j| format!("{:.4}", dij[(row, j)])).collect::<Vec<_>>());
        println!("DIJ[{}, 75..80] (lower TE area): {:?}", row,
                 (75..80).map(|j| format!("{:.4}", dij[(row, j)])).collect::<Vec<_>>());

        // Row sums separated by surface
        let upper_sum: f64 = (0..stag_idx).map(|j| dij[(row, j)]).sum();
        let lower_sum: f64 = (stag_idx..n).map(|j| dij[(row, j)]).sum();
        println!("\nRow {} sums:", row);
        println!("  Upper surface (0..{}): sum = {:.6}", stag_idx, upper_sum);
        println!("  Lower surface ({}..{}): sum = {:.6}", stag_idx, n, lower_sum);

        // Now simulate the DIJ coupling with typical mass values
        println!("\n=== Simulated DIJ coupling for row {} ===", row);

        // Typical initial mass = dstar * u ≈ 0.001 * 1.0 = 0.001 to 0.01
        let typical_mass = 0.005;  // Typical BL mass defect
        let mut mass = vec![0.0; n];

        // Set mass for all BL panels (excluding stagnation)
        for j in 0..stag_idx {
            if j != stag_idx - 1 {  // Skip stagnation-adjacent
                mass[j] = typical_mass;
            }
        }
        for j in stag_idx+1..n {
            mass[j] = typical_mass;
        }

        // Compute DUI with VTI signs (for row 78, lower surface)
        let vti_i = -1.0;  // Lower surface
        let mut dui = 0.0;
        let mut upper_contrib = 0.0;
        let mut lower_contrib = 0.0;

        for j in 0..n {
            let vti_j = if j < stag_idx { 1.0 } else { -1.0 };
            let contrib = -vti_i * vti_j * dij[(row, j)] * mass[j];
            dui += contrib;
            if j < stag_idx {
                upper_contrib += contrib;
            } else {
                lower_contrib += contrib;
            }
        }

        let qinv_row = qinv[row].abs();
        println!("QINV[{}] = {:.6}", row, qinv_row);
        println!("With typical mass = {:.4} at all stations:", typical_mass);
        println!("  Upper contribution: {:.6e}", upper_contrib);
        println!("  Lower contribution: {:.6e}", lower_contrib);
        println!("  Total DUI: {:.6e}", dui);
        println!("  UNEW = QINV + DUI = {:.6}", qinv_row + dui);

        // Check what mass values would make DUI reasonable
        println!("\n=== Sanity check ===");
        let max_dui = 0.1;  // Reasonable max DUI
        let needed_mass = max_dui / (upper_sum.abs() + lower_sum.abs());
        println!("For max DUI = {:.3}, max mass should be ~{:.6}", max_dui, needed_mass);

        // Also check a row in the middle of the lower surface
        let mid_row = (stag_idx + n) / 2;  // Middle of lower surface
        let upper_sum_mid: f64 = (0..stag_idx).map(|j| dij[(mid_row, j)]).sum();
        let lower_sum_mid: f64 = (stag_idx..n).map(|j| dij[(mid_row, j)]).sum();
        println!("\nRow {} (middle lower surface):", mid_row);
        println!("  Upper surface sum: {:.6}", upper_sum_mid);
        println!("  Lower surface sum: {:.6}", lower_sum_mid);

    } else {
        println!("DIJ matrix not computed!");
    }
}
