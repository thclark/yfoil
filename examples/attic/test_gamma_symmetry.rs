//! Debug test for gamma symmetry at LE

use yfoil::geometry::{naca_4digit, repanel_by_curvature, panel_foil, PaneConfig};
use yfoil::panel::solve_inviscid;

/// Find stagnation point (copy of the logic from viscal.rs)
fn find_stagnation_point(airfoil: &yfoil::geometry::PanelledFoil, velocity: &[f64]) -> usize {
    let le_idx = airfoil.le_index;
    let search_range = (airfoil.n / 4).max(5);

    let start = le_idx.saturating_sub(search_range);
    let end = (le_idx + search_range).min(airfoil.n);

    // Look for sign change in velocity
    for i in start..end.saturating_sub(1) {
        if velocity[i] > 0.0 && velocity[i + 1] <= 0.0 {
            return i + 1;
        }
    }

    // Fallback: find minimum velocity magnitude
    let mut min_vel = f64::MAX;
    let mut stag_idx = le_idx;

    for i in start..end {
        let vel_mag = velocity[i].abs();
        if vel_mag < min_vel {
            min_vel = vel_mag;
            stag_idx = i;
        }
    }

    stag_idx
}

fn main() {
    // Create NACA 0012 with XFOIL-style paneling
    let geom = naca_4digit("0012", 160).expect("Failed to create airfoil");

    // Check raw geometry before repaneling
    println!("=== Raw NACA geometry near LE (before repaneling) ===");
    let le_idx = geom.x_c.iter().enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(geom.x_c.len() / 2);
    println!("LE index in raw geometry: {}", le_idx);
    for i in le_idx.saturating_sub(3)..=(le_idx + 3).min(geom.x_c.len() - 1) {
        let marker = if i == le_idx { " <-- LE" } else { "" };
        println!("  Point[{}]: x={:12.9}, y={:12.9}{}",
                 i, geom.x_c[i], geom.y_c[i], marker);
    }

    let config = PaneConfig::default();
    let repaneled = repanel_by_curvature(&geom, 160, &config);
    let airfoil = panel_foil(&repaneled);

    println!("Airfoil: {} panels, LE at index {}", airfoil.n, airfoil.le_index);

    // Solve inviscid at alpha=0
    let inviscid = solve_inviscid(&airfoil);
    let gamma = inviscid.velocity_at_nodes(0.0);

    // Find stagnation point
    let stag_idx = find_stagnation_point(&airfoil, &gamma);
    println!("Stagnation index: {} (0-based)", stag_idx);

    // Print gamma values near stagnation
    println!("\nGamma values near stagnation:");
    let start = stag_idx.saturating_sub(5);
    let end = (stag_idx + 6).min(airfoil.n);
    for i in start..end {
        let marker = if i == stag_idx { " <-- stag_idx" }
                     else if i == stag_idx - 1 { " <-- upper[1]" }
                     else { "" };
        println!("  panel[{}]: gamma={:+.10}, x={:.6}, y={:.6}{}",
                 i, gamma[i], airfoil.x[i], airfoil.y[i], marker);
    }

    // Check symmetry: panels stag_idx-1 and stag_idx should have
    // nearly equal magnitudes for symmetric airfoil at alpha=0
    let upper_gamma = gamma[stag_idx - 1];
    let lower_gamma = gamma[stag_idx];
    println!("\nSymmetry check:");
    println!("  Upper (stag_idx-1={}): gamma={:+.10}, |gamma|={:.10}",
             stag_idx - 1, upper_gamma, upper_gamma.abs());
    println!("  Lower (stag_idx={}):   gamma={:+.10}, |gamma|={:.10}",
             stag_idx, lower_gamma, lower_gamma.abs());
    println!("  Ratio of magnitudes: {:.6}",
             upper_gamma.abs() / lower_gamma.abs().max(1e-10));

    // What XFOIL expects
    println!("\nXFOIL expects UINV(IBL=2) = 0.1461460999 on both surfaces");
    println!("YFoil computed:");
    println!("  Upper: |gamma|={:.10}", upper_gamma.abs());
    println!("  Lower: |gamma|={:.10}", lower_gamma.abs());

    // Check geometry near LE - compare with XFOIL's coordinates
    println!("\n=== Geometry near LE ===");
    println!("XFOIL panel coordinates for NACA 0012, 160 panels:");
    println!("  Panel 80 (1-based): should be symmetric upper LE node");
    println!("  Panel 81 (1-based): should be symmetric lower LE node");
    println!("\nYFoil panel coordinates (0-based):");
    for i in (stag_idx.saturating_sub(3))..=(stag_idx + 3).min(airfoil.n - 1) {
        let marker = if i == stag_idx { " <-- stag_idx" }
                     else if i == stag_idx - 1 { " <-- upper[1]" }
                     else { "" };
        println!("  Panel[{}] (1-based {}): x={:12.9}, y={:12.9}{}",
                 i, i + 1, airfoil.x[i], airfoil.y[i], marker);
    }

    // Compute arc lengths from LE
    let le_arc = airfoil.sle;
    println!("\n=== Arc lengths from LE ===");
    println!("LE arc length (sle) = {:.10}", le_arc);
    for i in (stag_idx.saturating_sub(3))..=(stag_idx + 3).min(airfoil.n - 1) {
        let arc_from_le = (airfoil.s[i] - le_arc).abs();
        let marker = if i == stag_idx { " <-- stag_idx" }
                     else if i == stag_idx - 1 { " <-- upper[1]" }
                     else { "" };
        println!("  Panel[{}]: s={:.10}, |s-sle|={:.10}{}",
                 i, airfoil.s[i], arc_from_le, marker);
    }
}
