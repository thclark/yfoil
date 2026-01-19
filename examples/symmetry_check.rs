use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::forces::integrate_forces;
use yfoil::bl::FlowConditions;
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    println!("=== Panel Setup ===");
    println!("Total nodes: {}", airfoil.n);
    println!("LE index: {}", airfoil.le_index);

    let le = airfoil.le_index;
    let n = airfoil.n;

    // Check if we have equal panels on each side
    let upper_panels = le; // panels 0 to le-1
    let lower_panels = n - le - 1; // panels le to n-2
    println!("Upper surface panels: {}", upper_panels);
    println!("Lower surface panels: {}", lower_panels);

    // Solve inviscid
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let gam = solution.gamma_at_alpha(0.0);

    // Check gamma (node values) at TE
    println!("\n=== Gamma (node values) at TE ===");
    println!("gam[0] = {:+.6} (upper TE node)", gam[0]);
    println!("gam[1] = {:+.6} (first upper interior)", gam[1]);
    println!("gam[{}] = {:+.6} (last lower interior)", n-2, gam[n-2]);
    println!("gam[{}] = {:+.6} (lower TE node)", n-1, gam[n-1]);
    println!("Kutta check: gam[0] + gam[{}] = {:.2e}", n-1, gam[0] + gam[n-1]);

    // Check velocity symmetry
    println!("\n=== Panel Velocity at α=0° (midpoints) ===");
    println!("Near TE:");
    for i in 0..5 {
        println!("  vel[{}]={:+.6}", i, vel[i]);
    }
    println!("  ...");
    for i in (n-5)..n {
        println!("  vel[{}]={:+.6}", i, vel[i]);
    }

    // Check forces
    let coeffs = integrate_forces(&airfoil, &vel, 0.0, 0.0);

    println!("\n=== Inviscid Forces at α=0° ===");
    println!("CL = {:+.6} (should be 0)", coeffs.cl);
    println!("CM = {:+.6}", coeffs.cm);
    println!("CDp = {:+.6}", coeffs.cdp);

    // Now check viscous
    println!("\n=== Viscous Solution at α=0° ===");
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();
    let result = solve_viscous(&airfoil, 0.0, &cond, &config);

    println!("CL = {:+.6} (should be 0)", result.cl);
    println!("CD = {:+.6}", result.cd);
    println!("Converged: {}, Iterations: {}", result.converged, result.iterations);

    // Check stagnation point
    use yfoil::solver::find_stagnation_point;
    let stag = find_stagnation_point(&airfoil, &vel);
    println!("\n=== Stagnation Point ===");
    println!("Stagnation index: {} (LE index: {})", stag, le);
    println!("Stagnation location: ({:.6}, {:.6})", airfoil.x[stag], airfoil.y[stag]);
    println!("Velocity at stag: {:.6}", vel[stag]);

    // Velocity near LE
    println!("\nVelocity near LE:");
    for i in (le-3)..=(le+3) {
        if i < n {
            println!("  vel[{}] = {:+.10} (|vel| = {:.10})", i, vel[i], vel[i].abs());
        }
    }

    // Check dq_source symmetry
    println!("\n=== Source Velocity Correction ===");
    for i in 0..5 {
        println!("  dq[{}]={:+.2e}, dq[{}]={:+.2e}, sum={:.2e}",
            i, result.dq_source[i],
            n-1-i, result.dq_source[n-1-i],
            result.dq_source[i] + result.dq_source[n-1-i]);
    }

    // Check TE velocity used for BL
    println!("\n=== TE Velocity for BL ===");
    println!("Upper TE (vel[0])  = {:.6} (panel 0->1)", vel[0]);
    println!("Lower TE (vel[159]) = {:.6} (panel 159->160)", vel[159]);
    println!("Lower TE (vel[160]) = {:.6} (TE gap panel 160->0)", vel[160]);

    // Check BL edge velocities that would be extracted
    use yfoil::solver::{extract_upper_surface, extract_lower_surface, solve_boundary_layer};
    use yfoil::bl::{NewtonConfig, WakeConfig};
    let stag = find_stagnation_point(&airfoil, &vel);
    let (_, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let (_, _, s_lower, ue_lower) = extract_lower_surface(&airfoil, &vel, stag);

    println!("\n=== Extracted BL Edge Velocities ===");
    println!("Upper surface: {} stations, s_max={:.6}", ue_upper.len(), s_upper.last().unwrap_or(&0.0));
    println!("Lower surface: {} stations, s_max={:.6}", ue_lower.len(), s_lower.last().unwrap_or(&0.0));

    println!("\nFirst 5 stations (from LE):");
    for i in 0..5.min(ue_upper.len()) {
        println!("  Upper[{}]: ue={:.6}, Lower[{}]: ue={:.6}, diff={:.2e}",
            i, ue_upper[i], i, ue_lower[i], ue_upper[i] - ue_lower[i]);
    }

    println!("\nLast 5 stations (near TE):");
    let nu = ue_upper.len();
    let nl = ue_lower.len();
    for i in 0..5.min(nu) {
        let ui = nu - 1 - i;
        let li = nl - 1 - i;
        println!("  Upper[{}]: ue={:.6}, Lower[{}]: ue={:.6}, diff={:.2e}",
            ui, ue_upper[ui], li, ue_lower[li], ue_upper[ui] - ue_lower[li]);
    }

    // Check arc length symmetry
    println!("\n=== Arc Length Symmetry ===");
    for i in 0..5.min(nu) {
        println!("  s_upper[{}]={:.6}, s_lower[{}]={:.6}, diff={:.2e}",
            i, s_upper[i], i, s_lower[i], s_upper[i] - s_lower[i]);
    }
    println!("  s_upper[TE]={:.6}, s_lower[TE]={:.6}, diff={:.2e}",
        s_upper[nu-1], s_lower[nl-1], s_upper[nu-1] - s_lower[nl-1]);

    // Solve BL and check results symmetry - using march_newton directly to track n_amp
    let newton_config = NewtonConfig::default();
    let wake_config = WakeConfig::default();
    let gamma_total = 0.0; // At alpha=0 for symmetric airfoil

    // March BL directly to see n_amp evolution
    use yfoil::bl::march_newton;
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &newton_config);
    let bl_lower = march_newton(&ue_lower, &s_lower, &cond, &newton_config);

    println!("\n=== Amplification Factor Evolution ===");
    println!("Stations around transition (n_amp approaching ncrit={}):", cond.ncrit);
    println!("{:>5} {:>12} {:>12} {:>12} {:>10} {:>10} {:>10} {:>10} {:>6} {:>8}",
             "Stn", "n_upper", "n_lower", "n_diff", "Hk_u", "Hk_l", "θ_u", "θ_l", "iter", "resid");
    for j in 45..70.min(bl_upper.len()) {
        if j < bl_lower.len() {
            let u = &bl_upper[j];
            let l = &bl_lower[j];
            let marker = if u.n_amp >= cond.ncrit || l.n_amp >= cond.ncrit { " <--" } else { "" };
            println!("{:>5} {:>12.6} {:>12.6} {:>12.2e} {:>10.4} {:>10.4} {:>10.2e} {:>10.2e} {:>6} {:>8.1e}{}",
                     j, u.n_amp, l.n_amp, u.n_amp - l.n_amp,
                     u.hk, l.hk, u.theta, l.theta, u.iterations, u.residual, marker);
        }
    }

    // Check Rtheta values around divergence point
    println!("\n=== Rtheta Evolution (controls amplification onset) ===");
    println!("{:>5} {:>12} {:>12} {:>12}", "Stn", "Rt_upper", "Rt_lower", "Rt_diff");
    for j in 45..70.min(bl_upper.len()) {
        if j < bl_lower.len() {
            let rt_u = ue_upper[j] * bl_upper[j].theta / cond.nu;
            let rt_l = ue_lower[j] * bl_lower[j].theta / cond.nu;
            println!("{:>5} {:>12.2} {:>12.2} {:>12.2e}", j, rt_u, rt_l, rt_u - rt_l);
        }
    }

    // Also get full BL solution for rest of analysis
    let bl = solve_boundary_layer(&airfoil, &vel, stag, 0.0, gamma_total, &cond, &newton_config, &wake_config);

    println!("\n=== BL Solution Symmetry ===");
    println!("Upper BL stations: {}", bl.upper.len());
    println!("Lower BL stations: {}", bl.lower.len());

    if !bl.upper.is_empty() && !bl.lower.is_empty() {
        println!("\nBL evolution (every 10 stations):");
        println!("{:>5} {:>12} {:>12} {:>12} {:>12} {:>10}",
            "Stn", "θ_upper", "θ_lower", "θ_diff", "H_upper", "H_lower");
        let step = 10;
        for j in (0..bl.upper.len()).step_by(step) {
            if j < bl.lower.len() {
                let u = &bl.upper[j];
                let l = &bl.lower[j];
                println!("{:>5} {:>12.4e} {:>12.4e} {:>12.2e} {:>10.4} {:>10.4}",
                    j, u.theta, l.theta, u.theta - l.theta, u.h, l.h);
            }
        }
        // Always show last station
        let u_te = bl.upper.last().unwrap();
        let l_te = bl.lower.last().unwrap();
        println!("{:>5} {:>12.4e} {:>12.4e} {:>12.2e} {:>10.4} {:>10.4}",
            bl.upper.len()-1, u_te.theta, l_te.theta,
            u_te.theta - l_te.theta, u_te.h, l_te.h);

        // Check edge velocity profile around problematic region
        println!("\nEdge velocity around station 50-60:");
        for j in 48..65.min(ue_upper.len()) {
            println!("  Station {}: ue_upper={:.4}, ue_lower={:.4}, diff={:.2e}",
                j, ue_upper[j], ue_lower[j], ue_upper[j] - ue_lower[j]);
        }

        // Check for velocity discontinuities
        println!("\nVelocity gradient (due/ds):");
        for j in 48..65.min(ue_upper.len()-1) {
            let ds = s_upper[j+1] - s_upper[j];
            if ds > 1e-10 {
                let due_upper = (ue_upper[j+1] - ue_upper[j]) / ds;
                let due_lower = (ue_lower[j+1] - ue_lower[j]) / ds;
                println!("  Station {}: due/ds_upper={:+.4}, due/ds_lower={:+.4}",
                    j, due_upper, due_lower);
            }
        }
    }
}
