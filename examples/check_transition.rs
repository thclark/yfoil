use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::bl::{FlowConditions, NewtonConfig};
use yfoil::solver::{find_stagnation_point, extract_upper_surface};
use yfoil::bl::march_newton;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);
    let stag = find_stagnation_point(&airfoil, &vel);

    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &vel, stag);
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl_upper = march_newton(&ue_upper, &s_upper, &cond, &config);

    println!("=== Transition Location Analysis ===\n");
    println!("Flow conditions: Re = 1e6, Ncrit = 9\n");

    // Find transition
    let mut trans_idx = None;
    for (i, r) in bl_upper.iter().enumerate() {
        if r.n_amp >= 9.0 && trans_idx.is_none() {
            trans_idx = Some(i);
        }
    }

    if let Some(idx) = trans_idx {
        let x_trans = x_upper[idx];
        let s_trans = s_upper[idx];
        let s_max = s_upper.last().unwrap();
        println!("Transition detected:");
        println!("  Station: {}", idx);
        println!("  x/c = {:.4}", x_trans);
        println!("  s/s_max = {:.4}", s_trans / s_max);
        println!("  n_amp = {:.4}", bl_upper[idx].n_amp);
        println!("  Hk at transition = {:.4}", bl_upper[idx].hk);
    }

    // Show Hk evolution approaching transition
    println!("\n=== Hk Evolution (stations 40-65) ===");
    println!("{:>5} {:>10} {:>10} {:>10} {:>10}", "Stn", "x/c", "Hk", "n_amp", "Rtheta");
    for i in 40..65.min(bl_upper.len()) {
        let rt = ue_upper[i] * bl_upper[i].theta / cond.nu;
        println!("{:>5} {:>10.4} {:>10.4} {:>10.4} {:>10.1}",
            i, x_upper[i], bl_upper[i].hk, bl_upper[i].n_amp, rt);
    }

    // Compare expected XFOIL transition
    println!("\n=== Expected XFOIL Results ===");
    println!("XFOIL NACA 0012, Re=1e6, α=0°:");
    println!("  Expected transition: x/c ≈ 0.45-0.50");
    println!("  Expected CD ≈ 0.0065-0.0070");
    println!();
    println!("Our results:");
    println!("  Transition x/c = {:.4}", x_upper[trans_idx.unwrap_or(0)]);

    // Check when Hk first exceeds 2.5 (where amplification becomes significant)
    println!("\n=== Hk > 2.5 Analysis ===");
    for (i, r) in bl_upper.iter().enumerate() {
        if r.hk > 2.5 {
            println!("  Hk first exceeds 2.5 at station {}, x/c = {:.4}", i, x_upper[i]);
            break;
        }
    }

    // Check when amplification rate becomes significant
    println!("\n=== Amplification Rate Analysis ===");
    use yfoil::bl::{hkin, dampl};
    println!("{:>5} {:>10} {:>10} {:>10} {:>10}", "Stn", "x/c", "Hk", "ax", "dn/ds");
    for i in 30..60.min(bl_upper.len()) {
        let rt = ue_upper[i] * bl_upper[i].theta / cond.nu;
        let (hk, _, _) = hkin(bl_upper[i].h, cond.msq);
        let (ax, _, _, _) = dampl(hk, bl_upper[i].theta, rt);
        if ax > 0.0 {
            let ds = if i > 0 { s_upper[i] - s_upper[i-1] } else { 0.0 };
            let dn_ds = ax;
            println!("{:>5} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
                i, x_upper[i], hk, ax, dn_ds);
        }
    }
}
