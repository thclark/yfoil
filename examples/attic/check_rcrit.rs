use yfoil::bl::{amplification_rate, rtheta_crit};

fn main() {
    println!("=== Critical Rtheta vs Hk ===\n");
    println!("At low Hk, Rcrit should be low to allow amplification to start early.\n");
    println!("{:>6} {:>12} {:>12}", "Hk", "log10(Rcrit)", "Rcrit");
    println!("{:-<34}", "");

    for hk in [2.2_f64, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.5, 4.0] {
        let grcrit = rtheta_crit(hk);
        let rcrit = 10.0_f64.powf(grcrit);
        println!("{:>6.2} {:>12.4} {:>12.1}", hk, grcrit, rcrit);
    }

    // Check actual XFOIL values from xblsys.f DAMPL subroutine
    // The critical Rtheta formula in XFOIL is:
    // GRCRIT = 2.492*(HMI)**0.43 + 0.7*(TANH(14*HMI - 9.24) + 1)
    // where HMI = 1/(HK-1)
    println!("\n=== XFOIL DAMPL Critical Values ===");
    println!("From xblsys.f, same formula. Let's verify:");
    for hk in [2.3_f64, 2.5, 2.7, 3.0] {
        let hmi: f64 = 1.0 / (hk - 1.0);
        let aa = 2.492 * hmi.powf(0.43);
        let bb = (14.0 * hmi - 9.24).tanh();
        let grcrit = aa + 0.7 * (bb + 1.0);
        let rcrit = 10.0_f64.powf(grcrit);
        println!("Hk={:.1}: hmi={:.4}, aa={:.4}, bb={:.4}, log10(Rcrit)={:.4}, Rcrit={:.1}",
            hk, hmi, aa, bb, grcrit, rcrit);
    }

    // Now let's check the amplification rate at some actual BL conditions
    println!("\n=== Amplification Rate at Realistic Conditions ===");
    println!("{:>6} {:>8} {:>12} {:>12}", "Hk", "Rθ", "ax (ours)", "ax*ds@ds=0.02");

    let theta = 0.0005; // Typical momentum thickness
    for (hk, rt) in [(2.4_f64, 300.0), (2.5, 400.0), (2.6, 500.0), (2.7, 600.0), (2.8, 700.0), (3.0, 800.0)] {
        let (ax, _, _, _) = amplification_rate(hk, theta, rt);
        let dn = ax * 0.02; // typical ds
        println!("{:>6.1} {:>8.0} {:>12.4} {:>12.4}", hk, rt, ax, dn);
    }

    // What we need for transition at x/c=0.5:
    // If we have ~20 stations from where amplification starts (say x/c=0.1) to x/c=0.5
    // And ds ≈ 0.02 per station
    // Total ds ≈ 0.4
    // To get n_amp = 9, we need average ax ≈ 9/0.4 = 22.5
    println!("\n=== Required Average Amplification Rate ===");
    println!("For transition at x/c ≈ 0.5:");
    println!("  Assume ~20 stations from amplification onset to transition");
    println!("  Average ds ≈ 0.02");
    println!("  Total distance ≈ 0.4 (arc length)");
    println!("  Required total n_amp = 9");
    println!("  Required average ax = 9/0.4 = 22.5");
    println!();
    println!("Our ax values at stations 30-40 range from 1.1 to 9.5");
    println!("This is much lower than the required ~22.5!");
}
