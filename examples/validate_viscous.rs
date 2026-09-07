//! Smoke check of the viscous analysis against a few hand-transcribed XFOIL screen values
//! (4–5 significant figures). The gates are the instrumented-dump tests under `tests/`; this
//! only shows that the CLI-level path runs and lands in the right place.

use yfoil::geometry::{naca_4digit, panel_foil};
use yfoil::solver::analysis::{analyse, FlowConditions};

fn main() {
    let geometry = naca_4digit("0012", 160).expect("NACA 0012");
    let airfoil = panel_foil(&geometry);
    let spec = FlowConditions {
        re: Some(1.0e6),
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    };

    println!("NACA 0012, Re = 1e6, M = 0, Ncrit = 9 (N = 160)");
    println!(
        "{:>6} {:>10} {:>10} {:>10} {:>10} {:>8} {:>8} {:>5} {:>4}",
        "alpha", "CL", "CD", "CDf", "CM", "xtr_u", "xtr_l", "iter", "conv"
    );
    for alpha_deg in [0.0_f64, 2.0, 4.0, 6.0] {
        let p = analyse(&airfoil, alpha_deg.to_radians(), &spec);
        println!(
            "{:>6.2} {:>10.6} {:>10.6} {:>10.6} {:>10.6} {:>8.4} {:>8.4} {:>5} {:>4}",
            alpha_deg,
            p.cl,
            p.cd,
            p.cd_friction,
            p.cm,
            p.transition_upper[0],
            p.transition_lower[0],
            p.iterations,
            if p.converged { "yes" } else { "NO" }
        );
    }
    println!("XFOIL (screen, 4–5 s.f.) at alpha = 0: CL 0.0000  CD 0.00541 (smoke reference only)");
}
