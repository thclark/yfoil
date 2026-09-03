//! Compare yfoil BL values with XFOIL at each station
//!
//! Outputs the same format as instrumented XFOIL for direct comparison

use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::{solve_viscous, ViscalConfig};

fn main() {
    println!("# IBL IS TURB X S THETA DSTAR HK UE CF RT AMPL");

    let geom = naca_4digit("0012", 160).unwrap(); // XFOIL uses 160 panels by default
    let airfoil = create_paneled_airfoil(&geom);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = ViscalConfig::default();

    let alpha = 0.0_f64.to_radians();
    let result = solve_viscous(&airfoil, alpha, &cond, &config);

    // Output upper surface (side 1 in XFOIL)
    let mut turb = false;
    for (i, st) in result.bl.upper.iter().enumerate() {
        let ibl = i + 2; // XFOIL starts at IBL=2
        let x = result.bl.x_upper.get(i).copied().unwrap_or(0.0);
        let s = result.bl.s_upper.get(i).copied().unwrap_or(0.0);
        let rt = st.ue * st.theta / cond.nu;

        // Check transition
        if !turb && st.n_amp >= cond.ncrit {
            turb = true;
        }

        let turb_flag = if turb { "T" } else { "F" };

        println!(
            "{:4} {:1} {} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E}",
            ibl, 1, turb_flag,
            x, s, st.theta, st.dstar, st.hk, st.ue, st.cf, rt, st.n_amp
        );
    }

    // Output lower surface (side 2 in XFOIL)
    turb = false;
    for (i, st) in result.bl.lower.iter().enumerate() {
        let ibl = i + 2;
        let x = result.bl.x_lower.get(i).copied().unwrap_or(0.0);
        let s = result.bl.s_lower.get(i).copied().unwrap_or(0.0);
        let rt = st.ue * st.theta / cond.nu;

        if !turb && st.n_amp >= cond.ncrit {
            turb = true;
        }

        let turb_flag = if turb { "T" } else { "F" };

        println!(
            "{:4} {:1} {} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E} {:14.6E}",
            ibl, 2, turb_flag,
            x, s, st.theta, st.dstar, st.hk, st.ue, st.cf, rt, st.n_amp
        );
    }

    eprintln!("\nSummary:");
    eprintln!("  Converged: {}", result.converged);
    eprintln!("  CL: {:.6}", result.cl);
    eprintln!("  CD: {:.6}", result.cd);
    eprintln!("  xtr_upper: {:.4}", result.xtr_upper);
    eprintln!("  xtr_lower: {:.4}", result.xtr_lower);
}
