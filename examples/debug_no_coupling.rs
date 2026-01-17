//! Run BL with inviscid velocities only (no coupling)

use yfoil::bl::{FlowConditions, march_newton, NewtonConfig};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    let inviscid = solve_inviscid(&airfoil);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();

    // Get node-based inviscid velocities
    let qinv = inviscid.velocity_at_nodes(0.0);

    // Extract upper surface from stagnation to TE
    let le_idx = airfoil.le_index;
    let mut x = Vec::new();
    let mut s = Vec::new();
    let mut ue = Vec::new();
    let mut arc_len = 0.0;

    for j in 0..=le_idx {
        let i = le_idx - j;
        x.push(airfoil.x[i]);
        s.push(arc_len);
        ue.push(qinv[i].abs());

        if i > 0 {
            let dx = airfoil.x[i-1] - airfoil.x[i];
            let dy = airfoil.y[i-1] - airfoil.y[i];
            arc_len += (dx*dx + dy*dy).sqrt();
        }
    }

    // March BL with inviscid velocities
    let bl = march_newton(&ue, &s, &cond, &config);

    println!("=== BL with inviscid Ue (no coupling) ===");
    println!("Station  X          Ue(inv)    theta      dstar      H");
    for (j, st) in bl.iter().take(15).enumerate() {
        println!("{:4}   {:10.6}  {:10.6}  {:10.3e}  {:10.3e}  {:6.3}",
                 j+2, x[j], ue[j], st.theta, st.dstar, st.hk);
    }

    // TE values
    if let Some(te) = bl.last() {
        println!("\nTE: theta={:.3e}, dstar={:.3e}, Ue={:.3}",
                 te.theta, te.dstar, ue.last().unwrap());
    }
}
