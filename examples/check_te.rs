use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;

fn main() {
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);

    let n = airfoil.x.len();
    println!("TE Geometry Check:");
    println!("==================");
    println!("Number of nodes: {}", n);
    println!();
    println!("First node (upper TE):  ({:.8}, {:.8})", airfoil.x[0], airfoil.y[0]);
    println!("Last node (lower TE):   ({:.8}, {:.8})", airfoil.x[n-1], airfoil.y[n-1]);
    println!();

    let dxte = airfoil.x[0] - airfoil.x[n-1];
    let dyte = airfoil.y[0] - airfoil.y[n-1];
    let dste = (dxte*dxte + dyte*dyte).sqrt();
    let chord = 1.0;  // NACA airfoil normalized to chord = 1

    println!("TE gap: dx={:.2e}, dy={:.2e}, distance={:.2e}", dxte, dyte, dste);
    println!("Gap/chord ratio: {:.2e}", dste / chord);
    println!("XFOIL SHARP threshold: 0.0001 * chord = {:.2e}", 0.0001 * chord);
    println!("Is SHARP: {}", dste < 0.0001 * chord);
    println!();

    // Check the TE panel itself
    println!("TE Panel (panel N-1, connecting node {} to node 0):", n-1);
    let dx_panel = airfoil.x[0] - airfoil.x[n-1];
    let dy_panel = airfoil.y[0] - airfoil.y[n-1];
    let panel_len = (dx_panel*dx_panel + dy_panel*dy_panel).sqrt();
    println!("  Length: {:.6}", panel_len);

    // Check velocity at TE
    let solution = solve_inviscid(&airfoil);
    let vel = solution.velocity_at_alpha(0.0);

    println!("\nVelocity at nodes near TE:");
    println!("  vel[0] (upper TE)    = {:.6}", vel[0]);
    println!("  vel[1]               = {:.6}", vel[1]);
    println!("  vel[{}]              = {:.6}", n-2, vel[n-2]);
    println!("  vel[{}] (lower TE)   = {:.6}", n-1, vel[n-1]);

    // Cp at TE
    let cp0 = 1.0 - vel[0]*vel[0];
    let cpn = 1.0 - vel[n-1]*vel[n-1];
    println!("\nCp at TE nodes:");
    println!("  Cp[0]   = {:.2}", cp0);
    println!("  Cp[{}] = {:.2}", n-1, cpn);
}
