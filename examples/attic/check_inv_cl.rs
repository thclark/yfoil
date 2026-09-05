use yfoil::bl::FlowConditions;
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::solver::solve_inviscid_only;

fn main() {
    let geom = naca_4digit("0012", 160).expect("Failed");
    let airfoil = create_paneled_airfoil(&geom);

    let alpha_deg = 2.0f64;
    let cond = FlowConditions::new(1e6, 0.0, 9.0, 1.0);

    let result = solve_inviscid_only(&airfoil, alpha_deg.to_radians(), cond.mach);

    println!("Inviscid-only at alpha={}deg:", alpha_deg);
    println!("CL = {:.6}", result.cl);
    println!("CM = {:.6}", result.cm);

    // XFOIL inviscid at alpha=2 is approximately CL=0.229
    println!("\nXFOIL inviscid CL ≈ 0.229");
    println!("Difference: {:.1}%", 100.0 * (result.cl - 0.229) / 0.229);
}
