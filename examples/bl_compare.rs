//! Compare YFoil BL march with XFOIL debug output

use yfoil::bl::{FlowConditions, NewtonConfig, march_newton, hkin};
use yfoil::geometry::{create_paneled_airfoil, naca_4digit};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{find_stagnation_point, extract_upper_surface};

fn main() {
    // Generate same geometry as XFOIL
    let geom = naca_4digit("0012", 160).unwrap();
    let airfoil = create_paneled_airfoil(&geom);
    
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);
    
    let (x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);
    
    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);
    let config = NewtonConfig::default();
    let bl = march_newton(&ue_upper, &s_upper, &cond, &config);
    
    // Output in XFOIL-like format
    println!("# IBL IS TURB X S THETA DSTAR HK UE CF RT AMPL");
    
    for (i, r) in bl.iter().enumerate() {
        if i >= x_upper.len() { break; }
        
        let turb = if r.n_amp >= cond.ncrit { 'T' } else { 'F' };
        let rt = r.ue * r.theta / cond.nu;
        
        println!("{:4} 1 {}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}  {:12.6E}",
                 i+2, turb,
                 x_upper[i],
                 s_upper[i],
                 r.theta,
                 r.dstar,
                 r.hk,
                 r.ue,
                 r.cf,
                 rt,
                 r.n_amp);
    }
}
