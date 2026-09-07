//! Boundary layer solver
//!
//! This module implements the integral boundary layer equations
//! with eN transition prediction.
//!
//! # Components
//! - `closure` - Closure relations (Cf, H, Hs, etc.)
//! - `params` - BLPAR constants, flow regime, global BL parameters (/V_VAR/)
//! - `station` - station state COM1/COM2 with BLPRV/BLKIN/BLVAR, BLMID, DSLIM
//! - `transition` - TRCHEK2, DAMPL/DAMPL2, AXSET
//! - `difference` - BLDIF/TRDIF local Newton system (/V_SYS/)
//! - `system` - re-exports of the four above
//! - `blsolv` - BLSOLV block solver for coupled Newton system

pub mod blsolv;
pub mod blsys;
mod closure;
pub mod difference;
pub mod gauss;
pub mod mrchdu;
pub mod mrchue;
pub mod params;
pub mod station;
pub mod system;
pub mod transition;

pub use blsolv::{solve_newton_system, solve_newton_system_traced, BlsolvTrace, NewtonDeltas, NewtonSystem};
pub use closure::{
    cdiss_laminar, cdiss_wake, cf_laminar, cf_turbulent, hk_from_h, hstar_laminar, hstar_turbulent, hstarstar, Closure,
};
