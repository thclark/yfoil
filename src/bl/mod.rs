//! Boundary layer solver
//!
//! This module implements the integral boundary layer equations
//! with eN transition prediction.
//!
//! # Components
//! - `state` - BL variable definitions
//! - `closure` - Closure relations (Cf, H, Hs, etc.)
//! - `transition` - eN method for laminar-turbulent transition

//! - `system` - BL Newton system data structures (XFOIL-compatible)
//! - `blsolv` - BLSOLV block solver for coupled Newton system

pub mod blsolv;
pub mod blsys;
mod closure;
pub mod gauss;
pub mod mrchdu;
pub mod mrchue;
mod state;
pub mod system;

pub use blsolv::*;
pub use closure::*;
pub use state::*;
