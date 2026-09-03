//! Boundary layer solver
//!
//! This module implements the integral boundary layer equations
//! with eN transition prediction.
//!
//! # Components
//! - `state` - BL variable definitions
//! - `closure` - Closure relations (Cf, H, Hs, etc.)
//! - `transition` - eN method for laminar-turbulent transition

//! - `newton` - Newton iteration system
//! - `system` - BL Newton system data structures (XFOIL-compatible)
//! - `blsolv` - BLSOLV block solver for coupled Newton system
//! - `wake` - Wake boundary layer model

pub mod blsolv;
pub mod blsys;
mod closure;
pub mod gauss;
pub mod march_legacy;
pub mod mrchdu;
pub mod mrchue;
mod newton;
mod state;
pub mod system;
mod wake;

pub use blsolv::*;
pub use closure::*;
pub use newton::*;
pub use state::*;
pub use wake::*;
