//! Boundary layer solver
//!
//! This module implements the integral boundary layer equations
//! with eN transition prediction.
//!
//! # Components
//! - `state` - BL variable definitions
//! - `closure` - Closure relations (Cf, H, Hs, etc.)
//! - `transition` - eN method for laminar-turbulent transition
//! - `march` - BL marching algorithms (direct and mixed modes)
//! - `newton` - Newton iteration system
//! - `system` - BL Newton system data structures (XFOIL-compatible)
//! - `blsolv` - BLSOLV block solver for coupled Newton system
//! - `wake` - Wake boundary layer model

pub mod blsolv;
mod closure;
pub mod gauss;
mod march;
pub mod mrchdu;
mod newton;
mod state;
pub mod system;
mod wake;

pub use blsolv::*;
pub use closure::*;
pub use march::*;
pub use newton::*;
pub use state::*;
pub use wake::*;
