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
//! - `wake` - Wake boundary layer model

mod closure;
mod march;
mod newton;
mod state;
mod wake;

pub use closure::*;
pub use march::*;
pub use newton::*;
pub use state::*;
pub use wake::*;
