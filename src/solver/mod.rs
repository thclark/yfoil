//! Viscous-inviscid coupling solver
//!
//! This module implements the coupled viscous-inviscid iteration
//! (VISCAL) that iterates between the panel method and boundary
//! layer solver until convergence.
//!
//! # Components
//! - `viscal` - Main coupling iteration loop
//! - `polar` - Alpha sweep for polar generation

mod polar;
mod viscal;

pub use polar::*;
pub use viscal::*;
