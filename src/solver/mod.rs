//! Viscous-inviscid coupling solver
//!
//! This module implements the coupled viscous-inviscid iteration
//! (VISCAL) that iterates between the panel method and boundary
//! layer solver until convergence.
//!
//! # Components
//! - `viscal` - Main coupling iteration loop
//! - `setbl` - SETBL Newton system setup (XFOIL-compatible)
//! - `polar` - Alpha sweep for polar generation

pub mod blstate;
pub mod pointers;
mod polar;
pub mod psilin;
pub mod setbl;
pub mod velocity;
mod viscal;
pub mod xywake;

pub use polar::*;
pub use setbl::*;
pub use viscal::*;
