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
pub mod ggcalc;
pub mod ludcmp;
pub mod pointers;
mod polar;
pub mod psilin;
pub mod qdcalc;
pub mod setbl;
pub mod setbl_legacy;
pub mod update;
pub mod velocity;
mod viscal;
pub mod xywake;

pub use polar::*;
pub use setbl_legacy::*;
pub use viscal::*;
