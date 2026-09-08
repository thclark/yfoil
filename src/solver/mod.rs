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

pub mod analysis;
pub mod blstate;
pub mod clcalc;
pub mod ggcalc;
pub mod ludcmp;
pub mod pointers;
pub mod psilin;
pub mod qdcalc;
pub mod setbl;
pub mod specal;
pub mod update;
pub mod velocity;
pub mod viscal;
pub mod xywake;

pub use analysis::{
    analyse, compute_polar, compute_polar_with, FlowConditions, PointResult, PolarConfig, PolarResult, Session,
};
