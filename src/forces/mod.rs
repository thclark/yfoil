//! Aerodynamic force calculation
//!
//! This module computes pressure coefficients and integrated
//! aerodynamic forces (CL, CD, CM).
//!
//! # Components
//! - `pressure` - Cp calculation with compressibility corrections
//! - `integrate` - Force and moment integration

pub mod integrate;
pub mod pressure;

pub use integrate::*;
pub use pressure::*;
