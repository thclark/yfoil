//! Inviscid panel method
//!
//! This module implements the vortex panel method for computing
//! the inviscid flow around an airfoil.
//!
//! # Components
//! - `influence` - Influence coefficient calculations
//! - `solver` - Linear system assembly and solution

pub mod influence;
pub mod solver;

pub use influence::*;
pub use solver::*;
