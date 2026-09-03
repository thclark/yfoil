//! Inviscid panel method
//!
//! This module implements the vortex panel method for computing
//! the inviscid flow around an airfoil.
//!
//! # Components
//! - `influence` - Influence coefficient calculations
//! - `solver` - Linear system assembly and solution
//! - `wake` - Wake panel generation and influence

pub mod influence;
pub mod solver;
pub mod wake;

pub use influence::*;
pub use solver::*;
pub use wake::*;
