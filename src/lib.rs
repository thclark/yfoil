//! YFoil - Rust-based aerofoil analysis tool
//!
//! YFoil reproduces the core analysis functionality of XFOIL in Rust.
//! It provides viscous/inviscid analysis of 2D airfoil sections using:
//! - Panel method for inviscid flow
//! - Integral boundary layer solver with eN transition
//! - Coupled viscous-inviscid iteration
//!
//! # Example
//!
//! ```ignore
//! use yfoil::geometry::read_geometry_from_file;
//!
//! let geometry = read_geometry_from_file("airfoil.json")?;
//! ```

pub mod error;
pub mod geometry;

// Placeholder modules (to be implemented)
pub mod bl;
pub mod forces;
pub mod output;
pub mod panel;
pub mod solver;

pub use error::{Result, YfoilError};
