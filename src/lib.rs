//! yFoil - Rust-based aerofoil analysis tool
//!
//! yFoil reproduces the core analysis functionality of XFOIL in Rust.
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

#![deny(rustdoc::broken_intra_doc_links)]
#![forbid(unsafe_code)]

pub mod error;
pub mod geometry;

pub mod bl;
pub mod output;
pub mod solver;

pub use error::{Result, YfoilError};
