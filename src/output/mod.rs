//! Output and visualization
//!
//! This module handles result serialization and optional plotting.
//!
//! # Components
//! - `results` - Result structs (JSON serializable)
//! - `plot` - Optional plotting with plotters library (feature-gated)

mod results;

#[cfg(feature = "plotting")]
mod plot;

pub use results::*;

#[cfg(feature = "plotting")]
pub use plot::*;
