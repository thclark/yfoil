//! Output and visualization
//!
//! This module handles result serialization and optional plotting.
//!
//! # Components
//! - `results` - Result structs (JSON serializable)
//! - `foil` - Geometry, wake and per-station BL output of an operating point
//! - `svg` / `foil_plot` - Precision-SVG canvas and the foil plot (feature-gated)
//! - `plot` - Optional plotting with plotters library (feature-gated)

mod foil;
mod results;

#[cfg(feature = "plotting")]
mod foil_plot;
#[cfg(feature = "plotting")]
mod plot;
#[cfg(feature = "plotting")]
mod svg;

pub use foil::*;
pub use results::*;

#[cfg(feature = "plotting")]
pub use foil_plot::*;
#[cfg(feature = "plotting")]
pub use plot::*;
