//! Airfoil geometry handling
//!
//! This module provides:
//! - Airfoil coordinate representation and validation
//! - Cubic spline interpolation
//! - Panel distribution and repaneling
//! - NACA airfoil generation
//! - File I/O for JSON and .dat formats

mod airfoil;
mod io;
mod naca;
mod panel;
mod spline;

pub use airfoil::{Geometry, InvalidGeometryError, PaneledAirfoil};
pub use io::{
    read_dat_file, read_geometry_auto, read_geometry_from_file, write_dat_file,
    write_geometry_to_json, GeometryReadError,
};
pub use naca::{naca_4digit, naca_5digit, NacaError};
pub use panel::{create_paneled_airfoil, repanel};
pub use spline::{d2val, deval, seval, spline};
