//! (c) yFoil functionality.
//!
//! yFoil-only functionality with no XFOIL equivalent: CLI, I/O, geometry generators, output
//! records, ids — against its specification, external references or regression snapshots.
//!
//! Categories and rules: `docs/conventions/testing.md`.

#[path = "../common/fixtures/mod.rs"]
mod fixtures;
#[path = "../common/utilities/mod.rs"]
mod utilities;

mod cli_geometry;
mod cli_plot;
mod fixed_cl;
mod geometry;
mod geometry_errors;
mod ids;
mod naca456;
mod output;
mod panelling;
mod polar_driver;
mod repanel_cosine;
mod validity_record;
