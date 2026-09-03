//! Airfoil geometry handling
//!
//! This module provides functions and data structures for working with 2D airfoil
//! geometry, including coordinate representation, spline interpolation, panel
//! distribution, NACA generation, and file I/O.
//!
//! # CLI Usage
//!
//! The geometry module is exposed through the `yfoil geom` subcommand, which provides
//! five operations:
//!
//! ## `yfoil geom convert` - Format Conversion
//!
//! Convert between geometry file formats (JSON and .dat):
//!
//! ```text
//! yfoil geom convert <input> [OPTIONS]
//!
//! Arguments:
//!   <input>              Input file path (.json or .dat)
//!
//! Options:
//!   --to <FORMAT>        Output format: json, dat [default: json]
//!   -o, --output <PATH>  Output file path [default: input with new extension]
//!   --name <NAME>        Airfoil name for .dat output [default: "Airfoil"]
//! ```
//!
//! Examples:
//! ```text
//! # Convert .dat to JSON
//! yfoil geom convert airfoil.dat --to json
//!
//! # Convert JSON to .dat with custom name
//! yfoil geom convert airfoil.json --to dat --name "My Airfoil"
//! ```
//!
//! ## `yfoil geom naca` - NACA Airfoil Generation
//!
//! Generate standard NACA 4-digit or 5-digit airfoil profiles:
//!
//! ```text
//! yfoil geom naca <spec> [OPTIONS]
//!
//! Arguments:
//!   <spec>               NACA designation (e.g., "0012", "4412", "23015")
//!
//! Options:
//!   -n, --panels <N>     Number of panels [default: 160]
//!   --to <FORMAT>        Output format: json, dat [default: json]
//!   -o, --output <PATH>  Output file path [default: naca<spec>.<format>]
//! ```
//!
//! ### NACA 4-Digit Format (e.g., "4412")
//! - 1st digit: Maximum camber as percentage of chord (0-9)
//! - 2nd digit: Position of maximum camber in tenths of chord (0-9)
//! - 3rd-4th digits: Maximum thickness as percentage of chord (00-99)
//!
//! ### NACA 5-Digit Format (e.g., "23012")
//! - 1st digit: Design lift coefficient × (20/3)
//! - 2nd digit: Position of maximum camber × 20
//! - 3rd digit: 0 = standard camber, 1 = reflex camber
//! - 4th-5th digits: Maximum thickness as percentage of chord
//!
//! Examples:
//! ```text
//! # Generate NACA 0012 with 160 panels
//! yfoil geom naca 0012
//!
//! # Generate NACA 4412 with 200 panels as .dat
//! yfoil geom naca 4412 -n 200 --to dat
//!
//! # Generate NACA 23015 (5-digit)
//! yfoil geom naca 23015 -o my_airfoil.json
//! ```
//!
//! ## `yfoil geom repanel` - Redistribute Panel Points
//!
//! Redistribute panel points on an existing geometry:
//!
//! ```text
//! yfoil geom repanel <input> [OPTIONS]
//!
//! Arguments:
//!   <input>              Input file path (.json or .dat)
//!
//! Options:
//!   -n, --panels <N>     Target number of panels [default: 160]
//!   --method <METHOD>    Repaneling method [default: xfoil]
//!   --le-ratio <RATIO>   LE/TE panel density ratio (cosine method only) [default: 0.15]
//!   -o, --output <PATH>  Output file path [default: <input>_repaneled.json]
//! ```
//!
//! ### Repaneling Methods
//!
//! **`xfoil` (default)**: Uses XFOIL's curvature-based PANE algorithm (PANGEN subroutine).
//! Distributes panels based on local surface curvature, placing more panels in
//! high-curvature regions (leading edge) and fewer in low-curvature regions (mid-chord).
//! This produces panel distributions that match XFOIL exactly.
//!
//! **`cosine`**: Uses modified cosine spacing with a configurable LE/TE density ratio.
//! The `--le-ratio` parameter controls panel clustering:
//! - Values < 1.0: Finer panels at LE, coarser at TE (recommended)
//! - Value = 1.0: Symmetric cosine spacing
//! - Values > 1.0: Finer panels at TE, coarser at LE
//!
//! Examples:
//! ```text
//! # Repanel using XFOIL's PANE algorithm (default)
//! yfoil geom repanel airfoil.dat -n 180
//!
//! # Repanel using modified cosine spacing
//! yfoil geom repanel airfoil.json --method cosine --le-ratio 0.2
//!
//! # Explicit XFOIL method
//! yfoil geom repanel airfoil.dat --method xfoil -n 200 -o repaneled.json
//! ```
//!
//! ## `yfoil geom info` - Display Geometry Information
//!
//! Display statistics about an airfoil geometry:
//!
//! ```text
//! yfoil geom info <input> [OPTIONS]
//!
//! Arguments:
//!   <input>              Input file path (.json or .dat)
//!
//! Options:
//!   -o, --output <PATH>  Output JSON file path (if not specified, prints summary to stdout)
//! ```
//!
//! ### Summary Output (stdout)
//!
//! When no `--output` is specified, prints a human-readable summary:
//! - Number of points and chord length
//! - X/Y coordinate ranges
//! - Maximum thickness
//! - Trailing edge gap and whether TE is sharp
//! - Reference point for moments
//! - Leading edge index and arc length
//! - Total arc length around the airfoil
//! - Maximum curvature
//! - First/last point coordinates (trailing edge)
//!
//! ### JSON Output (file)
//!
//! When `--output` is specified, writes a comprehensive JSON file containing:
//! - All summary statistics
//! - Full distributions: coordinates (x, y), arc length (s), curvature,
//!   panel angles, and normal vectors (nx, ny) at each node
//!
//! Examples:
//! ```text
//! # Print summary to stdout
//! yfoil geom info naca0012.json
//!
//! # Write full info to JSON file
//! yfoil geom info naca0012.json -o naca0012_info.json
//! ```
//!
//! ## `yfoil geom plot` - Visualize Geometry (Feature-gated)
//!
//! Generate an interactive HTML plot of the airfoil shape:
//!
//! ```text
//! yfoil geom plot <input> [OPTIONS]
//!
//! Arguments:
//!   <input>              Input file path (.json or .dat)
//!
//! Options:
//!   -o, --output <PATH>  Output HTML file [default: geometry.html]
//! ```
//!
//! Note: Requires the `plotting` feature flag:
//! ```text
//! cargo build --features plotting
//! ```
//!
//! # File Formats
//!
//! ## JSON Format
//!
//! The native YFoil format stores geometry as a JSON object:
//!
//! ```json
//! {
//!   "reference": [0.25, 0.0],
//!   "x_c": [1.0, 0.8, 0.5, 0.2, 0.0, 0.2, 0.5, 0.8, 1.0],
//!   "y_c": [0.0, 0.02, 0.04, 0.03, 0.0, -0.03, -0.04, -0.02, 0.0]
//! }
//! ```
//!
//! ## Selig .dat Format
//!
//! Standard XFOIL/Selig format with airfoil name on first line:
//!
//! ```text
//! NACA 0012
//!   1.000000   0.001260
//!   0.950000   0.008000
//!   ...
//! ```
//!
//! Both Selig (wrap-around) and Lednicer (separate upper/lower surfaces) formats
//! are supported for reading.
//!
//! # Coordinate Conventions
//!
//! - **Normalization**: All coordinates are normalized by chord (x/c, y/c)
//! - **Point ordering**: TE → upper surface → LE → lower surface → TE
//! - **Reference point**: Default [0.25, 0.0] for quarter-chord moment calculations
//! - **Trailing edge**: Located at x ≈ 1.0, validated to be at both first and last points
//!
//! # Data Structures
//!
//! - [`Geometry`]: Raw airfoil coordinates from input files
//! - [`PaneledAirfoil`]: Fully processed geometry ready for aerodynamic analysis,
//!   including arc length parameterization, spline derivatives, normal vectors,
//!   and panel angles
//!
//! # Key Functions
//!
//! - [`create_paneled_airfoil`]: Convert raw geometry to analysis-ready form
//! - [`naca_4digit`], [`naca_5digit`]: Generate NACA airfoil profiles
//! - [`repanel_xfoil`]: XFOIL's curvature-based PANE algorithm (default repaneling method)
//! - [`repanel_cosine`]: Modified cosine spacing repaneling (alternative method)
//! - [`spline`], [`seval`], [`deval`], [`d2val`]: Cubic spline interpolation

mod airfoil;
mod io;
mod naca;
mod panel;
mod spline;

pub use airfoil::{Geometry, InvalidGeometryError, PaneledAirfoil};
pub use io::{
    read_dat_file, read_geometry_auto, read_geometry_from_file, write_dat_file, write_geometry_to_json,
    GeometryReadError,
};
pub use naca::{naca_4digit, naca_5digit, NacaError};
pub use panel::{create_paneled_airfoil, repanel_cosine, repanel_xfoil, PaneConfig};
pub use spline::{d2val, deval, seval, spline};
