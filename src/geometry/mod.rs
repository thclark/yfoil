//! Aerofoil geometry handling
//!
//! This module provides functions and data structures for working with 2D aerofoil
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
//!   --name <NAME>        Aerofoil name for .dat output [default: "Aerofoil"]
//! ```
//!
//! Examples:
//! ```text
//! # Convert .dat to JSON
//! yfoil geom convert aerofoil.dat --to json
//!
//! # Convert JSON to .dat with custom name
//! yfoil geom convert aerofoil.json --to dat --name "My Aerofoil"
//! ```
//!
//! ## `yfoil geom naca` - NACA Aerofoil Generation
//!
//! Generate standard NACA 4-digit or 5-digit aerofoil profiles:
//!
//! ```text
//! yfoil geom naca <spec> [OPTIONS]
//!
//! Arguments:
//!   <spec>               NACA designation (e.g., "0012", "4412", "23015")
//!
//! Options:
//!   -n, --panels <N>     Number of panels [default: 160]
//!   --thickness <T>      perpendicular (the NACA definition) | vertical (XFOIL's NACA4/NACA5) [default: perpendicular]
//!   --sharp              Close the trailing edge
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
//! ## `yfoil geom repanel` - Redistribute panel points
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
//!   --method <METHOD>    Repanelling method: curvature | cosine [default: curvature]
//!   --te-le-ratio <R>    TE/LE panel density ratio, XFOIL's CTERAT (cosine method only) [default: 0.15]
//!   -o, --output <PATH>  Output file path [default: <input>_repanelled.json]
//! ```
//!
//! ### Repanelling Methods
//!
//! **`curvature` (default)**: Uses XFOIL's curvature-based PANE algorithm (PANGEN subroutine).
//! Distributes panels based on local surface curvature, placing more panels in
//! high-curvature regions (leading edge) and fewer in low-curvature regions (mid-chord).
//! This produces panel distributions that match XFOIL exactly.
//!
//! **`cosine`**: Uses modified cosine spacing with a configurable TE/LE density ratio.
//! The `--te-le-ratio` parameter controls panel clustering:
//! - Values < 1.0: Finer panels at TE, coarser at LE
//! - Value = 1.0: Symmetric cosine spacing
//! - Values > 1.0: Finer panels at LE, coarser at TE
//!
//! Examples:
//! ```text
//! # Repanel using XFOIL's PANE algorithm (default)
//! yfoil geom repanel aerofoil.dat -n 180
//!
//! # Repanel using modified cosine spacing
//! yfoil geom repanel aerofoil.json --method cosine --te-le-ratio 0.2
//!
//! # Explicit XFOIL method
//! yfoil geom repanel aerofoil.dat --method curvature -n 200 -o repanelled.json
//! ```
//!
//! ## `yfoil geom info` - Display Geometry Information
//!
//! Display statistics about an aerofoil geometry:
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
//! - Total arc length around the aerofoil
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
//! # File Formats
//!
//! ## JSON Format
//!
//! The native yFoil format stores geometry as a JSON object:
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
//! Standard XFOIL/Selig format with aerofoil name on first line:
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
//! - [`Geometry`]: Raw aerofoil coordinates from input files
//! - [`PanelledFoil`]: Fully processed geometry ready for aerodynamic analysis,
//!   including arc length parameterization, spline derivatives, normal vectors,
//!   and panel angles
//!
//! # Key Functions
//!
//! - [`panel_foil`]: Convert raw geometry to analysis-ready form
//! - [`naca_4digit`], [`naca_5digit`]: generate NACA aerofoil profiles
//! - [`repanel_by_curvature`]: XFOIL's curvature-based PANE algorithm (default repanelling method)
//! - [`repanel_cosine`]: modified cosine spacing (alternative method)
//! - [`spline_derivatives`], [`spline_value`], [`spline_slope`], [`spline_second_derivative`]: cubic spline interpolation

mod airfoil;
mod io;
mod naca;
mod panel;
mod spline;

pub use airfoil::{Geometry, InvalidGeometryError, PanelledFoil};
pub use io::{
    read_dat_file, read_geometry_auto, read_geometry_from_file, write_dat_file, write_geometry_to_json,
    GeometryReadError,
};
pub use naca::{
    naca_4digit, naca_4digit_vertical, naca_5digit, naca_5digit_vertical, NacaError, Thickness, XFOIL_NACA_NSIDE,
};
pub use panel::{
    arc_coordinate, curvature, find_le, panel_foil, repanel_by_curvature, repanel_cosine, solve_tridiagonal,
    spline_segmented, PaneConfig,
};
pub use spline::{spline_derivatives, spline_second_derivative, spline_slope, spline_value};
