//! Aerofoil geometry handling
//!
//! This module provides functions and data structures for working with 2D aerofoil
//! geometry, including coordinate representation, spline interpolation, panel
//! distribution, NACA generation, and file I/O.
//!
//! # CLI Usage
//!
//! The geometry module is exposed through the `yfoil geometry` subcommand (alias `geom`),
//! which provides five operations, each printing its full help when run without arguments:
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
//! Generate a NACA section from its designation — 4-digit (2412), 4-digit modified (0012-34),
//! 5-digit (23012, 23112 reflex), 16-series (16-212), 6-series (63-415) and 6A-series (64A010),
//! see [`series`] — panelled with XFOIL's PANGEN by default:
//!
//! ```text
//! yfoil geom naca <spec> [OPTIONS]
//!
//! Options:
//!   --to <FORMAT>               Output format: json, dat [default: json]
//!   --thickness <T>             perpendicular (the NACA definition, default) | vertical (XFOIL's NACA4/NACA5
//!                               model, always PANGEN-panelled; 4- and 5-digit only)
//!   --a <A>                     Extent of uniform loading of a 6-/16-series mean line, 0..1 [default: 1.0]
//!   -n, --method, PPAR flags, --sharp, --te-gap, --te-blend, --panelling   as `repanel` below
//!   -o, --output <PATH>         Output file path [default: naca<spec>.<format>]
//! ```
//!
//! ## `yfoil geom karman-trefftz` - Kármán–Trefftz Section
//!
//! The conformal map of a circle through ζ = 1 (`--x-centre`, `--y-centre`, `--te-angle`), an
//! analytic section with an exact potential-flow solution; the same panelling options as
//! `repanel`, `pangen` by default.
//!
//! ## `yfoil geom repanel` - Redistribute panel nodes
//!
//! XFOIL's `PANE`/`PPAR` on a loaded geometry (`.json` or `.dat`); the output is JSON only,
//! `<input stem>_repanelled.json` beside the input by default, with the panelling recorded under
//! `generator.panelling`.
//!
//! ```text
//! yfoil geom repanel <input> [OPTIONS]
//!
//! Panelling (both methods):
//!   -n, --panels <N>              Number of panel nodes, XFOIL's NPAN [default: 160]
//!   --method <pangen|cosine>      pangen: XFOIL's PANGEN (default); cosine: yFoil's own, no XFOIL equivalent
//!   --panelling <FILE>            The whole panelling from a JSON file in the record's shape (exclusive with the flags)
//! PANGEN parameters (--method pangen only; an error with cosine):
//!   --curvature-bunching <P>      CVPAR, PPAR menu P [default 1.0]
//!   --te-curvature-ratio <T>      CTERAT, PPAR menu T: fictitious TE curvature / LE curvature [default 0.15]
//!   --refined-curvature-ratio <R> CTRRAT, PPAR menu R [default 0.2]
//!   --refine-upper <X1,X2>        XSREF1, XSREF2, PPAR menu XT (off by default)
//!   --refine-lower <X1,X2>        XPREF1, XPREF2, PPAR menu XB (off by default)
//! Cosine parameters (--method cosine only; an error with pangen):
//!   --cosine-te-bias <B>          1 plain cosine; < 1 coarser at the TE, finer at the LE; > 1 finer at the TE [default 0.15]
//! Trailing edge (both methods, applied after panelling):
//!   --sharp                       Close the trailing edge (exclusive with --te-gap)
//!   --te-gap <GAP> [--te-blend F] XFOIL's TGAP on the panelled nodes
//! ```
//!
//! ### Methods
//!
//! **`pangen` (default)**: XFOIL's PANGEN, line for line ([`repanel_by_curvature`], gated against the
//! `pangen_*` fixtures). The curvature along the splined input is smoothed, a fictitious curvature
//! is added at the trailing edge (CTERAT) and inside the refinement windows (CTRRAT), and the nodes
//! are placed so that `(1 + 6·CVPAR·curvature)·Δs` is equal on every panel.
//!
//! **`cosine`**: yFoil's own arc-length cosine spacing ([`repanel_cosine`]), no XFOIL equivalent,
//! its parameter warped by a power law set by `--cosine-te-bias`: 1 is a plain cosine, below 1
//! coarser at the trailing edge and finer at the leading edge, above 1 finer at the trailing
//! edge. It writes N + 1 nodes (historic behaviour, frozen by `tests/repanel_cosine_tests.rs`).
//! On the generators `cosine` is the analytic sampling at cosine chord stations and has no bias.
//!
//! Examples:
//! ```text
//! yfoil geom repanel aerofoil.dat -n 180
//! yfoil geom repanel aerofoil.dat -n 200 --curvature-bunching 1.5 --refine-upper 0.2,0.4
//! yfoil geom repanel aerofoil.json --method cosine --cosine-te-bias 0.2
//! yfoil geom repanel aerofoil.dat --panelling panelling.json
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
//!   "cm_ref": [0.25, 0.0],
//!   "x": [1.0, 0.8, 0.5, 0.2, 0.0, 0.2, 0.5, 0.8, 1.0],
//!   "y": [0.0, 0.02, 0.04, 0.03, 0.0, -0.03, -0.04, -0.02, 0.0],
//!   "generator": { "...": "how it was generated and panelled, when yFoil made it" }
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
//! - [`repanel`] with a [`PanelConfig`]: the one entry point — method ([`PangenConfig`] for XFOIL's
//!   PANGEN, [`CosineConfig`] for yFoil's cosine), trailing-edge treatment, provenance record
//! - [`repanel_by_curvature`]: XFOIL's PANGEN itself; [`repanel_cosine`]: yFoil's cosine spacing itself
//! - [`spline_derivatives`], [`spline_value`], [`spline_slope`], [`spline_second_derivative`]: cubic spline interpolation

mod airfoil;
mod io;
mod naca;
mod panel;
pub mod series;
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
    arc_coordinate, curvature, find_le, panel_foil, repanel, repanel_by_curvature, repanel_cosine, set_te_gap,
    solve_tridiagonal, spline_segmented, CosineConfig, PanelConfig, PanelMethod, PangenConfig, RepanelError, TeGap,
    PANGEN_BUFFER_NODES,
};
pub use series::{KarmanTrefftz, KarmanTrefftzError, MeanLine, Section, Series, SixSeriesFamily, ThicknessForm};
pub use spline::{spline_derivatives, spline_second_derivative, spline_slope, spline_value};
