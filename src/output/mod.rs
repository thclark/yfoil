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

pub use foil::{
    BlQuantity, BoundaryLayerOutput, Closures, FoilNodes, LaggedClosures, Primaries, SeparationKind, SeparationMarker,
    SideStations, StagnationMarker, TransitionMarker, WakeNodes,
};
pub use results::{
    AnalysisOutput, GeometryDistributions, GeometryInfo, GeometrySummary, PolarOutput, PolarPoint, PolarSummary,
    SurfaceDistributions,
};

#[cfg(feature = "plotting")]
pub use foil_plot::{
    check_same_geometry, plot_foil, plot_foil_png, plot_foil_svg, point_style, render_foil_svg, scale_factors,
    surface_offset, wake_band, CurveStyle, DashPattern, DesignPoint, FoilPlotConfig, ImageFormat, MarkerSet,
    OffsetScale, PanelStyle, Pt,
};
#[cfg(feature = "plotting")]
pub use plot::{
    parse_xfoil_cp_file, parse_xfoil_dat_file, parse_xfoil_dump_file, parse_xfoil_polar_file,
    parse_xfoil_polar_file_with_cm, parse_yfoil_polar_file, parse_yfoil_polar_file_with_cm, plot_analysis_png,
    plot_analysis_svg, plot_bl_comparison_svg, plot_cp_ue_comparison_svg, plot_geometry_comparison_svg,
    plot_polar_3panel_svg, plot_polar_comparison_svg, plot_polars_png, plot_polars_svg, stitch_polars,
    stitch_polars_with_cm, AnalysisPlotConfig, GeometryComparisonPlotConfig, Marker, PlotError, PolarData,
    PolarDataWithCm, PolarPlotConfig, PolarSeries, PolarsPlotConfig, XfoilCp, XfoilDump, YfoilBLDist,
};
