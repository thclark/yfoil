//! Output and visualization
//!
//! This module handles result serialization and optional plotting.
//!
//! # Components
//! - `results` - Result structs (JSON serializable)
//! - `foil` - Geometry, wake and per-station BL output of an operating point
//! - `figures` - Shared conventions of yFoil's ironlab figures (feature-gated)
//! - `foil_plot` / `plot` - The `yfoil plot` figures (feature-gated)

mod foil;
mod results;

#[cfg(feature = "plotting")]
pub mod figures;
#[cfg(feature = "plotting")]
mod foil_plot;
#[cfg(feature = "plotting")]
mod plot;

pub use foil::{
    BlQuantity, BoundaryLayerOutput, Closures, FoilNodes, LaggedClosures, Primaries, SeparationKind, SeparationMarker,
    SideStations, StagnationMarker, TransitionMarker, WakeNodes,
};
pub use results::{
    AnalysisOutput, Diagnostics, GeometryDistributions, GeometryInfo, GeometrySummary, PointRecord, PointStatus,
    PointValues, PolarOutput, PolarSummary, Reason, SurfaceDistributions,
};

#[cfg(feature = "plotting")]
pub use foil_plot::{
    check_same_geometry, equal_aspect_limits, foil_figure, plot_foil, point_style, quantity_symbol, scale_factors,
    surface_offset, wake_band, CurveStyle, DesignPoint, FoilPlotConfig, MarkerSet, OffsetScale, PanelStyle, Pt,
    DASH_CYCLE, PLOT_AREA_MM,
};
#[cfg(feature = "plotting")]
pub use plot::{
    analysis_figure, plot_analysis, plot_polars, polars_figure, save_figure, AnalysisPlotConfig, PlotError,
    PolarSeries, PolarsPlotConfig, SeriesMarker, STUDY,
};
