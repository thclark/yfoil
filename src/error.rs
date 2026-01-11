//! Unified error types for yfoil

use thiserror::Error;

use crate::geometry::{GeometryReadError, InvalidGeometryError};

/// Top-level error type for yfoil operations
#[derive(Error, Debug)]
pub enum YfoilError {
    /// Geometry-related errors
    #[error("Geometry error: {0}")]
    Geometry(#[from] GeometryReadError),

    /// Invalid geometry
    #[error("Invalid geometry: {0}")]
    InvalidGeometry(#[from] InvalidGeometryError),

    /// Solver convergence failure
    #[error("Solver failed to converge after {iterations} iterations (residual: {residual:.2e})")]
    Convergence { iterations: usize, residual: f64 },

    /// Invalid flow conditions
    #[error("Invalid flow conditions: {message}")]
    InvalidConditions { message: String },

    /// Numerical instability
    #[error("Numerical error: {message}")]
    Numerical { message: String },

    /// Matrix singularity
    #[error("Singular matrix in {context}")]
    SingularMatrix { context: String },

    /// I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// JSON error
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
}

/// Result type alias for yfoil operations
pub type Result<T> = std::result::Result<T, YfoilError>;
