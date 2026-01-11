//! Airfoil coordinate representation and validation

use serde::{Deserialize, Serialize};

const MIN_PANELS: usize = 100;
const MAX_PANELS: usize = 250;

/// Raw airfoil geometry from input file
///
/// Coordinates are normalized by chord (x/c, y/c).
/// Points ordered: TE -> upper surface -> LE -> lower surface -> TE
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Geometry {
    /// Reference point for moment calculation [x/c, y/c]
    pub reference: [f64; 2],
    /// X-coordinates normalized by chord
    pub x_c: Vec<f64>,
    /// Y-coordinates normalized by chord
    pub y_c: Vec<f64>,
}

/// Paneled airfoil ready for aerodynamic analysis
///
/// This struct contains all geometric quantities needed by the panel method
/// and boundary layer solver.
#[derive(Debug, Clone)]
pub struct PaneledAirfoil {
    /// Panel node x-coordinates (TE -> upper -> LE -> lower -> TE)
    pub x: Vec<f64>,
    /// Panel node y-coordinates
    pub y: Vec<f64>,
    /// Arc length parameter (spline parameter)
    pub s: Vec<f64>,
    /// Spline derivatives dx/ds
    pub xp: Vec<f64>,
    /// Spline derivatives dy/ds
    pub yp: Vec<f64>,
    /// Normal vector x-components (pointing outward)
    pub nx: Vec<f64>,
    /// Normal vector y-components (pointing outward)
    pub ny: Vec<f64>,
    /// Panel angles (angle of panel tangent from horizontal)
    pub apanel: Vec<f64>,
    /// Number of panel nodes
    pub n: usize,
    /// Leading edge arc length parameter
    pub sle: f64,
    /// Leading edge node index
    pub le_index: usize,
    /// Chord length
    pub chord: f64,
    /// Whether trailing edge is sharp (zero thickness)
    pub sharp_te: bool,
    /// Reference point for moment calculation
    pub reference: [f64; 2],
}

#[derive(thiserror::Error, Debug)]
pub enum InvalidGeometryError {
    #[error("The x_c and y_c arrays must be the same size! Currently {0} and {1} elements.")]
    MismatchedDimension(usize, usize),

    #[error("The number of panels is too few! Minimum is {0}.")]
    TooFewPanels(usize),

    #[error("The number of panels is too many! Maximum is {0}.")]
    TooManyPanels(usize),

    #[error("The maximum value of x_c is < 0.95 or > 1.05, suggesting the input geometry is not a normalised aerofoil")]
    MaxXExtent,

    #[error("The minimum value in x_c array is < -0.05 or > 0.05, suggesting the input geometry is not a normalised aerofoil")]
    MinXExtent,

    #[error("At least one value in y_c is < -1.0 or > 1.0. YFoil is not intended for bluff bodies!")]
    MaxYExtent,

    #[error("All points are either above or below y_c=0. Perhaps your y_c has a nonzero y offset?")]
    OffsetY,

    #[error("First and last points in the x_c arrays are not near the trailing edge")]
    NotStartingAtTrailingEdge,
}

impl Geometry {
    /// Validate the geometry for basic sanity checks
    pub fn validate(&self) -> Result<(), InvalidGeometryError> {
        use InvalidGeometryError::*;

        let nx = self.x_c.len();
        let ny = self.y_c.len();

        // Find the maximum and minimum values
        let max_x = self.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let max_y = self.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_x = self.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
        let min_y = self.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        let first_x = *self.x_c.first().unwrap_or(&0.0);
        let last_x = *self.x_c.last().unwrap_or(&0.0);

        // Extremely basic panel quantity and aerofoil location / normalisation checks
        if nx != ny {
            Err(MismatchedDimension(nx, ny))
        } else if nx <= MIN_PANELS {
            Err(TooFewPanels(MIN_PANELS))
        } else if nx > MAX_PANELS + 1 {
            Err(TooManyPanels(MAX_PANELS))
        } else if max_x < 0.95 || max_x > 1.05 {
            Err(MaxXExtent)
        } else if min_x < -0.05 || min_x > 0.05 {
            Err(MinXExtent)
        } else if max_y > 1.0 || min_y < -1.0 {
            Err(MaxYExtent)
        } else if max_y < 0.0 || min_y > 0.0 {
            Err(OffsetY)
        } else if first_x < 0.95 || last_x < 0.95 {
            Err(NotStartingAtTrailingEdge)
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometry_validation_mismatched_dimensions() {
        let geom = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![0.0; 150],
            y_c: vec![0.0; 151],
        };
        assert!(matches!(
            geom.validate(),
            Err(InvalidGeometryError::MismatchedDimension(150, 151))
        ));
    }

    #[test]
    fn test_geometry_validation_too_few_panels() {
        let geom = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![0.0; 50],
            y_c: vec![0.0; 50],
        };
        assert!(matches!(
            geom.validate(),
            Err(InvalidGeometryError::TooFewPanels(_))
        ));
    }
}
