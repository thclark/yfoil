//! Airfoil coordinate representation and validation

use serde::{Deserialize, Serialize};

// XFOIL itself accepts any NB > 1 (ABCOPY) up to IQX-5; 100 was an arbitrary product floor that
// blocked small validation and fuzz cases. 20 is the smallest count PANGEN handles sensibly.
const MIN_PANELS: usize = 20;
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

    #[error(
        "The maximum value of x_c is < 0.95 or > 1.05, suggesting the input geometry is not a normalised aerofoil"
    )]
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
    /// Sharpen the trailing edge by moving first and last points to coincide.
    ///
    /// This creates a sharp TE geometry where nodes 0 and n-1 are at the same
    /// location, which is appropriate for airfoils with zero TE thickness.
    /// The sharp TE requires special handling in the panel method solver.
    ///
    /// # Returns
    /// A new geometry with the TE nodes moved to their midpoint.
    pub fn sharpen(&self) -> Geometry {
        let n = self.x_c.len();
        if n < 2 {
            return self.clone();
        }

        // Calculate midpoint of current TE
        let x_te = 0.5 * (self.x_c[0] + self.x_c[n - 1]);
        let y_te = 0.5 * (self.y_c[0] + self.y_c[n - 1]);

        let mut x_c = self.x_c.clone();
        let mut y_c = self.y_c.clone();

        // Move both TE points to coincide at the midpoint
        x_c[0] = x_te;
        y_c[0] = y_te;
        x_c[n - 1] = x_te;
        y_c[n - 1] = y_te;

        Geometry {
            reference: self.reference,
            x_c,
            y_c,
        }
    }

    /// Open the trailing edge by creating a small gap.
    ///
    /// This creates a blunt TE geometry where nodes 0 and n-1 have a small
    /// vertical separation, which is appropriate for airfoils with finite TE
    /// thickness. The blunt TE uses standard flow tangency equations in the
    /// panel method (no curvature extrapolation).
    ///
    /// # Arguments
    /// * `gap` - The vertical gap to create at TE (typical value: 0.002 for 0.2% chord)
    ///
    /// # Returns
    /// A new geometry with a blunt trailing edge.
    pub fn blunten(&self, gap: f64) -> Geometry {
        let n = self.x_c.len();
        if n < 2 {
            return self.clone();
        }

        // Calculate midpoint of current TE
        let x_te = 0.5 * (self.x_c[0] + self.x_c[n - 1]);
        let y_te = 0.5 * (self.y_c[0] + self.y_c[n - 1]);

        let mut x_c = self.x_c.clone();
        let mut y_c = self.y_c.clone();

        // Separate the TE points vertically (node 0 is upper surface, node n-1 is lower)
        // For standard airfoil ordering: TE -> upper -> LE -> lower -> TE
        x_c[0] = x_te;
        y_c[0] = y_te + gap / 2.0; // Upper surface TE moves up
        x_c[n - 1] = x_te;
        y_c[n - 1] = y_te - gap / 2.0; // Lower surface TE moves down

        Geometry {
            reference: self.reference,
            x_c,
            y_c,
        }
    }

    /// Check if the trailing edge is sharp (closed).
    ///
    /// Returns true if the TE gap is less than 0.01% of chord.
    pub fn is_sharp_te(&self) -> bool {
        let n = self.x_c.len();
        if n < 2 {
            return false;
        }
        let te_gap = ((self.x_c[0] - self.x_c[n - 1]).powi(2) + (self.y_c[0] - self.y_c[n - 1]).powi(2)).sqrt();
        // Use same threshold as panel.rs: 0.0001 * chord
        // For normalized airfoil, chord ≈ 1.0
        te_gap < 0.0001
    }

    /// Get the trailing edge gap as a fraction of chord.
    pub fn te_gap(&self) -> f64 {
        let n = self.x_c.len();
        if n < 2 {
            return 0.0;
        }
        ((self.x_c[0] - self.x_c[n - 1]).powi(2) + (self.y_c[0] - self.y_c[n - 1]).powi(2)).sqrt()
    }

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
    use approx::assert_relative_eq;

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
            x_c: vec![0.0; 10],
            y_c: vec![0.0; 10],
        };
        assert!(matches!(geom.validate(), Err(InvalidGeometryError::TooFewPanels(_))));
    }

    #[test]
    fn test_sharpen_closes_te_gap() {
        // Create a simple geometry with open TE
        let geom = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],       // TE at x=1, LE at x=0
            y_c: vec![0.01, 0.05, 0.0, -0.05, -0.01], // Gap of 0.02 at TE
        };

        assert!(!geom.is_sharp_te());
        assert_relative_eq!(geom.te_gap(), 0.02, epsilon = 1e-10);

        let sharpened = geom.sharpen();

        assert!(sharpened.is_sharp_te());
        assert_relative_eq!(sharpened.te_gap(), 0.0, epsilon = 1e-10);

        // TE should be at midpoint
        assert_relative_eq!(sharpened.x_c[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(sharpened.y_c[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(sharpened.x_c[sharpened.x_c.len() - 1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(sharpened.y_c[sharpened.y_c.len() - 1], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_blunten_opens_te_gap() {
        // Create a simple geometry with closed TE
        let geom = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y_c: vec![0.0, 0.05, 0.0, -0.05, 0.0], // Closed TE
        };

        assert!(geom.is_sharp_te());
        assert_relative_eq!(geom.te_gap(), 0.0, epsilon = 1e-10);

        let blunted = geom.blunten(0.02);

        assert!(!blunted.is_sharp_te());
        assert_relative_eq!(blunted.te_gap(), 0.02, epsilon = 1e-10);

        // Upper TE should move up, lower should move down
        assert_relative_eq!(blunted.y_c[0], 0.01, epsilon = 1e-10);
        assert_relative_eq!(blunted.y_c[blunted.y_c.len() - 1], -0.01, epsilon = 1e-10);
    }

    #[test]
    fn test_sharpen_then_blunten_preserves_midpoint() {
        let geom = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y_c: vec![0.01, 0.05, 0.0, -0.05, -0.01],
        };

        let sharpened = geom.sharpen();
        let blunted = sharpened.blunten(0.02);

        // The midpoint y should be at 0 (original gap midpoint)
        let mid_y = 0.5 * (blunted.y_c[0] + blunted.y_c[blunted.y_c.len() - 1]);
        assert_relative_eq!(mid_y, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_is_sharp_te_threshold() {
        // Gap of 0.00005 (below threshold of 0.0001)
        let sharp = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y_c: vec![0.000025, 0.05, 0.0, -0.05, -0.000025],
        };
        assert!(sharp.is_sharp_te());

        // Gap of 0.0002 (above threshold of 0.0001)
        let blunt = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y_c: vec![0.0001, 0.05, 0.0, -0.05, -0.0001],
        };
        assert!(!blunt.is_sharp_te());
    }
}
