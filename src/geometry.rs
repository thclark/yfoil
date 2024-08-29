use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

const MIN_PANELS: usize = 100;
const MAX_PANELS: usize = 250;

#[derive(Deserialize, Debug)]
pub struct Geometry {
    pub reference: [f64; 2], // 2-element vector allocated onto stack - see https://stackoverflow.com/a/30263497/3556110
    pub x_c: Vec<f64>,
    pub y_c: Vec<f64>,
}


#[derive(thiserror::Error, Debug)]
pub enum GeometryReadError {
    #[error("Failed to read the geometry file (missing or corrupted?)")]
    FileReadError(#[from] std::io::Error),

    #[error("Failed to parse geometry file contents (is it a valid JSON file?)")]
    FileParseError(#[from] serde_json::Error),

    #[error(transparent)]
    InvalidGeometryError(#[from] InvalidGeometryError),

}

pub fn read_geometry_from_file<P: AsRef<Path>>(path: P) -> Result<Geometry, GeometryReadError> {

    use GeometryReadError::*;

    // Open the file in read-only mode with buffered read
    let file = File::open(path).map_err(FileReadError)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Geometry`.
    let geometry = serde_json::from_reader(reader).map_err(FileParseError)?;

    // Validate the geometry by passing a reference
    validate_geometry(&geometry).map_err(InvalidGeometryError)?;

    // Return the `Geometry`
    Ok(geometry)
}

#[derive(thiserror::Error, Debug)]
pub enum InvalidGeometryError {
    #[error("The x_c and y_c arrays must be the same size! Currently {0} and {1} elements.")]
    MismatchedDimensionError(usize, usize),
    #[error("The number of panels is too few! Minimum is {0}.")]
    TooFewPanelsError(usize),
    #[error("The number of panels is too many! Maximum is {0}.")]
    TooManyPanelsError(usize),
    #[error("The maximum value of x_c is < 0.95 or > 1.05, suggesting the input geometry is not a normalised aerofoil (max x_c should approximately coincide with the trailing edge, [x=1,y=0])")]
    MaxXExtentError(),
    #[error("The minimum value in x_c array is < -0.05 or > 0.05, suggesting the input geometry is not a normalised aerofoil (min x_c should approximately coincide with the leading edge, [x=0,y=0])")]
    MinXExtentError(),
    #[error("At least one value in y_c is < -1.0 or > 1.0. YFoil is not intended for bluff bodies! Typically the aerofoil width should be less than half its chord length. Perhaps y_c is not normalised correctly?")]
    MaxYExtentError(),
    #[error("All points are either above or below y_c=0. YFoil is not suitable for sharp, extreme camber shapes like you might get in jet stator rows! Typically the aerofoil's camber line should deviate from its chord by a small amount. Perhaps your y_c has a nonzero y offset, or the aerofoil thickness isn't correct?")]
    OffsetYError(),
    #[error("First and last points in the x_c arrays are not near the trailing edge. Points should be ordered from TE (Trailing Edge), over the top in the -x direction, through LE (Leading Edge) in the -y direction, around the underside in +ve x direction to the TE.")]
    BeginAndEndAtTrailingEdgeError(),
}


fn validate_geometry(geometry: &Geometry) -> Result<(), InvalidGeometryError> {

    use InvalidGeometryError::*;

    let nx = geometry.x_c.len();
    let ny = geometry.y_c.len();

    // Find the maximum and minimum values
    let max_x = geometry.x_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let max_y = geometry.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_x = geometry.x_c.iter().cloned().fold(f64::INFINITY, f64::min);
    let min_y = geometry.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

    let first_x = *geometry.x_c.first().unwrap();
    let last_x = *geometry.x_c.last().unwrap();

    println!("Maximum values x_c, y_c: {}, {}", max_x, max_y);
    println!("Minimum values x_c, y_c: {}, {}", min_x, min_y);

    // Extremely basic panel quantity and aerofoil location / normalisation checks
    if nx != ny {
        Err(MismatchedDimensionError(nx, ny))
    } else if nx <= MIN_PANELS {
        Err(TooFewPanelsError(MIN_PANELS))
    } else if nx > MAX_PANELS + 1 {
        Err(TooManyPanelsError(MAX_PANELS))
    } else if max_x < 0.95 || max_x > 1.05 {
        Err(MaxXExtentError())
    } else if min_x < -0.05 || min_x > 0.05 {
        Err(MinXExtentError())
    } else if max_y > 1.0 || min_y < -1.0 {
        Err(MaxYExtentError())
    } else if max_y < 0.0 || min_y > 0.0 {
        Err(OffsetYError())
    } else if first_x < 0.95 || last_x < 0.95 {
        Err(BeginAndEndAtTrailingEdgeError())
    } else {
        Ok(())
    }
}
