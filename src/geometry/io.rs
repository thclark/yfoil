//! File I/O for airfoil geometry (JSON and .dat formats)

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use super::airfoil::{Geometry, InvalidGeometryError};

#[derive(thiserror::Error, Debug)]
pub enum GeometryReadError {
    #[error("Failed to read the geometry file: {0}")]
    FileRead(#[from] std::io::Error),

    #[error("Failed to parse geometry file contents: {0}")]
    FileParse(#[from] serde_json::Error),

    #[error("Invalid .dat file format: {0}")]
    DatParse(String),

    #[error(transparent)]
    InvalidGeometry(#[from] InvalidGeometryError),
}

/// Read airfoil geometry from a JSON file
///
/// The JSON file should contain:
/// - `reference`: [x, y] reference point for moments
/// - `x_c`: array of x/c coordinates
/// - `y_c`: array of y/c coordinates
///
/// Points should be ordered: TE -> upper surface -> LE -> lower surface -> TE
pub fn read_geometry_from_file<P: AsRef<Path>>(path: P) -> Result<Geometry, GeometryReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let geometry: Geometry = serde_json::from_reader(reader)?;

    // Validate geometry
    geometry.validate()?;

    Ok(geometry)
}

/// Read airfoil geometry from a Selig/XFOIL .dat file
///
/// Format: First line is the airfoil name (optional), subsequent lines are x y coordinates.
/// Coordinates can be separated by spaces or tabs.
///
/// Supports both Selig format (single wrap-around) and Lednicer format (separate upper/lower).
/// For Lednicer format, coordinates are combined into wrap-around order.
pub fn read_dat_file<P: AsRef<Path>>(path: P) -> Result<(String, Geometry), GeometryReadError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // First line is the airfoil name
    let name = lines
        .next()
        .ok_or_else(|| GeometryReadError::DatParse("Empty file".to_string()))??
        .trim()
        .to_string();

    let mut x_coords = Vec::new();
    let mut y_coords = Vec::new();
    let mut line_num = 1;

    for line_result in lines {
        line_num += 1;
        let line = line_result?;
        let trimmed = line.trim();

        // Skip empty lines
        if trimmed.is_empty() {
            continue;
        }

        // Parse x y coordinate pair
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 2 {
            let x: f64 = parts[0].parse().map_err(|_| {
                GeometryReadError::DatParse(format!(
                    "Invalid x coordinate '{}' on line {}",
                    parts[0], line_num
                ))
            })?;

            let y: f64 = parts[1].parse().map_err(|_| {
                GeometryReadError::DatParse(format!(
                    "Invalid y coordinate '{}' on line {}",
                    parts[1], line_num
                ))
            })?;

            x_coords.push(x);
            y_coords.push(y);
        }
        // Skip lines with fewer than 2 numbers (might be headers or separators)
    }

    if x_coords.is_empty() {
        return Err(GeometryReadError::DatParse(
            "No coordinates found in file".to_string(),
        ));
    }

    // Detect Lednicer format: if x values go 0->1 then 0->1 again
    // In this case, we need to reverse the upper surface and combine
    let geometry = if is_lednicer_format(&x_coords) {
        convert_lednicer_to_selig(&x_coords, &y_coords)?
    } else {
        Geometry {
            reference: [0.25, 0.0], // Default quarter-chord reference
            x_c: x_coords,
            y_c: y_coords,
        }
    };

    Ok((name, geometry))
}

/// Check if coordinates are in Lednicer format (LE->TE twice)
fn is_lednicer_format(x: &[f64]) -> bool {
    if x.len() < 6 {
        return false;
    }

    // Lednicer format: x goes from ~0 to ~1 (upper surface), then 0 to 1 again (lower)
    // Key distinction from Selig format:
    // - Lednicer STARTS near LE (x < 0.1) and first increases toward TE
    // - Selig STARTS near TE (x > 0.8) and first decreases toward LE
    //
    // Additionally, in Lednicer format we see x reach ~1 then jump back to ~0
    // In Selig format, x goes 1 -> 0 -> 1 smoothly (just passes through LE)

    // Check starting point: Selig starts near TE, Lednicer starts near LE
    let starts_near_le = x[0] < 0.2;
    let starts_near_te = x[0] > 0.8;

    if starts_near_te {
        // Selig format - x decreases from TE towards LE
        return false;
    }

    if !starts_near_le {
        // Neither clear format - default to Selig
        return false;
    }

    // Starts near LE - could be Lednicer
    // Check if x increases towards TE then restarts near LE
    for i in 1..x.len() {
        // Look for the restart: x was near TE, now near LE
        if x[i - 1] > 0.8 && x[i] < 0.1 {
            return true;
        }
    }

    false
}

/// Convert Lednicer format to Selig (wrap-around) format
fn convert_lednicer_to_selig(x: &[f64], y: &[f64]) -> Result<Geometry, GeometryReadError> {
    // Find the split point (where upper ends and lower begins)
    let mut split_idx = 0;
    for i in 1..x.len() {
        if x[i] < x[i - 1] - 0.3 {
            split_idx = i;
            break;
        }
    }

    if split_idx == 0 {
        return Err(GeometryReadError::DatParse(
            "Could not find split point in Lednicer format".to_string(),
        ));
    }

    // Upper surface: indices 0..split_idx (LE to TE), needs to be reversed
    // Lower surface: indices split_idx..end (LE to TE), keep as is
    // Combined: TE (upper) -> LE -> TE (lower)

    let mut x_c = Vec::with_capacity(x.len());
    let mut y_c = Vec::with_capacity(y.len());

    // Upper surface reversed (TE to LE)
    for i in (0..split_idx).rev() {
        x_c.push(x[i]);
        y_c.push(y[i]);
    }

    // Lower surface (LE to TE), skip first point if it duplicates LE
    let start_lower = if (x[split_idx] - x[split_idx - 1]).abs() < 0.001 {
        split_idx + 1
    } else {
        split_idx
    };

    for i in start_lower..x.len() {
        x_c.push(x[i]);
        y_c.push(y[i]);
    }

    Ok(Geometry {
        reference: [0.25, 0.0],
        x_c,
        y_c,
    })
}

/// Write airfoil geometry to a JSON file
pub fn write_geometry_to_json<P: AsRef<Path>>(
    geometry: &Geometry,
    path: P,
) -> Result<(), std::io::Error> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, geometry)?;
    Ok(())
}

/// Write airfoil geometry to a Selig/XFOIL .dat file
///
/// Format: First line is the airfoil name, subsequent lines are x y coordinates.
pub fn write_dat_file<P: AsRef<Path>>(
    geometry: &Geometry,
    name: &str,
    path: P,
) -> Result<(), std::io::Error> {
    let mut file = File::create(path)?;

    // Write name
    writeln!(file, "{}", name)?;

    // Write coordinates with full double precision for numerical consistency
    for i in 0..geometry.x_c.len() {
        writeln!(file, " {:22.16}  {:22.16}", geometry.x_c[i], geometry.y_c[i])?;
    }

    Ok(())
}

/// Detect file format from extension and read geometry
pub fn read_geometry_auto<P: AsRef<Path>>(path: P) -> Result<Geometry, GeometryReadError> {
    let path_ref = path.as_ref();
    let ext = path_ref
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "json" => read_geometry_from_file(path),
        "dat" => {
            let (_, geom) = read_dat_file(path)?;
            Ok(geom)
        }
        _ => {
            // Try JSON first, then .dat
            if let Ok(geom) = read_geometry_from_file(&path) {
                Ok(geom)
            } else {
                let (_, geom) = read_dat_file(path)?;
                Ok(geom)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_write_dat_roundtrip() {
        // Create a simple geometry
        let original = Geometry {
            reference: [0.25, 0.0],
            x_c: vec![1.0, 0.5, 0.0, 0.5, 1.0],
            y_c: vec![0.0, 0.05, 0.0, -0.05, 0.0],
        };

        // Write to temp file
        let temp_file = NamedTempFile::new().unwrap();
        write_dat_file(&original, "Test Airfoil", temp_file.path()).unwrap();

        // Read back
        let (name, loaded) = read_dat_file(temp_file.path()).unwrap();

        assert_eq!(name, "Test Airfoil");
        assert_eq!(loaded.x_c.len(), original.x_c.len());

        for i in 0..original.x_c.len() {
            assert!((loaded.x_c[i] - original.x_c[i]).abs() < 1e-5);
            assert!((loaded.y_c[i] - original.y_c[i]).abs() < 1e-5);
        }
    }

    #[test]
    fn test_read_selig_format() {
        // Create a Selig format file (wrap-around order)
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "NACA 0012").unwrap();
        writeln!(temp_file, "  1.000000   0.000000").unwrap();
        writeln!(temp_file, "  0.500000   0.050000").unwrap();
        writeln!(temp_file, "  0.000000   0.000000").unwrap();
        writeln!(temp_file, "  0.500000  -0.050000").unwrap();
        writeln!(temp_file, "  1.000000   0.000000").unwrap();
        temp_file.flush().unwrap();

        let (name, geom) = read_dat_file(temp_file.path()).unwrap();

        assert_eq!(name, "NACA 0012");
        assert_eq!(geom.x_c.len(), 5);
        assert!((geom.x_c[0] - 1.0).abs() < 1e-6);
        assert!((geom.x_c[2] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_skip_empty_lines() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "Test").unwrap();
        writeln!(temp_file, "").unwrap();
        writeln!(temp_file, "  1.0  0.0").unwrap();
        writeln!(temp_file, "").unwrap();
        writeln!(temp_file, "  0.0  0.0").unwrap();
        temp_file.flush().unwrap();

        let (_, geom) = read_dat_file(temp_file.path()).unwrap();

        assert_eq!(geom.x_c.len(), 2);
    }

    #[test]
    fn test_is_lednicer_format() {
        // Selig format: x goes 1 -> 0 -> 1
        let selig_x = vec![1.0, 0.5, 0.0, 0.5, 1.0];
        assert!(!is_lednicer_format(&selig_x));

        // Lednicer format: x goes 0 -> 1, then 0 -> 1 again
        let lednicer_x = vec![0.0, 0.5, 1.0, 0.0, 0.5, 1.0];
        assert!(is_lednicer_format(&lednicer_x));
    }

    #[test]
    fn test_read_xfoil_fortran_notation() {
        // XFOIL outputs coordinates in Fortran scientific notation
        // e.g., "0.1260000E-02" instead of "0.00126"
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "NACA 0012").unwrap();
        writeln!(temp_file, "    1.000000      0.1260000E-02").unwrap();
        writeln!(temp_file, "   0.9916796      0.2421450E-02").unwrap();
        writeln!(temp_file, "   0.9803692      0.3981369E-02").unwrap();
        writeln!(temp_file, "   0.0000260      0.9100000E-03").unwrap();
        writeln!(temp_file, "   0.9803692     -0.3981369E-02").unwrap();
        writeln!(temp_file, "   0.9916796     -0.2421450E-02").unwrap();
        writeln!(temp_file, "    1.000000     -0.1260000E-02").unwrap();
        temp_file.flush().unwrap();

        let (name, geom) = read_dat_file(temp_file.path()).unwrap();

        assert_eq!(name, "NACA 0012");
        assert_eq!(geom.x_c.len(), 7);

        // Check first point (TE upper)
        assert!(
            (geom.x_c[0] - 1.0).abs() < 1e-6,
            "x[0] = {} (expected 1.0)",
            geom.x_c[0]
        );
        assert!((geom.y_c[0] - 0.00126).abs() < 1e-8);

        // Check second point
        assert!((geom.x_c[1] - 0.9916796).abs() < 1e-6);
        assert!((geom.y_c[1] - 0.00242145).abs() < 1e-8);

        // Check LE point
        assert!((geom.x_c[3] - 0.000026).abs() < 1e-8);
        assert!((geom.y_c[3] - 0.00091).abs() < 1e-8);

        // Check TE lower (last point)
        assert!((geom.x_c[6] - 1.0).abs() < 1e-6);
        assert!((geom.y_c[6] - (-0.00126)).abs() < 1e-8);
    }

    #[test]
    fn test_xfoil_paneled_geometry_properties() {
        // Test that XFOIL-paneled coordinates produce correct geometric properties
        // This simulates loading a 160-panel NACA 0012 from XFOIL
        use crate::geometry::panel::create_paneled_airfoil;

        // Create a representative subset of XFOIL-paneled NACA 0012
        // (Full 160 panels would be too verbose, using 20 points for test)
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "NACA 0012 (XFOIL paneled)").unwrap();
        // Upper surface from TE to LE
        writeln!(temp_file, "  1.000000   0.001260").unwrap();
        writeln!(temp_file, "  0.950000   0.008000").unwrap();
        writeln!(temp_file, "  0.850000   0.017500").unwrap();
        writeln!(temp_file, "  0.700000   0.030000").unwrap();
        writeln!(temp_file, "  0.500000   0.042000").unwrap();
        writeln!(temp_file, "  0.300000   0.047000").unwrap();
        writeln!(temp_file, "  0.150000   0.044000").unwrap();
        writeln!(temp_file, "  0.050000   0.030500").unwrap();
        writeln!(temp_file, "  0.010000   0.015200").unwrap();
        writeln!(temp_file, "  0.000000   0.000000").unwrap();
        // Lower surface from LE to TE
        writeln!(temp_file, "  0.010000  -0.015200").unwrap();
        writeln!(temp_file, "  0.050000  -0.030500").unwrap();
        writeln!(temp_file, "  0.150000  -0.044000").unwrap();
        writeln!(temp_file, "  0.300000  -0.047000").unwrap();
        writeln!(temp_file, "  0.500000  -0.042000").unwrap();
        writeln!(temp_file, "  0.700000  -0.030000").unwrap();
        writeln!(temp_file, "  0.850000  -0.017500").unwrap();
        writeln!(temp_file, "  0.950000  -0.008000").unwrap();
        writeln!(temp_file, "  1.000000  -0.001260").unwrap();
        temp_file.flush().unwrap();

        let (name, geom) = read_dat_file(temp_file.path()).unwrap();
        assert_eq!(name, "NACA 0012 (XFOIL paneled)");
        assert_eq!(geom.x_c.len(), 19);

        // Create paneled airfoil and verify properties
        let paneled = create_paneled_airfoil(&geom);

        // Chord should be approximately 1.0
        assert!(
            (paneled.chord - 1.0).abs() < 0.01,
            "Chord should be ~1.0, got {}",
            paneled.chord
        );

        // LE should be at index 9 (the 0,0 point)
        assert_eq!(paneled.le_index, 9, "LE should be at index 9");

        // Arc lengths should be monotonically increasing
        for i in 1..paneled.s.len() {
            assert!(
                paneled.s[i] > paneled.s[i - 1],
                "Arc length should be monotonic at index {}",
                i
            );
        }

        // Normal vectors should have unit length
        for i in 0..paneled.n {
            let mag = (paneled.nx[i].powi(2) + paneled.ny[i].powi(2)).sqrt();
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "Normal at {} should be unit length",
                i
            );
        }
    }
}
