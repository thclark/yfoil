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

    // In Lednicer format, x goes from ~0 to ~1, then restarts from ~0 to ~1
    // Look for a restart (large decrease in x where next point is near LE)
    for i in 1..x.len() {
        // Large jump backwards (x decreases by more than 0.3)
        // And the new point is near LE (x < 0.1)
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

    // Write coordinates
    for i in 0..geometry.x_c.len() {
        writeln!(file, " {:10.6}  {:10.6}", geometry.x_c[i], geometry.y_c[i])?;
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
}
