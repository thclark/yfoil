use serde::Deserialize;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Deserialize)]
pub struct Geometry {
    pub reference: [f64; 2], // 2-element vector allocated onto stack - see https://stackoverflow.com/a/30263497/3556110
    pub x_c: Vec<f64>,
    pub y_c: Vec<f64>,
}

pub fn read_geometry_from_file<P: AsRef<Path>>(path: P) -> Result<Geometry, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Geometry`.
    let geometry = serde_json::from_reader(reader)?;

    // Return the `Geometry`
    Ok(geometry)
}
