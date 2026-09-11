//! NACA aerofoil generators: the 4- and 5-digit entry points kept for the fixture pipeline and
//! the CLI, and XFOIL's own `NACA4`/`NACA5` model (thickness applied vertically) for comparison.
//! Every NACA family is defined in [`super::series`]; the perpendicular-thickness generators here
//! delegate to it.

use super::airfoil::Geometry;
use super::panel::{repanel_by_curvature, PangenConfig};
use super::series::Section;

/// How a NACA section's thickness distribution is applied to its camber line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum Thickness {
    /// Perpendicular to the camber line: the NACA definition (yFoil's spacing)
    #[default]
    Perpendicular,
    /// Vertical, as XFOIL's NACA4/NACA5 apply it (`naca.f:62`), on XFOIL's 245-point buffer
    /// and then panelled with PANGEN to the requested count
    Vertical,
}

/// Translates XFOIL's `NACA4`.
///
/// Generate a NACA 4-digit airfoil
///
/// # Arguments
/// * `designation` - 4-digit string (e.g., "0012", "4412")
/// * `n_panels` - Number of panel points to generate
///
/// * `thickness` - perpendicular to the camber line (the NACA definition) or vertical (XFOIL's)
///
/// # Returns
/// Geometry with cosine-spaced points from TE around to TE
///
/// # Format
/// - First digit: maximum camber as percentage of chord
/// - Second digit: position of maximum camber in tenths of chord
/// - Last two digits: maximum thickness as percentage of chord
#[doc(alias = "NACA4")]
pub fn naca_4digit(designation: &str, n_panels: usize, thickness: Thickness) -> Result<Geometry, NacaError> {
    if thickness == Thickness::Vertical {
        let buffer = naca_4digit_vertical(designation)?;
        return Ok(repanel_by_curvature(&buffer, n_panels, &PangenConfig::default()));
    }
    if designation.len() != 4 {
        return Err(NacaError::InvalidDesignation(
            "NACA 4-digit designation must be exactly 4 characters".to_string(),
        ));
    }
    Ok(Section::from_designation(designation)?.geometry(n_panels))
}

/// Translates XFOIL's `NACA5`.
///
/// Generate a NACA 5-digit airfoil
///
/// # Arguments
/// * `designation` - 5-digit string (e.g., "23012", "23015")
/// * `n_panels` - Number of panel points to generate
///
/// # Format
/// - First digit: design lift coefficient * (2/3) * 10
/// - Second digit: position of max camber / 2 (in % chord / 10)
/// - Third digit: 0 = normal camber, 1 = reflex camber
/// - Last two digits: maximum thickness as percentage of chord
///
/// Common examples: 23012, 23015, 24112 (reflex)
#[doc(alias = "NACA5")]
pub fn naca_5digit(designation: &str, n_panels: usize, thickness: Thickness) -> Result<Geometry, NacaError> {
    if thickness == Thickness::Vertical {
        let buffer = naca_5digit_vertical(designation)?;
        return Ok(repanel_by_curvature(&buffer, n_panels, &PangenConfig::default()));
    }
    if designation.len() != 5 {
        return Err(NacaError::InvalidDesignation(
            "NACA 5-digit designation must be exactly 5 characters".to_string(),
        ));
    }
    Ok(Section::from_designation(designation)?.geometry(n_panels))
}

#[derive(thiserror::Error, Debug)]
pub enum NacaError {
    #[error("Invalid NACA designation: {0}")]
    InvalidDesignation(String),
}

/// XFOIL's own NACA generator (`naca.f`, `NACA4`/`NACA5`, called by the `NACA` command with
/// NSIDE = IQX/3 = 123 points per side). Reproduced as XFOIL does it, including the thickness
/// applied *vertically* (`YB = YC ± YT`) rather than perpendicular to the camber line, which
/// is CLAUDE.md's documented divergence from the NACA definition. The result is the 245-point
/// buffer airfoil XFOIL splines and then repanels with PANGEN (`repanel_by_curvature`); it takes no
/// panel count.
pub const XFOIL_NACA_NSIDE: usize = 123;

fn xfoil_naca_xx(nside: usize) -> Vec<f64> {
    // TE point bunching parameter
    let an: f64 = 1.5;
    let anp = an + 1.0;
    (1..=nside)
        .map(|i| {
            let frac = (i - 1) as f64 / (nside - 1) as f64;
            if i == nside {
                1.0
            } else {
                1.0 - anp * frac * (1.0 - frac).powf(an) - (1.0 - frac).powf(anp)
            }
        })
        .collect()
}

fn xfoil_naca_yt(xx: f64, t: f64) -> f64 {
    let x2 = xx * xx;
    let x3 = x2 * xx;
    let x4 = x2 * x2;
    (0.29690 * xx.sqrt() - 0.12600 * xx - 0.35160 * x2 + 0.28430 * x3 - 0.10150 * x4) * t / 0.20
}

fn xfoil_naca_assemble(xx: &[f64], yt: &[f64], yc: &[f64]) -> Geometry {
    let nside = xx.len();
    let mut x_c = Vec::with_capacity(2 * nside - 1);
    let mut y_c = Vec::with_capacity(2 * nside - 1);
    for i in (0..nside).rev() {
        x_c.push(xx[i]);
        y_c.push(yc[i] + yt[i]);
    }
    for i in 1..nside {
        x_c.push(xx[i]);
        y_c.push(yc[i] - yt[i]);
    }
    Geometry {
        cm_ref: [0.25, 0.0],
        x: x_c,
        y: y_c,
        generator: None,
    }
}

/// The provenance record of XFOIL's own NACA model: the divergence from the NACA definition
/// (thickness applied vertically, `docs/xfoil-known-issues.md` §6.1) is stated in the file
fn xfoil_naca_record(series: &str, designation: &str, buffer: &str) -> serde_json::Value {
    serde_json::json!({
        "yfoil": env!("CARGO_PKG_VERSION"),
        "series": series,
        "designation": format!("NACA {designation}"),
        "thickness_applied": "vertical",
        "buffer": buffer,
        "buffer_nodes": 2 * XFOIL_NACA_NSIDE - 1,
        "sharp_te": false,
        "references": ["jacobs1933", "jacobs1935"],
    })
}

/// `NACA4` as XFOIL runs it (vertical thickness, AN = 1.5 spacing, 2·NSIDE − 1 = 245 points).
#[doc(alias = "NACA4")]
pub fn naca_4digit_vertical(designation: &str) -> Result<Geometry, NacaError> {
    if designation.len() != 4 || !designation.chars().all(|c| c.is_ascii_digit()) {
        return Err(NacaError::InvalidDesignation(
            "NACA 4-digit designation must be exactly 4 digits".to_string(),
        ));
    }
    let ides: i64 = designation.parse().unwrap();
    let n4 = ides / 1000;
    let n3 = (ides - n4 * 1000) / 100;
    let n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;
    let m = n4 as f64 / 100.0;
    let p = n3 as f64 / 10.0;
    let t = (n2 * 10 + n1) as f64 / 100.0;

    let xx = xfoil_naca_xx(XFOIL_NACA_NSIDE);
    let yt: Vec<f64> = xx.iter().map(|&x| xfoil_naca_yt(x, t)).collect();
    let yc: Vec<f64> = xx
        .iter()
        .map(|&x| {
            if x < p {
                m / (p * p) * (2.0 * p * x - x * x)
            } else {
                m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x)
            }
        })
        .collect();
    let mut g = xfoil_naca_assemble(&xx, &yt, &yc);
    g.generator = Some(xfoil_naca_record("naca_4_digit", designation, "NACA4"));
    Ok(g)
}

/// `NACA5` as XFOIL runs it (210xx … 250xx camber lines by its M/C table, vertical thickness).
#[doc(alias = "NACA5")]
pub fn naca_5digit_vertical(designation: &str) -> Result<Geometry, NacaError> {
    if designation.len() != 5 || !designation.chars().all(|c| c.is_ascii_digit()) {
        return Err(NacaError::InvalidDesignation(
            "NACA 5-digit designation must be exactly 5 digits".to_string(),
        ));
    }
    let ides: i64 = designation.parse().unwrap();
    let n5 = ides / 10000;
    let n4 = (ides - n5 * 10000) / 1000;
    let n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    let n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    let n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;
    let n543 = 100 * n5 + 10 * n4 + n3;
    let (m, c) = match n543 {
        210 => (0.0580, 361.4),
        220 => (0.1260, 51.64),
        230 => (0.2025, 15.957),
        240 => (0.2900, 6.643),
        250 => (0.3910, 3.230),
        _ => {
            return Err(NacaError::InvalidDesignation(
                "Illegal 5-digit designation: first three digits must be 210, 220, ... 250".to_string(),
            ))
        }
    };
    let t = (n2 * 10 + n1) as f64 / 100.0;

    let xx = xfoil_naca_xx(XFOIL_NACA_NSIDE);
    let yt: Vec<f64> = xx.iter().map(|&x| xfoil_naca_yt(x, t)).collect();
    let yc: Vec<f64> = xx
        .iter()
        .map(|&x| {
            if x < m {
                (c / 6.0) * (x * x * x - 3.0 * m * (x * x) + m * m * (3.0 - m) * x)
            } else {
                (c / 6.0) * (m * m * m) * (1.0 - x)
            }
        })
        .collect();
    let mut g = xfoil_naca_assemble(&xx, &yt, &yc);
    g.generator = Some(xfoil_naca_record("naca_5_digit", designation, "NACA5"));
    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_naca_0012_generation() {
        let geom = naca_4digit("0012", 100, Thickness::Perpendicular).unwrap();

        // Should produce exactly the requested number of points
        assert_eq!(geom.x.len(), 100);

        // Should start and end at trailing edge (x=1.0)
        assert!((geom.x[0] - 1.0).abs() < 0.01);
        assert!((geom.x.last().unwrap() - 1.0).abs() < 0.01);

        // Symmetric airfoil: check that max thickness is approximately 12%
        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

        // Max thickness should be approximately 0.12 * 0.5 = 0.06 (half-thickness)
        assert_relative_eq!(max_y, 0.06, epsilon = 0.01);
        assert_relative_eq!(min_y, -0.06, epsilon = 0.01);
    }

    #[test]
    fn test_naca_4412_cambered() {
        let geom = naca_4digit("4412", 100, Thickness::Perpendicular).unwrap();

        // Cambered airfoil: upper surface should have more positive y values
        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

        // Upper surface should be thicker due to camber
        assert!(max_y > 0.08);
        assert!(min_y > -0.07);
    }

    #[test]
    fn test_thickness_distribution() {
        use super::super::series::four_digit_half_thickness as thickness_distribution;
        // At x=0.3 for 12% thick airfoil, thickness should be approximately maximum
        let yt = thickness_distribution(0.3, 0.12);
        assert!(yt > 0.05 && yt < 0.07);

        // At x=0, thickness should be 0
        let yt_le = thickness_distribution(0.0, 0.12);
        assert_relative_eq!(yt_le, 0.0, epsilon = 1e-10);

        // At x=1, thickness should be approximately 0.00126 (blunt TE, XFOIL-compatible)
        let yt_te = thickness_distribution(1.0, 0.12);
        assert!((yt_te - 0.00126).abs() < 0.0001);
    }

    #[test]
    fn test_naca_23012() {
        let geom = naca_5digit("23012", 100, Thickness::Perpendicular).unwrap();

        // The shared cosine spacing: exactly the requested count, both sides
        assert_eq!(geom.x.len(), 100);

        // Should start and end near trailing edge
        assert!(geom.x[0] > 0.95);
        assert!(geom.x.last().unwrap() > &0.95);

        // NACA 23012 is cambered - upper surface should be higher
        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

        // Check for positive camber (max_y should be larger in magnitude than |min_y|)
        assert!(max_y > -min_y * 0.8, "Expected positive camber");

        // Thickness should be approximately 12%
        let thickness = max_y - min_y;
        assert!(thickness > 0.10 && thickness < 0.14);
    }

    #[test]
    fn test_naca_23015() {
        let geom = naca_5digit("23015", 100, Thickness::Perpendicular).unwrap();

        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

        // Thickness should be approximately 15%
        let thickness = max_y - min_y;
        assert!(thickness > 0.13 && thickness < 0.17);
    }

    #[test]
    fn test_naca_5digit_different_cl() {
        // NACA 13012 has lower Cl than 23012
        let geom_low = naca_5digit("13012", 100, Thickness::Perpendicular).unwrap();
        let geom_high = naca_5digit("23012", 100, Thickness::Perpendicular).unwrap();

        // Higher Cl should have more camber (larger y values on upper surface)
        let max_y_low = geom_low.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let max_y_high = geom_high.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        assert!(max_y_high > max_y_low, "Higher Cl design should have more camber");
    }

    #[test]
    fn test_naca_5digit_symmetric() {
        // NACA 00012 - symmetric 5-digit (Cl = 0)
        let geom = naca_5digit("00012", 100, Thickness::Perpendicular).unwrap();

        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

        // Should be symmetric
        assert_relative_eq!(max_y, -min_y, epsilon = 0.001);
    }

    #[test]
    fn test_naca_5digit_invalid() {
        // Invalid length
        assert!(naca_5digit("2301", 100, Thickness::Perpendicular).is_err());
        assert!(naca_5digit("230123", 100, Thickness::Perpendicular).is_err());

        // Invalid camber position for non-reflex
        assert!(naca_5digit("26012", 100, Thickness::Perpendicular).is_err()); // p = 0.30 not in table
    }

    #[test]
    fn test_naca_5digit_reflex() {
        // NACA 23112 - reflex camber (third digit = 1)
        let geom = naca_5digit("23112", 100, Thickness::Perpendicular).unwrap();

        // Should generate valid geometry
        assert_eq!(geom.x.len(), 100);

        // Reflex camber should still have positive camber but trailing edge should curve up
        let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_y > 0.0);
    }
}
