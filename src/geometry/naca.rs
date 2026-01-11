//! NACA airfoil generators
//!
//! Generate standard NACA airfoil profiles from their designation numbers.

use super::airfoil::Geometry;

/// Generate a NACA 4-digit airfoil
///
/// # Arguments
/// * `designation` - 4-digit string (e.g., "0012", "4412")
/// * `n_panels` - Number of panel points to generate
///
/// # Returns
/// Geometry with cosine-spaced points from TE around to TE
///
/// # Format
/// - First digit: maximum camber as percentage of chord
/// - Second digit: position of maximum camber in tenths of chord
/// - Last two digits: maximum thickness as percentage of chord
pub fn naca_4digit(designation: &str, n_panels: usize) -> Result<Geometry, NacaError> {
    if designation.len() != 4 {
        return Err(NacaError::InvalidDesignation(
            "NACA 4-digit designation must be exactly 4 characters".to_string(),
        ));
    }

    // Parse designation
    let m = designation[0..1]
        .parse::<f64>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid camber digit".to_string()))?
        / 100.0; // max camber
    let p = designation[1..2]
        .parse::<f64>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid camber position digit".to_string()))?
        / 10.0; // camber position
    let t = designation[2..4]
        .parse::<f64>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid thickness digits".to_string()))?
        / 100.0; // thickness

    // Generate cosine-spaced x coordinates
    let n_half = n_panels / 2;
    let mut x_upper = Vec::with_capacity(n_half + 1);
    let mut y_upper = Vec::with_capacity(n_half + 1);
    let mut x_lower = Vec::with_capacity(n_half + 1);
    let mut y_lower = Vec::with_capacity(n_half + 1);

    for i in 0..=n_half {
        let beta = std::f64::consts::PI * (i as f64) / (n_half as f64);
        let x = 0.5 * (1.0 - beta.cos());

        // Thickness distribution (modified for closed TE)
        let yt = thickness_distribution(x, t);

        // Camber line and its derivative
        let (yc, dyc_dx) = camber_line(x, m, p);

        // Angle of camber line
        let theta = dyc_dx.atan();

        // Upper and lower surface coordinates
        x_upper.push(x - yt * theta.sin());
        y_upper.push(yc + yt * theta.cos());
        x_lower.push(x + yt * theta.sin());
        y_lower.push(yc - yt * theta.cos());
    }

    // Combine: TE -> upper (reverse) -> LE -> lower -> TE
    // Upper surface goes from TE (x=1) to LE (x=0), so reverse it
    // Lower surface goes from LE (x=0) to TE (x=1)
    let mut x_c = Vec::with_capacity(2 * n_half + 1);
    let mut y_c = Vec::with_capacity(2 * n_half + 1);

    // Upper surface from TE to LE (reverse order, skip last point which is LE)
    for i in (1..=n_half).rev() {
        x_c.push(x_upper[i]);
        y_c.push(y_upper[i]);
    }

    // Lower surface from LE to TE (skip first point which is LE, already included)
    for i in 0..=n_half {
        x_c.push(x_lower[i]);
        y_c.push(y_lower[i]);
    }

    Ok(Geometry {
        reference: [0.25, 0.0], // Quarter chord
        x_c,
        y_c,
    })
}

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
pub fn naca_5digit(designation: &str, n_panels: usize) -> Result<Geometry, NacaError> {
    if designation.len() != 5 {
        return Err(NacaError::InvalidDesignation(
            "NACA 5-digit designation must be exactly 5 characters".to_string(),
        ));
    }

    // Parse designation
    let first = designation[0..1]
        .parse::<u32>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid first digit".to_string()))?;
    let second = designation[1..2]
        .parse::<u32>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid second digit".to_string()))?;
    let third = designation[2..3]
        .parse::<u32>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid third digit".to_string()))?;
    let thickness = designation[3..5]
        .parse::<f64>()
        .map_err(|_| NacaError::InvalidDesignation("Invalid thickness digits".to_string()))?
        / 100.0;

    // Design lift coefficient
    let cl = (first as f64) * 0.15; // Cl = first_digit * 3/20

    // Position of maximum camber
    let p = (second as f64) * 0.05; // p = second_digit / 20

    // Check for reflex camber
    let reflex = third == 1;

    if p == 0.0 && cl != 0.0 {
        return Err(NacaError::InvalidDesignation(
            "Invalid camber position (second digit cannot be 0 with non-zero lift)".to_string(),
        ));
    }

    // Get camber line coefficients
    let (r, k1, k2_k1) = get_5digit_coefficients(cl, p, reflex)?;

    // Generate cosine-spaced x coordinates
    let n_half = n_panels / 2;
    let mut x_upper = Vec::with_capacity(n_half + 1);
    let mut y_upper = Vec::with_capacity(n_half + 1);
    let mut x_lower = Vec::with_capacity(n_half + 1);
    let mut y_lower = Vec::with_capacity(n_half + 1);

    for i in 0..=n_half {
        let beta = std::f64::consts::PI * (i as f64) / (n_half as f64);
        let x = 0.5 * (1.0 - beta.cos());

        // Thickness distribution (same as 4-digit, modified for closed TE)
        let yt = thickness_distribution(x, thickness);

        // Camber line and its derivative
        let (yc, dyc_dx) = camber_line_5digit(x, r, k1, k2_k1, reflex);

        // Angle of camber line
        let theta = dyc_dx.atan();

        // Upper and lower surface coordinates
        x_upper.push(x - yt * theta.sin());
        y_upper.push(yc + yt * theta.cos());
        x_lower.push(x + yt * theta.sin());
        y_lower.push(yc - yt * theta.cos());
    }

    // Combine: TE -> upper (reverse) -> LE -> lower -> TE
    let mut x_c = Vec::with_capacity(2 * n_half + 1);
    let mut y_c = Vec::with_capacity(2 * n_half + 1);

    // Upper surface from TE to LE (reverse order, skip last point which is LE)
    for i in (1..=n_half).rev() {
        x_c.push(x_upper[i]);
        y_c.push(y_upper[i]);
    }

    // Lower surface from LE to TE
    for i in 0..=n_half {
        x_c.push(x_lower[i]);
        y_c.push(y_lower[i]);
    }

    Ok(Geometry {
        reference: [0.25, 0.0],
        x_c,
        y_c,
    })
}

/// Get coefficients for NACA 5-digit mean camber line
///
/// Returns (r, k1, k2/k1) where r is the position where camber meets the
/// straight section, k1 is the camber multiplier, and k2/k1 is the ratio
/// for reflex cambers.
fn get_5digit_coefficients(cl: f64, p: f64, reflex: bool) -> Result<(f64, f64, f64), NacaError> {
    // Standard 5-digit camber line coefficients
    // These are tabulated values for specific p positions
    // p = 0.05, 0.10, 0.15, 0.20, 0.25

    if cl == 0.0 {
        // Symmetric airfoil
        return Ok((0.0, 0.0, 0.0));
    }

    // Lookup table for standard (non-reflex) 5-digit cambers
    // Format: (p, r, k1) - k1 is for Cl = 0.3 (first digit = 2)
    let standard_coeffs = [
        (0.05, 0.0580, 361.400),
        (0.10, 0.1260, 51.640),
        (0.15, 0.2025, 15.957),
        (0.20, 0.2900, 6.643),
        (0.25, 0.3910, 3.230),
    ];

    // Lookup table for reflex 5-digit cambers
    let reflex_coeffs = [
        (0.10, 0.1300, 51.990, 0.000764),
        (0.15, 0.2170, 15.793, 0.00677),
        (0.20, 0.3180, 6.520, 0.0303),
        (0.25, 0.4410, 3.191, 0.1355),
    ];

    // Find closest p value and interpolate if needed
    let (r, k1_base, k2_k1) = if reflex {
        // Find matching reflex coefficients
        let mut found = None;
        for &(pi, ri, k1i, k2_k1i) in &reflex_coeffs {
            if (pi - p).abs() < 0.001 {
                found = Some((ri, k1i, k2_k1i));
                break;
            }
        }
        found.ok_or_else(|| {
            NacaError::InvalidDesignation(format!(
                "Reflex camber position {} not supported. Use 0.10, 0.15, 0.20, or 0.25",
                p
            ))
        })?
    } else {
        // Find matching standard coefficients
        let mut found = None;
        for &(pi, ri, k1i) in &standard_coeffs {
            if (pi - p).abs() < 0.001 {
                found = Some((ri, k1i, 0.0));
                break;
            }
        }
        found.ok_or_else(|| {
            NacaError::InvalidDesignation(format!(
                "Camber position {} not supported. Use 0.05, 0.10, 0.15, 0.20, or 0.25",
                p
            ))
        })?
    };

    // Scale k1 for actual Cl (base values are for Cl = 0.3)
    let k1 = k1_base * (cl / 0.3);

    Ok((r, k1, k2_k1))
}

/// NACA 5-digit mean camber line
///
/// Returns (y_c, dy_c/dx) at given x/c coordinate
fn camber_line_5digit(x: f64, r: f64, k1: f64, k2_k1: f64, reflex: bool) -> (f64, f64) {
    if k1 == 0.0 {
        return (0.0, 0.0);
    }

    if reflex {
        // Reflex camber line (three regions)
        let k2 = k1 * k2_k1;
        if x < r {
            let yc = (k1 / 6.0)
                * (x.powi(3) - 3.0 * r * x.powi(2) + r.powi(2) * (3.0 - r) * x);
            let dyc = (k1 / 6.0)
                * (3.0 * x.powi(2) - 6.0 * r * x + r.powi(2) * (3.0 - r));
            (yc, dyc)
        } else {
            let yc = (k1 * r.powi(3) / 6.0) * (1.0 - x)
                - (k2 / 6.0)
                    * (x.powi(3) - 3.0 * r * x.powi(2)
                        + 3.0 * r.powi(2) * x
                        - r.powi(3));
            let dyc = -(k1 * r.powi(3) / 6.0)
                - (k2 / 6.0) * (3.0 * x.powi(2) - 6.0 * r * x + 3.0 * r.powi(2));
            (yc, dyc)
        }
    } else {
        // Standard camber line (two regions)
        if x < r {
            let yc = (k1 / 6.0)
                * (x.powi(3) - 3.0 * r * x.powi(2) + r.powi(2) * (3.0 - r) * x);
            let dyc = (k1 / 6.0)
                * (3.0 * x.powi(2) - 6.0 * r * x + r.powi(2) * (3.0 - r));
            (yc, dyc)
        } else {
            let yc = (k1 * r.powi(3) / 6.0) * (1.0 - x);
            let dyc = -(k1 * r.powi(3) / 6.0);
            (yc, dyc)
        }
    }
}

/// NACA 4-digit thickness distribution
///
/// Returns half-thickness at given x/c coordinate
fn thickness_distribution(x: f64, t: f64) -> f64 {
    // Standard NACA 4-digit thickness equation
    // Modified last coefficient for closed trailing edge (0.1036 instead of 0.1015)
    let yt = 5.0
        * t
        * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3)
            - 0.1036 * x.powi(4));
    yt
}

/// NACA 4-digit mean camber line
///
/// Returns (y_c, dy_c/dx) at given x/c coordinate
fn camber_line(x: f64, m: f64, p: f64) -> (f64, f64) {
    if m == 0.0 || p == 0.0 {
        // Symmetric airfoil
        return (0.0, 0.0);
    }

    let (yc, dyc) = if x < p {
        let yc = m / (p * p) * (2.0 * p * x - x * x);
        let dyc = 2.0 * m / (p * p) * (p - x);
        (yc, dyc)
    } else {
        let yc = m / ((1.0 - p).powi(2)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
        let dyc = 2.0 * m / ((1.0 - p).powi(2)) * (p - x);
        (yc, dyc)
    };

    (yc, dyc)
}

#[derive(thiserror::Error, Debug)]
pub enum NacaError {
    #[error("Invalid NACA designation: {0}")]
    InvalidDesignation(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_naca_0012_symmetric() {
        let geom = naca_4digit("0012", 100).unwrap();

        // Should have points on both sides
        assert!(geom.x_c.len() > 100);

        // Should start and end near trailing edge
        assert!(geom.x_c[0] > 0.95);
        assert!(geom.x_c.last().unwrap() > &0.95);

        // Symmetric airfoil: check that max thickness is approximately 12%
        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        // Max thickness should be approximately 0.12 * 0.5 = 0.06 (half-thickness)
        assert_relative_eq!(max_y, 0.06, epsilon = 0.01);
        assert_relative_eq!(min_y, -0.06, epsilon = 0.01);
    }

    #[test]
    fn test_naca_4412_cambered() {
        let geom = naca_4digit("4412", 100).unwrap();

        // Cambered airfoil: upper surface should have more positive y values
        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        // Upper surface should be thicker due to camber
        assert!(max_y > 0.08);
        assert!(min_y > -0.07);
    }

    #[test]
    fn test_thickness_distribution() {
        // At x=0.3 for 12% thick airfoil, thickness should be approximately maximum
        let yt = thickness_distribution(0.3, 0.12);
        assert!(yt > 0.05 && yt < 0.07);

        // At x=0, thickness should be 0
        let yt_le = thickness_distribution(0.0, 0.12);
        assert_relative_eq!(yt_le, 0.0, epsilon = 1e-10);

        // At x=1, thickness should be approximately 0 (closed TE)
        let yt_te = thickness_distribution(1.0, 0.12);
        assert!(yt_te.abs() < 0.002);
    }

    #[test]
    fn test_naca_23012() {
        let geom = naca_5digit("23012", 100).unwrap();

        // Should have points on both sides
        assert!(geom.x_c.len() > 100);

        // Should start and end near trailing edge
        assert!(geom.x_c[0] > 0.95);
        assert!(geom.x_c.last().unwrap() > &0.95);

        // NACA 23012 is cambered - upper surface should be higher
        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        // Check for positive camber (max_y should be larger in magnitude than |min_y|)
        assert!(max_y > -min_y * 0.8, "Expected positive camber");

        // Thickness should be approximately 12%
        let thickness = max_y - min_y;
        assert!(thickness > 0.10 && thickness < 0.14);
    }

    #[test]
    fn test_naca_23015() {
        let geom = naca_5digit("23015", 100).unwrap();

        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        // Thickness should be approximately 15%
        let thickness = max_y - min_y;
        assert!(thickness > 0.13 && thickness < 0.17);
    }

    #[test]
    fn test_naca_5digit_different_cl() {
        // NACA 13012 has lower Cl than 23012
        let geom_low = naca_5digit("13012", 100).unwrap();
        let geom_high = naca_5digit("23012", 100).unwrap();

        // Higher Cl should have more camber (larger y values on upper surface)
        let max_y_low = geom_low.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let max_y_high = geom_high.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        assert!(
            max_y_high > max_y_low,
            "Higher Cl design should have more camber"
        );
    }

    #[test]
    fn test_naca_5digit_symmetric() {
        // NACA 00012 - symmetric 5-digit (Cl = 0)
        let geom = naca_5digit("00012", 100).unwrap();

        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_y = geom.y_c.iter().cloned().fold(f64::INFINITY, f64::min);

        // Should be symmetric
        assert_relative_eq!(max_y, -min_y, epsilon = 0.001);
    }

    #[test]
    fn test_naca_5digit_invalid() {
        // Invalid length
        assert!(naca_5digit("2301", 100).is_err());
        assert!(naca_5digit("230123", 100).is_err());

        // Invalid camber position for non-reflex
        assert!(naca_5digit("26012", 100).is_err()); // p = 0.30 not in table
    }

    #[test]
    fn test_naca_5digit_reflex() {
        // NACA 23112 - reflex camber (third digit = 1)
        let geom = naca_5digit("23112", 100).unwrap();

        // Should generate valid geometry
        assert!(geom.x_c.len() > 100);

        // Reflex camber should still have positive camber but trailing edge should curve up
        let max_y = geom.y_c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_y > 0.0);
    }
}
