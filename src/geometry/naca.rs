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
    // To produce exactly n_panels points total, we need n_panels/2 points per surface.
    // XFOIL's PANE command:
    // - Places nodes at exactly x=1.0 (TE)
    // - Straddles the LE (no node at exactly x=0)
    let n_half = n_panels / 2;
    let mut x_upper = Vec::with_capacity(n_half);
    let mut y_upper = Vec::with_capacity(n_half);
    let mut x_lower = Vec::with_capacity(n_half);
    let mut y_lower = Vec::with_capacity(n_half);

    for i in 0..n_half {
        // Cosine distribution that places node at x=1 (TE) and straddles x=0 (LE).
        // After reversal, the ordering will be: TE (x=1) -> near-LE (x≈0).
        // XFOIL places station 1 exactly at x=1.0.
        //
        // Without half-cell offset: beta = π*i/(n_half-1) gives x=0 at i=0 and x=1 at i=n_half-1.
        // But we want to avoid x=0 exactly, so we use a slight offset at the LE end only:
        // - i=0: beta = 0.5π/(n_half-1+0.5), x ≈ 0.0001 (near LE)
        // - i=n_half-1: beta = π, x = 1.0 (exact TE)
        let beta = std::f64::consts::PI * (i as f64 + 0.5) / (n_half as f64 - 0.5);
        let beta = beta.min(std::f64::consts::PI); // Cap at π for last point
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

    // Combine: TE -> upper (reverse) -> near-LE upper -> near-LE lower -> lower -> TE
    // With half-cell offset, we have two near-LE points that straddle the actual LE:
    //   upper[0] at (x_small, +y) and lower[0] at (x_small, -y)
    // Both points must be included (like XFOIL's PANE does).
    // Total: n_half + n_half = n_panels points
    let mut x_c = Vec::with_capacity(n_panels);
    let mut y_c = Vec::with_capacity(n_panels);

    // Upper surface from TE to near-LE (n_half points: indices n_half-1, n_half-2, ..., 0)
    for i in (0..n_half).rev() {
        x_c.push(x_upper[i]);
        y_c.push(y_upper[i]);
    }

    // Lower surface from near-LE to TE (n_half points: indices 0, 1, ..., n_half-1)
    for i in 0..n_half {
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
        // Use half-cell offset to avoid putting a node at exactly x=0 (LE).
        // XFOIL's PANE command creates panels that straddle the LE, not pass through it.
        // With a node at exact x=0, gamma=0 there and the BL fails to converge.
        let beta = std::f64::consts::PI * (i as f64 + 0.5) / (n_half as f64 + 1.0);
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

    // Combine: TE -> upper (reverse) -> near-LE upper -> near-LE lower -> lower -> TE
    // With half-cell offset, we have two near-LE points that straddle the actual LE:
    //   upper[0] at (x_small, +y) and lower[0] at (x_small, -y)
    // Both points must be included (like XFOIL's PANE does).
    let mut x_c = Vec::with_capacity(2 * n_half + 2);
    let mut y_c = Vec::with_capacity(2 * n_half + 2);

    // Upper surface from TE to near-LE (include ALL points including i=0)
    for i in (0..=n_half).rev() {
        x_c.push(x_upper[i]);
        y_c.push(y_upper[i]);
    }

    // Lower surface from near-LE to TE (include ALL points including i=0)
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
            let yc = (k1 / 6.0) * (x.powi(3) - 3.0 * r * x.powi(2) + r.powi(2) * (3.0 - r) * x);
            let dyc = (k1 / 6.0) * (3.0 * x.powi(2) - 6.0 * r * x + r.powi(2) * (3.0 - r));
            (yc, dyc)
        } else {
            let yc = (k1 * r.powi(3) / 6.0) * (1.0 - x)
                - (k2 / 6.0) * (x.powi(3) - 3.0 * r * x.powi(2) + 3.0 * r.powi(2) * x - r.powi(3));
            let dyc = -(k1 * r.powi(3) / 6.0) - (k2 / 6.0) * (3.0 * x.powi(2) - 6.0 * r * x + 3.0 * r.powi(2));
            (yc, dyc)
        }
    } else {
        // Standard camber line (two regions)
        if x < r {
            let yc = (k1 / 6.0) * (x.powi(3) - 3.0 * r * x.powi(2) + r.powi(2) * (3.0 - r) * x);
            let dyc = (k1 / 6.0) * (3.0 * x.powi(2) - 6.0 * r * x + r.powi(2) * (3.0 - r));
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
    // Original coefficient -0.1015 gives blunt trailing edge (XFOIL default)
    // For NACA 0012: half-thickness at TE = 0.00126, gap = 0.00252
    5.0 * t * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3) - 0.1015 * x.powi(4))
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
        reference: [0.25, 0.0],
        x_c,
        y_c,
    }
}

/// `NACA4` as XFOIL runs it (vertical thickness, AN = 1.5 spacing, 2·NSIDE − 1 = 245 points).
pub fn naca_4digit_xfoil(designation: &str) -> Result<Geometry, NacaError> {
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
    Ok(xfoil_naca_assemble(&xx, &yt, &yc))
}

/// `NACA5` as XFOIL runs it (210xx … 250xx camber lines by its M/C table, vertical thickness).
pub fn naca_5digit_xfoil(designation: &str) -> Result<Geometry, NacaError> {
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
    Ok(xfoil_naca_assemble(&xx, &yt, &yc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_naca_0012_symmetric() {
        let geom = naca_4digit("0012", 100).unwrap();

        // Should produce exactly the requested number of points
        assert_eq!(geom.x_c.len(), 100);

        // Should start and end at trailing edge (x=1.0)
        assert!((geom.x_c[0] - 1.0).abs() < 0.01);
        assert!((geom.x_c.last().unwrap() - 1.0).abs() < 0.01);

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

        // At x=1, thickness should be approximately 0.00126 (blunt TE, XFOIL-compatible)
        let yt_te = thickness_distribution(1.0, 0.12);
        assert!((yt_te - 0.00126).abs() < 0.0001);
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

        assert!(max_y_high > max_y_low, "Higher Cl design should have more camber");
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
