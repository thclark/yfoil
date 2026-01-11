//! Pressure coefficient calculation with compressibility corrections
//!
//! Provides functions for computing pressure coefficients from surface
//! velocities, including the Karman-Tsien compressibility correction
//! for subsonic flow.

/// Calculate incompressible pressure coefficient from surface velocity
///
/// Cp = 1 - (q/q_∞)²
///
/// For unit freestream q_∞ = 1, this simplifies to:
/// Cp = 1 - q²
#[inline]
pub fn cp_incompressible(velocity: f64) -> f64 {
    1.0 - velocity * velocity
}

/// Apply Karman-Tsien compressibility correction to pressure coefficient
///
/// The Karman-Tsien rule corrects incompressible Cp for subsonic compressibility:
///
/// Cp_corrected = Cp_inc / (β + (M²/(1+β)) * (Cp_inc/2))
///
/// where β = √(1 - M²)
///
/// Valid for M < 1.0 (subsonic flow)
///
/// # Arguments
/// * `cp_inc` - Incompressible pressure coefficient
/// * `mach` - Freestream Mach number
///
/// # Returns
/// Compressibility-corrected pressure coefficient
pub fn cp_karman_tsien(cp_inc: f64, mach: f64) -> f64 {
    if mach < 0.001 {
        // Effectively incompressible
        return cp_inc;
    }

    if mach >= 1.0 {
        // Beyond validity of Karman-Tsien - return uncorrected
        // (In practice, should not be using this method for M >= 1)
        return cp_inc;
    }

    let m2 = mach * mach;
    let beta = (1.0 - m2).sqrt();

    // Karman-Tsien formula
    cp_inc / (beta + m2 / (1.0 + beta) * cp_inc / 2.0)
}

/// Calculate pressure coefficient array from surface velocities
///
/// # Arguments
/// * `velocities` - Surface velocities (γ values from panel method)
/// * `mach` - Freestream Mach number (0 for incompressible)
///
/// # Returns
/// Vector of pressure coefficients at each point
pub fn calculate_cp(velocities: &[f64], mach: f64) -> Vec<f64> {
    velocities
        .iter()
        .map(|&v| {
            let cp_inc = cp_incompressible(v);
            cp_karman_tsien(cp_inc, mach)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_cp_incompressible_stagnation() {
        // At stagnation point, v=0, Cp=1
        assert_relative_eq!(cp_incompressible(0.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cp_incompressible_freestream() {
        // At freestream velocity, v=1, Cp=0
        assert_relative_eq!(cp_incompressible(1.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_cp_incompressible_accelerated() {
        // Accelerated flow, v>1, Cp<0 (suction)
        let cp = cp_incompressible(1.2);
        assert!(cp < 0.0);
        assert_relative_eq!(cp, 1.0 - 1.44, epsilon = 1e-10);
    }

    #[test]
    fn test_karman_tsien_zero_mach() {
        // At M=0, correction should return original Cp
        let cp_inc = -0.5;
        assert_relative_eq!(cp_karman_tsien(cp_inc, 0.0), cp_inc, epsilon = 1e-10);
    }

    #[test]
    fn test_karman_tsien_increases_magnitude() {
        // Compressibility should increase magnitude of negative Cp (suction)
        let cp_inc = -0.5;
        let cp_corr = cp_karman_tsien(cp_inc, 0.5);

        // Corrected Cp should be more negative (higher suction)
        assert!(cp_corr < cp_inc);
    }

    #[test]
    fn test_karman_tsien_positive_cp() {
        // For positive Cp (compression), correction increases value
        let cp_inc = 0.5;
        let cp_corr = cp_karman_tsien(cp_inc, 0.5);

        // At M=0.5, beta = sqrt(1-0.25) = sqrt(0.75) ≈ 0.866
        // Denominator < 1, so Cp_corr > Cp_inc
        assert!(cp_corr > cp_inc);
    }

    #[test]
    fn test_calculate_cp_array() {
        let velocities = vec![0.0, 1.0, 1.2, 0.8];
        let cp = calculate_cp(&velocities, 0.0);

        assert_eq!(cp.len(), 4);
        assert_relative_eq!(cp[0], 1.0, epsilon = 1e-10); // stagnation
        assert_relative_eq!(cp[1], 0.0, epsilon = 1e-10); // freestream
        assert!(cp[2] < 0.0); // accelerated
        assert!(cp[3] > 0.0); // decelerated
    }
}
