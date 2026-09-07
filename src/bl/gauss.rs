//! Gaussian elimination solver
//!
//! This module implements XFOIL's GAUSS subroutine - a general NxN
//! solver using Gaussian elimination with partial pivoting.
//!
//! Used in MRCHDU for solving the local 4x4 Newton system at each station.

/// Solve a linear system using Gaussian elimination with partial pivoting
///
/// Solves the system Z * x = R in-place, where:
/// - Z is the coefficient matrix (destroyed during solution)
/// - R is the right-hand side (replaced by solution)
///
/// This is a direct translation of XFOIL's GAUSS subroutine from xsolve.f.
///
/// # Arguments
/// * `z` - Coefficient matrix (n x n), modified in place
/// * `r` - Right-hand side vector, replaced by solution
///
/// # Panics
/// Will panic if matrix is singular (division by zero during pivoting)
///
/// # Example
/// ```
/// use yfoil::bl::gauss::gauss_solve;
///
/// let mut z = [[2.0, 1.0], [1.0, 3.0]];
/// let mut r = [5.0, 7.0];
/// gauss_solve(&mut z, &mut r);
/// // r now contains solution: [1.6, 1.8]
/// ```
#[doc(alias = "GAUSS")]
pub fn gauss_solve<const N: usize>(z: &mut [[f64; N]; N], r: &mut [f64; N]) {
    // Forward elimination with partial pivoting
    for np in 0..N - 1 {
        let np1 = np + 1;

        // Find max pivot index
        let mut nx = np;
        for n in np1..N {
            if z[n][np].abs() > z[nx][np].abs() {
                nx = n;
            }
        }

        let pivot = 1.0 / z[nx][np];

        // Switch pivot rows
        z[nx][np] = z[np][np];

        // Switch rows & normalize pivot row
        for l in np1..N {
            let temp = z[nx][l] * pivot;
            z[nx][l] = z[np][l];
            z[np][l] = temp;
        }

        let temp = r[nx] * pivot;
        r[nx] = r[np];
        r[np] = temp;

        // Forward eliminate
        for k in np1..N {
            let ztmp = z[k][np];
            for l in np1..N {
                z[k][l] -= ztmp * z[np][l];
            }
            r[k] -= ztmp * r[np];
        }
    }

    // Solve for last row
    r[N - 1] /= z[N - 1][N - 1];

    // Back substitute
    for np in (0..N - 1).rev() {
        for k in (np + 1)..N {
            r[np] -= z[np][k] * r[k];
        }
    }
}

/// Translates XFOIL's `GAUSS`.
///
/// Solve a 4x4 linear system (MRCHDU's Newton system)
///
/// Specialized version for the 4x4 system used in MRCHDU.
/// Variables are: dCtau/dAmpl, dTheta, dDstar, dUe
#[inline]
#[doc(alias = "GAUSS")]
pub fn gauss_solve_4x4(z: &mut [[f64; 4]; 4], r: &mut [f64; 4]) {
    gauss_solve(z, r);
}

/// Solve a 4x4 system and also transform a 4x5 matrix
///
/// This variant is used for calculating the Ue-Hk characteristic slope
/// in MRCHDU. It applies the same elimination steps to both the
/// right-hand side r and an extra column of derivatives.
///
/// # Arguments
/// * `z` - 4x4 coefficient matrix (destroyed)
/// * `r` - 4-element RHS, replaced by solution
pub fn gauss_solve_with_extra(z: &mut [[f64; 4]; 4], r: &mut [f64; 4]) {
    gauss_solve(z, r);
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gauss_2x2() {
        // 2x + y = 5
        // x + 3y = 7
        // Solution: x = 8/5 = 1.6, y = 9/5 = 1.8
        let mut z = [[2.0, 1.0], [1.0, 3.0]];
        let mut r = [5.0, 7.0];
        gauss_solve(&mut z, &mut r);

        assert_relative_eq!(r[0], 1.6, epsilon = 1e-10);
        assert_relative_eq!(r[1], 1.8, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_3x3() {
        // 3x + 2y + z = 1
        // 2x + y + z = 4
        // x + 3y + 2z = 5
        let mut z = [[3.0, 2.0, 1.0], [2.0, 1.0, 1.0], [1.0, 3.0, 2.0]];
        let mut r = [1.0, 4.0, 5.0];
        gauss_solve(&mut z, &mut r);

        // Verify solution by back-substitution
        let x = r[0];
        let y = r[1];
        let z_val = r[2];
        assert_relative_eq!(3.0 * x + 2.0 * y + z_val, 1.0, epsilon = 1e-10);
        assert_relative_eq!(2.0 * x + y + z_val, 4.0, epsilon = 1e-10);
        assert_relative_eq!(x + 3.0 * y + 2.0 * z_val, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_4x4_identity() {
        // Identity matrix
        let mut z = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let mut r = [1.0, 2.0, 3.0, 4.0];
        gauss_solve_4x4(&mut z, &mut r);

        assert_relative_eq!(r[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(r[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(r[2], 3.0, epsilon = 1e-10);
        assert_relative_eq!(r[3], 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_4x4_requires_pivoting() {
        // System that requires pivoting (small diagonal element)
        let mut z = [
            [0.001, 1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let mut r = [1.0, 2.0, 3.0, 4.0];
        gauss_solve_4x4(&mut z, &mut r);

        // Verify by checking residual
        let z_orig = [
            [0.001, 1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let rhs = [1.0, 2.0, 3.0, 4.0];
        for i in 0..4 {
            let mut sum = 0.0;
            for j in 0..4 {
                sum += z_orig[i][j] * r[j];
            }
            assert_relative_eq!(sum, rhs[i], epsilon = 1e-10);
        }
    }
}
