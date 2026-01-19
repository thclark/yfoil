//! BLSOLV: Custom block solver for coupled viscous-inviscid Newton system
//!
//! This is a direct translation of XFOIL's BLSOLV subroutine (xsolve.f).
//! The algorithm solves a block-tridiagonal system with dense mass defect coupling.
//!
//! ## System Structure
//!
//! ```text
//!   A  |  |  .  |  |  .  |    d       R       S
//!   B  A  |  .  |  |  .  |    d       R       S
//!   |  B  A  .  |  |  .  |    d       R       S
//!   .  .  .  .  |  |  .  |    d   =   R - dRe S
//!   |  |  |  B  A  |  .  |    d       R       S
//!   |  Z  |  |  B  A  .  |    d       R       S
//!   .  .  .  .  .  .  .  |    d       R       S
//!   |  |  |  |  |  |  B  A    d       R       S
//! ```
//!
//! Where:
//! - A, B, Z: 3x2 blocks containing linearized BL equation coefficients
//! - |: 3x1 vectors containing mass defect influence coefficients on Ue
//! - d: 3x1 unknown vectors (Newton deltas for Ctau/Ampl, Theta, Mass)
//! - R: 3x1 residual vectors
//! - S: 3x1 Re influence vectors

/// Input matrices for BLSOLV
#[derive(Debug, Clone)]
pub struct BlsolvInput {
    /// Number of BL stations
    pub nsys: usize,
    /// Diagonal blocks VA[iv][k][l] - 3 equations, 2 columns
    pub va: Vec<[[f64; 2]; 3]>,
    /// Sub-diagonal blocks VB[iv][k][l]
    pub vb: Vec<[[f64; 2]; 3]>,
    /// RHS/solution VDEL[iv][k][l] - column 0 is residual, column 1 is Re sensitivity
    pub vdel: Vec<[[f64; 2]; 3]>,
    /// Mass defect coupling VM[iv][j][k] - coupling from station j to station iv, equation k
    pub vm: Vec<Vec<[f64; 3]>>,
    /// TE coupling block VZ[k][l] (optional, only used at trailing edge)
    pub vz: [[f64; 2]; 3],
    /// Index where upper surface ends at TE (0-based)
    pub ivte1: Option<usize>,
    /// System index for start of wake (lower surface TE + 1)
    pub ivz: Option<usize>,
    /// Acceleration parameter for sparse elimination
    pub vaccel: f64,
    /// Total arc length S(N) - S(1) for VACC scaling (optional)
    /// If provided, VACC2 and VACC3 are scaled by 2.0 / arc_length
    pub arc_length: Option<f64>,
}

impl Default for BlsolvInput {
    fn default() -> Self {
        Self {
            nsys: 0,
            va: Vec::new(),
            vb: Vec::new(),
            vdel: Vec::new(),
            vm: Vec::new(),
            vz: [[0.0; 2]; 3],
            ivte1: None,
            ivz: None,
            vaccel: 0.01,
            arc_length: None,
        }
    }
}

/// Solve the coupled BL Newton system using XFOIL's BLSOLV algorithm
///
/// This is a direct translation of XFOIL's BLSOLV subroutine.
/// The algorithm performs block Gaussian elimination with special handling
/// for the dense mass defect coupling (VM matrix).
///
/// # Arguments
/// * `input` - Mutable reference to the input matrices (modified in place during solve)
///
/// # Returns
/// The solution is stored in input.vdel after the solve completes.
pub fn blsolv(input: &mut BlsolvInput) {
    let nsys = input.nsys;
    if nsys == 0 {
        return;
    }

    // Compute acceleration thresholds
    // XFOIL: VACC1 = VACCEL, VACC2 = VACC3 = VACCEL * 2.0 / (S(N) - S(1))
    let vacc1 = input.vaccel;
    let vacc_scale = match input.arc_length {
        Some(arc_len) if arc_len > 0.0 => 2.0 / arc_len,
        _ => 1.0, // Default to no scaling if arc length not provided
    };
    let vacc2 = input.vaccel * vacc_scale;
    let vacc3 = input.vaccel * vacc_scale;

    // Forward sweep: IV = 0 to NSYS-1
    for iv in 0..nsys {
        let ivp = iv + 1;

        // ====== Invert VA(IV) block ======

        // Normalize first row by VA(1,1)
        let pivot = 1.0 / input.va[iv][0][0];
        input.va[iv][0][1] *= pivot;
        for l in iv..nsys {
            input.vm[iv][l][0] *= pivot;
        }
        input.vdel[iv][0][0] *= pivot;
        input.vdel[iv][0][1] *= pivot;

        // Eliminate lower first column in VA block (rows 2,3)
        for k in 1..3 {
            let vtmp = input.va[iv][k][0];
            input.va[iv][k][1] -= vtmp * input.va[iv][0][1];
            for l in iv..nsys {
                input.vm[iv][l][k] -= vtmp * input.vm[iv][l][0];
            }
            input.vdel[iv][k][0] -= vtmp * input.vdel[iv][0][0];
            input.vdel[iv][k][1] -= vtmp * input.vdel[iv][0][1];
        }

        // Normalize second row by VA(2,2)
        let pivot = 1.0 / input.va[iv][1][1];
        for l in iv..nsys {
            input.vm[iv][l][1] *= pivot;
        }
        input.vdel[iv][1][0] *= pivot;
        input.vdel[iv][1][1] *= pivot;

        // Eliminate lower second column in VA block (row 3)
        let vtmp = input.va[iv][2][1];
        for l in iv..nsys {
            input.vm[iv][l][2] -= vtmp * input.vm[iv][l][1];
        }
        input.vdel[iv][2][0] -= vtmp * input.vdel[iv][1][0];
        input.vdel[iv][2][1] -= vtmp * input.vdel[iv][1][1];

        // Normalize third row by VM(3,IV,IV) - the diagonal mass coupling
        let pivot = 1.0 / input.vm[iv][iv][2];
        for l in ivp..nsys {
            input.vm[iv][l][2] *= pivot;
        }
        input.vdel[iv][2][0] *= pivot;
        input.vdel[iv][2][1] *= pivot;

        // Eliminate upper third column in VA block (rows 1,2)
        let vtmp1 = input.vm[iv][iv][0];
        let vtmp2 = input.vm[iv][iv][1];
        for l in ivp..nsys {
            input.vm[iv][l][0] -= vtmp1 * input.vm[iv][l][2];
            input.vm[iv][l][1] -= vtmp2 * input.vm[iv][l][2];
        }
        input.vdel[iv][0][0] -= vtmp1 * input.vdel[iv][2][0];
        input.vdel[iv][1][0] -= vtmp2 * input.vdel[iv][2][0];
        input.vdel[iv][0][1] -= vtmp1 * input.vdel[iv][2][1];
        input.vdel[iv][1][1] -= vtmp2 * input.vdel[iv][2][1];

        // Eliminate upper second column in VA block (row 1)
        let vtmp = input.va[iv][0][1];
        for l in ivp..nsys {
            input.vm[iv][l][0] -= vtmp * input.vm[iv][l][1];
        }
        input.vdel[iv][0][0] -= vtmp * input.vdel[iv][1][0];
        input.vdel[iv][0][1] -= vtmp * input.vdel[iv][1][1];

        if iv == nsys - 1 {
            continue;
        }

        // ====== Eliminate VB(IV+1) block, rows 1 -> 3 ======
        for k in 0..3 {
            let vtmp1 = input.vb[ivp][k][0];
            let vtmp2 = input.vb[ivp][k][1];
            let vtmp3 = input.vm[ivp][iv][k];
            for l in ivp..nsys {
                input.vm[ivp][l][k] -= vtmp1 * input.vm[iv][l][0]
                    + vtmp2 * input.vm[iv][l][1]
                    + vtmp3 * input.vm[iv][l][2];
            }
            input.vdel[ivp][k][0] -= vtmp1 * input.vdel[iv][0][0]
                + vtmp2 * input.vdel[iv][1][0]
                + vtmp3 * input.vdel[iv][2][0];
            input.vdel[ivp][k][1] -= vtmp1 * input.vdel[iv][0][1]
                + vtmp2 * input.vdel[iv][1][1]
                + vtmp3 * input.vdel[iv][2][1];
        }

        // Handle VZ block at trailing edge (coupling from upper to lower surface)
        if let (Some(ivte1), Some(ivz)) = (input.ivte1, input.ivz) {
            if iv == ivte1 {
                for k in 0..3 {
                    let vtmp1 = input.vz[k][0];
                    let vtmp2 = input.vz[k][1];
                    for l in ivp..nsys {
                        input.vm[ivz][l][k] -=
                            vtmp1 * input.vm[iv][l][0] + vtmp2 * input.vm[iv][l][1];
                    }
                    input.vdel[ivz][k][0] -=
                        vtmp1 * input.vdel[iv][0][0] + vtmp2 * input.vdel[iv][1][0];
                    input.vdel[ivz][k][1] -=
                        vtmp1 * input.vdel[iv][0][1] + vtmp2 * input.vdel[iv][1][1];
                }
            }
        }

        if ivp == nsys - 1 {
            continue;
        }

        // ====== Eliminate lower VM column (sparse elimination) ======
        for kv in (iv + 2)..nsys {
            let vtmp1 = input.vm[kv][iv][0];
            let vtmp2 = input.vm[kv][iv][1];
            let vtmp3 = input.vm[kv][iv][2];

            if vtmp1.abs() > vacc1 {
                for l in ivp..nsys {
                    input.vm[kv][l][0] -= vtmp1 * input.vm[iv][l][2];
                }
                input.vdel[kv][0][0] -= vtmp1 * input.vdel[iv][2][0];
                input.vdel[kv][0][1] -= vtmp1 * input.vdel[iv][2][1];
            }

            if vtmp2.abs() > vacc2 {
                for l in ivp..nsys {
                    input.vm[kv][l][1] -= vtmp2 * input.vm[iv][l][2];
                }
                input.vdel[kv][1][0] -= vtmp2 * input.vdel[iv][2][0];
                input.vdel[kv][1][1] -= vtmp2 * input.vdel[iv][2][1];
            }

            if vtmp3.abs() > vacc3 {
                for l in ivp..nsys {
                    input.vm[kv][l][2] -= vtmp3 * input.vm[iv][l][2];
                }
                input.vdel[kv][2][0] -= vtmp3 * input.vdel[iv][2][0];
                input.vdel[kv][2][1] -= vtmp3 * input.vdel[iv][2][1];
            }
        }
    }

    // Backward sweep: IV = NSYS-1 down to 1
    for iv in (1..nsys).rev() {
        // Eliminate upper VM columns
        let vtmp = input.vdel[iv][2][0];
        for kv in (0..iv).rev() {
            input.vdel[kv][0][0] -= input.vm[kv][iv][0] * vtmp;
            input.vdel[kv][1][0] -= input.vm[kv][iv][1] * vtmp;
            input.vdel[kv][2][0] -= input.vm[kv][iv][2] * vtmp;
        }

        let vtmp = input.vdel[iv][2][1];
        for kv in (0..iv).rev() {
            input.vdel[kv][0][1] -= input.vm[kv][iv][0] * vtmp;
            input.vdel[kv][1][1] -= input.vm[kv][iv][1] * vtmp;
            input.vdel[kv][2][1] -= input.vm[kv][iv][2] * vtmp;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blsolv_decoupled_diagonal() {
        // Test with a decoupled diagonal system where each station is independent
        // This tests the basic elimination without inter-station coupling
        let nsys = 3;
        let mut input = BlsolvInput {
            nsys,
            va: vec![[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]; nsys],
            vb: vec![[[0.0; 2]; 3]; nsys],
            vdel: vec![[[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]; nsys],
            vm: vec![vec![[0.0, 0.0, 0.0]; nsys]; nsys],
            vz: [[0.0; 2]; 3],
            ivte1: None,
            ivz: None,
            vaccel: 0.01,
            arc_length: None,
        };

        // Set diagonal VM entries for the third equation
        for iv in 0..nsys {
            input.vm[iv][iv][2] = 1.0;
        }

        blsolv(&mut input);

        // Solution should be [1, 2, 3] at each station
        for iv in 0..nsys {
            assert!(
                (input.vdel[iv][0][0] - 1.0).abs() < 1e-10,
                "Station {} row 0: expected 1.0, got {}",
                iv,
                input.vdel[iv][0][0]
            );
            assert!(
                (input.vdel[iv][1][0] - 2.0).abs() < 1e-10,
                "Station {} row 1: expected 2.0, got {}",
                iv,
                input.vdel[iv][1][0]
            );
            assert!(
                (input.vdel[iv][2][0] - 3.0).abs() < 1e-10,
                "Station {} row 2: expected 3.0, got {}",
                iv,
                input.vdel[iv][2][0]
            );
        }
    }

    #[test]
    fn test_blsolv_with_coupling() {
        // Test with a system that has off-diagonal VM coupling
        let nsys = 2;
        let mut input = BlsolvInput {
            nsys,
            va: vec![[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]; nsys],
            vb: vec![[[0.0; 2]; 3]; nsys],
            vdel: vec![[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]; nsys],
            vm: vec![vec![[0.0, 0.0, 0.0]; nsys]; nsys],
            vz: [[0.0; 2]; 3],
            ivte1: None,
            ivz: None,
            vaccel: 0.01,
            arc_length: None,
        };

        // Set diagonal VM entries (required for the third equation)
        input.vm[0][0][2] = 1.0;
        input.vm[1][1][2] = 1.0;
        // Add some coupling from station 0 to station 1
        input.vm[1][0] = [0.1, 0.1, 0.1];

        blsolv(&mut input);

        // Solution should satisfy the system equations
        // This is a basic sanity check - actual values depend on the coupling
        assert!(input.vdel[0][0][0].is_finite(), "Station 0 row 0 is NaN");
        assert!(input.vdel[0][1][0].is_finite(), "Station 0 row 1 is NaN");
        assert!(input.vdel[0][2][0].is_finite(), "Station 0 row 2 is NaN");
        assert!(input.vdel[1][0][0].is_finite(), "Station 1 row 0 is NaN");
        assert!(input.vdel[1][1][0].is_finite(), "Station 1 row 1 is NaN");
        assert!(input.vdel[1][2][0].is_finite(), "Station 1 row 2 is NaN");
    }
}
