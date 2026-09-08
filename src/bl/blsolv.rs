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

/// Mirrors XFOIL's `VA`, `VB`, `VM`, `VZ`, `VDEL`.
///
/// Input matrices for BLSOLV
#[derive(Debug, Clone)]
#[doc(alias = "VA")]
#[doc(alias = "VB")]
#[doc(alias = "VM")]
#[doc(alias = "VZ")]
#[doc(alias = "VDEL")]
pub struct NewtonSystem {
    /// Number of BL stations
    pub n_rows: usize,
    /// Diagonal blocks `VA[iv][k][l]` - 3 equations, 2 columns
    pub diagonal: Vec<[[f64; 2]; 3]>,
    /// Sub-diagonal blocks `VB[iv][k][l]`
    pub subdiagonal: Vec<[[f64; 2]; 3]>,
    /// RHS/solution `VDEL[iv][k][l]` - column 0 is residual, column 1 is Re sensitivity
    pub rhs: Vec<[[f64; 2]; 3]>,
    /// Mass defect coupling `VM[iv][j][k]` - coupling from station j to station iv, equation k
    pub mass_influence: Vec<Vec<[f64; 3]>>,
    /// TE coupling block `VZ[k][l]` (optional, only used at trailing edge)
    pub te_block: [[f64; 2]; 3],
    /// Index where upper surface ends at TE (0-based)
    pub i_te_row_upper: Option<usize>,
    /// System index for start of wake (lower surface TE + 1)
    pub i_wake_row: Option<usize>,
    /// Acceleration parameter for sparse elimination
    pub elimination_threshold: f64,
    /// Total arc length S(N) - S(1) for VACC scaling (optional)
    /// If provided, VACC2 and VACC3 are scaled by 2.0 / arc_length
    pub s_total: Option<f64>,
}

impl Default for NewtonSystem {
    fn default() -> Self {
        Self {
            n_rows: 0,
            diagonal: Vec::new(),
            subdiagonal: Vec::new(),
            rhs: Vec::new(),
            mass_influence: Vec::new(),
            te_block: [[0.0; 2]; 3],
            i_te_row_upper: None,
            i_wake_row: None,
            elimination_threshold: 0.01,
            s_total: None,
        }
    }
}

/// Optional trace of BLSOLV's branch decisions and intermediate state, mirroring the
/// instrumentation in the reference build (xsolve.f) so the two can be diffed station by
/// station (CLAUDE.md Rule 3). Zero cost when not requested.
#[derive(Debug, Default, Clone)]
pub struct BlsolvTrace {
    /// Every sparse-elimination comparison `|VTMP| > VACC`: (iv, kv, k, |vtmp|, vacc, taken).
    pub skips: Vec<(usize, usize, usize, f64, f64, bool)>,
    /// Full VDEL after the forward sweep, before back-substitution.
    pub rhs_after_forward: Vec<[[f64; 2]; 3]>,
}

impl BlsolvTrace {
    /// The comparison closest to its threshold, as a relative margin `(|vtmp| - vacc) / vacc`.
    pub fn tightest_skip_margin(&self) -> Option<(usize, usize, usize, f64)> {
        self.skips
            .iter()
            .map(|&(i_row, kv, k, v, vacc, _)| (i_row, kv, k, (v - vacc) / vacc))
            .min_by(|a, b| a.3.abs().partial_cmp(&b.3.abs()).unwrap())
    }
}

/// The Newton deltas produced by [`solve_newton_system`].
///
/// This is a separate type on purpose: XFOIL's UPDATE aliases `UNEW` onto `VA` and `QNEW`
/// onto `VB` via EQUIVALENCE (xbl.f), so after BLSOLV the factored VA/VB/VM blocks are dead.
/// Consuming the [`NewtonSystem`] makes it impossible to read them by accident.
#[derive(Debug, Clone)]
pub struct NewtonDeltas {
    /// Number of system rows
    pub n_rows: usize,
    /// `VDEL[iv][k][l]`: column 0 is the Newton delta, column 1 the Re/alpha sensitivity
    pub deltas: Vec<[[f64; 2]; 3]>,
}

/// Solve the coupled BL Newton system using XFOIL's BLSOLV algorithm
///
/// This is a direct translation of XFOIL's BLSOLV subroutine (xsolve.f). The algorithm
/// performs block Gaussian elimination with special handling for the dense mass defect
/// coupling (VM matrix). Verified bit-identical to XFOIL on the tracked reference fixture
/// (all three calls, both columns) — see tests/xfoil_blsolv_tests.rs.
#[doc(alias = "BLSOLV")]
pub fn solve_newton_system(input: NewtonSystem) -> NewtonDeltas {
    solve_newton_system_traced(input, None)
}

/// `blsolv` with an optional [`BlsolvTrace`] collector.
pub fn solve_newton_system_traced(input: NewtonSystem, mut trace: Option<&mut BlsolvTrace>) -> NewtonDeltas {
    let mut input = input;
    let nsys = input.n_rows;
    if nsys == 0 {
        return NewtonDeltas {
            n_rows: nsys,
            deltas: input.rhs,
        };
    }

    // Compute acceleration thresholds
    // XFOIL: VACC1 = VACCEL, VACC2 = VACC3 = VACCEL * 2.0 / (S(N) - S(1))
    // Association must match the Fortran exactly — (VACCEL*2.0)/(S(N)-S(1)) — because these
    // thresholds gate branches; a 1-ULP difference in VACC2 can flip a skip decision.
    let vacc1 = input.elimination_threshold;
    let (vacc2, vacc3) = match input.s_total {
        Some(arc_len) if arc_len > 0.0 => (
            input.elimination_threshold * 2.0 / arc_len,
            input.elimination_threshold * 2.0 / arc_len,
        ),
        _ => (input.elimination_threshold, input.elimination_threshold),
    };

    // Forward sweep: IV = 0 to NSYS-1
    for i_row in 0..nsys {
        let ivp = i_row + 1;

        // ====== Invert VA(IV) block ======

        // Normalize first row by VA(1,1)
        let pivot = 1.0 / input.diagonal[i_row][0][0];
        input.diagonal[i_row][0][1] *= pivot;
        for l in i_row..nsys {
            input.mass_influence[i_row][l][0] *= pivot;
        }
        input.rhs[i_row][0][0] *= pivot;
        input.rhs[i_row][0][1] *= pivot;

        // Eliminate lower first column in VA block (rows 2,3)
        for k in 1..3 {
            let vtmp = input.diagonal[i_row][k][0];
            input.diagonal[i_row][k][1] -= vtmp * input.diagonal[i_row][0][1];
            for l in i_row..nsys {
                input.mass_influence[i_row][l][k] -= vtmp * input.mass_influence[i_row][l][0];
            }
            input.rhs[i_row][k][0] -= vtmp * input.rhs[i_row][0][0];
            input.rhs[i_row][k][1] -= vtmp * input.rhs[i_row][0][1];
        }

        // Normalize second row by VA(2,2)
        let pivot = 1.0 / input.diagonal[i_row][1][1];
        for l in i_row..nsys {
            input.mass_influence[i_row][l][1] *= pivot;
        }
        input.rhs[i_row][1][0] *= pivot;
        input.rhs[i_row][1][1] *= pivot;

        // Eliminate lower second column in VA block (row 3)
        let vtmp = input.diagonal[i_row][2][1];
        for l in i_row..nsys {
            input.mass_influence[i_row][l][2] -= vtmp * input.mass_influence[i_row][l][1];
        }
        input.rhs[i_row][2][0] -= vtmp * input.rhs[i_row][1][0];
        input.rhs[i_row][2][1] -= vtmp * input.rhs[i_row][1][1];

        // Normalize third row by VM(3,IV,IV) - the diagonal mass coupling
        let pivot = 1.0 / input.mass_influence[i_row][i_row][2];
        for l in ivp..nsys {
            input.mass_influence[i_row][l][2] *= pivot;
        }
        input.rhs[i_row][2][0] *= pivot;
        input.rhs[i_row][2][1] *= pivot;

        // Eliminate upper third column in VA block (rows 1,2)
        let vtmp1 = input.mass_influence[i_row][i_row][0];
        let vtmp2 = input.mass_influence[i_row][i_row][1];
        for l in ivp..nsys {
            input.mass_influence[i_row][l][0] -= vtmp1 * input.mass_influence[i_row][l][2];
            input.mass_influence[i_row][l][1] -= vtmp2 * input.mass_influence[i_row][l][2];
        }
        input.rhs[i_row][0][0] -= vtmp1 * input.rhs[i_row][2][0];
        input.rhs[i_row][1][0] -= vtmp2 * input.rhs[i_row][2][0];
        input.rhs[i_row][0][1] -= vtmp1 * input.rhs[i_row][2][1];
        input.rhs[i_row][1][1] -= vtmp2 * input.rhs[i_row][2][1];

        // Eliminate upper second column in VA block (row 1)
        let vtmp = input.diagonal[i_row][0][1];
        for l in ivp..nsys {
            input.mass_influence[i_row][l][0] -= vtmp * input.mass_influence[i_row][l][1];
        }
        input.rhs[i_row][0][0] -= vtmp * input.rhs[i_row][1][0];
        input.rhs[i_row][0][1] -= vtmp * input.rhs[i_row][1][1];

        if i_row == nsys - 1 {
            continue;
        }

        // ====== Eliminate VB(IV+1) block, rows 1 -> 3 ======
        for k in 0..3 {
            let vtmp1 = input.subdiagonal[ivp][k][0];
            let vtmp2 = input.subdiagonal[ivp][k][1];
            let vtmp3 = input.mass_influence[ivp][i_row][k];
            for l in ivp..nsys {
                input.mass_influence[ivp][l][k] -= vtmp1 * input.mass_influence[i_row][l][0]
                    + vtmp2 * input.mass_influence[i_row][l][1]
                    + vtmp3 * input.mass_influence[i_row][l][2];
            }
            input.rhs[ivp][k][0] -=
                vtmp1 * input.rhs[i_row][0][0] + vtmp2 * input.rhs[i_row][1][0] + vtmp3 * input.rhs[i_row][2][0];
            input.rhs[ivp][k][1] -=
                vtmp1 * input.rhs[i_row][0][1] + vtmp2 * input.rhs[i_row][1][1] + vtmp3 * input.rhs[i_row][2][1];
        }

        // Handle VZ block at trailing edge (coupling from upper to lower surface)
        if let (Some(ivte1), Some(ivz)) = (input.i_te_row_upper, input.i_wake_row) {
            if i_row == ivte1 {
                for k in 0..3 {
                    let vtmp1 = input.te_block[k][0];
                    let vtmp2 = input.te_block[k][1];
                    for l in ivp..nsys {
                        input.mass_influence[ivz][l][k] -=
                            vtmp1 * input.mass_influence[i_row][l][0] + vtmp2 * input.mass_influence[i_row][l][1];
                    }
                    input.rhs[ivz][k][0] -= vtmp1 * input.rhs[i_row][0][0] + vtmp2 * input.rhs[i_row][1][0];
                    input.rhs[ivz][k][1] -= vtmp1 * input.rhs[i_row][0][1] + vtmp2 * input.rhs[i_row][1][1];
                }
            }
        }

        if ivp == nsys - 1 {
            continue;
        }

        // ====== Eliminate lower VM column (sparse elimination) ======
        for kv in (i_row + 2)..nsys {
            let vtmp1 = input.mass_influence[kv][i_row][0];
            let vtmp2 = input.mass_influence[kv][i_row][1];
            let vtmp3 = input.mass_influence[kv][i_row][2];
            if let Some(t) = trace.as_mut() {
                t.skips.push((i_row, kv, 0, vtmp1.abs(), vacc1, vtmp1.abs() > vacc1));
                t.skips.push((i_row, kv, 1, vtmp2.abs(), vacc2, vtmp2.abs() > vacc2));
                t.skips.push((i_row, kv, 2, vtmp3.abs(), vacc3, vtmp3.abs() > vacc3));
            }

            if vtmp1.abs() > vacc1 {
                for l in ivp..nsys {
                    input.mass_influence[kv][l][0] -= vtmp1 * input.mass_influence[i_row][l][2];
                }
                input.rhs[kv][0][0] -= vtmp1 * input.rhs[i_row][2][0];
                input.rhs[kv][0][1] -= vtmp1 * input.rhs[i_row][2][1];
            }

            if vtmp2.abs() > vacc2 {
                for l in ivp..nsys {
                    input.mass_influence[kv][l][1] -= vtmp2 * input.mass_influence[i_row][l][2];
                }
                input.rhs[kv][1][0] -= vtmp2 * input.rhs[i_row][2][0];
                input.rhs[kv][1][1] -= vtmp2 * input.rhs[i_row][2][1];
            }

            if vtmp3.abs() > vacc3 {
                for l in ivp..nsys {
                    input.mass_influence[kv][l][2] -= vtmp3 * input.mass_influence[i_row][l][2];
                }
                input.rhs[kv][2][0] -= vtmp3 * input.rhs[i_row][2][0];
                input.rhs[kv][2][1] -= vtmp3 * input.rhs[i_row][2][1];
            }
        }
    }

    if let Some(t) = trace.as_mut() {
        t.rhs_after_forward = input.rhs.clone();
    }

    // Backward sweep: IV = NSYS-1 down to 1
    for i_row in (1..nsys).rev() {
        // Eliminate upper VM columns
        let vtmp = input.rhs[i_row][2][0];
        for kv in (0..i_row).rev() {
            input.rhs[kv][0][0] -= input.mass_influence[kv][i_row][0] * vtmp;
            input.rhs[kv][1][0] -= input.mass_influence[kv][i_row][1] * vtmp;
            input.rhs[kv][2][0] -= input.mass_influence[kv][i_row][2] * vtmp;
        }

        let vtmp = input.rhs[i_row][2][1];
        for kv in (0..i_row).rev() {
            input.rhs[kv][0][1] -= input.mass_influence[kv][i_row][0] * vtmp;
            input.rhs[kv][1][1] -= input.mass_influence[kv][i_row][1] * vtmp;
            input.rhs[kv][2][1] -= input.mass_influence[kv][i_row][2] * vtmp;
        }
    }

    NewtonDeltas {
        n_rows: nsys,
        deltas: input.rhs,
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
        let mut input = NewtonSystem {
            n_rows: nsys,
            diagonal: vec![[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]; nsys],
            subdiagonal: vec![[[0.0; 2]; 3]; nsys],
            rhs: vec![[[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]]; nsys],
            mass_influence: vec![vec![[0.0, 0.0, 0.0]; nsys]; nsys],
            te_block: [[0.0; 2]; 3],
            i_te_row_upper: None,
            i_wake_row: None,
            elimination_threshold: 0.01,
            s_total: None,
        };

        // Set diagonal VM entries for the third equation
        for i_row in 0..nsys {
            input.mass_influence[i_row][i_row][2] = 1.0;
        }

        let input = solve_newton_system(input);

        // Solution should be [1, 2, 3] at each station
        for i_row in 0..nsys {
            assert!(
                (input.deltas[i_row][0][0] - 1.0).abs() < 1e-10,
                "Station {} row 0: expected 1.0, got {}",
                i_row,
                input.deltas[i_row][0][0]
            );
            assert!(
                (input.deltas[i_row][1][0] - 2.0).abs() < 1e-10,
                "Station {} row 1: expected 2.0, got {}",
                i_row,
                input.deltas[i_row][1][0]
            );
            assert!(
                (input.deltas[i_row][2][0] - 3.0).abs() < 1e-10,
                "Station {} row 2: expected 3.0, got {}",
                i_row,
                input.deltas[i_row][2][0]
            );
        }
    }

    #[test]
    fn test_blsolv_with_coupling() {
        // Test with a system that has off-diagonal VM coupling
        let nsys = 2;
        let mut input = NewtonSystem {
            n_rows: nsys,
            diagonal: vec![[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]; nsys],
            subdiagonal: vec![[[0.0; 2]; 3]; nsys],
            rhs: vec![[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]]; nsys],
            mass_influence: vec![vec![[0.0, 0.0, 0.0]; nsys]; nsys],
            te_block: [[0.0; 2]; 3],
            i_te_row_upper: None,
            i_wake_row: None,
            elimination_threshold: 0.01,
            s_total: None,
        };

        // Set diagonal VM entries (required for the third equation)
        input.mass_influence[0][0][2] = 1.0;
        input.mass_influence[1][1][2] = 1.0;
        // Add some coupling from station 0 to station 1
        input.mass_influence[1][0] = [0.1, 0.1, 0.1];

        let input = solve_newton_system(input);

        // Solution should satisfy the system equations
        // This is a basic sanity check - actual values depend on the coupling
        assert!(input.deltas[0][0][0].is_finite(), "Station 0 row 0 is NaN");
        assert!(input.deltas[0][1][0].is_finite(), "Station 0 row 1 is NaN");
        assert!(input.deltas[0][2][0].is_finite(), "Station 0 row 2 is NaN");
        assert!(input.deltas[1][0][0].is_finite(), "Station 1 row 0 is NaN");
        assert!(input.deltas[1][1][0].is_finite(), "Station 1 row 1 is NaN");
        assert!(input.deltas[1][2][0].is_finite(), "Station 1 row 2 is NaN");
    }
}
