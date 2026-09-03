//! The only tolerance constants used by XFOIL-equivalence tests (CLAUDE.md Rule 1).
//!
//! Derived from the measured noise floor (`scripts/noise-floor.sh`, recorded in
//! `docs/validation/noise-floor.md`, 2026-09-03): perturbing every panel coordinate by +1 ULP on the
//! tracked reference case, with identical branch trace and iteration count, moves
//!   - geometry-derived quantities by ≤ 5e-14 relative,
//!   - inviscid GAM/QINV by ≤ 6e-11 (median 5e-14),
//!   - BL state after MRCHUE/MRCHDU by ≤ 2e-11,
//!   - per-iteration CL/CD/RMSBL by ≤ 2.2e-10 (median 1e-12).
//! Element-wise relative spread of DIJ (≤ 5e-6) and SETBL Jacobian entries (≤ 9e-4) is dominated by
//! near-zero entries: those matrices must be gated with a row-scaled `scale`, never bare relative.

/// Pure closure functions and other single-expression translations.
pub const TOL_PURE: f64 = 1e-12;
/// Linear solves (BLSOLV, LUDCMP/BAKSUB, DIJ assembly).
pub const TOL_LINALG: f64 = 1e-11;
/// State after one Newton iteration of a nonlinear solver, given identical branch trace.
pub const TOL_SOLVER: f64 = 1e-10;

/// The one error metric: `|a − b| ≤ tol · max(|a|, |b|, scale)`. Bare relative error is
/// undefined at CL≈0, VDEL≈0, laminar CTAU≈0; `scale` is the physical scale of the variable.
pub fn within(a: f64, b: f64, tol: f64, scale: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(scale)
}

/// Assert `within`, reporting the values, the achieved ratio and the tolerance on failure.
pub fn assert_within(a: f64, b: f64, tol: f64, scale: f64, what: &str) {
    let denom = a.abs().max(b.abs()).max(scale);
    let err = (a - b).abs() / denom;
    assert!(
        err <= tol,
        "{what}: yfoil={a:.17e} xfoil={b:.17e} err={err:.3e} > tol={tol:.1e} (scale {scale:.1e})"
    );
}
