//! The only tolerance constants used by XFOIL-equivalence tests (CLAUDE.md Rule 1).
//!
//! Derived from the measured noise floor (`scripts/noise-floor.sh`, recorded in
//! `docs/validation/noise-floor.md`, 2026-09-03). Perturbing every panel coordinate by +1 ULP on
//! the tracked reference case, with identical branch trace and iteration count, moves
//! geometry-derived quantities by at most 5e-14 relative, inviscid GAM/QINV by 6e-11 (median
//! 5e-14), BL state after MRCHUE/MRCHDU by 2e-11, and per-iteration CL/CD/RMSBL by 2.2e-10
//! (median 1e-12). Element-wise relative spread of DIJ (up to 5e-6) and SETBL Jacobian entries
//! (up to 9e-4) is dominated by near-zero entries: those matrices must be gated with a
//! row-scaled `scale`, never bare relative error.

/// Pure closure functions and other single-expression translations.
pub const TOL_PURE: f64 = 1e-12;
/// Linear solves (BLSOLV, LUDCMP/BAKSUB, DIJ assembly).
pub const TOL_LINALG: f64 = 1e-11;
/// State after one Newton iteration of a nonlinear solver, given identical branch trace.
pub const TOL_SOLVER: f64 = 1e-10;
/// Per-iteration *transients* (RLX, RMSBL, unconverged CL/CD) inside a multi-point sequence.
/// Measured 2026-09-04 on the tracked polar script (11 VISCAL calls, 62 iterations): a 1-ULP
/// geometry perturbation moves the reference's own intermediate RLX/RMSBL/CL by up to 5.1e-10
/// absolute while the converged per-point CL/CD/CM move by ≤ 8e-14 and XTR by ≤ 1e-12. Converged
/// points therefore stay at `TOL_SOLVER`; transients are gated at floor × 2.
pub const TOL_TRANSIENT: f64 = 1e-9;
/// One SETBL/UPDATE step replayed from XFOIL's dumped state on a host other than the one the
/// fixture was generated on (`tests/utilities/host.rs`). Both codes take their transcendentals
/// from the host libm, and Apple libSystem and glibc differ by 1 ULP on 0.1 % (`exp`, `ln`,
/// `powf`) to 18 % (`tanh`) of inputs (measured 2026-09-11, 20 000 inputs per function).
/// Measured on the same date from identical dumped inputs, the yFoil step on glibc 2.39 differs
/// from the yFoil step on libSystem 1345.120.2 by 1.2e-10 and 3.6e-10 relative in RMXBL (the
/// largest Newton delta, at the polar break points' near-singular iterations) — while the same
/// step on the fixture's own host reproduces XFOIL to ≤ 2e-12. Cross-host gates are that
/// measured spread × ~3; same-host stays at `TOL_SOLVER`.
pub const TOL_CROSS_HOST: f64 = 1e-9;
/// Safety factor applied to a case's own measured 1-ULP floor (`noise_floor.json`, written by the
/// fixture pipeline from the +1-ULP twin run): a value is accepted when it is within the base
/// tolerance *or* within `FLOOR_FACTOR` × the reference's own spread of that value. Chaotic
/// trajectories (separated flow, sharp-TE transients) are gated by what the reference can
/// reproduce itself, never by an asserted number.
pub const FLOOR_FACTOR: f64 = 4.0;
/// Hypersensitivity threshold for the third outcome. When the reference's own +1-ULP twin moves a
/// per-iteration value by more than this (absolute), the reference cannot reproduce itself there;
/// a yFoil run that matched every earlier iteration within `FLOOR_FACTOR` × floor and departs at
/// such an iteration is classified **threshold-straddling** (reported, not passed or failed),
/// provided the one-step replay from XFOIL's exact state at that iteration matches.
pub const STRADDLE_FLOOR: f64 = 1e-6;

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

// ---------------------------------------------------------------------------------------------
// Geometry generators against the NASA/PDAS `naca456` reference (tests/fixtures/naca456/,
// `scripts/naca456-fixtures.sh`). These are not XFOIL equivalence gates: naca456 is an
// independent implementation of the NACA definitions, and its own precision sets the floor.
// ---------------------------------------------------------------------------------------------

/// Closed-form families (4-digit, 4-digit modified, 2-, 3- and 3-reflex mean lines): the same
/// polynomials evaluated in a different association; expected agreement is round-off.
pub const TOL_NACA456_CLOSED_FORM: f64 = 1e-12;
/// The 6-series and 6A mean lines: naca456 evaluates them with `PI = 3.141592654` (a 10-digit
/// literal, 1.1e-10 relative from π), which yFoil does not reproduce; the gate is 10 × that.
pub const TOL_NACA456_SIX_SERIES_MEAN_LINE: f64 = 1e-9;
/// The 6-series thickness forms are looked up in naca456 by inverting the arc-length spline
/// with Brent's method at `TOL = 1E-6` on the arc coordinate, and the ordinate is reported at
/// the station reached, not the one requested; its ordinate is therefore in error by up to
/// `NACA456_ROOT_TOL × |dy_t/dx|`. yFoil inverts to round-off, so the gate on `y_t` is that
/// slope-scaled bound plus the mapping's own noise (`3.14159265` for π in the φ grid: 1.1e-9
/// relative, through a unit-scale mapping and a spline).
pub const NACA456_ROOT_TOL: f64 = 1e-6;
pub const TOL_NACA456_SIX_SERIES_THICKNESS: f64 = 1e-8;
