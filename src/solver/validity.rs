//! The validity record: the evidence for whether a reported result is valid.
//!
//! Every field here is **write-only inside the solver**. Nothing in `src/solver/` or `src/bl/`
//! may read one or branch on one: the record exists so that the output layer can classify a
//! result without the numerics ever depending on the classification (CLAUDE.md Rule 2).
//! `ci/no-deviations.sh` enforces that mechanically, so a read added here fails CI rather than
//! silently turning an observation into a control decision.
//!
//! Nothing is inferred from the size or the smoothness of a returned number. On the worked case
//! (`docs/xfoil-known-issues.md` §7.7) α = −1° reports CL = −52.09 and α = −5° reports +0.128,
//! and the two are equally broken — the tame-looking one is tame precisely because it is *more*
//! deeply out of domain.
//!
//! Three kinds of thing can be wrong with a solved point, and they are classified by *what it
//! would take to know* rather than by severity (`docs/guide/validity.md`, "How failures are
//! classified"). Only the first two are recorded here, both exact from a single run and both
//! derived from tests XFOIL itself performs and then discards:
//!
//! - **Class A, domain violation.** A formula was evaluated outside its analytic domain. The
//!   Kármán–Tsien denominator `β + BFAC·Cp_inc` has a pole at `q/Q∞ = √(1 + 2β(1+β)/M∞²)`;
//!   past it Cp has changed sign rather than merely lost accuracy. CPCALC tests exactly this
//!   into its `DENNEG` flag, prints `Local speed too large`, and returns the inverted Cp anyway;
//!   CLCALC forms the same denominator inline and does not warn at all.
//! - **Class B, iteration exhaustion.** A Newton ran to its cap without meeting its tolerance, so
//!   the value returned is not a solution to the problem that was posed.
//!
//! **Class C, conditioning** — a converged, in-domain point that a 1-ULP perturbation moves by
//! O(1) — cannot be seen from one run at all, because it is a property of the problem near that
//! point rather than of the arithmetic performed. It is not recorded here and does not make a
//! point invalid; the twins every validation run carries measure it. So a `Valid` verdict from this record is a
//! necessary condition for trusting a number, not a sufficient one.
//!
//! Each field describes **the state currently stored**, not the history of the point: a writer
//! overwrites, so after a point is solved the record refers to the values that point reports.

/// Observational record attached to `SolverState`, written by the routines whose formulas and
/// iterations it describes and never read back by them. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValidityRecord {
    /// Smallest Kármán–Tsien denominator `β + BFAC·Cp_inc` over the nodes of the CLCALC
    /// integration that produced the stored CL, CM and CDP. At or below zero the pole has been
    /// crossed and the forces are out of domain. `INFINITY` before any evaluation, and always
    /// `1` at M = 0, where `β = 1` and `BFAC = 0`.
    pub karman_tsien_margin_forces: f64,
    /// The same denominator, smallest over the nodes of the stored Cp arrays — CPCALC's own
    /// `DENNEG` test, kept as the margin rather than collapsed to its sign.
    pub karman_tsien_margin_pressure: f64,
    /// SPECAL's CL(M) Newton used all 20 iterations without reaching `|DCLM| ≤ 1e-6`
    /// (`SPECAL:  Minf convergence failed`). Only reachable for MATYP ≠ 1, since a fixed Mach
    /// converges on the first step.
    pub mach_cl_newton_exhausted: bool,
    /// MRCL floored the lift coefficient (`CLA = MAX(CLS, 1e-6)`) when it set the stored Mach
    /// and Reynolds number: the CL it was handed was at or below zero, so the CL-dependent type
    /// had no real value to work from and a substitute was used instead.
    pub cl_floored: bool,
    /// MRCL limited the stored Mach to 0.99 (`CL too low for chosen Mach(CL) dependence`),
    /// which also zeroes `M_CLS` and so flattens the Newton that called it.
    pub mach_limited: bool,
    /// MRCL limited the stored Reynolds number to 100 × REINF1
    /// (`CL too low for chosen Re(CL) dependence`).
    pub re_limited: bool,
}

impl Default for ValidityRecord {
    fn default() -> Self {
        Self {
            karman_tsien_margin_forces: f64::INFINITY,
            karman_tsien_margin_pressure: f64::INFINITY,
            mach_cl_newton_exhausted: false,
            cl_floored: false,
            mach_limited: false,
            re_limited: false,
        }
    }
}

impl ValidityRecord {
    /// A formula that produced a stored value was evaluated outside its domain, so that value is
    /// wrong rather than imprecise. Read by the output layer only.
    pub fn out_of_domain(&self) -> bool {
        self.karman_tsien_margin_forces <= 0.0 || self.karman_tsien_margin_pressure <= 0.0
    }

    /// An iteration that produced a stored value ran out without meeting its tolerance, so that
    /// value is not a solution to the problem posed. Read by the output layer only.
    ///
    /// Covers only what this record carries; VISCAL's own convergence is `PointResult::converged`
    /// and SPECCL's is `PointResult::inviscid_cl_iterations`.
    pub fn iteration_exhausted(&self) -> bool {
        self.mach_cl_newton_exhausted
    }

    /// The flow conditions the point was solved at are not the ones that were asked for: MRCL
    /// substituted a floored CL or a limited Mach/Re because the request had no real solution.
    pub fn conditions_substituted(&self) -> bool {
        self.cl_floored || self.mach_limited || self.re_limited
    }
}
