# Reference noise floor

**What this measures.** How much XFOIL's own outputs move when its *input* changes by the smallest
possible amount — every panel coordinate perturbed by +1 ULP (`math.nextafter`) — with everything else
identical. This is the floor below which no comparison between YFoil and XFOIL can be meaningful, and the
tolerances in `tests/utilities/tolerances.rs` are derived from it (CLAUDE.md Rule 1).

**How.** `scripts/noise-floor.sh` takes the tracked reference case
(`tests/fixtures/xfoil/naca0012_n60_a2_re1e6`: NACA 0012, N=60, α=2°, Re=1e6, M=0, Ncrit=9, ITER 20),
runs the instrumented double-precision reference twice (base and perturbed), and compares every
`ES24.16` dump line-by-line. Metric per value: `|a−b|` and, for `max(|a|,|b|) > 1e-8`, `|a−b| / max(|a|,|b|)`.
Both runs converged in the same number of iterations (5) with the same branch trace, so this is the
*continuous* floor — branch flips are a separate phenomenon (plan addendum R3).

**Result (2026-09-03, arm64 macOS, gfortran 14.2, `-O -fdefault-real-8 -ffp-contract=off`).**

| file | lines | misaligned | max abs diff | max rel (abs>1e-8) | p99 rel | median rel | stage |
|---|---|---|---|---|---|---|---|
| xfoil_panels.dat | 71 | 0 | 7.30e-15 | 5.25e-14 | 2.45e-14 | 1.45e-16 | geometry after ABCOPY |
| xfoil_inviscid.dat | 129 | 0 | 2.73e-12 | 5.78e-11 | 7.02e-12 | 4.99e-14 | inviscid GAM/QINV |
| xfoil_dij.dat | 5628 | 0 | 9.78e-09 | 5.20e-06 | 2.40e-07 | 1.67e-11 | DIJ matrix |
| mrchdu_input.dat | 85 | 0 | 1.59e-11 | 2.16e-11 | 8.22e-12 | 1.72e-13 | BL state after MRCHUE init |
| mrchdu_output.dat | 79 | 0 | 1.59e-11 | 2.16e-11 | 8.22e-12 | 1.73e-13 | BL state after MRCHDU #1 |
| blsolv_input.dat | 17802 | 0 | 1.62e-04 | 8.85e-04 | 1.33e-06 | 4.87e-11 | SETBL VA/VB/VM/VDEL (calls 1..3) |
| blsolv_output.dat | 231 | 0 | 1.23e-10 | 3.18e-08 | 1.74e-08 | 5.10e-11 | BLSOLV VDEL out |
| xfoil_update.dat | 770 | 0 | 2.00e-12 | 3.18e-08 | 9.02e-10 | 3.77e-13 | UPDATE DUI/UNEW |
| xfoil_iter_simple.dat | 90 | 0 | 5.98e-11 | 2.16e-10 | 4.83e-11 | 9.10e-13 | per-iteration CL/CD/RMSBL |

base: iterations=5 converged=True

ulp: iterations=5 converged=True

**Interpretation.**

- Geometry-derived quantities (`S`, normals, `APANEL`) sit at ≤ 5e-14: a 1-ULP input amplified ~200×
  through spline and arc-length evaluation. This is the floor for anything gated at "pure function"
  precision.
- The inviscid solve moves by ≤ 6e-11 (median 5e-14) and the BL state after MRCHUE/MRCHDU by ≤ 2e-11.
- Per-iteration **CL, CD and RMSBL move by up to 2.2e-10** (median 1e-12). So a "1e-10 on forces"
  acceptance is *at* the floor, not comfortably above it; the S9 gate of 1e-9 carries roughly a 5×
  safety factor and should not be tightened without re-measuring on the case in question.
- **Element-wise relative spread is not a usable metric for matrices.** DIJ entries move by up to
  5e-6 and SETBL Jacobian entries by up to 9e-4 *relative*, but the medians are 2e-11 and 5e-11: the
  tails are near-zero entries where relative error is meaningless. DIJ and VA/VB/VM/VDEL must be gated
  with the scaled metric `|a−b| ≤ tol · max(|a|, |b|, scale_row)`, with `scale_row` the row's largest
  magnitude.

**What it does not cover.** One case, one perturbation direction, M=0, blunt TE. Re-run on every case
added to `fixtures/cases.toml` (the branch-coverage cases especially) before trusting a tolerance there,
and on any other host before comparing fixtures across hosts (libm differences, addendum R2).
