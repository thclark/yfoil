# Tolerances

How close counts as equal. A translation of an iterative solver cannot be judged by "agrees to a
few digits": the two codes associate the same operations differently, so they differ in the last
bits, and those differences grow by exactly as much as the reference amplifies them itself. So
every tolerance is **derived from the reference's own behaviour**, not asserted, and a
comparison has two parts — the branch trace exact, and the values within a named tolerance by one
metric. If values differ by more than that with an identical branch trace, it is a bug: 1 %,
0.1 % and 1e-6 alike.

## The noise floor, in brief

The floor is how far the reference's own outputs move when its *input* moves by the least possible
amount — every panel coordinate jogged by one ULP — with everything else identical. No comparison
between the two codes can be sharper than that. It is measured per stage and per variable by
`scripts/noise-floor.sh`, and every measurement is recorded on
[the reference's noise floor](the-reference-noise-floor.md). Measured 2026-09-03 and remeasured on
the correctly rounded reference on 2026-10-15: forces and RMSBL move by up to 2.2e-10 per iteration
(1.55e-10 since), the boundary-layer state by 2e-11, the inviscid solution by 6e-11.

## Named tolerance constants

Four constants, in `tests/common/utilities/tolerances.rs`, are the only tolerances any
equivalence test uses; there are no ad-hoc literals anywhere else.

| Constant | Value | For | Derived from |
| --- | --- | --- | --- |
| `TOL_PURE` | 1e-12 | pure closure functions and other single-expression translations | geometry-derived quantities move by ≤ 5e-14 |
| `TOL_LINALG` | 1e-11 | linear solves (BLSOLV, LUDCMP/BAKSUB, DIJ assembly) | the BL state after a solve moves by ≤ 2e-11 |
| `TOL_SOLVER` | 1e-10 | the state after one Newton iteration, and a converged point, given an identical branch trace | forces and RMSBL move by up to 2.2e-10 per iteration: this tolerance is *at* the floor |
| `TOL_TRANSIENT` | 1e-9 | per-iteration transients (RLX, RMSBL, unconverged CL/CD) inside a multi-point sequence | measured 2026-09-04 on the tracked polar script: transients move by up to 5.1e-10 while converged points move by ≤ 8e-14 |

Expect pure closure functions at about 1e-14, linear solves at 1e-11 and converged Newton state at
1e-10. A single replayed step may also be held to four times its own measured sensitivity, when
that is larger ([step-by-step replay](../execution-equivalence/step-by-step-replay.md)). The same
file holds the naca456 constants, whose basis is naca456's own precision
([how the naca456 tables were produced](../../engineering-validation/generated-sections-against-naca456.md#how-the-naca456-tables-were-produced)).

## The error metric

```
|a − b| ≤ tol · max(|a|, |b|, scale_v)
```

with `scale_v` a physical scale per variable. **Bare relative error is never used**: it is
undefined at CL ≈ 0, VDEL ≈ 0 and laminar CTAU ≈ 0, all of which occur in ordinary cases. Matrices
(DIJ, the SETBL Jacobian) are compared row-scaled, `scale` the row's largest magnitude: their
element-wise relative spread is dominated by near-zero entries, up to 5e-6 for DIJ and 9e-4 for
the Jacobian against medians of 2e-11 and 5e-11.

## Divergent comparisons and ill-conditioned solutions

XFOIL contains about 130 exact real equalities and 140 real inequality branches: RLX clamps,
BLSOLV's sparsity skips, `AMPL > ACRIT`, `SHARP`, `RMSBL < 1e-4`, MRCHUE's `DIRECT` switch, the
Hk clamps. Where the operands of one are within rounding of each other, a last-bit difference
between the codes falls on opposite sides, the two take different routes, and the results differ
by O(1). That is a **divergent comparison** ([terminology](../../contributing/conventions/terminology.md#divergent-comparison)):
it says nothing about the translation, so no test is one. So:

- the instrumented reference logs its decisions — `ITRAN`, the limiting station, skip counts, the
  iteration count, `SHARP`, `DIRECT` — and those are compared **exactly**;
- for an iterative solver the iteration count must be identical, and RMSBL is compared per
  iteration.

A divergent comparison is a consequence of an **ill-conditioned solution** — stall, separation,
the Kármán–Tsien pole — and those solutions *are* tested, because most rarely-taken branches exist
for them: one step at a time, each replayed from XFOIL's exact state, so that no step inherits the
drift a whole run amplifies ([execution equivalence](../execution-equivalence/index.md)). Whether
the reference is ill-conditioned at a point is measured on the reference alone, by its
[twins](../../contributing/conventions/terminology.md#twins): five reruns with every panel
coordinate jogged by one ULP.
