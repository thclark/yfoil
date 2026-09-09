# Reference noise floor

**What this measures.** How much XFOIL's own outputs move when its *input* changes by the smallest
possible amount — every panel coordinate perturbed by +1 ULP (`math.nextafter`) — with everything else
identical. This is the floor below which no comparison between yFoil and XFOIL can be meaningful, and the
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
added to `xtask/fixtures-config/cases.toml` (the branch-coverage cases especially) before trusting a tolerance there,
and on any other host before comparing fixtures across hosts (libm differences, addendum R2).


## 2026-09-04 — current dumps, and the polar script

Re-measured after stages S5–S10 replaced the dump set (`scripts/noise-floor.sh`, same method:
every panel coordinate +1 ULP, instrumented XFOIL, line-aligned ES24.16 comparison).

### Single point (`naca0012_n60_a2_re1e6`, ALFA 2)

| file | lines | misaligned | max abs diff | max rel (abs>1e-8) | p99 rel | median rel | stage |
|---|---|---|---|---|---|---|---|
| xfoil_panels.dat | 71 | 0 | 7.30e-15 | 5.25e-14 | 2.45e-14 | 1.45e-16 | geometry after ABCOPY |
| xfoil_inviscid.dat | 129 | 0 | 2.73e-12 | 5.78e-11 | 7.02e-12 | 4.99e-14 | inviscid GAM/QINV |
| xfoil_dij.dat | 5628 | 0 | 9.78e-09 | 5.20e-06 | 2.40e-07 | 1.67e-11 | DIJ matrix |
| mrchdu_input_1.dat | 105 | 0 | 1.59e-11 | 2.16e-11 | 8.22e-12 | 1.68e-13 | BL state after MRCHUE init |
| mrchdu_output_1.dat | 160 | 0 | 1.59e-11 | 3.90e-11 | 8.49e-12 | 2.79e-13 | BL state after MRCHDU #1 |
| blsolv_input.dat | 17802 | 0 | 1.62e-04 | 8.85e-04 | 1.33e-06 | 4.87e-11 | SETBL VA/VB/VM/VDEL (calls 1..3) |
| blsolv_output.dat | 231 | 0 | 1.23e-10 | 3.18e-08 | 1.74e-08 | 5.10e-11 | BLSOLV VDEL out (calls 1..3) |
| update_output_1.dat | 198 | 0 | 5.98e-11 | 1.10e-10 | 1.87e-11 | 5.46e-13 | UPDATE state/RLX/RMSBL (call 1) |
| viscal_iter.dat | 120 | 0 | 6.46e-12 | 2.16e-10 | 1.68e-10 | 2.15e-13 | per-iteration CL/CD/RMSBL/RLX (VISCAL call 1) |
| viscal_final.dat | 169 | 0 | 1.43e-11 | 3.24e-10 | 6.67e-12 | 8.27e-14 | converged point (VISCAL call 1) |
| viscal_points.dat | 21 | 0 | 2.78e-13 | 2.16e-10 | 5.50e-12 | 1.01e-13 | per-point CL/CD/CM/XTR |

Both runs: 5 iterations, converged. Same floor as 2026-09-03.

### Polar script (`naca0012_n60_polar_re1e6`: ALFA 0 / ASEQ 1 5 1 / INIT / ALFA -1 / ASEQ -2 -5 -1)

| file | lines | misaligned | max abs diff | max rel (abs>1e-8) | p99 rel | median rel | stage |
|---|---|---|---|---|---|---|---|
| mrchdu_output_1.dat | 160 | 0 | 1.45e-11 | 2.41e-09 | 8.60e-11 | 5.65e-13 | BL state after MRCHDU #1 |
| blsolv_output.dat | 231 | 0 | 3.00e-10 | 7.53e-08 | 1.18e-08 | 2.87e-10 | BLSOLV VDEL out (calls 1..3) |
| update_output_1.dat | 198 | 0 | 6.09e-11 | 6.77e-09 | 8.34e-11 | 1.33e-12 | UPDATE state/RLX/RMSBL (call 1) |
| viscal_iter.dat | 144 | 0 | 1.02e-10 | 3.84e-07 | 3.77e-07 | 2.38e-12 | per-iteration record (VISCAL call 1) |
| viscal_final.dat | 169 | 0 | 3.98e-11 | 2.95e-10 | 1.75e-11 | 8.54e-14 | converged point (VISCAL call 1) |
| viscal_iters_all.dat | 73 | 0 | **5.12e-10** | 3.77e-07 | 1.12e-08 | 9.41e-13 | per-iteration record, all 11 VISCAL calls |
| viscal_points.dat | 231 | 0 | 1.19e-12 | 5.66e-09 | 2.31e-09 | 1.84e-13 | per-point CL/CD/CM/XTR (all calls) |

Per point, base vs +1 ULP (identical iteration counts and converged flags on all 11 points):

| call | alpha | iterations | \|dCL\| | \|dCD\| | \|dCM\| | \|dXTR1\| |
|---|---|---|---|---|---|---|
| 1 | 0 | 6 | 3.9e-14 | 2.0e-15 | 9.9e-15 | 9.8e-13 |
| 2 | 1 | 6 | 6.6e-15 | 2.7e-15 | 2.7e-15 | 5.3e-13 |
| 3 | 2 | 6 | 1.3e-14 | 1.1e-15 | 2.5e-16 | 2.3e-13 |
| 4 | 3 | 5 | 3.6e-15 | 1.4e-15 | 1.7e-15 | 3.3e-13 |
| 5 | 4 | 5 | 1.7e-14 | 1.3e-15 | 4.7e-15 | 2.0e-13 |
| 6 | 5 | 6 | 7.8e-14 | 2.0e-16 | 1.4e-14 | 5.5e-14 |
| 7 | −1 | 6 | 4.6e-14 | 2.0e-15 | 8.9e-15 | 5.4e-14 |
| 8 | −2 | 6 | 4.0e-14 | 5.4e-15 | 4.9e-15 | 4.2e-14 |
| 9 | −3 | 5 | 4.6e-14 | 3.1e-15 | 5.1e-15 | 1.1e-13 |
| 10 | −4 | 5 | 5.1e-14 | 3.0e-17 | 5.5e-15 | 9.6e-14 |
| 11 | −5 | 6 | 3.8e-14 | 5.4e-17 | 9.6e-15 | 6.9e-15 |

**Reading.** Converged points are fixed points of the Newton iteration: the 1-ULP perturbation
moves them by ≤ 8e-14 (forces) and ≤ 1e-12 (XTR) — `TOL_SOLVER` (1e-10) is generous for them.
The *path* to convergence is not a fixed point: intermediate RLX/RMSBL/CL move by up to
5.1e-10 absolute (relative spreads of 4e-7 are on values ~1e-7). Hence `TOL_TRANSIENT = 1e-9`
(floor × 2) for per-iteration comparisons inside multi-point sequences, with iteration counts,
converged flags, IST and ITRAN still compared exactly. yFoil's polar replay meets both: worst
per-iteration transient difference 1.5e-10 (RLX, point 4 iteration 2), all converged points
within `TOL_SOLVER`.


## 2026-09-04 — per-case floors (mechanised) and the coverage cases

`cargo xtask fixtures` now runs the +1-ULP twin of every case and records `noise_floor.json`
next to the fixtures: per VISCAL call the spread of the converged point, and per iteration the
spread of RMSBL/RLX/CL/CD/CM/ALFA/MINF/REINF (plus RMXBL and whether UPDATE's reported limiter
flipped, on call 1). Gates use `max(tol·scale, FLOOR_FACTOR · floor)` with `FLOOR_FACTOR = 4`.

| case | branch trace under 1 ULP | worst point spread | worst per-iteration spread | yFoil outcome |
|---|---|---|---|---|
| naca0012_n60_a2_re1e6 | identical | 2.3e-13 | 6.5e-12 | match |
| naca0012_n60_polar_re1e6 (11 calls) | identical | 1.2e-12 | 5.1e-10 | match |
| naca0012_n60_cl03_re1e6 (CL 0.3) | identical | 3.3e-13 | 2.8e-10 | match |
| naca0012_n60_a2_re1e6_type2 (TYPE 2) | identical | 3.0e-8 | 3.8e-6 | match |
| naca0012_n60_sharp_a2_re1e6 | identical | 2.6e-9 | 1.0 (reported limiter flips at a tie, iteration 3) | match (tie confirmed by the twin) |
| naca0012_n60_a2_re1e6_m03 (M 0.3) | identical | 8.7e-14 | 2.8e-11 | match |
| naca0012_n60_a12_re1e6 (12°, unconverged) | identical | 1.5e-3 | 4.6e-2 | **threshold-straddling at iteration 19** |
| naca0012_n60_a4_re1e5 (Re 1e5) | identical | 5.4e-14 | 8.5e-12 | match |
| naca0012_n60_a2_re1e6_xtr03 (XTR 0.3) | identical | 1.2e-14 | 8.8e-13 | match |

**The 12° case, in full.** Through iteration 17 yFoil differs from the reference by 1.5–2× the
reference's own 1-ULP spread at every iteration (e.g. 1.09e-8 vs 6.7e-9 on RMSBL at 17): it
behaves as a ~2-ULP perturbation of XFOIL, which is what a translation with identical branch
trace and libm should look like. At iteration 18 the reference becomes hypersensitive (its own
spread jumps from 6.7e-9 to 5.1e-5 on RMSBL and 7.7e-4 on RLX); yFoil is still at 1.7× the
floor there. At iteration 19 the runs part (RMSBL 1.87 vs 0.54). Replaying iterations 18 and
19 from XFOIL's dumped state (`dump_calls = [18, 19]`) reproduces XFOIL's own next state within
the floor (iteration 19: every BL array within 1.9e-13, RLX/RMSBL identical, same limiter),
so the step map is faithful and the departure is accumulated 2-ULP-level differences crossing
a threshold the reference itself cannot hold to better than 3.5e-4. Classified
threshold-straddling (`STRADDLE_FLOOR = 1e-6`), reported, not passed.

**Ties in UPDATE's bookkeeping.** VMXBL/IMXBL record which normalised change was largest; at
the sharp-TE case's iteration 3 the θ and δ* changes at the side-2 similarity station are equal
to 7 digits and the label flips — in the twin as well. RLX is a minimum over every variable and
is unaffected; such flips are accepted only when |RMXBL| agrees within the recorded floor or the
twin flips too.
