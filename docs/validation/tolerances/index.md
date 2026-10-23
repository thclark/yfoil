# Tolerances

--8<-- "validation/tolerances/tables/run.md"

!!! info "Used by the test suite"

    The four tolerance constants every test is held to — `TOL_PURE`, `TOL_LINALG`, `TOL_SOLVER` and
    `TOL_TRANSIENT`, in `tests/common/utilities/tolerances.rs` — were set by hand from this study's
    findings. No test reads the study's output, so rerunning it changes nothing until the constants
    are revisited by hand. How the tests use them is in [Testing](../../development/testing/README.md).

**What this measures.** How much XFOIL's own outputs move when its *input* changes by the smallest
possible amount — every panel coordinate perturbed by +1 ULP (`math.nextafter`) — with everything else
identical. This is the floor below which no comparison between yFoil and XFOIL can be meaningful
(CLAUDE.md Rule 1).

**How.** The tolerances study, `studies/tolerances/tolerances.sh`, takes the tracked reference case
(`tests/fixtures/xfoil/naca0012_n60_a2_re1e6`: NACA 0012, N=60, α=2°, Re=1e6, M=0, Ncrit=9, ITER 20),
runs the instrumented double-precision reference twice (base and perturbed), and compares every
`ES24.16` dump line-by-line. Metric per value: `|a−b|` and, for `max(|a|,|b|) > 1e-8`, `|a−b| / max(|a|,|b|)`.
Both runs converged in the same number of iterations (5) with the same branch trace, so this is the
*continuous* floor — branch flips are a separate phenomenon (plan addendum R3).

## Current measurement

The latest published run of the study (`cargo xtask publish tolerances`). The study writes these tables
into its run folder, and this page includes them.

!!! note
    These tables were produced by an earlier version of the script, whose list of dumps differs from
    the current one's; they are replaced at the study's next run.

--8<-- "validation/tolerances/tables/floor.md"

--8<-- "validation/tolerances/tables/iterations.md"

## In this section

| Page | What it shows |
|---|---|
| [Measurement history](measurement-history.md) | Every earlier measurement of the floor, as written at the time |
