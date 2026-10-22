# Whole-run comparison

A whole run is a run test only if it takes XFOIL's route through the whole OPER script with every
value in tolerance to the end. There is no third outcome: a run whose route differs from XFOIL's
anywhere is a divergent comparison, not a run test, and such a case is tested by single steps
instead ([step-by-step replay](step-by-step-replay.md)).

## What a run test compares

`tests/execution/run.rs` drives yFoil through the case's own OPER script — every `ALFA`, `CL`,
`ASEQ` point and `INIT`, in one session — on the case's panels, and compares, from
`viscal_points.dat` and `viscal_iters_all.dat`:

- per VISCAL call: the iteration count and convergence, **exactly**;
- per iteration: IST and ITRAN exactly, and RMSBL, RLX, CL, CD, CM, ALFA, MINF and REINF within
  `TOL_SOLVER` for a run's first call and `TOL_TRANSIENT` for the calls after it, which carry the
  previous point's state;
- the converged point — ALFA, CL, CM, CD, CDF, CDP, the transition points, MINF, REINF, with IST
  and ITRAN exactly — within `TOL_SOLVER`.

An inviscid run compares every SPECAL and SPECCL point. The tolerances hold on every host, since
both codes take their transcendentals from the same correctly rounded functions
([correctly rounded maths](../the-reference-xfoil-build/correctly-rounded-maths-in-both-codes.md)).

## The run tests

`tests/execution/runs.rs` holds one test per case marked `run = true` (or `run_through = N`) in
`cases.toml`: the CI reference case, the polar, the fixed-CL and TYPE 2/TYPE 3 cases, M = 0.3,
laminar separation at Re 1e5, forced transition, a trip in the natural-transition interval and one
upstream of the stagnation point, a repeated alpha, `DAMP`, the inviscid Kármán–Tsien and TYPE 2 cases, the analytic section, and the
first calls of the [polar break points](polar-break-points.md).

On the reference case `naca0012_n60_a2_re1e6` the run is also checked from the first stage up,
everything yFoil's own:

| Test | From | To |
| --- | --- | --- |
| `tests/execution/prologue_pointers.rs` | yFoil's panels | TECALC → QISET → STFIND → IBLPAN → XICALC → IBLSYS → UICALC, against `xfoil_pointers.dat` / `xfoil_uinv.dat` |
| `tests/execution/prologue_dij.rs` | yFoil's panels | GGCALC → QISET → SPECAL's GAM → XYWAKE → QWCALC → QDCALC, against QINVU and DIJ |
| `tests/execution/analysis.rs` | yFoil's panels (`panels.json`) | the inviscid solve, wake, pointers, DIJ and the VISCAL loop to the converged point (`viscal_final.dat`, `viscal_iter.dat`); everything upstream of the Newton loop is yFoil's, so this is the propagated floor: the wake and DIJ at about 5e-11 bound the final forces |

The polar's own run test, with its reinitialisation, is
[polars and reinitialisation](polars-and-reinitialisation.md).
