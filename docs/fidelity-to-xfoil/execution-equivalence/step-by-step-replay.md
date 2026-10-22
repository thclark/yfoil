# Step-by-step replay

The replay harness seeds yFoil with XFOIL's complete state at the start of one step and requires
XFOIL's state at its end. That isolates one iteration of one call from everything upstream, which
is how an ill-conditioned solution is tested at all: a step replayed from XFOIL's exact state is
far less likely to diverge than a whole run, so no step inherits the drift a whole run amplifies.
It is also the method for finding a translation bug
([finding a translation bug](../../contributing/finding-a-translation-bug.md)).

## One step

A step is one VISCAL iteration — iteration `i` of call `c`, global SETBL call `k` — or one
inviscid call. `tests/execution/step.rs` loads the inviscid side from the panels and the BL side
from XFOIL's dump entering that SETBL call (`mrchdu_input_<k>.dat`), runs SETBL → BLSOLV → UPDATE
and the rest of the iteration, and compares with XFOIL's state after it (`update_output_<k>.dat`,
`viscal_points.dat`):

- the step's **decisions exactly**: IST, ITRAN, RLX, convergence and the iteration count;
- each **value** within the named tolerance (`TOL_SOLVER`, on every host) **or within four times
  the step's own sensitivity**, whichever is larger.

## The step's own sensitivity

The sensitivity is measured in the test: the same step is replayed from jogged copies of its
inputs — every panel coordinate, from which yFoil builds the inviscid side, and every BL entry of
XFOIL's state moved by −1, 0 or +1 ULP — and the furthest a jog moves a value is how much the step
amplifies a last-bit difference. An ill-conditioned solution amplifies within one step; the test
allows for exactly that amplification and no more. That is how an ill-conditioned solution is
tested without a tolerance asserted for it.

## The step tests

| Test | Steps |
| --- | --- |
| `tests/execution/steps.rs` | every iteration of the VISCAL calls a case names in `step_calls` — the ill-conditioned solutions: the sharp trailing edge, α = 12° separation, the TYPE 2 clamp, fixed CL far from alpha, the stalled first march at Re 3e6, deep stall, Re 1e4 with `DAMP`, and the [polar break points](polar-break-points.md) |
| `tests/execution/branches.rs` | the fewest steps that take every branch any step takes, each one whose route yFoil is observed to follow ([which test covers each branch](../branch-coverage/which-test-covers-each-branch.md)) |
| `tests/execution/events.rs` | the non-finite probes, one rarely-taken branch at a time, compared with XFOIL's own branch event ([non-finite probes](../branch-coverage/non-finite-probes.md)) |
| `tests/execution/viscal.rs` | the VISCAL loop on the reference case, seeded entering the first SETBL and run to convergence, every iteration's RMSBL, RLX, forces, IST and ITRAN compared |
| `tests/execution/update.rs` | UPDATE as the last step of VISCAL iterations 1 to 3 on the reference case, with its diagnostics (RLX, RMSBL, the largest change and where) |
