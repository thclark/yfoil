# Execution equivalence

[Subroutine equivalence](../subroutine-equivalence/index.md) shows each routine right in
isolation. This section shows them right *together*: yFoil executing what XFOIL executed — a whole
OPER session, or one iteration of one — taking XFOIL's route through the code, with its values
within tolerance. These are the `execution` tests (`tests/execution/`).

A nonlinear iterative solver cannot be judged from the end. Small early differences compound, and
where the solution is ill-conditioned the reference cannot even reproduce itself under a 1-ULP
change of its input ([tolerances](../tolerances/index.md#divergent-comparisons-and-ill-conditioned-solutions)).
So there are two ways to compare an execution, and each case uses the one its conditioning allows:

- a **whole run**, for a run that takes XFOIL's route to the end with every value in tolerance;
- a **step**, one iteration replayed from XFOIL's exact dumped state, for an ill-conditioned
  solution — so that no step inherits the drift a whole run amplifies.

## How each case is tested

`xtask/fixtures-config/cases.toml` records, beside every reference case, how it is tested and why
(`# how tested:`), in fields the tests and the fixture generator both read:

- `run = true` — the whole OPER script is compared with XFOIL's run
  ([whole-run comparison](whole-run-comparison.md)); `run_through = N` compares the first N calls
  of a run whose later calls are tested as steps.
- `step_calls = [...]` — every iteration of these VISCAL calls is replayed as a single step from
  XFOIL's exact state ([step-by-step replay](step-by-step-replay.md)).
- `dump_calls` with `events = true` — the non-finite probes, one step per rarely-taken branch,
  compared with XFOIL's branch events ([non-finite probes](../branch-coverage/non-finite-probes.md)).

`tests/apparatus/case_tests.rs` checks that every case so marked has its test and every such test
names a case so marked, so the record and the tests cannot drift apart.

## Pages of execution equivalence

| Page | What it covers |
| --- | --- |
| [Whole-run comparison](whole-run-comparison.md) | Well-conditioned runs compared from the panels to the last iteration |
| [Step-by-step replay](step-by-step-replay.md) | One iteration seeded from XFOIL's exact state, and the step's own sensitivity |
| [Polars and reinitialisation](polars-and-reinitialisation.md) | The sweep state machine, and the re-solved 0° that seeds the downward leg |
| [Polar break points](polar-break-points.md) | Sweeps past maximum lift: a whole run up to the break, then steps |
| [Validity flags against XFOIL's events](validity-flags-against-xfoils-events.md) | yFoil's invalid-point flags checked against the reference's recorded branch events |

Which steps are replayed so that every branch XFOIL takes is tested is
[branch coverage](../branch-coverage/index.md). The same comparison, applied to whole sweeps under
normal conditions and read as polars, is the engineer's view under Engineering validation:
[agreement with XFOIL under normal conditions](../../engineering-validation/agreement-with-xfoil-under-normal-conditions.md).
That page demonstrates; the tests on this section check.
