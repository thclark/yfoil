# Finding a translation bug

VISCAL is a nonlinear iterative solver, and **it cannot be debugged from the end.** Small early
differences compound; the final CD or CL says nothing about where a bug is. Fixing one bug may make
the end result *worse*, because several bugs can cancel. If you are looking at final forces to
understand a difference, go back to the second step below.

## When something is wrong

1. **Stop.** Do not tweak parameters, and do not add a workaround.
2. **Isolate the smallest unit that could be wrong**: a subroutine, a statement, one operation.
3. **Instrument XFOIL there**, rebuild the reference and capture a fixture
   ([instrumentation patches and logging](../fidelity-to-xfoil/the-reference-xfoil-build/instrumentation-patches-and-logging.md#adding-new-instrumentation)).
4. **Write the test against the fixture** ([the test suite](the-test-suite.md)).
5. **Fix yFoil to match**, and check that the branch trace matches too.
6. **Advance** to the next computation.

Progress is measured by fixture matches, never by an end-to-end number moving in the right
direction.

## Forward-stepping through XFOIL's execution

1. Instrument XFOIL at the exact point under investigation, as a tracked patch.
2. Rebuild the reference (`cargo xtask xfoil-build --verify`, which also proves the new
   instrumentation inert).
3. Run both codes on identical inputs at full precision: `ES24.16` on the Fortran side, `{:.17e}`
   on the Rust side, on panels yFoil generated and XFOIL loaded
   ([the bitwise geometry handoff](../fidelity-to-xfoil/fixtures.md#the-bitwise-geometry-handoff)).
4. Compare at that point. Values outside tolerance with an identical branch trace are a bug: fix
   it. A different branch at a threshold whose operands are within rounding is a divergent
   comparison and says nothing about the translation
   ([tolerances](../fidelity-to-xfoil/tolerances/index.md#divergent-comparisons-and-ill-conditioned-solutions)).
5. Advance to the next computation **in XFOIL's execution order** — not panel order, not module
   order.

Prefer the **replay harness**: seed yFoil with XFOIL's complete state at iteration *k* and require
XFOIL's state at *k + 1* ([step-by-step replay](../fidelity-to-xfoil/execution-equivalence/step-by-step-replay.md)).
That isolates one iteration of one subroutine from everything upstream. For a single subroutine,
capture its inputs and outputs and test it in isolation
([subroutine equivalence](../fidelity-to-xfoil/subroutine-equivalence/index.md)). To find the Rust
that translates a given Fortran routine, search `cargo doc` by the Fortran name (every translating
function carries `#[doc(alias = "...")]`) or the [translation map](../fidelity-to-xfoil/translation-map/index.md).

## Afterwards

If the bug was XFOIL's — a quirk yFoil must replicate, a dated fix, a dead branch — register it in
[the issue register](../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md). If a
variable or routine was added or renamed, update
[variable and routine names](../fidelity-to-xfoil/translation-map/variable-and-routine-names.md).
