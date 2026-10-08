# Testing

What each test is for, how tests are organised, and where their fixtures come from. Terms such
as *ill-conditioned*, *twin* and *noise floor* are defined in [terminology.md](terminology.md).

## Every test has exactly one purpose

Every test performs one, and only one, of the following functions.

| | Category | Function | Compared against | Typical shape | Binary |
|---|---|---|---|---|---|
| **a** | Subroutine equivalence | Verify that a well-encapsulated subroutine/function of yFoil produces identical-within-tolerance output to an equivalent subroutine in XFOIL. Necessary to test the numerical correctness of the basic input blocks. | XFOIL's dumped inputs → outputs | a pure call | `subroutine` |
| **b** | Execution equivalence | Verify that an execution of yFoil (e.g. for an iteration from a given starting point) either matches XFOIL's, or returns failure where XFOIL is ill-conditioned. This is where branch coverage is key, because we need to check all conditions that XFOIL will have encountered. It's also where we have to run a study in order to generate test fixtures, because in some cases the only way we can know if XFOIL is ill-conditioned is by looking at the outputs of perturbed runs. This has required many cases (the branch coverage and twins studies, for us to *discover* which cases actually fulfil all these requirements and what the particular tolerance should be in each case) but should be able to be distilled down to particular narrow test cases replicating XFOIL's behaviour on that branch. | XFOIL step dumps; run records for well-conditioned whole runs | one step seeded from XFOIL's state; short whole runs | `execution` |
| **c** | yFoil functionality | Verify that a well-encapsulated functionality of yFoil not equivalent to anything in XFOIL succeeds with correct results or fails in the expected way (e.g. file writers/parsers, geometry routines, CLI commands). This is standard unit/integration testing for the yFoil program. | the specification; external non-XFOIL references (PDAS naca456); regression snapshots | ordinary unit/integration tests | `application` |
| **d** | Known XFOIL weaknesses | Verify (and define) deterministic behaviour of yFoil in all special cases (other than the routine branch coverage above) where XFOIL's behaviour is noted to be incorrect or ill-conditioned (tracked in [`docs/xfoil-known-issues.md`](../xfoil-known-issues.md)). Ideally this is a correct solution, where case flags have been given that allow solutions to vary from XFOIL; otherwise it is a well-defined failure per the [solution validity](../guide/validity.md) contract. This is what allows us to be confident that we're handling the worst weaknesses and corner cases found in XFOIL along the way. | the known-issue entry; XFOIL's branch events | one test module per known-issue section | `known_issues` |
| **e** | Physical invariants | Verify that yFoil satisfies properties that hold independently of XFOIL: a symmetric aerofoil at α = 0 gives CL = CM = 0, a mirrored aerofoil at −α gives the mirrored solution, and a flat plate gives the Blasius solution. Two codes can share a misunderstanding; these tests check the physics rather than the translation. | physics | yFoil only | `invariants` |
| **f** | Reference integrity | Verify the reference and its fixtures rather than yFoil: the bitwise geometry handoff, inert instrumentation, recorded provenance, and the reference agreeing with itself (e.g. the polar's two 0° solves are identical). If these fail, no comparison in a–b can be trusted. | the reference itself | checks on fixtures | `apparatus` |

## Rules

1. **One category per test.** A test that does two of these jobs is two tests.
2. **A test never accepts a divergence.** For (b), the branch decisions must be identical to
   XFOIL's, and each value must agree within one of the named tolerances in
   `tests/common/utilities/tolerances.rs` — or, for a single step, within four times that step's own
   sensitivity, whichever is larger. The sensitivity is measured in the test: the same step is
   replayed from three copies of its inputs (the panels and XFOIL's dumped state) jogged by −1, 0 or
   +1 ULP, and the furthest a jog moves a value is how much the step amplifies a last-bit
   difference. That is how an ill-conditioned solution is tested without a tolerance asserted for
   it. There is no third outcome in a test.
3. **A (b) test is never a divergent comparison.** Its purpose is to prove that yFoil takes the
   same route through the code as XFOIL, so XFOIL and yFoil must take the same branches through the
   step it replays (`terminology.md`). That is a property of the comparison, decided by observing
   both codes' routes, not by how sensitive the solution is.
4. **Ill-conditioned solutions are tested, not avoided.** Most of XFOIL's rarely-taken branches
   exist for ill-conditioned solutions (stall, separation, the Kármán–Tsien pole), and branch
   coverage needs them. A step replayed from XFOIL's exact state is far less likely to diverge than
   a whole run, so such cases are tested as single steps wherever possible.
5. **Divergence belongs to the studies.** Explaining that an XFOIL–yFoil difference in a polar is
   a divergent comparison, not an error, is the job of the studies and of the validation reports —
   not of a test. No test reads twin data.
6. **Studies demonstrate; tests gate.** A study shows an external reader how a condition, regime
   or behaviour arises and whether yFoil follows the reference through it; a test gates correct,
   consistent functionality as minimally and quickly as it can. So one behaviour may need both: a
   minimal test (a step or an event) that gates it, and a study run that shows it arising — a
   `series` case for a section or regime under normal conditions, a `known-issues` case for one of
   XFOIL's known issues (`docs/xfoil-known-issues.md`). A run whose only purpose is a known issue is
   a `known-issues` case, never a `series` one. Every validation run — each `branch-coverage`,
   `non-finite` and `series` case — carries its own five 1-ULP twins (`twins = true`), at its own
   panelling, which is the only way to know whether XFOIL was ill-conditioned at a given point; a
   twin of a different case, or of the same section at another panel count, says nothing about it.
7. **A test asserts something.** A diagnostic that only prints is not a test; it goes in
   `examples/attic/`.
8. **A missing fixture fails the test** (`require_fixture`); it is never skipped.

## Layout

One test binary per category, so the category of a test is the directory it is in:

```
tests/
  subroutine/main.rs      (a)
  execution/main.rs       (b)
  application/main.rs     (c)
  known_issues/main.rs    (d)   one module per known-issues section
  invariants/main.rs      (e)
  apparatus/main.rs       (f)
  common/                 shared helpers and fixture loaders
```

`cargo test --test execution` runs one category. Each module starts with a `//!` header that
states its category, what is compared with what, which fixture family it reads, and the command
that regenerates that family.

Unit tests inside `src/` are limited to (c), (e) and the mathematics of generic routines
(splines, Gauss elimination, LU). A comparison with XFOIL belongs in `tests/subroutine/`, against
a fixture with a manifest, never against literals typed into the source.

## How each case is tested

`xtask/fixtures-config/cases.toml` records, beside every reference case, how it is tested and why
(`# how tested:`), in fields the tests and the fixture generator both read:

- `run = true` — the whole OPER script is compared with XFOIL's run (`tests/execution/run.rs`):
  iteration counts, convergence, IST and ITRAN exactly, every iteration's values and the converged
  point within the named tolerances. Only a run that takes XFOIL's route with values in tolerance to
  the end is a run test; `run_through = N` compares the first N calls of a run whose later calls
  are tested as steps. An inviscid run compares every SPECAL/SPECCL point.
- `step_calls = [...]` — every iteration of these VISCAL calls is replayed as a single step from
  XFOIL's exact state (`tests/execution/step.rs`, `steps.rs`). Ill-conditioned solutions are tested
  this way, so that no step inherits the drift a whole run amplifies.
- `dump_calls` with `events = true` — the non-finite probes, one step per rarely-taken branch,
  compared with XFOIL's branch events (`tests/execution/events.rs`).

## How branch coverage is established

1. **Which step takes which branch.** `cargo xtask steps` reruns each candidate case on the
   gcov build of the reference, cut off after each step, and attributes to a step the branches whose
   counts rose over the previous step's prefix (`target/coverage/steps.json`).
2. **The route.** `cargo xtask route` reruns each candidate case on the instrumented reference
   with every SETBL call dumped, replays every candidate step in a build of the `execution` tests
   with LLVM source coverage, and compares, step by step, the call count of every translated
   subroutine (the functions carrying its `#[doc(alias)]`) with XFOIL's for the same step, which
   step 1 also records. The replay's profile windows are lined up with XFOIL's prefix-run
   differences as `tests/execution/route.rs` describes. yFoil makes no call at some of XFOIL's call
   sites: work it does inline (SETBL's COMSET), and calls whose result nothing reads, which it does
   not translate (SETBL's trailing-edge BLMID(3), the BLMID calls of the march fallbacks). Each
   such site is listed with its reason in `xtask/fixtures-config/route-map.toml`, and the calls
   XFOIL made from it in the step — gcov counts them per line — are taken off XFOIL's counts,
   with the calls the dropped routine makes in turn, before comparing. yFoil's code is not
   changed to suit the measurement. Steps whose counts differ are **divergent** and
   are listed in `xtask/fixtures-config/route.toml` with the subroutines that differ, beside
   divergences found by other means (`found`, with the evidence).
3. **The step cover** (`terminology.md`). The fewest steps taking every branch that any step whose route agrees takes
   (`xtask/fixtures-config/step-cover.toml`), each replayed by one test in
   `tests/execution/branches.rs`; the fixture generator adds the dumps and files those steps need.
   Every step is eligible, ill-conditioned solutions included, except the divergent ones; the
   twins' conditioning of each step is recorded beside it for information.
4. **The record.** `docs/validation/branch-gating.md` (generated) lists every candidate case and
   how it is tested, every chosen step with its branches, and the branches no test gates.

## What "gated" means

A branch is **gated** when a test replays a step in which XFOIL took the branch, compares the step's
decisions exactly and its values as above, **and** yFoil's own route through that step is observed
to equal XFOIL's (`cargo xtask route`). The route is observed at the granularity of subroutine call
counts: every translated subroutine is called the same number of times in the step by both codes.
That detects a divergent route wherever it changes what is called — a different transition
station, a different number of TRCHEK2 or station-Newton iterates, a fallback taken or not — and
the step's compared decisions (IST, ITRAN, the iteration count, RLX, convergence) cover the
decisions those reach. A branch inside a subroutine that changes no call count and no compared
decision is gated only through the step's values; observing it directly needs a line map between
the Fortran and the Rust, which is not built yet.

Measured 2026-10-06 over 870 steps of 36 candidate cases. The first measurement found 6 divergent
steps, every one a different number of TRCHEK2 Newton iterates; they were one yFoil bug, not a
knife edge (`docs/xfoil-known-issues.md` §7.12: AMI was not set from AMPL2 after a forced
transition, which changed no converged value and so no value test), and after the fix every
measured step's route agrees. One step remains excluded, found from the branch events rather than
the call counts (MRCHDU's fallback at a different station with the same counts, `route.toml`).
Every branch any candidate step takes is taken by a step whose route agrees, so none is left
ungated by the exclusion.

## Status

This convention was adopted on 2026-10-06 and every test lives in its category binary. One exception
remains: `subroutine/pane_legacy.rs` compares against panels XFOIL generated itself (CLAUDE.md Rule 4);
PANGEN is gated against yFoil-generated panels by `subroutine/pangen.rs`.
