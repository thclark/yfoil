# Testing

What each test is for, how tests are organised, and where their fixtures come from. Terms such
as *ill-conditioned*, *twin* and *noise floor* are defined in [terminology.md](terminology.md).

## Every test has exactly one purpose

Every test performs one, and only one, of the following functions.

| | Category | Function | Compared against | Typical shape | Binary |
|---|---|---|---|---|---|
| **a** | Subroutine equivalence | Verify that a well-encapsulated subroutine/function of yFoil produces identical-within-tolerance output to an equivalent subroutine in XFOIL. Necessary to test the numerical correctness of the basic input blocks. | XFOIL's dumped inputs → outputs | a pure call | `subroutine` |
| **b** | Execution equivalence | Verify that an execution of yFoil (e.g. for an iteration from a given starting point) either matches XFOIL's, or returns failure where XFOIL is ill-conditioned. This is where branch coverage is key, because we need to check all conditions that XFOIL will have encountered. It's also where we have to run a study in order to generate test fixtures, because in some cases the only way we can know if XFOIL is ill-conditioned is by looking at the outputs of perturbed runs. This has required many cases (the branch coverage and twins studies, for us to *discover* which cases actually fulfil all these requirements and what the particular tolerance should be in each case) but should be able to be distilled down to particular narrow test cases replicating XFOIL's behaviour on that branch. | XFOIL step dumps; run records for well-conditioned whole runs | one step seeded from XFOIL's state; short whole runs | `execution` |
| **c** | yFoil functionality | Verify that a well-encapsulated functionality of yFoil not equivalent to anything in XFOIL succeeds with correct results or fails in the expected way (e.g. file writers/parsers, geometry routines, CLI commands). This is standard unit/integration testing for the yFoil program. | the specification; external non-XFOIL references (PDAS naca456); regression snapshots | ordinary unit/integration tests | `yfoil` |
| **d** | Known XFOIL weaknesses | Verify (and define) deterministic behaviour of yFoil in all special cases (other than the routine branch coverage above) where XFOIL's behaviour is noted to be incorrect or ill-conditioned (tracked in [`docs/xfoil-known-issues.md`](../xfoil-known-issues.md)). Ideally this is a correct solution, where case flags have been given that allow solutions to vary from XFOIL; otherwise it is a well-defined failure per the [solution validity](../guide/validity.md) contract. This is what allows us to be confident that we're handling the worst weaknesses and corner cases found in XFOIL along the way. | the known-issue entry; XFOIL's branch events | one test module per known-issue section | `known_issues` |
| **e** | Physical invariants | Verify that yFoil satisfies properties that hold independently of XFOIL: a symmetric aerofoil at α = 0 gives CL = CM = 0, a mirrored aerofoil at −α gives the mirrored solution, and a flat plate gives the Blasius solution. Two codes can share a misunderstanding; these tests check the physics rather than the translation. | physics | yFoil only | `invariants` |
| **f** | Reference integrity | Verify the reference and its fixtures rather than yFoil: the bitwise geometry handoff, inert instrumentation, recorded provenance, and the reference agreeing with itself (e.g. the polar's two 0° solves are identical). If these fail, no comparison in a–b can be trusted. | the reference itself | checks on fixtures | `apparatus` |

## Rules

1. **One category per test.** A test that does two of these jobs is two tests.
2. **A test never accepts a divergence.** A test passes when values agree within one of the
   named tolerances in `tests/utilities/tolerances.rs` and, for (b), the branch decisions are
   identical to XFOIL's. There is no third outcome in a test.
3. **(b) test cases are well-conditioned by construction.** The purpose of a (b) test is to prove
   that yFoil takes the same code path as XFOIL, so the step it replays must be one where a
   one-ULP change cannot flip a decision. A step that is ill-conditioned is not a test case.
4. **The studies select the test cases; the tests do not read study data.** Whether a step is
   well-conditioned is decided once, when its fixture is chosen, from the twins: every twin must
   reproduce the reference's branch trace through the step, and the twins' spread at the step
   must be far below the tolerance. A branch reachable only through ill-conditioned steps is
   reported by the coverage study, and no test is written for it.
5. **Ill-conditioned behaviour belongs to the studies.** Explaining that an XFOIL–yFoil
   difference is a rounding-induced change of branch, not an error, is the job of the twins and
   sensitivity studies and of the validation reports — not of a test.
6. **A test asserts something.** A diagnostic that only prints is not a test; it goes in
   `examples/attic/`.
7. **A missing fixture fails the test** (`require_fixture`); it is never skipped.

## Layout

One test binary per category, so the category of a test is the directory it is in:

```
tests/
  subroutine/main.rs      (a)
  execution/main.rs       (b)
  yfoil/main.rs           (c)
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

## Status

This convention was adopted on 2026-10-06. The existing flat `tests/*_tests.rs` files are being
migrated into the category binaries; until that is complete, some files still mix categories and
some whole-run comparisons still read the noise floor (`tests/utilities/records.rs`).
