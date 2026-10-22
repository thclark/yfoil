# The test suite

What each test is for, how tests are organised, and what every test must do. Terms such
as *ill-conditioned*, *twin* and *noise floor* are defined in the
[terminology](conventions/terminology.md). What the tests show is under Fidelity to XFOIL and
Engineering validation; this page is the rules for writing them.

## The six test purposes

Every test performs one, and only one, of the following functions.

| | Category | Function | Compared against | Typical shape | Binary |
|---|---|---|---|---|---|
| **a** | Subroutine equivalence | Verify that a well-encapsulated subroutine/function of yFoil produces identical-within-tolerance output to an equivalent subroutine in XFOIL. Necessary to test the numerical correctness of the basic input blocks. | XFOIL's dumped inputs → outputs | a pure call | `subroutine` |
| **b** | Execution equivalence | Verify that an execution of yFoil (e.g. for an iteration from a given starting point) either matches XFOIL's, or returns failure where XFOIL is ill-conditioned. This is where branch coverage is key, because we need to check all conditions that XFOIL will have encountered. It's also where we have to run a study in order to generate test fixtures, because in some cases the only way we can know if XFOIL is ill-conditioned is by looking at the outputs of perturbed runs. This has required many cases (the branch coverage and twins studies, for us to *discover* which cases actually fulfil all these requirements and what the particular tolerance should be in each case) but should be able to be distilled down to particular narrow test cases replicating XFOIL's behaviour on that branch. | XFOIL step dumps; run records for well-conditioned whole runs | one step seeded from XFOIL's state; short whole runs | `execution` |
| **c** | yFoil functionality | Verify that a well-encapsulated functionality of yFoil not equivalent to anything in XFOIL succeeds with correct results or fails in the expected way (e.g. file writers/parsers, geometry routines, CLI commands). This is standard unit/integration testing for the yFoil program. | the specification; external non-XFOIL references (PDAS naca456); regression snapshots | ordinary unit/integration tests | `application` |
| **d** | Known XFOIL weaknesses | Verify (and define) deterministic behaviour of yFoil in all special cases (other than the routine branch coverage above) where XFOIL's behaviour is noted to be incorrect or ill-conditioned (tracked in [the issue register](../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md)). Ideally this is a correct solution, where case flags have been given that allow solutions to vary from XFOIL; otherwise it is a well-defined failure per the [validity contract](../using-yfoil/whether-to-trust-a-result.md#the-validity-contract). This is what allows us to be confident that we're handling the worst weaknesses and corner cases found in XFOIL along the way. | the known-issue entry; XFOIL's branch events | one test module per known-issue section | `known_issues` |
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
   step it replays ([terminology](conventions/terminology.md#divergent-comparison)). That is a property of the comparison, decided by observing
   both codes' routes, not by how sensitive the solution is.
4. **Ill-conditioned solutions are tested, not avoided.** Most of XFOIL's rarely-taken branches
   exist for ill-conditioned solutions (stall, separation, the Kármán–Tsien pole), and branch
   coverage needs them. A step replayed from XFOIL's exact state is far less likely to diverge than
   a whole run, so such cases are tested as single steps wherever possible.
5. **Divergence belongs to the studies.** Explaining that an XFOIL–yFoil difference in a polar is
   a divergent comparison, not an error, is the job of the studies and of the validation reports —
   not of a test. No test reads twin data.
6. **Studies demonstrate; tests check.** A study shows an external reader how a condition, regime
   or behaviour arises and whether yFoil follows the reference through it; a test checks correct,
   consistent functionality as minimally and quickly as it can. So one behaviour may need both: a
   minimal test (a step or an event) that tests it, and a study run that shows it arising — a
   `series` case for a section or regime under normal conditions, a `known-issues` case for one of
   XFOIL's known issues ([the issue register](../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md)). A run whose only purpose is a known issue is
   a `known-issues` case, never a `series` one. Every study run — each `branch-coverage`,
   `non-finite` and `series` case — carries its own five 1-ULP twins (`twins = true`), at its own
   panelling, which is the only way to know whether XFOIL was ill-conditioned at a given point; a
   twin of a different case, or of the same section at another panel count, says nothing about it.
7. **A test asserts something.** A diagnostic that only prints is not a test, and is not
   committed: it is written under `.tmp/` and discarded.
8. **A missing fixture fails the test** (`require_fixture`); it is never skipped.

## Where the evidence is shown

The reports the tests and studies produce live in the tab of the question they answer, not in a
bucket of their own: whether yFoil is fit for engineering work — the generators against naca456,
the physical invariants, agreement with XFOIL on the series cases under normal conditions — is
[Engineering validation](../engineering-validation/index.md); whether it is a faithful translation
— the reference, the fixtures, the tolerances, subroutine and execution equivalence, branch
coverage, the known-issues runs — is [Fidelity to XFOIL](../fidelity-to-xfoil/index.md). Every
report is generated ([regenerating the reports](regenerating-the-reports.md)).

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

How each reference case is tested — as a whole run, step by step, or as events — is recorded
beside it in `xtask/fixtures-config/cases.toml` and described in
[execution equivalence](../fidelity-to-xfoil/execution-equivalence/index.md#how-each-case-is-tested);
how full branch coverage is established, and what "tested" means for a branch, is
[branch coverage](../fidelity-to-xfoil/branch-coverage/index.md).

## Testing yFoil's own features

The `application` tests (c) check what XFOIL does not have, against its specification or an
external reference rather than against XFOIL: the CLI end to end through the binary
(`cli_geometry.rs`, `cli_plot.rs`), geometry generation, I/O, repanelling and input validation
(`geometry.rs`, `geometry_errors.rs`, `panelling.rs`), yFoil's own cosine repanelling frozen as a
regression snapshot (`repanel_cosine.rs`), the series generators against naca456
(`naca456.rs`, shown in
[generated sections against naca456](../engineering-validation/generated-sections-against-naca456.md)),
the operating-point commands' contracts (`fixed_cl.rs`), the polar driver (`polar_driver.rs`), the
output and validity records (`output.rs`, `validity_record.rs`) and the content-addressed point
ids (`ids.rs`).

## Missing fixtures and ignored tests

- **A missing fixture fails the test** (`require_fixture`); it is never `eprintln!("Skipping")`
  and green ([tracked and regenerated fixtures](../fidelity-to-xfoil/fixtures.md#tracked-and-regenerated-fixtures)).
- **An `#[ignore]` is tagged with the stage that unblocks it** and listed in `tests/IGNORED.txt`;
  the `ignore-drift` CI job (`ci/ignore-drift.sh`) fails when the two disagree, and the `unit` job
  runs with zero ignores.
- **A test with no assertion is not a test.** A diagnostic that only prints is written under
  `.tmp/` and discarded, never committed.

## Status of the convention

This convention was adopted on 2026-10-06 and every test lives in its category binary. One exception
remains: `subroutine/pane_legacy.rs` compares against panels XFOIL generated itself (CLAUDE.md Rule 4);
PANGEN is tested against yFoil-generated panels by `subroutine/pangen.rs`.
