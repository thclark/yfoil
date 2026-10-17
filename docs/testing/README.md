# Testing

How the test suite tracks yFoil against XFOIL, for software developers and for aerodynamicists deep into
XFOIL's edge cases. The engineering assessment of yFoil against XFOIL on normal cases is
[Validation](../validation/README.md); what each test is for is `docs/conventions/testing.md`.

Every page here is generated from data the fixture pipeline produced (the reference's dumps under
`tests/fixtures/xfoil/<case>/` or `target/fixtures/<case>/`), and can be regenerated at any time from
tracked inputs. The tests themselves are the evidence (`docs/conventions/testing.md`); the study pages
measure the reference and compare whole runs, and add nothing the tests rely on.

Figures are [ironlab](https://github.com/thclark/ironlab) figures, built in Rust by each study from the numbers it
computes (`docs/conventions/figures.md`). They are embedded live in these pages: pan, zoom and toggle lines from the
legend. Each carries parameters (case, conditions, outputs) and labels, so every study's figures can be loaded
together into the ironlab viewer and sorted, filtered and grouped. A study's run folder also holds a PDF of each. The branch-coverage study's branch map and case matrix are the exception: they are SVGs drawn by
matplotlib (`scripts/branch-coverage/plot.py` through `scripts/figures/render.sh`, which needs
[`uv`](https://docs.astral.sh/uv/) installed).

| Page | What it shows | Regenerate with |
|---|---|---|
| [coverage.md](coverage.md) | Rule 6: gcov branch completeness of the reference over every tracked case, with the open/unreachable/loop-entry classification | `cargo xtask coverage` |
| [branch-tests.md](branch-tests.md) | Which test covers each branch: every candidate case and how it is tested (run or single steps), which step takes which branch, the fewest steps that take them all (each replayed by a test), and the branches no test covers | `cargo xtask steps`, `cargo xtask route`, then `cargo xtask steps --from-json` |
| [branch-coverage/README.md](branch-coverage/README.md) | Rule 6, sharpened: the minimal finite set of cases that takes every finite-reachable branch of the analysis path (panelling excluded), each compared with its instrumented reference by the studies' whole-run floor comparison with difference metrics; the branch map and the case × subroutine matrix; the minimal non-finite set, whose branches are each tested as an event (`tests/execution/events.rs`). A study: the tests cover branches by single steps instead ([branch-tests.md](branch-tests.md)). `README.tex` in the run folder carries the same tables for the paper | `cargo xtask fixtures --group branch-coverage` and `--group non-finite` (each runs its twins), then `cargo run --release -p branch-coverage -- --docs` |
| [known-issues/README.md](known-issues/README.md) | Runs that show XFOIL's known issues arising (`docs/xfoil-known-issues.md`): the dependence on the panel count (each case at N = 60, 160 and 240), the plot-label hang, and the sweeps whose purpose is one issue; yFoil compared point by point as for the series cases | `cargo xtask fixtures --group known-issues`, then `cargo run --release -p known-issues -- --docs` |
| [noise-floor.md](noise-floor.md) | Rule 1: the reference's own 1-ULP sensitivity per stage and variable, from which `tests/common/utilities/tolerances.rs` is derived | `scripts/noise-floor.sh` |
| [geometry/README.md](geometry/README.md) | Stage G: XFOIL-model NACA4/5 and PANGEN tested bitwise against `pangen_*` dumps; the bitwise LOAD handoff that every solver case asserts | `cargo xtask fixtures --case pangen_…`, `cargo test --test subroutine pangen` |
| [subroutines/README.md](subroutines/README.md) | BL closure relations (HKIN, HSL, HST, CFL, CFT, DIL, DAMPL) against the instrumented reference's own calls, dumped at 17 significant digits, in `tests/fixtures/subroutines/` | `cargo run --bin generate_subroutine_validation` |
| `scripts/input-sensitivity/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the figures (`.fig` and PDF), JSON summaries and a LaTeX index) | The reference's own input sensitivity: NACA 0012 swept 0°–25° with node coordinates (1 ULP … 1e-7), panel count and alpha step perturbed; seven per-alpha quantities per family, base case in front (`scripts/input-sensitivity/src/perturb.rs` for the node-perturbation method) | `cargo run --release -p input-sensitivity` |
| `scripts/xfoil-instrumentation-check/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the 2×2 figure (`.fig` and PDF), `summary.json` and a LaTeX index) | The instrumentation is inert on a polar through stall: pristine and instrumented DP builds run one identical OPER script (XFOIL's own `NACA 0012`, `PPAR N 160`, Re 1e6, 0°–25° by 0.5°, `ALFA` + `DUMP` per point) and are compared byte-for-byte on XFOIL's standard `PACC` polar and `DUMP` files only, at XFOIL's own output precision | `cargo run --release -p xfoil-instrumentation-check` |

**Polar break points (divergent comparisons)** — two tracked cases drive a 0.5° `ASEQ` sweep at `ITER 100`
past CL_max, where the Newton iteration stops converging and both codes wander for the full 105 iterations
per point. The solutions there are ill-conditioned (`docs/conventions/terminology.md`), and how each case is
tested is recorded beside it in `xtask/fixtures-config/cases.toml`: the calls before the break as a whole run
(`tests/execution/runs.rs`), the break call one step at a time from XFOIL's dumped state
(`tests/execution/steps.rs`):

| Case (`tests/fixtures/xfoil/…`) | OPER script | Reference behaviour | Tested as |
|---|---|---|---|
| `naca0012_n160_polar_up22_re1e6_iter100` | `ALFA 0 / ASEQ 0.5 22 0.5` | 20.5°–22° unconverged (on the host-libm reference, before 2026-10-15, 22° converged on the separated branch in 7 iterations) | run through call 41; every iteration of call 42 as a step |
| `naca4412_n160_polar_down16_re1e6_iter100` | `ALFA 0 / INIT / ALFA 0 / ASEQ -0.5 -16 -0.5` | −14.5° converges on the separated branch (CL ≈ +0.10) in 53 iterations, its 1-ULP twins do not (four wander for 105, one converges elsewhere); −15° to −16° follow on that branch (CL ≈ +0.05 at −16°). On the host-libm reference −14.5° to −15.5° were unconverged and −16° converged in 39 | run through call 30; every iteration of call 31 as a step |

In both cases the reference's own 1-ULP twins part from the reference, and flip its stagnation and
transition stations, before yFoil does, and every step replayed from XFOIL's dumped state at the break
reproduces the step. Which attempt converges is therefore decided at the noise floor — changing the
reference's libm moved it in both cases; a polar
that "carries on" past the break in one code and halts in the other is a divergent comparison, not a
translation bug, and no test compares the two codes' own marches through it.

The reference's behaviour at the break does not depend on the host. Both codes take their
transcendentals from CORE-MATH's correctly rounded functions, not the host libm, and the reference
reproduces every tracked record byte for byte on glibc (`scripts/cross-host.sh`), so a test holds every
value to the same tolerance everywhere. Until 2026-10-15 both took them from the host libm, Apple's
libSystem and glibc disagreeing by 1 ULP on 0.1 % (`exp`, `ln`, `pow`) to 18 % (`tanh`) of inputs, and
the break was host-dependent: on glibc the fourth attempt's chance convergence did not survive and XFOIL
failed 22° and −16° as well (`docs/cross-platform-quirks.md`).
The same cases wake the MRCHUE inverse wake march, the MRCHDU iterate fallback, BLVAR's Us and Hk
clamps and TRCHEK2's iteration cap that every other case left open ([coverage.md](coverage.md)). Two
XFOIL quirks surfaced here are registered in `docs/xfoil-known-issues.md` (§2.8, §4): MASS is never
written on side 1's wake slots, and ASEQ's sequence-plot label loops forever on a non-finite CL/CM even
with graphics off, so a reference sweep that blows up past the break hangs instead of halting.
