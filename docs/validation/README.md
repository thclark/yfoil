# Validation

Every page in this directory is generated from data the fixture pipeline produced, and can be
regenerated at any time from tracked inputs. Nothing here is a baseline: the reference is the
instrumented double-precision XFOIL 6.99 build (CLAUDE.md, *Reference implementation*), its dumps
live under `tests/fixtures/xfoil/<case>/` (tracked) or `target/fixtures/<case>/` (`--big` cases),
and the reports read those directories. Only Markdown and SVG are committed under `docs/validation/`;
`.gitignore` refuses everything else.

Figures are drawn by matplotlib from each study's JSON outputs (`scripts/<study>/plot.py`, presentation only;
`scripts/figures/style.py` holds the publication conventions, `scripts/figures/render.sh` pins the version and runs it
through uv, so running any of the studies needs [`uv`](https://docs.astral.sh/uv/) installed). Every number in a figure, including axis extents, is computed by the Rust study and is in its run folder,
so the plots can be redrawn by anyone from the JSON alone.

| Page | What it shows | Regenerate with |
|---|---|---|
| [coverage.md](coverage.md) | Rule 6: gcov branch completeness of the reference over every tracked case, with the open/unreachable/loop-entry classification | `cargo xtask coverage` |
| [branch-gating.md](branch-gating.md) | How full branch coverage is gated: every candidate case and how it is tested (run or single steps), which step takes which branch, the fewest steps that take them all (each replayed by a test), and the branches no test gates | `cargo xtask steps`, `cargo xtask route`, then `cargo xtask steps --from-json` |
| [branch-coverage/README.md](branch-coverage/README.md) | Rule 6, sharpened: the minimal set of cases that takes every finite-reachable branch of the analysis path (panelling excluded), each compared with its instrumented reference by the studies' whole-run floor comparison with difference metrics; the branch map and the case × subroutine matrix; the branches reached only in non-finite reference runs, each gated as an event (`tests/execution/events.rs`). A study: the tests gate branches by single steps instead ([branch-gating.md](branch-gating.md)). `README.tex` in the run folder carries the same tables for the paper | `cargo xtask fixtures --group branch-coverage` and `--group non-finite` (each runs its twins), then `cargo run --release -p branch-coverage -- --docs` |
| [series-cases/README.md](series-cases/README.md) | Runs of the aerofoil series and conditions of interest — full ±30° polars by XFOIL's polar procedure, inviscid alpha sweeps and fixed-CL sweeps — with yFoil driven the same way and compared point by point through stall and into the non-finite region; each run with its own five 1-ULP twins, so the page also says where XFOIL itself was ill-conditioned | `cargo xtask fixtures --group series` (runs the twins too), then `cargo run --release -p series-cases -- --docs` |
| [noise-floor.md](noise-floor.md) | Rule 1: the reference's own 1-ULP sensitivity per stage and variable, from which `tests/common/utilities/tolerances.rs` is derived | `scripts/noise-floor.sh` |
| [geometry/README.md](geometry/README.md) | Stage G: XFOIL-model NACA4/5 and PANGEN gated bitwise against `pangen_*` dumps; the bitwise LOAD handoff that every solver case asserts | `cargo xtask fixtures --case pangen_…`, `cargo test --test subroutine pangen` |
| `scripts/input-sensitivity/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the SVG sheet, JSON summaries and a LaTeX index) | The reference's own input sensitivity: NACA 0012 swept 0°–25° with node coordinates (1 ULP … 1e-7), panel count and alpha step perturbed; seven per-alpha quantities per family, base case in front (`scripts/input-sensitivity/src/perturb.rs` for the node-perturbation method) | `cargo run --release -p input-sensitivity` |
| `scripts/xfoil-instrumentation-check/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the 2×2 SVG, `summary.json` and a LaTeX index) | The instrumentation is inert on a polar through stall: pristine and instrumented DP builds run one identical OPER script (XFOIL's own `NACA 0012`, `PPAR N 160`, Re 1e6, 0°–25° by 0.5°, `ALFA` + `DUMP` per point) and are compared byte-for-byte on XFOIL's standard `PACC` polar and `DUMP` files only, at XFOIL's own output precision | `cargo run --release -p xfoil-instrumentation-check` |
| [subroutines/README.md](subroutines/README.md) | BL closure relations (HKIN, HSL, HST, CFL, CFT, DIL, DAMPL) against instrumented E24.16 samples in `tests/fixtures/subroutines/` | `cargo run --bin generate_subroutine_validation` |
| [aerofoil-series/README.md](aerofoil-series/README.md) | The selection of aerofoils the equivalence cases run on (one per family, with why), every generator family drawn one panel per parameter, and the NACA generators against the public-domain NASA/PDAS `naca456` ordinates (`tests/fixtures/naca456/`, written by `scripts/naca456-fixtures.sh`; gated by `tests/application/naca456.rs`) | `cargo run --release -p aerofoil-series -- --docs` |

**Analysis (polar) validation** — the ±15° N=160 polar with per-alpha CL/CD/CM/XTR and BL-distribution
difference tables against the instrumented per-alpha dumps of the `naca0012_n160_polar_re1e6` case — is
the remaining validation stage and will land here as a generated `analysis/` section. It is deliberately
absent until then: an earlier version of this directory carried XFOIL's formatted `DUMP`/`CPWR`/`PACC`
files (F9.5, F11.5, F10.5 — about 1e-5) as committed reference data, which cannot support the project's
tolerances (CLAUDE.md Rule 4), and was removed from history on 2026-09-04.

**Polar break points (divergent comparisons)** — two tracked cases drive a 0.5° `ASEQ` sweep at `ITER 100`
past CL_max, where the Newton iteration stops converging and both codes wander for the full 105 iterations
per point. The solutions there are ill-conditioned (`docs/conventions/terminology.md`), and how each case is
tested is recorded beside it in `xtask/fixtures-config/cases.toml`: the calls before the break as a whole run
(`tests/execution/runs.rs`), the break call one step at a time from XFOIL's dumped state
(`tests/execution/steps.rs`):

| Case (`tests/fixtures/xfoil/…`) | OPER script | Reference behaviour | Tested as |
|---|---|---|---|
| `naca0012_n160_polar_up22_re1e6_iter100` | `ALFA 0 / ASEQ 0.5 22 0.5` | 20.5°–21.5° unconverged; 22° converges on the separated branch (CL ≈ 0.35) in 7 iterations, its 1-ULP twins do not | run through call 41; every iteration of call 42 as a step |
| `naca4412_n160_polar_down16_re1e6_iter100` | `ALFA 0 / INIT / ALFA 0 / ASEQ -0.5 -16 -0.5` | −14.5°–−15.5° unconverged; −16° converges on the separated branch (CL ≈ +0.05) in 39 iterations, its twins do not | run through call 30; every iteration of call 31 as a step |

In both cases the reference's own 1-ULP twins part from the reference, and flip its stagnation and
transition stations, before yFoil does, and every step replayed from XFOIL's dumped state at the break
reproduces the step. Whether the fourth attempt converges is therefore decided at the noise floor; a polar
that "carries on" past the break in one code and halts in the other is a divergent comparison, not a
translation bug, and no test compares the two codes' own marches through it.

The reference's behaviour holds on the host that generated the fixtures (arm64 macOS; each `manifest.json` records
it). Both codes take `EXP`/`LOG`/`**`/`ATAN2`/`SIN`/`COS`/`TANH` from the host C library, and Apple's
libSystem and glibc disagree by exactly 1 ULP on 0.1 % (`exp`, `ln`, `pow`) to 18 % (`tanh`) of inputs
(measured 2026-09-11, 20 000 inputs per function; glibc's results are the same on x86_64 and aarch64). The
same two cases regenerated on glibc show what that does to the reference itself: every converged point
before the break moves by ≤ 1.4e-11 in CL, the unconverged wanderings by O(1), and the fourth attempt's
chance convergence does not survive — on glibc XFOIL fails 22° and −16° as well (its converged points
agree with the macOS run to 1.4e-11 and 7.5e-11 in CL). Run against the macOS fixtures on
glibc, yFoil matches every converged call at the same tolerances, parts inside the break call at a
different iteration, and its one-step replays from XFOIL's dumped state
differ by 3.6e-10 and 1.2e-10 relative in the largest Newton delta where the same-host replays differ by
≤ 2e-12 — the libm's 1-ULP differences through one near-singular Newton step, gated by `TOL_CROSS_HOST`
(`tests/common/utilities/tolerances.rs`); `tests/common/utilities/host.rs` holds a test to `TOL_SOLVER` on
the fixture's host and to `TOL_CROSS_HOST` on any other.
The same cases wake the MRCHUE inverse wake march, the MRCHDU extrapolation fallback, BLVAR's Us and Hk
clamps and TRCHEK2's iteration cap that every other case left open ([coverage.md](coverage.md)). Two
XFOIL quirks surfaced here are registered in `docs/xfoil-known-issues.md` (§2.8, §4): MASS is never
written on side 1's wake slots, and ASEQ's sequence-plot label loops forever on a non-finite CL/CM even
with graphics off, so a reference sweep that blows up past the break hangs instead of halting.

The gates themselves are the tests (`docs/conventions/testing.md`); the study pages above measure the
reference and compare whole sweeps, and add no evidence the tests rely on.
