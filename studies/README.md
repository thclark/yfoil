# Studies

A study demonstrates: it shows an outside reader how a condition, regime or behaviour arises, and
whether yFoil follows the reference through it, and writes its findings as a report. A test checks:
it asserts correct, consistent functionality as minimally and quickly as it can, and lives in
`tests/`. One behaviour may need both (`docs/development/conventions/testing.md`, "Studies demonstrate; tests
check"). No test reads a study's output.

Each study writes a dated run folder under its own `runs/` (gitignored), and never writes under
`docs/`. A run folder holds `publish/`, the run's tables and figures laid out as the docs hold them
(`study-support/publish.rs`); `cargo xtask publish <study> [--run <UTC>]` copies one run's bundle into
`docs/validation/<study>/`, whose pages are written by hand
(`docs/development/conventions/documentation.md`). Run every command from the repository root.

| Study | What it shows | Run with | Published to |
|---|---|---|---|
| `aerofoil-series/` | The selection of aerofoils the series and test cases run on, every generator family drawn one panel per parameter, and the NACA generators against the NASA/PDAS `naca456` ordinates | `cargo run --release -p aerofoil-series` | `docs/validation/aerofoil-series/` (`cargo xtask publish aerofoil-series`) |
| `series-cases/` | The aerofoil series and flow regimes of interest under normal conditions (±30° polars, an inviscid sweep, a fixed-CL sweep at N = 160, ITER 100), yFoil against the reference point by point, with each run's 1-ULP twins showing where XFOIL itself was ill-conditioned | `cargo xtask fixtures --group series`, then `cargo run --release -p series-cases` | `docs/validation/series-cases/` (`cargo xtask publish series-cases`) |
| `branch-coverage/` | The minimal set of cases that takes every finite-reachable branch of the analysis path, each compared with its reference run; the branch map and case × subroutine matrix; the minimal non-finite set | `cargo xtask fixtures --group branch-coverage` and `--group non-finite`, then `cargo run --release -p branch-coverage` | `docs/validation/branch-coverage/` (`cargo xtask publish branch-coverage`) |
| `known-issues/` | Runs that show XFOIL's known issues arising: the dependence on the panel count (N = 60, 160, 240), the plot-label hang, and the sweeps whose purpose is one issue | `cargo xtask fixtures --group known-issues`, then `cargo run --release -p known-issues` | `docs/validation/known-issues/` (`cargo xtask publish known-issues`) |
| `tolerances/` | The reference's noise floor under a +1 ULP perturbation of the reference case, per stage and variable; the constants in `tests/common/utilities/tolerances.rs` were set from its findings, by hand | `studies/tolerances/tolerances.sh` | `docs/validation/tolerances/` (`cargo xtask publish tolerances`) |
| `input-sensitivity/` | The reference's own input sensitivity: NACA 0012 swept 0°–25° with node coordinates, panel count and alpha step perturbed | `cargo run --release -p input-sensitivity` | its run folder only |
| `xfoil-instrumentation-check/` | The instrumentation is inert on a polar through stall: pristine and instrumented builds compared byte for byte on XFOIL's own `PACC` and `DUMP` output | `cargo run --release -p xfoil-instrumentation-check` | its run folder only |

Shared support, not studies: `study-support/records.rs` is the studies' whole-run floor comparison
(used by `branch-coverage` and `series-cases`), `study-support/publish.rs` writes a run's `publish/`, and `figures/` (`render.sh`, `style.py`) draws the
branch-coverage study's two matplotlib figures (needs [`uv`](https://docs.astral.sh/uv/)).
