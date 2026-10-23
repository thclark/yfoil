# Validation

The engineering assessment: whether yFoil stacks up sensibly against XFOIL on the aerofoil series and
flow regimes an aerodynamicist would actually run. The results read as polars, transition locations and
boundary layers, and what they mean is plain to anyone who uses XFOIL. How the test suite tracks yFoil
against XFOIL statement by statement — branch coverage, the minimal finite and non-finite case sets,
the known-issues and pathological runs — is for software developers, and
lives in [Testing](../testing/README.md).

Every page here is generated from data the fixture pipeline produced, and can be regenerated at any time
from tracked inputs. Nothing here is a baseline: the reference is the instrumented double-precision XFOIL
6.99 build (CLAUDE.md, *Reference implementation*), its dumps live under `tests/fixtures/xfoil/<case>/`
(tracked) or `target/fixtures/<case>/` (`--big` cases), and the reports read those directories. Only
Markdown, ironlab figures (`.fig`) and the branch-coverage SVGs are committed under `docs/`.

Figures are [ironlab](https://github.com/thclark/ironlab) figures, built in Rust by each study from the numbers it
computes (`docs/conventions/figures.md`). They are embedded live in these pages: pan, zoom and toggle lines from the
legend. Each carries parameters (case, conditions, outputs) and labels, so every study's figures can be loaded
together into the ironlab viewer and sorted, filtered and grouped. A study's run folder also holds a PDF of each.

| Page | What it shows | Regenerate with |
|---|---|---|
| [Series cases](series-cases.md) | Runs of the aerofoil series and flow regimes of interest under normal conditions — full ±30° polars by XFOIL's polar procedure, an inviscid alpha sweep and a fixed-CL sweep, at N = 160, ITER 100 — with yFoil driven the same way and compared point by point; each run with its own five 1-ULP twins, so the page also says where XFOIL itself was ill-conditioned | `cargo xtask fixtures --group series` (runs the twins too), then `cargo run --release -p series-cases -- --docs` |
| [Aerofoil series](aerofoil-series.md) | The selection of aerofoils the series and the test cases run on (one per family, with why), every generator family drawn one panel per parameter, and the NACA generators against the public-domain NASA/PDAS `naca456` ordinates (`tests/fixtures/naca456/`, written by `scripts/naca456-fixtures.sh`; tested by `tests/application/naca456.rs`) | `cargo run --release -p aerofoil-series -- --docs` |
| [Tolerances](tolerances.md) | Rule 1: the reference's own noise floor — its 1-ULP sensitivity per stage and variable — from which the constants in `tests/common/utilities/tolerances.rs` were set by hand. The tables are fragments the study writes into `tables/`, included into the page | `studies/tolerances/tolerances.sh` |

**Analysis (polar) validation** — the ±15° N=160 polar with per-alpha CL/CD/CM/XTR and BL-distribution
difference tables against the instrumented per-alpha dumps of the `naca0012_n160_polar_re1e6` case — is
the remaining validation stage and will land here as a generated `analysis/` section. It is deliberately
absent until then: an earlier version of this directory carried XFOIL's formatted `DUMP`/`CPWR`/`PACC`
files (F9.5, F11.5, F10.5 — about 1e-5) as committed reference data, which cannot support the project's
tolerances (CLAUDE.md Rule 4), and was removed from history on 2026-09-04.
