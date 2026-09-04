# Validation

Every page in this directory is generated from data the fixture pipeline produced, and can be
regenerated at any time from tracked inputs. Nothing here is a baseline: the reference is the
instrumented double-precision XFOIL 6.99 build (CLAUDE.md, *Reference implementation*), its dumps
live under `tests/fixtures/xfoil/<case>/` (tracked) or `target/fixtures/<case>/` (`--big` cases),
and the reports read those directories. Only Markdown and SVG are committed under `docs/validation/`;
`.gitignore` refuses everything else.

| Page | What it shows | Regenerate with |
|---|---|---|
| [coverage.md](coverage.md) | Rule 6: gcov branch completeness of the reference over every tracked case, with the open/unreachable/loop-entry classification | `cargo xtask coverage` |
| [noise-floor.md](noise-floor.md) | Rule 1: the reference's own 1-ULP sensitivity per stage and variable, from which `tests/utilities/tolerances.rs` is derived | `scripts/noise-floor.sh` |
| [geometry/README.md](geometry/README.md) | Stage G: XFOIL-model NACA4/5 and PANGEN gated bitwise against `pangen_*` dumps; the bitwise LOAD handoff that every solver case asserts | `cargo xtask fixtures --case pangen_…`, `cargo test --test xfoil_pangen_tests` |
| [subroutines/README.md](subroutines/README.md) | BL closure relations (HKIN, HSL, HST, CFL, CFT, DIL, DAMPL) against instrumented E24.16 samples in `tests/fixtures/subroutines/` | `cargo run --bin generate_subroutine_validation` |

**Analysis (polar) validation** — the ±15° N=160 polar with per-alpha CL/CD/CM/XTR and BL-distribution
difference tables against the instrumented per-alpha dumps of the `naca0012_n160_polar_re1e6` case — is
the remaining validation stage and will land here as a generated `analysis/` section. It is deliberately
absent until then: an earlier version of this directory carried XFOIL's formatted `DUMP`/`CPWR`/`PACC`
files (F9.5, F11.5, F10.5 — about 1e-5) as committed reference data, which cannot support the project's
tolerances (CLAUDE.md Rule 4), and was removed from history on 2026-09-04.

The gates themselves are the tests (`tests/xfoil_*_tests.rs`); these pages tabulate and plot what the
tests assert, they do not add evidence of their own.
