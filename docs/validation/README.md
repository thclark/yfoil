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
| `scripts/xfoil-sensitivity/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the SVG sheet, JSON summaries and a LaTeX index) | The reference's own input sensitivity: NACA 0012 swept 0°–25° with node coordinates (1 ULP … 1e-7), panel count and alpha step perturbed; seven per-alpha quantities per family, base case in front (`scripts/xfoil-sensitivity/src/perturb.rs` for the node-perturbation method) | `cargo run --release -p xfoil-sensitivity` |
| `scripts/xfoil-instrumentation-check/` (not a page: each run writes a dated folder under its `runs/`, gitignored, with the 2×2 SVG, `summary.json` and a LaTeX index) | The instrumentation is inert on a polar through stall: pristine and instrumented DP builds run one identical OPER script (XFOIL's own `NACA 0012`, `PPAR N 160`, Re 1e6, 0°–25° by 0.5°, `ALFA` + `DUMP` per point) and are compared byte-for-byte on XFOIL's standard `PACC` polar and `DUMP` files only, at XFOIL's own output precision | `cargo run --release -p xfoil-instrumentation-check` |
| [subroutines/README.md](subroutines/README.md) | BL closure relations (HKIN, HSL, HST, CFL, CFT, DIL, DAMPL) against instrumented E24.16 samples in `tests/fixtures/subroutines/` | `cargo run --bin generate_subroutine_validation` |

**Analysis (polar) validation** — the ±15° N=160 polar with per-alpha CL/CD/CM/XTR and BL-distribution
difference tables against the instrumented per-alpha dumps of the `naca0012_n160_polar_re1e6` case — is
the remaining validation stage and will land here as a generated `analysis/` section. It is deliberately
absent until then: an earlier version of this directory carried XFOIL's formatted `DUMP`/`CPWR`/`PACC`
files (F9.5, F11.5, F10.5 — about 1e-5) as committed reference data, which cannot support the project's
tolerances (CLAUDE.md Rule 4), and was removed from history on 2026-09-04.

**Polar break points (threshold-straddling)** — two tracked cases drive a 0.5° `ASEQ` sweep at `ITER 100`
past CL_max, where the Newton iteration stops converging and both codes wander for the full 105 iterations
per point. They are the worked examples of Rule 1's third outcome for a multi-point sequence, gated by
`tests/xfoil_polar_break_tests.rs`:

| Case (`tests/fixtures/xfoil/…`) | OPER script | Reference behaviour | Outcome |
|---|---|---|---|
| `naca0012_n160_polar_up22_re1e6_iter100` | `ALFA 0 / ASEQ 0.5 22 0.5` | 20.5°–21.5° unconverged; 22° converges on the separated branch (CL ≈ 0.35) in 7 iterations, its +1-ULP twin does not | calls 1–41 match; call 42 straddles at iteration 32; yFoil fails 22° and the sweep halts (NSEQEX = 4) |
| `naca4412_n160_polar_down16_re1e6_iter100` | `ALFA 0 / INIT / ALFA 0 / ASEQ -0.5 -16 -0.5` | −14.5°–−15.5° unconverged; −16° converges on the separated branch (CL ≈ +0.05) in 39 iterations, its twin does not | calls 1–30 match; call 31 straddles at iteration 16 (IST/ITRAN); yFoil converges −16° in 17 iterations on the same branch |

In both cases the reference's own +1-ULP twin parts from the reference, and flips its stagnation and
transition stations, before yFoil does, and the one-step replays from XFOIL's dumped state at the parting
iterations (`mrchdu_input_<k>.dat` → `update_output_<k>.dat`) reproduce the step to 2e-12 in every array.
Whether the fourth attempt converges is therefore decided at the noise floor; a polar that "carries on"
past the break in one code and halts in the other is not a gate difference and not a translation bug.
The same cases wake the MRCHUE inverse wake march, the MRCHDU extrapolation fallback, BLVAR's Us and Hk
clamps and TRCHEK2's iteration cap that every other case left open ([coverage.md](coverage.md)). Two
XFOIL quirks surfaced here are registered in `docs/xfoil-known-issues.md` (§2.8, §4): MASS is never
written on side 1's wake slots, and ASEQ's sequence-plot label loops forever on a non-finite CL/CM even
with graphics off, so a reference sweep that blows up past the break hangs instead of halting.

The gates themselves are the tests (`tests/xfoil_*_tests.rs`); these pages tabulate and plot what the
tests assert, they do not add evidence of their own.
