# Test fixtures: what they are and how to reproduce them

Every tracked file under `tests/fixtures/` is read by a test (`cargo xtask fixtures --audit` fails
otherwise), and every family has one producer. The tests that read each family name it, with its
producer, in their module header (`docs/conventions/testing.md`).

| Family | What it holds | Producer | Needs | Read by |
|---|---|---|---|---|
| `xfoil/<case>/` | One instrumented XFOIL run on yFoil's panels: its inputs (`manifest.json`, `xfoil.inp`, `panels.json`, `panels.dat`) and the dumps its tests read | `cargo xtask fixtures [--case NAME]` from `xtask/fixtures-config/cases.toml` (and the branch cover in `step-cover.toml`) | the instrumented DP reference (`cargo xtask xfoil-build`, built on demand), gfortran | `tests/subroutine/`, `tests/execution/`, `tests/known_issues/`, `tests/apparatus/` |
| `subroutines/<closure>/<case>_NNN.json`, `subroutines/manifest.json` | Distinct calls of the BL closures (HKIN … BLVAR), inputs to outputs | `cargo xtask fixtures --case closures_naca0012_n60_re1e6 --case closures_naca0012_n60_re1e6_m03` (`closures = true`, `xtask/src/closures.rs`) | as above | `tests/subroutine/closures.rs` |
| `naca456/` | The NASA/PDAS `naca456` ordinates for the NACA series (an external reference, no XFOIL) | `scripts/naca456-fixtures.sh` | gfortran, network on the first run (the PDAS source, pinned by SHA-256) | `tests/application/naca456.rs` |
| `repanel_cosine/` | yFoil's own cosine repanelling, frozen as regression snapshots | `UPDATE_SNAPSHOT=1 cargo test --test application repanel_cosine` | nothing | `tests/application/repanel_cosine.rs` |
| `naca0012/panels*.json` | Legacy: panels XFOIL generated itself with PANE (contrary to CLAUDE.md Rule 4) | hand-made, 2026-01; no producer | — | `tests/subroutine/pane_legacy.rs` |
| `naca0012_reference.dat`, `aerofoil_with_invalid_last_point.json` | Hand-written inputs | hand-written | — | `tests/application/geometry.rs`, `geometry_errors.rs` |

## Reproducing the XFOIL families

```bash
cargo xtask xfoil-build --verify          # the reference, with the inertness proof
cargo xtask fixtures                       # every tracked case; --case NAME for one
cargo xtask fixtures --verify              # regenerate and assert byte-identity (same host)
cargo xtask fixtures --audit               # every tracked file is read by some test
```

A case's entry in `cases.toml` states how it is tested (`# how tested:`, `run`, `run_through`,
`step_calls`) and so which files are tracked. The branch cover is recomputed by

```bash
cargo xtask twins --case NAME ...          # the five seeded 1-ULP twins (study data, untracked)
cargo xtask steps                          # which step takes which branch, and its subroutine calls
cargo xtask route                          # yFoil's route through every step against XFOIL's; writes route.toml
cargo xtask steps --from-json              # the cover from the route-agreeing steps; writes step-cover.toml
```

and `cargo xtask fixtures` then adds the dumps the chosen steps need. `cargo xtask route` needs
`rustup component add llvm-tools`; its reference runs, with every SETBL call dumped, live in
`target/route/` and are not tracked. Fixture data is a property of
the host it was generated on (manifest `host:` line); on another host the tests use
`TOL_CROSS_HOST`.

## Study data (untracked)

The studies (`scripts/branch-coverage`, `scripts/series-cases`, `scripts/known-issues`,
`scripts/input-sensitivity`)
read the reference's full runs from `target/fixtures/<case>/`, never from `tests/fixtures/`, and no
test reads what they produce. Their cases carry `group = "..."` in `cases.toml`; the ones only a study
uses are `track = false`.

```bash
cargo xtask fixtures --group series        # the series cases, each with its twins (twins = true)
cargo xtask fixtures --group known-issues  # the known-issues cases (references only)
cargo xtask fixtures --group branch-coverage  # likewise; --group non-finite for the probes
cargo xtask twins --case NAME ...           # rerun a case's five twins and noise_floor.json alone
cargo run --release -p branch-coverage -- --docs
cargo run --release -p series-cases -- --docs
cargo run --release -p known-issues -- --docs
```

Regenerating a case keeps its twins when its panels and script are unchanged.

