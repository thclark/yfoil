# Test fixtures: what they are and how to reproduce them

Every tracked file under `tests/fixtures/` is read by a test (`cargo xtask fixtures --audit` fails
otherwise), and every family has one producer. The tests that read each family name it, with its
producer, in their module header (`docs/development/conventions/testing.md`).

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
`target/route/` and are not tracked. Fixture data does not depend on the host it was generated
on: both codes take their maths from CORE-MATH's correctly rounded functions, not the host libm, so
the tests hold every value to the same tolerance on every host. That the reference reproduces the
fixtures bit for bit elsewhere is checked by

```bash
scripts/cross-host.sh    # needs docker: reruns the reference on glibc on every tracked run case's
                         # own panels.dat and xfoil.inp; fails unless every record is byte-identical
```

## Study data (untracked)

The studies (`studies/branch-coverage`, `studies/series-cases`, `studies/known-issues`,
`studies/input-sensitivity`)
read the reference's full runs from `target/fixtures/<case>/`, never from `tests/fixtures/`, and no
test reads what they produce. Their cases carry `group = "..."` in `cases.toml`; the ones only a study
uses are `track = false`.

```bash
cargo xtask fixtures --group series        # the series cases, each with its twins (twins = true)
cargo xtask fixtures --group known-issues  # the known-issues cases (references only)
cargo xtask fixtures --group branch-coverage  # likewise; --group non-finite for the probes
cargo xtask twins --case NAME ...           # rerun a case's five twins and noise_floor.json alone
cargo run --release -p branch-coverage      # each run writes studies/<study>/runs/<UTC>/
cargo run --release -p series-cases
cargo run --release -p known-issues
cargo xtask publish series-cases            # copy a run's tables and figures into docs/validation/
```

Regenerating a case keeps its twins when its panels and script are unchanged.

## Pathological runs on the host libm (untracked)

Before 2026-10-15 the reference took its maths from the host libm, and four known-issues sweeps
hung in XFOIL's plot label (known issues §4). On CORE-MATH they do not hang. Their
`pathological_…_host_libm` copies in `cases.toml` recreate the hangs. Each copy carries
`reference = "host-libm"`, which runs it on the instrumented build linked against the host libm.
The 64A010 copy also carries `panels = "xtask/fixtures-config/panels/…_b448f8d4.json"`, the panels
yFoil generated before its own generators moved to CORE-MATH; with today's panels the 64A010 sweep
does not hang on either reference.

```bash
cargo xtask xfoil-build --host-libm        # target/xfoil-ref/host-libm/ (built on demand too)
cargo xtask fixtures --case pathological_naca16-212_n160_polar30_re1e6_m07_iter100_host_libm
cargo xtask fixtures --group pathological  # every pathological copy, the host-libm ones included
```

A run that hangs is ended by the watchdog and its manifest carries `truncated_by_watchdog`. The
copies' descriptions give where each one hung (alpha, leg, VISCAL and SETBL call). The host-libm
build is not a reference and no test or study reads its runs. They depend on the host: the
descriptions were recorded on macOS 15.6.1 (arm64, libSystem), and
`target/xfoil-ref/host-libm/manifest.txt` names the libm a rerun used. On glibc the same run may
take another route. `cargo xtask fixtures` refuses `reference` and `panels` on any case that is not
an untracked `pathological` copy.

