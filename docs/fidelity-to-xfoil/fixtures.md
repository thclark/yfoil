---
icon: lucide/database
---

# Fixtures

A fixture is one captured run of the [reference build](the-reference-xfoil-build/index.md): its
inputs and the full-precision dumps a test reads. Every equivalence test on this tab compares yFoil
with a fixture, never with a number typed into the source, and never with XFOIL's formatted
output. This page says how cases are defined, how both codes are made to run on the same panels,
what each fixture records about its own making, and which fixtures are committed. The runbook for
every fixture family, with the commands that reproduce it, is
[`tests/fixtures/README.md`](https://github.com/thclark/yfoil/blob/main/tests/fixtures/README.md).

## Case definitions and groups

`cargo xtask fixtures [--case NAME] [--group NAME] [--verify] [--big]` reads
`xtask/fixtures-config/cases.toml`, builds the reference if needed, generates the panels **with
yFoil**, runs the instrumented XFOIL on them, asserts the bitwise geometry handoff, and keeps the
case under `tests/fixtures/xfoil/<case>/` (`track = true`) or `target/fixtures/<case>/`.

A case's entry says what to run:

| Option | Meaning |
| --- | --- |
| `alphas`, `alphas_after_reinit` | OPER `ALFA` points, with an `INIT` between the two lists |
| `polar = true` | drive with `ALFA a0 / ASEQ a1 aN da`, as the polar procedure does; `ASEQ` points get ITMAX + 5 iterations |
| `cls` | OPER `CL` points |
| `matyp` | OPER `TYPE n` |
| `inviscid = true` | no `VISC`: each `ALFA` is one SPECAL, each `CL` one SPECCL |
| `events = true` | keep `events.dat`, the reference's branch events |
| `keep = [...]` | track these work-directory files besides the inputs |
| `closures = true` | cut the reference's closure log into `tests/fixtures/subroutines/` |
| `twins = true` | also run the case's five seeded 1-ULP twins and write its `noise_floor.json` (study data, under `target/fixtures/<case>/ulp<seed>/`) |
| `stop_on_nonfinite`, `stop_after_setbl = N` | stop the reference on the first SETBL call whose state is non-finite, or after N calls; never allowed on a case that tests or measures branches |

and how it is tested — `run`, `run_through`, `step_calls`, `dump_calls` — which is
[execution equivalence](execution-equivalence/index.md#how-each-case-is-tested). `group = "name"`
names the study a case belongs to: `branch-coverage` (the minimal set), `non-finite` (its NaN
probes), `series` (the series cases), `known-issues` (the runs that show XFOIL's known issues) and
`pathological` (copies of cases whose runs are pathologically slow or hang, run only by name).

Every reference run, twins included, is watched: killed after 120 s without output (XFOIL's
plot-label hang on a non-finite value, [known issues §4](xfoils-known-issues/the-issue-register.md))
and kept as truncated. Where a run ends is never decided by the clock.

## The bitwise geometry handoff

Different panels make any comparison meaningless, so **yFoil generates, XFOIL consumes**:

``` sh
yfoil geometry naca 0012 -n 160 --method cosine -o geometry.json
yfoil geometry convert geometry.json --to dat -o geometry.dat
# XFOIL:  PLOP / G F / LOAD geometry.dat / ...      never NACA, never PANE, never PPAR
```

This works because `LOAD` does not transform coordinates (`LNORM` is false by default, and
`ABCOPY` copies buffer to current with no repanelling). The fixture pipeline **asserts** it on every
case of every run: the `.dat` is written at 17 significant digits, XFOIL dumps `X/Y` at `ES24.16`
immediately after `ABCOPY`, and the two must be bitwise equal before any fixture from that run is
accepted. It has held on every case since 2026-09-03. `tests/apparatus/inputs.rs` checks each
tracked case's inputs against each other, and `tests/apparatus/tgap_handoff.rs` checks the same
handoff from inside XFOIL's `TGAP`.

One consequence: yFoil's panelling only has to produce *good* panels, not XFOIL-identical ones.
That yFoil can also reproduce XFOIL's own `NACA` and `PANGEN` bitwise is an optional check,
[XFOIL's own geometry routines](subroutine-equivalence/xfoils-own-geometry-routines.md). And
XFOIL's `LOAD` reverses a clockwise input, so yFoil always writes counter-clockwise
([known issues §6.2](xfoils-known-issues/the-issue-register.md)).

## Manifests and fingerprints

Every fixture directory carries a `manifest.json`: the case definition, the SHA-256 of the
sources, the patch series' SHA, `gfortran --version`, the FFLAGS, the host triple and the maths
library, and `inputs_fingerprint` — the fingerprint of the reference stage that produced it
([the reference build](the-reference-xfoil-build/index.md#layout-and-how-it-is-built)). A fixture
whose fingerprint is not the current build's is regenerated.

## Tracked and regenerated fixtures

`.gitignore` excludes generated files by default (`*.dat`, `*.log`, `fort.*`, `*.bl`, `.tmp/`)
and **allowlists** what is deliberately committed (`tests/fixtures/**`, `xfoil/third-party/**`):
committing a fixture is a decision.

- A **tracked** case is its inputs plus exactly the files its tests read (`keep`, and what `run`,
  `step_calls` and the branch cover add). `cargo xtask fixtures --audit` fails on any tracked file
  no test opens.
- A case no test reads is **untracked** (`track = false`) and regenerated into
  `target/fixtures/<case>/`; the studies read only those, and no test reads what a study produces.
- `cargo xtask fixtures --verify` regenerates every tracked case and asserts it byte-identical to
  what is committed, on any host (the manifest, which names the host, aside).
- A missing fixture **fails** its test (`require_fixture`); it is never skipped
  ([the test suite](../contributing/the-test-suite.md)).

The tracked CI reference case is `naca0012_n60_a2_re1e6`.
