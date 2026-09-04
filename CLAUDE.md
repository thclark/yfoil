# YFoil — XFOIL in Rust

## Committing

When making commits in git, NEVER attribute Claude (yourself) as a contributor. Reasons: 
1. A contributor is a human who takes responsibility for the code; LLMs (you) cannot do this. 
2. We (the Open-Source community) built the code that was used to train Claude (you), and never got any credit or compensation for that. Not attributing Anthropic/Claude to our outputs built with this tool is consistent with Anthropic's own practice.
3. Years of Mark Drela's work went into building XFOIL; to claim credit for somthing that's bit-for-bit translated is crass at best.

## STOP. READ THIS FIRST.

YFoil is a **translation** of XFOIL 6.99 into Rust. Not a reimplementation, not "something like XFOIL". The
rules below are non-negotiable; the nuance in Rule 1 exists to make them *provable*, not to soften them.

### Rule 1: Numerical equivalence — defined, measured, then enforced

**The reference** is XFOIL 6.99 compiled in **double precision** (`-fdefault-real-8`, `libPlt_gDP.a`) with
pinned flags (`-ffp-contract=off`, one fixed optimisation level) on a recorded toolchain. Stock single-precision
XFOIL agrees with this build only to ~1e-7; when this document says "XFOIL" it means the DP reference build.

**Tolerances are derived, not asserted.** The reference's own noise floor is measured by `scripts/noise-floor.sh`
(1-ULP input perturbation, per stage, per variable; results in `docs/validation/noise-floor.md`) and every
tolerance is that floor times a safety factor. Measured 2026-09-03: forces/RMSBL move ≤ 2.2e-10 per iteration,
BL state ≤ 2e-11, inviscid ≤ 6e-11 — so 1e-10 on forces is *at* the floor, and matrix entries (DIJ, SETBL
Jacobian) can only be gated with a row-scaled metric. Four named constants live in `tests/utilities/tolerances.rs` (`TOL_PURE`, `TOL_LINALG`, `TOL_SOLVER`, and
`TOL_TRANSIENT` for per-iteration transients inside multi-point sequences — measured 2026-09-04, floor 5.1e-10);
no ad-hoc literals anywhere else. Expect pure closure functions at ~1e-14, linear solves at ~1e-11, converged Newton
state at ~1e-10. Transcendentals (`**`, `EXP`, `LOG`, `ATAN2`) come from the host libm in *both* codes: bit-identity
is a same-host property, cross-host is an ULP budget.

**The error metric** is `|a − b| ≤ tol · max(|a|, |b|, scale_v)` with a physical per-variable scale. Bare relative
error is undefined at CL≈0, VDEL≈0 and laminar CTAU≈0 and must not be used.

**A comparison passes only when the branch trace is identical AND values are within tolerance.** XFOIL contains
~130 exact real equalities and ~140 real inequality branches (RLX clamps, BLSOLV sparsity skips, `AMPL>ACRIT`,
`SHARP`, `RMSBL<1e-4`, the MRCHUE `DIRECT` switch, Hk clamps…). A 1-ULP difference at a threshold flips a branch
and produces O(1) differences that are **not translation bugs**. So:

- Instrumented XFOIL logs its decisions (`ITRAN`, limiting station, skip counts, iteration count, `SHARP`,
  `DIRECT`). Those are compared exactly.
- For iterative solvers the iteration count must be identical and per-iteration `RMSBL` is compared.
- When a branch flips, show both inputs lie within the noise floor of the threshold and record the case as
  **threshold-straddling** — a third outcome, reported separately, never silently passed or failed.
- **Mechanised per case:** `cargo xtask fixtures` runs every case twice — as generated and with every panel
  coordinate +1 ULP — and writes `noise_floor.json` (the reference's own spread of every recorded value, and
  whether its branch trace survived). Tests gate a value at `max(tol · scale, FLOOR_FACTOR · floor)`
  (`tests/utilities/records.rs`); a twin that flips its own branch trace, or a run that matches every
  iteration until one where the reference moves by more than `STRADDLE_FLOOR` under 1 ULP, is classified
  threshold-straddling — and the one-step replay from XFOIL's dumped state at that iteration (`dump_calls`)
  is the evidence that the step itself is faithful. The 12° NACA 0012 case is the worked example.

Under those conditions: if values differ by more than tolerance, **it is a bug**. 1%, 0.1%, 1e-6 — all bugs. There
is no "acceptable engineering tolerance" in this project.

### Rule 2: NEVER modify the algorithm

When something doesn't work the temptation is to "fix" it. **This is always wrong.** Do not:

- add relaxation factors, iteration limits, convergence bail-outs, clamps, limits or safeguards not in XFOIL;
- use "simplified" versions of XFOIL's algorithms or approximate formulas;
- change signs or coefficients to "make it work";
- disable a coupling term "temporarily".

Any such change masks the real bug and creates new ones. CI runs a `no-deviations` grep gate for the tell-tale
vocabulary (`disable_`, `TEMPORARILY`, `simplified`, `causes divergence`, `coupling_scale`, bare `break` inside a
Newton loop). If XFOIL does X, YFoil does X — including XFOIL's own bugs and quirks.

**Deliberate divergences from XFOIL are allowed only when documented here.** Currently:

| Divergence | Why | Where |
|---|---|---|
| NACA 4/5-digit thickness applied perpendicular to the camber line | XFOIL's `NACA4` (`naca.f:62`) applies it vertically, which is not the NACA definition. Irrelevant for validation because YFoil generates the panels and XFOIL consumes them (Rule 4). XFOIL's variant is available as `--naca-model xfoil`. | `src/geometry/naca.rs` |

Adding a row to that table requires the same evidence standard as a bug report against XFOIL.

### Rule 3: Debug by forward-stepping through XFOIL's execution

VISCAL is a nonlinear iterative solver. **You cannot debug it from the end.** Small early differences compound;
final CD/CL tell you nothing about where the bug is. Fixing one bug may make the end result *worse* — multiple bugs
cancel.

The only valid approach:

1. Instrument XFOIL at the exact point under investigation (tracked patch series, see *Reference implementation*).
2. Rebuild the reference.
3. Run both codes on identical inputs at full precision (`E24.16` on the Fortran side, `{:.17e}` on the Rust side).
4. Compare at that point. Values outside tolerance with an identical branch trace = bug. Fix it.
5. Advance to the next computation **in XFOIL's execution order** (not panel order, not module order).

Prefer the **replay harness**: seed YFoil with XFOIL's complete state at iteration *k* and require XFOIL's state at
*k+1*. That isolates one iteration of one subroutine from everything upstream. For a single subroutine, capture its
inputs and outputs and test it in isolation.

Progress is measured by fixture matches, never by an end-to-end number moving in the right direction.

### Rule 4: Identical geometry — YFoil generates, XFOIL consumes

Different panels make any comparison meaningless. The workflow, and the only one:

```bash
yfoil geometry naca 0012 -n 160 -o geometry.json      # or repanel: yfoil geometry repanel --method xfoil
yfoil geometry convert geometry.json --to dat -o geometry.dat
# XFOIL:  PLOP / G F / LOAD geometry.dat / ...      never NACA, never PANE, never PPAR
```

This works because `LOAD` does not transform coordinates (`LNORM=.FALSE.` by default; `ABCOPY` copies
buffer→current with no repaneling). The fixture pipeline **asserts** it: `.dat` is written at 17 significant
digits, XFOIL dumps `X/Y` at `E24.16` immediately after `ABCOPY`, and the two must be bitwise equal before any
fixture from that run is accepted.

Consequences:
- YFoil's PANE port only has to produce *good* panels, not XFOIL-identical ones. XFOIL-fidelity of `NACA4`/`PANGEN`
  is an optional quality goal, off the critical path.
- **Never gate on XFOIL's formatted output files.** `PSAVE` is `G15.7`, `.pol` is `F9.4/F10.5`, `DUMP` is
  `F9.5/F10.6`, `CPWR` is `F11.5`. They cannot support anything below ~1e-5. Everything that gates comes from
  instrumented `E24.16` dumps.

### Rule 5: What to do when something is wrong

1. STOP. Do not tweak parameters. Do not add workarounds.
2. Isolate the smallest unit that could be wrong — a subroutine, a statement, one operation.
3. Instrument XFOIL there; rebuild; capture a fixture.
4. Write the test against the fixture.
5. Fix YFoil to match. Check the branch trace matches too.
6. Advance.

If you are looking at final CD/CL to understand a bug, you are doing it wrong. Go back to step 2.

### Rule 6: "Numerically exact" is only claimable over branches that have been exercised

The fixture set is complete when `gcov` on the reference build shows every reachable branch in the translated
subroutines hit at least once. NACA 0012/4412 at Re=1e6, M=0, Ncrit=9 leaves the sharp-TE path, all compressibility
sensitivities, laminar separation, the MRCHDU fallback, forced transition, `MATYP≠1` and RLX limiting **dead**.
Cases that exercise each of those are part of the validation set, not extras. The parameter/fuzz harness optimises
branch coverage; it reports per-variable ULP distributions, iteration-count and branch-flip counts, and
threshold-straddling cases — not a single worst-case number.

**The measurement is mechanised.** `scripts/xfoil-build.sh --gcov` builds the pristine DP reference with
`-fprofile-arcs -ftest-coverage`; `cargo xtask coverage` runs every tracked case through it from its tracked
`xfoil.inp`/`panels.dat`, reads the counters back with gcov and writes `docs/validation/coverage.md`.
`xtask/fixtures-config/coverage.toml` names the translated subroutine set (77 subroutines, 10 files) and carries the
annotations: `[[unreachable]]` entries with a class (*structural*, *mode*, *guard*, *compiler*) and a reason,
`[[note]]` entries recording how an open branch can be reached, `[[dead]]` entries for subroutines XFOIL never
calls. Every never-taken branch that is not annotated is **open**; annotations that stop matching are reported
as stale. Measured 2026-09-04 over 20 cases: 1320 branches, 1125 taken, 90 open (every one with a reach note,
8 of them Newton-loop iteration caps), 84 annotated unreachable, 21 DO-loop zero-trip edges. The first
measurement found one translation gap — OPER `DAMP` (IDAMPV=1, `DAMPL2`) was reachable and not translated — now
closed and gated by the `naca0012_n60_a2_re1e6_damp` case.

XFOIL-independent invariants are also required, because two codes can share a misunderstanding: symmetric airfoil at
α=0 → CL=CM=0 to the noise floor; mirrored airfoil at −α; Blasius flat plate.

### Rule 7: Fixtures are intentional, tracked, and never silently skipped

- `.gitignore` excludes generated files by default (`*.dat`, `*.log`, `fort.*`, `*.bl`, `.tmp/`) and **allowlists**
  what is deliberately committed (`!tests/fixtures/**`, `!xfoil/third-party/**`). Committing a fixture is a decision.
- Tracked fixture budget ~25 MB; full-resolution cases live behind `cargo xtask fixtures --big`.
- Every fixture directory carries a `manifest.json`: source SHA256, patch-series SHA, `gfortran --version`,
  FFLAGS, host triple, libm/OS.
- A missing fixture **fails** the test (`require_fixture`), never `eprintln!("Skipping")` + green.
- `#[ignore]` is tagged with its unblocking stage and tracked in `tests/IGNORED.txt`; CI fails on drift.
- A test with no assertion is not a test. Diagnostic scripts go in `examples/attic/`.

## Project intent

It's XFoil, but in Rust, with unit and integration tests and easily-automated I/O. YFoil reproduces the core
analysis functionality of XFOIL (Mark Drela, MIT) exactly. The name follows the tradition: "y" comes after "x",
but the true choice is a little more esoteric. "Why"?

**Licence:** GPL, because XFOIL is vendored in-tree.

**Scope:** full OPER *analysis* parity — fixed-alpha and fixed-CL (`SPECAL`/`SPECCL`/`MRCL`), forced transition
(`XSTRIP`), polars, Kármán–Tsien compressibility. **Not** in scope: inverse design, interactive mode, multi-element,
flap hinge moments, `FCPMIN`, XFOIL-format output files, 3D effects.

## Architecture

XFOIL uses Fortran COMMON blocks; YFoil passes state explicitly. Pure functions for closures, influence
coefficients and splines; explicit data passing everywhere; plotting behind a cargo feature.

**Inside `src/solver/` and `src/bl/`, the BL state mirrors XFOIL's data model exactly**: two sides with the wake
appended to side 2 (`NBL(2) = IBLTE(2) + NW`), explicit `IPAN/VTI/ISYS/IBLTE/NBL/ITRAN`, and 1-based station
arrays with a dummy slot 0, so translated lines read identically to the Fortran and can be diffed by eye. There is
no third "wake surface". Elsewhere in the crate, idiomatic Rust.

```
geometry/   - Airfoil coordinates, splines, paneling, NACA generation
panel/      - Inviscid panel method, influence matrices, wake geometry, DIJ
bl/         - BL closures, BLDIF/TRDIF/TESYS, MRCHUE, MRCHDU, BLSOLV
solver/     - Pointer layer (IBLPAN/XICALC/IBLSYS/STFIND/STMOVE), velocity layer (UESET/QVFUE/GAMQV…),
              PSILIN/GGCALC/XYWAKE/QDCALC, SETBL, UPDATE, CPCALC/CLCALC/CDCALC (clcalc.rs), VISCAL,
              SPECAL (specal.rs), the OPER session and polar driver (analysis.rs)
output/     - Results serialisation, optional plotting
```

The variable mapping in `docs/xfoil-reference/xfoil-to-yfoil-mapping.md` must be kept current whenever a struct
field or state variable is added, renamed or re-represented.

## CLI

```
yfoil geometry - convert | naca | repanel | info | plot        (alias: geom)
yfoil analyze  - single operating point (--alpha or --cl)
yfoil polar    - alpha sweep (state machine: 0°→max, reinitialise, 0°→min, stitched ascending);
                 --label names the curve, -o writes the polar JSON
yfoil plot     - analysis <file> | polar <file>... (feature-gated; several polars overlay for comparison)
```

Analysis commands accept JSON geometry only; use `yfoil geometry convert` for `.dat`. Typical session:

```bash
yfoil geometry naca 4412 -n 160 -o naca4412.json
yfoil polar naca4412.json --alpha-min -5 --alpha-max 15 --alpha-step 0.5 -r 1e6 --iterations 20 \
    --label "NACA 4412" -o naca4412_polar.json
yfoil plot polar --title "NACA 0012 vs 4412" naca0012_polar.json naca4412_polar.json -o compare.svg
```

## Reference implementation

Never look above `/Users/thc29/source/thclark/yfoil/`.

**Layout:**

```
xfoil/third-party/xfoil-6.99/     - pristine XFOIL 6.99 source (MIT tarball, sha256 5c025064…), tracked,
                                    with a per-file checksum manifest (xfoil/third-party/xfoil-6.99.sha256)
xfoil/instrumentation/build/      - build-config patches applied to BOTH reference builds (DP, -ffp-contract=off,
                                    -ffixed-line-length-none, macOS X11 paths)
xfoil/instrumentation/instrument/ - instrumentation patches (xlog harness + per-file dumps), instrumented build only
target/xfoil-ref/{pristine,instrumented,snan}/ - build output; never edit in place
scripts/xfoil-build.sh            - the build; `cargo xtask xfoil-build [--verify] [--snan]` wraps it
```

Two proofs are part of the reference build and re-run nightly:

1. **Inert instrumentation** — pristine-DP and instrumented-DP produce byte-identical `cp.dat`, `bl.dat`
   and OPER summaries on the smoke case (`--verify`). Proven 2026-09-03; the patch-series build is also
   byte-identical to the previous in-place instrumented binary.
2. **No uninitialised reads** — a `-finit-real=snan -ffpe-trap=invalid` build runs the smoke case clean
   (`--snan`). Proven 2026-09-03. The only raised (untrapped) flag is `IEEE_DIVIDE_BY_ZERO` from
   `PLTINI` (`xplots.f:37`, plot-scale setup that runs even with graphics off) — outside the solver.

All instrumentation writes into the **working directory** (never `/tmp`) at **`ES24.16`** — 17 significant
figures, the minimum that round-trips an f64. The original `E24.16`/`E18.10`/`E12.4` formats were 16, 10 and
4 significant figures respectively and were the reason earlier fixtures were not bit-exact. Prefer the
`xlog.f` JSON harness (`LOGENTER/LOGINT/LOGREAL/LOGBOOL/LOGARR1/LOGARR2`); per-file `.dat` dumps are legacy
and are migrated stage by stage.

Two toolchain facts that cost a ULP each until fixed, now enforced in code: `.dat` files are written with
`{:.17e}` (fixed-point `{:22.16}` lost digits for small `y`), and `serde_json` is built with
`float_roundtrip` (its default parser is best-effort).

Key Fortran files: `xpanel.f` (PSILIN, GGCALC, QDCALC, XYWAKE, IBLPAN, XICALC, STFIND, UESET…), `xbl.f` (SETBL,
IBLSYS, MRCHUE, MRCHDU, UPDATE, DSLIM), `xblsys.f` (BLDIF, TRDIF, TESYS, closures), `xoper.f` (VISCAL, SPECAL,
SPECCL), `xsolve.f` (BLSOLV, GAUSS, LUDCMP/BAKSUB), `xfoil.f` (CLCALC, CDCALC, MRCL, COMSET, PANGEN), `xutils.f`
(SETEXP — port it character-for-character; its 1e-5 Newton tolerance is load-bearing).

Non-interactive XFOIL scripts always start with `PLOP` / `G F` to disable graphics. Defaults that must be pinned
explicitly in every equivalence run because the two codes' defaults differ: `ITER` (XFOIL `ITMAX=20`), `VACCEL`
(0.01), `WAKLEN` (1.0), `TYPE` (MATYP/RETYP; 1 unless the case says otherwise — `MRCL` is gated by the TYPE 2 case).

## Fixture pipeline

`cargo xtask fixtures [--case NAME] [--verify] [--big]` reads `xtask/fixtures-config/cases.toml`, builds the reference if
needed, generates panels **with YFoil**, runs instrumented XFOIL via `LOAD`, asserts the bitwise geometry
handoff, and keeps the raw dumps plus `manifest.json` under `tests/fixtures/xfoil/<case>/` (`track = true`,
budget 8 MB) or `target/fixtures/<case>/`. Case options: `alphas`, `alphas_after_reinit` (INIT between), `polar = true`
(drive with `ALFA a0 / ASEQ a1 aN da` like the polar procedure; ASEQ points get ITMAX+5), `cls` (OPER `CL x`
points), `matyp` (OPER `TYPE n`), `minimal = true` (keep only the `viscal_*.dat` records — for coverage cases). The tracked CI reference case is
`naca0012_n60_a2_re1e6` (`tests/fixtures/mod.rs::REF_CASE`). Stage-specific JSON parsers are added as each
plan stage lands. `--verify` regenerates and asserts byte-identity with what is tracked (same host; cross-host is an ULP
budget). `cargo xtask coverage [--big] [--rebuild]` is the Rule 6 measurement over the same cases. Per-case directories, so a sweep of thousands of runs is just more directories and differencing is a
directory walk.

## Testing

```
tests/
├── utilities/          - shared helpers; tolerances.rs holds the only tolerance constants
├── fixtures/           - tracked XFOIL fixtures (allowlisted in .gitignore) + loaders
├── cli_*_tests.rs      - CLI behaviour
├── integration_*_tests.rs
└── xfoil_*_tests.rs    - equivalence against XFOIL fixtures
```

All test files end in `_tests.rs`. CI jobs: `lint`, `unit` (zero ignores), `fixtures`, `ignore-drift`,
`no-deviations`, `examples`, and nightly `xfoil-parity` (rebuild reference, regenerate, compare).

**Validation reports are generated, never hand-fed.** `docs/validation/` holds only Markdown and SVG; every
number in it is derived from a fixture directory (`tests/fixtures/xfoil/<case>/` or `target/fixtures/<case>/`
for `--big` cases) or from `tests/fixtures/subroutines/`, at the time the report is generated. No XFOIL dump,
`.dat`, `.pol`, `DUMP`/`CPWR` output or other reference data is ever committed under `docs/` — `.gitignore` enforces
it — and a report that cannot be regenerated from tracked inputs plus `cargo xtask fixtures` is not evidence.
Current generators: `cargo xtask coverage` (coverage.md), `scripts/noise-floor.sh` (noise-floor.md),
`generate_subroutine_validation` (subroutines/). The analysis report (BL-distribution difference tables and
plots for the ±15° N=160 polar) is generated from the `--big` fixture case and is the remaining validation stage.

## Conventions

- Coordinates normalised by chord; panel ordering TE → upper → LE → lower → TE, counter-clockwise (XFOIL's `LOAD`
  reverses clockwise input, which would break the bitwise handoff).
- Angles in radians internally, degrees in the CLI. SI units.
- Integer powers written as explicit multiplications where bit-exactness matters (`x*x*x`, not `powi(3)`) to match
  gfortran's inline expansion.
- Temporary files go in `.tmp/`, never `/tmp`.

## Polar sweep procedure

A polar is a state machine, not N independent solves. Model XFOIL's flags (`LWAKE, LIPAN, LBLINI, LADIJ, LWDIJ,
LVCONV, LGAMU, LQAIJ`) and `AWAKE/AVISC/MVISC` explicitly; the previous alpha's BL is the next one's initial
condition. Sweep 0°→+max, reinitialise, −1°→−min, stitch ascending.

XFOIL equivalent:

```
PLOP
G F

LOAD geometry.dat
OPER
VISC 1000000
ITER 20
PACC
polar_pos.txt

ALFA 0
ASEQ 1 15 1
PACC

INIT
PACC
polar_neg.txt

ALFA -1
ASEQ -2 -15 -1

QUIT
```

The `.pol` file is for humans; the gate is the instrumented per-alpha dump.

## Validation set

Minimum cases, each tabulated and plotted with BL-variable difference tables:

- NACA 0012, α=0 — symmetry, base viscous case, XFOIL-independent invariants
- NACA 0012, 0°→+1° and 0°→−1° — initialisation from a previous solution, sign checks
- NACA 4412, full ±15° polar — beyond the limits of convergence
- Plus the branch-coverage cases of Rule 6: sharpened TE, M=0.3, high-α separated, low-Re laminar separation,
  `XSTRIP`, `MATYP≠1`

BL distributions are extracted at 0°, ±5°, ±10°, ±15°, but every intermediate angle is computed so initialisation
matches XFOIL.
