# YFoil - Viscous and inviscid 2d Aerofoil Analyses

It's XFoil, but in Rust, with unit and integration tests, and an easily-automated I/O.

## Project Intent

YFoil reproduces the core analysis functionality of XFOIL (created by Mark Drela at MIT in the 1980s) in Rust.

The name follows the tradition: "y" comes after "x", but the true choice is a little more esoteric. "Why

## What YFoil Does

YFoil is intended to replicate XFoil as near to exactly as possible. Features include:

- Generation of aerofoil geometry
- Panel method for 2D airfoil sections
    - Inviscid solver
    - Integral boundary layer solver with eN transition prediction
    - Coupled viscous-inviscid iteration
- Polar sweep generation (from 0° up, then 0° down)
- Compressibility corrections (Karman-Tsien for subsonic Mach)

## What YFoil Does NOT Do

- No inverse design (Qspec manipulation for aerofoil design)
- No interactive mode (CLI is one-shot commands only)
- No multi-element airfoils
- No 3D effects

## Architecture Principles

### Decoupled Design (Unlike XFOIL)

XFOIL uses Fortran COMMON blocks for global state, tightly coupling all modules. YFoil explicitly avoids this:

1. **Pure Functions**: Closure relations, influence coefficients, splines are pure functions with no side effects
2. **Explicit Data Passing**: Each function takes inputs and returns outputs, no global state
3. **Trait-Based Extension**: Transition models, closures can be swapped via traits
4. **Feature-Gated Optionals**: Plotting behind cargo feature flag

### Module Responsibilities

```
geometry/   - Airfoil coordinates, splines, paneling, NACA generation
panel/      - Inviscid panel method, influence matrices, wake
bl/         - Boundary layer solver, closures, transition, marching
solver/     - Viscous-inviscid coupling (VISCAL), polar sweeps
forces/     - Pressure coefficients, CL/CD/CM integration
output/     - Results serialization, optional plotting
```

## CLI Structure

```
yfoil geom     - Geometry operations to create a set of panels (convert, naca, repanel, smooth)
yfoil analyze  - Single operating point analysis for a given paneling
yfoil polar    - Alpha sweep
yfoil plot     - Visualization (feature-gated)
```

Analysis commands accept JSON geometry only. Use `yfoil geom convert` to convert .dat files.

## Reference Implementation

Everything needed is in this repository folder. Never look above `/Users/thc29/source/thclark/yfoil/`.

XFOIL 6.99 source code and compiled binary are in `xfoil/xfoil6.99/`. Key source files in `src/`:

- `xpanel.f` - Panel method (PSILIN, influence coefficients)
- `xbl.f` - BL solver (SETBL, MRCHUE, MRCHDU)
- `xblsys.f` - BL Newton system, transition (DAMPL, AXSET)
- `xoper.f` - VISCAL coupling loop
- `XFOIL.INC` - Data structure definitions

The compiled XFOIL binary is located at `xfoil/xfoil6.99/bin/xfoil`. Always use this binary for generating reference
data - do NOT search for xfoil elsewhere.

When running XFOIL non-interactively (via script/stdin), always start the script with:

```
PLOP
G F
```

This disables graphics mode to prevent "Cannot open display" errors.

## Variable Mapping Documentation

The file `docs/xfoil-reference/xfoil-to-yfoil-mapping.md` contains a comprehensive mapping between XFOIL Fortran
variables/common blocks and their YFoil Rust equivalents.

**Keep this mapping updated** when:

- Adding or renaming struct fields in YFoil
- Changing how XFOIL variables are represented in YFoil
- Adding new BL or solver state variables
- Modifying closure result structures

This document is essential for debugging discrepancies between XFOIL and YFoil behavior.

## Development Approach

Build each module to completion with full tests before proceeding to the next. Order:

1. Geometry (splines, paneling, NACA)
2. Inviscid panel method
3. Boundary layer solver
4. Viscous-inviscid coupling
5. Polar sweeps and output

## CRITICAL: Numerical Precision Requirements

**YFoil MUST produce numerically identical results to XFOIL.** This is not negotiable.

- 5-10% error is UNACCEPTABLE - this indicates a bug, not acceptable tolerance
- The only acceptable differences are floating-point rounding errors (typically < 1e-10)
- There should be NO differences in calculation logic between YFoil and XFOIL
- Every formula, coefficient, and algorithm must match XFOIL exactly
- When in doubt, instrument XFOIL to verify the exact values being computed

**NO WORKAROUNDS. NO APPROXIMATIONS.**

If something doesn't match XFOIL output, it is a BUG. The objective is not to replicate "something like XFOIL" - it is
to TRANSLATE XFOIL IDENTICALLY. Do not:

- Adjust relaxation factors to "make it stable"
- Limit iterations to "avoid divergence"
- Change signs to "make it work"
- Use simplified approaches instead of the actual XFOIL algorithms

If XFOIL uses a Newton system, implement the Newton system. If XFOIL uses a specific formula, use that exact formula.
Any deviation from XFOIL's actual implementation is wrong and must be fixed, not worked around.

**Debugging Approach:**

1. Instrument XFOIL Fortran source with WRITE statements to output intermediate values
2. Recompile the instrumented XFOIL binary
3. Compare yfoil output at each step against instrumented XFOIL output
4. Values must match to machine precision (typically 12+ significant figures)

**Instrumentation Process:**

1. Add WRITE statements to the relevant XFOIL subroutine
2. Rebuild XFOIL: `cd xfoil/xfoil6.99 && make clean && make`
3. Run instrumented XFOIL with identical inputs
4. Parse output and compare against yfoil values
5. Any discrepancy > 1e-10 relative error indicates a bug to fix

## HARD RULE: Identical Geometry for XFOIL vs YFoil Comparisons

**When comparing XFOIL and YFoil results, you MUST use the EXACT SAME panel coordinates.**

This is non-negotiable. Different paneling produces different results, making any comparison meaningless.

**Required workflow for any XFOIL vs YFoil comparison:**

1. Generate geometry with YFoil: `yfoil geom naca 0012 --paneler xfoil -o geometry.json`
2. Export to .dat format: `yfoil geom convert geometry.json -o geometry.dat`
3. Run XFOIL with that .dat file: `LOAD geometry.dat` (do NOT use `NACA 0012` then `PANE`)
4. Run YFoil with the same geometry.json

**NEVER do this:**

- Use `NACA 0012` + `PANE` in XFOIL while using `naca_4digit("0012", N)` in YFoil
- Compare results without verifying panel coordinates match exactly
- Assume similar panel counts mean identical geometry

**Verification:** Before any comparison, print the first few panel coordinates from both XFOIL and YFoil and confirm
they match to machine precision.

## Testing Strategy

- Tests for pure functions (closures, splines, influence)
- Tests comparing against XFOIL results with EXACT numerical matching
    - Run instrumented XFOIL to get input and output values to create fixtures
- Regression tests with stored known-good outputs
- Test airfoils: NACA 0012 (symmetric), NACA 4412 (cambered)
- Instrumented comparison tests at each solver stage

### Test File Naming Convention

All test files in `tests/` must be suffixed with `_tests.rs` so they are easily identifiable from filename tabs.

**Structure:**

```
tests/
├── utilities/                        - Shared test utilities (not run as tests)
│   └── mod.rs
├── fixtures/                         - XFOIL-generated test fixture data
│   └── mod.rs                        - Fixture loading utilities
├── cli_*_tests.rs                    - CLI command tests (yfoil geom, analyze, etc.)
├── integration_*_tests.rs            - Library module integration tests
├── xfoil_*_tests.rs                  - XFOIL validation tests (exact numeric comparison)
└── ...
```

**Naming patterns:**

- `cli_` prefix: Tests for CLI commands (e.g., `cli_geom_tests.rs`, `cli_analyze_tests.rs`)
- `integration_` prefix: Integration tests for library modules (e.g., `integration_geometry_tests.rs`)
- `xfoil_` prefix: Tests that validate against XFOIL output (exact numeric matching)
- `_errors` suffix: Tests specifically for error handling

## Conventions

- Coordinates normalized by chord (x/c, y/c)
- Panel ordering: TE → upper surface → LE → lower surface → TE
- Angles in radians internally, degrees in CLI
- SI units throughout

## Development Workflow

For temporary files and debug scripts, use the `.tmp/` directory in the repo root instead of `/tmp`. This avoids
permission issues and keeps debug artifacts with the project.

## Polar Sweep Procedure

**CRITICAL: The correct order to run a polar sweep is:**

1. **Start at α = 0°** - this provides a good initial condition
2. **Sweep upward** in increments of +1° from 0° to +15° (or desired max)
3. **Reinitialize** the solution at α = -1° - this prevents non-converged high-alpha results from polluting the negative sweep
4. **Sweep downward** in increments of -1° from -1° to -15° (or desired min)
5. **Stitch results** in ascending order from -15° to +15° for clean output

This procedure ensures:
- Each point uses the previous converged solution as initialization (good for convergence)
- The negative alpha sweep starts fresh, not from a potentially diverged high-alpha state
- Final output appears as a single continuous polar

**For XFOIL (non-interactive):**
```
PLOP
G F

LOAD geometry.dat
PCOP
OPER
VISC 1000000
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

Then stitch `polar_neg.txt` (reversed) + α=0 from `polar_pos.txt` + rest of `polar_pos.txt`.

## Full validation against XFOIL

There should be a set of validations maintained in the final documentation (these cases may also be useful for debugging
purposes).

- NACA0012 at 0 degree angle of attack (for paneling, inviscid solution and base case viscous solutions which check
  symmetry in the solution)
- NACA0012 at 0 then 1 degree angle of attack (for inviscid and viscous solutions and polars, where the subsequent value
  of the polar at 1 degree is initialised from the first solution at 0 degrees.
- NACA0012 at 0 and then -1 degree angle of attack (to test the same in reverse angles to check signs)
- NACA4412 at 0-15 degrees, then reinitialised at 0, then -1--15 as a complete polar to test beyond limits of divergence

Each of these cases should be tabulated and plotted, including differences between all boundary layer variable
distributions.

All xfoil runs should be run with the same default number of iterations as yfoil has.

For every validation, test comparison output or debugging run of xfoil, it's imperative to use the same panels -
now we have the geometry working, do that by generating the geometry then repaneling it using the xfoil-based paneler (
this can be done with yfoil geom).

### Validation Tool Locations

Validation generators are binaries in `src/bin/`, NOT examples. These tools must be re-run whenever YFoil changes to
verify correctness.

**Source code:**
```
src/bin/
├── generate_geometry_validation.rs    - Geometry validation (paneling)
├── generate_analysis_validation.rs    - Analysis validation (viscous solver)
├── generate_subroutine_validation.rs  - Subroutine validation (BL closures)
└── ...
```

**Running validation generators:**
```bash
# Geometry validation
cargo run --bin generate_geometry_validation --features plotting

# Analysis validation (regenerate YFoil data and plots)
cargo run --bin generate_analysis_validation --features plotting

# Analysis validation (also regenerate XFOIL baseline data)
cargo run --bin generate_analysis_validation --features plotting -- --run-xfoil

# Subroutine validation (BL closure relations)
cargo run --bin generate_subroutine_validation
```

**XFOIL scripts and assets:**
```
docs/validation/assets/
├── geometry/
│   ├── xfoil/           - XFOIL reference .dat files
│   ├── yfoil/           - YFoil generated .json files
│   └── plots/           - Comparison SVG plots
└── analysis/
    ├── geometry/        - Shared geometry files for analysis
    ├── scripts/         - XFOIL .xfoil scripts (sweep scripts only)
    ├── xfoil/           - XFOIL output (polars, BL dumps, Cp files)
    ├── yfoil/           - YFoil output
    └── plots/           - Comparison SVG plots
```

**Generated documentation:**
```
docs/validation/
├── geometry/
│   └── README.md        - Geometry validation report
├── analysis/
│   ├── README.md        - Analysis validation overview
│   ├── naca0012.md      - NACA 0012 detailed results
│   └── naca4412.md      - NACA 4412 detailed results
└── subroutines/
    ├── README.md        - Subroutine validation overview
    ├── closure.md       - Closure functions (HKIN, HSL, etc.)
    └── transition.md    - Transition (DAMPL)
```

**Subroutine validation fixtures** (used by tests and validation):
```
tests/fixtures/subroutines/
├── hkin/     - Kinematic shape factor fixtures
├── hsl/      - Laminar energy shape factor fixtures
├── hst/      - Turbulent energy shape factor fixtures
├── cfl/      - Laminar skin friction fixtures
├── cft/      - Turbulent skin friction fixtures
├── dil/      - Laminar dissipation fixtures
└── dampl/    - Amplification rate fixtures
```

These fixtures were generated from instrumented XFOIL (see `docs/validation/subroutines/README.md`
for details on the instrumentation added to `xfoil/xfoil6.99/src/xbl.f`).

### Analysis Validation Methodology

Analysis validation uses proper polar sweeps to ensure solutions are correctly initialized:

1. **Positive sweep**: 0° → 1° → 2° → ... → 15° (each angle initialized from previous)
2. **Reinitialize** at α = -1°
3. **Negative sweep**: -1° → -2° → ... → -15° (each angle initialized from previous)

BL distributions are extracted only at validation angles (0°, ±5°, ±10°, ±15°), but ALL intermediate angles are computed
to ensure proper initialization. This applies to both XFOIL scripts and YFoil code.                                                                              
