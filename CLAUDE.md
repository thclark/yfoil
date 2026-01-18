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

## Testing Strategy

- Unit tests for pure functions (closures, splines, influence)
- Integration tests comparing against XFOIL results with EXACT numerical matching
    - Run instrumented XFOIL to get input and output values to create fixtures
- Regression tests with stored known-good outputs
- Test airfoils: NACA 0012 (symmetric), NACA 4412 (cambered)
- Instrumented comparison tests at each solver stage

## Conventions

- Coordinates normalized by chord (x/c, y/c)
- Panel ordering: TE → upper surface → LE → lower surface → TE
- Angles in radians internally, degrees in CLI
- SI units throughout

## Development Workflow

- For temporary files and debug scripts, use the `.tmp/` directory in the repo root instead of `/tmp`. This avoids
  permission issues and keeps debug artifacts with the project.
