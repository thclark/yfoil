# YFoil - Viscous and inviscid 2d Aerofoil Analyses

It's XFoil, but in Rust, with unit and integration tests, and an easily-automated I/O.

## Project Intent

YFoil reproduces the core analysis functionality of XFOIL (created by Mark Drela at MIT in the 1980s) in Rust.

The name follows the tradition: "y" comes after "x", but the true choice is a little more esoteric. "Why

## What YFoil Does

- Viscous/inviscid analysis of 2D airfoil sections
- Panel method for inviscid flow with vortex distribution
- Integral boundary layer solver with eN transition prediction
- Coupled viscous-inviscid iteration (VISCAL)
- Polar sweep generation (from 0° up, then 0° down)
- Compressibility corrections (Karman-Tsien for subsonic Mach)

## What YFoil Does NOT Do

- No inverse design (Qspec manipulation for airfoil design)
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

### Data Flow

```
Geometry → Panel Method → Inviscid Solution
                ↓
         BL Solver ← Edge Velocities
                ↓
         Mass Defect → Source Distribution
                ↓
         Updated Velocities (iterate until converged)
                ↓
         Force Integration → CL, CD, CM
```

## CLI Structure

```
yfoil geom     - Geometry operations (convert, naca, repanel, smooth)
yfoil analyze  - Single operating point analysis
yfoil polar    - Alpha sweep
yfoil plot     - Visualization (feature-gated)
```

Analysis commands accept JSON geometry only. Use `yfoil geom convert` to convert .dat files.

## Reference Implementation

XFOIL 6.99 source code is in `xfoil/xfoil6.99/src/`. Key files:

- `xpanel.f` - Panel method (PSILIN, influence coefficients)
- `xbl.f` - BL solver (SETBL, MRCHUE, MRCHDU)
- `xblsys.f` - BL Newton system, transition (DAMPL, AXSET)
- `xoper.f` - VISCAL coupling loop
- `XFOIL.INC` - Data structure definitions

## Development Approach

Build each module to completion with full tests before proceeding to the next. Order:

1. Geometry (splines, paneling, NACA)
2. Inviscid panel method
3. Boundary layer solver
4. Viscous-inviscid coupling
5. Polar sweeps and output

## Testing Strategy

- Unit tests for pure functions (closures, splines, influence)
- Integration tests comparing against XFOIL results
- Regression tests with stored known-good outputs
- Test airfoils: NACA 0012 (symmetric), NACA 4412 (cambered)

## Conventions

- Coordinates normalized by chord (x/c, y/c)
- Panel ordering: TE → upper surface → LE → lower surface → TE
- Angles in radians internally, degrees in CLI
- SI units throughout
