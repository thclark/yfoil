# XFOIL Reference Documentation

This directory contains detailed documentation of XFOIL 6.99 for the purpose of:

1. Validating yFoil produces numerically identical results
2. Understanding the algorithm flow for debugging

## Overview

XFOIL (created by Mark Drela at MIT, 1980s) performs viscous-inviscid analysis of 2D airfoil sections using:

- Panel method with vortex distribution for inviscid flow
- Integral boundary layer solver with e^N transition prediction
- Coupled viscous-inviscid iteration (VISCAL)

## Directory Structure

```
xfoil-reference/
├── README.md                 # This file
├── architecture.md           # Systems diagram description
├── architecture.mmd          # Mermaid diagram source
├── architecture.pdf          # Generated PDF (see below)
├── common-blocks/
│   ├── XFOIL.INC.md         # Main data structures
│   ├── XBL.INC.md           # BL local variables
│   └── BLPAR.INC.md         # BL closure parameters
└── modules/
    ├── xoper.md             # VISCAL coupling loop
    ├── xblsys.md            # BL Newton system
    ├── xbl.md               # BL marching
    ├── xpanel.md            # Panel method
    ├── spline.md            # Spline utilities
    ├── xgeom.md             # Geometry operations
    └── naca.md              # NACA airfoil generation
```

## XFOIL Source Files

| File       | Purpose                      | Key Subroutines                     |
|------------|------------------------------|-------------------------------------|
| `xoper.f`  | Main operations and coupling | VISCAL, SPECAL, SPECCL              |
| `xblsys.f` | BL Newton system             | BLPRV, BLKIN, BLVAR, BLSYS, TRCHEK2 |
| `xbl.f`    | BL marching                  | SETBL, MRCHUE, MRCHDU, UPDATE       |
| `xpanel.f` | Panel method                 | PSILIN, QDCALC, UICALC, QVFUE       |
| `spline.f` | Cubic splines                | SPLINE, SEVAL, SPLIND               |
| `xgeom.f`  | Geometry                     | PANGEN, LEFIND, TECALC              |
| `naca.f`   | NACA generation              | NACA4, NACA5                        |

## Data Flow

```
Input Geometry (X, Y coordinates)
        │
        ▼
┌───────────────────┐
│   PANGEN          │  Repanel airfoil
│   (xgeom.f)       │
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   SPECAL          │  Specify alpha
│   (xoper.f)       │
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   Panel Method    │  Solve inviscid flow
│   (xpanel.f)      │  → GAM, QINV
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   VISCAL          │  Viscous-inviscid coupling
│   (xoper.f)       │
│   ┌─────────────┐ │
│   │ SETBL       │ │  Initialize BL
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ MRCHUE      │ │  March BL equations
│   │ (xbl.f)     │ │  → DSTR, THET
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ UPDATE      │ │  Mass defect → source
│   │ (xbl.f)     │ │  → SIG
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ QDCALC      │ │  Update Qvis
│   └─────────────┘ │
│         ↑         │
│         └─────────│  Iterate until converged
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   Force Coeffs    │  CL, CD, CM
└───────────────────┘
```

## yFoil equivalents

The full XFOIL → yFoil name mapping is [`xfoil-to-yfoil-mapping.md`](xfoil-to-yfoil-mapping.md);
the naming rules are [`docs/conventions/naming.md`](../conventions/naming.md). Where the
translated routines live:

| XFOIL             | yFoil                                                  |
|-------------------|--------------------------------------------------------|
| `xoper.f:VISCAL`  | `src/solver/viscal.rs::solve_viscous`                  |
| `xoper.f:SPECAL`  | `src/solver/specal.rs::solve_inviscid_at_alpha`        |
| `xbl.f:SETBL`     | `src/solver/setbl.rs::assemble_newton_system`          |
| `xbl.f:UPDATE`    | `src/solver/update.rs::apply_newton_update`            |
| `xbl.f:MRCHUE`    | `src/bl/mrchue.rs::march_direct`                       |
| `xbl.f:MRCHDU`    | `src/bl/mrchdu.rs::march_prescribed_dstar`             |
| `xblsys.f:BLKIN`  | `src/bl/station.rs::StationState::set_kinematic_variables` |
| `xblsys.f:BLVAR`  | `src/bl/station.rs::StationState::set_closure_variables`   |
| `xblsys.f:BLSYS`  | `src/bl/blsys.rs::assemble_interval_system`            |
| `xblsys.f:BLDIF`  | `src/bl/difference.rs::IntervalSystem::assemble_interval_equations` |
| `xblsys.f:TRCHEK2`| `src/bl/transition.rs::check_transition`               |
| `xsolve.f:BLSOLV` | `src/bl/blsolv.rs::solve_newton_system`                |
| `xpanel.f:PSILIN` | `src/solver/psilin.rs::panel_influence`                |
| `xpanel.f:GGCALC` | `src/solver/ggcalc.rs::build_inviscid_system`          |

Every translating function carries `#[doc(alias = "XFOIL NAME")]`, so `cargo doc` search by the
Fortran name finds it.

## Generating Architecture PDF

To generate the architecture diagram PDF from the Mermaid source:

```bash
# Install mermaid-cli if not present
npm install -g @mermaid-js/mermaid-cli

# Generate PDF
cd docs/xfoil-reference
mmdc -i architecture.mmd -o architecture.pdf -b white
```

Alternatively, you can use the online Mermaid Live Editor at https://mermaid.live/ to render and export the diagram.

## Usage

This documentation serves as a reference when:

1. Debugging mismatches between yFoil and XFOIL
2. Understanding algorithmic details for implementation
3. Writing test fixtures and validation cases
4. Preparing paper content describing the method
