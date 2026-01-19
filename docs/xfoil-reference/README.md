# XFOIL Reference Documentation

This directory contains detailed documentation of XFOIL 6.99 for the purpose of:

1. Validating YFoil produces numerically identical results
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

## YFoil Equivalents

| XFOIL             | YFoil                    |
|-------------------|--------------------------|
| `xoper.f:VISCAL`  | `src/solver/viscal.rs`   |
| `xblsys.f:BLKIN`  | `src/bl/closure.rs`      |
| `xblsys.f:BLVAR`  | `src/bl/secondary.rs`    |
| `xblsys.f:BLSYS`  | `src/bl/system.rs`       |
| `xbl.f:MRCHUE`    | `src/bl/march.rs`        |
| `xpanel.f:PSILIN` | `src/panel/influence.rs` |

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

1. Debugging mismatches between YFoil and XFOIL
2. Understanding algorithmic details for implementation
3. Writing test fixtures and validation cases
4. Preparing paper content describing the method
