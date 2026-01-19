# YFoil Validation Documentation

This directory contains validation documentation comparing YFoil results against XFOIL reference data.
The goal is to verify that YFoil produces numerically identical results to XFOIL at machine precision.

## Validation Sections

### [Geometry Validation](geometry/README.md)

Validates YFoil's geometry generation (NACA airfoil generation and XFOIL-style paneling) against
XFOIL reference coordinates.

- **Test Cases:** NACA 0012 and NACA 4412 with 81 and 160 panels
- **Results:** See [geometry/README.md](geometry/README.md) for current metrics

### [Analysis Validation](analysis/README.md)

*(Future)* Validates YFoil's aerodynamic analysis results against XFOIL.

## Regenerating Validation

To regenerate validation after code changes:

```bash
cargo run --bin generate_validation --features plotting
```

This command:
1. Regenerates YFoil geometry for all test cases
2. Recomputes comparison metrics against XFOIL baseline
3. Regenerates comparison SVG plots
4. Updates the markdown reports

**Note:** XFOIL reference data is committed as the baseline and is NOT regenerated.
If XFOIL reference data needs updating, run the scripts in `assets/geometry/scripts/`.

## Directory Structure

```
docs/validation/
├── README.md                              # This index page
├── geometry/
│   └── README.md                          # Geometry validation report
├── analysis/
│   └── README.md                          # Analysis validation (future)
└── assets/
    └── geometry/
        ├── xfoil/                         # XFOIL reference .dat files
        ├── yfoil/                         # YFoil generated .json files
        ├── plots/                         # Comparison SVG plots
        └── scripts/                       # XFOIL generation scripts
```

## Design Principles

1. **XFOIL Baseline**: Reference data is generated once by XFOIL and committed.
   This provides a stable baseline for comparison.

2. **Machine Precision**: YFoil must match XFOIL to floating-point precision (< 1e-12).
   Larger errors indicate bugs, not acceptable tolerance.

3. **Visual Verification**: SVG plots allow humans to visually verify that
   geometries are coincident.

4. **Automated Updates**: The `generate_validation` binary automatically updates
   all reports and plots when run.
