# YFoil Subroutine Validation Report

This report compares YFoil's boundary layer functions against XFOIL reference values.

## Summary

| Subroutine | Description | Test Cases | Passed | Failed | Max Relative Error | Status |
|------------|-------------|------------|--------|--------|-------------------|--------|
| HKIN | Kinematic shape factor | 200 | 200 | 0 | 3.66e-16 | ✓ |
| CFL | Laminar skin friction | 200 | 200 | 0 | 2.63e-15 | ✓ |
| HSL | Laminar energy shape factor | 200 | 200 | 0 | 2.77e-14 | ✓ |
| DIL | Laminar dissipation | 200 | 200 | 0 | 6.00e-15 | ✓ |
| HST | Turbulent energy shape factor | 200 | 200 | 0 | 6.51e-15 | ✓ |
| CFT | Turbulent skin friction | 200 | 200 | 0 | 2.14e-15 | ✓ |
| DAMPL | Amplification rate | 200 | 200 | 0 | 2.11e-13 | ✓ |
| BLKIN | Kinematic BL variables | 200 | 200 | 0 | <1e-12 | ✓ |
| BLVAR | All secondary BL variables | 200 | 200 | 0 | <1e-10 | ✓ |
| TRCHEK2 | Transition detection | 6 | 6 | 0 | N/A | ✓ |

## Tolerance

Tests use the following tolerances:
- **Closure functions (HKIN, CFL, HSL, DIL, HST, CFT, DAMPL)**: 1e-12 relative tolerance (12 significant figures)
- **BLKIN**: 1e-12 relative tolerance
- **BLVAR**: 1e-10 relative tolerance (slightly relaxed due to chain of computations)
- **TRCHEK2**: Structural validation (fixtures loaded successfully)

## Detailed Reports

- [Closure Functions (HKIN, CFL, HSL, DIL, CFT, HST)](closure.md)
- [Transition (DAMPL)](transition.md)

## Test Execution

Run all subroutine tests with:
```bash
cargo test --test subroutine_tests
```
