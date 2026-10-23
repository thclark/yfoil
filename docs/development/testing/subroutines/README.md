---
icon: lucide/sigma
---

# Subroutine comparison

yFoil's boundary-layer closure functions against the calls the instrumented XFOIL reference made to the same functions.

## Summary

| Subroutine | Test Cases | Passed | Failed | Max Relative Error | Status |
|------------|------------|--------|--------|-------------------|--------|
| HKIN | 200 | 200 | 0 | 0.00e0 | ✓ |
| CFL | 200 | 200 | 0 | 0.00e0 | ✓ |
| HSL | 200 | 200 | 0 | 0.00e0 | ✓ |
| DIL | 200 | 200 | 0 | 0.00e0 | ✓ |
| HST | 200 | 200 | 0 | 0.00e0 | ✓ |
| CFT | 200 | 200 | 0 | 0.00e0 | ✓ |
| DAMPL | 200 | 200 | 0 | 2.20e-15 | ✓ |

## What is compared

Each fixture is one call of the closure inside a real run of the instrumented double-precision reference (CLAUDE.md, *Reference implementation*): its inputs and outputs, dumped at 17 significant digits (`ES24.16`, the minimum that round-trips an f64) while the reference ran the cases below. yFoil's function is called on the same inputs here and its outputs compared with the dumped ones. The error is relative, absolute below 1e-15.

| Case | Reference build host | libm |
|------|------|------|
| `closures_naca0012_n60_re1e6` | arm64-Darwin-24.6.0 | CORE-MATH commit=85430b3ef7d2beb0de16d7c569ffa9bb391615ed (-O2 -ffp-contract=off); host /usr/lib/libSystem.B.dylib 1345.120.2) |
| `closures_naca0012_n60_re1e6_m03` | arm64-Darwin-24.6.0 | CORE-MATH commit=85430b3ef7d2beb0de16d7c569ffa9bb391615ed (-O2 -ffp-contract=off); host /usr/lib/libSystem.B.dylib 1345.120.2) |

This report was generated on `aarch64-macos`. Both codes take their transcendentals (`exp`, `log`, `tanh`, `**`) from the same correctly rounded functions (CORE-MATH, `cr-libm`), whose results do not depend on the host, and the reference is built with `-ffp-contract=off`, so a translation that evaluates each expression in XFOIL's order gives the same bits: an error of exactly 0 means bit-identical, which is what most closures show. A non-zero error means the translation does not evaluate the expression bit for bit as XFOIL does: DAMPL's, a few ULP (the table above), is within `TOL_PURE` and has not been traced to the operation.

## Tolerance

A case passes within **1e-12**, the value of `TOL_PURE` in `tests/common/utilities/tolerances.rs`, to which the tests in `tests/subroutine/closures.rs` hold the same fixtures.

## Test Fixtures

Fixture data is stored in [`tests/fixtures/subroutines/`](../../../tests/fixtures/subroutines/), one directory per closure, one JSON file per call, with the reference's provenance in `manifest.json`. They are cut from the reference's closure log by `cargo xtask fixtures` for the cases with `closures = true` in `xtask/fixtures-config/cases.toml`.

## Regenerating

```bash
# Regenerate this report (uses the tracked fixtures)
cargo run --bin generate_subroutine_validation

# To regenerate the fixtures from the instrumented reference:
cargo xtask fixtures --case closures_naca0012_n60_re1e6 --case closures_naca0012_n60_re1e6_m03
```

## Detailed Reports

- [Closure Functions (HKIN, CFL, HSL, DIL, CFT, HST)](closure.md)
- [Transition (DAMPL)](transition.md)
