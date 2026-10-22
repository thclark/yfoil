# Proofs that the reference is sound

A comparison against the reference is only as good as the reference. Three proofs are part of
building it, and the `reference-build` workflow re-runs them on every change to `xfoil/` or the
build script (and on demand); a fourth check reruns it on another host. If any of these fails, no
comparison on this tab can be trusted.

## Proof that instrumentation changes nothing

The instrumentation writes; it may not change control flow or any stored value. `cargo xtask
xfoil-build --verify` runs the smoke case on the pristine and the instrumented double-precision
builds and requires byte-identical `cp.dat`, `bl.dat` and OPER summaries. Proven 2026-09-03; the
patch-series build is also byte-identical to the earlier in-place instrumented binary.

The smoke case is one converged point. `scripts/xfoil-instrumentation-check` extends the proof
through stall: both builds run one identical OPER script (XFOIL's own `NACA 0012`, `PPAR N 160`,
Re 1e6, 0° to 25° by 0.5°, `ALFA` and `DUMP` per point) and are compared byte for byte on XFOIL's
standard `PACC` polar and `DUMP` files, at XFOIL's own output precision. It writes a dated run
folder (a 2×2 figure, `summary.json`), not a page:

``` sh
cargo run --release -p xfoil-instrumentation-check
```

## Proof of no uninitialised reads

`cargo xtask xfoil-build --snan` builds the reference with `-finit-real=snan
-ffpe-trap=invalid`, so any real read before it is written traps, and runs the smoke case clean.
Proven 2026-09-03. The only floating-point flag raised (untrapped) is `IEEE_DIVIDE_BY_ZERO` from
`PLTINI` (`xplots.f:37`), plot-scale setup that runs even with graphics off, outside the solver.

## Proof of no host maths

Every build fails if any maths function is still resolved from the host libm (`nm -u`), whatever
the compiler lowered an intrinsic to: GCC turns a `SIN`/`COS` pair of one argument into `sincos` on
glibc and into `cexp(i·x)` on Darwin, and both are shimmed onto CORE-MATH. The only exceptions are
MDES's `cabs`/`clog` (inverse design, out of scope) and `sqrt`, which IEEE 754 requires correctly
rounded. Why this matters is [correctly rounded maths in both
codes](correctly-rounded-maths-in-both-codes.md).

## Reproduction on another host

`scripts/cross-host.sh` (`cargo xtask cross-host`, in a Linux container) builds the reference on
glibc (`scripts/cross-host/Dockerfile`), copies the tree into `target/cross-host-tree/` and
requires the reference to reproduce **every tracked record byte for byte**. The
`reference-build` workflow checks the same for the CI reference case, `naca0012_n60_a2_re1e6`, on
Linux. Since both codes moved to correctly rounded maths, any difference between hosts is a bug.
