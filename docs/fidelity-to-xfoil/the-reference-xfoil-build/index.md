# The reference XFOIL build

Every comparison on this tab is against one build of XFOIL, *the reference*: XFOIL 6.99 compiled
in double precision with pinned flags, linked against correctly rounded maths, and instrumented to
dump its state at full precision. When this documentation says "XFOIL" it means this build, not
the program you download. This section says what the build is, how it is made, and how it is
proved sound — three proofs re-run whenever the reference changes.

## Double precision and pinned compiler flags

Stock XFOIL is single precision, and agrees with a double-precision build only to about 1e-7 —
far coarser than the differences a translation has to be judged by. So the reference is XFOIL 6.99
compiled with `-fdefault-real-8` and linked against the double-precision plot library
(`libPlt_gDP.a`), with `-ffp-contract=off` so that the compiler fuses no multiply-add that the
source does not write, at one fixed optimisation level, on a recorded toolchain (gfortran 14.2).
The gfortran Makefile XFOIL ships has its bounds and floating-point checks commented out; the
reference's own checks are the [proofs](proofs-that-the-reference-is-sound.md).

## Layout and how it is built

The reference is **never modified in place and never built in-tree**.

| Directory | What it holds |
|---|---|
| `xfoil/third-party/xfoil-6.99/` | the pristine MIT tarball (sha256 `5c025064…`), tracked, with a per-file checksum manifest (`xfoil/third-party/xfoil-6.99.sha256`) |
| `xfoil/instrumentation/build/` | build-config patches applied to **both** builds: double precision, `-ffp-contract=off`, `-ffixed-line-length-none`, macOS X11 paths |
| `xfoil/instrumentation/instrument/` | instrumentation patches, the **instrumented build only**: the `xlog` harness and the per-subroutine dumps |
| `third-party/core-math/` | the correctly rounded maths functions both codes use ([correctly rounded maths in both codes](correctly-rounded-maths-in-both-codes.md)) |
| `target/xfoil-ref/{pristine,instrumented,snan}/` | build output; never edited in place |

`scripts/xfoil-build.sh` (wrapped by `cargo xtask xfoil-build [--verify] [--snan]`) stages each
tree under `target/xfoil-ref/`, applies the series listed in `xfoil/instrumentation/series.build`
and `series.instrument`, and builds. Every stage records a `fingerprint`, the hash of everything
the build depends on (sources, patch series, maths shim, build script); every `cargo xtask` that
runs the reference rebuilds a stage whose fingerprint is not current, and the fixture manifests
carry it.

## Pages of the reference build

| Page | What it covers |
| --- | --- |
| [Correctly rounded maths in both codes](correctly-rounded-maths-in-both-codes.md) | Why neither code takes its transcendentals from the host, and the host differences that forced it |
| [Instrumentation patches and logging](instrumentation-patches-and-logging.md) | What the instrumented build dumps, in what format, and how to add a dump |
| [Proofs that the reference is sound](proofs-that-the-reference-is-sound.md) | That the instrumentation changes nothing, that nothing reads an uninitialised value, that no maths comes from the host, and that another host reproduces it |

Next in this tab: the runs of this build that the tests read, [fixtures](../fixtures.md).
