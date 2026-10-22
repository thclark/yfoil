# Correctly rounded maths in both codes

Both codes take every transcendental — `**`, `EXP`, `LOG`, `LOG10`, `TANH`, `SIN`, `COS`, `ATAN`,
`ATAN2`, `ASIN` — from CORE-MATH's correctly rounded functions, compiled from one vendored copy
(`third-party/core-math`, upstream commit in `core-math.ORIGIN`) with the same flags. The
reference links them in place of the host libm through `scripts/xfoil-build/core_math_libm.c`;
yFoil calls them through the `cr-libm` crate, and `clippy.toml` disallows `std`'s `f64` methods,
which call the host libm. A correctly rounded result is a property of its argument alone, so it is
the same on every host, and the same in both codes: a fixture holds on every host at the same
tolerance, and there is no cross-host allowance anywhere. A function not there yet is added to
`cr-libm` and to the reference's shim together.

## Why the host's maths had to go

Until 2026-10-15 both codes used the host's libm, and the September investigation of a CI failure
on Linux (every test passed on the arm64 macOS fixture host) found that this was the whole story.
It was neither over-tight tolerances nor a translation difference: Rust's `f64` methods call the
platform C library exactly as gfortran's intrinsics do, and Apple's libSystem and glibc disagree by
one ULP on some inputs. Everywhere well-conditioned that was absorbed by the tolerances; the four
failing tests were the [polar break points](../execution-equivalence/polar-break-points.md), which
probe the region where the reference cannot reproduce itself under a 1-ULP perturbation.

Measured in an Ubuntu 24.04 container (glibc 2.39):

| Measurement | Result |
|---|---|
| Probe, 20 000 inputs per function | libSystem 1345.120.2 vs glibc 2.39 differ by 1 ULP on 0.1 % (`exp`, `ln`, `powf`) to 18 % (`tanh`) of inputs; `tanh` occasionally 2 ULP; `sqrt` identical (IEEE correctly rounded) |
| arm64 glibc container vs x86_64 GitHub runner | failing numbers bit-for-bit identical, so it is the libm, not the CPU or the compiler |
| yFoil's cosine-spaced NACA panels, macOS vs glibc | 4 of 160 nodes differ, by up to 13 ULP |
| XFOIL regenerated on glibc, converged polar points | CL moves by at most 1.4e-11 (NACA 0012) and 7.5e-11 (NACA 4412) |
| XFOIL regenerated on glibc, unconverged post-CL_max wanderings | CL moves by O(1); 22° (0012) and −16° (4412) no longer converge at all |
| One-step replay from identical dumped XFOIL state, yFoil on macOS vs glibc | 3.6e-10 and 1.2e-10 relative in the largest Newton delta; same-host the replay reproduces XFOIL to at most 2e-12 |

So the pins in those tests — the iteration at which the runs part, and a chance convergence at the
fourth attempt — were properties of the macOS libm, not of the translation. For a month the tests
carried a cross-host tolerance and a host check; both were removed when the maths was made the
same everywhere.

## What it took

| Finding | Consequence |
|---|---|
| Linking CORE-MATH in place of `exp`, `log`, `pow`… was not enough. GCC fuses a `SIN`/`COS` pair of one argument: into `sincos` on glibc, and into `cexp(i·x)` on Darwin, which has no `sincos` | Both are shimmed. The Darwin `cexp` leak alone moved one inviscid TYPE 2 case (`naca64a010_n60_inviscid_type2_m03`, the `CL 8` point) by 1 ULP in `SIN(ALFA)` inside CLCALC, located by forward-stepping a probe through SPECCL on both hosts |
| A shim can fail silently | Every build now fails if any maths function is left to the host libm ([proof of no host maths](proofs-that-the-reference-is-sound.md#proof-of-no-host-maths)) |
| The cross-host container reused a reference binary built a week earlier, and `cargo xtask fixtures` reused any existing binary | Each build stage records a fingerprint of everything it depends on, and every xtask that runs the reference rebuilds a stage whose fingerprint is not current ([the reference build](index.md#layout-and-how-it-is-built)) |
| CORE-MATH compiled by Apple clang and by GCC | Identical hashes of 5 million outputs per function, all 11 functions; `cr_sincos` equals `cr_sin`/`cr_cos` on 20 million inputs |
| The smoke case | `cp.dat`, `bl.dat` and the console output byte-identical between arm64 macOS and aarch64 glibc |

Changing the libm moved the reference's own trajectories wherever they are ill-conditioned: the
non-finite probes' rarely taken branches now fall at other SETBL calls (`UPDATE_RDN4` in the 64A010
TYPE 2 M 0.3 case moved from call 7 to call 17), which is why the fixtures, the step cover and the
coverage measurement were regenerated together. That the reference now reproduces itself on
another host is checked by [reproduction on another host](proofs-that-the-reference-is-sound.md#reproduction-on-another-host).
