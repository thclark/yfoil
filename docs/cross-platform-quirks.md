# Cross-platform quirks

Superseded by the resolution below: since 2026-10-15 neither code takes its maths from the host. The
notes above it record why that was needed.

Raw notes from the 2026-09-11 investigation of the Build and Test workflow failing on GitHub (Linux runner)
while every test passed on the fixture host (arm64 macOS). To be tidied into a coherent overview of
cross-platform issues. Branch `fix-ci-cross-host`, PR #9.

## The question

Is the failure a problem of the fixtures having too tight a tolerance (in which case the solution is to
re-run the analysis of the subroutine differences on different platforms), or of the computation in yFoil
being different (in which case: why does it depend on an OS-specific library)?

## The answer

It is neither over-tight tolerances nor a translation difference. yFoil's computation does depend on an
OS-specific library, and so does XFOIL's. Rust's `f64` methods for `exp`, `ln`, `powf`, `atan2`, `sin`,
`cos` and `tanh` call the platform C library, exactly as gfortran's intrinsics do. Apple's libSystem and
glibc disagree by exactly 1 ULP on some inputs. The fixtures were generated on a Mac and CI runs on Linux,
so on CI every transcendental in both the reference trajectory and yFoil's differs by an ULP here and
there. Everywhere well-conditioned that is absorbed by the existing tolerances, which is why every other
fixture test passes on Linux. The four failing tests are the polar break cases
(`tests/execution/polar_break.rs`), which deliberately probe the region where the reference cannot
reproduce itself under a 1-ULP perturbation, and there the ULPs get amplified.

## Evidence

All from an Ubuntu 24.04 container (glibc 2.39), Dockerfile, pipeline script and libm probe in
`.tmp/linux/` (untracked).

| Measurement | Result |
|---|---|
| Probe, 20 000 inputs per function | libSystem 1345.120.2 vs glibc 2.39 differ by 1 ULP on 0.1 % (`exp`, `ln`, `powf`) to 18 % (`tanh`) of inputs; `tanh` occasionally 2 ULP; `sqrt` identical (IEEE correctly rounded) |
| arm64 glibc container vs x86_64 GitHub runner | failing numbers bit-for-bit identical, so it is the libm, not the CPU or the compiler |
| yFoil's cosine-spaced NACA panels, macOS vs glibc | 4 of 160 nodes differ, by up to 13 ULP |
| XFOIL regenerated on glibc, converged polar points | CL moves by at most 1.4e-11 (NACA 0012) and 7.5e-11 (NACA 4412) |
| XFOIL regenerated on glibc, unconverged post-CL_max wanderings | CL moves by O(1); 22° (0012) and −16° (4412) no longer converge at all |
| One-step replay from identical dumped XFOIL state, yFoil on macOS vs glibc | 3.6e-10 and 1.2e-10 relative in the largest Newton delta (RMXBL); same-host the replay reproduces XFOIL to at most 2e-12 |

So the pins in those tests — the iteration at which the runs part (32 and 16 on macOS; 23 and 22 on
glibc), and the fourth attempt's chance convergence — were properties of the macOS libm, not of the
translation. yFoil on Linux with Linux-generated fixtures matched every converged call and passed every
replay value, failing only on those macOS-specific pins.

## What the branch changed

- `tests/common/utilities/host.rs` compares each fixture's recorded host (`manifest.json`, `host:` line) with
  the running one. Same-host pins are asserted on the fixture's host and reported on any other, never
  skipped (CLAUDE.md Rule 1's third outcome). Tested on every host: the converged calls before the break,
  the straddle classification by the reference's own floor, the replay within `TOL_CROSS_HOST`
  (`tests/common/utilities/tolerances.rs`, 1e-9 from the measured 3.6e-10), and the polar driver's bookkeeping
  against the sweep's own outcomes.
- The straddle classifier in `tests/common/utilities/records.rs` had a real gap: a run parting in the relaxation
  factor (RLX) at a hypersensitive iteration was failed rather than classified, where one parting in RMSBL
  was classified. Fixed.
- The replay (`tests/common/utilities/replay.rs`) asserted UPDATE's reported limiter letter (VMXBL) exactly; at
  the similarity station it is a rounding-decided tie (dn2 == dn3), now reported as the sweep comparison
  already did.
- `scripts/xfoil-build.sh` used BSD-only `sed -i ''`, which is why the former nightly `xfoil-parity` job died on
  GNU sed at the snan stage. Fixed with `sed -i.bak … && rm`.
- The 4412 twin run regenerated on glibc hung at exit in XFOIL's sequence-plot label loop on a non-finite
  CL (known-issues §4), even with graphics off; all 34 calls were on disk. `cargo xtask fixtures` now
  keeps such a run: its watchdog ends a reference that stops producing output and keeps what it wrote.
- CLAUDE.md, `docs/testing/README.md` (*Polar break points*) and `docs/xfoil-known-issues.md` §7.4
  record the measurement.

## Resolution: one correctly rounded libm for both codes (2026-10-15)

The host dependence is removed at its source rather than budgeted for. Both codes now take every
transcendental from CORE-MATH's correctly rounded functions, compiled from one vendored copy
(`third-party/core-math`, upstream commit in `core-math.ORIGIN`) with the same flags: the reference is
linked against them through `scripts/xfoil-build/core_math_libm.c`, yFoil calls them through the
`cr-libm` crate, and `clippy.toml` disallows `std`'s `f64` methods. A correctly rounded result is a
property of its argument alone, so it is the same on every host, and the same in both codes.

What it took, found on the way:

| Finding | Consequence |
|---|---|
| Linking CORE-MATH in place of `exp`, `log`, `pow`… was not enough. GCC fuses a `SIN`/`COS` pair of one argument: into `sincos` on glibc, and into `cexp(i·x)` on Darwin, which has no `sincos` | Both are shimmed. The Darwin `cexp` leak alone moved one inviscid TYPE 2 case (`naca64a010_n60_inviscid_type2_m03`, the `CL 8` point) by 1 ULP in `SIN(ALFA)` inside CLCALC, located by forward-stepping a probe through SPECCL on both hosts |
| A shim can fail silently | Every build now fails if `nm -u` shows any maths function left to the host libm, except MDES's `cabs`/`clog` (inverse design, out of scope) and `sqrt`, which IEEE 754 requires correctly rounded |
| The cross-host container reused a reference binary built a week earlier, and `cargo xtask fixtures` reused any existing binary | Each build stage records a fingerprint of everything it depends on; every xtask that runs the reference rebuilds a stage whose fingerprint is not current, and the fixture manifests carry it (`inputs_fingerprint`). Twins are kept on regeneration only if the reference is unchanged too |
| CORE-MATH compiled by Apple clang and by GCC | Identical hashes of 5 million outputs per function, all 11 functions; `cr_sincos` equals `cr_sin`/`cr_cos` on 20 million inputs |
| The smoke case | `cp.dat`, `bl.dat` and the console output byte-identical between arm64 macOS and aarch64 glibc |

Changing the libm moved the reference's own trajectories wherever they are ill-conditioned: the
non-finite probes' rarely taken branches now fall at other SETBL calls (`UPDATE_RDN4` in the 64A010
TYPE 2 M 0.3 case moved from call 7 to call 17), which is why the fixtures, the step cover and the
coverage measurement are regenerated together.

With the host's libm out of both codes, the cross-host tolerances (`TOL_CROSS_HOST`,
`CROSS_HOST_FACTOR`), the measured spread (`tests/fixtures/cross-host/spread.json`) and the host check
(`tests/common/utilities/host.rs`) described above were removed: every test holds its values to the same
tolerance on every host. `scripts/cross-host.sh` now checks that the reference reproduces every tracked
record byte for byte on glibc, and the `reference-build` workflow that the CI reference case regenerates
byte for byte on Linux. Whether glibc's own versions agree is no longer a question for either code.

## Reproducing

```bash
# image: ubuntu:24.04 + build-essential gfortran libx11-dev libfontconfig1-dev pkg-config rustup
docker build -t yfoil-ci:ubuntu24 .tmp/linux
rsync -a --delete --exclude target --exclude .tmp --exclude site ./ .tmp/linux/yfoil/   # never bind-mount the live checkout: `cargo xtask fixtures` overwrites tracked fixture directories
docker run --rm -v "$PWD/.tmp/linux:/host" yfoil-ci:ubuntu24 /host/pipeline.sh
```
