# Cross-platform quirks

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
(`tests/xfoil_polar_break_tests.rs`), which deliberately probe the region where the reference cannot
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

- `tests/utilities/host.rs` compares each fixture's recorded host (`manifest.json`, `host:` line) with
  the running one. Same-host pins are asserted on the fixture's host and reported on any other, never
  skipped (CLAUDE.md Rule 1's third outcome). Gated on every host: the converged calls before the break,
  the straddle classification by the reference's own floor, the replay within `TOL_CROSS_HOST`
  (`tests/utilities/tolerances.rs`, 1e-9 from the measured 3.6e-10), and the polar driver's bookkeeping
  against the sweep's own outcomes.
- The straddle classifier in `tests/utilities/records.rs` had a real gap: a run parting in the relaxation
  factor (RLX) at a hypersensitive iteration was failed rather than classified, where one parting in RMSBL
  was classified. Fixed.
- The replay (`tests/utilities/replay.rs`) asserted UPDATE's reported limiter letter (VMXBL) exactly; at
  the similarity station it is a rounding-decided tie (dn2 == dn3), now reported as the sweep comparison
  already did.
- `scripts/xfoil-build.sh` used BSD-only `sed -i ''`, which is why the nightly `xfoil-parity` job died on
  GNU sed at the snan stage. Fixed with `sed -i.bak … && rm`.
- CLAUDE.md, `docs/validation/README.md` (*Polar break points*) and `docs/xfoil-known-issues.md` §7.4
  record the measurement.

## Still open

- The nightly `xfoil-parity` job's `cargo xtask fixtures --verify` asserts byte-identity with the tracked
  fixtures, which is same-host by design, so it will not pass on a Linux runner against macOS-generated
  fixtures. Making that step an ULP-budget comparison is a separate change.
- The 4412 twin run regenerated on glibc hung at exit in XFOIL's sequence-plot label loop on a non-finite
  CL (known-issues §4), even with graphics off; all 34 calls were on disk. `cargo xtask fixtures` treats
  the killed run as a failure.
- Whether glibc's results are stable across its own versions (2.36 vs 2.39) was not measured; a same-OS
  libm drift would surface as a `--verify` byte-identity failure on the nightly job.

## Reproducing

```bash
# image: ubuntu:24.04 + build-essential gfortran libx11-dev libfontconfig1-dev pkg-config rustup
docker build -t yfoil-ci:ubuntu24 .tmp/linux
rsync -a --delete --exclude target --exclude .tmp --exclude site ./ .tmp/linux/yfoil/   # never bind-mount the live checkout: `cargo xtask fixtures` overwrites tracked fixture directories
docker run --rm -v "$PWD/.tmp/linux:/host" yfoil-ci:ubuntu24 /host/pipeline.sh
```
