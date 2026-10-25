# NACA 64A010, TYPE 2

Shows: [§7.13 Whether a sweep converges, and what it converges to, can depend on the panel count](../category-7/7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md); [§4 Silent recovery](../category-4/index.md).

--8<-- "validation/known-issues/tables/cases-naca-64a010-type-2.md"

!!! note "Stub — to be written"

    What this run shows, and where and why yFoil and XFOIL part.

## N = 60

`naca64a010_n60_polar30_re1e6_type2` — Known issues 7.13 at N = 60: NACA 64A010, viscous TYPE 2, M 0 (MRCL 100xRe clamp at CL ~ 0; 5 of 14 points converge).

<ironlab-figure src="../figures/polars-naca64a010_n60_polar30_re1e6_type2.fig" format="fig" name="polars-naca64a010_n60_polar30_re1e6_type2" alt="naca64a010_n60_polar30_re1e6_type2: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n60_polar30_re1e6_type2: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

## N = 160

`naca64a010_n160_polar30_re1e6_type2_iter100` — Known issues 7.13 at N = 160: NACA 64A010, viscous TYPE 2, M 0 (on CORE-MATH 12-15 deg fail and the downward leg halts at -4 deg, no hang; on the host-libm reference (before 2026-10-15) with the panels of then, 11 of 13 converge, CD is infinite at 12 deg and it hangs, 4: pathological_naca64a010_n160_polar30_re1e6_type2_iter100_host_libm).

<ironlab-figure src="../figures/polars-naca64a010_n160_polar30_re1e6_type2_iter100.fig" format="fig" name="polars-naca64a010_n160_polar30_re1e6_type2_iter100" alt="naca64a010_n160_polar30_re1e6_type2_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n160_polar30_re1e6_type2_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

## N = 240

`naca64a010_n240_polar30_re1e6_type2_iter200` — Known issues 7.13 at N = 240: NACA 64A010, viscous TYPE 2, M 0 (no point converges).

<ironlab-figure src="../figures/polars-naca64a010_n240_polar30_re1e6_type2_iter200.fig" format="fig" name="polars-naca64a010_n240_polar30_re1e6_type2_iter200" alt="naca64a010_n240_polar30_re1e6_type2_iter200: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n240_polar30_re1e6_type2_iter200: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

## Reproducing the hang (N = 160)

On the host-libm reference, the reference used before 2026-10-15, this sweep hangs in XFOIL's plot
label. On the CORE-MATH reference it does not. The bundle
`xtask/fixtures-config/reproduce/pathological_naca64a010_n160_polar30_re1e6_type2_iter100_host_libm/` holds everything needed to repeat the hang without yFoil:
the `panels.dat` XFOIL loads (the panels yFoil generated before its generators moved to CORE-MATH; today's panels do not hang on either reference), the `xfoil.inp` that drives it, and the `panels.json`
that the fixture pipeline checks the geometry handoff against.

!!! info "Verified 2026-10-25"

    yFoil 0.1.0 <!-- commit: fill in after merge -->, host-libm reference `inputs_fingerprint`
    `a72447693d15adef02ff948d7ee42dc1d4abba50835c9e9e9b05764e3e0846b0`
    (`target/xfoil-ref/host-libm/manifest.txt`), GNU Fortran 14.2.0 (Homebrew GCC), on macOS 15.6.1
    (24G90), arm64, Apple `libSystem.B.dylib` 1345.120.2. Run from the committed bundle alone; the
    pristine and the instrumented host-libm builds stop at the same point.

From the repository root (macOS; on Linux, install gfortran and the X11 development library instead):

```bash
brew install gcc libx11                 # gfortran, and the X11 library XFOIL's plot library links
scripts/xfoil-build.sh --host-libm      # target/xfoil-ref/host-libm/{pristine,instrumented}/
mkdir -p .tmp/hang && cp xtask/fixtures-config/reproduce/pathological_naca64a010_n160_polar30_re1e6_type2_iter100_host_libm/{panels.dat,xfoil.inp} .tmp/hang/
cd .tmp/hang && GFORTRAN_UNBUFFERED_ALL=y ../../target/xfoil-ref/host-libm/pristine/bin/xfoil < xfoil.inp
```

(`GFORTRAN_UNBUFFERED_ALL=y` only stops gfortran buffering the output, so the last lines are
visible when the run stops.) Within a minute the output stops at alpha = 12° on the upward leg, at
point number 13 of the run (VISCAL call 13, SETBL call 464). That point runs to its
iteration limit, and its last lines are

```text
       a = 12.000      CL =  2.7442
      Cm = -0.5457     CD = Infinity   =>   CDf =  0.00245    CDp = Infinity
 VISCAL:  Convergence failed
```

then one line holding `Infinity`. XFOIL never returns: it stays at 100 % CPU in the label's
digit loop, so stop it with Ctrl-C. That CD is `Infinity` is what marks this as the plot-label
hang. Alternatively, with yFoil's tooling, run
`cargo xtask fixtures --case pathological_naca64a010_n160_polar30_re1e6_type2_iter100_host_libm`. It runs the instrumented build on the same bundle, the
watchdog ends the run after 120 s without output, and
`target/fixtures/pathological_naca64a010_n160_polar30_re1e6_type2_iter100_host_libm/manifest.json` then carries `"truncated_by_watchdog": true`. The last record
of `viscal_iters_all.dat` in that directory has CD = `Infinity`.

!!! warning "Depends on the platform's maths library"

    XFOIL's own plot-label loop on a non-finite CD causes the hang (§4.17), not the maths library.
    What the library decides is whether the run gets there. The run is ill-conditioned, so one-ULP
    differences between maths libraries put it on another path. On another OS or libm (glibc, or
    the CORE-MATH reference yFoil is validated against), the same case may not hang, or may hang at
    another operating point. Three known-issues sweeps show the same thing on the CORE-MATH
    reference, where they hang at other angles: `naca0012_n240_polar30_re1e6_iter200` (24° there,
    22° on the host libm), `naca4412_n240_polar30_re1e6_m06_iter200` (21°) and
    `naca64a010_n160_polar30_re1e6_m03_type2_iter100` (13°).
