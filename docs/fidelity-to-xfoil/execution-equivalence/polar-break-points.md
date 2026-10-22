# Polar break points

Past maximum lift a sweep stops converging, and where exactly it breaks — which attempt converges,
on which branch — is decided at the reference's noise floor. These sweeps are where a polar that
"carries on" in one code and halts in the other looks most like a translation bug and is not one.
This page records the two tracked cases, how each is tested, and why.

## The two break-point cases

Two tracked cases drive a 0.5° `ASEQ` sweep at `ITER 100` past CL_max, where the Newton iteration
stops converging and both codes wander for the full 105 iterations per point. The solutions there
are ill-conditioned ([terminology](../../contributing/conventions/terminology.md#ill-conditioned-solution)), and how each case is
tested is recorded beside it in `xtask/fixtures-config/cases.toml`: the calls before the break as a whole run
([whole-run comparison](whole-run-comparison.md)), the break call one step at a time from XFOIL's
dumped state ([step-by-step replay](step-by-step-replay.md)):

| Case (`tests/fixtures/xfoil/…`) | OPER script | Reference behaviour | Tested as |
|---|---|---|---|
| `naca0012_n160_polar_up22_re1e6_iter100` | `ALFA 0 / ASEQ 0.5 22 0.5` | 20.5°–22° unconverged (on the host-libm reference, before 2026-10-15, 22° converged on the separated branch in 7 iterations) | run through call 41; every iteration of call 42 as a step |
| `naca4412_n160_polar_down16_re1e6_iter100` | `ALFA 0 / INIT / ALFA 0 / ASEQ -0.5 -16 -0.5` | −14.5° converges on the separated branch (CL ≈ +0.10) in 53 iterations, its 1-ULP twins do not (four wander for 105, one converges elsewhere); −15° to −16° follow on that branch (CL ≈ +0.05 at −16°). On the host-libm reference −14.5° to −15.5° were unconverged and −16° converged in 39 | run through call 30; every iteration of call 31 as a step |

## Why the break is not a translation difference

In both cases the reference's own 1-ULP twins part from the reference, and flip its stagnation and
transition stations, before yFoil does, and every step replayed from XFOIL's dumped state at the break
reproduces the step. Which attempt converges is therefore decided at the noise floor — changing the
reference's libm moved it in both cases; a polar
that "carries on" past the break in one code and halts in the other is a divergent comparison, not a
translation bug, and no test compares the two codes' own marches through it.

The reference's behaviour at the break does not depend on the host. Both codes take their
transcendentals from CORE-MATH's correctly rounded functions, not the host libm, and the reference
reproduces every tracked record byte for byte on glibc (`scripts/cross-host.sh`), so a test holds every
value to the same tolerance everywhere. Until 2026-10-15 both took them from the host libm, Apple's
libSystem and glibc disagreeing by 1 ULP on 0.1 % (`exp`, `ln`, `pow`) to 18 % (`tanh`) of inputs, and
the break was host-dependent: on glibc the fourth attempt's chance convergence did not survive and XFOIL
failed 22° and −16° as well ([correctly rounded maths in both codes](../the-reference-xfoil-build/correctly-rounded-maths-in-both-codes.md)).
The same cases wake the MRCHUE inverse wake march, the MRCHDU iterate fallback, BLVAR's Us and Hk
clamps and TRCHEK2's iteration cap that every other case left open
([branches the reference takes](../branch-coverage/branches-the-reference-takes.md)). Two XFOIL quirks
surfaced here are registered in [the issue register](../xfoils-known-issues/the-issue-register.md) (§2.8, §4): MASS is never
written on side 1's wake slots, and ASEQ's sequence-plot label loops forever on a non-finite CL/CM even
with graphics off, so a reference sweep that blows up past the break hangs instead of halting.
