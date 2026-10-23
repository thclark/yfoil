# Known issues

--8<-- "validation/known-issues/tables/run.md"

## What this study is for

XFOIL has known issues: places where its answer depends on the panel count, where it continues after
a failure, where it hangs, and where a mode drives it somewhere it cannot recover from. They are all
listed in [Known XFOIL issues](../../user-guide/reference/known-xfoil-issues.md). This study runs cases
that show those issues arising, and checks whether yFoil does the same thing. yFoil is a translation,
so it should reproduce XFOIL's issues rather than fix them. The study is to grow until it shows every
documented issue.

Each run is compared exactly as the [series cases](../series-cases/index.md#how-a-run-is-compared) are,
and its figure reads the same way ([reading the figures](../series-cases/index.md#reading-the-figures)).
Unlike the series cases, these runs sit where XFOIL is often ill-conditioned or non-finite, so partings
are expected. What matters is whether XFOIL's own twins part there too.

## Main results

!!! note "Stub — to be written"

    The headline per issue: which runs yFoil follows to the end, and where the partings are XFOIL's
    own ill-conditioning.

--8<-- "validation/known-issues/tables/cases.md"

## Effect of accumulated error on run parting

A run [parts](../../development/conventions/terminology.md#run-parting) where yFoil stops matching
XFOIL, either because a branch decision differs or because a value leaves tolerance; *parted by* in
the table says which. These runs have no 1-ULP twins, so the floor column is zero and each ratio is against the named
tolerance alone; the table cannot say whether XFOIL is ill-conditioned at a point. It shows each point
where yFoil parts and the four points of the leg before it, with the worst array of yFoil's final state
against XFOIL's, the difference there, and the difference as a ratio of its tolerance. On some runs the
ratio is above 1 at points whose iterations all agree with XFOIL's, and on the closed-trailing-edge
NACA 63-415 it grows from point to point along the downward leg.

!!! note "Stub — to be written"

    What this means for the known-issues runs.

--8<-- "validation/known-issues/tables/final-state.md"

## The issues and the runs that show them

Every issue in [Known XFOIL issues](../../user-guide/reference/known-xfoil-issues.md), with the runs of this study that show it. Categories 3, 4, 8 and 9 are registers, a numbered table or list on one page, so each is one row.

| | Issue | Handling | Runs |
|---|---|---|---|
| [1.1](category-1/1-1-dump-prints-one-step-stale-closure-quantities-h-cf-k-tau-di.md) | [DUMP prints one-step-stale closure quantities (H*, Cf, K, tau, Di)](category-1/1-1-dump-prints-one-step-stale-closure-quantities-h-cf-k-tau-di.md) | Overcome (output only) | — |
| [1.2](category-1/1-2-dump-s-ue-vinf-is-signed-by-gam.md) | [DUMP's `Ue/Vinf` is signed by GAM](category-1/1-2-dump-s-ue-vinf-is-signed-by-gam.md) | Overcome (output only) | — |
| [1.3](category-1/1-3-dump-writes-h-h-1-where-thet-0.md) | [DUMP writes `H = H* = 1` where THET = 0](category-1/1-3-dump-writes-h-h-1-where-thet-0.md) | Replicated | — |
| [1.4](category-1/1-4-cpdisp-s-displacement-surface-has-no-scale-factor-and-splits-the-wake-by-te-values.md) | [CPDISP's displacement surface has no scale factor and splits the wake δ* by TE values](category-1/1-4-cpdisp-s-displacement-surface-has-no-scale-factor-and-splits-the-wake-by-te-values.md) | Replicated, extended | — |
| [1.5](category-1/1-5-formatted-output-precision.md) | [Formatted output precision](category-1/1-5-formatted-output-precision.md) | Overcome (instrumentation) | — |
| [2.1](category-2/2-1-hvrat-is-never-assigned-on-the-analysis-path.md) | [HVRAT is never assigned on the analysis path](category-2/2-1-hvrat-is-never-assigned-on-the-analysis-path.md) | Replicated | — |
| [2.2](category-2/2-2-mrchdu-locals-persist-across-sides.md) | [MRCHDU locals persist across sides](category-2/2-2-mrchdu-locals-persist-across-sides.md) | Replicated | — |
| [2.3](category-2/2-3-trdif-s-transition-point-station-reuses-station-2-common.md) | [TRDIF's transition-point station reuses station-2 COMMON](category-2/2-3-trdif-s-transition-point-station-reuses-station-2-common.md) | Replicated | — |
| [2.4](category-2/2-4-blvar-clamps-hk2-in-common-without-touching-its-derivatives.md) | [BLVAR clamps HK2 in COMMON without touching its derivatives](category-2/2-4-blvar-clamps-hk2-in-common-without-touching-its-derivatives.md) | Replicated | — |
| [2.5](category-2/2-5-trchek2-leaves-ampl2-xt-and-xt-as-the-loop-left-them.md) | [TRCHEK2 leaves AMPL2, XT and XT_* as the loop left them](category-2/2-5-trchek2-leaves-ampl2-xt-and-xt-as-the-loop-left-them.md) | Replicated | — |
| [2.6](category-2/2-6-vsrez-4-is-stale-at-mrchue-s-inner-loop-entry.md) | [VSREZ(4) is stale at MRCHUE's inner-loop entry](category-2/2-6-vsrez-4-is-stale-at-mrchue-s-inner-loop-entry.md) | Replicated (instrumentation note) | — |
| [2.7](category-2/2-7-radle-stays-zero-for-a-flat-le.md) | [RADLE stays zero for a flat LE](category-2/2-7-radle-stays-zero-for-a-flat-le.md) | Replicated | — |
| [2.8](category-2/2-8-mass-is-never-written-on-side-1-s-wake-slots.md) | [MASS is never written on side 1's wake slots](category-2/2-8-mass-is-never-written-on-side-1-s-wake-slots.md) | Replicated (harness note) | — |
| [2.9](category-2/2-9-the-closure-log-s-s2-at-laminar-blvar-calls-is-not-a-value-blvar-uses.md) | [The closure log's `s2` at laminar BLVAR calls is not a value BLVAR uses](category-2/2-9-the-closure-log-s-s2-at-laminar-blvar-calls-is-not-a-value-blvar-uses.md) | Out of scope | — |
| [3](category-3/index.md) | [Drela's fixes](category-3/index.md) | see the category | — |
| [4](category-4/index.md) | [Silent recovery](category-4/index.md) | see the category | [NACA 16-212 at M 0.7](runs/naca-16-212-at-m-0-7.md), [NACA 64A010, TYPE 2](runs/naca-64a010-type-2.md), [NACA 4412 inviscid at M 0.7](runs/naca-4412-inviscid-at-m-0-7.md), [NACA 0012, Re 1e6, N = 240](runs/naca-0012-re-1e6-n-240.md), [NACA 4412, Re 1e6, N = 240](runs/naca-4412-re-1e6-n-240.md), [NACA 4412, Re 3e6, N = 240](runs/naca-4412-re-3e6-n-240.md), [NACA 4412 at M 0.6, N = 240](runs/naca-4412-at-m-0-6-n-240.md), [NACA 63-415 at M 0.5](runs/naca-63-415-at-m-0-5.md), [NACA 64A010, TYPE 2 at M 0.3](runs/naca-64a010-type-2-at-m-0-3.md), [NACA 23012 with a trip upstream of stagnation](runs/naca-23012-with-a-trip-upstream-of-stagnation.md) |
| [5.1](category-5/5-1-subroutines-never-called.md) | [Subroutines never called](category-5/5-1-subroutines-never-called.md) | Not translated | — |
| [5.2](category-5/5-2-limage-is-permanently-false.md) | [LIMAGE is permanently false](category-5/5-2-limage-is-permanently-false.md) | Out of scope | — |
| [5.3](category-5/5-3-stale-duplicate-source-files.md) | [Stale duplicate source files](category-5/5-3-stale-duplicate-source-files.md) | Out of scope | — |
| [5.4](category-5/5-4-matyp-3-is-unreachable.md) | [MATYP = 3 is unreachable](category-5/5-4-matyp-3-is-unreachable.md) | Replicated | — |
| [5.5](category-5/5-5-oper-damp-is-reachable-and-undocumented.md) | [OPER `DAMP` is reachable and undocumented](category-5/5-5-oper-damp-is-reachable-and-undocumented.md) | Replicated | — |
| [5.6](category-5/5-6-structurally-dead-branches.md) | [Structurally dead branches](category-5/5-6-structurally-dead-branches.md) | Out of scope | — |
| [6.1](category-6/6-1-naca4-applies-thickness-vertically.md) | [NACA4 applies thickness vertically](category-6/6-1-naca4-applies-thickness-vertically.md) | Overcome (the one live divergence) | — |
| [6.2](category-6/6-2-load-reverses-clockwise-input.md) | [LOAD reverses clockwise input](category-6/6-2-load-reverses-clockwise-input.md) | Constraint | — |
| [6.3](category-6/6-3-chord-and-apanel-are-not-what-their-names-suggest.md) | [CHORD and APANEL are not what their names suggest](category-6/6-3-chord-and-apanel-are-not-what-their-names-suggest.md) | Replicated | — |
| [6.4](category-6/6-4-nw-n-12-10int-waklen.md) | [NW = N/12 + 10·INT(WAKLEN)](category-6/6-4-nw-n-12-10int-waklen.md) | Replicated | — |
| [6.5](category-6/6-5-tgap-on-a-sharp-trailing-edge-does-not-deliver-the-requested-gap.md) | [TGAP on a sharp trailing edge does not deliver the requested gap](category-6/6-5-tgap-on-a-sharp-trailing-edge-does-not-deliver-the-requested-gap.md) | Replicated | — |
| [6.6](category-6/6-6-naca-generator-constants-the-ordinate-program-s-precision.md) | [NACA generator constants: the ordinate program's precision](category-6/6-6-naca-generator-constants-the-ordinate-program-s-precision.md) | Reference note | — |
| [7.1](category-7/7-1-setexp-s-1e-5-newton-tolerance.md) | [SETEXP's 1e-5 Newton tolerance](category-7/7-1-setexp-s-1e-5-newton-tolerance.md) | Replicated | — |
| [7.2](category-7/7-2-blsolv-s-vacc2-vacc3-association.md) | [BLSOLV's VACC2/VACC3 association](category-7/7-2-blsolv-s-vacc2-vacc3-association.md) | Replicated | — |
| [7.3](category-7/7-3-karmantsien-tk-must-be-formed-as-comset-forms-it.md) | [Kármán–Tsien TK must be formed as COMSET forms it](category-7/7-3-karmantsien-tk-must-be-formed-as-comset-forms-it.md) | Overcome (yFoil bug, fixed) | — |
| [7.4](category-7/7-4-integer-powers.md) | [Integer powers](category-7/7-4-integer-powers.md) | Replicated (gfortran) | — |
| [7.5](category-7/7-5-equivalence-aliasing-of-unew-qnew-onto-va-vb.md) | [EQUIVALENCE aliasing of UNEW/QNEW onto VA/VB](category-7/7-5-equivalence-aliasing-of-unew-qnew-onto-va-vb.md) | Overcome (structure only) | — |
| [7.6](category-7/7-6-single-precision-and-fpe-checks.md) | [Single precision and FPE checks](category-7/7-6-single-precision-and-fpe-checks.md) | Reference-build facts | — |
| [7.7](category-7/7-7-type-2-on-a-section-that-can-reach-negative-cl-drives-itself-onto-the-karmantsien-pole.md) | [TYPE 2 on a section that can reach negative CL drives itself onto the Kármán–Tsien pole](category-7/7-7-type-2-on-a-section-that-can-reach-negative-cl-drives-itself-onto-the-karmantsien-pole.md) | Replicated | [NACA 64A010 inviscid, TYPE 2 at M 0.3](runs/naca-64a010-inviscid-type-2-at-m-0-3.md), [NACA 64A010, TYPE 2 at M 0.3](runs/naca-64a010-type-2-at-m-0-3.md) |
| [7.8](category-7/7-8-update-s-largest-change-report-is-a-console-label-that-rounding-can-decide.md) | [UPDATE's largest-change report is a console label that rounding can decide](category-7/7-8-update-s-largest-change-report-is-a-console-label-that-rounding-can-decide.md) | Replicated | — |
| [7.9](category-7/7-9-mrchdu-s-station-newton-at-a-transition-station-can-amplify-its-seed-a-thousandfold.md) | [MRCHDU's station Newton at a transition station can amplify its seed a thousandfold](category-7/7-9-mrchdu-s-station-newton-at-a-transition-station-can-amplify-its-seed-a-thousandfold.md) | Replicated | — |
| [7.10](category-7/7-10-the-fallback-tests-are-written-so-that-a-nan-residual-takes-the-fallback.md) | [The fallback tests are written so that a NaN residual takes the fallback](category-7/7-10-the-fallback-tests-are-written-so-that-a-nan-residual-takes-the-fallback.md) | Replicated | — |
| [7.11](category-7/7-11-mrchue-s-failure-path-calls-blvar-2-and-then-blvar-3-at-a-wake-station.md) | [MRCHUE's failure path calls BLVAR(2) and then BLVAR(3) at a wake station](category-7/7-11-mrchue-s-failure-path-calls-blvar-2-and-then-blvar-3-at-a-wake-station.md) | Replicated | — |
| [7.12](category-7/7-12-ami-ampl2-after-a-forced-transition.md) | [AMI = AMPL2 after a forced transition](category-7/7-12-ami-ampl2-after-a-forced-transition.md) | Overcome (yFoil bug, fixed) | — |
| [7.13](category-7/7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md) | [Whether a sweep converges, and what it converges to, can depend on the panel count](category-7/7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md) | Replicated | [NACA 16-212 at M 0.7](runs/naca-16-212-at-m-0-7.md), [NACA 64A010, TYPE 2](runs/naca-64a010-type-2.md), [NACA 4412 inviscid at M 0.7](runs/naca-4412-inviscid-at-m-0-7.md), [NACA 64A010 inviscid, TYPE 2 at M 0.3](runs/naca-64a010-inviscid-type-2-at-m-0-3.md) |
| [8](category-8/index.md) | [Open items](category-8/index.md) | see the category | [NACA 63-415, closed trailing edge](runs/naca-63-415-closed-trailing-edge.md), [NACA 64A010, closed trailing edge](runs/naca-64a010-closed-trailing-edge.md) |
| [9](category-9/index.md) | [Pathological cases](category-9/index.md) | see the category | [NACA 23012 with a trip upstream of stagnation](runs/naca-23012-with-a-trip-upstream-of-stagnation.md) |
