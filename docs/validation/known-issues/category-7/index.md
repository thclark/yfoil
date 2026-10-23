---
icon: lucide/binary
---

# 7. Numerics

Numerical details a reproduction has to copy exactly, and behaviours they produce.

| | Issue | Description | Handling |
|---|---|---|---|
| [7.1](7-1-setexp-s-1e-5-newton-tolerance.md) | [SETEXP's 1e-5 Newton tolerance](7-1-setexp-s-1e-5-newton-tolerance.md) | SETEXP's loose Newton tolerance shapes every wake node, so it has to be kept. | Replicated |
| [7.2](7-2-blsolv-s-vacc2-vacc3-association.md) | [BLSOLV's VACC2/VACC3 association](7-2-blsolv-s-vacc2-vacc3-association.md) | The thresholds behind BLSOLV's sparsity skips have to be formed exactly as XFOIL forms them. | Replicated |
| [7.3](7-3-karmantsien-tk-must-be-formed-as-comset-forms-it.md) | [Kármán–Tsien TK must be formed as COMSET forms it](7-3-karmantsien-tk-must-be-formed-as-comset-forms-it.md) | The Kármán–Tsien factor has to be computed in COMSET's form, which differs in the last bits from an algebraically equal one. | Overcome (yFoil bug, fixed) |
| [7.4](7-4-integer-powers.md) | [Integer powers](7-4-integer-powers.md) | Integer powers have to be multiplied out, as gfortran does, to match bit for bit. | Replicated (gfortran) |
| [7.5](7-5-equivalence-aliasing-of-unew-qnew-onto-va-vb.md) | [EQUIVALENCE aliasing of UNEW/QNEW onto VA/VB](7-5-equivalence-aliasing-of-unew-qnew-onto-va-vb.md) | UPDATE's arrays share storage with BLSOLV's through EQUIVALENCE. | Overcome (structure only) |
| [7.6](7-6-single-precision-and-fpe-checks.md) | [Single precision and FPE checks](7-6-single-precision-and-fpe-checks.md) | Facts about the reference build: its precision and floating-point checks. | Reference-build facts |
| [7.7](7-7-type-2-on-a-section-that-can-reach-negative-cl-drives-itself-onto-the-karmantsien-pole.md) | [TYPE 2 on a section that can reach negative CL drives itself onto the Kármán–Tsien pole](7-7-type-2-on-a-section-that-can-reach-negative-cl-drives-itself-onto-the-karmantsien-pole.md) | On a section that can reach negative CL, TYPE 2 drives the Mach number onto the Kármán–Tsien pole. | Replicated |
| [7.8](7-8-update-s-largest-change-report-is-a-console-label-that-rounding-can-decide.md) | [UPDATE's largest-change report is a console label that rounding can decide](7-8-update-s-largest-change-report-is-a-console-label-that-rounding-can-decide.md) | UPDATE's largest-change report is a console label that rounding can decide. | Replicated |
| [7.9](7-9-mrchdu-s-station-newton-at-a-transition-station-can-amplify-its-seed-a-thousandfold.md) | [MRCHDU's station Newton at a transition station can amplify its seed a thousandfold](7-9-mrchdu-s-station-newton-at-a-transition-station-can-amplify-its-seed-a-thousandfold.md) | At a transition station, MRCHDU's station Newton can amplify a tiny difference a thousandfold before converging. | Replicated |
| [7.10](7-10-the-fallback-tests-are-written-so-that-a-nan-residual-takes-the-fallback.md) | [The fallback tests are written so that a NaN residual takes the fallback](7-10-the-fallback-tests-are-written-so-that-a-nan-residual-takes-the-fallback.md) | A NaN residual in MRCHUE or MRCHDU takes the iterate fallback. | Replicated |
| [7.11](7-11-mrchue-s-failure-path-calls-blvar-2-and-then-blvar-3-at-a-wake-station.md) | [MRCHUE's failure path calls BLVAR(2) and then BLVAR(3) at a wake station](7-11-mrchue-s-failure-path-calls-blvar-2-and-then-blvar-3-at-a-wake-station.md) | MRCHUE's failure path at a wake station applies both its clamps in turn. | Replicated |
| [7.12](7-12-ami-ampl2-after-a-forced-transition.md) | [AMI = AMPL2 after a forced transition](7-12-ami-ampl2-after-a-forced-transition.md) | After a forced transition, the amplification seed has to be set as XFOIL sets it. | Overcome (yFoil bug, fixed) |
| [7.13](7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md) | [Whether a sweep converges, and what it converges to, can depend on the panel count](7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md) | Whether a sweep converges, and what it converges to, can change with the number of panels. | Replicated |
