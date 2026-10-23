---
icon: lucide/database-zap
---

# 2. Persistent state

Variables XFOIL leaves unset, stale or carried over between calls, which a translation has to reproduce.

| | Issue | Description | Handling |
|---|---|---|---|
| [2.1](2-1-hvrat-is-never-assigned-on-the-analysis-path.md) | [HVRAT is never assigned on the analysis path](2-1-hvrat-is-never-assigned-on-the-analysis-path.md) | The viscosity-law constant HVRAT is set only by plotting routines, so the analysis runs with zero. | Replicated |
| [2.2](2-2-mrchdu-locals-persist-across-sides.md) | [MRCHDU locals persist across sides](2-2-mrchdu-locals-persist-across-sides.md) | Several MRCHDU variables carry over from side 1 into side 2. | Replicated |
| [2.3](2-3-trdif-s-transition-point-station-reuses-station-2-common.md) | [TRDIF's transition-point station reuses station-2 COMMON](2-3-trdif-s-transition-point-station-reuses-station-2-common.md) | At the transition point TRDIF overwrites station 2's values without recomputing their derivatives. | Replicated |
| [2.4](2-4-blvar-clamps-hk2-in-common-without-touching-its-derivatives.md) | [BLVAR clamps HK2 in COMMON without touching its derivatives](2-4-blvar-clamps-hk2-in-common-without-touching-its-derivatives.md) | BLVAR clamps HK2 but leaves its sensitivities as they were. | Replicated |
| [2.5](2-5-trchek2-leaves-ampl2-xt-and-xt-as-the-loop-left-them.md) | [TRCHEK2 leaves AMPL2, XT and XT_* as the loop left them](2-5-trchek2-leaves-ampl2-xt-and-xt-as-the-loop-left-them.md) | When its Newton iteration fails, TRCHEK2 leaves the transition quantities where the loop left them. | Replicated |
| [2.6](2-6-vsrez-4-is-stale-at-mrchue-s-inner-loop-entry.md) | [VSREZ(4) is stale at MRCHUE's inner-loop entry](2-6-vsrez-4-is-stale-at-mrchue-s-inner-loop-entry.md) | One of MRCHUE's residuals still holds the previous iteration's value when its inner loop starts. | Replicated (instrumentation note) |
| [2.7](2-7-radle-stays-zero-for-a-flat-le.md) | [RADLE stays zero for a flat LE](2-7-radle-stays-zero-for-a-flat-le.md) | A leading edge flatter than a threshold leaves the leading-edge radius at zero. | Replicated |
| [2.8](2-8-mass-is-never-written-on-side-1-s-wake-slots.md) | [MASS is never written on side 1's wake slots](2-8-mass-is-never-written-on-side-1-s-wake-slots.md) | MASS is never set on side 1's wake stations. | Replicated (harness note) |
| [2.9](2-9-the-closure-log-s-s2-at-laminar-blvar-calls-is-not-a-value-blvar-uses.md) | [The closure log's `s2` at laminar BLVAR calls is not a value BLVAR uses](2-9-the-closure-log-s-s2-at-laminar-blvar-calls-is-not-a-value-blvar-uses.md) | A value in the reference's closure log at laminar BLVAR calls is not one BLVAR uses. | Out of scope |
