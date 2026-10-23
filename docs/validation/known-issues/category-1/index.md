---
icon: lucide/file-output
---

# 1. Output

Where XFOIL's printed and plotted output differs from the solver's state, or is too coarse to compare against.

| | Issue | Description | Handling |
|---|---|---|---|
| [1.1](1-1-dump-prints-one-step-stale-closure-quantities-h-cf-k-tau-di.md) | [DUMP prints one-step-stale closure quantities (H*, Cf, K, tau, Di)](1-1-dump-prints-one-step-stale-closure-quantities-h-cf-k-tau-di.md) | DUMP prints closure quantities evaluated one Newton step behind the final state. | Overcome (output only) |
| [1.2](1-2-dump-s-ue-vinf-is-signed-by-gam.md) | [DUMP's `Ue/Vinf` is signed by GAM](1-2-dump-s-ue-vinf-is-signed-by-gam.md) | DUMP signs the edge velocity by the vortex strength, so the lower surface prints a negative Ue. | Overcome (output only) |
| [1.3](1-3-dump-writes-h-h-1-where-thet-0.md) | [DUMP writes `H = H* = 1` where THET = 0](1-3-dump-writes-h-h-1-where-thet-0.md) | Where the momentum thickness is zero, DUMP writes H and H* as 1. | Replicated |
| [1.4](1-4-cpdisp-s-displacement-surface-has-no-scale-factor-and-splits-the-wake-by-te-values.md) | [CPDISP's displacement surface has no scale factor and splits the wake δ* by TE values](1-4-cpdisp-s-displacement-surface-has-no-scale-factor-and-splits-the-wake-by-te-values.md) | The plotted displacement surface is drawn at true scale, and the wake δ* is split using trailing-edge values. | Replicated, extended |
| [1.5](1-5-formatted-output-precision.md) | [Formatted output precision](1-5-formatted-output-precision.md) | XFOIL's output formats carry too few digits to support the project's tolerances. | Overcome (instrumentation) |
