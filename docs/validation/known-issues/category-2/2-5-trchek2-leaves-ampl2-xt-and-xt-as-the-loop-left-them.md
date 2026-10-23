# 2.5 TRCHEK2 leaves AMPL2, XT and XT_* as the loop left them — Replicated

After 30 Newton iterations TRCHEK2 prints `N2 convergence failed.` and continues (`xblsys.f:277`);
`AX <= 0` and non-convergence fall through to the transition tests rather than returning; the
returned `AMPL2` is the iterated value and may exceed Ncrit; with no transition, `XT = X2` and the
`XT_*` sensitivities are left stale (`src/bl/system.rs:136-142, 294`, `src/bl/mrchue.rs:120`).
