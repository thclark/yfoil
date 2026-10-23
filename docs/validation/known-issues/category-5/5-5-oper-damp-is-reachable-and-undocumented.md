# 5.5 OPER `DAMP` is reachable and undocumented — Replicated

`IDAMP` toggles the modified envelope method `DAMPL2` (`+0.1·exp(−20·HMI)` term). It is absent from the
OPER menu text but reachable; the first coverage measurement found it untranslated and it is now
ported and tested (`FlowConditions.amplification_model`, case `naca0012_n60_a2_re1e6_damp`).
