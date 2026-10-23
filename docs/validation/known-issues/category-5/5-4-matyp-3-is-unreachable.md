# 5.4 MATYP = 3 is unreachable — Replicated

The `TYPE` command maps TYPE 3 to `MATYP = 1, RETYP = 3` (`xoper.f:357-366`); nothing ever sets
`MATYP = 3`, so MRCL's third Mach branch (`xfoil.f:812, 817`) and SPECAL's `MINF_CLM = 0` branch
(`xoper.f:2784`) are dead. yFoil's `FlowConditions { mach_cl_dependence: Fixed, re_cl_dependence: InverseCl }` reproduces TYPE 3
(`tests/execution/coverage.rs`).
