# 2.1 HVRAT is never assigned on the analysis path — Replicated

`HVRAT` (Sutherland-constant ratio in the `REYBL` viscosity law) is set to 0.35 only by the plotting
routines (`blplot.f:72`, `dplot.f:158, 413`). In a non-plotting run it keeps the static zero of
uninitialised COMMON and the viscosity law reduces to `HERAT**1.5`. `src/bl/system.rs:847-851` uses
`hvrat = 0.0`; using 0.35 gave `REYBL = 1e6 − 1 ULP` on the reference case where XFOIL gives 1e6. Two
unit tests whose expected values assumed 0.35 are permanently `#[ignore]`d (`tests/IGNORED.txt`).
Note: a run that has *plotted* a boundary layer before analysing continues with HVRAT = 0.35, so
XFOIL's own results depend on whether VPLO was visited. The DP reference build never plots.
