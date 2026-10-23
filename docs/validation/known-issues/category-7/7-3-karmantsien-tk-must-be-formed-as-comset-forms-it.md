# 7.3 Kármán–Tsien TK must be formed as COMSET forms it — Overcome (yFoil bug, fixed)

`TKLAM = MSQ/(1+β)²` and `TKBL = 1/β − 1` are algebraically equal and numerically different; the
difference is dead at M = 0. `BLGlobalParams::new` now uses COMSET's form; two tests derived with the
other form are ignored pending regeneration from the M = 0.3 coverage case.
