# Polars and reinitialisation

A polar is a state machine, not N independent solves: the boundary layer converged at one angle is
the next one's initial condition ([sweeping a polar](../../using-yfoil/commands/sweeping-a-polar.md#how-the-sweep-runs)).
So yFoil models XFOIL's state flags explicitly — `LWAKE`, `LIPAN`, `LBLINI`, `LADIJ`, `LWDIJ`,
`LVCONV`, `LGAMU`, `LQAIJ` and `AWAKE`/`AVISC`/`MVISC` — and a sweep is one persistent
`Session` driven through the same commands XFOIL's OPER receives.

## The polar run test

The tracked polar case `naca0012_n60_polar_re1e6` runs XFOIL's OPER script
`ALFA 0 / ASEQ 1 5 1 / INIT / ALFA 0 / ASEQ -1 -5 -1` on yFoil's panels; `viscal_points.dat`
records every VISCAL call (ITMAX passed, `LBLINI` on entry, iterations, `LVCONV`, CL, CM, CD, CDF,
CDP, XOCTR, IST, ITRAN). `tests/execution/polar.rs` replays the same sequence through one
`Session`, and every record must match: identical iteration counts, convergence, IST and ITRAN;
converged forces within `TOL_SOLVER`; per-iteration transients within `TOL_TRANSIENT`.

## The re-solved 0°

OPER cannot save a boundary-layer state, so the downward leg starts by reinitialising (`INIT`) and
re-solving 0° with a fresh MRCHUE march. That point seeds the leg and is not a polar point. Two
facts make that sound, and each is tested:

- **in XFOIL**, the re-solved 0° (VISCAL call 7 of the polar case) is the first 0° again,
  textually identical `ES24.16` records — `tests/apparatus/polar_reseed.rs`;
- **in yFoil**, the polar driver's re-solved 0° is its first 0° again —
  `tests/application/polar_driver.rs`.

That identity is also what lets yFoil's content-addressed ids cite the first solve for the
seeding one ([output data](../../using-yfoil/output-data.md#point-identifiers-and-the-initialisation-chain)).
