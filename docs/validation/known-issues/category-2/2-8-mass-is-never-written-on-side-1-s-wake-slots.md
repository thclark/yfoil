# 2.8 MASS is never written on side 1's wake slots — Replicated (harness note)

UPDATE's final loop equates the "upper wake arrays" `(IBLTE(1)+K, 1)` to side 2's wake for CTAU,
THET, DSTR, UEDG, TAU, DIS, CTQ, DELT and TSTR (`xbl.f:1546-1557`) but not MASS, and nothing else
writes `MASS(IBL, 1)` beyond `NBL(1)` (MRCHUE/MRCHDU, UPDATE, DSSET and the SETBL/QVFUE sums all
run to `NBL(IS)`, and `NBL(1) = IBLTE(1)`). Those slots hold whatever an earlier stagnation-point
position left there and are read by nothing. yFoil's `mass_defect[1][..]` beyond `n_stations[1]` is
the same dead storage with a different history, so the replay harness (`tests/common/utilities/replay.rs`)
does not compare it; it surfaced in the polar break-point cases, where IST moves between iterations
and the instrumented `update_output_<k>.dat` dumps every row up to `NBL(1) + NW`.
