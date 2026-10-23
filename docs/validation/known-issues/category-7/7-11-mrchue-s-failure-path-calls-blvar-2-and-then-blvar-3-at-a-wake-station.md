# 7.11 MRCHUE's failure path calls BLVAR(2) and then BLVAR(3) at a wake station — Replicated

After a failed station Newton, MRCHUE recomputes the closures at label 109 with `BLVAR(1)` or
`BLVAR(2)` *and then* `BLVAR(3)` when the station is in the wake (`xbl.f:820-823`); MRCHDU does
the same (`xbl.f:1142-1144`). Each call clamps `HK2` in COMMON (§2.4), so a failed wake station
leaves with the turbulent floor of 1.05 applied before the wake closures are formed, and the
next station's "1" quantities carry that clamp. The effect is visible whenever
1.00005 < Hk < 1.05 at a failed wake station — the first march of the NACA 4412 at M 0.6, α 16°
fails at wake station 40 with Hk = 1.011. yFoil keeps the call sequence in both marches.
**Test:** `tests/execution/events.rs::test_mrchue_iterate_fallback_events` (that
case's first march, station 41 side 2 is where a single call would differ) and
`tests/subroutine/mrchdu.rs` (MRCHDU's failure path on the reference case).
