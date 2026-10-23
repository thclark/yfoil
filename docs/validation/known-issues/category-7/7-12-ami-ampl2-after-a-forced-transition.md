# 7.12 AMI = AMPL2 after a forced transition — Overcome (yFoil bug, fixed)

MRCHUE, MRCHDU and SETBL set `AMI = AMPL2` after every `CALL TRCHEK` (`xbl.f:225`, `:624`, `:815`,
`:968`, `:1137`), whatever TRCHEK2 decided: after a forced transition (`XIFORC` in the interval — with
the default `XSTRIP = 1` that is the trailing-edge interval of a side still laminar there) AMPL2 is the
N2 Newton's converged value. yFoil's `TransitionCheck::Forced` carried no `ampl2`, so its callers kept
the previous AMI. The converged state is unaffected — AMI at a transition station feeds only the
starting value of the next TRCHEK2 call at that station — so no value test saw it; it changed the
number of TRCHEK2 iterates, which `cargo xtask route` found on 2026-10-06 (NACA 0012 at 7° of the
0 → 30° sweep, lower side, station 29: AMPL2 8.645 against 8.830). Fixed by carrying `ampl2` in
`Forced`. **Test:** `tests/subroutine/mrchdu.rs::
test_mrchdu_station_newton_replays_from_xfoil_state_at_the_7deg_point` (AMPL2 at every station
iterate, both sides).
