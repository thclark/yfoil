# 7.9 MRCHDU's station Newton at a transition station can amplify its seed a thousandfold — Replicated

The station Newton of MRCHDU (`xbl.f`, `DO 100 ITBL=1, 25`, `DEPS = 5.0E-6`) does not always contract: at a
transition station its residual DMAX can grow for many iterates before converging (on the NACA 0012 at Re 1e6,
7° of the 0 → 30° sweep, iteration 3, station 8 of the upper side: 18 iterates, DMAX 0.015 at the sixth rising
to 0.185 at the fifteenth). Through such a phase two runs that enter the station 1e-12 apart leave it 1e-9
apart, by a constant factor per iterate, on the same trajectory (the same count, DMAX to seven digits). From
*identical* inputs (yFoil's march seeded with the reference's dumped entering state) the same 18 iterates agree
to 8e-15 relative — the solver is the same, and the amplification acts on whatever seed it is given, rounding
included. The reference's own 1-ULP twins are amplified by it too, by factors that depend on the perturbation:
each independent random-sign ±1-ULP pattern of the panels leaves that station 1e-9 to 7e-9 apart, RLX at the
next iteration 1e-9 to 1e-8 from the reference, while an all-+1-ULP twin — a translation of the aerofoil —
leaves it only 1e-10 apart; yFoil enters it 2e-12 apart (its accumulated last-digit rounding, the size of the
twins' spread everywhere) and leaves it 1e-9 apart — at the bottom of the reference's own range. This is the
measurement that replaced the pipeline's single +1-ULP twin with five seeded −1/0/+1-ULP twins (CLAUDE.md
Rule 1). TRCHEK2 (`DAEPS = 5.0E-5`, the transition location's own Newton) is bit-exact from the
same inputs. It is what UPDATE's RLX
then reports whenever that station is the one that limits the step (RLX = DLO / DN there, one station,
undiluted). yFoil follows the same trajectory. **Test:** `tests/subroutine/mrchdu.rs::
test_mrchdu_station_newton_replays_from_xfoil_state_at_the_7deg_point` (every iterate of the station Newton
from XFOIL's exact state at `TOL_PURE`), `tests/subroutine/transition.rs` (TRCHEK2's every iterate from the
reference's dumped inputs, `trchek2_<k>.dat`), and `tests/execution/steps.rs` (every iteration of the 7° call
replayed from XFOIL's state, each within four times the step's own measured sensitivity). Each code's own
march through the call — the amplification acting on each code's own last-bit seed — is a divergent
comparison in the sense of `docs/development/conventions/terminology.md` and belongs to the series-cases study. The
instrumented reference dumps `trchek2_<k>.dat` and `mrchdu_trace_<k>.dat` for every `dump_calls` SETBL call.
