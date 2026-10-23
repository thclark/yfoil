# 7.10 The fallback tests are written so that a NaN residual takes the fallback — Replicated

MRCHUE and MRCHDU decide between keeping a failed station's last iterate and the iterate
fallback with `IF(DMAX .LE. 0.1) GO TO 109` (`xbl.f:785, 1109`): the fallback is taken
whenever the residual is *not* at or below 0.1, which includes a non-finite one. That is how a
station whose closures went non-finite (§4, TRCHEK2 row) beyond the third is set from the station upstream with finite values
and the march carries on. The same holds for every comparison in the reference: NaN makes it
false, and the branch that follows is whichever arm the false case leads to. A reproduction must
therefore keep each comparison's own operator and operand order — `if !(dmax <= 0.1)`, never the
complement `if dmax > 0.1`, which agrees for finite values and disagrees for NaN. yFoil does
(`src/bl/mrchue.rs`, `src/bl/mrchdu.rs`; the negated-comparison lint is allowed for this reason,
`Cargo.toml`). **Test:** `tests/execution/events.rs::test_mrchue_iterate_fallback_events`
and `::test_mrchdu_iterate_fallback_events_call_1` (the fallback stations of every non-finite
probe, from XFOIL's own `events.dat`, and the arrays after them, NaN where XFOIL has NaN).
