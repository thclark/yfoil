# 2.6 VSREZ(4) is stale at MRCHUE's inner-loop entry — Replicated (instrumentation note)

The fourth residual is assigned only in the direct/inverse branch, so a trace written before it holds
the previous iteration's value; the test excludes it (`tests/subroutine/mrchue.rs`).
