# 1.3 DUMP writes `H = H* = 1` where THET = 0 — Replicated

`xoper.f:1972-1978` (`IF(TH.EQ.0.0) THEN H = 1.0; HS = 1.0`). `src/output/foil.rs` does the same for
`stored.hs_dump` and for an unsolved station's live values.
