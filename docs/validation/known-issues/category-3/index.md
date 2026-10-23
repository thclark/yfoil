---
icon: lucide/wrench
---

# 3. Drela's fixes

| | Where | Comment | yFoil |
|---|---|---|---|
| 3.1 | `xbl.f:934` (MRCHDU) | `fixed BUG   MD 7 June 99` — whether `CTAU` holds an amplification (laminar) or a shear coefficient (turbulent) at stations upstream of the previous ITRAN; `CTI <= 0 → 0.03` | Replicated verbatim, `src/bl/mrchdu.rs:103` |
| 3.2 | `xbl.f:742`, `:1066` | `added Ue clamp   MD  3 Apr 03` — Ue enters DMAX in the under-relaxation of both marches | Replicated, `src/bl/mrchue.rs:271`, `src/bl/mrchdu.rs:281` |
| 3.3 | `xbl.f:83` | `initialize BL by marching with Ue (fudge at separation)` — MRCHUE's direct/inverse switch at Hk limits | Replicated, `src/solver/setbl.rs:117` |
| 3.4 | `xpanel.f:1509` (XICALC, TE-gap cubic) | `ccc DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)  !!! BUG 2/2/95` — replaced by the `CROSP` form | Replicated (the live, fixed form), `src/solver/pointers.rs:143-147` |
| 3.5 | `xblsys.f:2458-2465` (HST) | `fudge HS slightly to make sure HS -> 2 as HK -> 1 (unnecessary with new correlation)` — left commented out | Replicated (omitted, as live code), `src/bl/closure.rs` |
| 3.6 | `xblsys.f:944` | `CCC CALL DIT(...)` — the dissipation closure `DIT` is dead | Not translated; `[[dead]]` in `coverage.toml` |
| 3.7 | `xfoil.f:1938` (PANGEN) | `fudge equations adjacent to TE to get TE panel length ratio RTF` | Replicated, `src/geometry/panel.rs:307` |
| 3.8 | `xfoil.f:1837` (PANGEN) | `temporarily used for more reliable convergence` — IPFAC = 5 oversampling that was never removed | Replicated, `repanel_by_curvature` |
| 3.9 | `xblsys.f:1990` (DAMPL) | `NEW VERSION. March 1991 (latest bug fix July 93)`; other closures dated 1991–94 | Shipped versions ported |
| 3.10 | `sort.f:244-247` | `Modified 4/24/01 HHY ... cures a bug for sharp LE foils where there were 3 LE points` | Out of scope (design path) |
| 3.11 | `xtcam.f:1289` | `ccc TOL = 1.0E-3*(S(N)-S(1))  ! Bad bug -- was losing x=1.0 point` | Out of scope (design path) |
