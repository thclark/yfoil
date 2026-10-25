---
icon: lucide/volume-x
---

# 4. Silent recovery

XFOIL prints and carries on in every case below; none aborts, none is flagged in the result. yFoil
replicates each (with the message as a code comment) so that its results and branch traces match.

| | Where | Behaviour | yFoil |
|---|---|---|---|
| 4.1 | `xfoil.f:796, 801` (MRCL) | `Illegal Re(CL)/Mach(CL) dependence trigger. Setting fixed` — RETYP/MATYP outside 1..3 silently reset to 1 | `src/solver/setbl.rs:39, 43`; unreachable through TYPE (§5.4) |
| 4.2 | `xfoil.f:845, 856` (MRCL) | `CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99`; Re limited to 100× REINF1 | `setbl.rs:77, 86`; the Mach limit and CL floor are checked against XFOIL's events by `tests/known_issues/kt_pole_7_7.rs` on `naca64a010_n60_inviscid_type2_m03_am23_a0_a19`; reach notes in `coverage.toml` |
| 4.3 | `xfoil.f` (CPCALC, CLCALC) | the Kármán–Tsien denominator `DEN = β + BFAC·Cp_inc` goes non-positive once any panel's speed passes `q/Q∞ = √(1 + 2β(1+β)/M∞²)`, so Cp does not merely lose accuracy — it passes through a **pole** and changes sign. CPCALC prints `Local speed too large. Compressibility corrections invalid.` once per call and returns the inverted Cp anyway; **CLCALC forms the same denominator inline and does not warn at all**, so CL/CM/CDP are silently wrong, frequently in sign. At M = 0.3 the threshold is q/Q∞ = 6.513 (a NACA 64A010 at N = 240 crosses it between 19° and 20°); at the M = 0.99 that MRCL's clamp produces it is 1.153, which that section exceeds by α = ±1°. See §7.7 for how TYPE 2 walks into this on its own | Replicated (`clcalc.rs:26` warns as CPCALC does, `:57-69` is the unguarded CLCALC form) |
| 4.4 | `xutils.f:41-50` (SETEXP) | 100-iteration Newton on the spacing ratio to `|dRatio| < 1e-5`; `Convergence failed. Continuing anyway ...` | `src/solver/xywake.rs:33-49`; the 1e-5 tolerance is load-bearing (§7.1) |
| 4.5 | `xoper.f` (VISCAL, SPECAL, SPECCL) | `Convergence failed` after ITMAX / 20 / 12 iterations; state kept | `viscal.rs:187`, `specal.rs:110, 172` |
| 4.6 | `xbl.f:1474` (UPDATE) | under-relaxation whenever a step would change a BL variable by more than 50 % (DLO) | `src/solver/update.rs` |
| 4.7 | `xbl.f:1515-1517` (UPDATE) | `eliminate absurd transients`: CTAU capped at 0.25 for turbulent stations only | `update.rs:296` |
| 4.8 | `xbl.f:1537` (UPDATE) | `make sure there are no "islands" of negative Ue` | `update.rs:311`; reached only in non-finite reference runs so far (`naca4412_n60_a16_re1e6_m06_iter60`) |
| 4.9 | `xbl.f:785-827`, `:1111-1149` | MRCHUE/MRCHDU iterate fallback: when the residual is not ≤ 0.1 (NaN included) the iterate is rejected (XFOIL's comment calls it `garbage`) and, beyond station 3, the station is set from the one upstream — θ, δ* scaled by √(ξ/ξ_prev) on the aerofoil, the TE values at the first wake station, δ* relaxed towards θ over 10 δ* further down the wake, Cτ 0.05 at transition or copied, Ue averaged from the neighbours (MRCHUE) or copied (MRCHDU); at stations 1–3 the rejected iterate is kept. Events `MRCHUE_ITERATE_FALLBACK`, `MRCHDU_ITERATE_FALLBACK` | Translated. MRCHUE's fallback is exercised at airfoil, wake-start and wake stations by `naca4412_n60_a18_re3e6_iter40` (compared: match); MRCHDU's only ever fired in runs the reference had already driven non-finite (next row), reach notes in `coverage.toml` |
| 4.10 | `xblsys.f` (BLVAR) | laminar Cf used in the turbulent branch at unreasonably small Rθ | `src/bl/system.rs:1287, 1408` |
| 4.11 | `xblsys.f` (BLVAR) | `DE` capped at `12·θ` with all four sensitivities zeroed — a Jacobian discontinuity | `system.rs:1458-1465` |
| 4.12 | `xblsys.f` (BLDIF) | similarity-station amplification row is a dummy `VS2(1,1) = 1` pivot | `system.rs:1744-1746` |
| 4.13 | `xblsys.f` (BLDIF) | forms `UQ_RTA`, `UQ_T1..UQ_RE` and never uses them | omitted where they feed nothing, `system.rs:1830` |
| 4.14 | `xpanel.f:1737` (UESET) | `tweak Ue so it's not zero, in case stag. point is right on node` (UEPS = 1e-7) | `pointers.rs:325` |
| 4.15 | `xpanel.f:1385` (STFIND) | `tweak stagnation point if it falls right on a node (very unlikely)` — needs `GAM(I) == 0` bitwise | `pointers.rs:47`; a divergent comparison wherever it is reached (`docs/development/conventions/terminology.md`) |
| 4.16 | `xbl.f` (SETBL) | `SETBL: Xtr???  n1 n2:` diagnostic when the transition interval and ITRAN disagree | `setbl.rs:272` |
| 4.17 | `plotlib/plt_font.f:249` (PLNUMBABS, via ASEQ → SEQPLT at `xoper.f:719` and via ALFA → CPX → COEFPL) | plot labels are formatted even with graphics off (`PLOP / G F`), and PLNUMBABS's digit-extraction loop never terminates on a non-finite value. Reached by both the sequence-plot label and the Cp-plot coefficient label, so every OPER point whose CL/CD has become Infinity hangs XFOIL at 100 % CPU (NACA 0012, `ITER 100`, past ~21° on either leg; stacks sampled 2026-09-10). The `input-sensitivity` driver and the fixture runner detect the stall (stdout stops growing) and kills the run | Out of scope (plot library); the sequence data up to the hang is intact. yFoil has no plot label and halts the sequence normally (`compute_polar`) |
| 4.18 | `xblsys.f:396-425` (TRCHEK2) | the transition-point Newton can drive `XT` onto `X2` (the interval end); the station's Newton then fails with `Res = NaN`, the iterate fallback above sets the station from the one upstream with finite values (beyond station 3) and the run carries on (`MRCHUE: Convergence failed at 37 side 1 Res = NaN`, `x: 0.88282 0.88572 0.88572 N: 0.683 9.000 NaN`). Seen on the NACA 63-415 (closed TE) from α = 10° and on stalled compressible points (`docs/validation/branch-coverage/`, the `non-finite` cases of `cases.toml`) | Replicated at the event level: the fallback and every other branch these runs reach are tested one call at a time by `tests/execution/events.rs` (NaN where XFOIL has NaN); the whole run is not, because the reference's own 1-ULP twins wander by O(10²) on it. Such runs are excluded from the validation set and the branches only they reach are reported as *non-finite only* (§8); the fallback's NaN behaviour is §7.10 |

!!! note "Which reference reaches the plot-label hang (4.17)"

    Whether a sweep reaches a non-finite CL or CD, and so the hang, depends on the maths library,
    because the runs that reach it are ill-conditioned. Four known-issues sweeps hung on the
    host-libm reference used before 2026-10-15 and do not hang on the CORE-MATH reference:
    `naca4412_n240_polar30_re1e6_iter200` (30°), `naca4412_n240_polar30_re3e6_iter200` (−24°),
    `naca16-212_n160_polar30_re1e6_m07_iter100` (−18°) and, with the panels yFoil generated then,
    `naca64a010_n160_polar30_re1e6_type2_iter100` (12°). Their `pathological_…_host_libm` copies in
    `cases.toml` recreate those hangs on `scripts/xfoil-build.sh --host-libm`. Because that build
    takes its maths from the host, it is not a reference, and what it does depends on the host
    (`tests/fixtures/README.md`). On CORE-MATH, `naca0012_n240_polar30_re1e6_iter200`,
    `naca4412_n240_polar30_re1e6_m06_iter200` and `naca64a010_n160_polar30_re1e6_m03_type2_iter100`
    still hang (checked 2026-10-25).

## Runs that show it

- [NACA 16-212 at M 0.7](../runs/naca-16-212-at-m-0-7.md)
- [NACA 64A010, TYPE 2](../runs/naca-64a010-type-2.md)
- [NACA 4412 inviscid at M 0.7](../runs/naca-4412-inviscid-at-m-0-7.md)
- [NACA 0012, Re 1e6, N = 240](../runs/naca-0012-re-1e6-n-240.md)
- [NACA 4412, Re 1e6, N = 240](../runs/naca-4412-re-1e6-n-240.md)
- [NACA 4412, Re 3e6, N = 240](../runs/naca-4412-re-3e6-n-240.md)
- [NACA 4412 at M 0.6, N = 240](../runs/naca-4412-at-m-0-6-n-240.md)
- [NACA 63-415 at M 0.5](../runs/naca-63-415-at-m-0-5.md)
- [NACA 64A010, TYPE 2 at M 0.3](../runs/naca-64a010-type-2-at-m-0-3.md)
- [NACA 23012 with a trip upstream of stagnation](../runs/naca-23012-with-a-trip-upstream-of-stagnation.md)
