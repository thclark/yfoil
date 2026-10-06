---
icon: lucide/triangle-alert
---

# Known XFOIL issues

XFOIL is the reference (CLAUDE.md Rule 1) and its algorithm is reproduced including its own bugs and
quirks (Rule 2). This page is the register of everything found while translating it that a reader
comparing the two codes, or reading XFOIL's output, needs to know about. Each entry says what XFOIL
does, cites both sources, and states yFoil's handling:

- **Replicated** — yFoil does the same thing, and the behaviour is fixture-gated where a gate exists.
- **Overcome** — yFoil deliberately differs. Only one such divergence affects results (the NACA
  generator); the others are output-side or instrumentation-side and leave the solver untouched.
- **Out of scope** — code XFOIL itself never executes on the analysis path, or design/plotting
  features yFoil does not translate. Listed so nobody looks for them.

Line numbers refer to `xfoil/third-party/xfoil-6.99/src/` and to yFoil at the commit that last
touched this page. Every "unreachable" claim below is also carried, with its reason, in
`xtask/fixtures-config/coverage.toml` and measured by `cargo xtask coverage`
(`docs/validation/coverage.md`).

---

## 1. Output quantities

### 1.1 DUMP prints one-step-stale closure quantities (H*, Cf, K, tau, Di) — Overcome (output only)

Verified in the source, not inferred:

- SETBL calls MRCHDU at the top of **every** Newton iteration (`xbl.f:93`, matched at
  `src/solver/setbl.rs:125`). MRCHDU re-marches and stores `THET, DSTR, UEDG, CTAU, MASS, TAU, DIS,
  CTQ, DELT, TSTR` (`xbl.f:1157-1167`). SETBL then evaluates the closures on that state
  (BLPRV/BLKIN/BLVAR, `xbl.f:219-220, 470`) and stores `TAU, DIS, CTQ, DELT, USLP` (`xbl.f:277-282`,
  matched at `setbl.rs:316-320`). Nothing writes TSTR in SETBL. Only after that do BLSOLV and UPDATE
  apply the Newton correction to `THET, DSTR, UEDG, CTAU, MASS`.
- So at the end of VISCAL the primaries are post-final-correction, while `TAU, DIS, CTQ, DELT, USLP,
  TSTR` describe the state *entering* the final iteration: they lag by exactly the last Newton
  correction.
- **Consequence for XFOIL's own output.** DUMP (`xoper.f:1955-1995`) prints `Ue, Dstar, Theta, H, HK,
  P, m` from current values but `Cf = TAU/(½Q²)`, `H* = TSTR/THET`, `K = TSTR·Ue³`, `tau`, `Di` from
  the lagged arrays. `H*` and `K` are the worst case: a lagged numerator over a current denominator.
  VPLO's `CF`, `CD` and `DELT` plots use `TAU, DIS, DELT` (`blplot.f:399, 499, 573`) and are lagged
  too; `H` is recomputed live with HKIN (`blplot.f:1385`), and the `TSTR/THET` version of that plot is
  commented out (`blplot.f:1387`), consistent with Drela knowing. `DT, DB, UE, N, CT` are live.
- Magnitude: the final correction, bounded by the convergence test (RMSBL < 1e-4 on the normalised
  corrections), so ~1e-4 relative or better on a converged point; unbounded on an unconverged one.
  On the reference case (NACA 0012, N=60, α=2°, Re=1e6) live and DUMP `H*` differ by ~1e-8.

**yFoil.** The stored arrays are reproduced exactly (`tests/subroutine/mrchdu.rs`,
`tests/execution/update.rs`, the SETBL gates); the solver is untouched. In the analysis JSON
(`BlSideOutput`, `src/output/foil.rs`) the canonical columns `hs, cf, cdis, delta, ctq, uslp` are
computed *live* by running XFOIL's own BLPRV → BLKIN → BLVAR on the converged primaries, and the lagged
arrays are emitted verbatim under `stored` (with `hs_dump`, `cf_dump` as DUMP prints them). A
side-by-side of live `hs` against DUMP's `H*` will differ at the ~RMSBL level by construction;
`tests/application/output.rs` bounds it. Column table: `docs/xfoil-reference.md`,
"Output" section.

### 1.2 DUMP's `Ue/Vinf` is signed by GAM — Overcome (output only)

DUMP forms `UE = (GAM(I)/QINF)(1-TKLAM)/(1-TKLAM(GAM(I)/QINF)²)` (`xoper.f:1987`), so the lower surface
prints negative Ue. yFoil's `ue` column applies the same Kármán–Tsien transformation to UEDG (BLPRV's
`U2/QINF`), which is unsigned on both sides. `|DUMP Ue| == ue` on the surface.

### 1.3 DUMP writes `H = H* = 1` where THET = 0 — Replicated

`xoper.f:1972-1978` (`IF(TH.EQ.0.0) THEN H = 1.0; HS = 1.0`). `src/output/foil.rs` does the same for
`stored.hs_dump` and for an unsolved station's live values.

### 1.4 CPDISP's displacement surface has no scale factor and splits the wake δ* by TE values — Replicated, extended

`xplots.f:699-754`: the displacement surface is `X + NX·DSTR` at true geometric scale; the wake δ* is
split into upper/lower fractions `DSF1 = (DSTR(IBLTE(1),1) + ½ANTE)/DSTR(IBLTE(2)+1,2)` (and DSF2
likewise) with a 0.5/0.5 fallback when the first wake δ* is exactly zero, and the wake normal points
to the *lower* side so the upper edge is `X − N·DSTR·DSF1`. `src/output/foil_plot.rs` follows this
exactly for δ*, and additionally (yFoil's own): a per-quantity scale factor, the same construction for
other quantities with a 0.5/0.5 split, and markers for stagnation, transition and separation, which
XFOIL never draws on the airfoil plot. The separation/reattachment markers are derived from the sign of
the live Cf; XFOIL reports no separation location anywhere.

### 1.5 Formatted output precision — Overcome (instrumentation)

`PSAVE` is `G15.7`, `.pol` is `F9.4/F10.5`, `DUMP` is `F9.5/F10.6`, `CPWR` is `F11.5`: none can
support a comparison below ~1e-5, and XFOIL's stock debug formats `E24.16/E18.10/E12.4` carry 16, 10
and 4 significant figures. The instrumented reference build writes `ES24.16` (17 figures, the minimum
that round-trips an f64); the instrumentation is proven inert (byte-identical pristine vs instrumented
outputs, `scripts/xfoil-build.sh --verify`). Nothing gates on XFOIL's formatted files (CLAUDE.md
Rule 4). An old validation table that showed 4.9e-4 to 1.6e-3 RMS "paneling error" was measuring
7-digit PSAVE output plus the two NACA definitions of §6.1 (`docs/validation/geometry/README.md`).

---

## 2. Uninitialised, stale and COMMON-persistent state

### 2.1 HVRAT is never assigned on the analysis path — Replicated

`HVRAT` (Sutherland-constant ratio in the `REYBL` viscosity law) is set to 0.35 only by the plotting
routines (`blplot.f:72`, `dplot.f:158, 413`). In a non-plotting run it keeps the static zero of
uninitialised COMMON and the viscosity law reduces to `HERAT**1.5`. `src/bl/system.rs:847-851` uses
`hvrat = 0.0`; using 0.35 gave `REYBL = 1e6 − 1 ULP` on the reference case where XFOIL gives 1e6. Two
unit tests whose expected values assumed 0.35 are permanently `#[ignore]`d (`tests/IGNORED.txt`).
Note: a run that has *plotted* a boundary layer before analysing continues with HVRAT = 0.35, so
XFOIL's own results depend on whether VPLO was visited. The DP reference build never plots.

### 2.2 MRCHDU locals persist across sides — Replicated

`AMI, SENS/SENNEW, UEREF/HKREF, CTE/TTE/DTE` are Fortran locals never re-initialised between side 1
and side 2, and `COM1/COM2/XT` are COMMON (`src/bl/mrchdu.rs:64-65`). A per-side reset diverges from
XFOIL. The same persistence in MRCHUE is trace-faithful only: side 2's similarity station overwrites
COM1 before anything reads it (`src/bl/mrchue.rs:54-56`).

### 2.3 TRDIF's transition-point station reuses station-2 COMMON — Replicated

`xblsys.f:355` "temporarily set '2' variables from 'T' for BLKIN": XFOIL overwrites
`X2/T2/D2/U2/AMPL2/S2` in place without calling BLPRV, so `U2_UEI`, `U2_MS` and `DW2` at the transition
point are still station 2's (`src/bl/system.rs:2191-2193`).

### 2.4 BLVAR clamps HK2 in COMMON without touching its derivatives — Replicated

`xblsys.f:798-799`: `HK2 = MAX(HK2, 1.00005)` in the wake, `MAX(HK2, 1.05)` otherwise, while
`HK2_U/T/D/MS` stay as BLKIN set them. The clamped value then feeds BLMID, the next station's HK1
(after `COM1 = COM2`) and MRCHUE's HTARG (`src/bl/system.rs:1171-1176`). The wake `US2` is likewise
clamped (`xblsys.f:833, 843`). The order of the repeated BLVAR calls in MRCHDU matters only because of
this in-place clamp (`src/bl/mrchdu.rs:385-387`).

### 2.5 TRCHEK2 leaves AMPL2, XT and XT_* as the loop left them — Replicated

After 30 Newton iterations TRCHEK2 prints `N2 convergence failed.` and continues (`xblsys.f:277`);
`AX <= 0` and non-convergence fall through to the transition tests rather than returning; the
returned `AMPL2` is the iterated value and may exceed Ncrit; with no transition, `XT = X2` and the
`XT_*` sensitivities are left stale (`src/bl/system.rs:136-142, 294`, `src/bl/mrchue.rs:120`).

### 2.6 VSREZ(4) is stale at MRCHUE's inner-loop entry — Replicated (instrumentation note)

The fourth residual is assigned only in the direct/inverse branch, so a trace written before it holds
the previous iteration's value; the gate excludes it (`tests/subroutine/mrchue.rs`).

### 2.7 RADLE stays zero for a flat LE — Replicated

`xgeom.f:350`: LE curvature below `0.001/(S(N)-S(1))` leaves `RADLE = 0` (flat plate, wedge). Open
branch with a reach note in `coverage.toml`.

### 2.8 MASS is never written on side 1's wake slots — Replicated (harness note)

UPDATE's final loop equates the "upper wake arrays" `(IBLTE(1)+K, 1)` to side 2's wake for CTAU,
THET, DSTR, UEDG, TAU, DIS, CTQ, DELT and TSTR (`xbl.f:1546-1557`) but not MASS, and nothing else
writes `MASS(IBL, 1)` beyond `NBL(1)` (MRCHUE/MRCHDU, UPDATE, DSSET and the SETBL/QVFUE sums all
run to `NBL(IS)`, and `NBL(1) = IBLTE(1)`). Those slots hold whatever an earlier stagnation-point
position left there and are read by nothing. yFoil's `mass_defect[1][..]` beyond `n_stations[1]` is
the same dead storage with a different history, so the replay harness (`tests/common/utilities/replay.rs`)
does not compare it; it surfaced in the polar break-point cases, where IST moves between iterations
and the instrumented `update_output_<k>.dat` dumps every row up to `NBL(1) + NW`.

---

### 2.9 The closure log's `s2` at laminar BLVAR calls is not a value BLVAR uses — Out of scope

The reference's closure log (`xfoil_subroutine_log.dat`, cut into `tests/fixtures/subroutines/` by
`closures = true`) records BLVAR's inputs as they stand in COM2. At some laminar calls (`ityp = 1`,
the stations near the leading edge) `S2` holds a subnormal (about 3e-314) that changed, with every
output unchanged, when the instrumentation's dump list grew from 32 to 1000 entries — a value that
depends on memory layout, so not one XFOIL computed for that station. BLVAR(1) does not read `S2`
(the laminar closures take no shear stress), so the outputs and the gate are unaffected; the input
is regenerated with the fixtures and `--verify` holds it for a given build. Noted 2026-10-06.

## 3. Drela's own dated fixes and flagged defects in the shipped source

| Where | Comment | yFoil |
|---|---|---|
| `xbl.f:934` (MRCHDU) | `fixed BUG   MD 7 June 99` — whether `CTAU` holds an amplification (laminar) or a shear coefficient (turbulent) at stations upstream of the previous ITRAN; `CTI <= 0 → 0.03` | Replicated verbatim, `src/bl/mrchdu.rs:103` |
| `xbl.f:742`, `:1066` | `added Ue clamp   MD  3 Apr 03` — Ue enters DMAX in the under-relaxation of both marches | Replicated, `src/bl/mrchue.rs:271`, `src/bl/mrchdu.rs:281` |
| `xbl.f:83` | `initialize BL by marching with Ue (fudge at separation)` — MRCHUE's direct/inverse switch at Hk limits | Replicated, `src/solver/setbl.rs:117` |
| `xpanel.f:1509` (XICALC, TE-gap cubic) | `ccc DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)  !!! BUG 2/2/95` — replaced by the `CROSP` form | Replicated (the live, fixed form), `src/solver/pointers.rs:143-147` |
| `xblsys.f:2458-2465` (HST) | `fudge HS slightly to make sure HS -> 2 as HK -> 1 (unnecessary with new correlation)` — left commented out | Replicated (omitted, as live code), `src/bl/closure.rs` |
| `xblsys.f:944` | `CCC CALL DIT(...)` — the dissipation closure `DIT` is dead | Not translated; `[[dead]]` in `coverage.toml` |
| `xfoil.f:1938` (PANGEN) | `fudge equations adjacent to TE to get TE panel length ratio RTF` | Replicated, `src/geometry/panel.rs:307` |
| `xfoil.f:1837` (PANGEN) | `temporarily used for more reliable convergence` — IPFAC = 5 oversampling that was never removed | Replicated, `repanel_by_curvature` |
| `xblsys.f:1990` (DAMPL) | `NEW VERSION. March 1991 (latest bug fix July 93)`; other closures dated 1991–94 | Shipped versions ported |
| `sort.f:244-247` | `Modified 4/24/01 HHY ... cures a bug for sharp LE foils where there were 3 LE points` | Out of scope (design path) |
| `xtcam.f:1289` | `ccc TOL = 1.0E-3*(S(N)-S(1))  ! Bad bug -- was losing x=1.0 point` | Out of scope (design path) |

---

## 4. Silent self-corrections and continue-on-failure

XFOIL prints and carries on in every case below; none aborts, none is flagged in the result. yFoil
replicates each (with the message as a code comment) so that its results and branch traces match.

| Where | Behaviour | yFoil |
|---|---|---|
| `xfoil.f:796, 801` (MRCL) | `Illegal Re(CL)/Mach(CL) dependence trigger. Setting fixed` — RETYP/MATYP outside 1..3 silently reset to 1 | `src/solver/setbl.rs:39, 43`; unreachable through TYPE (§5.4) |
| `xfoil.f:845, 856` (MRCL) | `CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99`; Re limited to 100× REINF1 | `setbl.rs:77, 86`; the Mach limit and CL floor are checked against XFOIL's events by `tests/known_issues/kt_pole_7_7.rs` on `naca64a010_n60_inviscid_sweep30_type2_m03`; reach notes in `coverage.toml` |
| `xfoil.f` (CPCALC, CLCALC) | the Kármán–Tsien denominator `DEN = β + BFAC·Cp_inc` goes non-positive once any panel's speed passes `q/Q∞ = √(1 + 2β(1+β)/M∞²)`, so Cp does not merely lose accuracy — it passes through a **pole** and changes sign. CPCALC prints `Local speed too large. Compressibility corrections invalid.` once per call and returns the inverted Cp anyway; **CLCALC forms the same denominator inline and does not warn at all**, so CL/CM/CDP are silently wrong, frequently in sign. At M = 0.3 the threshold is q/Q∞ = 6.513 (a NACA 64A010 at N = 240 crosses it between 19° and 20°); at the M = 0.99 that MRCL's clamp produces it is 1.153, which that section exceeds by α = ±1°. See §7.7 for how TYPE 2 walks into this on its own | Replicated (`clcalc.rs:26` warns as CPCALC does, `:57-69` is the unguarded CLCALC form) |
| `xutils.f:41-50` (SETEXP) | 100-iteration Newton on the spacing ratio to `|dRatio| < 1e-5`; `Convergence failed. Continuing anyway ...` | `src/solver/xywake.rs:33-49`; the 1e-5 tolerance is load-bearing (§7.1) |
| `xoper.f` (VISCAL, SPECAL, SPECCL) | `Convergence failed` after ITMAX / 20 / 12 iterations; state kept | `viscal.rs:187`, `specal.rs:110, 172` |
| `xbl.f:1474` (UPDATE) | under-relaxation whenever a step would change a BL variable by more than 50 % (DLO) | `src/solver/update.rs` |
| `xbl.f:1515-1517` (UPDATE) | `eliminate absurd transients`: CTAU capped at 0.25 for turbulent stations only | `update.rs:296` |
| `xbl.f:1537` (UPDATE) | `make sure there are no "islands" of negative Ue` | `update.rs:311`; reached only in non-finite reference runs so far (`naca4412_n60_a16_re1e6_m06_iter60`) |
| `xbl.f:785-827`, `:1111-1149` | MRCHUE/MRCHDU non-convergence fallbacks: extrapolate from the previous station when the residual > 0.1 (`garbage solution`) | Translated. MRCHUE's fallback is exercised at airfoil, wake-start and wake stations by `naca4412_n60_a18_re3e6_iter40` (gated: match); MRCHDU's only ever fired in runs the reference had already driven non-finite (next row), reach notes in `coverage.toml` |
| `xblsys.f` (BLVAR) | laminar Cf used in the turbulent branch at unreasonably small Rθ | `src/bl/system.rs:1287, 1408` |
| `xblsys.f` (BLVAR) | `DE` capped at `12·θ` with all four sensitivities zeroed — a Jacobian discontinuity | `system.rs:1458-1465` |
| `xblsys.f` (BLDIF) | similarity-station amplification row is a dummy `VS2(1,1) = 1` pivot | `system.rs:1744-1746` |
| `xblsys.f` (BLDIF) | forms `UQ_RTA`, `UQ_T1..UQ_RE` and never uses them | omitted where they feed nothing, `system.rs:1830` |
| `xpanel.f:1737` (UESET) | `tweak Ue so it's not zero, in case stag. point is right on node` (UEPS = 1e-7) | `pointers.rs:325` |
| `xpanel.f:1385` (STFIND) | `tweak stagnation point if it falls right on a node (very unlikely)` — needs `GAM(I) == 0` bitwise | `pointers.rs:47`; a divergent comparison wherever it is reached (`docs/conventions/terminology.md`) |
| `xbl.f` (SETBL) | `SETBL: Xtr???  n1 n2:` diagnostic when the transition interval and ITRAN disagree | `setbl.rs:272` |
| `plotlib/plt_font.f:249` (PLNUMBABS, via ASEQ → SEQPLT at `xoper.f:719` and via ALFA → CPX → COEFPL) | plot labels are formatted even with graphics off (`PLOP / G F`), and PLNUMBABS's digit-extraction loop never terminates on a non-finite value. Reached by both the sequence-plot label and the Cp-plot coefficient label, so every OPER point whose CL/CD has become Infinity hangs XFOIL at 100 % CPU (NACA 0012, `ITER 100`, past ~21° on either leg; stacks sampled 2026-09-10). The `xfoil-sensitivity` driver detects the stall (stdout stops growing) and kills the run | Out of scope (plot library); the sequence data up to the hang is intact. yFoil has no plot label and halts the sequence normally (`compute_polar`) |

---
| `xblsys.f:396-425` (TRCHEK2) | the transition-point Newton can drive `XT` onto `X2` (the interval end); the station's Newton then fails with `Res = NaN`, the fallback above extrapolates finite values over it and the run carries on (`MRCHUE: Convergence failed at 37 side 1 Res = NaN`, `x: 0.88282 0.88572 0.88572 N: 0.683 9.000 NaN`). Seen on the NACA 63-415 (closed TE) from α = 10° and on stalled compressible points (`docs/validation/branch-coverage/`, the `non-finite` cases of `cases.toml`) | Replicated at the event level: the fallback and every other branch these runs reach are gated one call at a time by `tests/execution/events.rs` (NaN where XFOIL has NaN); the whole run is not, because the reference's own 1-ULP twins wander by O(10²) on it. Such runs are excluded from the validation set and the branches only they reach are reported as *non-finite only* (§8); the fallback's NaN behaviour is §7.10 |

## 5. Dead and unreachable code in 6.99

### 5.1 Subroutines never called — Not translated

`UECALC` and `DSSET` (`xpanel.f`) are never called anywhere in 6.99, and `DIT`'s only call site is
commented out (`xblsys.f:944`). None of the three is translated: they were carried for completeness
until 2026-09-13 and removed once the branch-coverage study confirmed no case can reach them, so the
translated set of `xtask/fixtures-config/coverage.toml` no longer names them.

### 5.2 LIMAGE is permanently false — Out of scope

The ground-image toggle is commented out (`xoper.f:405-410`), so PSILIN's image-airfoil block
(`xpanel.f:488` and 12 dependent lines) is dead. Not translated (`src/solver/psilin.rs:6-7`).
GEOLIN (geometric sensitivities) is inverse-design only and also not translated.

### 5.3 Stale duplicate source files — Out of scope

`dplot1.f` is an older snapshot of `dplot.f` (missing `FBLGET`); `xqdes1.f`, `modify1/2/3.f`,
`frplot0.f` and `xpanel.new` are the same pattern. None is in the build (`bin/Makefile` links
`dplot.o` only) and none is in `coverage.toml`'s file list.

### 5.4 MATYP = 3 is unreachable — Replicated

The `TYPE` command maps TYPE 3 to `MATYP = 1, RETYP = 3` (`xoper.f:357-366`); nothing ever sets
`MATYP = 3`, so MRCL's third Mach branch (`xfoil.f:812, 817`) and SPECAL's `MINF_CLM = 0` branch
(`xoper.f:2784`) are dead. yFoil's `FlowConditions { mach_cl_dependence: Fixed, re_cl_dependence: InverseCl }` reproduces TYPE 3
(`tests/execution/coverage.rs`).

### 5.5 OPER `DAMP` is reachable and undocumented — Replicated

`IDAMP` toggles the modified envelope method `DAMPL2` (`+0.1·exp(−20·HMI)` term). It is absent from the
OPER menu text but reachable; the first coverage measurement found it untranslated and it is now
ported and gated (`FlowConditions.amplification_model`, case `naca0012_n60_a2_re1e6_damp`).

### 5.6 Structurally dead branches — Out of scope

Coincident/duplicate consecutive nodes (ABCOPY deletes them on LOAD, PANGEN never produces them):
`xpanel.f:29, 66, 77, 81, 180, 861, 867, 875`, `spline.f:542-543`. SPLIND's `999`/specified-slope end
conditions (only SEGSPL's `-999/-999` is reached on the analysis path): `spline.f:96, 101, 113, 117`.
Array-bound `STOP`s that yFoil, having no fixed dimensions, cannot hit. Interactive prompts (`NITER = 0`,
`NACA` without a designation). Flap hinge moments (`LFLAP`, `LBFLAP`). `LIPAN` cleared with `LBLINI`
kept (`xqdes.f:497` only). The impossible flag combinations `LADIJ = F ∧ LWDIJ = T` and
`LQAIJ = F ∧ LGAMU = T`. Each is an `[[unreachable]]` entry with a class and reason in `coverage.toml`;
84 in total.

---

## 6. Geometry

### 6.1 NACA4 applies thickness vertically — Overcome (the one live divergence)

`naca.f:62`: `YB(IB) = YC(I) + YT(I)`, i.e. thickness added vertically rather than perpendicular to
the camber line, which is not the NACA definition (irrelevant for symmetric sections). yFoil's
generators lay the thickness perpendicular to the mean line by default, whatever the panelling
method (`yfoil geometry naca …`, `Section`). XFOIL's model is `--thickness vertical`
(`naca_4digit_vertical`, `naca_5digit_vertical`): XFOIL's own 245-point NACA4/NACA5 buffer, always
PANGEN-panelled, and it exists only to replicate the output of XFOIL's `NACA` command bitwise
(`tests/subroutine/pangen.rs`); the geometry file records `thickness_applied: "vertical"`. This
never enters a comparison because yFoil generates the panels and XFOIL consumes them (Rule 4).
Recorded in CLAUDE.md's divergence table, the only row.

### 6.2 LOAD reverses clockwise input — Constraint

`xfoil.f:1251-1260`: a negative signed area flips the node order. yFoil writes counter-clockwise
(TE → upper → LE → lower → TE) so the bitwise handoff holds; `LNORM` defaults to false
(`xfoil.f:495`) so LOAD does not rescale.

### 6.3 CHORD and APANEL are not what their names suggest — Replicated

`CHORD` (GEOPAR) is the distance from the LE point *on the spline* to the TE midpoint, 0.99995721 for
the 60-panel NACA 0012 rather than 1; `APANEL` (APCALC) is the panel-*normal* angle, 3π/2 from a
node-tangent angle. Both were yFoil misreadings found by the S3 gate; both are now XFOIL's, and
`CHORD` is what SETEXP scales the wake by.

### 6.4 NW = N/12 + 10·INT(WAKLEN) — Replicated

Integer truncation in both terms (`xpanel.f:1280`); `src/solver/analysis.rs:90`.

### 6.5 TGAP on a sharp trailing edge does not deliver the requested gap — Replicated

`xgdes.f:1238-1244` (TGAP): when the buffer airfoil's trailing-edge gap is zero, the direction
along which the two surfaces are moved apart is `(DXU, DYU) = (−½(YBP(NB) − YBP(1)), ½(XBP(NB) −
XBP(1)))`, the mean of the two end tangents of the spline. That vector is unit length only for a
cusp: for a trailing-edge angle τ it has length cos(τ/2) (times the small stretch of the
spline-parameter tangents), and the surfaces end up `cos(τ/2)` × the requested gap apart. With a
gap already open the direction is the existing gap's unit vector and the requested gap is
delivered exactly. Measured on the tracked `naca63-415_n160_tgap` case (τ ≈ 7°): `TGAP 0.002`
gives 2.0012e-3. yFoil's `set_te_gap` reproduces the moved nodes bitwise
(`tests/subroutine/tgap.rs`); the resulting gap is asserted against the requested one only when a
gap existed. The gap TGAP produces is what `TECALC` then measures, so the `SHARP` decision and
the base-drag treatment see the delivered gap, not the requested one.

### 6.6 NACA generator constants: the ordinate program's precision — Reference note

Not XFOIL: the NASA/PDAS `naca456` program that the geometry generators are gated against
(`tests/fixtures/naca456/`) carries `PI = 3.141592654` in its 6-series mean lines and
`3.14159265` in the φ grid of the 6-series mapping, and inverts its arc-length spline with Brent's
method at a 1e-6 tolerance, reporting the ordinate at the station reached. yFoil uses π and inverts
to round-off; the tolerances in `tests/common/utilities/tolerances.rs` are derived from those three facts.
naca456 also reports the aft slope of the 4-digit modified thickness form as dy/d(1 − x).

---

## 7. Numerical conventions that a reproduction must keep

### 7.1 SETEXP's 1e-5 Newton tolerance — Replicated

Tightening it moves every wake node; ported character for character (`src/solver/xywake.rs:7-9`).

### 7.2 BLSOLV's VACC2/VACC3 association — Replicated

`VACC2 = VACC3 = (VACCEL*2.0)/(S(N)-S(1))` gates the sparse-elimination skips; a 1-ULP difference in the
threshold flips a branch (`src/bl/blsolv.rs:122-124`). A historic "~1 % BLSOLV error" was a yFoil test
hard-coding `S(N)−S(1) = 2.0` instead of the fixture's 2.0387 (`tests/subroutine/blsolv.rs`);
BLSOLV is bit-identical.

### 7.3 Kármán–Tsien TK must be formed as COMSET forms it — Overcome (yFoil bug, fixed)

`TKLAM = MSQ/(1+β)²` and `TKBL = 1/β − 1` are algebraically equal and numerically different; the
difference is dead at M = 0. `BLGlobalParams::new` now uses COMSET's form; two tests derived with the
other form are ignored pending regeneration from the M = 0.3 coverage case.

### 7.4 Integer powers — Replicated (gfortran)

gfortran expands `X**2`, `X**3` inline as multiplications; yFoil writes `x*x*x` where bit-exactness
matters rather than `powi` (CLAUDE.md conventions). Transcendentals come from the host libm in both
codes, so bit-identity is a same-host property: Apple libSystem and glibc differ by 1 ULP on 0.1 % (`exp`,
`ln`, `pow`) to 18 % (`tanh`) of inputs (measured 2026-09-11), and XFOIL's own post-CL_max wanderings
differ by O(1) between the two while its converged points move by ≤ 1.4e-11 (`docs/validation/README.md`,
*Polar break points*; `tests/common/utilities/host.rs`).

### 7.5 EQUIVALENCE aliasing of UNEW/QNEW onto VA/VB — Overcome (structure only)

UPDATE's `UNEW`/`QNEW` alias the factored `VA`/`VB` blocks (`xbl.f`), which are dead after BLSOLV. yFoil
uses separate locals and a consuming `BlsolvInput`, same numerics, so the alias cannot be read by
accident (`src/bl/blsolv.rs:92-94`, `src/solver/update.rs:6-8`).

### 7.6 Single precision and FPE checks — Reference-build facts

Stock XFOIL is single precision and agrees with the double-precision reference only to ~1e-7. The
gfortran Makefile ships with the bounds/FPE checks commented out (`bin/Makefile_gfortran:120`;
`bin/Makefile` has them on). The DP reference built with `-finit-real=snan -ffpe-trap=invalid` runs the
smoke case clean; the only raised (untrapped) flag is `IEEE_DIVIDE_BY_ZERO` from `PLTINI`
(`xplots.f:37`), plot-scale setup that runs even with graphics off, outside the solver.

---

### 7.7 TYPE 2 on a section that can reach negative CL drives itself onto the Kármán–Tsien pole — Replicated

**What TYPE 2 asks XFOIL to do.** TYPE 2 (`MATYP = 2`) is the fixed-lift type: it models an aircraft in
steady level flight, holding weight and chord fixed and letting the speed vary, so that the Mach number
follows the lift as `M∞ = M₁/√CL` (the three types are defined in
[the XFOIL reference](xfoil-reference.md#the-oper-type-command-and-the-meaning-of-matyp-and-retyp)). It therefore does
not fix the Mach number. It asks for a state in which two conditions hold at once: the Mach number follows
the lift, `M∞ = M₁/√CL`, and the lift is whatever the aerofoil produces at that Mach, `CL = CL(α, M∞)`. SPECAL (`xoper.f`) solves this by Newton
iteration on a variable `CLM`, starting from `CLM = 1.0` on every call. That is the whole mechanism, so
the quality of the answer depends entirely on whether the fixed-point problem is well posed at the given α.

**Why the equation is well posed above zero lift and hostile below it.** The boundary is set by the
aerofoil's own compressibility correction. CLCALC forms Cp through the Kármán–Tsien denominator
`DEN = β + BFAC·Cp_inc`, which crosses zero when the local speed reaches `q/Q∞ = √(1 + 2β(1+β)/M∞²)`.
That is a **pole, not a loss of accuracy**: below it Cp is a smooth function of Mach, above it Cp has
changed sign (§4, CPCALC/CLCALC row).

At positive lift the coupling is **self-stabilising** — more lift lowers the Mach number, which weakens the
compressibility correction, which moves the denominator further from zero. At negative lift it is
**self-trapping**, in four steps:

1. `CL ≤ 0` puts the case outside the type's domain of definition — the speed at which a fixed weight is
   carried at negative lift is imaginary — so MRCL's floor `CLA = MAX(CLS, 1e-6)` stands in for it, and M∞
   is pushed to its *maximum* and clamped at 0.99.
2. At M∞ = 0.99 the pole sits at `q/Q∞ = 1.153`, which a 10 % section exceeds by α = ±1°.
3. Past the pole the denominator is negative, so **CL comes back positive at negative α**.
4. A positive CL is perfectly admissible to MRCL, so the fixed-point equation acquires a spurious root.
   Because `CL(α, M)` has a pole it takes every value in a neighbourhood of it, so a root always exists
   *adjacent to the pole*. SPECAL converges onto it in a handful of iterations and reports success.

**Method.** Measured on `naca64a010_n240_inviscid_polar30_type2_m03` (2026-09-17) by exact offline replay of
SPECAL's Newton. The replay is legitimate because on the inviscid path
`GAM(I) = cos α · GAMU(I,1) + sin α · GAMU(I,2)` is Mach-independent — Mach enters only through `BETA` and
`BFAC` inside CLCALC — so `CL(α, M)` can be re-evaluated for any M from the GAM already dumped at that α,
and the whole iteration replayed. The replay reproduces all 61 tracked points with zero relative
disagreement, so the figures below are XFOIL's own arithmetic.

**Every negative-α point is wrong**, not only the conspicuous ones, and most of them converge:

| α | XFOIL reports | truth (M→0) | Newton | min DEN |
|---:|---:|---:|---|---:|
| −30° | **+1.111** | −3.390 | converged, 6 it | −0.949 |
| −20° | **+1.040** | −2.318 | converged, 5 it | −0.0023 |
| −10° | **+0.327** | −1.176 | converged, 4 it | −0.0158 |
| −5° | **+0.128** | −0.590 | converged, 13 it | −0.361 |
| −1° | **−52.086** | −0.118 | failed, 20 it | −0.030 |

The Newton fails outright (`Minf convergence failed`) at only 6 of the 61 points (α = −11, −8, −3, −2, −1, 0).
The converged points are the worse case, because they carry no warning at all: XFOIL reports positive lift on
a symmetric aerofoil at −30° and prints nothing. α = 0 looks innocent only because CL = 0 is pinned by
symmetry whatever the Mach; its Newton failed just as hard, emitting 186 Mach-limit messages. α = +20° and
above is the same failure reached from the other side, once the suction peak crosses the M ≈ 0.3 threshold
of `q/Q∞ = 6.513`.

**The size of the anomaly tracks distance from the pole, not severity of the failure.** At α = −20° the
Newton converged to `M∞ = 0.294270`, where the minimum denominator is −0.0023 — essentially *on* the pole:

```
M = 0.294000   DEN = -0.000435   CL =   +59.89
M = 0.294270   DEN = -0.002316   CL =    +1.07   <- XFOIL's converged answer
M = 0.295000   DEN = -0.007411   CL =  -137.84
```

A change of 1e-3 in M∞ moves CL from +60 to −138, and XFOIL happens to land where it reads +1.04. By the same
token α = −1° reports −52.09 while α = −5° reports +0.128: the first is spectacular only because its
denominator is −0.030 rather than solidly negative, and the second is tamer for being *more* deeply inverted.
Both are equally broken. A polar plotted from this case therefore looks locally plausible exactly where it is
not, and a departure metric built on the size of the anomaly will rank these points in the wrong order.

**These points are ill-conditioned pointwise, not by accumulation** (`docs/conventions/terminology.md`). They reproduce bitwise,
and forward (−30→+30) and reverse (+30→−30) sweeps agree at all 61 points, because SPECAL resets `CLM = 1.0`
at every call and GAMU is α-independent: nothing carries between points on the inviscid path. Their wide
1-ULP twin envelope is pointwise ill-conditioning — each sits next to a pole, so a 1-ULP perturbation
genuinely does move it by O(1). The case is consequently a useful control for sensitivity studies: a wide
envelope that is provably non-cumulative. It says nothing either way about accumulation on the viscous
cases, where BL state really is carried forward.

**yFoil.** Replicated: `set_mach_re_from_cl` (`src/solver/setbl.rs:40-82`) carries both MRCL guards and
`clcalc.rs:57-69` the unguarded denominator, so yFoil converges onto the same spurious roots — and
then flags them, rather than reporting the numbers as if they were sound. This case is the worked
example behind [Solution validity](guide/validity.md), which catalogues every reason yFoil withholds
a result and what each one means. **Gate:** `tests/known_issues/kt_pole_7_7.rs` sweeps the same section at N = 60
(`naca64a010_n60_inviscid_sweep30_type2_m03`) and checks yFoil's validity record against XFOIL's own events
point by point. The measurements above were made at N = 240 (`naca64a010_n240_inviscid_polar30_type2_m03`,
defined in `cases.toml` and regenerated by `cargo xtask fixtures --case` with that name; untracked). The
negative leg of such a sweep must not be read as a polar, and it is excluded from any accuracy claim.

---

### 7.8 UPDATE's largest-change report is a console label that rounding can decide — Replicated

UPDATE records the single largest change it is about to make — its size RMXBL, the variable
VMXBL (`n`/`C`, `T`, `D`, `U`) and the station IMXBL, ISMXBL (`xbl.f:1433-1471`) — and VISCAL only
prints them (`xoper.f:3001`, `:3003`); nothing reads them back. Convergence is `RMSBL < EPS1`, the
root-mean-square change over every variable at every station (`xoper.f:3009`), and the relaxation
factor comes from the changes themselves. Two properties make the label a poor quantity to
compare: at the similarity station the relative changes of θ and δ* are mathematically equal, so
which one is reported is decided by the last bit of rounding; and for Ue the stored size is the raw
change `DUEDG` (`:1468`) while the candidates are ranked by the normalised `DN4`, so the reported
size changes units with the winner. **yFoil** replicates both (`update.rs`, `residual_max`). The
step tests compare the size's magnitude and the solution arrays, not the label; the label is
compared exactly on the well-conditioned reference case (`tests/execution/update.rs`). Not
corrected: yFoil outputs every distribution in full, so nothing depends on the console label.

### 7.9 MRCHDU's station Newton at a transition station can amplify its seed a thousandfold — Replicated

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
undiluted). yFoil follows the same trajectory. **Gate:** `tests/subroutine/mrchdu.rs::
test_mrchdu_station_newton_replays_from_xfoil_state_at_the_7deg_point` (every iterate of the station Newton
from XFOIL's exact state at `TOL_PURE`), `tests/subroutine/transition.rs` (TRCHEK2's every iterate from the
reference's dumped inputs, `trchek2_<k>.dat`), and `tests/execution/steps.rs` (every iteration of the 7° call
replayed from XFOIL's state, each within four times the step's own measured sensitivity). Each code's own
march through the call — the amplification acting on each code's own last-bit seed — is a divergent
comparison in the sense of `docs/conventions/terminology.md` and belongs to the branch-case-polars study. The
instrumented reference dumps `trchek2_<k>.dat` and `mrchdu_trace_<k>.dat` for every `dump_calls` SETBL call.

### 7.10 The fallback tests are written so that a NaN residual takes the fallback — Replicated

MRCHUE and MRCHDU decide between keeping a failed station's last iterate and the garbage
extrapolation with `IF(DMAX .LE. 0.1) GO TO 109` (`xbl.f:785, 1109`): the fallback is taken
whenever the residual is *not* at or below 0.1, which includes a non-finite one. That is how a
station whose closures went non-finite (§4, TRCHEK2 row) is extrapolated over with finite values
and the march carries on. The same holds for every comparison in the reference: NaN makes it
false, and the branch that follows is whichever arm the false case leads to. A reproduction must
therefore keep each comparison's own operator and operand order — `if !(dmax <= 0.1)`, never the
complement `if dmax > 0.1`, which agrees for finite values and disagrees for NaN. yFoil does
(`src/bl/mrchue.rs`, `src/bl/mrchdu.rs`; the negated-comparison lint is allowed for this reason,
`Cargo.toml`). **Gate:** `tests/execution/events.rs::test_mrchue_garbage_extrapolation_events`
and `::test_mrchdu_garbage_extrapolation_events_call_1` (the fallback stations of every non-finite
probe, from XFOIL's own `events.dat`, and the arrays after them, NaN where XFOIL has NaN).

### 7.11 MRCHUE's failure path calls BLVAR(2) and then BLVAR(3) at a wake station — Replicated

After a failed station Newton, MRCHUE recomputes the closures at label 109 with `BLVAR(1)` or
`BLVAR(2)` *and then* `BLVAR(3)` when the station is in the wake (`xbl.f:820-823`); MRCHDU does
the same (`xbl.f:1142-1144`). Each call clamps `HK2` in COMMON (§2.4), so a failed wake station
leaves with the turbulent floor of 1.05 applied before the wake closures are formed, and the
next station's "1" quantities carry that clamp. The effect is visible whenever
1.00005 < Hk < 1.05 at a failed wake station — the first march of the NACA 4412 at M 0.6, α 16°
fails at wake station 40 with Hk = 1.011. yFoil keeps the call sequence in both marches.
**Gate:** `tests/execution/events.rs::test_mrchue_garbage_extrapolation_events` (that
case's first march, station 41 side 2 is where a single call would differ) and
`tests/subroutine/mrchdu.rs` (MRCHDU's failure path on the reference case).

### 7.12 AMI = AMPL2 after a forced transition — Overcome (yFoil bug, fixed)

MRCHUE, MRCHDU and SETBL set `AMI = AMPL2` after every `CALL TRCHEK` (`xbl.f:225`, `:624`, `:815`,
`:968`, `:1137`), whatever TRCHEK2 decided: after a forced transition (`XIFORC` in the interval — with
the default `XSTRIP = 1` that is the trailing-edge interval of a side still laminar there) AMPL2 is the
N2 Newton's converged value. yFoil's `TransitionCheck::Forced` carried no `ampl2`, so its callers kept
the previous AMI. The converged state is unaffected — AMI at a transition station feeds only the
starting value of the next TRCHEK2 call at that station — so no value test saw it; it changed the
number of TRCHEK2 iterates, which `cargo xtask route` found on 2026-10-06 (NACA 0012 at 7° of the
0 → 30° sweep, lower side, station 29: AMPL2 8.645 against 8.830). Fixed by carrying `ampl2` in
`Forced`. **Gate:** `tests/subroutine/mrchdu.rs::
test_mrchdu_station_newton_replays_from_xfoil_state_at_the_7deg_point` (AMPL2 at every station
iterate, both sides).

## 8. Open items in yFoil's own tooling

- `src/bin/generate_subroutine_validation.rs:691-695` still emits an `E24.16` instrumentation
  template (16 significant figures) where the rule is `ES24.16`. Harmless today because nothing
  gates on that generator's output, but it contradicts §1.5.
- Two `#[ignore]`d compressible closure unit tests (§2.1, §7.3) await regeneration from the M = 0.3
  coverage case (`tests/IGNORED.txt`).

- About twenty unit tests in `src/` (`bl/transition.rs`, `bl/difference.rs`, `bl/station.rs`,
  `bl/params.rs`) compare closures with hard-coded values of unknown provenance, three of them
  `#[ignore]`d (`tests/IGNORED.txt`). The same closures are gated at `TOL_PURE` against the reference's
  own log, incompressible and at M = 0.3, by `tests/subroutine/closures.rs`; the unit tests are to be
  moved onto those fixtures (`docs/conventions/testing.md`, Layout).
- `tests/subroutine/pane_legacy.rs` compares the curvature repanelling with panels XFOIL generated
  itself (`tests/fixtures/naca0012/panels*.json`, contrary to CLAUDE.md Rule 4); PANGEN is gated against
  yFoil-generated panels by `tests/subroutine/pangen.rs`.
- One step of the branch candidates is a divergent comparison
  (`xtask/fixtures-config/route.toml`: NACA 63-415 M 0.5, SETBL 2, MRCHDU's garbage stations differ);
  its branches are gated by other steps. The step tests compare a step's outputs, not the route
  through it; observing the route is the next stage (`docs/conventions/testing.md`, "What gated means").
