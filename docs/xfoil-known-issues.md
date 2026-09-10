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

**yFoil.** The stored arrays are reproduced exactly (`tests/xfoil_mrchdu_tests.rs:60`,
`tests/xfoil_update_tests.rs:75`, the SETBL gates); the solver is untouched. In the analysis JSON
(`BlSideOutput`, `src/output/foil.rs`) the canonical columns `hs, cf, cdis, delta, ctq, uslp` are
computed *live* by running XFOIL's own BLPRV → BLKIN → BLVAR on the converged primaries, and the lagged
arrays are emitted verbatim under `stored` (with `hs_dump`, `cf_dump` as DUMP prints them). A
side-by-side of live `hs` against DUMP's `H*` will differ at the ~RMSBL level by construction;
`tests/foil_output_tests.rs` bounds it. Column table: `docs/xfoil-reference.md`,
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
the previous iteration's value; the gate excludes it (`tests/xfoil_mrchue_tests.rs:123`).

### 2.7 RADLE stays zero for a flat LE — Replicated

`xgeom.f:350`: LE curvature below `0.001/(S(N)-S(1))` leaves `RADLE = 0` (flat plate, wedge). Open
branch with a reach note in `coverage.toml`.

### 2.8 MASS is never written on side 1's wake slots — Replicated (harness note)

UPDATE's final loop equates the "upper wake arrays" `(IBLTE(1)+K, 1)` to side 2's wake for CTAU,
THET, DSTR, UEDG, TAU, DIS, CTQ, DELT and TSTR (`xbl.f:1546-1557`) but not MASS, and nothing else
writes `MASS(IBL, 1)` beyond `NBL(1)` (MRCHUE/MRCHDU, UPDATE, DSSET and the SETBL/QVFUE sums all
run to `NBL(IS)`, and `NBL(1) = IBLTE(1)`). Those slots hold whatever an earlier stagnation-point
position left there and are read by nothing. yFoil's `mass_defect[1][..]` beyond `n_stations[1]` is
the same dead storage with a different history, so the replay harness (`tests/utilities/replay.rs`)
does not compare it; it surfaced in the polar break-point cases, where IST moves between iterations
and the instrumented `update_output_<k>.dat` dumps every row up to `NBL(1) + NW`.

---

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
| `xfoil.f:1837` (PANGEN) | `temporarily used for more reliable convergence` — IPFAC = 5 oversampling that was never removed | Replicated, `repanel_xfoil` |
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
| `xfoil.f:845, 856` (MRCL) | `CL too low for chosen Mach(CL) dependence — artificially limiting Mach to 0.99`; Re limited to 100× REINF1 | `setbl.rs:77, 86`; reach notes in `coverage.toml` |
| `xutils.f:41-50` (SETEXP) | 100-iteration Newton on the spacing ratio to `|dRatio| < 1e-5`; `Convergence failed. Continuing anyway ...` | `src/solver/xywake.rs:33-49`; the 1e-5 tolerance is load-bearing (§7.1) |
| `xoper.f` (VISCAL, SPECAL, SPECCL) | `Convergence failed` after ITMAX / 20 / 12 iterations; state kept | `viscal.rs:187`, `specal.rs:110, 172` |
| `xbl.f:1474` (UPDATE) | under-relaxation whenever a step would change a BL variable by more than 50 % (DLO) | `src/solver/update.rs` |
| `xbl.f:1515-1517` (UPDATE) | `eliminate absurd transients`: CTAU capped at 0.25 for turbulent stations only | `update.rs:296` |
| `xbl.f:1537` (UPDATE) | `make sure there are no "islands" of negative Ue` | `update.rs:311` |
| `xbl.f:785-827`, `:1111-1149` | MRCHUE/MRCHDU non-convergence fallbacks: extrapolate from the previous station when the residual > 0.1 (`garbage solution`) | Translated; still unexercised (every failure so far had DMAX ≤ 0.1), reach notes in `coverage.toml` |
| `xblsys.f` (BLVAR) | laminar Cf used in the turbulent branch at unreasonably small Rθ | `src/bl/system.rs:1287, 1408` |
| `xblsys.f` (BLVAR) | `DE` capped at `12·θ` with all four sensitivities zeroed — a Jacobian discontinuity | `system.rs:1458-1465` |
| `xblsys.f` (BLDIF) | similarity-station amplification row is a dummy `VS2(1,1) = 1` pivot | `system.rs:1744-1746` |
| `xblsys.f` (BLDIF) | forms `UQ_RTA`, `UQ_T1..UQ_RE` and never uses them | omitted where they feed nothing, `system.rs:1830` |
| `xpanel.f:1737` (UESET) | `tweak Ue so it's not zero, in case stag. point is right on node` (UEPS = 1e-7) | `pointers.rs:325` |
| `xpanel.f:1385` (STFIND) | `tweak stagnation point if it falls right on a node (very unlikely)` — needs `GAM(I) == 0` bitwise | `pointers.rs:47`; threshold-straddling territory |
| `xbl.f` (SETBL) | `SETBL: Xtr???  n1 n2:` diagnostic when the transition interval and ITRAN disagree | `setbl.rs:272` |
| `plotlib/plt_font.f:249` (PLNUMBABS, via ASEQ → SEQPLT at `xoper.f:719` and via ALFA → CPX → COEFPL) | plot labels are formatted even with graphics off (`PLOP / G F`), and PLNUMBABS's digit-extraction loop never terminates on a non-finite value. Reached by both the sequence-plot label and the Cp-plot coefficient label, so every OPER point whose CL/CD has become Infinity hangs XFOIL at 100 % CPU (NACA 0012, `ITER 100`, past ~21° on either leg; stacks sampled 2026-09-10). The `xfoil-sensitivity` driver detects the stall (stdout stops growing) and kills the run | Out of scope (plot library); the sequence data up to the hang is intact. yFoil has no plot label and halts the sequence normally (`compute_polar`) |

---

## 5. Dead and unreachable code in 6.99

### 5.1 Subroutines never called — Out of scope / translated for completeness

`UECALC` and `DSSET` are never called anywhere in 6.99 (translated anyway,
`src/solver/velocity.rs:29, 77`); `DIT`'s only call site is commented out (`xblsys.f:944`, not
translated). `docs/validation/coverage.md` lists them as the three never-called subroutines.

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
(`tests/xfoil_coverage_tests.rs:290`).

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
default `naca_4digit`/`naca_5digit` use the NACA definition; `--naca-model xfoil` reproduces XFOIL's
generator bitwise. This never enters a comparison because yFoil generates the panels and XFOIL
consumes them (Rule 4). Recorded in CLAUDE.md's divergence table, the only row.

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

---

## 7. Numerical conventions that a reproduction must keep

### 7.1 SETEXP's 1e-5 Newton tolerance — Replicated

Tightening it moves every wake node; ported character for character (`src/solver/xywake.rs:7-9`).

### 7.2 BLSOLV's VACC2/VACC3 association — Replicated

`VACC2 = VACC3 = (VACCEL*2.0)/(S(N)-S(1))` gates the sparse-elimination skips; a 1-ULP difference in the
threshold flips a branch (`src/bl/blsolv.rs:122-124`). A historic "~1 % BLSOLV error" was a yFoil test
hard-coding `S(N)−S(1) = 2.0` instead of the fixture's 2.0387 (`tests/xfoil_blsolv_tests.rs:9-12`);
BLSOLV is bit-identical.

### 7.3 Kármán–Tsien TK must be formed as COMSET forms it — Overcome (yFoil bug, fixed)

`TKLAM = MSQ/(1+β)²` and `TKBL = 1/β − 1` are algebraically equal and numerically different; the
difference is dead at M = 0. `BLGlobalParams::new` now uses COMSET's form; two tests derived with the
other form are ignored pending regeneration from the M = 0.3 coverage case.

### 7.4 Integer powers — Replicated (gfortran)

gfortran expands `X**2`, `X**3` inline as multiplications; yFoil writes `x*x*x` where bit-exactness
matters rather than `powi` (CLAUDE.md conventions). Transcendentals come from the host libm in both
codes, so bit-identity is a same-host property.

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

## 8. Open items in yFoil's own tooling

- `src/bin/generate_subroutine_validation.rs:691-695` still emits an `E24.16` instrumentation
  template (16 significant figures) where the rule is `ES24.16`. Harmless today because nothing
  gates on that generator's output, but it contradicts §1.5.
- Two `#[ignore]`d compressible closure unit tests (§2.1, §7.3) await regeneration from the M = 0.3
  coverage case (`tests/IGNORED.txt`).
