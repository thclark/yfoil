# XFOIL to YFoil Variable Mapping

This document maps XFOIL Fortran variables and common blocks to their YFoil Rust equivalents.

**Important**: Keep this document updated when modifying struct fields or variable names in either codebase.

---

## 1. Geometry Variables

**XFOIL**: Common block CR05, integer indices in CI04
**YFoil**: `PaneledAirfoil` struct in `src/geometry/airfoil.rs`

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `X(I)` | CR05 | `coords[i].0` | `PaneledAirfoil.coords` |
| `Y(I)` | CR05 | `coords[i].1` | `PaneledAirfoil.coords` |
| `S(I)` | CR05 | `s[i]` | `PaneledAirfoil.s` |
| `NX(I)` | CR05 | `normals[i].0` | `PaneledAirfoil.normals` |
| `NY(I)` | CR05 | `normals[i].1` | `PaneledAirfoil.normals` |
| `APANEL(I)` | CR05 | `theta[i]` | `PaneledAirfoil.theta` |
| `N` | CI04 | `n()` method | `PaneledAirfoil.coords.len()` |
| `SHARP` | CL01 | `sharp` | `PaneledAirfoil.sharp` |
| `XLE, YLE` | CR05 | `le` | `PaneledAirfoil.le` |
| `XTE, YTE` | CR05 | `te` | `PaneledAirfoil.te` |
| `SLE` | CR05 | `s_le` | `PaneledAirfoil.s_le` |
| `CHORD` | CR05 | `chord` | `PaneledAirfoil.chord` |

**Index convention**: XFOIL uses 1-based indexing; YFoil uses 0-based. Convert with `i_yfoil = i_xfoil - 1`.

---

## 2. Inviscid Solution Variables

**XFOIL**: Common blocks CR03 (matrices), CR04 (velocities/Cp), CR06 (circulation)
**YFoil**: `InviscidSolution` struct in `src/panel/solver.rs`

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `AIJ(I,J)` | CR03 | `aij[(i,j)]` | `InviscidSolution.aij` |
| `DIJ(I,J)` | CR03 | `dij[(i,j)]` | `InviscidSolution.dij` |
| `GAM(I)` | CR06 | `gam[i]` | `InviscidSolution.gam` |
| `GAMU(I,1)` | CR06 | `gamu_alpha[i]` | `InviscidSolution.gamu_alpha` |
| `GAMU(I,2)` | CR06 | `gamu_beta[i]` | `InviscidSolution.gamu_beta` |
| `SIG(I)` | CR06 | `sig[i]` | `InviscidSolution.sig` |
| `QINV(I)` | CR04 | `qinv[i]` | `InviscidSolution.qinv` |
| `CPI(I)` | CR04 | `cpi[i]` | `InviscidSolution.cpi` |
| `ALFA` | CR09 | `alpha` | `InviscidSolution.alpha` |
| `CL` | CR09 | `cl` | `InviscidSolution.cl` |
| `CM` | CR09 | `cm` | `InviscidSolution.cm` |

---

## 3. Viscous Solution Variables

**XFOIL**: CR04 (viscous velocities), CR09 (force coefficients)
**YFoil**: `ViscousResult` struct in `src/solver/viscal.rs`

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `QVIS(I)` | CR04 | `qvis[i]` | `ViscousResult.qvis` |
| `CPV(I)` | CR04 | `cpv[i]` | `ViscousResult.cpv` |
| `CD` | CR09 | `cd` | `ViscousResult.cd` |
| `CDF` | CR09 | `cdf` | `ViscousResult.cdf` |
| `CDP` | CR09 | `cdp` | `ViscousResult.cdp` |
| `RMSBL` | local | `rms_bl` | `ViscousResult.rms_bl` |
| `RMXBL` | local | `max_bl` | `ViscousResult.max_bl` |

---

## 4. Boundary Layer State Variables

**XFOIL**: Common block CR15 (BL arrays), CI05 (BL indices)
**YFoil**: `BLStation` in `src/bl/state.rs`, `BLSide` in `src/bl/state.rs`

### Per-Station Variables (indexed by IBL, IS in XFOIL)

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `XSSI(IBL,IS)` | CR15 | `xssi` | `BLStation.xssi` |
| `UEDG(IBL,IS)` | CR15 | `uedg` | `BLStation.uedg` |
| `DSTR(IBL,IS)` | CR15 | `dstr` | `BLStation.dstr` |
| `THET(IBL,IS)` | CR15 | `thet` | `BLStation.thet` |
| `CTAU(IBL,IS)` | CR15 | `ctau` | `BLStation.ctau` |
| `MASS(IBL,IS)` | CR15 | `mass` | `BLStation.mass` |
| `TAU(IBL,IS)` | CR15 | `tau` | `BLStation.tau` |
| `DIS(IBL,IS)` | CR15 | `dis` | `BLStation.dis` |
| `CTQ(IBL,IS)` | CR15 | `ctq` | `BLStation.ctq` |
| `TSTR(IBL,IS)` | CR15 | `tstr` | `BLStation.tstr` |
| `DELT(IBL,IS)` | CR15 | `delt` | `BLStation.delt` |
| `ENTR(IBL,IS)` | CR15 | `amplification` | `BLStation.amplification` |

### Side-Level Variables

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `NBL(IS)` | CI05 | `stations.len()` | `BLSide.stations` |
| `ITRAN(IS)` | CI05 | `transition_index` | `BLSide.transition_index` |
| `IBLTE(IS)` | CI05 | `i_te` | `BLSide.i_te` |
| `XOCTR(IS)` | local | `transition_x` | `BLSide.transition_x` |
| `XSSITR(IS)` | local | `transition_xssi` | `BLSide.transition_xssi` |

### Global BL Variables

| XFOIL Variable | XFOIL Location | YFoil Variable | YFoil Location |
|----------------|----------------|----------------|----------------|
| `IST` | CI05 | `i_stag` | `BLState.i_stag` |
| `SIMI` | CL01 | (similarity flag) | `BLState.use_similarity` |
| `TRAN` | CL01 | (transition flag) | computed from `transition_index` |

---

## 5. BL Index Mapping Arrays

**XFOIL**: CI05, local arrays in xbl.f
**YFoil**: Fields in `BLStation` or computed on-the-fly

| XFOIL Variable | Purpose | YFoil Equivalent |
|----------------|---------|------------------|
| `IPAN(IBL,IS)` | BL station → panel index | `BLStation.ipan` |
| `ISYS(IBL,IS)` | BL station → Newton system row | Computed in system assembly |
| `VTI(IBL,IS)` | Velocity sign (+1 upper, -1 lower) | `BLStation.vti` |

---

## 6. BL Closure Variables

**XFOIL**: Local variables in BLVAR, HKIN, HSL/HST, CFL/CFT, DIL/DIT
**YFoil**: `ClosureResult` in `src/bl/closure.rs`, `BLStationState` in `src/bl/system.rs`

### Primary Closure Outputs

| XFOIL Variable | XFOIL Subroutine | YFoil Variable | YFoil Location |
|----------------|------------------|----------------|----------------|
| `HK` | HKIN | `hk` | `ClosureResult.hk` |
| `HS` | HSL/HST | `hs` | `ClosureResult.hs` |
| `CF` | CFL/CFT | `cf` | `ClosureResult.cf` |
| `CD` | DIL/DIT | `cd` | `ClosureResult.cd` |
| `US` | USL | `us` | `ClosureResult.us` |

### Derived Variables

| XFOIL Variable | XFOIL Context | YFoil Variable | YFoil Location |
|----------------|---------------|----------------|----------------|
| `RT` | BLVAR (Rθ) | `rt` | `BLStationState.rt` |
| `MSQ` | BLVAR (M²) | `msq` | `BLStationState.msq` |
| `AMPL` | DAMPL | `amplification` | computed in `transition.rs` |

### Closure Derivatives

YFoil and XFOIL use the same derivative structure with chain rule conversion:

**Step 1 - Closure functions** return derivatives w.r.t. intermediate variables:
- `hk_h`, `hk_msq` from `hkin()`
- `hs_hk`, `hs_rt`, `hs_msq` from `hs_lam()`/`hs_turb()`
- `cf_hk`, `cf_rt`, `cf_msq` from `cf_lam()`/`cf_turb()`
- `di_hk`, `di_rt` from `di_lam()`

**Step 2 - BLKIN equivalent** (`BLStationState::blkin()`) computes intermediate derivatives:
- `hk_u`, `hk_t`, `hk_d` = ∂Hk/∂U, ∂Hk/∂θ, ∂Hk/∂δ*
- `rt_u`, `rt_t` = ∂Rθ/∂U, ∂Rθ/∂θ
- `msq_u` = ∂M²/∂U

**Step 3 - BLVAR equivalent** (`BLStationState::blvar()`) applies chain rule:
- `hs_u = hs_hk * hk_u + hs_rt * rt_u + hs_msq * msq_u`
- `cf_u = cf_hk * hk_u + cf_rt * rt_u + cf_msq * msq_u`
- etc.

This matches XFOIL's BLKIN + BLVAR approach exactly.

---

## 7. BL System Variables

**XFOIL**: Local arrays in BLSYS
**YFoil**: `BLSystem` struct (or equivalent matrices)

| XFOIL Variable | Dimensions | YFoil Variable | YFoil Location |
|----------------|------------|----------------|----------------|
| `VS1(4,5)` | 4×5 | `vs1` | BL system matrix |
| `VS2(4,5)` | 4×5 | `vs2` | BL system matrix |
| `VSREZ(4)` | 4 | `vsrez` | BL residual vector |
| `VSM(4)` | 4 | `vsm` | Mass equation coeffs |
| `VSR(4)` | 4 | `vsr` | Shape equation coeffs |
| `VSX(4)` | 4 | `vsx` | Ctau equation coeffs |

---

## 8. Solver Configuration

**XFOIL**: Various globals in XFOIL.INC
**YFoil**: `ViscalConfig` in `src/solver/viscal.rs`

| XFOIL Variable | YFoil Variable | Purpose |
|----------------|----------------|---------|
| `ACRIT` | `n_crit` | Critical amplification factor |
| `VACCEL` | `vaccel` | BL solution acceleration |
| `REINF` | `re` | Reynolds number |
| `MINF` | `mach` | Mach number |
| `RLXBL` | (internal) | Under-relaxation factor |

---

## Notes on Structural Differences

### Wake Handling
- **XFOIL**: Wake stations are appended to lower surface (side 2) after `IBLTE(2)`
- **YFoil**: Wake is stored separately in `BLState.wake: Vec<BLStation>`

### Transition Representation
- **XFOIL**: `ITRAN(IS)` is an integer; 0 means no transition yet
- **YFoil**: `transition_index` is `Option<usize>`; `None` means no transition

### Matrix Storage
- **XFOIL**: Dense 2D arrays with fixed maximum dimensions
- **YFoil**: `nalgebra::DMatrix` with dynamic sizing

## BL state (`src/solver/blstate.rs`) — stages S2 onward

`BlState` mirrors the BL COMMON blocks one-to-one and keeps XFOIL's indexing: **1-based with a
dummy slot 0**, sides `is = 1, 2`, wake appended to side 2 (`NBL(2) = IBLTE(2) + NW`). The field
names *are* the Fortran names, lower-cased, so the table is only the exceptions and shapes.

| XFOIL | YFoil | Notes |
|---|---|---|
| `N`, `NW` | `st.n`, `st.nw` | |
| `X(I)`, `Y(I)`, `S(I)` for `I = 1..N+NW` | `st.x[i]`, `st.y[i]`, `st.s[i]` | wake nodes `n+1..=n+nw` |
| `XP(I)`, `YP(I)` (airfoil spline derivatives) | `st.xp[i]`, `st.yp[i]` | `1..=n` |
| `GAM(I)`, `GAM_A(I)` | `st.gam[i]`, `st.gam_a[i]` | |
| `QINVU(I,1)`, `QINVU(I,2)` | `st.qinvu[1][i]`, `st.qinvu[2][i]` | |
| `QINV(I)`, `QINV_A(I)`, `QVIS(I)` | `st.qinv[i]`, `st.qinv_a[i]`, `st.qvis[i]` | |
| `CHORD, SLE, XLE, YLE, XTE, YTE` | same names | |
| `ANTE, ASTE, DSTE, SHARP` (TECALC) | same names | `pointers::tecalc` |
| `IST, SST, SST_GO, SST_GP` (STFIND) | same names | `pointers::stfind` |
| `NBL(IS)`, `IBLTE(IS)`, `ITRAN(IS)`, `NSYS` | `st.nbl[is]`, `st.iblte[is]`, `st.itran[is]`, `st.nsys` | `[_; 3]`, index 0 unused |
| `IPAN(IBL,IS)`, `VTI(IBL,IS)`, `ISYS(IBL,IS)` | `st.ipan[is][ibl]`, `st.vti[is][ibl]`, `st.isys[is][ibl]` | `pointers::{iblpan, iblsys}` |
| `XSSI, UEDG, UINV, UINV_A, MASS, THET, DSTR, CTAU, DELT, TSTR, USLP, GUXQ, GUXD, TAU, DIS, CTQ (IBL,IS)` | `st.<name>[is][ibl]` | |
| `WGAP(IW)` | `st.wgap[iw]` | `1..=nw`, set by `pointers::xicalc` |
| `XSTRIP(IS)` | `st.xstrip[is]` | `XIFORC` is returned by `pointers::xifset(&st, is)` |

Velocity-layer subroutines (`src/solver/velocity.rs`): `QISET → qiset`, `UICALC → uicalc`,
`UECALC → uecalc`, `QVFUE → qvfue`, `GAMQV → gamqv`, `UESET → ueset` (takes the
`(N+NW)×(N+NW)` DIJ, 0-based storage), `DSSET → dsset`.

The legacy `SetblState` in `src/solver/setbl.rs` (two 0-based sides, no wake stations) is
superseded by `BlState` and is deleted with the legacy VISCAL (stage S9).

Inviscid / wake subroutines on `BlState` (stages S3–S4): `PSILIN → solver::psilin::psilin`
(returns `Psilin { psi, psi_ni, qtan1, qtan2, qtanm, dzdg, dqdg, dzdm, dqdm, z_qinf, z_alfa }`),
`PSWLIN → solver::qdcalc::pswlin`, `SETEXP/XYWAKE/QWCALC → solver::xywake`, `GGCALC → solver::ggcalc::ggcalc`
(returns `InviscidSystem { aij: LuFactors, bij, ladij }`), `LUDCMP/BAKSUB → solver::ludcmp`,
`QDCALC → solver::qdcalc::qdcalc` (fills `st.dij`, 1-based (N+NW)²), `ATANC → solver::ggcalc::atanc`.
`PI/HOPI/QOPI` are computed as XFOIL's INIT does (`4*atan(1)`), see `psilin::pi_consts`.

BL march subroutines (stage S5): `MRCHUE → bl::mrchue::mrchue(&mut st, &params, acrit, trace)`
(fills `THET/DSTR/CTAU/UEDG/MASS/TAU/DIS/CTQ/DELT/TSTR`, sets `ITRAN`; `MrchueTrace` mirrors the
`xfoil_newton_trace.dat` records: `StationIter { ampl, primary, kinematic, closure, residual, vs2,
solution, dmax, rlx, updated, converged }`), `BLSYS → bl::blsys::blsys(sys, s1, s2, IntervalFlags
{simi, tran, turb, wake}, trans, acrit, params)`, `TESYS → bl::blsys::tesys`, `TRCHEK2 →
bl::system::trchek` (returns `TransitionResult::{NoTransition{ampl2}, FreeTransition{location,
ampl2}, ForcedTransition{location}}` with `TransitionLocation` = `XT` and its `XT_*`
sensitivities; `ampl2` is the iterated `AMPL2`, which may exceed `AMCRIT` exactly as XFOIL leaves
it), `BLVAR/BLKIN/BLPRV → BLStationState::{blvar, blkin, blprv}`, `BLDIF → BLLocalSystem::bldif`,
`TRDIF → BLLocalSystem::trdif`, `BLMID → MidpointCf::compute`, `DILW → bl::closure::dilw`.
`COM1/COM2` are `s1: BLStationState` / `s2: BLStationState`; the `COM1 = COM2` copies are
`s1 = s2.clone()`. XFOIL quirks reproduced: `HVRAT` is never assigned on the analysis path
(`BLGlobalParams.hvrat = 0.0`); `BLVAR` clamps `HK2` in COMMON without recomputing its
derivatives (`blvar` writes the clamped `hk` back); `BLDIF` forms `UQ_T1..UQ_RE` but never uses them.

Stage S6: `MRCHDU → bl::mrchdu::mrchdu(&mut st, &params, acrit, trace)` (mixed-mode march on the
Ue–Hk characteristic: `SENSWT`, `UEREF/HKREF`, the `ITROLD` re-laminarisation/re-turbulisation
logic, the 25-iteration Newton with `DEPS = 5e-6`, `DSLIM`, and the extrapolation fallback;
`MrchduTrace` mirrors `xfoil_mrchdu_trace.dat`). XFOIL's `COM1`, `COM2` and `XT` COMMON state
persists across MRCHUE/MRCHDU/SETBL calls, so they live on `BlState` as `com1`, `com2`, `xt` and
each march takes them out and puts them back (`std::mem::take`). `XSSITR(IS)`, `TFORCE(IS)` →
`st.xssitr[is]`, `st.tforce[is]`. The pre-S6 station-at-a-time march is `bl::march_legacy`
(used only by the legacy `SetblState` path; both go with the S7 SETBL rewrite).

Stage S7: `SETBL → solver::setbl::setbl(&mut st) -> SetblResult { sys: BlsolvInput, params,
re_clmr, msq_clmr, dule }` — MRCL/COMSET/parameter setup, MRCHUE (if `!st.lblini`), MRCHDU, the
USAV/UESET swap, ULE/UTE sensitivities from DIJ, the assembly sweep with the full VM chain rule,
the VDEL Re/Mach column, the VZ block, TAU/DIS/CTQ/DELT/USLP, XOCTR/YOCTR/TINDEX. `MRCL →
solver::setbl::mrcl`. The XFOIL.INC controls it reads live on `BlState`: `LALFA, CL, CLSPEC,
MATYP, RETYP, MINF1/REINF1, MINF/REINF, LBLINI, ACRIT(IS), VACCEL, GAMMA`. `XT` and the `XT_*`
sensitivities (TRCHEK2's COMMON outputs) are `st.trloc: TransitionLocation`, taken out and put
back by every march like `com1`/`com2`. `IDAMPV` is pinned at 0. `BLGlobalParams::new` forms
TKLAM exactly as COMSET (was `1/beta - 1`; dead at M = 0). The legacy `SetblState` path is
`solver::setbl_legacy` (with `bl::march_legacy`, `bl::wake`), used only by the legacy VISCAL and
deleted with the S9 rewrite.

Stage S8: `UPDATE → solver::update::update(&mut st, &vdel, minf_cl) -> UpdateResult { rlx,
rmsbl, rmxbl, vmxbl, imxbl, ismxbl, dac, clnew, cl_a, cl_ms, cl_ac }` — UNEW/U_AC from DIJ and
the solved VDEL, QNEW/Q_AC, the Kármán–Tsien CLNEW integral with CL_A/CL_MS/CL_AC, DAC for CL
(LALFA) or alpha, the two-pass RLX search (DHI = 1.5, DLO = −0.5, DCLMAX/DALMAX), RMSBL/RMXBL,
the under-relaxed update with the CTAU ≤ 0.25 clamp for IBL ≥ ITRAN, DSLIM, MASS = DSTR·UEDG,
the negative-Ue island fix-up and the wake array equating. `UNEW/U_AC/QNEW/Q_AC` are locals (XFOIL
EQUIVALENCEs them onto VA/VB; `blsolv` consumes its input, so nothing can read VA/VB after the
solve). `CL`/`ALFA` updated in place on `BlState`.

Stage S9: `VISCAL → solver::viscal::viscal(&mut st, sys, niter, waklen, trace) -> bool`
(prologue XYWAKE/QWCALC/QISET/STFIND/IBLPAN/XICALC/IBLSYS/UICALC/QDCALC guarded by
`LWAKE/LIPAN/LBLINI/LWDIJ`, then SETBL → BLSOLV → UPDATE → MRCL+COMSET | QISET+UICALC → QVFUE →
GAMQV → STMOVE → CLCALC → CDCALC to `RMSBL < EPS1`; `ViscalIter` mirrors `viscal_iter.dat`).
`STMOVE → solver::pointers::stmove`, `CPCALC/CLCALC/CDCALC/COMSET → solver::clcalc::{cpcalc,
clcalc, cdcalc, comset}`, `SPECAL → solver::specal::specal` and OPER's `ALFA` command →
`specal::alfa_command` (SPECAL, then the AWAKE/AVISC/MVISC invalidations), OPER `INIT` →
`analysis::Session::init` (LBLINI toggle). The XFOIL.INC flags/outputs on `BlState`: `LWAKE,
LIPAN, LWDIJ, LVISC, LVCONV, AWAKE, AVISC, MVISC, TKLAM, TKL_MSQ, MINF_CL, REINF_CL, CM, CDP, CD,
CDF, CL_ALF, CL_MSQ, XCMREF, YCMREF, CPI, CPV`; `LADIJ` stays on `InviscidSystem`; `LGAMU/LQAIJ`
are the presence of the `InviscidSystem`. XYWAKE sets `LWAKE/AWAKE/LWDIJ`, IBLPAN `LIPAN`,
QDCALC `LWDIJ`, exactly where XFOIL does. `solver::analysis::{FlowSpec, Session, analyze,
compute_polar}` is the OPER driver (one persistent `BlState` per session, as XFOIL's COMMON).
The legacy `SetblState`/station-at-a-time/`bl::wake`/`bl::newton`/`forces` implementations are
deleted; `panel::solve_inviscid` remains only for the pre-S3 inviscid tests (removed with
stage G).

Stage S10: the polar is XFOIL's OPER script run through one persistent `Session`
(`analysis::compute_polar`): `ALFA 0` (`Session::alfa`, VISCAL(ITMAX)), `ASEQ step alpha_max step`
(`Session::aseq` → `specal::aseq_point` — invalidations, SPECAL — then VISCAL(ITMAX + 5), as ASEQ
does), `INIT` (`Session::init`, the LBLINI toggle that also clears LIPAN), `ALFA -step`, `ASEQ`
down to alpha_min; each ASEQ halts after `NSEQEX` (4) consecutive non-converged points, PACC keeps
only converged points, and the result is stitched ascending. `QDCALC` keeps the airfoil DIJ block
across calls (LADIJ) and refreshes only the wake part when XYWAKE has moved the wake (LWDIJ).

Stage S11: `SPECCL → solver::specal::speccl` (MRCL(CLSPEC)+COMSET, then the 20-iteration Newton on
ALFA with CL_ALF to |DALFA| ≤ 1e-6, TECALC, QISET, Cp), OPER's `CL` command →
`specal::cl_command` (LALFA = .FALSE., ALFA = 0 as the initial guess, SPECCL, invalidations) →
`Session::cl`; VISCAL then runs its `QISET+UICALC` branch and UPDATE moves ALFA (DAC with
DALMAX/DALMIN). `MRCL` with MATYP/RETYP = 2 (Re, M ∝ 1/√CL, re-evaluated from CL after every
UPDATE) is gated by the `naca0012_n60_a2_re1e6_type2` case (`FlowSpec { matyp, retyp }`; the
pipeline emits OPER `TYPE n`). Minimal cases (`minimal = true`) keep only the `viscal_*.dat`
records (~50 KB each) so coverage cases stay cheap to track.

Coverage cases and the replay harness: `cases.toml` options `airfoil = "naca4:0012:sharp"` (`geom
naca --sharp`, `Geometry::sharpen`), `xtr = [xu, xl]` (OPER `VPAR`/`XTR`), `dump_calls = [k, ...]`
(the SETBL/UPDATE per-call dumps for arbitrary iterations; `RDDUMP` reads `dump_calls.txt`), and
the +1-ULP twin run behind every case (`noise_floor.json`, `xtask::ulp_twin`). `tests/utilities/
records.rs` is the VISCAL-level checker (exact branch trace, floor-derived value tolerances, the
threshold-straddling outcome); `tests/xfoil_coverage_tests.rs::replay_iteration` seeds a
`Session` with XFOIL's dumped state at call k (pointers rebuilt from IST/SST) and runs one
SETBL → BLSOLV → UPDATE against `update_output_k.dat`.

Stage G (optional, off the solver's critical path): `NACA4/NACA5 (naca.f) → geometry::naca::{naca_4digit_xfoil,
naca_5digit_xfoil}` (NSIDE = IQX/3 = 123, AN = 1.5 spacing, thickness applied *vertically* — the
documented divergence from the NACA definition that `naca_4digit` does not share; 245-point buffer,
no panel count), `PANGEN → geometry::panel::repanel_xfoil` line for line (IPFAC = 5, RDSTE = 0.667,
the LE-adjacent smoothing equations on both neighbours, exact-equality corner/LE tests, TRISOL
split at a sharp LE, corner insertion), with `SCALC/SEGSPL/CURV/LEFIND/TRISOL →
geometry::panel::{scalc, segspl, curv, lefind, trisol}` as XFOIL has them (`create_paneled_airfoil`
uses the same `lefind`). CLI: `yfoil geometry naca --naca-model xfoil`. The pipeline's geometry-only
cases (`airfoil = "xfoil-naca:dddd"`, `geometry_only = true`) run `NACA dddd / PPAR / N n` and dump
`xfoil_pangen.dat` from PANGEN. `panel::solve_inviscid` (the pre-S3 inviscid solver) is deleted;
GGCALC/PSILIN on `BlState` are the inviscid solve.

## OPER DAMP (2026-09-04)

| XFOIL | YFoil | Notes |
|---|---|---|
| `IDAMP` (OPER `DAMP` toggle) | `BlState.idamp`, `FlowSpec.idamp` | 0/false = envelope e^n (`DAMPL`), 1/true = modified envelope (`DAMPL2`) |
| `IDAMPV` (SETBL: `IDAMPV = IDAMP`) | `BLGlobalParams.idampv` | read by `axset` (AXSET's `IF(IDAMPV.EQ.0)`), reached through `trchek` and `bldif` |
| `DAMPL2` | `bl::system::dampl2` | verbatim; gated end-to-end by `naca0012_n60_a2_re1e6_damp` |

## Output: BLDUMP ↔ `BlSideOutput` (2026-09-04)

`yfoil analyze -o` / `yfoil polar --distributions` write `AnalysisOutput.boundary_layer`
(`src/output/foil.rs`): three struct-of-arrays sides (`upper` = side 1 to IBLTE, `lower` = side 2
to IBLTE, `wake` = side 2 past IBLTE), each column one entry per station.

| XFOIL DUMP column (`xoper.f:1955-1995`) | XFOIL source | `BlSideOutput` column | How YFoil forms it |
|---|---|---|---|
| `s`, `x`, `y` | `S(I), X(I), Y(I)` | `x`, `y` (`xssi` is XSSI, not S) | verbatim node coordinates |
| `Ue/Vinf` | `(GAM/QINF)(1-TKLAM)/(1-TKLAM(GAM/QINF)²)` — signed by GAM, negative on the lower side | `ue` | BLPRV's `U2/QINF` on UEDG: the same transformation, unsigned |
| `Dstar`, `Theta` | `DSTR, THET` (post-UPDATE) | `dstr`, `thet` | verbatim |
| `H` | `DSTR/THET` | `h` | BLKIN on the primaries (wake: without WGAP, as BLKIN) |
| `HK` | `HKIN(H, AMSQ)` live | `hk` | BLKIN's HKIN, unclamped |
| `Cf` | `TAU/(½QINF²)`, TAU **lagged** | `cf` (live), `stored.cf_dump` (DUMP's) | live: BLVAR's `CF2` scaled `R2·U2²/QINF²` |
| `H*` | `TSTR/THET`, TSTR **lagged** | `hs` (live), `stored.hs_dump` (DUMP's) | live: BLVAR's `HS2` |
| `CDIS` (`Di` in the wide format) | `DIS/QINF³`, DIS **lagged** | `cdis` (live), `stored.dis` | live: `R2·U2³·DI2·HS2·½/QINF³` |
| — | `DELT, CTQ, USLP` (VPLO) **lagged** | `delta`, `ctq`, `uslp` (live), `stored.delt/ctq/uslp` | live: BLVAR's `DE2, CQ2, 1.6/(1+US2)` |
| — | `CTAU, MASS` | `ctau`, `mass` | verbatim |
| — | `CPV(I)` | `cp` | verbatim |

**Why "lagged".** SETBL calls MRCHDU at the top of every Newton iteration (`xbl.f:93`); MRCHDU
and SETBL store `TAU, DIS, CTQ, DELT, USLP, TSTR` from the state entering that iteration
(`xbl.f:277-282, 1157-1167`), then BLSOLV/UPDATE correct `THET, DSTR, UEDG, CTAU, MASS` and nothing
refreshes the closure arrays. So DUMP's `Cf`, `H*`, `K`, `tau`, `Di` and VPLO's `CF`, `CD`, `DELT`
are one Newton correction behind `Dstar`, `Theta`, `Ue`, `H`, `HK`, `N`, `CT` (`H*` and `K` divide a
lagged TSTR by a current THET). The mismatch is bounded by the convergence test (RMSBL < 1e-4) on a
converged point and unbounded on an unconverged one. It is an inconsistency in XFOIL's own output,
not a solver error; `blplot.f:1387` has the `TSTR/THET` plot commented out in favour of a live
HKIN. YFoil stores the same lagged arrays (fixture-gated: `tests/xfoil_mrchdu_tests.rs`,
`tests/xfoil_update_tests.rs`) and emits them under `stored`; its canonical columns are the closures
re-evaluated live on the converged primaries. A side-by-side of live `hs` against DUMP's `H*`
therefore differs at the ~RMSBL level (`tests/foil_output_tests.rs` bounds it).

**Markers** (`BoundaryLayerOutput`): `stagnation` = `IST, SST` and the spline point at SST;
`transition[is]` = `ITRAN, TFORCE, XOCTR, YOCTR` and the spline point at `SST ∓ XSSITR`;
`wake_split` = CPDISP's `DSF1, DSF2` (`xplots.f:714-721`); `derived_separation` is YFoil's (XFOIL
reports no separation location): the sign change of the live `cf` between consecutive surface
stations, linear in S, on the spline.
