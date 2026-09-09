# Naming

How variables, fields, functions, types, JSON keys and CLI flags are named in yFoil, and why.
This document holds the *rules*. The *application* of the rules, the table of every XFOIL
name and the yFoil name it became, lives in
[`docs/xfoil-reference/xfoil-to-yfoil-mapping.md`](../xfoil-reference/xfoil-to-yfoil-mapping.md);
read that when you need to know what a particular XFOIL variable is called here, and read this
when you need to name something new.

## Why the names were remapped

yFoil is a line-for-line translation of XFOIL. The first translation kept XFOIL's identifiers,
lower-cased: `thet`, `dstr`, `tstr`, `sst_go`, `vsrez`, `iblte`. Those names come from an era
when identifiers were six characters and every byte counted; they are opaque to anyone who does
not have `XFOIL.INC` open, and several of them are wrong about what they hold (`CTAU` holds
√Cτ, `M2` is a Mach number squared, `RMSBL` is a convergence metric). We are moving away from
XFOIL's names, not towards them. The translation stays provable because every XFOIL name is
recorded once, in the mapping table and in a `#[doc(alias)]` on the translating function, and
never inside an identifier.

The opposite failure is also a failure: spelling a symbol out into English
(`kinetic_energy_thickness` for θ*) is no clearer to a reader of Drela's papers than `tstr`
was. Symbols from the equation system are names.

## The rules

1. **Symbols from the equation system are names.** `theta`, `thetastar`, `dstar`, `h`, `hk`,
   `hstar`, `hstarstar`, `cf`, `sqrtctau`, `cdiss`, `us`, `ue`, `q`, `qinf`, `psi`, `xi`,
   `gamma`, `sigma`, `dij`, `cp`, `cl`, `cd`, `cm`, `re`, `retheta`, `mach`, `alpha`, `ncrit`,
   `tau`, `rho`, `nu`, `delta`. A symbol is never spelled out into English. Drela's papers use
   q for surface speed and q∞ for the freestream, so `q_inviscid`, `q_viscous` and `qinf` are
   readable to anyone with the paper open.

2. **One symbol is one glued token.** Subscripts and operators that belong to the symbol in
   the equations are written without underscores: `retheta` (Reθ), `hstar` (H*), `hstarstar`
   (H**), `thetastar` (θ*), `sqrtctau` (Cτ^½), `sqrtctaueq` (Cτ_eq^½), `machsqd` (M∞²). A
   square says `sqd`, never `sq`, so it cannot be read as a square root.

3. **Underscore has exactly five meanings.**
   - a kind prefix: `n_` the length of a collection, `i_` an index;
   - a coordinate prefix: `x_`, `y_`, `s_`, `xi_` = the value of that coordinate at a named
     point: `x_le`, `y_te`, `s_le`, `s_stagnation`, `x_trip`, `x_transition`, `xi_transition`;
   - the derivative marker `_d_<token>` (see *Derivative tokens*);
   - a frame suffix naming the index system (see *Index systems*): `_station1`, `_station2`,
     `_node`, `_row`;
   - an English qualifier word: `cd_friction`, `machsqd_edge`, `ue_inviscid`,
     `ue_compressible`, `rho_stagnation`, `wake_built`.

   Counters are plain words, not `n_`: `iterations`, `max_iterations`.

4. **Opaque abbreviations are expanded** to the symbol or the word: `thet` → `theta`,
   `dstr` → `dstar`, `tstr` → `thetastar`, `uedg` → `ue`, `xssi` → `xi`, `vsrez` → `residual`,
   `lblini` → `bl_initialised`, `SCCON` → `LAG_CONSTANT`.

5. **Every count and index names its index system.** `n_foil_nodes`, `i_te_station`,
   `i_stagnation_node`. An index array that is stored per element of another frame is named
   for the frame it points *into*; the frame it is indexed *by* is given by the container:
   `state.i_node[side][i_station]` is the panel node of a BL station.

6. **Output coefficients keep their names**, so what refers to them is named from them:
   `cm_ref_x`, `cl_d_alpha`, `cd_friction`, `ldratio`.

7. **British English** in identifiers, keys, commands and docs: `analyse`, `panelled`,
   `initialised`, `normalised`. "Aerofoil" is the word in prose; in identifiers and keys it is
   `foil` (`PanelledFoil`, `FoilNodes`, `n_foil_nodes`, `panel_foil`, JSON `foil`), matching
   the `plot foil` command.

8. **The panel-node frame has two named subsets**, foil nodes 1..N and wake nodes N+1..N+NW,
   and every count, type and key says which: `n_foil_nodes` / `n_wake_nodes`, `FoilNodes` /
   `WakeNodes`, `geometry` / `geometry.wake`. A bare `_node` index (`i_le_node`,
   `i_stagnation_node`, `i_node`) is an index into the whole frame.

9. **Coordinates are chord-normalised everywhere** (`Geometry`, `PanelledFoil`, every output
   block), so they are plainly `x`, `y`, `s`; the normalisation is documented on the type, not
   encoded in the name (`x_c` would read as "x at point c" under rule 3).

10. **Type names** are descriptive and unique across the crate; the first doc line names the
    XFOIL COMMON block or subroutine they mirror.

11. **Functions translating an XFOIL subroutine** get a descriptive name plus
    `#[doc(alias = "SETBL")]` so rustdoc search by the Fortran name finds them, and the doc
    comment states the equivalent XFOIL subroutine, or subroutines when functionality was
    merged or split (`panel_foil` ≡ SCALC + SEGSPL + LEFIND + TECALC + NCALC + APCALC;
    `Session::alpha` ≡ OPER `ALFA` = SPECAL + VISCAL). Module files keep the routine name
    (`setbl.rs`) because the fixture pipeline, `coverage.toml`, the tests and the docs address
    subroutines by that name.

12. **Inside a translated body**, short locals bound at the top of the function from the named
    fields are allowed (`let hk1 = station1.hk;`) so formula lines stay diffable against the
    Fortran. That binding block is the function's legend.

13. **JSON keys are the variable names.** A key that is not a one-to-one print of a variable
    (a pair, a unit conversion, a derived ratio) is listed in the mapping document's
    "JSON ↔ variable" table. CLI flags are user vocabulary, with no aliases.

## Index systems

XFOIL runs four discretisations at once. Every count, index, coordinate and station label is
assigned to exactly one of them, and the suffix says which.

| System | What is indexed | Count | Suffix | Coordinate | XFOIL |
|---|---|---|---|---|---|
| **buffer geometry** | input coordinates before panelling | `n_points` | `_point` | — | XB, YB, NB |
| **panel nodes** | foil nodes 1..N then wake nodes N+1..N+NW | `n_foil_nodes`, `n_wake_nodes` | `_node` | `s` (foil arc coordinate) | X, Y, S, N, NW, I |
| **BL stations** | per side: stagnation → TE (→ wake on side 2) | `n_stations[side]` | `_station` | `xi` (BL arc coordinate from the stagnation point) | XSSI, NBL, IBL, IS |
| **Newton rows** | one row per BL station on both sides, in system order | `n_rows` | `_row` | — | NSYS, IV |

A named location carries its system's suffix: `i_stagnation_node` (IST), `i_le_node`,
`i_te_station` (IBLTE), `i_transition_station` (ITRAN), `i_te_row_upper` (IVTE1),
`i_wake_row` (IVZ). Per-station arrays that point into another frame are `i_node` (IPAN) and
`i_row` (ISYS), indexed `[side][i_station]`. XFOIL's "1"/"2" station pair of an interval is
the BL-station frame, so it is written `_station1` / `_station2`: `sqrtctau_d_ue_station1` is
∂(Cτ^½)/∂Ue at the upstream station of the interval. Loop variables are `i_node`, `i_station`,
`i_row` and `side`.

## Derivative tokens

A sensitivity is `<name>_d_<token>`, read "derivative of *name* with respect to *token*". The
token is the name of the variable differentiated against, so it needs no separate legend, but
XFOIL's suffixes map onto it as follows.

| Token | With respect to | XFOIL suffix |
|---|---|---|
| `ue` | compressible edge velocity Ue | `_U`, `_U1`, `_U2` |
| `uei` | incompressible edge velocity | `_UEI` |
| `theta` | momentum thickness θ | `_T` |
| `dstar` | displacement thickness δ* | `_D` |
| `sqrtctau` | Cτ^½, the lag variable | `_S` |
| `ampl` | amplification factor N (`n` would read as a count) | `_A`, `_A1`, `_A2` |
| `xi` | BL arc coordinate ξ | `_X`, `_X1`, `_X2` |
| `machsqd` | freestream Mach number squared M∞² | `_MS`, `_MSQ` |
| `re` | Reynolds number | `_RE` |
| `h` | shape factor H | (HKIN) |
| `hk` | kinematic shape factor Hk | `_HK` |
| `retheta` | momentum-thickness Reynolds number Rθ | `_RT` |
| `alpha` | angle of attack | `_A`, `_ALF`, `_ALFA` |
| `cl` | lift coefficient | `_CL`, `_CLMR` |
| `gamma` | vortex sheet strength γ | `_G`, `_GO`, `_GP` |
| `sigma` | source strength σ | `_M` (mass) |
| `mass` | mass defect m | `_M` |
| `x_trip` | forced-transition location | `_XF` |
| `free` | the free operating variable (see below) | `_AC` |
| `qinf` | freestream speed | `_QINF` |

A station-frame suffix follows the token when XFOIL's suffix carries one:
`XT_T1` → `xi_transition_d_theta_station1`, `CFM_U2` → `cf_d_ue_station2`.

## Interpretation notes

These are the XFOIL names whose meaning had to be established from the equations before they
could be named. The evidence is recorded so the choice is not re-litigated.

**`AC` → `free_variable`.** `AC` appears only in UPDATE (`DAC`, `U_AC`, `Q_AC`, `CL_AC`) and
as the second right-hand-side column of `VDEL` that SETBL fills and BLSOLV solves. In all
three it is the same thing: the operating-point variable that is *not* prescribed. With α
fixed (`LALFA`) it is CL, because Re and Mach depend on CL through MRCL and CL must move with
the BL solution; with CL fixed it is α. One concept, one name, documented as "α or CL,
whichever is not prescribed"; token `_d_free`.

**`CTAU`, `S1`/`S2` → `sqrtctau`.** Cτ is the maximum shear-stress coefficient τ_max/(ρUe²).
The digit in `S2` is the station index (COM2 is `X2 U2 T2 D2 S2 AMPL2`, all "at station 2").
`S` itself is Cτ^½, shown three ways in the Fortran: the equilibrium value is
`CQ2 = SQRT(CTCON·HS2·(HK2−1)·HKC²/((1−US2)·H2·HK2²))`; the shear-lag residual is
`SCC·(CQA − SA·ALD)·Δξ − 2·δ·ln(S2/S1) + …`, which is Drela's
(δ/Cτ)·dCτ/dξ = 5.6·(Cτ_eq^½ − Cτ^½) + … written for Cτ^½ (since (δ/Cτ)·dCτ/dξ =
2δ·d ln Cτ^½/dξ); and the dissipation `DI = (½·CF·US + ST²·(1−US))·2/HS` uses ST² = Cτ.
`XFOIL.INC` itself says "CTAU sqrt(max shear coefficient)". The equilibrium value `CQ`/`CTQ`
is `sqrtctaueq`. In laminar stations the same array holds the amplification N; that overload
is XFOIL's data model and is documented on the field, not hidden in the name.

**`US` → `us`; `USLP` → `us_plot_scale`.** Us is Drela's equivalent normalised wall-slip
velocity Us/Ue, a closure variable (`US2 = 0.5·HS2·(1 − (HK2−1)/(GBCON·H2))`, capped) used in
the dissipation and equilibrium-shear closures. `USLP = 1.6/(1+Us)` is computed in SETBL only
for XFOIL's velocity-profile plot; nothing in the solver reads it, and the outputs print `us`.

**`RMSBL` → `residual`.** XFOIL's convergence metric is the rms of the normalised Newton
changes. It is called `residual` throughout (state, iteration record, JSON) because that is
what a reader expects a convergence metric to be called; `RMXBL` is `residual_max`. The true
equation residuals of one interval (`VSREZ`) live only inside `IntervalSystem.residual`.

**`GAMMA` → `gamma_gas`; `GAM` → `gamma`.** Both are γ in Drela's papers: the ratio of specific
heats and the vortex-sheet strength. The sheet strength is referenced throughout the panel
method and in every `_d_gamma` sensitivity, so it keeps the bare `gamma`; the gas constant
(1.4 for air, used by the compressibility relations in COMSET, BLKIN and UPDATE) is qualified.

**`FlowConditions` versus `PointResult`.** The inputs shared by every point of a polar (Re,
Mach, Ncrit, wake length, trip, iteration limit, CL-dependence modes) are the *flow
conditions*. An operating point in the aerodynamic sense is those plus α or CL. What the
solver returns for one point (forces, transition, convergence) is a *result*. The JSON carries
`conditions` and `results` in both the analysis and the polar files, and the inviscid analysis
file is a strict subset of the viscous one.

**The BLPAR closure constants.** `SCCON = 5.6` is the constant in the shear-lag equation
(`LAG_CONSTANT`); `GACON`, `GBCON`, `GCCON` are the G–β locus G = A·√(1+B·β) + C/(H·Rθ·√(Cf/2))
(`GBETA_LOCUS_A`, `GBETA_LOCUS_B`, `GBETA_LOCUS_WALL`); `DLCON` is the wake/wall
dissipation-length ratio (`WAKE_DISSIPATION_LENGTH_RATIO`); `CTRCON`, `CTRCEX` set the Cτ^½ at
transition as 1.8·exp(−3.3/(Hk−1)) times the equilibrium value
(`TRANSITION_SQRTCTAU_FACTOR`, `TRANSITION_SQRTCTAU_EXPONENT`); `DUXCON` weights the
pressure-gradient term of the lag equation (`LAG_PRESSURE_GRADIENT_WEIGHT`); `CTCON` is the
coefficient in the equilibrium Cτ^½ closure (`SQRTCTAUEQ_COEFFICIENT`); `CFFAC` scales the
turbulent Cf correlation (`CF_TURBULENT_FACTOR`).

## Related notes

- [`docs/xfoil-reference/xfoil-to-yfoil-mapping.md`](../xfoil-reference/xfoil-to-yfoil-mapping.md)
  — the remap table: every XFOIL name, the yFoil name, and what it is
- [[git-commits-and-versioning]] — commit message form for the rename commits (`REF:`)
