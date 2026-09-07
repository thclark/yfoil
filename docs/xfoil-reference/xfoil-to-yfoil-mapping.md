# XFOIL → YFoil mapping

Every XFOIL name (COMMON variable, subroutine, local of note) and the YFoil name it became, with
what the quantity is. The *rules* behind these names — symbols as names, the underscore
meanings, the four index systems, the derivative tokens, and the interpretation notes on `AC`,
`CTAU`, `US`, `RMSBL` and `GAMMA` — are in
[`docs/conventions/naming.md`](../conventions/naming.md); this file is their application and is
kept current whenever a field, variable or function is added or renamed (CLAUDE.md, Architecture).

Indexing is unchanged from XFOIL: 1-based, side 1 = upper and 2 = lower, the wake appended to
side 2 (`NBL(2) = IBLTE(2) + NW`), station arrays `[side][i_station]`. "—" in the XFOIL column
means no Fortran counterpart. XFOIL's own dumps under `tests/fixtures/xfoil/` keep the Fortran
names by design; everything YFoil writes uses the YFoil names.

## Table 1: `BlState` → `SolverState` (`src/solver/blstate.rs`), panel-node section

| XFOIL | YFoil | What it is |
|---|---|---|
| N | `n_foil_nodes` | number of aerofoil panel nodes |
| NW | `n_wake_nodes` | number of wake nodes |
| X, Y | `x`, `y` | node coordinates, aerofoil then wake, chord-normalised |
| S | `s` | arc coordinate along the aerofoil (spline parameter) |
| XP, YP | `dxds`, `dyds` | spline derivatives dx/ds, dy/ds at aerofoil nodes |
| NX, NY | `normal_x`, `normal_y` | outward unit normal components (`nx` would read as a count) |
| APANEL | `panel_angle` | panel angle, counter-clockwise positive |
| SIG | `sigma` | mass-defect source strength per node |
| QINF | `qinf` | freestream speed q∞ (1.0) |
| ALFA | `alpha` | angle of attack, radians |
| GAM | `gamma` | surface vortex sheet strength (= tangential velocity) |
| GAM_A | `gamma_d_alpha` | |
| QINVU(.,1..2) | `q_inviscid_basis` | inviscid surface speed at α = 0° and 90° |
| QINV | `q_inviscid` | inviscid surface speed at current α |
| QINV_A | `q_inviscid_d_alpha` | |
| QVIS | `q_viscous` | surface speed including source influence |
| CHORD | `chord` | |
| SLE | `s_le` | value of s at the leading edge |
| XLE, YLE | `x_le`, `y_le` | leading-edge point |
| XTE, YTE | `x_te`, `y_te` | trailing-edge midpoint |
| ANTE | `te_thickness_normal` | TE thickness projected perpendicular to the TE bisector |
| ASTE | `te_thickness_parallel` | TE thickness projected along the bisector |
| DSTE | `te_gap` | trailing-edge gap length |
| SHARP | `sharp_te` | TE gap below 1e-4·chord |
| IST | `i_stagnation_node` | stagnation point lies on the panel between nodes IST and IST+1 |
| SST | `s_stagnation` | value of s at the stagnation point |
| SST_GO | `s_stagnation_d_gamma_node0` | dSST/dγ at node IST |
| SST_GP | `s_stagnation_d_gamma_node1` | dSST/dγ at node IST+1 |
| DIJ | `dij` | dQtan(i)/dσ(j), dense (N+NW)² |
| WGAP | `wake_gap` | dead-air thickness inside the wake behind a blunt TE |
| XCMREF, YCMREF | `cm_ref_x`, `cm_ref_y` | moment reference point |
| CPI, CPV | `cp_inviscid`, `cp_viscous` | pressure coefficient per node |

## Table 2: `SolverState`, BL-station section (`[side][i_station]`)

| XFOIL | YFoil | What it is |
|---|---|---|
| NBL(IS) | `n_stations` | last station index on each side (side 2 includes the wake) |
| IBLTE(IS) | `i_te_station` | station at the trailing edge |
| ITRAN(IS) | `i_transition_station` | station of the transition interval |
| NSYS | `n_rows` | rows in the BL Newton system |
| IPAN(IBL,IS) | `i_node` | panel node of each station |
| VTI(IBL,IS) | `velocity_sign` | ±1 between panel tangential velocity and BL edge velocity |
| ISYS(IBL,IS) | `i_row` | Newton row of each station |
| XSSI | `xi` | BL arc coordinate ξ from the stagnation point |
| UEDG | `ue` | edge velocity |
| UINV | `ue_inviscid` | edge velocity without source influence |
| UINV_A | `ue_inviscid_d_alpha` | |
| MASS | `mass_defect` | m = Ue·δ* |
| THET | `theta` | momentum thickness |
| DSTR | `dstar` | displacement thickness |
| CTAU | `sqrtctau` | Cτ^½ at turbulent/wake stations; amplification N at laminar stations (documented overload) |
| DELT | `delta` | boundary-layer thickness (plotting) |
| TSTR | `thetastar` | kinetic-energy thickness θ* = H*·θ |
| USLP | `us_plot_scale` | 1.6/(1+Us), XFOIL's profile-plot scale (never read by the solver) |
| GUXQ, GUXD | *delete* | never assigned on the analysis path; commented out even in blplot.f |
| TAU | `tau` | wall shear stress ½ρUe²Cf (plotting) |
| DIS | `dissipation` | ½ρUe³·CD·H* (plotting) |
| CTQ | `sqrtctaueq` | equilibrium Cτ^½ |
| XSTRIP(IS) | `x_trip` | forced-transition x/c per side; ≥ 1 means free |
| XSSITR(IS) | `xi_transition` | ξ of actual transition |
| TFORCE(IS) | `transition_forced` | |
| XOCTR, YOCTR | `x_transition`, `y_transition` | actual transition point, chord fractions |
| TINDEX(IS) | `transition_node_fraction` | fractional panel-node position of transition (plotting) |
| COM1 | `station1` | the upstream station of the current interval, persisted between calls |
| COM2 | `station2` | the current station |
| XT block | `transition` | transition location and sensitivities (table 6) |

## Table 3: `SolverState`, flow conditions, flags, forces

| XFOIL | YFoil | What it is |
|---|---|---|
| MINF1 | `mach_cl1` | freestream Mach at CL = 1 (the user's input) |
| REINF1 | `re_cl1` | Reynolds number at CL = 1 (the user's input) |
| MINF | `mach` | Mach at the current CL |
| REINF | `re` | Reynolds at the current CL |
| MATYP | `mach_cl_dependence` | enum `Fixed`, `InverseSqrtCl` (was 1/2) |
| RETYP | `re_cl_dependence` | enum `Fixed`, `InverseSqrtCl`, `InverseCl` (was 1/2/3) |
| IDAMP | `amplification_model` | enum `Envelope` (DAMPL), `ModifiedEnvelope` (DAMPL2) |
| MINF_CL, REINF_CL | `mach_d_cl`, `re_d_cl` | |
| LALFA | `alpha_specified` | true: α fixed, CL solved; false: CL fixed, α solved |
| CLSPEC | `cl_specified` | target CL when `!alpha_specified` |
| ACRIT(IS) | `ncrit` | log critical amplification ratio per side |
| VACCEL | `elimination_threshold` | BLSOLV skips off-diagonal entries below this |
| GAMMA | `gamma_gas` | ratio of specific heats Cp/Cv (1.4 for air); see the GAMMA note |
| TKLAM | `karman_tsien` | λ = M²/(1+√(1−M²))² |
| TKL_MSQ | `karman_tsien_d_machsqd` | |
| LWAKE | `wake_built` | wake geometry exists |
| LIPAN | `pointers_built` | station→node and station→row maps exist |
| LBLINI | `bl_initialised` | BL arrays have been marched once |
| LWDIJ | `dij_wake_built` | wake columns of DIJ exist |
| LVISC | `viscous` | viscous analysis requested |
| LVCONV | `converged` | a converged BL solution exists |
| AWAKE | `alpha_wake` | α the wake geometry was built for |
| AVISC, MVISC | `alpha_converged`, `mach_converged` | α and Mach of the converged BL solution |
| CL, CM, CD | same | |
| CDF, CDP | `cd_friction`, `cd_pressure` | |
| CL_ALF, CL_MSQ | `cl_d_alpha`, `cl_d_machsqd` | |

## Table 4: `BLStationState` → `StationState` (`/V_VAR2/`, XBL.INC)

The struct is instantiated as `station1` and `station2`, so the frame suffix lives on the
instance (`station2.sqrtctau`), not on the fields. Sensitivities follow `<name>_d_<token>`;
base names are listed with the sensitivities XFOIL carries.

| XFOIL | YFoil | What it is |
|---|---|---|
| X2 | `xi` | BL arc coordinate at the station |
| U2, U2_UEI, U2_MS | `ue`, `ue_d_uei`, `ue_d_machsqd` | Kármán–Tsien-corrected edge velocity |
| T2 | `theta` | |
| D2 | `dstar` | δ* excluding the wake gap |
| S2 | `sqrtctau` | Cτ^½, the turbulent lag variable |
| AMPL2 | `ampl` | amplification factor N (same token as the derivative suffix) |
| DW2 | `wake_gap` | wake-gap part of δ* |
| H2, H2_T2, H2_D2 | `h`, `h_d_theta`, `h_d_dstar` | H = δ*/θ |
| M2, M2_U2, M2_MS | `machsqd_edge`, `_d_ue`, `_d_machsqd` | edge Mach number squared (the token is freestream M∞²) |
| R2 (+_U2 _MS) | `rho`, … | edge density / stagnation density |
| V2 (+_U2 _MS _RE) | `nu`, … | edge kinematic viscosity ÷ (U∞·c), i.e. 1/Re at edge conditions |
| HK2 (+_U2 _T2 _D2 _MS) | `hk`, … | kinematic shape factor |
| RT2 (+_U2 _T2 _MS _RE) | `retheta`, … | Rθ |
| HC2 (+_U2 _T2 _D2 _MS) | `hstarstar`, … | density-thickness shape parameter H** |
| HS2 (+_U2 _T2 _D2 _MS _RE) | `hstar`, … | energy shape factor H* |
| US2 (+…) | `us`, … | equivalent normalised wall-slip velocity Us/Ue |
| CQ2 (+…) | `sqrtctaueq`, … | equilibrium Cτ^½ |
| CF2 (+…) | `cf`, … | |
| DI2 (+_U2 _T2 _D2 _S2 _MS _RE) | `cdiss`, …, `cdiss_d_sqrtctau` | dissipation coefficient in XFOIL's form 2·CD/H* |
| DE2 (+_U2 _T2 _D2 _MS) | `delta`, … | δ from Green's correlation |
| — | `mass_defect` | Ue·δ* (YFoil-only cached copy) |

## Table 5: `BLGlobalParams` → `FlowParameters` (`/V_VAR/`) and the closure constants (`/BLPAR/`, set in BLPINI)

| XFOIL | YFoil | What it is |
|---|---|---|
| IDAMPV | `amplification_model` | copy of IDAMP for the BL routines |
| QINFBL | `qinf` | |
| TKBL, TKBL_MS | `karman_tsien`, `karman_tsien_d_machsqd` | |
| RSTBL, RSTBL_MS | `rho_stagnation`, `rho_stagnation_d_machsqd` | ρ0/ρ∞ |
| HSTINV, HSTINV_MS | `h_stagnation_inv`, `h_stagnation_inv_d_machsqd` | 1/h0 in freestream units |
| REYBL, REYBL_MS, REYBL_RE | `re`, `re_d_machsqd`, `re_d_re` | Reynolds number on freestream density and viscosity |
| GAMBL, GM1BL | `gamma_gas`, `gamma_gas_m1` | γ, γ−1 |
| HVRAT | `sutherland_ratio` | Sutherland constant / freestream temperature (0 on the analysis path) |
| SCCON = 5.6 | `LAG_CONSTANT` | the 5.6 in the shear-lag equation (δ/Cτ)·dCτ/dξ = 5.6·(Cτ_eq^½ − Cτ^½) + … |
| GACON = 6.7 | `GBETA_LOCUS_A` | G–β locus: G = A·√(1 + B·β) + C/(H·Rθ·√(Cf/2)) |
| GBCON = 0.75 | `GBETA_LOCUS_B` | |
| GCCON = 18.0 | `GBETA_LOCUS_WALL` | the wall term C of the locus |
| DLCON = 0.9 | `WAKE_DISSIPATION_LENGTH_RATIO` | wake/wall dissipation-length ratio Lo/L (applied in the wake) |
| CTRCON = 1.8, CTRCEX = 3.3 | `TRANSITION_SQRTCTAU_FACTOR`, `TRANSITION_SQRTCTAU_EXPONENT` | Cτ^½ at transition = 1.8·exp(−3.3/(Hk−1)) × equilibrium value (TRDIF) |
| DUXCON = 1.0 | `LAG_PRESSURE_GRADIENT_WEIGHT` | weight on the (Ue-gradient − equilibrium gradient) term of the lag equation |
| CTCON = 0.5/(A²·B) | `SQRTCTAUEQ_COEFFICIENT` | coefficient in the equilibrium Cτ^½ closure (BLVAR) |
| CFFAC = 1.0 | `CF_TURBULENT_FACTOR` | multiplier on the turbulent Cf correlation (CFT) |

## Table 6: interval-level types (`/V_SYS/`, `/V_VARA/`, `/V_INT/`)

| XFOIL | YFoil | What it is |
|---|---|---|
| /V_SYS/ | `IntervalSystem` | the 4×5 linearised system for one interval |
| VS1, VS2 | `jacobian_station1`, `jacobian_station2` | ∂residual/∂(Cτ^½, θ, δ*, Ue, ξ) at each station |
| VSREZ | `residual` | the interval's equation residuals |
| VSR, VSM, VSX | `residual_d_re`, `residual_d_machsqd`, `residual_d_xi` | |
| /V_INT/ | `IntervalFlags` | |
| SIMI, TRAN, TURB, WAKE | `similarity`, `transition`, `turbulent`, `wake` | |
| ITYP | `FlowRegime` {`Laminar`, `Turbulent`, `Wake`} | closure-set selector |
| XT block | `Transition` | |
| XT | `xi_transition` | ξ of transition (the same quantity `SolverState.xi_transition[side]` stores per side after the sweep) |
| XT_A1, XT_X1, XT_T1, XT_D1, XT_U1 | `xi_transition_d_ampl_station1`, `xi_transition_d_xi_station1`, `xi_transition_d_theta_station1`, `xi_transition_d_dstar_station1`, `xi_transition_d_ue_station1` | |
| XT_X2, XT_T2, XT_D2, XT_U2 | `xi_transition_d_xi_station2`, … | |
| XT_MS, XT_RE, XT_XF | `xi_transition_d_machsqd`, `xi_transition_d_re`, `xi_transition_d_x_trip` | |
| TRCHEK2 outcome | `TransitionCheck` {`None`, `Free`, `Forced`}; payload `location` → `transition` | |
| CFM block | `MidpointCf` | Cf at the interval midpoint |
| CFM, CFM_MS, CFM_RE, CFM_U1 … | `cf`, `cf_d_machsqd`, `cf_d_re`, `cf_d_ue_station1`, … | |
| UPW block | `Upwinding` | |
| UPW, UPW_U1 … UPW_MS | `weight`, `weight_d_ue_station1`, …, `weight_d_machsqd` | |
| AX, AX_HK, AX_TH, AX_RT | `AmplificationRate` {`rate`, `rate_d_hk`, `rate_d_theta`, `rate_d_retheta`} | dN/dξ from DAMPL |
| AXSET outputs | `IntervalAmplificationRate` {`rate`, `rate_d_hk_station1`, `rate_d_theta_station1`, `rate_d_retheta_station1`, `rate_d_ampl_station1`, …} | averaged dN/dξ over an interval |

## Table 7: global Newton system and iteration records

| XFOIL | YFoil | What it is |
|---|---|---|
| VA/VB/VDEL/VM/VZ | `NewtonSystem` | the block system BLSOLV consumes by value |
| NSYS | `n_rows` | |
| VA | `diagonal` | 3×2 diagonal blocks per row |
| VB | `subdiagonal` | 3×2 blocks coupling to the previous row |
| VZ | `te_block` | block coupling the first wake row to the upper-surface TE row |
| VM | `mass_influence` | dense 3-vectors ∂equations/∂mass defect of every row |
| VDEL | `rhs` | column 0 residual, column 1 ∂residual/∂free variable; solution after the solve |
| IVTE1 | `i_te_row_upper` | row of the upper-surface TE station |
| IVZ | `i_wake_row` | first wake row, where the TE block applies |
| VACCEL | `elimination_threshold` | |
| S(N)−S(1) | `s_total` | scales the threshold |
| — | `NewtonDeltas` {`deltas`} | |
| — | keep | instrumentation |
| SETBL outputs | `AssembledSystem` | |
| — | `newton`, `flow` | |
| RE_CLMR, MSQ_CLMR | `re_d_cl`, `machsqd_d_cl` | from MRCL |
| M_CLS | `mach_d_cl` | |
| DULE1, DULE2 | `ue_le_mismatch` | UEDG − USAV at the first station per side |
| UPDATE outputs | `UpdateSummary` | |
| RLX | `relaxation` | under-relaxation factor applied |
| RMSBL, RMXBL | `residual`, `residual_max` | rms and largest normalised Newton change (see the RMSBL note) |
| VMXBL | `residual_max_variable` | enum {`Ampl`, `Sqrtctau`, `Theta`, `Dstar`, `Ue`} (was `char`) |
| IMXBL, ISMXBL | `i_residual_max_station`, `residual_max_side` | |
| DAC | `free_variable_change` | Newton change in the free variable before relaxation |
| CLNEW, CL_A, CL_MS, CL_AC | `cl_new`, `cl_d_alpha`, `cl_d_machsqd`, `cl_d_free` | |
| U_AC, Q_AC (locals) | `ue_d_free`, `q_d_free` | |
| one VISCAL iteration | `IterationRecord` | fields as `UpdateSummary` plus `alpha`, `mach`, `re`, forces, `i_stagnation_node`, `s_stagnation`, `i_transition_station`, `x_transition`, `converged` |
| EPS1 | `CONVERGENCE_TOLERANCE` | RMSBL < 1e-4 |

## Table 8: inviscid system

| XFOIL | YFoil | What it is |
|---|---|---|
| AIJ, BIJ, LADIJ | {`aij_lu`, `bij`, `dij_foil_built`} | dψ/dγ (LU-factored), dγ/dσ, flag |
| LUDCMP output | {`n`, `lu`, `pivots`} | |
| PSILIN outputs | `PanelInfluence` | streamfunction and velocity influence at one point |
| PSI, PSI_NI | `psi`, `psi_d_n` | ψ and ∂ψ/∂n |
| QTAN1, QTAN2 | `qtan_alpha0`, `qtan_alpha90` | tangential velocity at α = 0°, 90° |
| QTANM | `qtan_sigma` | tangential velocity induced by the sources |
| Z_QINF, Z_ALFA | `psi_d_qinf`, `psi_d_alpha` | |
| DZDG, DQDG | `psi_d_gamma`, `qtan_d_gamma` | per panel |
| DZDM, DQDM | `psi_d_sigma`, `qtan_d_sigma` | |
| PSWLIN outputs | `WakeSourceInfluence` | same names |

## Table 9: session, flow conditions and results (`src/solver/analysis.rs`)

| XFOIL | YFoil | What it is |
|---|---|---|
| OPER settings | `FlowConditions` | inputs shared by every point of a polar; also the serialised `conditions` block |
| REINF1, MINF1, ACRIT | `re` (`None` for inviscid), `mach`, `ncrit` | |
| ITMAX | `max_iterations` | |
| WAKLEN | `wake_length` | chords |
| VACCEL | `elimination_threshold` | |
| XSTRIP | `x_trip` | |
| MATYP, RETYP, IDAMP | `mach_cl_dependence`, `re_cl_dependence`, `amplification_model` | enums |
| — | `Session` {`state`, `inviscid`, `conditions`}, private with accessors | |
| — | `PointResult` | the solved result of one point; also the serialised `results` block |
| ALFA | `alpha` | radians |
| CL, CM, CD, CDF, CDP, CL_ALF | `cl`, `cm`, `cd`, `cd_friction`, `cd_pressure`, `cl_d_alpha` | viscous-only fields are `Option` |
| XOCTR(1..2), YOCTR(1..2) | `transition_upper: [f64; 2]`, `transition_lower: [f64; 2]` | (x, y) of the transition point per side |
| ITRAN | `i_transition_station` | |
| RMSBL | `residual` | last iteration |
| — | `iterations`, `iteration_records` | |
| ASEQ | keep | degrees |
| NSEQEX | `max_consecutive_failures` | |
| — | {`results`, `failed_alphas`, `conditions`, `completed`} | |

## Table 10: functions (each carries `#[doc(alias = "XFOIL NAME")]`)

| XFOIL | YFoil | What it does |
|---|---|---|
| BLPRV | `set_primary_variables` | loads ξ, N/Cτ^½, θ, δ*, wake gap, Ue and applies Kármán–Tsien |
| BLKIN | `set_kinematic_variables` | H, Mₑ², ρ, ν, Hk, Rθ and sensitivities |
| BLVAR | `set_closure_variables` | H**, H*, Us, Cτ_eq^½, Cf, CD, δ for the regime |
| BLMID | `MidpointCf::compute` | |
| BLDIF | `assemble_interval_equations` | |
| BLDIF blocks | `shear_lag_equation`, `momentum_equation`, `shape_equation`, `upwinding` | doc: "part of BLDIF" |
| TRDIF | `assemble_transition_equations` | |
| BLSYS | `assemble_interval_system` | |
| TESYS | `assemble_te_system` | |
| TRCHEK2 | `check_transition` | doc: TRCHEK2; XFOIL's TRCHEK wrapper is not translated |
| DAMPL, DAMPL2 | `amplification_rate`, `amplification_rate_modified` | |
| DAMPL | *delete* | duplicate with a wrong doc header |
| AXSET | `interval_amplification_rate` | |
| DSLIM | `limit_dstar` | keeps Hk ≥ HKLIM |
| HKIN | `hk_from_h` | returns (Hk, dHk/dH, dHk/dM²) |
| CFL, HSL, DIL | `cf_laminar`, `hstar_laminar`, `cdiss_laminar` | |
| CFT, HST, HCT | `cf_turbulent`, `hstar_turbulent`, `hstarstar` | |
| DILW | `cdiss_wake` | |
| DIT, —, — | *delete* | unused |
| closure return | `Closure` {`value`, `value_d_hk`, `value_d_retheta`, `value_d_machsqd`} | |
| MRCHUE | `march_direct` | prescribed Ue, inverse step where separating |
| MRCHDU | `march_prescribed_dstar` | current Ue and δ*, to locate transition |
| SETBL | `assemble_newton_system` | |
| MRCL | `set_mach_re_from_cl` | |
| BLSOLV | `solve_newton_system` | |
| GAUSS | keep | |
| UPDATE | `apply_newton_update` | |
| VISCAL | `solve_viscous` | |
| SPECAL, SPECCL | `solve_inviscid_at_alpha`, `solve_inviscid_at_cl` | |
| OPER ALFA / CL / ASEQ | `Session::alpha`, `Session::cl`, `Session::sequence_point` | doc: OPER `ALFA` = SPECAL + VISCAL, etc. |
| COMSET | `set_compressibility` | |
| CPCALC, CLCALC, CDCALC | `compute_cp`, `compute_cl_cm`, `compute_cd` | |
| QISET | `set_q_inviscid` | |
| UICALC | `set_ue_inviscid` | |
| UECALC | `set_ue_from_q_viscous` | |
| QVFUE | `set_q_viscous_from_ue` | |
| GAMQV | `set_gamma_from_q_viscous` | |
| UESET | `set_ue_with_sources` | |
| DSSET | `set_dstar_from_mass` | |
| TECALC | `set_te_thickness` | |
| STFIND, STMOVE | `find_stagnation`, `move_stagnation` | |
| IBLPAN, XICALC, IBLSYS | `map_stations_to_nodes`, `set_station_xi`, `map_stations_to_rows` | |
| XIFSET | `xi_trip` | ξ of the trip on a side |
| SINVRT | `s_at_x` | inverts the spline for x |
| XYWAKE, QWCALC | `build_wake`, `set_wake_q_basis` | |
| SETEXP | `exponential_spacing` | |
| QDCALC, PSWLIN | `build_dij`, `wake_source_influence` | |
| PSILIN | `panel_influence` | |
| GGCALC | `build_inviscid_system` | |
| LUDCMP, BAKSUB | `lu_decompose`, `lu_back_substitute` | |
| ATANC | `continuous_atan2` | |
| SPLINE, SEVAL, DEVAL, D2VAL | `spline_derivatives`, `spline_value`, `spline_slope`, `spline_second_derivative` | |
| SEGSPL, CURV, LEFIND, SCALC | `spline_segmented`, `curvature`, `find_le`, `arc_coordinate` | `find_leading_edge` and `calculate_arc_length` duplicates merge into these |
| TRISOL | `solve_tridiagonal` (one, XFOIL argument order) | |
| NCALC, APCALC | `node_normals`, `panel_angles` | |
| PANGEN | `repanel_by_curvature` | |
| — | keep | doc: no XFOIL equivalent |
| NACA4/NACA5 | `naca_4digit`, `naca_5digit` with a `Thickness` {`Perpendicular`, `Vertical`} argument | the un-suffixed name currently holds the non-XFOIL algorithm |
| SCALC + SEGSPL + LEFIND + TECALC + NCALC + APCALC | `panel_foil` | doc lists all six |
| OPER ALFA (fresh session) | `analyse` | |

## Table 11: geometry and output types, JSON keys

| XFOIL | YFoil | What it is |
|---|---|---|
| XB, YB | `Geometry` {`x`, `y`, `cm_ref`} | buffer-geometry points, chord-normalised; `cm_ref: [x, y]` |
| /CR05/ | `PanelledFoil` {`x`, `y`, `s`, `dxds`, `dyds`, `normal_x`, `normal_y`, `panel_angle`, `n_foil_nodes`, `s_le`, `i_le_node`, `chord`, `sharp_te`, `cm_ref`} | |
| — | one `AnalysisOutput` {`foil`, `conditions`, `results`, `geometry`, `surface`, `boundary_layer: Option`} | inviscid output is the same shape with `boundary_layer` absent and `geometry.wake` absent |
| — | `PolarPoint` | the `results` record of one polar point |
| — | fold into `FlowConditions` | one flow-condition type instead of three |
| QINV/QVIS, CPI/CPV | `SurfaceDistributions` {`q`, `cp`} | per panel node; the inviscid or viscous pair according to `conditions.re` |
| GEOPAR THICK | `y_extent` | it is max_y − min_y |
| — | keep names | |
| — | `FoilNodes`, `WakeNodes` | fields as table 1 |
| BLDUMP columns | `SideStations` {`i_station`, `i_node`, `x`, `y`, `xi`, `cp`, `primaries`, `closures`, `lagged_closures: Option`} | JSON below |
| UEDG THET DSTR CTAU MASS | `Primaries` {`ue`, `theta`, `dstar`, `sqrtctau`, `mass_defect`} | the converged solver state |
| BLPRV→BLKIN→BLVAR on the primaries | `Closures` {`ue_compressible`, `h`, `hk`, `hstar`, `cf`, `cdiss`, `delta`, `sqrtctaueq`, `us`, `retheta`, `machsqd_edge`} | closures evaluated on the converged primaries |
| TAU DIS CTQ DELT USLP TSTR | `LaggedClosures` {`tau`, `dissipation`, `sqrtctaueq`, `delta`, `us_plot_scale`, `thetastar`, `hstar_dump`, `cf_dump`}, only with `--include-lagged-closures` | XFOIL's arrays left by the last march, one iterate behind the primaries; what XFOIL's DUMP prints |
| VPLO variables | {`Dstar`, `Theta`, `Delta`, `H`, `Hk`, `Hstar`, `Ue`, `Cf`, `Cdiss`, `Sqrtctau`, `Sqrtctaueq`, `Us`, `MassDefect`, `Cp`}; CLI tokens unchanged | |
| IST, SST | {`i_stagnation_node`, `s_stagnation`, `x`, `y`} | |
| XOCTR, YOCTR | {`i_station`, `forced`, `x_transition`, `y_transition`, `s_transition`}; the near-duplicate interpolated `x`, `y` are dropped | resolves the `x_c` collision with `Geometry` |
| — | delete; `output::PanelStyle`, `ImageFormat` derive `ValueEnum` | |

## Table 12: CLI (no aliases; a short flag only where it is the same letter under every subcommand)

| XFOIL | YFoil | Note |
|---|---|---|
| OPER ALFA | `analyse` | |
| ITER | `--max-iterations` | |
| VISC | keep | |
| VPAR N | `--ncrit` (no short flag) | `-n` is `--panels` under `geometry` |
| PANE N | keep | |
| CTERAT | `--te-le-ratio` | the value *is* TE/LE; doc corrected |
| PANGEN vs none |cosine` | `--method curvature\|cosine` | |
| NACA4 thickness |xfoil` | `--thickness perpendicular\|vertical` | |
| DUMP | keep | |
| — | `--include-lagged-closures` | adds `lagged_closures` to the boundary-layer output |

## Table 13: tests and xtask

| XFOIL | YFoil | Note |
|---|---|---|
| ACRIT | `ncrit` | |
| — | `FixtureStation`, `FixtureSide`; `delta_star` → `dstar` | |
| THET_TE1 … | `theta_te_station1` … | the tracked JSON fixture inputs are translated to the new keys (step 9) |
| REINF1, ITMAX | `re`, `max_iterations` | `xtask/fixtures-config/cases.toml` translated in the same step |

---

## JSON ↔ variable: keys that are not a one-to-one print of a variable

Every other key in the analysis, polar and geometry-info JSON is the name of the variable it prints. These are the exceptions, each
with the variables it is formed from.

| Key | Formed from | Why it differs |
|---|---|---|
| `alpha_deg` | `alpha` (radians) | unit conversion, stated in the key |
| `transition_upper`, `transition_lower` | `x_transition[side]`, `y_transition[side]` | a point printed as a pair, as `cm_ref` already is |
| `cm_ref` | `cm_ref_x`, `cm_ref_y` | same |
| `ldratio` | `cl / cd` | derived, not stored |
| `residual` (in `results`) | `IterationRecord::residual` of the last iteration | the per-iteration records are not printed |
| `surface.q`, `surface.cp`, per-station `cp` | `q_inviscid` / `cp_inviscid` when `conditions.re` is null, `q_viscous` / `cp_viscous` otherwise | the state holds both; the output holds the one the analysis produced |
| `closures.us` | `us` (BLVAR), replacing `uslp` = 1.6/(1+Us) | prints the closure variable rather than XFOIL's plot scale of it |
| `i_node` (per station) | `i_node[side][i_station]` | same name; the frame it is indexed by is the enclosing `upper`/`lower`/`wake` block |
| `summary.cl_max`, `alpha_at_cl_max`, `ldratio_max`, `cl_at_ldratio_max`, `cd0`, `n_converged`, `n_failed` | polar post-processing | derived summary, no state variable |
| `y_extent`, `x_range`, `y_range`, `max_curvature`, `first_point`, `last_point` | geometry post-processing | derived summary |
| `lagged_closures.hstar_dump`, `cf_dump` | XFOIL DUMP's `H*` = TSTR/THET and `Cf` = TAU/(½q∞²) | reproductions of XFOIL's printed columns, diagnostics only |

---


## JSON ↔ variable: keys that are not a one-to-one print of a variable

Every other key in the analysis, polar and geometry-info JSON is the name of the variable it prints. These are the exceptions, each
with the variables it is formed from.

| Key | Formed from | Why it differs |
|---|---|---|
| `alpha_deg` | `alpha` (radians) | unit conversion, stated in the key |
| `transition_upper`, `transition_lower` | `x_transition[side]`, `y_transition[side]` | a point printed as a pair, as `cm_ref` already is |
| `cm_ref` | `cm_ref_x`, `cm_ref_y` | same |
| `ldratio` | `cl / cd` | derived, not stored |
| `residual` (in `results`) | `IterationRecord::residual` of the last iteration | the per-iteration records are not printed |
| `surface.q`, `surface.cp`, per-station `cp` | `q_inviscid` / `cp_inviscid` when `conditions.re` is null, `q_viscous` / `cp_viscous` otherwise | the state holds both; the output holds the one the analysis produced |
| `closures.us` | `us` (BLVAR), replacing `uslp` = 1.6/(1+Us) | prints the closure variable rather than XFOIL's plot scale of it |
| `i_node` (per station) | `i_node[side][i_station]` | same name; the frame it is indexed by is the enclosing `upper`/`lower`/`wake` block |
| `summary.cl_max`, `alpha_at_cl_max`, `ldratio_max`, `cl_at_ldratio_max`, `cd0`, `n_converged`, `n_failed` | polar post-processing | derived summary, no state variable |
| `y_extent`, `x_range`, `y_range`, `max_curvature`, `first_point`, `last_point` | geometry post-processing | derived summary |
| `lagged_closures.hstar_dump`, `cf_dump` | XFOIL DUMP's `H*` = TSTR/THET and `Cf` = TAU/(½q∞²) | reproductions of XFOIL's printed columns, diagnostics only |

