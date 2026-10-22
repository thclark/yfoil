# XFOIL's COMMON blocks

XFOIL keeps its state in COMMON blocks declared in three include files. This page describes
them block by block; where each variable went in yFoil is [variable and routine
names](variable-and-routine-names.md).

## XFOIL.INC - Main XFOIL Data Structures

This is the primary COMMON block include file containing global state for XFOIL.

### Dimensioning Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `IQX` | 370 | Number of surface panel nodes + 6 |
| `IWX` | IQX/8+2 = 48 | Number of wake panel nodes |
| `IPX` | 5 | Number of Qspec distributions |
| `ISX` | 2 | Number of airfoil sides |
| `IBX` | 4*IQX = 1480 | Number of buffer airfoil nodes |
| `IZX` | IQX+IWX = 418 | Number of panel nodes (airfoil + wake) |
| `IVX` | IQX/2+IWX+50 = 283 | Number of BL nodes per side |
| `NAX` | 800 | Number of points in stored polar |
| `NPX` | 12 | Number of polars |
| `NFX` | 128 | Number of reference polar points |
| `NTX` | 2*IBX = 2960 | Number of thickness/camber points |

### Geometry Arrays (COMMON/CR05/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `X(IZX)` | (418) | x-coordinates (airfoil + wake) |
| `Y(IZX)` | (418) | y-coordinates (airfoil + wake) |
| `XP(IZX)` | (418) | dX/dS spline derivatives |
| `YP(IZX)` | (418) | dY/dS spline derivatives |
| `S(IZX)` | (418) | Arc length parameter |
| `SLE` | scalar | Arc length at leading edge |
| `XLE, YLE` | scalars | Leading edge coordinates |
| `XTE, YTE` | scalars | Trailing edge coordinates |
| `CHORD` | scalar | Chord length |
| `WGAP(IWX)` | (48) | Dead air thickness in wake |
| `WAKLEN` | scalar | Wake length / chord ratio |

### Panel Method Variables

#### Vortex/Source Strengths (COMMON/CR06/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `GAM(IQX)` | (370) | Surface vortex panel strength |
| `GAMU(IQX,2)` | (370,2) | GAM for alpha = 0, 90 deg |
| `GAM_A(IQX)` | (370) | dGAM/dALFA |
| `SIG(IZX)` | (418) | Mass defect source strength |
| `NX(IZX)` | (418) | Normal x-component |
| `NY(IZX)` | (418) | Normal y-component |
| `APANEL(IZX)` | (418) | Panel angle |

#### Influence Matrices (COMMON/CR03/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `AIJ(IQX,IQX)` | (370,370) | dPsi/dGam influence matrix |
| `DIJ(IZX,IZX)` | (418,418) | dQtan/dSig influence matrix |

#### Velocities (COMMON/CR04/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `QINV(IZX)` | (418) | Inviscid tangential velocity |
| `QVIS(IZX)` | (418) | Viscous tangential velocity |
| `CPI(IZX)` | (418) | Inviscid pressure coefficient |
| `CPV(IZX)` | (418) | Viscous pressure coefficient |
| `QINVU(IZX,2)` | (418,2) | QINV for alpha = 0, 90 deg |
| `QINV_A(IZX)` | (418) | dQINV/dalpha |

### Boundary Layer Arrays (COMMON/CR15/)

| Variable | Dimension | Description | Units |
|----------|-----------|-------------|-------|
| `XSSI(IVX,ISX)` | (283,2) | BL arc length coordinate | - |
| `UEDG(IVX,ISX)` | (283,2) | BL edge velocity | - |
| `UINV(IVX,ISX)` | (283,2) | Edge velocity without mass defect | - |
| `MASS(IVX,ISX)` | (283,2) | Mass defect = Ue * delta* | - |
| `THET(IVX,ISX)` | (283,2) | Momentum thickness theta | - |
| `DSTR(IVX,ISX)` | (283,2) | Displacement thickness delta* | - |
| `TSTR(IVX,ISX)` | (283,2) | Kinetic energy thickness theta* | - |
| `CTAU(IVX,ISX)` | (283,2) | sqrt(max shear coeff) or log(amp ratio) | - |
| `TAU(IVX,ISX)` | (283,2) | Wall shear stress | - |
| `DIS(IVX,ISX)` | (283,2) | Dissipation | - |
| `CTQ(IVX,ISX)` | (283,2) | sqrt(equilibrium shear coeff) | - |
| `VTI(IVX,ISX)` | (283,2) | +/-1 panel to BL conversion | - |

### BL Indexing (COMMON/CI05/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `IBLTE(ISX)` | (2) | BL index at trailing edge |
| `NBL(ISX)` | (2) | Maximum BL array index |
| `IPAN(IVX,ISX)` | (283,2) | Panel index for BL location |
| `ISYS(IVX,ISX)` | (283,2) | Newton system line number |
| `NSYS` | scalar | Total Newton system lines |
| `ITRAN(ISX)` | (2) | BL index of transition |

### Flow State Variables (COMMON/CR09/)

| Variable | Description |
|----------|-------------|
| `ADEG, ALFA` | Angle of attack (degrees, radians) |
| `AWAKE` | AoA for wake geometry (radians) |
| `AVISC` | AoA for BL solution (radians) |
| `MVISC` | Mach number for BL solution |
| `CL, CM, CD` | Lift, moment, drag coefficients |
| `CDP, CDF` | Pressure and friction drag |
| `CL_ALF` | dCL/dALFA |
| `CL_MSQ` | dCL/d(Minf^2) |
| `PSIO` | Streamfunction inside airfoil |
| `CIRC` | Circulation |
| `COSA, SINA` | cos(ALFA), sin(ALFA) |
| `QINF` | Freestream speed (defined as 1) |
| `GAMMA, GAMM1` | Cp/Cv, Cp/Cv - 1 |
| `MINF1` | Freestream Mach at CL=1 |
| `MINF` | Current Mach number |
| `REINF1` | Reynolds number at CL=1 |
| `REINF` | Current Reynolds number |

### Newton System (COMMON/VMAT/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `VA(3,2,IZX)` | (3,2,418) | Diagonal blocks |
| `VB(3,2,IZX)` | (3,2,418) | Off-diagonal blocks |
| `VDEL(3,2,IZX)` | (3,2,418) | Residual/solution vectors |
| `VM(3,IZX,IZX)` | (3,418,418) | Mass-influence vectors |
| `VZ(3,2)` | (3,2) | TE station block |

### Convergence (COMMON/CR17/)

| Variable | Description |
|----------|-------------|
| `RMSBL` | RMS change from Newton solution |
| `RMXBL` | Max change from Newton solution |
| `RLX` | Under-relaxation factor |
| `VACCEL` | Acceleration parameter |

### Logical Flags (COMMON/CL01/)

| Flag | Description |
|------|-------------|
| `SHARP` | TRUE if trailing edge is sharp |
| `LVISC` | TRUE if viscous mode active |
| `LALFA` | TRUE if alpha specified (vs CL) |
| `LWAKE` | TRUE if wake geometry calculated |
| `LBLINI` | TRUE if BL initialized |
| `LIPAN` | TRUE if IPAN pointers calculated |
| `LQAIJ` | TRUE if AIJ computed/factored |
| `LADIJ` | TRUE if DIJ computed for airfoil |
| `LWDIJ` | TRUE if DIJ computed for wake |
| `LVCONV` | TRUE if converged BL exists |

### Stagnation Point (COMMON/CR06/)

| Variable | Description |
|----------|-------------|
| `SST` | S value at stagnation point |
| `SST_GO` | dSST/dGAM(IST) |
| `SST_GP` | dSST/dGAM(IST+1) |
| `GAMTE` | Vortex strength across TE |
| `SIGTE` | Source strength across TE |
| `GAMTE_A` | dGAMTE/dALFA |
| `SIGTE_A` | dSIGTE/dALFA |
| `DSTE` | TE panel length |
| `ANTE, ASTE` | Projected TE thickness |

### Transition Parameters (COMMON/CR15/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `ACRIT(ISX)` | (2) | log(critical amplification ratio) |
| `XSTRIP(ISX)` | (2) | Transition trip x/c locations |
| `XOCTR(ISX)` | (2) | Actual transition x/c |
| `YOCTR(ISX)` | (2) | Actual transition y/c |
| `XSSITR(ISX)` | (2) | Actual transition xi locations |

### yFoil mapping of `XFOIL.INC`

See [variable and routine names](variable-and-routine-names.md) (tables 1–4 cover this
block) and the rules in [naming](../../contributing/conventions/naming.md).

## XBL.INC - Boundary Layer Local Variables

This include file contains local variables used during boundary layer calculations. These are "station" variables representing the BL state at a single location, with derivatives for the Newton solver.

### Overview of `XBL.INC`

The BL solver uses two stations:
- **Station 1**: Upstream station (known values)
- **Station 2**: Downstream station (being solved)

Each variable at a station includes derivatives with respect to primary unknowns for implicit Newton solution.

### Parameter of `XBL.INC`

| Parameter | Value | Description |
|-----------|-------|-------------|
| `NCOM` | 73 | Number of common variables per station |

### Primary BL Variables (Station 1)

| Variable | Description |
|----------|-------------|
| `X1` | Arc length coordinate |
| `U1` | Edge velocity Ue |
| `T1` | Momentum thickness theta |
| `D1` | Displacement thickness delta* |
| `S1` | Kinetic energy thickness theta* |
| `AMPL1` | Amplification factor (laminar) or Ctau (turbulent) |
| `DW1` | Wake gap thickness |

### Primary BL Variables (Station 2)

| Variable | Description |
|----------|-------------|
| `X2` | Arc length coordinate |
| `U2` | Edge velocity Ue |
| `T2` | Momentum thickness theta |
| `D2` | Displacement thickness delta* |
| `S2` | Kinetic energy thickness theta* |
| `AMPL2` | Amplification factor (laminar) or Ctau (turbulent) |
| `DW2` | Wake gap thickness |

### Derived Variables and Derivatives (Station 1)

#### Shape Factor H (COMMON/V_VAR1/)

| Variable | Description |
|----------|-------------|
| `H1` | Shape factor H = delta*/theta |
| `H1_T1` | dH/dtheta |
| `H1_D1` | dH/d(delta*) |

#### Mach Number Correction

| Variable | Description |
|----------|-------------|
| `M1` | Local Mach number |
| `M1_U1` | dM/dUe |
| `M1_MS` | dM/d(Minf^2) |

#### Density Ratio

| Variable | Description |
|----------|-------------|
| `R1` | rho_e / rho_inf |
| `R1_U1` | dR/dUe |
| `R1_MS` | dR/d(Minf^2) |

#### Kinematic Viscosity Ratio

| Variable | Description |
|----------|-------------|
| `V1` | nu / nu_inf |
| `V1_U1` | dV/dUe |
| `V1_MS` | dV/d(Minf^2) |
| `V1_RE` | dV/dRe |

#### Kinematic Shape Factor Hk

| Variable | Description |
|----------|-------------|
| `HK1` | Kinematic shape factor |
| `HK1_U1` | dHk/dUe |
| `HK1_T1` | dHk/dtheta |
| `HK1_D1` | dHk/d(delta*) |
| `HK1_MS` | dHk/d(Minf^2) |

#### Density Shape Factor Hs

| Variable | Description |
|----------|-------------|
| `HS1` | H* (theta*/theta) |
| `HS1_U1` | dHs/dUe |
| `HS1_T1` | dHs/dtheta |
| `HS1_D1` | dHs/d(delta*) |
| `HS1_MS` | dHs/d(Minf^2) |
| `HS1_RE` | dHs/dRe |

#### Thickness Shape Factor Hc

| Variable | Description |
|----------|-------------|
| `HC1` | Hc = (delta - delta*)/theta |
| `HC1_U1` | dHc/dUe |
| `HC1_T1` | dHc/dtheta |
| `HC1_D1` | dHc/d(delta*) |
| `HC1_MS` | dHc/d(Minf^2) |

#### Reynolds Number Based on Theta

| Variable | Description |
|----------|-------------|
| `RT1` | Re_theta |
| `RT1_U1` | dRt/dUe |
| `RT1_T1` | dRt/dtheta |
| `RT1_MS` | dRt/d(Minf^2) |
| `RT1_RE` | dRt/dRe |

#### Skin Friction Coefficient

| Variable | Description |
|----------|-------------|
| `CF1` | Cf/2 |
| `CF1_U1` | dCf/dUe |
| `CF1_T1` | dCf/dtheta |
| `CF1_D1` | dCf/d(delta*) |
| `CF1_MS` | dCf/d(Minf^2) |
| `CF1_RE` | dCf/dRe |

#### Dissipation Integral

| Variable | Description |
|----------|-------------|
| `DI1` | 2*CD/H* |
| `DI1_U1` | dDI/dUe |
| `DI1_T1` | dDI/dtheta |
| `DI1_D1` | dDI/d(delta*) |
| `DI1_S1` | dDI/d(theta*) |
| `DI1_MS` | dDI/d(Minf^2) |
| `DI1_RE` | dDI/dRe |

#### Normalized Shear Stress

| Variable | Description |
|----------|-------------|
| `US1` | tau_max / (rho * Ue^2) |
| `US1_U1` | dUs/dUe |
| `US1_T1` | dUs/dtheta |
| `US1_D1` | dUs/d(delta*) |
| `US1_MS` | dUs/d(Minf^2) |
| `US1_RE` | dUs/dRe |

#### Equilibrium Shear Stress

| Variable | Description |
|----------|-------------|
| `CQ1` | sqrt(Ctau_eq) |
| `CQ1_U1` | dCq/dUe |
| `CQ1_T1` | dCq/dtheta |
| `CQ1_D1` | dCq/d(delta*) |
| `CQ1_MS` | dCq/d(Minf^2) |
| `CQ1_RE` | dCq/dRe |

#### Energy Thickness Factor

| Variable | Description |
|----------|-------------|
| `DE1` | delta** (density thickness) |
| `DE1_U1` | dDe/dUe |
| `DE1_T1` | dDe/dtheta |
| `DE1_D1` | dDe/d(delta*) |
| `DE1_MS` | dDe/d(Minf^2) |

### Station 2 Variables

Station 2 variables follow the same pattern as Station 1 but with subscript `2`:
- `H2, H2_T2, H2_D2`
- `M2, M2_U2, M2_MS`
- etc.

### Mean Skin Friction (COMMON/V_VARA/)

| Variable | Description |
|----------|-------------|
| `CFM` | Mean Cf between stations |
| `CFM_MS` | dCfm/d(Minf^2) |
| `CFM_RE` | dCfm/dRe |
| `CFM_U1, CFM_T1, CFM_D1` | Station 1 derivatives |
| `CFM_U2, CFM_T2, CFM_D2` | Station 2 derivatives |

### Transition Variables (COMMON/V_VARA/)

| Variable | Description |
|----------|-------------|
| `XT` | Transition location |
| `XT_A1` | dXt/d(amp factor) |
| `XT_MS` | dXt/d(Minf^2) |
| `XT_RE` | dXt/dRe |
| `XT_XF` | dXt/d(forced trip) |
| `XT_X1, XT_T1, XT_D1, XT_U1` | Station 1 derivatives |
| `XT_X2, XT_T2, XT_D2, XT_U2` | Station 2 derivatives |

### Global BL Variables (COMMON/V_VAR/)

| Variable | Description |
|----------|-------------|
| `DWTE` | Wake thickness at TE |
| `QINFBL` | Reference velocity for BL |
| `TKBL, TKBL_MS` | Karman-Tsien parameter |
| `RSTBL, RSTBL_MS` | Sutherland constant ratio |
| `HSTINV, HSTINV_MS` | Inverse stagnation enthalpy |
| `REYBL, REYBL_MS, REYBL_RE` | Reynolds number |
| `GAMBL` | Gamma (Cp/Cv) |
| `GM1BL` | Gamma - 1 |
| `HVRAT` | Heat/viscous ratio |
| `BULE` | BL starting location |
| `XIFORC` | Forced transition location |
| `AMCRIT` | Critical amplification |

### Flow State Flags (COMMON/V_INT/)

| Variable | Description |
|----------|-------------|
| `SIMI` | TRUE if similar (initial) profile |
| `TRAN` | TRUE if in transition interval |
| `TURB` | TRUE if turbulent |
| `WAKE` | TRUE if in wake region |
| `TRFORC` | TRUE if transition forced |
| `TRFREE` | TRUE if natural transition |
| `IDAMPV` | Damping mode for e^n |

### Local Newton System (COMMON/V_SYS/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `VS1(4,5)` | (4,5) | Station 1 Jacobian |
| `VS2(4,5)` | (4,5) | Station 2 Jacobian |
| `VSREZ(4)` | (4) | Residual vector |
| `VSR(4)` | (4) | Reynolds number sensitivity |
| `VSM(4)` | (4) | Mach number sensitivity |
| `VSX(4)` | (4) | xi sensitivity |

### Storage Arrays (COMMON/V_SAV/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `C1SAV(NCOM)` | (73) | Saved station 1 values |
| `C2SAV(NCOM)` | (73) | Saved station 2 values |

### EQUIVALENCE

The COM1 and COM2 arrays are equivalenced to the station variables:
```fortran
EQUIVALENCE (X1,COM1(1)), (X2,COM2(1))
```

This allows copying all station variables with a single array operation.

### yFoil mapping of `XBL.INC`

See [variable and routine names](variable-and-routine-names.md) (tables 1–4 cover this
block) and the rules in [naming](../../contributing/conventions/naming.md).

## BLPAR.INC - Boundary Layer Closure Parameters

This include file contains the empirical closure constants used in XFOIL's integral boundary layer formulation.

### Parameters (COMMON/BLPAR/)

| Variable | Description | Value | Source |
|----------|-------------|-------|--------|
| `SCCON` | Shear coefficient lag constant | 5.6 | Based on Green's lag equation |
| `GACON` | G-beta locus constant | 6.7 | Empirical fit |
| `GBCON` | G-beta locus constant | 0.75 | Empirical fit |
| `GCCON` | G-beta wall term constant | 18.0 | Empirical fit |
| `DLCON` | Wall/wake dissipation length ratio Lo/L | 0.9 | Empirical |
| `CTRCON` | Ctau root constant | 1.8 | Transition modeling |
| `CTRCEX` | Ctau root exponent | 3.3 | Transition modeling |
| `DUXCON` | dUe/dx amplification factor constant | 1.0 | Transition modeling |
| `CTCON` | Ctau weighting coefficient | Derived | From G-beta constants |
| `CFFAC` | Skin friction factor | 1.0 | Calibration factor |

### Usage in BL Equations

#### Shear Stress Lag (SCCON)

The lag equation for the maximum shear stress coefficient Ctau:
```
dCtau/dxi = (Ctau_eq - Ctau) * Ue / (SCCON * delta)
```

Where:
- `Ctau_eq` is the equilibrium value from the G-beta relation
- `delta` is the boundary layer thickness
- `SCCON = 5.6` controls the lag rate

#### G-Beta Relation (GACON, GBCON, GCCON)

The equilibrium Ctau is determined from the G-beta locus:
```
G = GACON * sqrt(1.0 + GBCON * beta) + GCCON / (H * Re_theta * sqrt(Cf/2))
```

Where:
- `beta` is the Clauser parameter
- The last term is a wall correction

#### Dissipation Length Ratio (DLCON)

Relates wall and wake dissipation lengths:
```
Lo = DLCON * L
```

Used in the dissipation integral formulation.

#### Ctau Root Relations (CTRCON, CTRCEX)

For turbulent flow, the Ctau correlation:
```
Ctau^(1/CTRCEX) = CTRCON * ...
```

#### Amplification Factor (DUXCON)

Modifies the e^N transition criterion based on pressure gradient:
```
dn/dRe_theta = f(H) * (1.0 + DUXCON * dUe/dx * ...)
```

### Default Values (from xfoil.f initialization)

```fortran
SCCON  = 5.6
GACON  = 6.70
GBCON  = 0.75
GCCON  = 18.0
DLCON  = 0.9
CTRCON = 1.8
CTRCEX = 3.3
DUXCON = 1.0
CTCON  = 0.5/(GACON**2 * GBCON)
CFFAC  = 1.0
```

### Physical Basis of `BLPAR.INC`

These parameters originate from:

1. **Green's lag entrainment method** (SCCON)
   - Models the finite response time of turbulence to pressure gradients

2. **Bradshaw's structural parameter** (GACON, GBCON)
   - Relates shear stress to mean velocity profile shape

3. **Skin friction correlations** (CFFAC)
   - Empirical adjustment for friction prediction

4. **Abu-Ghannam & Shaw transition** (CTRCON, CTRCEX)
   - e^N method calibration constants

### yFoil Mapping of `BLPAR.INC`

These parameters are defined in `src/bl/params.rs` (`LAG_CONSTANT`, `GBETA_LOCUS_*`, `WAKE_DISSIPATION_LENGTH_RATIO`, `TRANSITION_SQRTCTAU_*`, `LAG_PRESSURE_GRADIENT_WEIGHT`, `SQRTCTAUEQ_COEFFICIENT`, `CF_TURBULENT_FACTOR`; see the mapping table 5):

| XFOIL | yFoil |
|-------|-------|
| `SCCON` | `SCCON` |
| `GACON` | `GACON` |
| `GBCON` | `GBCON` |
| `GCCON` | `GCCON` |
| `DLCON` | `DLCON` |
| `CTRCON` | `CTRCON` |
| `CTRCEX` | `CTRCEX` |
| `DUXCON` | `DUXCON` |
| `CTCON` | `CTCON` |
| `CFFAC` | `CFFAC` |

### Notes of `BLPAR.INC`

1. These values must match XFOIL exactly for numerical agreement
2. The constants are interdependent - changing one may require adjusting others
3. CFFAC allows post-hoc calibration without changing other parameters
