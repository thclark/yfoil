# XFOIL.INC - Main XFOIL Data Structures

This is the primary COMMON block include file containing global state for XFOIL.

## Dimensioning Parameters

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

## Geometry Arrays (COMMON/CR05/)

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

## Panel Method Variables

### Vortex/Source Strengths (COMMON/CR06/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `GAM(IQX)` | (370) | Surface vortex panel strength |
| `GAMU(IQX,2)` | (370,2) | GAM for alpha = 0, 90 deg |
| `GAM_A(IQX)` | (370) | dGAM/dALFA |
| `SIG(IZX)` | (418) | Mass defect source strength |
| `NX(IZX)` | (418) | Normal x-component |
| `NY(IZX)` | (418) | Normal y-component |
| `APANEL(IZX)` | (418) | Panel angle |

### Influence Matrices (COMMON/CR03/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `AIJ(IQX,IQX)` | (370,370) | dPsi/dGam influence matrix |
| `DIJ(IZX,IZX)` | (418,418) | dQtan/dSig influence matrix |

### Velocities (COMMON/CR04/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `QINV(IZX)` | (418) | Inviscid tangential velocity |
| `QVIS(IZX)` | (418) | Viscous tangential velocity |
| `CPI(IZX)` | (418) | Inviscid pressure coefficient |
| `CPV(IZX)` | (418) | Viscous pressure coefficient |
| `QINVU(IZX,2)` | (418,2) | QINV for alpha = 0, 90 deg |
| `QINV_A(IZX)` | (418) | dQINV/dalpha |

## Boundary Layer Arrays (COMMON/CR15/)

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

## BL Indexing (COMMON/CI05/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `IBLTE(ISX)` | (2) | BL index at trailing edge |
| `NBL(ISX)` | (2) | Maximum BL array index |
| `IPAN(IVX,ISX)` | (283,2) | Panel index for BL location |
| `ISYS(IVX,ISX)` | (283,2) | Newton system line number |
| `NSYS` | scalar | Total Newton system lines |
| `ITRAN(ISX)` | (2) | BL index of transition |

## Flow State Variables (COMMON/CR09/)

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

## Newton System (COMMON/VMAT/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `VA(3,2,IZX)` | (3,2,418) | Diagonal blocks |
| `VB(3,2,IZX)` | (3,2,418) | Off-diagonal blocks |
| `VDEL(3,2,IZX)` | (3,2,418) | Residual/solution vectors |
| `VM(3,IZX,IZX)` | (3,418,418) | Mass-influence vectors |
| `VZ(3,2)` | (3,2) | TE station block |

## Convergence (COMMON/CR17/)

| Variable | Description |
|----------|-------------|
| `RMSBL` | RMS change from Newton solution |
| `RMXBL` | Max change from Newton solution |
| `RLX` | Under-relaxation factor |
| `VACCEL` | Acceleration parameter |

## Logical Flags (COMMON/CL01/)

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

## Stagnation Point (COMMON/CR06/)

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

## Transition Parameters (COMMON/CR15/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `ACRIT(ISX)` | (2) | log(critical amplification ratio) |
| `XSTRIP(ISX)` | (2) | Transition trip x/c locations |
| `XOCTR(ISX)` | (2) | Actual transition x/c |
| `YOCTR(ISX)` | (2) | Actual transition y/c |
| `XSSITR(ISX)` | (2) | Actual transition xi locations |

## YFoil Mapping

| XFOIL Variable | YFoil Equivalent |
|---------------|------------------|
| `X, Y, S` | `Geometry.x, y, s` |
| `GAM` | `PanelSolution.gamma` |
| `SIG` | `BLState.mass_defect` |
| `QINV, QVIS` | `Velocities.q_inv, q_vis` |
| `DSTR, THET` | `BLStation.delta_star, theta` |
| `UEDG` | `BLStation.ue` |
| `CL, CD, CM` | `ForceCoefficients.cl, cd, cm` |
