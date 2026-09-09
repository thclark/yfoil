# XBL.INC - Boundary Layer Local Variables

This include file contains local variables used during boundary layer calculations. These are "station" variables representing the BL state at a single location, with derivatives for the Newton solver.

## Overview

The BL solver uses two stations:
- **Station 1**: Upstream station (known values)
- **Station 2**: Downstream station (being solved)

Each variable at a station includes derivatives with respect to primary unknowns for implicit Newton solution.

## Parameter

| Parameter | Value | Description |
|-----------|-------|-------------|
| `NCOM` | 73 | Number of common variables per station |

## Primary BL Variables (Station 1)

| Variable | Description |
|----------|-------------|
| `X1` | Arc length coordinate |
| `U1` | Edge velocity Ue |
| `T1` | Momentum thickness theta |
| `D1` | Displacement thickness delta* |
| `S1` | Kinetic energy thickness theta* |
| `AMPL1` | Amplification factor (laminar) or Ctau (turbulent) |
| `DW1` | Wake gap thickness |

## Primary BL Variables (Station 2)

| Variable | Description |
|----------|-------------|
| `X2` | Arc length coordinate |
| `U2` | Edge velocity Ue |
| `T2` | Momentum thickness theta |
| `D2` | Displacement thickness delta* |
| `S2` | Kinetic energy thickness theta* |
| `AMPL2` | Amplification factor (laminar) or Ctau (turbulent) |
| `DW2` | Wake gap thickness |

## Derived Variables and Derivatives (Station 1)

### Shape Factor H (COMMON/V_VAR1/)

| Variable | Description |
|----------|-------------|
| `H1` | Shape factor H = delta*/theta |
| `H1_T1` | dH/dtheta |
| `H1_D1` | dH/d(delta*) |

### Mach Number Correction

| Variable | Description |
|----------|-------------|
| `M1` | Local Mach number |
| `M1_U1` | dM/dUe |
| `M1_MS` | dM/d(Minf^2) |

### Density Ratio

| Variable | Description |
|----------|-------------|
| `R1` | rho_e / rho_inf |
| `R1_U1` | dR/dUe |
| `R1_MS` | dR/d(Minf^2) |

### Kinematic Viscosity Ratio

| Variable | Description |
|----------|-------------|
| `V1` | nu / nu_inf |
| `V1_U1` | dV/dUe |
| `V1_MS` | dV/d(Minf^2) |
| `V1_RE` | dV/dRe |

### Kinematic Shape Factor Hk

| Variable | Description |
|----------|-------------|
| `HK1` | Kinematic shape factor |
| `HK1_U1` | dHk/dUe |
| `HK1_T1` | dHk/dtheta |
| `HK1_D1` | dHk/d(delta*) |
| `HK1_MS` | dHk/d(Minf^2) |

### Density Shape Factor Hs

| Variable | Description |
|----------|-------------|
| `HS1` | H* (theta*/theta) |
| `HS1_U1` | dHs/dUe |
| `HS1_T1` | dHs/dtheta |
| `HS1_D1` | dHs/d(delta*) |
| `HS1_MS` | dHs/d(Minf^2) |
| `HS1_RE` | dHs/dRe |

### Thickness Shape Factor Hc

| Variable | Description |
|----------|-------------|
| `HC1` | Hc = (delta - delta*)/theta |
| `HC1_U1` | dHc/dUe |
| `HC1_T1` | dHc/dtheta |
| `HC1_D1` | dHc/d(delta*) |
| `HC1_MS` | dHc/d(Minf^2) |

### Reynolds Number Based on Theta

| Variable | Description |
|----------|-------------|
| `RT1` | Re_theta |
| `RT1_U1` | dRt/dUe |
| `RT1_T1` | dRt/dtheta |
| `RT1_MS` | dRt/d(Minf^2) |
| `RT1_RE` | dRt/dRe |

### Skin Friction Coefficient

| Variable | Description |
|----------|-------------|
| `CF1` | Cf/2 |
| `CF1_U1` | dCf/dUe |
| `CF1_T1` | dCf/dtheta |
| `CF1_D1` | dCf/d(delta*) |
| `CF1_MS` | dCf/d(Minf^2) |
| `CF1_RE` | dCf/dRe |

### Dissipation Integral

| Variable | Description |
|----------|-------------|
| `DI1` | 2*CD/H* |
| `DI1_U1` | dDI/dUe |
| `DI1_T1` | dDI/dtheta |
| `DI1_D1` | dDI/d(delta*) |
| `DI1_S1` | dDI/d(theta*) |
| `DI1_MS` | dDI/d(Minf^2) |
| `DI1_RE` | dDI/dRe |

### Normalized Shear Stress

| Variable | Description |
|----------|-------------|
| `US1` | tau_max / (rho * Ue^2) |
| `US1_U1` | dUs/dUe |
| `US1_T1` | dUs/dtheta |
| `US1_D1` | dUs/d(delta*) |
| `US1_MS` | dUs/d(Minf^2) |
| `US1_RE` | dUs/dRe |

### Equilibrium Shear Stress

| Variable | Description |
|----------|-------------|
| `CQ1` | sqrt(Ctau_eq) |
| `CQ1_U1` | dCq/dUe |
| `CQ1_T1` | dCq/dtheta |
| `CQ1_D1` | dCq/d(delta*) |
| `CQ1_MS` | dCq/d(Minf^2) |
| `CQ1_RE` | dCq/dRe |

### Energy Thickness Factor

| Variable | Description |
|----------|-------------|
| `DE1` | delta** (density thickness) |
| `DE1_U1` | dDe/dUe |
| `DE1_T1` | dDe/dtheta |
| `DE1_D1` | dDe/d(delta*) |
| `DE1_MS` | dDe/d(Minf^2) |

## Station 2 Variables

Station 2 variables follow the same pattern as Station 1 but with subscript `2`:
- `H2, H2_T2, H2_D2`
- `M2, M2_U2, M2_MS`
- etc.

## Mean Skin Friction (COMMON/V_VARA/)

| Variable | Description |
|----------|-------------|
| `CFM` | Mean Cf between stations |
| `CFM_MS` | dCfm/d(Minf^2) |
| `CFM_RE` | dCfm/dRe |
| `CFM_U1, CFM_T1, CFM_D1` | Station 1 derivatives |
| `CFM_U2, CFM_T2, CFM_D2` | Station 2 derivatives |

## Transition Variables (COMMON/V_VARA/)

| Variable | Description |
|----------|-------------|
| `XT` | Transition location |
| `XT_A1` | dXt/d(amp factor) |
| `XT_MS` | dXt/d(Minf^2) |
| `XT_RE` | dXt/dRe |
| `XT_XF` | dXt/d(forced trip) |
| `XT_X1, XT_T1, XT_D1, XT_U1` | Station 1 derivatives |
| `XT_X2, XT_T2, XT_D2, XT_U2` | Station 2 derivatives |

## Global BL Variables (COMMON/V_VAR/)

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

## Flow State Flags (COMMON/V_INT/)

| Variable | Description |
|----------|-------------|
| `SIMI` | TRUE if similar (initial) profile |
| `TRAN` | TRUE if in transition interval |
| `TURB` | TRUE if turbulent |
| `WAKE` | TRUE if in wake region |
| `TRFORC` | TRUE if transition forced |
| `TRFREE` | TRUE if natural transition |
| `IDAMPV` | Damping mode for e^n |

## Local Newton System (COMMON/V_SYS/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `VS1(4,5)` | (4,5) | Station 1 Jacobian |
| `VS2(4,5)` | (4,5) | Station 2 Jacobian |
| `VSREZ(4)` | (4) | Residual vector |
| `VSR(4)` | (4) | Reynolds number sensitivity |
| `VSM(4)` | (4) | Mach number sensitivity |
| `VSX(4)` | (4) | xi sensitivity |

## Storage Arrays (COMMON/V_SAV/)

| Variable | Dimension | Description |
|----------|-----------|-------------|
| `C1SAV(NCOM)` | (73) | Saved station 1 values |
| `C2SAV(NCOM)` | (73) | Saved station 2 values |

## EQUIVALENCE

The COM1 and COM2 arrays are equivalenced to the station variables:
```fortran
EQUIVALENCE (X1,COM1(1)), (X2,COM2(1))
```

This allows copying all station variables with a single array operation.

## yFoil mapping

See [`xfoil-to-yfoil-mapping.md`](../xfoil-to-yfoil-mapping.md) (tables 1–4 cover this
block) and the rules in [`docs/conventions/naming.md`](../../conventions/naming.md).

