# BLPAR.INC - Boundary Layer Closure Parameters

This include file contains the empirical closure constants used in XFOIL's integral boundary layer formulation.

## Parameters (COMMON/BLPAR/)

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

## Usage in BL Equations

### Shear Stress Lag (SCCON)

The lag equation for the maximum shear stress coefficient Ctau:
```
dCtau/dxi = (Ctau_eq - Ctau) * Ue / (SCCON * delta)
```

Where:
- `Ctau_eq` is the equilibrium value from the G-beta relation
- `delta` is the boundary layer thickness
- `SCCON = 5.6` controls the lag rate

### G-Beta Relation (GACON, GBCON, GCCON)

The equilibrium Ctau is determined from the G-beta locus:
```
G = GACON * sqrt(1.0 + GBCON * beta) + GCCON / (H * Re_theta * sqrt(Cf/2))
```

Where:
- `beta` is the Clauser parameter
- The last term is a wall correction

### Dissipation Length Ratio (DLCON)

Relates wall and wake dissipation lengths:
```
Lo = DLCON * L
```

Used in the dissipation integral formulation.

### Ctau Root Relations (CTRCON, CTRCEX)

For turbulent flow, the Ctau correlation:
```
Ctau^(1/CTRCEX) = CTRCON * ...
```

### Amplification Factor (DUXCON)

Modifies the e^N transition criterion based on pressure gradient:
```
dn/dRe_theta = f(H) * (1.0 + DUXCON * dUe/dx * ...)
```

## Default Values (from xfoil.f initialization)

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

## Physical Basis

These parameters originate from:

1. **Green's lag entrainment method** (SCCON)
   - Models the finite response time of turbulence to pressure gradients

2. **Bradshaw's structural parameter** (GACON, GBCON)
   - Relates shear stress to mean velocity profile shape

3. **Skin friction correlations** (CFFAC)
   - Empirical adjustment for friction prediction

4. **Abu-Ghannam & Shaw transition** (CTRCON, CTRCEX)
   - e^N method calibration constants

## YFoil Mapping

These parameters are defined in `src/bl/params.rs` (`LAG_CONSTANT`, `GBETA_LOCUS_*`, `WAKE_DISSIPATION_LENGTH_RATIO`, `TRANSITION_SQRTCTAU_*`, `LAG_PRESSURE_GRADIENT_WEIGHT`, `SQRTCTAUEQ_COEFFICIENT`, `CF_TURBULENT_FACTOR`; see the mapping table 5):

| XFOIL | YFoil |
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

## Notes

1. These values must match XFOIL exactly for numerical agreement
2. The constants are interdependent - changing one may require adjusting others
3. CFFAC allows post-hoc calibration without changing other parameters
