# XFOIL Architecture

This document describes the overall architecture of XFOIL's viscous-inviscid coupling method.

## Overview

XFOIL solves the 2D viscous airfoil problem using:

1. **Panel Method** - Linear vortex panels for inviscid flow
2. **Integral Boundary Layer** - 3-equation integral BL with e^N transition
3. **Viscous-Inviscid Coupling** - Semi-inverse Newton iteration

## Solution Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                        INPUT                                    │
│  Geometry: X, Y coordinates                                     │
│  Conditions: α (or CL), Re, M∞                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                     GEOMETRY SETUP                              │
│  PANGEN: Repanel airfoil to N panels                            │
│  SPLINE: Compute spline coefficients for X(S), Y(S)             │
│  TECALC: Calculate TE geometry (sharp or blunt)                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   INVISCID SOLUTION                             │
│  PSILIN: Compute vortex influence coefficients AIJ              │
│  SETUP: Factor AIJ matrix (LU decomposition)                    │
│  GAMCALC: Solve [AIJ]{γ} = {RHS} for vortex strengths           │
│  QINV: Tangential velocity from vortex distribution             │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      VISCAL                                     │
│  Main viscous-inviscid coupling loop                            │
└─────────────────────────────────────────────────────────────────┘
         │                    │                    │
         ▼                    ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│   WAKE SETUP    │  │   DIJ MATRIX    │  │   BL SETUP      │
│                 │  │                 │  │                 │
│ XYWAKE: Compute │  │ QDCALC: Source  │  │ STFIND: Locate  │
│ wake trajectory │  │ influence on    │  │ stagnation pt   │
│ from TE         │  │ tangential      │  │                 │
│                 │  │ velocity        │  │ IBLPAN: Set     │
│ QWCALC: Set     │  │                 │  │ BL→panel ptrs   │
│ wake velocities │  │ DIJ(i,j) =      │  │                 │
│                 │  │ ∂Qtan_i/∂σ_j    │  │ XICALC: Arc     │
└─────────────────┘  └─────────────────┘  │ length coords   │
                                          └─────────────────┘
                              │
                              ▼
         ┌────────────────────────────────────────┐
         │         NEWTON ITERATION LOOP          │
         │         (ITER = 1 to ITMAX)            │
         └────────────────────────────────────────┘
                              │
         ┌────────────────────┴────────────────────┐
         ▼                                         ▼
┌─────────────────────────┐           ┌─────────────────────────┐
│       SETBL             │           │       BLSOLV            │
│                         │           │                         │
│ Fill Newton system      │           │ Block tri-diagonal      │
│ for BL variables        │           │ solution with           │
│                         │           │ underrelaxation         │
│ For each BL station:    │           │                         │
│ ├─ BLPRV: Primary vars  │           │ Solve:                  │
│ ├─ BLKIN: Secondary     │           │ [VA VB] [Δ] = [VDEL]    │
│ ├─ BLVAR: Closures      │           │                         │
│ ├─ BLSYS: Local system  │           └─────────────────────────┘
│ └─ TRCHEK: Transition   │
└─────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        UPDATE                                   │
│  Apply Newton deltas with underrelaxation                       │
│  Update: CTAU, THET, DSTR, UEDG, MASS                           │
│  Enforce bounds: Hk > 1, Ctau > 0                               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    VELOCITY UPDATE                              │
│  QVFUE: Compute QVIS from mass defect                           │
│         QVIS_i = QINV_i + Σ_j DIJ(i,j)·σ_j                      │
│  GAMQV: Set γ from QVIS                                         │
│  STMOVE: Relocate stagnation point                              │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  FORCE CALCULATION                              │
│  CLCALC: Integrate Cp → CL, CM                                  │
│  CDCALC: CD from momentum deficit                               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │ RMSBL < 1e-4 ?  │
                    └─────────────────┘
                     │ No          │ Yes
                     │             │
         ┌───────────┘             └───────────┐
         ▼                                     ▼
    Next iteration                      ┌───────────┐
                                        │ CONVERGED │
                                        └───────────┘
                                               │
                                               ▼
┌─────────────────────────────────────────────────────────────────┐
│                          OUTPUT                                 │
│  CL, CD, CM coefficients                                        │
│  Cp distribution                                                │
│  BL parameters: δ*, θ, Cf, H                                    │
│  Transition locations                                           │
└─────────────────────────────────────────────────────────────────┘
```

## Boundary Layer Solution Detail

```
SETBL / MRCHDU Flow:

Station 1 (upstream)              Station 2 (downstream)
       │                                   │
       ▼                                   ▼
┌─────────────┐                    ┌─────────────┐
│   BLPRV     │                    │   BLPRV     │
│ Set primary │                    │ Set primary │
│ variables   │                    │ variables   │
└─────────────┘                    └─────────────┘
       │                                   │
       ▼                                   ▼
┌─────────────┐                    ┌─────────────┐
│   BLKIN     │                    │   BLKIN     │
│ M, ρ, Hk    │                    │ M, ρ, Hk    │
│ Re_θ        │                    │ Re_θ        │
└─────────────┘                    └─────────────┘
       │                                   │
       ▼                                   ▼
┌─────────────┐                    ┌─────────────┐
│   BLVAR     │                    │   BLVAR     │
│ H*, Cf, DI  │                    │ H*, Cf, DI  │
│ Closures    │                    │ Closures    │
└─────────────┘                    └─────────────┘
       │                                   │
       └──────────────┬────────────────────┘
                      ▼
              ┌─────────────┐
              │   BLSYS     │
              │ Assemble    │
              │ 4×4 system  │
              └─────────────┘
                      │
                      ▼
              ┌─────────────┐
              │  TRCHEK     │
              │ Transition? │
              └─────────────┘
                     │
           ┌─────────┴─────────┐
           ▼                   ▼
      [Laminar]           [Turbulent]
      Amplification       Shear stress
      equation            lag equation
```

## Key Data Structures

### Panel Method
```
AIJ(N,N)     Vortex influence matrix: dψ/dγ
DIJ(N+NW,N+NW)  Source influence matrix: dQtan/dσ
GAM(N)       Vortex strengths
SIG(N+NW)    Source strengths (mass defect)
QINV(N+NW)   Inviscid tangential velocity
QVIS(N+NW)   Viscous tangential velocity
```

### Boundary Layer
```
UEDG(IVX,2)  Edge velocity
DSTR(IVX,2)  Displacement thickness δ*
THET(IVX,2)  Momentum thickness θ
CTAU(IVX,2)  Shear coefficient or amplification
MASS(IVX,2)  Mass defect = Ue·δ*
```

### Newton System
```
VA(3,2,N)    Diagonal blocks
VB(3,2,N)    Off-diagonal blocks
VDEL(3,2,N)  Residual/solution vectors
VM(3,N,N)    Mass influence vectors
```

## Coupling Mechanism

The viscous-inviscid interaction is handled through:

1. **Mass defect equivalence**:
   ```
   σ = d(ρ·Ue·δ*)/ds
   ```
   BL displacement thickness creates equivalent source distribution.

2. **Source influence**:
   ```
   Qvis = Qinv + DIJ·σ
   ```
   Sources modify surface velocity.

3. **Edge velocity feedback**:
   ```
   Ue = |Qvis|
   ```
   Modified velocity becomes new BL edge condition.

4. **Newton iteration**:
   Simultaneous update of BL variables and flow field.

## Transition Model

The e^N method:
```
dN/ds = f(Hk, Re_θ)  [Drela-Giles correlation]

Transition when: N > N_crit (typically 9)
```

The amplification rate f depends on:
- Kinematic shape factor Hk
- Momentum thickness Reynolds number Re_θ
- Pressure gradient (implicit in Hk)
