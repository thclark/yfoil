---
icon: lucide/library
title: XFOIL reference
---

# XFOIL reference

Everything recorded about XFOIL 6.99 during the translation, on one page: the
architecture and data flow, the instrumentation, the variable and control-flow
map, the XFOIL → yFoil name mapping, and a walk through every translated module
and COMMON block.


## XFOIL Reference Documentation

This directory contains detailed documentation of XFOIL 6.99 for the purpose of:

1. Validating yFoil produces numerically identical results
2. Understanding the algorithm flow for debugging

### Overview

XFOIL (created by Mark Drela at MIT, 1980s) performs viscous-inviscid analysis of 2D airfoil sections using:

- Panel method with vortex distribution for inviscid flow
- Integral boundary layer solver with e^N transition prediction
- Coupled viscous-inviscid iteration (VISCAL)

### Directory Structure

```
xfoil-reference/
├── index.md                  # This page: every reference document, in one place
├── architecture.mmd          # Mermaid diagram source
└── architecture.pdf          # Generated PDF (see below)
```

The reference documents that used to be one file each — the architecture, the
instrumentation guide, the variable map, the XFOIL → yFoil mapping, the module
walkthroughs (`xoper`, `xblsys`, `xbl`, `xpanel`, `spline`) and the COMMON block
descriptions (`XFOIL.INC`, `XBL.INC`, `BLPAR.INC`) — are the sections of this
page.

### XFOIL Source Files

| File       | Purpose                      | Key Subroutines                     |
|------------|------------------------------|-------------------------------------|
| `xoper.f`  | Main operations and coupling | VISCAL, SPECAL, SPECCL              |
| `xblsys.f` | BL Newton system             | BLPRV, BLKIN, BLVAR, BLSYS, TRCHEK2 |
| `xbl.f`    | BL marching                  | SETBL, MRCHUE, MRCHDU, UPDATE       |
| `xpanel.f` | Panel method                 | PSILIN, QDCALC, UICALC, QVFUE       |
| `spline.f` | Cubic splines                | SPLINE, SEVAL, SPLIND               |
| `xfoil.f`  | Panelling                    | PANGEN, GETPAN (PPAR), TECALC       |
| `xgeom.f`  | Geometry                     | LEFIND, SCALC, SEGSPL               |
| `naca.f`   | NACA generation              | NACA4, NACA5                        |

### Data Flow

```
Input Geometry (X, Y coordinates)
        │
        ▼
┌───────────────────┐
│   PANGEN          │  Repanel airfoil
│   (xgeom.f)       │
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   SPECAL          │  Specify alpha
│   (xoper.f)       │
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   Panel Method    │  Solve inviscid flow
│   (xpanel.f)      │  → GAM, QINV
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   VISCAL          │  Viscous-inviscid coupling
│   (xoper.f)       │
│   ┌─────────────┐ │
│   │ SETBL       │ │  Initialize BL
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ MRCHUE      │ │  March BL equations
│   │ (xbl.f)     │ │  → DSTR, THET
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ UPDATE      │ │  Mass defect → source
│   │ (xbl.f)     │ │  → SIG
│   └─────────────┘ │
│   ┌─────────────┐ │
│   │ QDCALC      │ │  Update Qvis
│   └─────────────┘ │
│         ↑         │
│         └─────────│  Iterate until converged
└───────────────────┘
        │
        ▼
┌───────────────────┐
│   Force Coeffs    │  CL, CD, CM
└───────────────────┘
```

### yFoil equivalents

The full XFOIL → yFoil name mapping is [the mapping section below](#xfoil-yfoil-mapping);
the naming rules are [`docs/conventions/naming.md`](conventions/naming.md). Where the
translated routines live:

| XFOIL             | yFoil                                                  |
|-------------------|--------------------------------------------------------|
| `xoper.f:VISCAL`  | `src/solver/viscal.rs::solve_viscous`                  |
| `xoper.f:SPECAL`  | `src/solver/specal.rs::solve_inviscid_at_alpha`        |
| `xbl.f:SETBL`     | `src/solver/setbl.rs::assemble_newton_system`          |
| `xbl.f:UPDATE`    | `src/solver/update.rs::apply_newton_update`            |
| `xbl.f:MRCHUE`    | `src/bl/mrchue.rs::march_direct`                       |
| `xbl.f:MRCHDU`    | `src/bl/mrchdu.rs::march_prescribed_dstar`             |
| `xblsys.f:BLKIN`  | `src/bl/station.rs::StationState::set_kinematic_variables` |
| `xblsys.f:BLVAR`  | `src/bl/station.rs::StationState::set_closure_variables`   |
| `xblsys.f:BLSYS`  | `src/bl/blsys.rs::assemble_interval_system`            |
| `xblsys.f:BLDIF`  | `src/bl/difference.rs::IntervalSystem::assemble_interval_equations` |
| `xblsys.f:TRCHEK2`| `src/bl/transition.rs::check_transition`               |
| `xsolve.f:BLSOLV` | `src/bl/blsolv.rs::solve_newton_system`                |
| `xpanel.f:PSILIN` | `src/solver/psilin.rs::panel_influence`                |
| `xpanel.f:GGCALC` | `src/solver/ggcalc.rs::build_inviscid_system`          |

Every translating function carries `#[doc(alias = "XFOIL NAME")]`, so `cargo doc` search by the
Fortran name finds it.

### Generating Architecture PDF

To generate the architecture diagram PDF from the Mermaid source:

```bash
# Install mermaid-cli if not present
npm install -g @mermaid-js/mermaid-cli

# Generate PDF
cd docs/xfoil-reference
mmdc -i architecture.mmd -o architecture.pdf -b white
```

Alternatively, you can use the online Mermaid Live Editor at https://mermaid.live/ to render and export the diagram.

### Usage

This documentation serves as a reference when:

1. Debugging mismatches between yFoil and XFOIL
2. Understanding algorithmic details for implementation
3. Writing test fixtures and validation cases
4. Preparing paper content describing the method

## XFOIL Architecture

This document describes the overall architecture of XFOIL's viscous-inviscid coupling method.

### Overview

XFOIL solves the 2D viscous airfoil problem using:

1. **Panel Method** - Linear vortex panels for inviscid flow
2. **Integral Boundary Layer** - 3-equation integral BL with e^N transition
3. **Viscous-Inviscid Coupling** - Semi-inverse Newton iteration

### Solution Flow

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

### Boundary Layer Solution Detail

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

### Key Data Structures

#### Panel Method
```
AIJ(N,N)     Vortex influence matrix: dψ/dγ
DIJ(N+NW,N+NW)  Source influence matrix: dQtan/dσ
GAM(N)       Vortex strengths
SIG(N+NW)    Source strengths (mass defect)
QINV(N+NW)   Inviscid tangential velocity
QVIS(N+NW)   Viscous tangential velocity
```

#### Boundary Layer
```
UEDG(IVX,2)  Edge velocity
DSTR(IVX,2)  Displacement thickness δ*
THET(IVX,2)  Momentum thickness θ
CTAU(IVX,2)  Shear coefficient or amplification
MASS(IVX,2)  Mass defect = Ue·δ*
```

#### Newton System
```
VA(3,2,N)    Diagonal blocks
VB(3,2,N)    Off-diagonal blocks
VDEL(3,2,N)  Residual/solution vectors
VM(3,N,N)    Mass influence vectors
```

### Coupling Mechanism

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

### Transition Model

The e^N method:
```
dN/ds = f(Hk, Re_θ)  [Drela-Giles correlation]

Transition when: N > N_crit (typically 9)
```

The amplification rate f depends on:
- Kinematic shape factor Hk
- Momentum thickness Reynolds number Re_θ
- Pressure gradient (implicit in Hk)

## XFOIL Instrumentation Guide

This document describes how to instrument XFOIL to capture intermediate values for validation against yFoil.

### Current Instrumentation

The XFOIL source code in `xfoil/xfoil6.99/src/` has been modified with diagnostic output:

#### xoper.f - VISCAL Coupling

Location: `xoper.f:2982-3104`

Outputs:
- `/tmp/xfoil_dij.dat` - Source influence matrix DIJ
- `/tmp/xfoil_viscal_iter.dat` - Iteration-by-iteration BL state

#### xfoil.f - Panel Geometry

Location: `xfoil.f:2124`

Outputs:
- `/tmp/xfoil_panels.dat` - Panel coordinates and normals

### Output Format

#### DIJ Matrix (`/tmp/xfoil_dij.dat`)

```
=== DIJ MATRIX ===
N = 160
NW = 23
    1     1  1.234567890123456E+00
    1     2  2.345678901234567E-01
    ...
```

Format: `I J VALUE` where VALUE is in Fortran E24.16 format.

#### VISCAL Iteration Log (`/tmp/xfoil_viscal_iter.dat`)

```
=== VISCAL ITERATION LOG ===
--- ITERATION 1 ---
ALFA =  0.000000000000000E+00
CL =  1.234567890123456E-01
CD =  9.876543210987654E-03
CDF =  5.432109876543210E-03
CDP =  4.444433333222211E-03
CM = -1.234567890123456E-02
RMSBL =  1.234567890123456E-02
RMXBL =  9.876543210987654E-02
RLX =  1.000000000000000E+00
IST = 80
NBL(1) = 84
NBL(2) = 84
ITRAN(1) = 45
ITRAN(2) = 52
--- Upper surface BL (IS=1) ---
IBL, UEDG, DSTR, THET, MASS, CTAU
    1  1.234567890E+00  1.234567890E-03  ...
    2  ...
--- Lower surface BL (IS=2) ---
...
```

### Adding New Instrumentation

#### Step 1: Locate Subroutine

Find the subroutine in the XFOIL source:

```bash
grep -n "SUBROUTINE BLKIN" xfoil/xfoil6.99/src/xblsys.f
```

#### Step 2: Add WRITE Statements

Example for BLKIN output:

```fortran
C---- INSTRUMENTATION: Output BLKIN variables
      IF (IPRINT .GT. 0) THEN
        WRITE(IPRINT,'(A)') '=== BLKIN OUTPUT ==='
        WRITE(IPRINT,'(A,E24.16)') 'X2 =', X2
        WRITE(IPRINT,'(A,E24.16)') 'U2 =', U2
        WRITE(IPRINT,'(A,E24.16)') 'T2 =', T2
        WRITE(IPRINT,'(A,E24.16)') 'D2 =', D2
        WRITE(IPRINT,'(A,E24.16)') 'H2 =', H2
        WRITE(IPRINT,'(A,E24.16)') 'HK2 =', HK2
        WRITE(IPRINT,'(A,E24.16)') 'RT2 =', RT2
        WRITE(IPRINT,'(A,E24.16)') 'M2 =', M2
        WRITE(IPRINT,'(A,E24.16)') 'R2 =', R2
        WRITE(IPRINT,'(A,E24.16)') 'V2 =', V2
      ENDIF
```

#### Step 3: Rebuild XFOIL

```bash
cd xfoil/xfoil6.99/bin
make clean
make
```

#### Step 4: Run Test Case

```bash
./xfoil << EOF
PLOP
G F

NACA 0012
OPER
VISC 1e6
ITER 20
ALFA 2
EOF
```

#### Step 5: Parse Output

Use the fixture generation script to convert output to JSON.

### Key Subroutines to Instrument

#### Panel Method

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| PSILIN | xpanel.f | PSI, DZDG, DZDM per panel |
| QDCALC | xpanel.f | DIJ matrix, CIJ matrix |

#### BL System

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| BLPRV | xblsys.f | X2, U2, T2, D2, U2_UEI, U2_MS |
| BLKIN | xblsys.f | M2, R2, H2, HK2, RT2, V2 + derivatives |
| BLVAR | xblsys.f | HS2, CF2, DI2, US2, HC2, DE2 + derivatives |
| BLSYS | xblsys.f | VS1, VS2, VSREZ matrices |
| TRCHEK2 | xblsys.f | AMPL2, XT, transition derivatives |

#### BL Marching

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| SETBL | xbl.f | VA, VB, VDEL, VM matrices |
| MRCHUE | xbl.f | Per-station THET, DSTR, CTAU |
| UPDATE | xbl.f | RMSBL, RMXBL, RLX, changes |

### Precision Requirements

**All output must use E24.16 format** to capture full double-precision values:

```fortran
WRITE(LU,'(A,E24.16)') 'VAR =', VAR
```

This gives 16 significant figures, matching IEEE double precision.

### Comparison Tolerance

When comparing yFoil to XFOIL:

| Tolerance | Meaning |
|-----------|---------|
| < 1e-14 | Perfect match (roundoff only) |
| 1e-14 to 1e-10 | Acceptable (numerical precision) |
| 1e-10 to 1e-6 | Concerning (investigate) |
| > 1e-6 | Bug (must fix) |

Relative error calculation:
```
rel_error = |yfoil - xfoil| / max(|xfoil|, 1e-14)
```

### Existing Debug Files

The instrumented XFOIL writes to:
- `/tmp/xfoil_dij.dat`
- `/tmp/xfoil_viscal_iter.dat`
- `/tmp/xfoil_inviscid.dat`
- `/tmp/xfoil_panels.dat`
- `/tmp/xfoil_bl_debug.dat` (from xbl.f MRCHUE)

These files are overwritten on each run. Copy them before running another case.

## XFOIL Variable and Control Flow Map

This document maps the variables updated within each subroutine and common block
during a typical XFOIL run (NACA 0012, repanel, viscous solve at alpha=0, alpha=1).

Generated from instrumented XFOIL run on 2026-01-17.

### Overview of Key Common Blocks

| Common Block | File      | Key Variables                                    | Purpose                                    |
|--------------|-----------|--------------------------------------------------|--------------------------------------------|
| CR03         | XFOIL.INC | AIJ, DIJ                                         | Influence matrices (dPsi/dGam, dQtan/dSig) |
| CR04         | XFOIL.INC | QINV, QVIS, CPI, CPV, QINVU                      | Flow velocities and Cp                     |
| CR05         | XFOIL.INC | X, Y, XP, YP, S, SLE, XLE, YLE, XTE, YTE, CHORD  | Panel geometry                             |
| CR06         | XFOIL.INC | GAM, GAMU, GAM_A, SIG, NX, NY, APANEL, SST       | Vortex/source strengths                    |
| CR09         | XFOIL.INC | ALFA, CL, CM, CD, CDF, CDP, MINF, QINF           | Flow parameters                            |
| CR15         | XFOIL.INC | XSSI, UEDG, UINV, MASS, THET, DSTR, CTAU         | BL arrays (2 sides)                        |
| CI04         | XFOIL.INC | N, NB, NW, NPAN, IST, ITMAX                      | Panel/BL indices                           |
| CI05         | XFOIL.INC | IBLTE, NBL, IPAN, ISYS, NSYS, ITRAN              | BL indices                                 |
| CL01         | XFOIL.INC | LVISC, LALFA, LWAKE, LBLINI, LIPAN, LQAIJ, LADIJ | State flags                                |
| VMAT         | XFOIL.INC | VA, VB, VDEL, VM, VZ                             | BL Newton system matrices                  |
| COM1/COM2    | XBL.INC   | Local BL state (73 variables per side)           | BL marching state                          |
| BLPAR        | BLPAR.INC | SCCON, GACON, GBCON, GCCON, DLCON, CTRCON        | BL closure constants                       |

### Control Flow: NACA 0012 Generation

```
NACA(0012)
├─ NACA4(IDES=12)
│  ├─ Computes XX(i), YT(i), YC(i) thickness/camber distributions
│  └─ Generates XB(1:NB), YB(1:NB) buffer airfoil coords
│
├─ SCALC(XB, YB, SB, NB)
│  └─ Computes SB(i) = arc length parameter
│
├─ SEGSPL(XB, XBP, SB, NB)
│  └─ Computes XBP = dXB/dSB spline derivatives
│
├─ SEGSPL(YB, YBP, SB, NB)
│  └─ Computes YBP = dYB/dSB spline derivatives
│
├─ GEOPAR(...)
│  └─ Computes SBLE, CHORDB, AREAB, RADBLE, ANGBTE, THICKB, CAMBRB
│
└─ PANGEN(.TRUE.)
   └─ [see PANGEN below]
```

#### Variables Updated by NACA:

| Variable  | Common Block | Description                    |
|-----------|--------------|--------------------------------|
| XB(1:NB)  | CR14         | Buffer airfoil x-coords        |
| YB(1:NB)  | CR14         | Buffer airfoil y-coords        |
| SB(1:NB)  | CR14         | Buffer arc length parameter    |
| XBP(1:NB) | CR14         | dXB/dSB                        |
| YBP(1:NB) | CR14         | dYB/dSB                        |
| NB        | CI04         | Number of buffer points (=245) |
| SBLE      | CR14         | LE arc length on buffer        |
| CHORDB    | CR14         | Buffer chord                   |
| NAME      | CC01         | Airfoil name = "NACA0012"      |

### Control Flow: PANGEN (Paneling)

```
PANGEN(SHOPAR)
├─ SCALC(XB, YB, SB, NB)
├─ SEGSPL(XB, XBP, SB, NB)
├─ SEGSPL(YB, YBP, SB, NB)
├─ LEFIND(SBLE, ...)  → Finds LE arc length
├─ CURV(...) loop     → Computes curvature array W5
├─ TRISOL(...)        → Solves for smoothed curvature
├─ Node distribution iteration (curvature-based bunching)
├─ SEVAL(...) loop    → Extracts panel coords X, Y, S from splines
├─ SEGSPL(X, XP, S, N)
├─ SEGSPL(Y, YP, S, N)
├─ LEFIND(SLE, X, XP, Y, YP, S, N)
├─ SCALC(X, Y, S, N)
├─ NCALC(...)         → Computes normals NX, NY
└─ Output: N panels with X, Y, S, XP, YP, NX, NY, APANEL
```

#### Variables Updated by PANGEN:

| Variable    | Common Block | Description                      |
|-------------|--------------|----------------------------------|
| N           | CI04         | Number of panel nodes (=160)     |
| X(1:N)      | CR05         | Panel x-coordinates              |
| Y(1:N)      | CR05         | Panel y-coordinates              |
| S(1:N)      | CR05         | Panel arc length                 |
| XP(1:N)     | CR05         | dX/dS                            |
| YP(1:N)     | CR05         | dY/dS                            |
| NX(1:N)     | CR06         | Normal x-component               |
| NY(1:N)     | CR06         | Normal y-component               |
| APANEL(1:N) | CR06         | Panel angle                      |
| SLE         | CR05         | LE arc length                    |
| XLE, YLE    | CR05         | LE coordinates                   |
| XTE, YTE    | CR05         | TE coordinates                   |
| CHORD       | CR05         | LE-point-to-TE-midpoint distance (slightly under 1 for a NACA section) |
| SHARP       | CL01         | .TRUE. if sharp TE               |

### Control Flow: SPECAL (Inviscid Solve at Alpha)

```
SPECAL
├─ IF NOT LGAMU OR NOT LQAIJ:
│  └─ GGCALC
│     ├─ Build AIJ matrix (dPsi/dGam):
│     │  └─ PSILIN(I, X(I), Y(I), ...) for each panel
│     ├─ LUDCMP(AIJ, ...)  → LU factor AIJ
│     ├─ BAKSUB(AIJ, GAMU(:,1))  → Solve for alpha=0
│     ├─ BAKSUB(AIJ, GAMU(:,2))  → Solve for alpha=90
│     ├─ QINVU = GAMU
│     └─ Set LGAMU=.TRUE., LQAIJ=.TRUE.
│
├─ Compute GAM = cos(α)*GAMU(:,1) + sin(α)*GAMU(:,2)
├─ Compute GAM_A = -sin(α)*GAMU(:,1) + cos(α)*GAMU(:,2)
├─ TECALC → TE vortex/source GAMTE, SIGTE
├─ QISET → QINV from GAM
│
├─ Newton loop for Mach correction (if MATYP != 1):
│  ├─ MRCL(CLM, MINF_CLM, REINF_CLM)
│  ├─ COMSET → TKLAM, TKL_MSQ, CPSTAR, QSTAR
│  └─ CLCALC → CL, CM, CDP, CL_ALF, CL_MSQ
│
├─ CPCALC(N, QINV, QINF, MINF, CPI) → Cp from QINV
└─ Return with inviscid solution
```

#### Variables Updated by SPECAL:

| Variable       | Common Block | Description                  |
|----------------|--------------|------------------------------|
| GAMU(1:N+1,1)  | CR06         | Vorticity at alpha=0         |
| GAMU(1:N+1,2)  | CR06         | Vorticity at alpha=90        |
| GAM(1:N)       | CR06         | Vorticity at current alpha   |
| GAM_A(1:N)     | CR06         | dGAM/dAlpha                  |
| QINV(1:N)      | CR04         | Inviscid tangential velocity |
| QINVU(1:N,1:2) | CR04         | QINV at alpha=0, 90          |
| CPI(1:N)       | CR04         | Inviscid Cp                  |
| CL, CM, CDP    | CR09         | Lift, moment, pressure drag  |
| CL_ALF, CL_MSQ | CR09         | Derivatives                  |
| GAMTE, SIGTE   | CR06         | TE panel strengths           |
| AIJ(N+1,N+1)   | CR03         | Influence matrix (factored)  |
| LGAMU, LQAIJ   | CL01         | State flags = .TRUE.         |

### Control Flow: VISCAL (Viscous-Inviscid Coupling)

```
VISCAL(NITER1)
├─ IF NOT LWAKE:
│  └─ XYWAKE → Wake trajectory X, Y for N+1:N+NW
│
├─ QWCALC → Wake velocities QINVU(N+1:N+NW, 1:2)
├─ QISET → QINV from GAMU at current alpha
│
├─ IF NOT LIPAN:
│  ├─ IF LBLINI: GAMQV → GAM from QVIS
│  ├─ STFIND → SST, IST (stagnation point)
│  ├─ IBLPAN → IPAN(:,1:2) (BL→panel pointers)
│  ├─ XICALC → XSSI(:,1:2) (BL arc lengths)
│  └─ IBLSYS → ISYS(:,1:2), NSYS (BL→system pointers)
│
├─ UICALC → UINV from QINV
│
├─ IF NOT LBLINI:
│  └─ Initialize UEDG = UINV
│
├─ IF LVCONV:
│  ├─ QVFUE → QVIS from UEDG + mass defect
│  ├─ CPCALC → CPV, CPI
│  ├─ GAMQV → GAM from QVIS
│  ├─ CLCALC → CL, CM, CDP
│  └─ CDCALC → CD, CDF
│
├─ IF NOT LWDIJ OR NOT LADIJ:
│  └─ QDCALC → DIJ(N+NW, N+NW) source influence matrix
│
└─ Newton Iteration Loop (ITER = 1 to NITER):
   ├─ SETBL → Assemble BL Newton system
   │  └─ [see SETBL below]
   │
   ├─ BLSOLV → Solve BL system
   │  └─ Custom block-tridiagonal solver
   │
   ├─ UPDATE → Apply solution to BL variables
   │  └─ UEDG, DSTR, THET, CTAU with underrelaxation RLX
   │
   ├─ IF LALFA:
   │  ├─ MRCL(CL, ...) → Update MINF, REINF from CL
   │  └─ COMSET → Update compressibility params
   │ ELSE:
   │  ├─ QISET → Update QINV for new alpha
   │  └─ UICALC → Update UINV
   │
   ├─ QVFUE → QVIS(1:N+NW) from UEDG + DIJ*SIG
   ├─ GAMQV → GAM(1:N) from QVIS
   ├─ STMOVE → Update SST, IST
   ├─ CLCALC → CL, CM, CDP
   ├─ CDCALC → CD, CDF
   │
   └─ IF RMSBL < EPS1: LVCONV=.TRUE., EXIT
```

#### Variables Updated by VISCAL (per iteration):

| Variable     | Common Block | Description                 |
|--------------|--------------|-----------------------------|
| UEDG(IVX,2)  | CR15         | BL edge velocity            |
| DSTR(IVX,2)  | CR15         | Displacement thickness      |
| THET(IVX,2)  | CR15         | Momentum thickness          |
| CTAU(IVX,2)  | CR15         | Max shear / amplitude       |
| MASS(IVX,2)  | CR15         | Mass defect = UEDG*DSTR     |
| SIG(1:N+NW)  | CR06         | Source panel strengths      |
| QVIS(1:N+NW) | CR04         | Viscous tangential velocity |
| CPV(1:N+NW)  | CR04         | Viscous Cp                  |
| GAM(1:N)     | CR06         | Vorticity from QVIS         |
| SST          | CR06         | Stagnation arc length       |
| IST          | CI04         | Stagnation panel index      |
| CL, CM, CDP  | CR09         | Forces from GAM             |
| CD, CDF      | CR09         | Total and friction drag     |
| RMSBL, RMXBL | CR17         | Convergence metrics         |
| ITRAN(1:2)   | CI05         | Transition indices          |
| XOCTR(1:2)   | CR15         | Transition x/c              |

### Control Flow: SETBL (BL Newton System Assembly)

```
SETBL
├─ Set CLMR = CL or CLSPEC
├─ MRCL(CLMR, ...) → MINF, REINF
├─ COMSET → TKLAM, TKL_MSQ
├─ Set BL parameters: GAMBL, GM1BL, QINFBL, TKBL, etc.
│
├─ Initialize: RMSBL=0, RMXBL=0
│
└─ For IS = 1, 2 (upper/lower surfaces):
   ├─ IBLTE = IBLTE(IS)
   ├─ NBL = NBL(IS)
   │
   ├─ Stagnation point initialization:
   │  └─ THET(1,IS), DSTR(1,IS), CTAU(1,IS) = 0
   │
   └─ BL Marching (IBL = 2 to NBL):
      ├─ Set station indices: I=IPAN(IBL,IS), IM=IPAN(IBL-1,IS)
      │
      ├─ BLKIN(2) → T2, D2, U2, H2, HK2, RT2, M2 from THET, DSTR, UEDG
      │
      ├─ IF WAKE: Use wake closures
      │ ELIF TURB (IBL >= ITRAN):
      │  ├─ BLKIN(1) → Previous station T1, D1, U1, etc.
      │  ├─ HST(HK2, RT2, ...) → Shape parameter HS2
      │  ├─ CFT(HK2, RT2, ...) → Skin friction CF2
      │  ├─ DIT(HS2, US2, ...) → Dissipation DI2
      │  └─ BLVAR(2) → Set auxiliary variables
      │ ELSE LAMINAR:
      │  ├─ BLKIN(1)
      │  ├─ DAMPL(HK1, TH1, ...) → Amplification AX1
      │  ├─ DAMPL(HK2, TH2, ...) → Amplification AX2
      │  ├─ HSL(HK2, RT2, ...) → Shape parameter HS2
      │  ├─ CFL(HK2, RT2, ...) → Skin friction CF2
      │  ├─ DIL(HK2, RT2, ...) → Dissipation DI2
      │  └─ BLVAR(2)
      │
      ├─ BLSYS → Build local system coefficients
      │  └─ VSREZ, VS1, VS2 → Residuals and Jacobians
      │
      ├─ IF LAMINAR and AX >= AMCRIT:
      │  └─ TRCHEK2 → Check/set transition
      │
      ├─ Update global residual:
      │  ├─ RMSBL += (residuals)^2
      │  └─ RMXBL = max(RMXBL, max_residual)
      │
      └─ Build global Newton system blocks:
         ├─ VA(3,2,IBL) → Diagonal block
         ├─ VB(3,2,IBL) → Off-diagonal block
         ├─ VM(3,IBL,JBL) → Mass influence
         └─ VDEL(3,2,IBL) → RHS/residual
```

#### Closure Functions Called by SETBL:

| Function | Input          | Output                       | Purpose                 |
|----------|----------------|------------------------------|-------------------------|
| HKIN     | H, MSQ         | HK, HK_H, HK_MSQ             | Kinematic shape factor  |
| BLKIN    | station vars   | T2, D2, U2, H2, HK2, RT2, M2 | BL input variables      |
| HSL      | HK, RT, MSQ    | HS, HS_HK, HS_RT, HS_MSQ     | Laminar shape param     |
| HST      | HK, RT, MSQ    | HS, HS_HK, HS_RT, HS_MSQ     | Turbulent shape param   |
| CFL      | HK, RT, MSQ    | CF, CF_HK, CF_RT, CF_MSQ     | Laminar skin friction   |
| CFT      | HK, RT, MSQ    | CF, CF_HK, CF_RT, CF_MSQ     | Turbulent skin friction |
| DIL      | HK, RT         | DI, DI_HK, DI_RT             | Laminar dissipation     |
| DIT      | HS, US, CF, ST | DI, DI_*                     | Turbulent dissipation   |
| DAMPL    | HK, TH, RT     | AX, AX_HK, AX_TH, AX_RT      | e^N amplification       |
| BLVAR    | all            | auxiliary BL vars            | Combined BL state       |
| BLSYS    | all            | VSREZ, VS1, VS2              | Newton system coeffs    |
| TRCHEK2  | AX1, AX2, ...  | XT, ITRAN                    | Transition check        |

### Subroutine Call Frequency (from instrumented run)

| Subroutine | Calls  | Per-Iteration Average   |
|------------|--------|-------------------------|
| DAMPL      | 12,110 | ~867 (laminar stations) |
| HKIN       | 12,067 | ~862                    |
| CFL        | 7,916  | ~566 (laminar)          |
| BLKIN      | 7,711  | ~551                    |
| BLVAR      | 4,635  | ~331                    |
| BLSYS      | 4,327  | ~309                    |
| CFT        | 4,041  | ~289 (turbulent)        |
| DIL        | 4,003  | ~286                    |
| HSL        | 3,258  | ~233                    |
| HST        | 2,009  | ~144                    |
| TRCHEK2    | 90     | ~6 (transition checks)  |
| VISCAL     | 2      | 1 per alpha             |

### Key State Flags

| Flag   | Set By       | Purpose                              |
|--------|--------------|--------------------------------------|
| LGAMU  | GGCALC       | .TRUE. when GAMU arrays computed     |
| LQAIJ  | GGCALC       | .TRUE. when AIJ factored             |
| LWAKE  | XYWAKE       | .TRUE. when wake trajectory computed |
| LBLINI | SETBL/MRCHUE | .TRUE. when BL initialized           |
| LIPAN  | IBLPAN       | .TRUE. when BL→panel pointers set    |
| LADIJ  | QDCALC       | .TRUE. when DIJ computed for airfoil |
| LWDIJ  | QDCALC       | .TRUE. when DIJ computed for wake    |
| LVCONV | VISCAL       | .TRUE. when RMSBL < EPS1 (converged) |
| LALFA  | OPER         | .TRUE. if alpha specified (vs CL)    |
| LVISC  | OPER         | .TRUE. if viscous mode enabled       |

### Data Flow Summary

```
NACA/LOAD → XB, YB, NB (buffer airfoil)
    ↓
PANGEN → X, Y, S, N, NX, NY (panels)
    ↓
GGCALC → AIJ, GAMU (influence matrix, unit solutions)
    ↓
SPECAL → GAM, QINV, CPI, CL, CM (inviscid solution)
    ↓
VISCAL:
    XYWAKE → X, Y for wake panels
    QDCALC → DIJ (source influence)
    ↓
    SETBL (iterate):
        UEDG → BLKIN → HKIN, HSL/HST, CFL/CFT, DIL/DIT
                    ↓
               BLVAR → BLSYS
                    ↓
        BLSOLV → VDEL (solution)
                    ↓
        UPDATE → UEDG, DSTR, THET, CTAU
                    ↓
        QVFUE → QVIS = QINV + DIJ * MASS/UEDG
                    ↓
        GAMQV → GAM from QVIS
                    ↓
        CLCALC → CL, CM, CDP
        CDCALC → CD, CDF
    ↓
    Converged: CL, CD, CM, UEDG, DSTR, THET, etc.
```

## XFOIL → yFoil mapping

Every XFOIL name (COMMON variable, subroutine, local of note) and the yFoil name it became, with
what the quantity is. The *rules* behind these names — symbols as names, the underscore
meanings, the four index systems, the derivative tokens, and the interpretation notes on `AC`,
`CTAU`, `US`, `RMSBL` and `GAMMA` — are in
[`docs/conventions/naming.md`](conventions/naming.md); this file is their application and is
kept current whenever a field, variable or function is added or renamed (CLAUDE.md, Architecture).

Indexing is unchanged from XFOIL: 1-based, side 1 = upper and 2 = lower, the wake appended to
side 2 (`NBL(2) = IBLTE(2) + NW`), station arrays `[side][i_station]`. "—" in the XFOIL column
means no Fortran counterpart. XFOIL's own dumps under `tests/fixtures/xfoil/` keep the Fortran
names by design; everything yFoil writes uses the yFoil names.

### Table 1: `BlState` → `SolverState` (`src/solver/blstate.rs`), panel-node section

| XFOIL | yFoil | What it is |
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

### Table 2: `SolverState`, BL-station section (`[side][i_station]`)

| XFOIL | yFoil | What it is |
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

### Table 3: `SolverState`, flow conditions, flags, forces

| XFOIL | yFoil | What it is |
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

### Table 4: `BLStationState` → `StationState` (`/V_VAR2/`, XBL.INC)

The struct is instantiated as `station1` and `station2`, so the frame suffix lives on the
instance (`station2.sqrtctau`), not on the fields. Sensitivities follow `<name>_d_<token>`;
base names are listed with the sensitivities XFOIL carries.

| XFOIL | yFoil | What it is |
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
| — | `mass_defect` | Ue·δ* (yFoil-only cached copy) |

### Table 5: `BLGlobalParams` → `FlowParameters` (`/V_VAR/`) and the closure constants (`/BLPAR/`, set in BLPINI)

| XFOIL | yFoil | What it is |
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

### Table 6: interval-level types (`/V_SYS/`, `/V_VARA/`, `/V_INT/`)

| XFOIL | yFoil | What it is |
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

### Table 7: global Newton system and iteration records

| XFOIL | yFoil | What it is |
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

### Table 8: inviscid system

| XFOIL | yFoil | What it is |
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

### Table 9: session, flow conditions and results (`src/solver/analysis.rs`)

| XFOIL | yFoil | What it is |
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

### Table 10: functions (each carries `#[doc(alias = "XFOIL NAME")]`)

| XFOIL | yFoil | What it does |
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
| PANGEN | `repanel_by_curvature` (`PangenConfig`); `repanel` (`PanelConfig`) is PANE/PPAR with the trailing-edge treatment and the provenance record | |
| GETPAN (PPAR menu) | `PanelConfig` / `PangenConfig` fields, `yfoil geometry repanel` flags | see table 12 |
| TGAP | `set_te_gap` (`gap`, `blend` = XFOIL's DOC) | GDES trailing-edge gap; gated by `tests/xfoil_tgap_tests.rs` |
| — | keep | doc: no XFOIL equivalent |
| NACA4/NACA5 | `naca_4digit`, `naca_5digit` with a `Thickness` {`Perpendicular`, `Vertical`} argument | the un-suffixed name currently holds the non-XFOIL algorithm; every NACA family (4, 4M, 5, 16, 6, 6A) is `series::Section`, which has no XFOIL equivalent |
| SCALC + SEGSPL + LEFIND + TECALC + NCALC + APCALC | `panel_foil` | doc lists all six |
| OPER ALFA (fresh session) | `analyse` | |

### Table 11: geometry and output types, JSON keys

| XFOIL | yFoil | What it is |
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

### Table 12: CLI (no aliases; a short flag only where it is the same letter under every subcommand)

| XFOIL | yFoil | Note |
|---|---|---|
| OPER ALFA | `analyse` | |
| ITER | `--max-iterations` | |
| VISC | keep | |
| VPAR N | `--ncrit` (no short flag) | `-n` is `--panels` under `geometry` |
| PANE, PPAR | `geometry repanel --method pangen` (the default; also on `naca` and `karman-trefftz`) | `--method cosine` is yFoil's own, no XFOIL equivalent |
| NPAN (PPAR `N`) | `-n`, `--panels`; `PanelConfig::n_nodes` | |
| CVPAR (PPAR `P`) | `--curvature-bunching`; `PangenConfig::curvature_bunching` | attraction coefficient 6·CVPAR |
| CTERAT (PPAR `T`) | `--te-curvature-ratio`; `PangenConfig::te_curvature_ratio` | fictitious TE curvature / averaged LE curvature |
| CTRRAT (PPAR `R`) | `--refined-curvature-ratio`; `PangenConfig::refined_curvature_ratio` | fictitious curvature in the windows / LE curvature |
| XSREF1, XSREF2 (PPAR `XT`) | `--refine-upper X1,X2`; `PangenConfig::refine_upper` | `None` = XFOIL's 1.0 1.0 |
| XPREF1, XPREF2 (PPAR `XB`) | `--refine-lower X1,X2`; `PangenConfig::refine_lower` | |
| TGAP gap doc | `--te-gap`, `--te-blend`; `PanelConfig::te_gap` {`gap`, `blend`} | applied after panelling (XFOIL: to the buffer before PANE) |
| — | `--cosine-te-bias`; `CosineConfig::te_bias` | yFoil's cosine method only |
| — | `--panelling FILE` | the `generator.panelling` record as input |
| NACA4/NACA5 vertical thickness | `--thickness vertical` | XFOIL's model, always PANGEN-panelled; `perpendicular` (default) is the NACA definition |
| DUMP | keep | |
| — | `--include-lagged-closures` | adds `lagged_closures` to the boundary-layer output |

### Table 13: tests and xtask

| XFOIL | yFoil | Note |
|---|---|---|
| ACRIT | `ncrit` | |
| — | `FixtureStation`, `FixtureSide`; `delta_star` → `dstar` | |
| THET_TE1 … | `theta_te_station1` … | the tracked JSON fixture inputs are translated to the new keys (step 9) |
| REINF1, ITMAX | `re`, `max_iterations` | `xtask/fixtures-config/cases.toml` translated in the same step |

---

### JSON ↔ variable: keys that are not a one-to-one print of a variable

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

## spline.f - Spline Utilities

This module contains cubic spline routines used for smooth interpolation of geometry and flow variables.

### Key Subroutines

#### SPLINE (Line ~1)

**Purpose:** Computes cubic spline coefficients for a set of data points with natural (zero second derivative) end conditions.

**Input:**
| Variable | Description |
|----------|-------------|
| `X(N)` | Independent variable array |
| `S(N)` | Dependent variable array (arc length) |
| `N` | Number of data points |

**Output:**
| Variable | Description |
|----------|-------------|
| `XS(N)` | Spline coefficients dX/dS |

**Algorithm:**

Natural cubic spline with:
```
X''(S_1) = 0  (start)
X''(S_N) = 0  (end)
```

The spline interpolant between points i and i+1:
```
X(s) = a_i + b_i*(s-s_i) + c_i*(s-s_i)² + d_i*(s-s_i)³
```

Solved via tridiagonal system.

---

#### SPLIND (Line ~50)

**Purpose:** Computes cubic spline coefficients with specified end slopes.

**Input:**
| Variable | Description |
|----------|-------------|
| `X(N)` | Independent variable array |
| `S(N)` | Arc length array |
| `N` | Number of points |
| `XS1` | Slope dX/dS at start |
| `XS2` | Slope dX/dS at end |

**Output:**
| Variable | Description |
|----------|-------------|
| `XS(N)` | Spline coefficients |

**Algorithm:**

Clamped cubic spline with:
```
X'(S_1) = XS1  (start slope)
X'(S_N) = XS2  (end slope)
```

---

#### SEVAL (Line ~100)

**Purpose:** Evaluates spline value at a given parameter.

**Input:**
| Variable | Description |
|----------|-------------|
| `SS` | Parameter value |
| `X(N)` | Data array |
| `XS(N)` | Spline coefficients |
| `S(N)` | Parameter array |
| `N` | Number of points |

**Output:**
| Return value | Description |
|--------------|-------------|
| `SEVAL` | Interpolated X value |

**Algorithm:**

1. Binary search to find interval [i, i+1] containing SS
2. Evaluate cubic polynomial:
```
t = SS - S(i)
h = S(i+1) - S(i)

SEVAL = X(i) + XS(i)*t
      + [3*(X(i+1)-X(i))/h - 2*XS(i) - XS(i+1)] * t²/h
      + [2*(X(i)-X(i+1))/h + XS(i) + XS(i+1)] * t³/h²
```

---

#### DEVAL (Line ~150)

**Purpose:** Evaluates spline derivative at a given parameter.

**Input:**
| Variable | Description |
|----------|-------------|
| `SS` | Parameter value |
| `X(N)` | Data array |
| `XS(N)` | Spline coefficients |
| `S(N)` | Parameter array |
| `N` | Number of points |

**Output:**
| Return value | Description |
|--------------|-------------|
| `DEVAL` | Interpolated dX/dS value |

**Algorithm:**

Derivative of cubic polynomial:
```
DEVAL = XS(i)
      + 2*[3*(X(i+1)-X(i))/h - 2*XS(i) - XS(i+1)] * t/h
      + 3*[2*(X(i)-X(i+1))/h + XS(i) + XS(i+1)] * t²/h²
```

---

#### D2VAL (Line ~200)

**Purpose:** Evaluates spline second derivative at a given parameter.

**Output:**
| Return value | Description |
|--------------|-------------|
| `D2VAL` | Interpolated d²X/dS² value |

---

### Usage in XFOIL

#### Geometry Splines

Airfoil coordinates are stored with splines:
```fortran
CALL SPLINE(X, XP, S, N)     ! X coordinates
CALL SPLINE(Y, YP, S, N)     ! Y coordinates
```

To evaluate at arc length SS:
```fortran
XVAL = SEVAL(SS, X, XP, S, N)
YVAL = SEVAL(SS, Y, YP, S, N)
```

#### Curvature

Local curvature κ is:
```
κ = (X'*Y'' - Y'*X'') / (X'² + Y'²)^(3/2)
```

where primes denote derivatives with respect to arc length S.

#### Panel Normal

Normal vector at point I:
```
NX(I) = YP(I) / sqrt(XP(I)² + YP(I)²)
NY(I) = -XP(I) / sqrt(XP(I)² + YP(I)²)
```

---

### yFoil Equivalent

`src/geometry/spline.rs`

```rust
/// Cubic spline with natural end conditions
pub struct CubicSpline {
    s: Vec<f64>,      // Parameter values
    x: Vec<f64>,      // Data values
    xs: Vec<f64>,     // Spline coefficients (derivatives)
}

impl CubicSpline {
    /// Create natural cubic spline
    pub fn natural(s: &[f64], x: &[f64]) -> Self {
        let n = s.len();
        let mut xs = vec![0.0; n];

        // Set up tridiagonal system
        // [b_i  c_i     ] [xs_1  ]   [d_1]
        // [a_i  b_i  c_i] [xs_2  ] = [d_2]
        // [     a_i  b_i] [xs_n  ]   [d_n]

        // Natural BC: xs'' = 0 at ends

        // Thomas algorithm for tridiagonal solve
        // ...

        Self { s: s.to_vec(), x: x.to_vec(), xs }
    }

    /// Evaluate spline at parameter value
    pub fn eval(&self, ss: f64) -> f64 {
        // Binary search for interval
        let i = self.find_interval(ss);

        // Cubic interpolation
        let t = ss - self.s[i];
        let h = self.s[i + 1] - self.s[i];

        let a = self.x[i];
        let b = self.xs[i];
        let c = (3.0 * (self.x[i + 1] - self.x[i]) / h
                - 2.0 * self.xs[i] - self.xs[i + 1]) / h;
        let d = (2.0 * (self.x[i] - self.x[i + 1]) / h
                + self.xs[i] + self.xs[i + 1]) / (h * h);

        a + b * t + c * t * t + d * t * t * t
    }

    /// Evaluate derivative at parameter value
    pub fn deval(&self, ss: f64) -> f64 {
        let i = self.find_interval(ss);
        let t = ss - self.s[i];
        let h = self.s[i + 1] - self.s[i];

        let b = self.xs[i];
        let c = (3.0 * (self.x[i + 1] - self.x[i]) / h
                - 2.0 * self.xs[i] - self.xs[i + 1]) / h;
        let d = (2.0 * (self.x[i] - self.x[i + 1]) / h
                + self.xs[i] + self.xs[i + 1]) / (h * h);

        b + 2.0 * c * t + 3.0 * d * t * t
    }
}
```

## xpanel.f - Panel Method

This module contains the panel method routines for computing inviscid flow via vortex distribution.

### Key Subroutines

#### PSILIN (Line ~99)

**Purpose:** Calculates streamfunction ψ at a panel node due to freestream and all bound vorticity. Computes sensitivities for the Newton system.

**Called by:** `QDCALC`, inviscid solver setup

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `I` | Panel node index |
| `XI, YI` | Node coordinates |
| `NXI, NYI` | Node normal vector |
| `GAM(IQX)` | Vortex distribution |
| `SIG(IZX)` | Source distribution (if SIGLIN) |
| `ALFA` | Angle of attack |
| `QINF` | Freestream speed |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `PSI` | Streamfunction at node |
| `PSI_NI` | dψ/dn (normal derivative) |
| `DZDG(IQX)` | dψ/dγ for each panel |
| `DZDM(IZX)` | dψ/dσ if SIGLIN=true |
| `DZDN(IQX)` | dψ/dn if GEOLIN=true |
| `Z_ALFA` | dψ/dα |
| `Z_QINF` | dψ/dQ∞ |
| `Z_QDOF0-3` | Inverse design sensitivities |

**Algorithm:**

1. **Freestream contribution**:
```
ψ_∞ = -Q∞ * (X*sinα - Y*cosα)
dψ/dα = -Q∞ * (X*cosα + Y*sinα)
```

2. **Panel loop** (J = 1 to N):

   For each panel J with endpoints (X_J, Y_J) and (X_{J+1}, Y_{J+1}):

   a. **Compute panel geometry**:
   ```
   DSJ = sqrt((X_{J+1}-X_J)² + (Y_{J+1}-Y_J)²)  [panel length]
   SX = (X_{J+1}-X_J) / DSJ                      [tangent x]
   SY = (Y_{J+1}-Y_J) / DSJ                      [tangent y]
   ```

   b. **Transform to local coordinates**:
   ```
   X1 = (XI-X_J)*SX + (YI-Y_J)*SY               [start distance]
   X2 = (XI-X_{J+1})*SX + (YI-Y_{J+1})*SY       [end distance]
   YY = (XI-X_J)*(-SY) + (YI-Y_J)*SX            [normal distance]
   ```

   c. **Evaluate influence integrals**:
   ```
   R1² = X1² + YY²
   R2² = X2² + YY²

   If YY ≠ 0:
       Θ1 = atan2(YY, X1)
       Θ2 = atan2(YY, X2)
   Else:
       Θ1, Θ2 = π or 0 depending on sign of X

   G = 0.5 * ln(R1/R2)                          [log term]
   T = Θ1 - Θ2                                  [angle term]
   ```

   d. **Vortex panel influence** (linear distribution):
   ```
   ψ_J = (γ_J + γ_{J+1})/2 * YY*G / (2π)
       + (γ_{J+1} - γ_J) * (X2*T - X1*G - YY*ln(R2/R1)) / (2π*DSJ)
   ```

3. **Accumulate sensitivities**:
```
DZDG(J) = ∂ψ/∂γ_J
DZDN(I) = ∂ψ/∂n_I     (geometric sensitivity)
```

**Mathematical Background:**

The streamfunction due to a point vortex at (ξ, η):
```
ψ = γ/(2π) * ln(r)
```

where r = √((x-ξ)² + (y-η)²).

For a linear vortex panel, integrate analytically:
```
ψ = ∫ γ(s)/(2π) * ln(r(s)) ds
```

with γ(s) = γ_J + (γ_{J+1} - γ_J) * s/DSJ.

---

#### QDCALC (Line ~1161)

**Purpose:** Computes the source panel influence coefficient matrix DIJ, which gives the sensitivity of tangential velocity to source strength.

**Called by:** `VISCAL`

**Calls:**
- `PSILIN` - For influence coefficients
- `BAKSUB` - Backsubstitution with factored AIJ

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `AIJ(IQX,IQX)` | Factored vortex influence matrix |
| `AIJPIV(IQX)` | Pivot indices |
| `N` | Number of airfoil panels |
| `NW` | Number of wake panels |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `DIJ(IZX,IZX)` | dQ_tan/dσ matrix |
| `CIJ(IWX,IQX)` | dQ_tan/dγ at wake |
| `LADIJ` | TRUE when airfoil DIJ computed |
| `LWDIJ` | TRUE when wake DIJ computed |

**Algorithm:**

1. **For each source location J**:

   a. **Compute direct influence on streamfunction**:
   ```
   Call PSILIN(I, ..., SIGLIN=true) for each node I
   → DZDM(I) = dψ_I/dσ_J
   ```

   b. **Solve for induced vorticity change**:
   ```
   Kutta condition: ψ = constant on airfoil

   [AIJ] * [Δγ] = -[DZDM]

   Use BAKSUB with pre-factored AIJ to get Δγ/Δσ_J
   ```

   c. **Compute tangential velocity sensitivity**:
   ```
   DIJ(I,J) = direct source effect + Δγ induced effect

   Direct: ∂Q_tan/∂σ from source integral
   Induced: Σ_k (∂Q_tan/∂γ_k) * (∂γ_k/∂σ_J)
   ```

2. **For wake points**:

   The wake sees both direct source effect and indirect effect through airfoil circulation changes:
   ```
   DIJ_wake(I,J) = DIJ_direct(I,J) + CIJ(I,:) * DIJ_airfoil(:,J)
   ```

**Matrix Structure:**

The full DIJ matrix:
```
        Source at:
        1  2  ... N  N+1 ... N+NW
      ┌────────────────────────────┐
    1 │                            │
    2 │     Airfoil-Airfoil        │  LADIJ
  ... │        (N × N)             │
    N │                            │
  ────│────────────────────────────│
  N+1 │                            │
  ... │   Wake-Airfoil + Wake-Wake │  LWDIJ
N+NW │        (NW × N+NW)         │
      └────────────────────────────┘
```

**Usage in VISCAL:**

The DIJ matrix enables efficient coupling:
```
Q_vis(I) = Q_inv(I) + Σ_J DIJ(I,J) * σ_J
```

where σ_J is the mass defect source:
```
σ_J = d(ρ*Ue*δ*)/ds
```

---

#### QVFUE (Line ~1400)

**Purpose:** Calculates edge velocities QVIS from UEDG via the mass defect influence.

**Called by:** `VISCAL`

**Algorithm:**

```
1. Compute source distribution from mass defect:
   σ_I = (MASS_I - MASS_{I-1}) / (S_I - S_{I-1})

2. Apply influence matrix:
   QVIS_I = QINV_I + Σ_J DIJ(I,J) * σ_J
```

---

#### UICALC (Line ~1300)

**Purpose:** Sets inviscid BL edge velocity UINV from panel velocity QINV.

**Algorithm:**

```
For each BL station:
   UINV(IBL,IS) = |QINV(IPAN(IBL,IS))| * VTI(IBL,IS)
```

where VTI = ±1 accounts for direction convention.

---

### Influence Coefficient Theory

#### Vortex Panel Method

The airfoil surface is discretized into N panels. Each panel has:
- Linear vortex distribution γ(s)
- Endpoints shared with neighbors

The influence of panel J on point I:
```
ψ_IJ = ∫_panel_J γ(s)/(2π) * ln|r_I - r(s)| ds
```

For linear γ(s), this integrates to:
```
ψ_IJ = (γ_J * F1 + γ_{J+1} * F2) / (2π)
```

where F1, F2 are geometric functions involving ln(r) and arctan terms.

#### Source Panels for Viscous Coupling

Mass defect creates equivalent sources:
```
σ = d(m)/ds = d(ρ*Ue*δ*)/ds
```

The source influence on velocity:
```
ΔQ_tan = ∫ σ(s')/(2π) * sin(θ)/r ds'
```

#### Newton System Coupling

The global system couples:
1. Panel method: [AIJ][γ] = [RHS] (Kutta condition)
2. BL equations: f(θ, δ*, Ue) = 0
3. Mass defect: σ = g(Ue, δ*)
4. Velocity update: Ue_new = Ue_inv + DIJ * σ

---

### yFoil Equivalent

`src/solver/psilin.rs`, `src/solver/ggcalc.rs`, `src/solver/qdcalc.rs`, `src/solver/xywake.rs`, `src/solver/pointers.rs`, `src/solver/velocity.rs`

```rust
/// Calculate streamfunction influence coefficients
pub fn psilin(
    node: usize,
    geometry: &Geometry,
    gamma: &[f64],
    sigma: Option<&[f64]>,
) -> PsilinResult {
    let (xi, yi) = geometry.node(node);
    let (nxi, nyi) = geometry.normal(node);

    let mut psi = 0.0;
    let mut dzdg = vec![0.0; geometry.n_panels()];

    // Freestream contribution
    psi -= geometry.qinf * (xi * geometry.sina - yi * geometry.cosa);

    // Panel contributions
    for j in 0..geometry.n_panels() {
        let (x1, y1) = geometry.panel_start(j);
        let (x2, y2) = geometry.panel_end(j);

        // Local coordinates
        let dsj = geometry.panel_length(j);
        let sx = (x2 - x1) / dsj;
        let sy = (y2 - y1) / dsj;

        let xl1 = (xi - x1) * sx + (yi - y1) * sy;
        let xl2 = (xi - x2) * sx + (yi - y2) * sy;
        let yy = -(xi - x1) * sy + (yi - y1) * sx;

        // Influence integrals
        let (g, t) = panel_integrals(xl1, xl2, yy);

        // Accumulate
        let psi_j = vortex_influence(gamma[j], gamma[j + 1], yy, g, t, dsj);
        psi += psi_j;
        dzdg[j] += dpsi_dgamma_j(...);
    }

    PsilinResult { psi, dzdg, ... }
}

/// Compute source influence matrix
pub fn qdcalc(
    geometry: &Geometry,
    aij: &FactoredMatrix,
) -> SourceInfluenceMatrix {
    let mut dij = Matrix::zeros(geometry.n_total(), geometry.n_total());

    for j in 0..geometry.n_total() {
        // Direct influence
        let dzdm = source_influence(j, geometry);

        // Induced vorticity
        let dgamma = aij.solve(&(-dzdm));

        // Tangential velocity sensitivity
        for i in 0..geometry.n_total() {
            dij[(i, j)] = direct_qtan_sigma(i, j)
                + cij_row(i).dot(&dgamma);
        }
    }

    SourceInfluenceMatrix { dij }
}
```

## xbl.f - Boundary Layer Marching

This module contains the BL marching routines that solve the integral boundary layer equations station by station.

### Key Subroutines

#### SETBL (Line ~21)

**Purpose:** Sets up the BL Newton system coefficients for current BL variables and incorporates them into the global Newton system.

**Called by:** `VISCAL`

**Calls:**
- `COMSET` - Set compressibility parameters
- `MRCHUE` - March BL in direct mode (initialization)
- `MRCHDU` - March BL in direct/inverse mode
- `BLPRV` - Set primary BL variables
- `BLKIN` - Compute secondary variables
- `TRCHEK` - Check for transition

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `UEDG(IVX,ISX)` | Edge velocity array |
| `THET(IVX,ISX)` | Momentum thickness |
| `DSTR(IVX,ISX)` | Displacement thickness |
| `CTAU(IVX,ISX)` | Shear stress coefficient |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `VA, VB` | Newton system diagonal/off-diagonal blocks |
| `VDEL` | Residual vectors |
| `VM` | Mass influence vectors |
| `TAU, DIS` | Wall shear and dissipation (for output) |

**Algorithm:**

```
1. Set compressibility parameters via COMSET
   ├── TKLAM = Karman-Tsien parameter
   └── TK_MSQ = d(TKLAM)/d(M²)

2. Initialize BL if needed (LBLINI = false)
   └── Call MRCHUE for direct-mode march

3. March BL with current Ue (MRCHDU)
   └── Mixed direct/inverse mode

4. For each BL station on each side:
   │
   ├── Set primary variables (BLPRV)
   │   └── X2, U2, T2, D2, S2
   │
   ├── Compute secondary variables (BLKIN)
   │   └── M2, R2, HK2, RT2, etc.
   │
   ├── Check for transition (TRCHEK)
   │   └── Set ITRAN if transition occurs
   │
   └── Assemble 10×4 linearized system
       ├── VS1(4,5) - Station 1 Jacobian
       ├── VS2(4,5) - Station 2 Jacobian
       └── VSREZ(4) - Residuals

5. Store diagnostic quantities
   ├── TAU = wall shear stress
   └── DIS = dissipation
```

---

#### MRCHUE (Line ~542)

**Purpose:** Marches boundary layers and wake in direct mode using prescribed UEDG array. Handles separation via kinematic shape parameter extrapolation.

**Called by:** `SETBL` (for BL initialization)

**Calls:**
- `BLPRV` - Set primary variables
- `BLKIN` - Compute secondary variables
- `BLVAR` - Compute closure relations
- `TRCHEK` - Check transition
- `GAUSS` - Solve 4×4 system

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `UEDG(IVX,ISX)` | Prescribed edge velocity |
| `XSSI(IVX,ISX)` | Arc length coordinates |
| `AMCRIT` | Critical amplification factor |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `THET(IVX,ISX)` | Momentum thickness |
| `DSTR(IVX,ISX)` | Displacement thickness |
| `CTAU(IVX,ISX)` | Shear/amplification |
| `ITRAN(ISX)` | Transition location index |

**Algorithm:**

```
For each side (IS = 1, 2):
│
├── Initialize similarity station (IBL = 2)
│   └── θ² ≈ 0.45ν / (U·(5β+1)·Re)  [Thwaites]
│
└── March downstream (IBL = 3 to NBL):
    │
    ├── Set station 1 from previous station
    │
    ├── Newton iteration (up to 25 iterations):
    │   │
    │   ├── Set primary variables (BLPRV)
    │   │
    │   ├── Compute secondary variables (BLKIN)
    │   │
    │   ├── Check transition (TRCHEK)
    │   │   └── If N > Ncrit, set transition
    │   │
    │   ├── Assemble 4×4 direct-mode system
    │   │   └── dUe = 0 (prescribed)
    │   │
    │   ├── Solve via Gaussian elimination
    │   │
    │   └── Check Hk for separation
    │       └── If Hk > Hmax, switch to inverse mode
    │
    └── Handle separation/failure
        └── Extrapolate using power law
```

**Direct Mode System:**

In direct mode, edge velocity is prescribed:
```
dUe = 0
```

The 3×3 system solves for:
- `dθ` - momentum thickness change
- `dδ*` - displacement thickness change
- `dA` - amplification (laminar) or Ctau (turbulent)

---

#### MRCHDU (Line ~886)

**Purpose:** Marches BL in mixed direct/inverse mode tracking the Ue-Hk characteristic to avoid the Goldstein singularity at separation.

**Called by:** `SETBL` (for coupled iteration)

**Calls:**
- `BLPRV`, `BLKIN`, `BLVAR` - BL variable computation
- `TRCHEK` - Transition detection

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `UEDG(IVX,ISX)` | Current edge velocity |
| `THET, DSTR, CTAU` | Current BL state |
| `SENSWT = 1000` | Sensitivity weight |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `THET, DSTR, CTAU` | Updated BL state |
| `UEDG` | Updated edge velocity |

**Algorithm:**

```
For each side (IS = 1, 2):
│
└── March downstream (IBL = 2 to NBL):
    │
    ├── Set baseline Ue, Hk references
    │
    ├── For similarity station or first wake point:
    │   └── Prescribe Ue (direct mode)
    │
    └── For interior stations:
        │
        ├── Solve dHk system to find dUe/dHk slope
        │
        ├── Construct prescribed Ue-Hk line
        │   └── Quasi-normal to characteristic
        │
        └── Solve mixed system with weighted constraint
            └── SENSWT * (Ue - Ue_target) + (Hk - Hk_target) = 0
```

**Mixed Mode Formulation:**

The characteristic relation near separation becomes ill-conditioned:
```
dHk/dUe → 0  at separation
```

MRCHDU avoids this by:
1. Computing the local Ue-Hk characteristic slope
2. Imposing a constraint quasi-normal to the characteristic
3. Weighting between Ue and Hk prescription

---

#### UPDATE (Line ~1264)

**Purpose:** Adds Newton deltas to BL variables with underrelaxation to prevent excessive changes.

**Called by:** `VISCAL`

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `VDEL(3,2,IZX)` | Newton solution vector |
| `DALMAX = 0.5°` | Max alpha change |
| `DCLMAX = 0.5` | Max CL change |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `ALFA` or `CL` | Updated control variable |
| `CTAU, THET, DSTR` | Updated BL arrays |
| `UEDG, MASS` | Updated mass defect |
| `RMSBL, RMXBL` | RMS and max changes |
| `RLX` | Relaxation factor used |

**Algorithm:**

```
1. Calculate new Ue from mass defect updates
   └── QNEW = QOLD + DIJ · ΔSIG

2. Recompute CL from updated velocity
   └── CL_NEW from integrated pressure

3. Set ΔAC (alpha or CL change) from Newton solution

4. Apply underrelaxation
   │
   ├── If |ΔAC| > limit:
   │   └── RLX = limit / |ΔAC|
   │
   ├── If any δ* change > DMAX:
   │   └── RLX = min(RLX, DMAX / |Δδ*|)
   │
   └── Apply: Δ → RLX · Δ

5. Update BL variables
   ├── CTAU = CTAU + RLX · ΔCTAU
   ├── THET = THET + RLX · ΔTHET
   ├── DSTR = DSTR + RLX · ΔDSTR
   └── UEDG = UEDG + RLX · ΔUEDG

6. Compute mass defect
   └── MASS = UEDG · DSTR

7. Enforce physical bounds
   ├── Hk > 1.0
   └── Ctau > 0

8. Sync upper/lower wake arrays
```

**Underrelaxation:**

```
RLX = min(1.0,
          DALMAX / |ΔAlfa|,
          DCLMAX / |ΔCL|,
          DMAX / max(|Δδ*|))
```

---

### yFoil Equivalent

`src/bl/mrchue.rs`, `src/bl/mrchdu.rs`, `src/solver/setbl.rs`, `src/solver/update.rs`

```rust
pub fn march_bl_direct(
    state: &mut BLState,
    uedg: &[f64],
    params: &BLParams,
) -> Result<(), BLError> {
    // March from stagnation point outward
    for ibl in 2..state.nbl {
        // Initialize station
        let station1 = state.station(ibl - 1);

        // Newton iteration
        for _ in 0..25 {
            // Compute secondary variables
            let vars = blkin(&station1, &station2, params);

            // Check transition
            if vars.ampl > params.ncrit {
                state.set_transition(ibl);
            }

            // Solve local system
            let delta = solve_direct_system(&vars);

            // Update and check convergence
            station2.apply_delta(&delta);
            if delta.norm() < 1e-6 {
                break;
            }
        }

        state.set_station(ibl, station2);
    }

    Ok(())
}
```

## xblsys.f - Boundary Layer Newton System

This module contains the core BL closure relations and the Newton system assembly routines.

### Key Subroutines

#### BLPRV (Line ~701)

**Purpose:** Sets BL primary "2" variables from parameter list with compressibility correction.

**Called by:** `MRCHUE`, `MRCHDU`, `SETBL`

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `XSI` | Arc length position |
| `AMI` | Amplification (laminar) or Ctau (turbulent) |
| `THI` | Momentum thickness θ |
| `DSI` | Displacement thickness δ* |
| `DSWAKI` | Wake gap thickness |
| `UEI` | Incompressible edge velocity |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `X2` | Arc length |
| `AMPL2` | Amplification/Ctau |
| `T2` | θ |
| `D2` | δ* |
| `DW2` | Wake gap |
| `U2` | Compressible edge velocity |
| `U2_UEI` | dU2/dUei |
| `U2_MS` | dU2/d(M²) |

**Algorithm:**

Karman-Tsien compressibility correction:
```
U2 = UEI * (1 - TKLAM) / (1 - TKLAM * (UEI/Qinf)²)
```

Derivatives:
```
U2_UEI = (1 - TKLAM) / (1 - TKLAM * (UEI/Qinf)²)²
       + 2 * UEI² * TKLAM * (1 - TKLAM) / (Qinf² * (1 - TKLAM * (UEI/Qinf)²)²)

U2_MS = ... (similar chain rule)
```

---

#### BLKIN (Line ~725)

**Purpose:** Calculates turbulence-independent secondary "2" variables and their sensitivities.

**Called by:** `MRCHUE`, `MRCHDU`, `SETBL`

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `X2, T2, D2, U2` | Primary variables from BLPRV |
| `HSTINV` | Stagnation temperature ratio |
| `GAMBL` | Cp/Cv |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `M2` | Local Mach number |
| `M2_U2, M2_MS` | Mach derivatives |
| `R2` | Density ratio ρ_e/ρ_∞ |
| `R2_U2, R2_MS` | Density derivatives |
| `H2` | Shape factor H = δ*/θ |
| `H2_T2, H2_D2` | H derivatives |
| `HK2` | Kinematic shape factor |
| `HK2_U2, HK2_T2, HK2_D2, HK2_MS` | Hk derivatives |
| `RT2` | Reynolds number Re_θ |
| `RT2_U2, RT2_T2, RT2_MS, RT2_RE` | Re_θ derivatives |

**Algorithm:**

1. **Mach number** (isentropic relation):
```
M² = U² * (γ-1) * (M∞/Q∞)² / (1 - 0.5 * U² * (γ-1) * (M∞/Q∞)²)
```

2. **Density ratio** (isentropic):
```
TR = 1 - 0.5 * (γ-1) * M²      [Temperature ratio]
R = TR^(1/(γ-1))               [Density ratio]
```

3. **Shape factor**:
```
H = δ* / θ
```

4. **Kinematic shape factor** (via HKIN subroutine):
```
Hk = (H - 0.29*M²) / (1 + 0.113*M² + 0.0056*M⁴)    [Whitfield-Jameson]
```

5. **Reynolds number**:
```
V = μ/μ∞ = TR^1.5 * (1 + S∞)/(TR + S∞)  [Sutherland]
Re_θ = Re∞ * U * θ / V
```

---

#### BLVAR (Line ~784)

**Purpose:** Calculates all secondary "2" variables and sensitivities for laminar, turbulent, and wake regimes.

**Called by:** `MRCHUE`, `MRCHDU`, `BLSYS`

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `ITYP` | 1=laminar, 2=turbulent, 3=wake |
| `X2, U2, T2, D2, S2` | Primary variables |
| `HK2, RT2` | From BLKIN |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `HS2` | Energy shape factor H* = θ*/θ |
| `HS2_* ` | H* derivatives |
| `HC2` | Density shape factor |
| `HC2_*` | Hc derivatives |
| `US2` | Slip velocity |
| `US2_*` | Us derivatives |
| `CF2` | Skin friction Cf/2 |
| `CF2_*` | Cf derivatives |
| `DI2` | Dissipation 2CD/H* |
| `DI2_*` | DI derivatives |
| `DE2` | Density thickness δ** |
| `DE2_*` | δ** derivatives |

**Algorithm:**

1. **Enforce bounds**:
```
Hk ≥ 1.05  (airfoil)
Hk ≥ 1.00005  (wake)
```

2. **Energy shape factor H*** (via HSL/HST):

   Laminar (Falkner-Skan family):
   ```
   H* = f(Hk)  [tabulated correlation]
   ```

   Turbulent (Cebeci-Smith type):
   ```
   H* = 1.505 + 4/(Hk-1) + ...  [complex correlation]
   ```

3. **Density shape factor Hc** (via HCT):
```
Hc = (δ - δ*) / θ = f(Hk, M²)
```

4. **Slip velocity Us**:
```
Us = 0.5 * H* * (1 - (Hk-1)/(Gb*H))
```

5. **Skin friction Cf** (via CFL/CFT):

   Laminar (Blasius/Thwaites):
   ```
   Cf = f(Hk, Re_θ)
   ```

   Turbulent (Ludwieg-Tillmann type):
   ```
   Cf = 0.3 * exp(-1.33*Hk) / (log10(Re_θ))^(1.74+0.31*Hk) + ...
   ```

6. **Dissipation DI** (via DIL/DIT):
```
2*CD/H* = Cf*Us + 2*(Hk*Us/H*)²/(0.5 - a*(...))**2 + ...
```

---

#### BLSYS (Line ~583)

**Purpose:** Assembles the linearized 4-equation finite-difference system governing BL evolution over one interval.

**Called by:** `MRCHUE`, `MRCHDU` (within iteration)

**Input Variables:**
| Variable | Description |
|----------|-------------|
| Station 1, 2 variables | From BLPRV, BLKIN, BLVAR |
| `SIMI` | TRUE if similarity station |
| `WAKE` | TRUE if in wake |
| `TRAN` | TRUE if transition interval |
| `TURB` | TRUE if turbulent |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `VS1(4,5)` | Station 1 Jacobian |
| `VS2(4,5)` | Station 2 Jacobian |
| `VSREZ(4)` | Residual vector |
| `VSM(4)` | Mach sensitivity |
| `VSR(4)` | Reynolds sensitivity |
| `VSX(4)` | ξ sensitivity |

**System Structure:**

The 4 equations for BL marching:

1. **Momentum integral equation**:
```
dθ/dξ + (2 + H - M²) * θ/Ue * dUe/dξ = Cf/2
```

2. **Shape parameter equation** (kinetic energy):
```
d(θ·H*)/dξ + (1 + H* - M²) * θ·H*/Ue * dUe/dξ = 2CD
```

3. **Shear stress lag equation** (turbulent):
```
dCtau/dξ = (Ctau_eq - Ctau) * Ue / (SCCON * δ)
```

   Or **Amplification equation** (laminar):
```
dN/dξ = dn/dRe_θ * dRe_θ/dξ
```

4. **Closure constraint** (varies by mode):
   - Direct mode: `dUe = 0`
   - Inverse mode: `dHk = 0`
   - Mixed mode: `SENSWT * dUe + dHk = 0`

**Finite Difference Form:**

```
[VS1(i,j)] * [Δθ₁, Δδ*₁, ΔA₁, ΔUe₁]ᵀ +
[VS2(i,j)] * [Δθ₂, Δδ*₂, ΔA₂, ΔUe₂]ᵀ = [VSREZ(i)]
```

where:
- j=1: θ column
- j=2: δ* column
- j=3: A (amplification or Ctau) column
- j=4: Ue column
- j=5: RHS (residual)

---

#### TRCHEK2 (Line ~231)

**Purpose:** Checks if natural transition occurs in current X1→X2 interval using second-order implicit amplification equation.

**Called by:** `MRCHUE`, `MRCHDU`, `SETBL`

**Input Variables:**
| Variable | Description |
|----------|-------------|
| `X1, X2` | Interval bounds |
| `AMPL1` | Amplification at station 1 |
| `HK1, HK2` | Kinematic shape factors |
| `RT1, RT2` | Reynolds numbers Re_θ |
| `AMCRIT` | Critical amplification N_crit |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `AMPL2` | Amplification at station 2 |
| `XT` | Transition location (0-1 in interval) |
| `XT_A1, XT_MS, XT_RE` | XT sensitivities |
| `TRAN` | TRUE if transition in interval |

**Algorithm:**

1. **Calculate amplification rate** (via AXSET):
```
dN/dξ = f(Hk, Re_θ)

AX = 0.5 * (dN/dξ|₁ + dN/dξ|₂)  [trapezoidal average]
```

2. **Implicit amplification equation**:
```
N₂ = N₁ + AX * (X₂ - X₁)
```

3. **Check for transition**:
```
If N₂ > N_crit:
    Transition in interval
    XT = (N_crit - N₁) / (N₂ - N₁)
    N₂ = N_crit
Else:
    No transition
    XT = 1.0
```

4. **Newton iteration** for implicit N₂ (up to 30 iterations):
```
Residual = N₂ - N₁ - 0.5*(dN/dξ(N₂) + dN/dξ(N₁))*(X₂-X₁)
```

5. **Calculate sensitivities** via chain rule:
```
XT_A1 = -1 / (N₂ - N₁)
XT_H1 = -XT * dN/dξ_H1 * (X₂-X₁) / (N₂-N₁)
...
```

**Amplification Correlation (DAMPL/AXSET):**

The e^N method uses:
```
dN/dRe_θ = 0.028*(Hk - 1) - 0.0345*exp(-(3.87/(Hk - 1) - 2.52)²)
```

with Re_θ threshold:
```
Re_θ_crit = 10^(2.492/(Hk - 1)^0.43 + 0.7*(tanh(14*(Hk-1)-9.24)+1))
```

---

### Helper Subroutines

#### HKIN (Kinematic Shape Factor)

```
Hk = (H - 0.29*M²) / (1 + 0.113*M² + 0.0056*M⁴)
```

#### DAMPL (Amplification Rate)

The Drela-Giles envelope e^N correlation.

#### CFL/CFT (Skin Friction)

Laminar and turbulent Cf correlations.

#### HSL/HST (Energy Shape Factor)

H* correlations for laminar and turbulent flow.

---

### yFoil Equivalent

`src/bl/closure.rs`, `src/bl/station.rs`, `src/bl/transition.rs`, `src/bl/difference.rs`, `src/bl/blsys.rs`

```rust
/// Kinetic energy shape factor correlation
pub fn hstar(hk: f64, regime: BLRegime) -> (f64, f64) {
    match regime {
        BLRegime::Laminar => hstar_laminar(hk),
        BLRegime::Turbulent => hstar_turbulent(hk),
    }
}

/// Skin friction correlation
pub fn skin_friction(hk: f64, re_theta: f64, regime: BLRegime) -> CfResult {
    match regime {
        BLRegime::Laminar => cf_laminar(hk, re_theta),
        BLRegime::Turbulent => cf_turbulent(hk, re_theta),
    }
}

/// Assemble local BL Newton system
pub fn blsys(
    station1: &BLStation,
    station2: &BLStation,
    regime: BLRegime,
) -> (Matrix4x5, Matrix4x5, Vector4) {
    // Compute closure relations
    let vars1 = blvar(station1, regime);
    let vars2 = blvar(station2, regime);

    // Finite-difference discretization
    let (vs1, vs2, vsrez) = discretize_bl_equations(
        station1, station2, &vars1, &vars2
    );

    (vs1, vs2, vsrez)
}
```

## xoper.f - Main Operations and Viscous-Inviscid Coupling

This module contains the primary analysis routines including the VISCAL viscous-inviscid coupling loop.

### Key Subroutines

#### VISCAL (Line ~2907)

**Purpose:** Converges viscous operating point through Newton iteration.

**Called by:** SPECAL, SPECCL (when prescribing alpha or CL)

**Calls:**
- `XYWAKE` - Calculate wake trajectory
- `QWCALC` - Set wake velocities
- `QISET` - Set inviscid velocities
- `STFIND` - Locate stagnation point
- `IBLPAN` - Set BL→panel pointers
- `XICALC` - Calculate surface arc length
- `IBLSYS` - Set BL→system pointers
- `UICALC` - Set inviscid edge velocity
- `QDCALC` - Compute source influence matrix
- `SETBL` - Fill Newton system
- `BLSOLV` - Solve Newton system
- `UPDATE` - Update BL variables
- `QVFUE` - Calculate edge velocities from mass defect
- `GAMQV` - Set GAM from QVIS
- `STMOVE` - Relocate stagnation point
- `CLCALC` - Calculate lift coefficient
- `CDCALC` - Calculate drag coefficient

**Input Variables:**
| Variable | Source | Description |
|----------|--------|-------------|
| `NITER1` | Argument | Maximum iterations |
| `ALFA` | COMMON | Angle of attack |
| `GAM(IQX)` | COMMON | Vortex distribution |
| `X, Y, S` | COMMON | Geometry arrays |
| `REINF` | COMMON | Reynolds number |
| `MINF` | COMMON | Mach number |

**Output Variables:**
| Variable | Description |
|----------|-------------|
| `LVCONV` | TRUE if converged |
| `CL, CD, CM` | Force coefficients |
| `DSTR, THET` | BL thickness arrays |
| `UEDG` | Edge velocity |
| `SIG` | Mass defect source |
| `QVIS` | Viscous velocity |
| `CPV` | Viscous Cp |

**Algorithm:**

```
1. Initialize
   ├── Calculate wake trajectory if needed (XYWAKE)
   ├── Set wake velocities for alpha = 0, 90 (QWCALC)
   └── Set initial velocities (QISET)

2. Set up BL pointers (if not done)
   ├── Locate stagnation point (STFIND)
   ├── Set BL→panel pointers (IBLPAN)
   ├── Calculate arc length (XICALC)
   └── Set BL→system pointers (IBLSYS)

3. Set inviscid edge velocity (UICALC)
   └── UEDG = UINV if not initialized

4. Set up source influence matrix (QDCALC)
   └── DIJ(i,j) = dQtan_i/dSig_j

5. Newton iteration loop (ITER = 1 to NITER)
   │
   ├── Fill Newton system (SETBL)
   │   └── March BL, build Jacobian
   │
   ├── Solve Newton system (BLSOLV)
   │   └── Block elimination with underrelaxation
   │
   ├── Update BL variables (UPDATE)
   │   └── Apply deltas with limits
   │
   ├── Update flow field
   │   ├── Set new Mach, Re from CL (MRCL)
   │   ├── Set new velocities (QISET, UICALC) if CL specified
   │   ├── Calculate edge velocities (QVFUE)
   │   ├── Set GAM from QVIS (GAMQV)
   │   └── Relocate stagnation point (STMOVE)
   │
   ├── Calculate coefficients (CLCALC, CDCALC)
   │
   └── Check convergence
       └── RMSBL < EPS1 (1e-4) → converged

6. Finalize
   ├── Calculate Cp distributions
   └── Calculate hinge moment if flapped
```

**Convergence Criterion:**
```
RMS of Newton residuals < 1.0e-4
```

**Key Mathematical Relations:**

The coupling is achieved by:
1. BL mass defect creates equivalent source distribution:
   ```
   SIG_i = (MASS_i - MASS_{i-1}) / (S_i - S_{i-1})
   ```

2. Source distribution modifies surface velocity:
   ```
   QVIS_i = QINV_i + SUM_j(DIJ(i,j) * SIG_j)
   ```

3. Modified velocity becomes new edge velocity for BL:
   ```
   UEDG = |QVIS|
   ```

---

#### SPECAL (Line ~1405)

**Purpose:** Converges on specified angle of attack.

**Calls:** `VISCAL`

**Algorithm:**
1. Set `LALFA = .TRUE.` (alpha specified)
2. Set `ALFA` from input
3. Call `VISCAL` for viscous iteration

---

#### SPECCL (Line ~1500)

**Purpose:** Converges on specified lift coefficient.

**Calls:** `VISCAL`

**Algorithm:**
1. Set `LALFA = .FALSE.` (CL specified)
2. Set target CL
3. Call `VISCAL` which iterates on both BL and alpha

---

### yFoil Equivalent

`src/solver/viscal.rs`

The Rust implementation follows the same structure:

```rust
pub fn viscal(
    state: &mut SolverState,
    max_iter: usize,
) -> Result<ViscalResult, SolverError> {
    // Initialize wake if needed
    if !state.wake_computed {
        compute_wake(&mut state.wake, &state.panel);
    }

    // Set up BL pointers
    if !state.bl_initialized {
        find_stagnation_point(state);
        setup_bl_pointers(state);
    }

    // Compute DIJ matrix
    compute_dij_matrix(state);

    // Newton iteration
    for iter in 0..max_iter {
        // Fill Newton system
        let residual = setup_bl_system(state);

        // Solve and update
        solve_bl_system(state);
        update_bl_variables(state);

        // Update velocities
        compute_qvis(state);
        update_gamma(state);

        // Check convergence
        if residual.rms < 1e-4 {
            return Ok(ViscalResult::Converged { ... });
        }
    }

    Ok(ViscalResult::NotConverged { ... })
}
```

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

### yFoil mapping

See [the mapping section](#xfoil-yfoil-mapping) (tables 1–4 cover this
block) and the rules in [`docs/conventions/naming.md`](conventions/naming.md).

## XBL.INC - Boundary Layer Local Variables

This include file contains local variables used during boundary layer calculations. These are "station" variables representing the BL state at a single location, with derivatives for the Newton solver.

### Overview

The BL solver uses two stations:
- **Station 1**: Upstream station (known values)
- **Station 2**: Downstream station (being solved)

Each variable at a station includes derivatives with respect to primary unknowns for implicit Newton solution.

### Parameter

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

### yFoil mapping

See [the mapping section](#xfoil-yfoil-mapping) (tables 1–4 cover this
block) and the rules in [`docs/conventions/naming.md`](conventions/naming.md).

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

### Physical Basis

These parameters originate from:

1. **Green's lag entrainment method** (SCCON)
   - Models the finite response time of turbulence to pressure gradients

2. **Bradshaw's structural parameter** (GACON, GBCON)
   - Relates shear stress to mean velocity profile shape

3. **Skin friction correlations** (CFFAC)
   - Empirical adjustment for friction prediction

4. **Abu-Ghannam & Shaw transition** (CTRCON, CTRCEX)
   - e^N method calibration constants

### yFoil Mapping

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

### Notes

1. These values must match XFOIL exactly for numerical agreement
2. The constants are interdependent - changing one may require adjusting others
3. CFFAC allows post-hoc calibration without changing other parameters
