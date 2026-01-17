# XFOIL Variable and Control Flow Map

This document maps the variables updated within each subroutine and common block
during a typical XFOIL run (NACA 0012, repanel, viscous solve at alpha=0, alpha=1).

Generated from instrumented XFOIL run on 2026-01-17.

## Overview of Key Common Blocks

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

## Control Flow: NACA 0012 Generation

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

### Variables Updated by NACA:

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

## Control Flow: PANGEN (Paneling)

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
├─ NORM(...)          → Normalizes by chord
└─ Output: N panels with X, Y, S, XP, YP, NX, NY, APANEL
```

### Variables Updated by PANGEN:

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
| CHORD       | CR05         | Chord length (normalized to 1.0) |
| SHARP       | CL01         | .TRUE. if sharp TE               |

## Control Flow: SPECAL (Inviscid Solve at Alpha)

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

### Variables Updated by SPECAL:

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

## Control Flow: VISCAL (Viscous-Inviscid Coupling)

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

### Variables Updated by VISCAL (per iteration):

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

## Control Flow: SETBL (BL Newton System Assembly)

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

### Closure Functions Called by SETBL:

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

## Subroutine Call Frequency (from instrumented run)

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

## Key State Flags

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

## Data Flow Summary

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
