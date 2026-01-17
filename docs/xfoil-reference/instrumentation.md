# XFOIL Instrumentation Guide

This document describes how to instrument XFOIL to capture intermediate values for validation against YFoil.

## Current Instrumentation

The XFOIL source code in `xfoil/xfoil6.99/src/` has been modified with diagnostic output:

### xoper.f - VISCAL Coupling

Location: `xoper.f:2982-3104`

Outputs:
- `/tmp/xfoil_dij.dat` - Source influence matrix DIJ
- `/tmp/xfoil_viscal_iter.dat` - Iteration-by-iteration BL state

### xfoil.f - Panel Geometry

Location: `xfoil.f:2124`

Outputs:
- `/tmp/xfoil_panels.dat` - Panel coordinates and normals

## Output Format

### DIJ Matrix (`/tmp/xfoil_dij.dat`)

```
=== DIJ MATRIX ===
N = 160
NW = 23
    1     1  1.234567890123456E+00
    1     2  2.345678901234567E-01
    ...
```

Format: `I J VALUE` where VALUE is in Fortran E24.16 format.

### VISCAL Iteration Log (`/tmp/xfoil_viscal_iter.dat`)

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

## Adding New Instrumentation

### Step 1: Locate Subroutine

Find the subroutine in the XFOIL source:

```bash
grep -n "SUBROUTINE BLKIN" xfoil/xfoil6.99/src/xblsys.f
```

### Step 2: Add WRITE Statements

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

### Step 3: Rebuild XFOIL

```bash
cd xfoil/xfoil6.99/bin
make clean
make
```

### Step 4: Run Test Case

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

### Step 5: Parse Output

Use the fixture generation script to convert output to JSON.

## Key Subroutines to Instrument

### Panel Method

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| PSILIN | xpanel.f | PSI, DZDG, DZDM per panel |
| QDCALC | xpanel.f | DIJ matrix, CIJ matrix |

### BL System

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| BLPRV | xblsys.f | X2, U2, T2, D2, U2_UEI, U2_MS |
| BLKIN | xblsys.f | M2, R2, H2, HK2, RT2, V2 + derivatives |
| BLVAR | xblsys.f | HS2, CF2, DI2, US2, HC2, DE2 + derivatives |
| BLSYS | xblsys.f | VS1, VS2, VSREZ matrices |
| TRCHEK2 | xblsys.f | AMPL2, XT, transition derivatives |

### BL Marching

| Subroutine | File | Variables to Capture |
|------------|------|---------------------|
| SETBL | xbl.f | VA, VB, VDEL, VM matrices |
| MRCHUE | xbl.f | Per-station THET, DSTR, CTAU |
| UPDATE | xbl.f | RMSBL, RMXBL, RLX, changes |

## Precision Requirements

**All output must use E24.16 format** to capture full double-precision values:

```fortran
WRITE(LU,'(A,E24.16)') 'VAR =', VAR
```

This gives 16 significant figures, matching IEEE double precision.

## Comparison Tolerance

When comparing YFoil to XFOIL:

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

## Existing Debug Files

The instrumented XFOIL writes to:
- `/tmp/xfoil_dij.dat`
- `/tmp/xfoil_viscal_iter.dat`
- `/tmp/xfoil_inviscid.dat`
- `/tmp/xfoil_panels.dat`
- `/tmp/xfoil_bl_debug.dat` (from xbl.f MRCHUE)

These files are overwritten on each run. Copy them before running another case.
