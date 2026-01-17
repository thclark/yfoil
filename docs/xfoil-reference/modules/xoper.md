# xoper.f - Main Operations and Viscous-Inviscid Coupling

This module contains the primary analysis routines including the VISCAL viscous-inviscid coupling loop.

## Key Subroutines

### VISCAL (Line ~2907)

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

### SPECAL (Line ~1405)

**Purpose:** Converges on specified angle of attack.

**Calls:** `VISCAL`

**Algorithm:**
1. Set `LALFA = .TRUE.` (alpha specified)
2. Set `ALFA` from input
3. Call `VISCAL` for viscous iteration

---

### SPECCL (Line ~1500)

**Purpose:** Converges on specified lift coefficient.

**Calls:** `VISCAL`

**Algorithm:**
1. Set `LALFA = .FALSE.` (CL specified)
2. Set target CL
3. Call `VISCAL` which iterates on both BL and alpha

---

## YFoil Equivalent

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
