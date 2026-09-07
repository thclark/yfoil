# xbl.f - Boundary Layer Marching

This module contains the BL marching routines that solve the integral boundary layer equations station by station.

## Key Subroutines

### SETBL (Line ~21)

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

### MRCHUE (Line ~542)

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

### MRCHDU (Line ~886)

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

### UPDATE (Line ~1264)

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

## YFoil Equivalent

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
