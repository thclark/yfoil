# xblsys.f - Boundary Layer Newton System

This module contains the core BL closure relations and the Newton system assembly routines.

## Key Subroutines

### BLPRV (Line ~701)

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

### BLKIN (Line ~725)

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

### BLVAR (Line ~784)

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

### BLSYS (Line ~583)

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

### TRCHEK2 (Line ~231)

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

## Helper Subroutines

### HKIN (Kinematic Shape Factor)

```
Hk = (H - 0.29*M²) / (1 + 0.113*M² + 0.0056*M⁴)
```

### DAMPL (Amplification Rate)

The Drela-Giles envelope e^N correlation.

### CFL/CFT (Skin Friction)

Laminar and turbulent Cf correlations.

### HSL/HST (Energy Shape Factor)

H* correlations for laminar and turbulent flow.

---

## YFoil Equivalent

`src/bl/closure.rs`, `src/bl/system.rs`

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
