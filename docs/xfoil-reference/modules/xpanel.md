# xpanel.f - Panel Method

This module contains the panel method routines for computing inviscid flow via vortex distribution.

## Key Subroutines

### PSILIN (Line ~99)

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

### QDCALC (Line ~1161)

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

### QVFUE (Line ~1400)

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

### UICALC (Line ~1300)

**Purpose:** Sets inviscid BL edge velocity UINV from panel velocity QINV.

**Algorithm:**

```
For each BL station:
   UINV(IBL,IS) = |QINV(IPAN(IBL,IS))| * VTI(IBL,IS)
```

where VTI = ±1 accounts for direction convention.

---

## Influence Coefficient Theory

### Vortex Panel Method

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

### Source Panels for Viscous Coupling

Mass defect creates equivalent sources:
```
σ = d(m)/ds = d(ρ*Ue*δ*)/ds
```

The source influence on velocity:
```
ΔQ_tan = ∫ σ(s')/(2π) * sin(θ)/r ds'
```

### Newton System Coupling

The global system couples:
1. Panel method: [AIJ][γ] = [RHS] (Kutta condition)
2. BL equations: f(θ, δ*, Ue) = 0
3. Mass defect: σ = g(Ue, δ*)
4. Velocity update: Ue_new = Ue_inv + DIJ * σ

---

## YFoil Equivalent

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
