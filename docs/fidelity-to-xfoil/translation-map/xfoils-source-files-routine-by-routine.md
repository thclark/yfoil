# XFOIL's source files, routine by routine

A walk through the translated routines of each XFOIL source file: what each does, its inputs and
outputs, and where yFoil translates it. Line numbers are approximate.

## spline.f - Spline Utilities

This module contains cubic spline routines used for smooth interpolation of geometry and flow variables.

### Key subroutines of `spline.f`

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

### How XFOIL uses `spline.f`

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

### `spline.f` in yFoil

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

### Key subroutines of `xpanel.f`

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

### `xpanel.f` in yFoil

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

### Key subroutines of `xbl.f`

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

### `xbl.f` in yFoil

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

### Key subroutines of `xblsys.f`

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

### Helper subroutines of `xblsys.f`

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

### `xblsys.f` in yFoil

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

### Key subroutines of `xoper.f`

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

### `xoper.f` in yFoil

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
