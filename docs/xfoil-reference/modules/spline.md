# spline.f - Spline Utilities

This module contains cubic spline routines used for smooth interpolation of geometry and flow variables.

## Key Subroutines

### SPLINE (Line ~1)

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

### SPLIND (Line ~50)

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

### SEVAL (Line ~100)

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

### DEVAL (Line ~150)

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

### D2VAL (Line ~200)

**Purpose:** Evaluates spline second derivative at a given parameter.

**Output:**
| Return value | Description |
|--------------|-------------|
| `D2VAL` | Interpolated d²X/dS² value |

---

## Usage in XFOIL

### Geometry Splines

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

### Curvature

Local curvature κ is:
```
κ = (X'*Y'' - Y'*X'') / (X'² + Y'²)^(3/2)
```

where primes denote derivatives with respect to arc length S.

### Panel Normal

Normal vector at point I:
```
NX(I) = YP(I) / sqrt(XP(I)² + YP(I)²)
NY(I) = -XP(I) / sqrt(XP(I)² + YP(I)²)
```

---

## yFoil Equivalent

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
