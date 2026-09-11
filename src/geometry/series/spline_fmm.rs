//! The spline machinery of the NASA/PDAS NACA ordinate program (`splprocs.f90`, Carmichael 2001):
//! the cubic spline with Forsythe–Malcolm–Moler end conditions, Hermite evaluation of one segment
//! with derivatives, and the inversion of a monotone segment. Used by the 6-series thickness form,
//! which is defined by a table of mapped points rather than a formula.
//!
//! XFOIL's own splines (`geometry::spline`) are a different algorithm (SPLINE/SEVAL with its own
//! end conditions); the 6-series ordinates are defined by *this* interpolant, so it is kept
//! separate and ported as written.

/// Cubic spline slopes with FMM end conditions (`FMMspline`): the slope at every knot of the
/// spline through `(x, y)`, the end slopes chosen so that the third derivative is continuous at the
/// second and penultimate knots (a "not-a-knot" condition).
pub fn fmm_spline(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert_eq!(n, y.len(), "fmm_spline: x and y must have the same length");
    let mut yp = vec![0.0; n];
    if n < 2 {
        return yp;
    }
    let dx: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let delta: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / dx[i]).collect();
    if n == 2 {
        yp[0] = delta[0];
        yp[1] = delta[0];
        return yp;
    }
    let dd: Vec<f64> = (0..n - 2).map(|i| delta[i + 1] - delta[i]).collect();
    if n == 3 {
        let deriv2 = dd[0] / (x[2] - x[0]);
        let deriv1 = delta[0] - deriv2 * dx[0];
        yp[0] = deriv1;
        yp[1] = deriv1 + deriv2 * dx[0];
        yp[2] = deriv1 + deriv2 * (x[2] - x[0]);
        return yp;
    }

    let mut alpha = vec![0.0; n];
    let mut beta = vec![0.0; n];
    let mut sigma = vec![0.0; n];
    alpha[0] = -dx[0];
    for i in 1..n - 1 {
        alpha[i] = 2.0 * (dx[i - 1] + dx[i]);
    }
    for i in 1..n - 1 {
        alpha[i] -= dx[i - 1] * dx[i - 1] / alpha[i - 1];
    }
    alpha[n - 1] = -dx[n - 2] - dx[n - 2] * dx[n - 2] / alpha[n - 2];

    beta[0] = dd[1] / (x[3] - x[1]) - dd[0] / (x[2] - x[0]);
    beta[0] = beta[0] * dx[0] * dx[0] / (x[3] - x[0]);
    beta[1..n - 1].copy_from_slice(&dd[..n - 2]);
    beta[n - 1] = dd[n - 3] / (x[n - 1] - x[n - 3]) - dd[n - 4] / (x[n - 2] - x[n - 4]);
    beta[n - 1] = -beta[n - 1] * dx[n - 2] * dx[n - 2] / (x[n - 1] - x[n - 4]);
    for i in 1..n {
        beta[i] -= dx[i - 1] * beta[i - 1] / alpha[i - 1];
    }

    sigma[n - 1] = beta[n - 1] / alpha[n - 1];
    for i in (0..n - 1).rev() {
        sigma[i] = (beta[i] - dx[i] * sigma[i + 1]) / alpha[i];
    }
    for i in 0..n - 1 {
        yp[i] = delta[i] - dx[i] * (sigma[i] + sigma[i] + sigma[i + 1]);
    }
    yp[n - 1] = yp[n - 2] + dx[n - 2] * 3.0 * (sigma[n - 1] + sigma[n - 2]);
    yp
}

/// One Hermite segment of a spline: the cubic on `[a, b]` with values `fa`, `fb` and slopes
/// `fpa`, `fpb` at the ends (`EvaluateCubicAndDerivs`, mapped onto t ∈ [0, 1]).
#[derive(Debug, Clone, Copy)]
pub struct Segment {
    pub a: f64,
    pub b: f64,
    coef: [f64; 4],
}

impl Segment {
    pub fn new(a: f64, fa: f64, fpa: f64, b: f64, fb: f64, fpb: f64) -> Self {
        let rhs = [fa, fb, fpa * (b - a), fpb * (b - a)];
        let coef = [
            2.0 * rhs[0] - 2.0 * rhs[1] + rhs[2] + rhs[3],
            -3.0 * rhs[0] + 3.0 * rhs[1] - 2.0 * rhs[2] - rhs[3],
            rhs[2],
            rhs[0],
        ];
        Self { a, b, coef }
    }

    /// `(f, f', f'')` at `u`
    pub fn eval(&self, u: f64) -> (f64, f64, f64) {
        let c = &self.coef;
        let h = 1.0 / (self.b - self.a);
        let t = (u - self.a) * h;
        let f = c[3] + t * (c[2] + t * (c[1] + t * c[0]));
        let fp = h * (c[2] + t * (2.0 * c[1] + t * 3.0 * c[0]));
        let fpp = h * h * (2.0 * c[1] + t * 6.0 * c[0]);
        (f, fp, fpp)
    }
}

/// Index `k` with `x[k] <= u < x[k + 1]` in an increasing table (`Lookup`), clamped so that
/// `k + 1` is a valid index; `u` outside the table gives the end segment.
pub fn lookup(x: &[f64], u: f64) -> usize {
    let n = x.len();
    if u < x[0] {
        return 0;
    }
    if u >= x[n - 1] {
        return n - 2;
    }
    let (mut i, mut j) = (0usize, n - 1);
    while j > i + 1 {
        let k = (i + j) / 2;
        if u < x[k] {
            j = k;
        } else {
            i = k;
        }
    }
    i
}

/// A spline `f(x)` given by knots, values and slopes: `(f, f', f'')` at `u` (`PClookup`).
pub fn spline_eval(x: &[f64], f: &[f64], fp: &[f64], u: f64) -> (f64, f64, f64) {
    let k = lookup(x, u);
    let seg = Segment::new(x[k], f[k], fp[k], x[k + 1], f[k + 1], fp[k + 1]);
    let (v, vp, vpp) = seg.eval(u);
    // at a knot the spline takes the knot's value and slope by definition; the Hermite sum
    // reproduces them only to round-off
    if u == x[k] {
        return (f[k], fp[k], vpp);
    }
    if u == x[k + 1] {
        return (f[k + 1], fp[k + 1], vpp);
    }
    (v, vp, vpp)
}

/// The parameter `u` in `[a, b]` at which a monotone segment takes the value `target`. Newton
/// iteration safeguarded by bisection, run to the round-off of `u`; the segment must bracket
/// `target` (`SplineZero` solves the same problem with Brent's method at a 1e-6 tolerance — see
/// the tolerance note in `tests/utilities/tolerances.rs`).
pub fn segment_inverse(seg: &Segment, target: f64) -> f64 {
    let (mut lo, mut hi) = (seg.a, seg.b);
    let (flo, fhi) = (seg.eval(lo).0 - target, seg.eval(hi).0 - target);
    if flo == 0.0 {
        return lo;
    }
    if fhi == 0.0 {
        return hi;
    }
    let increasing = fhi > flo;
    let mut u = 0.5 * (lo + hi);
    for _ in 0..100 {
        let (f, fp, _) = seg.eval(u);
        let r = f - target;
        if r == 0.0 {
            return u;
        }
        if (r > 0.0) == increasing {
            hi = u;
        } else {
            lo = u;
        }
        let step = if fp != 0.0 { u - r / fp } else { f64::NAN };
        let next = if step.is_finite() && step > lo && step < hi {
            step
        } else {
            0.5 * (lo + hi)
        };
        if next == u || (hi - lo) <= f64::EPSILON * u.abs().max(1.0) {
            return next;
        }
        u = next;
    }
    u
}

/// Linear interpolation in a short table (`TableLookup` with order 1): the two knots bracketing
/// `u`, extrapolated linearly outside the table.
pub fn table_lookup_linear(x: &[f64], y: &[f64], u: f64) -> f64 {
    let n = x.len();
    let j = lookup(x, u).min(n - 2);
    let (x0, x1, y0, y1) = (x[j], x[j + 1], y[j], y[j + 1]);
    if u == x0 {
        return y0;
    }
    if u == x1 {
        return y1;
    }
    // Lagrange form, as InterpolatePolynomial writes it
    y0 * (u - x1) / (x0 - x1) + y1 * (u - x0) / (x1 - x0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmm_spline_reproduces_a_cubic() {
        // a cubic is reproduced exactly by a cubic spline with not-a-knot ends
        let x: Vec<f64> = (0..9).map(|i| i as f64 * 0.25).collect();
        let f = |t: f64| 1.0 + t - 2.0 * t * t + 0.5 * t * t * t;
        let fp = |t: f64| 1.0 - 4.0 * t + 1.5 * t * t;
        let y: Vec<f64> = x.iter().map(|&t| f(t)).collect();
        let yp = fmm_spline(&x, &y);
        for (i, &t) in x.iter().enumerate() {
            assert!((yp[i] - fp(t)).abs() < 1e-12, "slope at {t}: {} vs {}", yp[i], fp(t));
        }
        let (v, vp, _) = spline_eval(&x, &y, &yp, 0.6);
        assert!((v - f(0.6)).abs() < 1e-12);
        assert!((vp - fp(0.6)).abs() < 1e-12);
    }

    #[test]
    fn segment_inverse_hits_round_off() {
        let seg = Segment::new(0.0, 0.0, 0.1, 1.0, 1.0, 2.5);
        for target in [0.01, 0.3, 0.7, 0.999] {
            let u = segment_inverse(&seg, target);
            assert!((seg.eval(u).0 - target).abs() < 1e-15, "target {target}");
        }
    }

    #[test]
    fn linear_table_lookup_matches_the_knots() {
        let x = [0.05, 0.1, 0.15, 0.2, 0.25];
        let y = [0.0580, 0.126, 0.2025, 0.29, 0.391];
        for i in 0..5 {
            assert_eq!(table_lookup_linear(&x, &y, x[i]), y[i]);
        }
        assert!((table_lookup_linear(&x, &y, 0.125) - 0.5 * (0.126 + 0.2025)).abs() < 1e-15);
    }
}
