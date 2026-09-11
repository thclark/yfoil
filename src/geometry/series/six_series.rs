//! The NACA 6- and 6A-series basic thickness forms.
//!
//! These families have no closed form. Each is defined (Theodorsen 1931, Theodorsen and Garrick
//! 1933; Abbott, von Doenhoff and Stivers 1945) by the conformal transformation of a circle
//! through two tabulated functions of the circle angle φ — ε(φ), the angular distortion, and
//! ψ(φ), the radial one — scaled by a factor that sets the thickness ratio. The NASA ordinate
//! program (Ladson and Brooks 1974; Ladson, Brooks, Hill and Sproles 1996, TM-4741) carries
//! 201-point tables of ε and ψ for the 63, 64, 65, 66, 67, 63A, 64A and 65A families and maps
//! them as follows (`SetSixDigitPoints` in PDAS `naca456`, Carmichael 2001):
//!
//! ```text
//! z      = exp(ψ₀ + iφ)                      the circle
//! z'     = z · exp((ψ − ψ₀) − iε)            the near-circle
//! ζ      = z' + 1/z'                          the Joukowski transform
//! x + iy = (ζ₀ − ζ) / |ζ₂₀₀ − ζ₀|             normalised so the chord runs 0 → 1
//! ```
//!
//! with the scale factor applied to ε and ψ before mapping. The result is 201 points from the
//! leading edge (φ = 0) to the trailing edge (φ = π); ordinates at other stations come from the
//! arc-length spline through them (`spline_fmm`). TM-4741 states the program reproduces the
//! published ordinates to within 5 × 10⁻⁵ chord, and that the original ε/ψ graphs are lost: the
//! tables are the definition.

use super::six_series_tables as tables;
use super::spline_fmm::{fmm_spline, segment_inverse, spline_eval, Segment};

/// The eight tabulated families
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SixSeriesFamily {
    F63,
    F64,
    F65,
    F66,
    F67,
    F63A,
    F64A,
    F65A,
}

impl SixSeriesFamily {
    pub const ALL: [SixSeriesFamily; 8] = [
        Self::F63,
        Self::F64,
        Self::F65,
        Self::F66,
        Self::F67,
        Self::F63A,
        Self::F64A,
        Self::F65A,
    ];

    /// The designation prefix: "63", "64", …, "63A", "64A", "65A"
    pub fn label(self) -> &'static str {
        match self {
            Self::F63 => "63",
            Self::F64 => "64",
            Self::F65 => "65",
            Self::F66 => "66",
            Self::F67 => "67",
            Self::F63A => "63A",
            Self::F64A => "64A",
            Self::F65A => "65A",
        }
    }

    pub fn from_label(label: &str) -> Option<Self> {
        let upper = label.to_ascii_uppercase();
        Self::ALL.iter().copied().find(|f| f.label() == upper)
    }

    /// The A-series (Loftin 1948): straight-sided aft of 0.8 chord, used with the 6A mean line
    pub fn is_a_series(self) -> bool {
        matches!(self, Self::F63A | Self::F64A | Self::F65A)
    }

    fn tables(self) -> (&'static [f64; tables::N_PHI], &'static [f64; tables::N_PHI]) {
        match self {
            Self::F63 => (&tables::EPSILON_63, &tables::PSI_63),
            Self::F64 => (&tables::EPSILON_64, &tables::PSI_64),
            Self::F65 => (&tables::EPSILON_65, &tables::PSI_65),
            Self::F66 => (&tables::EPSILON_66, &tables::PSI_66),
            Self::F67 => (&tables::EPSILON_67, &tables::PSI_67),
            Self::F63A => (&tables::EPSILON_63A, &tables::PSI_63A),
            Self::F64A => (&tables::EPSILON_64A, &tables::PSI_64A),
            Self::F65A => (&tables::EPSILON_65A, &tables::PSI_65A),
        }
    }

    /// Coefficients of the scale-factor polynomial in t/c (`ScaleFactor`, Carmichael 2001):
    /// constant, linear, quadratic, cubic, quartic.
    fn scale_coefficients(self) -> [f64; 5] {
        match self {
            Self::F63 => [0.0, 8.1827699, 1.3776209, -0.092851684, 7.5942563],
            Self::F64 => [0.0, 4.6535511, 1.038063, -1.5041794, 4.7882784],
            Self::F65 => [0.0, 6.5718716, 0.49376292, 0.7319794, 1.9491474],
            Self::F66 => [0.0, 6.7581414, 0.19253769, 0.81282621, 0.85202897],
            Self::F67 => [0.0, 6.627289, 0.098965859, 0.96759774, 0.90537584],
            Self::F63A => [0.0, 8.1845925, 1.0492569, 1.31150930, 4.4515579],
            Self::F64A => [0.0, 8.2125018, 0.76855961, 1.4922345, 3.6130133],
            Self::F65A => [0.0, 8.2514822, 0.46569361, 1.50113018, 2.0908904],
        }
    }
}

/// The factor by which the tabulated ε and ψ are multiplied to give thickness ratio `t`
pub fn scale_factor(family: SixSeriesFamily, t: f64) -> f64 {
    let c = family.scale_coefficients();
    // Horner's rule, as `Polynomial` writes it
    let mut f = c[4];
    for &ck in c[..4].iter().rev() {
        f = f * t + ck;
    }
    f
}

/// A complex number, enough for the mapping
#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn exp(self) -> Self {
        let (s, c) = self.im.sin_cos();
        let e = self.re.exp();
        Self { re: e * c, im: e * s }
    }
    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
    fn recip(self) -> Self {
        let d = self.re * self.re + self.im * self.im;
        Self {
            re: self.re / d,
            im: -self.im / d,
        }
    }
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}

/// The 201 mapped points of the basic thickness form at thickness ratio `t`, leading edge to
/// trailing edge, as `(x, y_t)` on the upper surface (`SetSixDigitPoints`).
pub fn basic_thickness_form(family: SixSeriesFamily, t: f64) -> (Vec<f64>, Vec<f64>) {
    let (eps_table, psi_table) = family.tables();
    let n = tables::N_PHI;
    let sf = scale_factor(family, t);
    let psi0 = sf * psi_table[0];
    let mut zeta = Vec::with_capacity(n);
    for k in 0..n {
        let phi = std::f64::consts::PI * k as f64 / (n - 1) as f64;
        let eps = sf * eps_table[k];
        let psi = sf * psi_table[k];
        let z = Complex { re: psi0, im: phi }.exp();
        let zprime = z.mul(
            Complex {
                re: psi - psi0,
                im: -eps,
            }
            .exp(),
        );
        let zr = zprime.recip();
        zeta.push(Complex {
            re: zprime.re + zr.re,
            im: zprime.im + zr.im,
        });
    }
    let chord = Complex {
        re: zeta[n - 1].re - zeta[0].re,
        im: zeta[n - 1].im - zeta[0].im,
    }
    .abs();
    let mut x: Vec<f64> = zeta.iter().map(|z| (zeta[0].re - z.re) / chord).collect();
    let mut y: Vec<f64> = zeta.iter().map(|z| -(zeta[0].im - z.im) / chord).collect();
    // the trailing edge closes by construction (ε = ψ = 0 at φ = π); remove the round-off of
    // sin π so that the two surfaces meet exactly
    x[n - 1] = 1.0;
    y[n - 1] = 0.0;
    (x, y)
}

/// A 6-series thickness form ready to evaluate: the arc-length spline through the mapped points
/// (`Thickness6`). The spline is fitted around the whole symmetric section, trailing edge →
/// upper → leading edge → lower → trailing edge, and the lower half is inverted for `x`.
#[derive(Debug, Clone)]
pub struct SixSeriesForm {
    pub family: SixSeriesFamily,
    pub t: f64,
    /// Arc coordinate of the lower half, leading edge (index 0) to trailing edge
    s: Vec<f64>,
    x: Vec<f64>,
    /// Half thickness (positive) along the lower half
    y: Vec<f64>,
    x_d_s: Vec<f64>,
    y_d_s: Vec<f64>,
}

impl SixSeriesForm {
    pub fn new(family: SixSeriesFamily, t: f64) -> Self {
        let (xt, yt) = basic_thickness_form(family, t);
        let n = xt.len();
        // ParametrizeAirfoil: upper surface reversed (TE → LE), then the lower surface (LE → TE)
        let nn = 2 * n - 1;
        let mut x = Vec::with_capacity(nn);
        let mut y = Vec::with_capacity(nn);
        for k in (0..n).rev() {
            x.push(xt[k]);
            y.push(yt[k]);
        }
        for k in 1..n {
            x.push(xt[k]);
            y.push(-yt[k]);
        }
        let mut s = vec![0.0; nn];
        for k in 1..nn {
            let (dx, dy) = (x[k] - x[k - 1], y[k] - y[k - 1]);
            s[k] = s[k - 1] + (dx * dx + dy * dy).sqrt();
        }
        let xp = fmm_spline(&s, &x);
        let yp = fmm_spline(&s, &y);
        // the lower half, with the sign of y flipped so the thickness is positive
        let lo = n - 1;
        Self {
            family,
            t,
            s: s[lo..].to_vec(),
            x: x[lo..].to_vec(),
            y: y[lo..].iter().map(|v| -v).collect(),
            x_d_s: xp[lo..].to_vec(),
            y_d_s: yp[lo..].iter().map(|v| -v).collect(),
        }
    }

    /// The arc coordinate at which the form passes chord station `x`
    fn s_at(&self, x: f64) -> f64 {
        let n = self.x.len();
        if x <= self.x[0] {
            return self.s[0];
        }
        if x >= self.x[n - 1] {
            return self.s[n - 1];
        }
        let k = super::spline_fmm::lookup(&self.x, x);
        let seg = Segment::new(
            self.s[k],
            self.x[k],
            self.x_d_s[k],
            self.s[k + 1],
            self.x[k + 1],
            self.x_d_s[k + 1],
        );
        segment_inverse(&seg, x)
    }

    /// Half thickness `y_t` and its slope `dy_t/dx` at chord station `x` (the slope is 0 where
    /// `dx/ds` vanishes, at the leading edge)
    pub fn at(&self, x: f64) -> (f64, f64) {
        let s = self.s_at(x);
        let (_, x_d_s, _) = spline_eval(&self.s, &self.x, &self.x_d_s, s);
        let (y, y_d_s, _) = spline_eval(&self.s, &self.y, &self.y_d_s, s);
        let slope = if x_d_s == 0.0 { 0.0 } else { y_d_s / x_d_s };
        (y, slope)
    }

    /// The chord station actually reached for a requested `x` (a check on the inversion)
    pub fn x_at(&self, x: f64) -> f64 {
        let s = self.s_at(x);
        spline_eval(&self.s, &self.x, &self.x_d_s, s).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mapped_form_runs_from_le_to_te_on_the_upper_surface() {
        for family in SixSeriesFamily::ALL {
            let (x, y) = basic_thickness_form(family, 0.12);
            assert_eq!(x.len(), 201);
            assert!(x[0].abs() < 1e-15, "{}: x(0) = {}", family.label(), x[0]);
            assert!((x[200] - 1.0).abs() < 1e-14, "{}: x(π) = {}", family.label(), x[200]);
            assert!(y[200].abs() < 1e-14, "{}: y(π) = {}", family.label(), y[200]);
            for k in 1..201 {
                assert!(x[k] > x[k - 1], "{}: x not increasing at {k}", family.label());
                assert!(y[k] >= 0.0 || k == 200, "{}: negative ordinate at {k}", family.label());
            }
            let t_max = 2.0 * y.iter().cloned().fold(0.0_f64, f64::max);
            assert!(
                (t_max - 0.12).abs() < 2e-3,
                "{}: max thickness {t_max} for t = 0.12",
                family.label()
            );
        }
    }

    #[test]
    fn form_inverts_x_to_round_off() {
        let form = SixSeriesForm::new(SixSeriesFamily::F63, 0.15);
        for &x in &[1e-5, 0.005, 0.0125, 0.1, 0.3, 0.5, 0.9, 0.999] {
            assert!((form.x_at(x) - x).abs() < 1e-15, "x = {x}: reached {}", form.x_at(x));
        }
        assert_eq!(form.at(0.0).0, 0.0);
        assert!(form.at(1.0).0.abs() < 1e-14);
    }

    #[test]
    fn scale_factor_is_the_published_polynomial() {
        // t/c = 0.15 on the 63 family, Horner on the tabulated coefficients
        let t: f64 = 0.15;
        let expect = 8.1827699 * t + 1.3776209 * t * t - 0.092851684 * t * t * t + 7.5942563 * t * t * t * t;
        assert!((scale_factor(SixSeriesFamily::F63, t) - expect).abs() < 1e-15);
    }
}
