//! NACA basic thickness forms: the half-thickness `y_t(x)` of the symmetrical section, from which
//! every cambered section is built by adding a mean line (`mean_line.rs`).
//!
//! Three families are defined:
//!
//! - **4-digit** (Jacobs, Ward and Pinkerton 1933, NACA Report 460): the closed form
//!   `y_t = 5t (0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1015x⁴)`. The −0.1015 coefficient
//!   leaves a trailing-edge half-thickness of 0.0105t (XFOIL's `NACA` command uses the same);
//!   −0.1036 would close it.
//! - **4-digit modified** (Stack and von Doenhoff 1934, NACA Report 492): two cubics joined at
//!   the position of maximum thickness, with the leading-edge radius index `I` (6 is the plain
//!   4-digit radius, 0 is a sharp nose) and the chordwise position of maximum thickness `M`
//!   (tenths of chord) chosen independently. The 16-series (Stack 1943, NACA Report 763) is the
//!   member `I = 4, M = 5`. Coefficients as in NASA TM-4741 (`Thickness4M`).
//! - **6- and 6A-series** (`six_series.rs`): no closed form; tabulated conformal maps.

use super::six_series::{SixSeriesFamily, SixSeriesForm};
use serde_json::{json, Value};

/// A basic thickness form. `t` is the maximum thickness as a fraction of chord in every family.
#[derive(Debug, Clone)]
pub enum ThicknessForm {
    FourDigit {
        t: f64,
    },
    FourDigitModified {
        t: f64,
        /// Leading-edge radius index `I`: `r_le = 1.1019 (t I / 6)²`
        le_radius_index: f64,
        /// Chordwise position of maximum thickness (fraction of chord)
        x_max_thickness: f64,
    },
    SixSeries(SixSeriesForm),
}

/// The 4-digit half-thickness at `x` for thickness ratio `t`. This expression, in this
/// association, is what the tracked XFOIL fixtures were generated with (`cargo xtask fixtures
/// --verify` asserts the panels bitwise), so it is kept verbatim.
pub fn four_digit_half_thickness(x: f64, t: f64) -> f64 {
    5.0 * t * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x.powi(2) + 0.2843 * x.powi(3) - 0.1015 * x.powi(4))
}

/// `dy_t/dx` of the 4-digit form; infinite at the leading edge
pub fn four_digit_slope(x: f64, t: f64) -> f64 {
    if x == 0.0 {
        return f64::INFINITY;
    }
    5.0 * t * (0.5 * 0.2969 / x.sqrt() - 0.1260 - 2.0 * 0.3516 * x + 3.0 * 0.2843 * x * x - 4.0 * 0.1015 * x * x * x)
}

/// The coefficients of a 4-digit-modified form: `(a0, a1, a2, a3)` of the forward cubic-in-√x
/// and `(d0, d1, d2, d3)` of the aft cubic in `(1 − x)`, on the t/c = 0.2 basis
/// (TM-4741 p. 6; `Thickness4M`, `CalA1`–`CalA3`, `CalculateD1` in `naca456`).
pub fn four_digit_modified_coefficients(t: f64, le_radius_index: f64, x_max_thickness: f64) -> ([f64; 4], [f64; 4]) {
    let xmt = x_max_thickness;
    // d1, the trailing-edge half-angle parameter: a quartic fit through the tabulated values
    // (M = 0.2 → 0.200, 0.3 → 0.234, 0.4 → 0.315, 0.5 → 0.465, 0.6 → 0.700)
    let d1 = {
        let c = [3.48E-5, 2.3076628, -10.127712, 19.961478, -10.420597];
        let mut f = c[4];
        for &ck in c[..4].iter().rev() {
            f = f * xmt + ck;
        }
        f
    };
    // leading-edge radius (Abbott and von Doenhoff p. 117): the plain 4-digit 1.1019 t² at I = 6
    let rle = 1.1019 / 36.0 * (t * le_radius_index) * (t * le_radius_index);
    let a0 = (0.2 / t) * (rle + rle).sqrt();
    let a3 = {
        let omxmt = 1.0 - xmt;
        let v1 = 0.1 / (xmt * xmt * xmt);
        let v2 = (d1 * omxmt - 0.294) / (xmt * omxmt * omxmt);
        let v3 = (3.0 / 8.0) * a0 / xmt.powf(2.5);
        v1 + v2 - v3
    };
    let a2 = {
        let v1 = 0.1 / (xmt * xmt);
        let v2 = 0.5 * a0 / (xmt * xmt * xmt).sqrt();
        let v3 = 2.0 * a3 * xmt;
        -v1 + v2 - v3
    };
    let a1 = {
        let v1 = 0.5 * a0 / xmt.sqrt();
        let v2 = 2.0 * a2 * xmt;
        let v3 = 3.0 * a3 * xmt * xmt;
        -v1 - v2 - v3
    };
    let omxmt = 1.0 - xmt;
    let omxmsq = omxmt * omxmt;
    let d3 = ((3. * d1) - (0.588 / omxmt)) / (3. * omxmsq);
    let d2 = (-1.5 * omxmt * d3) - ((0.5 * d1) / omxmt);
    let d0 = 0.002;
    ([a0, a1, a2, a3], [d0, d1, d2, d3])
}

impl ThicknessForm {
    pub fn six_series(family: SixSeriesFamily, t: f64) -> Self {
        Self::SixSeries(SixSeriesForm::new(family, t))
    }

    /// Maximum thickness as a fraction of chord
    pub fn t(&self) -> f64 {
        match self {
            Self::FourDigit { t } | Self::FourDigitModified { t, .. } => *t,
            Self::SixSeries(f) => f.t,
        }
    }

    /// Whether the form closes at the trailing edge (the 6-series does; the 4-digit families keep
    /// a finite trailing-edge thickness)
    pub fn sharp_te(&self) -> bool {
        matches!(self, Self::SixSeries(_))
    }

    /// Half thickness `y_t` and slope `dy_t/dx` at chord station `x ∈ [0, 1]`
    pub fn at(&self, x: f64) -> (f64, f64) {
        match self {
            Self::FourDigit { t } => (four_digit_half_thickness(x, *t), four_digit_slope(x, *t)),
            Self::FourDigitModified {
                t,
                le_radius_index,
                x_max_thickness,
            } => {
                if x == 0.0 {
                    return (0.0, f64::INFINITY);
                }
                let (a, d) = four_digit_modified_coefficients(*t, *le_radius_index, *x_max_thickness);
                let (y, yp) = if x < *x_max_thickness {
                    let srx = x.sqrt();
                    (
                        a[0] * srx + x * (a[1] + x * (a[2] + x * a[3])),
                        0.5 * a[0] / srx + a[1] + x * (2.0 * a[2] + x * 3.0 * a[3]),
                    )
                } else {
                    let xx = 1.0 - x;
                    // the aft cubic is in (1 − x), so dy/dx is minus its derivative in xx
                    (
                        d[0] + xx * (d[1] + xx * (d[2] + xx * d[3])),
                        -(d[1] + xx * (2.0 * d[2] + xx * 3.0 * d[3])),
                    )
                };
                (5.0 * t * y, 5.0 * t * yp)
            }
            Self::SixSeries(f) => f.at(x),
        }
    }

    /// The provenance record of the form (the `thickness_form` entry of `Geometry.generator`)
    pub fn record(&self) -> Value {
        match self {
            Self::FourDigit { t } => json!({ "family": "4", "t": t }),
            Self::FourDigitModified {
                t,
                le_radius_index,
                x_max_thickness,
            } => json!({
                "family": "4M", "t": t, "le_radius_index": le_radius_index, "x_max_thickness": x_max_thickness
            }),
            Self::SixSeries(f) => json!({ "family": f.family.label(), "t": f.t }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn four_digit_modified_63_is_close_to_the_plain_four_digit() {
        // I = 6, M = 0.3 reproduces the 4-digit leading-edge radius and thickness position; the
        // ordinates agree to a few 1e-4 (the aft cubic differs from the quartic)
        let plain = ThicknessForm::FourDigit { t: 0.12 };
        let modified = ThicknessForm::FourDigitModified {
            t: 0.12,
            le_radius_index: 6.0,
            x_max_thickness: 0.3,
        };
        for &x in &[0.005, 0.05, 0.2, 0.3, 0.5, 0.8, 0.95] {
            let (yp, _) = plain.at(x);
            let (ym, _) = modified.at(x);
            assert!((yp - ym).abs() < 2e-3, "x = {x}: {yp} vs {ym}");
        }
        // the maximum thickness is t at x = M
        let (y, slope) = modified.at(0.3);
        assert!((2.0 * y - 0.12).abs() < 1e-12, "2 y_t(M) = {}", 2.0 * y);
        assert!(slope.abs() < 1e-12, "slope at M = {slope}");
    }

    #[test]
    fn four_digit_modified_slope_is_a_derivative() {
        let form = ThicknessForm::FourDigitModified {
            t: 0.12,
            le_radius_index: 4.0,
            x_max_thickness: 0.5,
        };
        for &x in &[0.1, 0.4, 0.6, 0.9] {
            let h = 1e-6;
            let fd = (form.at(x + h).0 - form.at(x - h).0) / (2.0 * h);
            assert!((fd - form.at(x).1).abs() < 1e-8, "x = {x}: fd {fd} vs {}", form.at(x).1);
        }
    }

    #[test]
    fn four_digit_matches_the_closed_form() {
        assert_eq!(four_digit_half_thickness(0.0, 0.12), 0.0);
        let x: f64 = 0.3;
        let expect =
            0.6 * (0.2969 * x.sqrt() - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
        assert!((four_digit_half_thickness(x, 0.12) - expect).abs() < 1e-16);
    }
}
