//! NACA mean lines: the camber line `y_c(x)` and its slope, added to a basic thickness form
//! perpendicular to the line (`Section::surface_at`).
//!
//! - **2-digit** (Report 460): two parabolic arcs meeting at the maximum camber `m` at `x = p`.
//! - **3-digit** (Jacobs and Pinkerton 1935, NACA Report 537; the 5-digit sections): a cubic
//!   forward of `x = r` and a straight line aft, the constants `r` and `k₁` tabulated against
//!   the position of maximum camber for a design lift coefficient of 0.3 and scaled linearly.
//! - **3-digit reflex** (Report 537): the same with a second cubic aft, its constant `k₂/k₁`
//!   chosen for zero pitching moment.
//! - **6-series** `a` mean line (Abbott and von Doenhoff, eq. 4.26–4.27): uniform load from the
//!   leading edge to `x = a`, falling linearly to zero at the trailing edge.
//! - **6A** modified mean line (Loftin 1948, NACA Report 903): the `a = 0.8` line, scaled by
//!   0.97948 and replaced by a straight line aft of 0.86 chord, which with the 6A thickness
//!   forms gives straight-sided sections aft of 0.8 chord.
//!
//! The formulas are those of NASA TM-4741 as coded in PDAS `naca456` (`MeanLine2`, `MeanLine3`,
//! `MeanLine3Reflex`, `MeanLine6`, `MeanLine6M`); the 2-digit line keeps the expression the
//! tracked XFOIL fixtures were generated with.

use super::spline_fmm::table_lookup_linear;
use serde_json::{json, Value};

/// A mean line. `p` is the chordwise position of maximum camber (fraction of chord), `cl` the
/// design lift coefficient, `a` the extent of uniform loading.
#[derive(Debug, Clone, PartialEq)]
pub enum MeanLine {
    Symmetric,
    TwoDigit { m: f64, p: f64 },
    ThreeDigit { cl: f64, p: f64 },
    ThreeDigitReflex { cl: f64, p: f64 },
    SixSeries { a: f64, cl: f64 },
    SixSeriesModified { cl: f64 },
}

/// The 2-digit line, `(y_c, dy_c/dx)`. Kept verbatim: the tracked XFOIL fixtures were generated
/// with this expression.
pub fn two_digit_mean_line(x: f64, m: f64, p: f64) -> (f64, f64) {
    if m == 0.0 || p == 0.0 {
        // Symmetric airfoil
        return (0.0, 0.0);
    }

    let (yc, dyc) = if x < p {
        let yc = m / (p * p) * (2.0 * p * x - x * x);
        let dyc = 2.0 * m / (p * p) * (p - x);
        (yc, dyc)
    } else {
        let yc = m / ((1.0 - p).powi(2)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
        let dyc = 2.0 * m / ((1.0 - p).powi(2)) * (p - x);
        (yc, dyc)
    };

    (yc, dyc)
}

/// `(r, k₁)` of the 3-digit line at position of maximum camber `p` (TM-4741 p. 8, `GetRk1`):
/// linear interpolation in the table for p = 0.05, 0.10, 0.15, 0.20, 0.25.
pub fn three_digit_constants(p: f64) -> (f64, f64) {
    const M: [f64; 5] = [0.05, 0.1, 0.15, 0.2, 0.25];
    const R: [f64; 5] = [0.0580, 0.126, 0.2025, 0.29, 0.391];
    const K1: [f64; 5] = [361.4, 51.64, 15.957, 6.643, 3.23];
    (table_lookup_linear(&M, &R, p), table_lookup_linear(&M, &K1, p))
}

/// `(r, k₁, k₂/k₁)` of the reflex line (`GetRk1k2`). `k₂/k₁ = (3(r − p)² − r³)/(1 − r)³` is
/// Report 537's expression; TM X-3284 and TM-4741 both misprint it (Carmichael's note).
pub fn three_digit_reflex_constants(p: f64) -> (f64, f64, f64) {
    const M: [f64; 4] = [0.1, 0.15, 0.2, 0.25];
    const R: [f64; 4] = [0.13, 0.217, 0.318, 0.441];
    const K1: [f64; 4] = [51.99, 15.793, 6.52, 3.191];
    let r = table_lookup_linear(&M, &R, p);
    let k1 = table_lookup_linear(&M, &K1, p);
    let k21 = (3.0 * (r - p) * (r - p) - r * r * r) / ((1.0 - r) * (1.0 - r) * (1.0 - r));
    (r, k1, k21)
}

/// Below this distance from either end the logarithmic terms of the 6-series lines are replaced
/// by their limits (the reference program's `EPS`)
const END_GUARD: f64 = 1e-7;

/// The 6-series `a` mean line at `x` for unit design lift coefficient
fn six_series_unit(a: f64, x: f64) -> (f64, f64) {
    use std::f64::consts::PI;
    let two_pi = 2.0 * PI;
    let oma = 1.0 - a;
    let omx = 1.0 - x;
    if oma.abs() < END_GUARD {
        // a = 1: uniform load over the whole chord
        if x < END_GUARD || omx < END_GUARD {
            return (0.0, 0.0);
        }
        let ym = -(omx * omx.ln() + x * x.ln()) * (0.25 / PI);
        let ymp = (omx.ln() - x.ln()) * (0.25 / PI);
        return (ym, ymp);
    }
    if x < END_GUARD || omx.abs() < END_GUARD {
        return (0.0, 0.0);
    }
    let (g, h) = if a.abs() < END_GUARD {
        (-0.25, -0.5)
    } else {
        let g = -(a * a * (0.5 * a.ln() - 0.25) + 0.25) / oma;
        let h = g + (0.5 * oma * oma * oma.ln() - 0.25 * oma * oma) / oma;
        (g, h)
    };
    let amx = a - x;
    let (term1, term1p) = if amx.abs() < END_GUARD {
        (0.0, 0.0)
    } else {
        (amx * amx * (2.0 * amx.abs().ln() - 1.0), -amx * amx.abs().ln())
    };
    let term2 = omx * omx * (1.0 - 2.0 * omx.ln());
    let term2p = omx * omx.ln();
    let ym = 0.25 * (term1 + term2) / oma - x * x.ln() + g - h * x;
    let ymp = (term1p + term2p) / oma - 1.0 - x.ln() - h;
    (ym / (two_pi * (a + 1.0)), ymp / (two_pi * (a + 1.0)))
}

impl MeanLine {
    pub fn is_symmetric(&self) -> bool {
        match self {
            Self::Symmetric => true,
            Self::TwoDigit { m, p } => *m == 0.0 || *p == 0.0,
            Self::ThreeDigit { cl, .. } | Self::ThreeDigitReflex { cl, .. } => *cl == 0.0,
            Self::SixSeries { cl, .. } | Self::SixSeriesModified { cl } => *cl == 0.0,
        }
    }

    /// Camber `y_c` and slope `dy_c/dx` at chord station `x ∈ [0, 1]`
    pub fn at(&self, x: f64) -> (f64, f64) {
        match *self {
            Self::Symmetric => (0.0, 0.0),
            Self::TwoDigit { m, p } => two_digit_mean_line(x, m, p),
            Self::ThreeDigit { cl, p } => {
                if cl == 0.0 {
                    return (0.0, 0.0);
                }
                let (r, k1) = three_digit_constants(p);
                let (ym, ymp) = if x < r {
                    (
                        x * (x * (x - 3.0 * r) + r * r * (3.0 - r)),
                        3.0 * x * (x - r - r) + r * r * (3.0 - r),
                    )
                } else {
                    (r * r * r * (1.0 - x), -r * r * r)
                };
                // the tabulated k₁ is for cl = 0.3 and the line carries a factor 1/6
                let scale = k1 * cl / 1.8;
                (scale * ym, scale * ymp)
            }
            Self::ThreeDigitReflex { cl, p } => {
                if cl == 0.0 {
                    return (0.0, 0.0);
                }
                let (r, k1, k21) = three_digit_reflex_constants(p);
                let r3 = r * r * r;
                let mr3 = (1.0 - r) * (1.0 - r) * (1.0 - r);
                let xr = x - r;
                let (ym, ymp) = if x < r {
                    (
                        xr * xr * xr - k21 * mr3 * x - x * r3 + r3,
                        3.0 * xr * xr - k21 * mr3 - r3,
                    )
                } else {
                    (
                        k21 * xr * xr * xr - k21 * mr3 * x - x * r3 + r3,
                        3.0 * k21 * xr * xr - k21 * mr3 - r3,
                    )
                };
                let scale = k1 * cl / 1.8;
                (scale * ym, scale * ymp)
            }
            Self::SixSeries { a, cl } => {
                if cl == 0.0 {
                    return (0.0, 0.0);
                }
                let (ym, ymp) = six_series_unit(a, x);
                (cl * ym, cl * ymp)
            }
            Self::SixSeriesModified { cl } => {
                if cl == 0.0 {
                    return (0.0, 0.0);
                }
                let te_slope = -0.24521 * cl;
                if x > 0.86 {
                    return (te_slope * (x - 1.0), te_slope);
                }
                if x < END_GUARD {
                    return (0.0, 0.0);
                }
                let (ym, ymp) = six_series_unit(0.8, x);
                (cl * 0.97948 * ym, cl * 0.97948 * ymp)
            }
        }
    }

    /// The provenance record of the line (the `mean_line` entry of `Geometry.generator`)
    pub fn record(&self) -> Value {
        match *self {
            Self::Symmetric => Value::Null,
            Self::TwoDigit { m, p } => json!({ "family": "2", "m": m, "p": p }),
            Self::ThreeDigit { cl, p } => json!({ "family": "3", "cl": cl, "p": p }),
            Self::ThreeDigitReflex { cl, p } => json!({ "family": "3R", "cl": cl, "p": p }),
            Self::SixSeries { a, cl } => json!({ "family": "6", "cl": cl, "a": a }),
            Self::SixSeriesModified { cl } => json!({ "family": "6A", "cl": cl }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finite_difference(line: &MeanLine, x: f64) -> f64 {
        let h = 1e-6;
        (line.at(x + h).0 - line.at(x - h).0) / (2.0 * h)
    }

    #[test]
    fn slopes_are_derivatives() {
        let lines = [
            MeanLine::TwoDigit { m: 0.04, p: 0.4 },
            MeanLine::ThreeDigit { cl: 0.3, p: 0.15 },
            MeanLine::ThreeDigitReflex { cl: 0.3, p: 0.15 },
            MeanLine::SixSeries { a: 1.0, cl: 0.4 },
            MeanLine::SixSeries { a: 0.5, cl: 0.4 },
            MeanLine::SixSeries { a: 0.0, cl: 0.4 },
            MeanLine::SixSeriesModified { cl: 0.4 },
        ];
        for line in &lines {
            for &x in &[0.05, 0.1, 0.2, 0.3, 0.45, 0.7, 0.9, 0.95] {
                let (_, slope) = line.at(x);
                let fd = finite_difference(line, x);
                assert!((slope - fd).abs() < 1e-7, "{line:?} at {x}: slope {slope} vs fd {fd}");
            }
        }
    }

    #[test]
    fn three_digit_reflex_constants_reproduce_the_published_k2_k1() {
        // TM-4741 table (as corrected by Carmichael): p = 0.10 → 0.000764, 0.15 → 0.00677,
        // 0.20 → 0.0303, 0.25 → 0.1355
        for (p, k21) in [(0.10, 0.000764), (0.15, 0.00677), (0.20, 0.0303), (0.25, 0.1355)] {
            let (_, _, got) = three_digit_reflex_constants(p);
            assert!((got - k21).abs() / k21 < 2e-3, "p = {p}: {got} vs {k21}");
        }
    }

    #[test]
    fn six_series_line_is_zero_at_both_ends_and_peaks_near_mid_chord() {
        let line = MeanLine::SixSeries { a: 1.0, cl: 1.0 };
        assert_eq!(line.at(0.0), (0.0, 0.0));
        assert_eq!(line.at(1.0), (0.0, 0.0));
        // Abbott and von Doenhoff: the a = 1 line has y_c/c = 0.05516 cl at x = 0.5
        assert!((line.at(0.5).0 - 0.05516).abs() < 1e-4, "{}", line.at(0.5).0);
    }

    #[test]
    fn six_a_line_is_straight_aft_of_0_86() {
        let line = MeanLine::SixSeriesModified { cl: 0.4 };
        let (y1, s1) = line.at(0.9);
        let (y2, s2) = line.at(0.95);
        assert_eq!(s1, s2);
        assert!((y2 - (y1 + s1 * 0.05)).abs() < 1e-15);
        assert_eq!(line.at(1.0).0, 0.0);
    }
}
