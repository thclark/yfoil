//! The Kármán–Trefftz section (von Kármán and Trefftz 1918), the Joukowski section (Joukowski
//! 1910) being the cusped member.
//!
//! A circle in the ζ-plane through ζ = 1 with centre `(x_centre, y_centre)` is mapped by
//!
//! ```text
//! z = n · [1 + w^n] / [1 − w^n],   w = (ζ − 1)/(ζ + 1),   n = 2 − τ/π
//! ```
//!
//! where τ is the trailing-edge angle: τ = 0 gives the Joukowski map z = ζ + 1/ζ and its cusp.
//! `x_centre < 0` sets the thickness, `y_centre > 0` the camber. The image is normalised so
//! that the trailing edge (the image of ζ = 1, which is z = n) sits at x = 1 and the leading edge
//! (the point of minimum x) at x = 0; the trailing edge has y = 0.
//!
//! The section is analytic: its potential-flow solution is known exactly, which is what makes it
//! a solver-independent check (CLAUDE.md Rule 6). The trailing edge is sharp, so it exercises the
//! `SHARP` branch of the panel method.

use super::super::airfoil::Geometry;
use super::super::panel::{
    apply_te_treatment, record_panelling, repanel, CosineConfig, PanelConfig, PanelMethod, RepanelError,
    PANGEN_BUFFER_NODES,
};
use super::cosine_stations;
use serde_json::{json, Value};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KarmanTrefftz {
    /// Circle centre in the ζ-plane; the circle passes through ζ = 1
    pub x_centre: f64,
    pub y_centre: f64,
    /// Trailing-edge angle τ in degrees, 0 ≤ τ < 180 (0 is the Joukowski cusp)
    pub te_angle_deg: f64,
}

#[derive(thiserror::Error, Debug)]
pub enum KarmanTrefftzError {
    #[error("trailing-edge angle must be in [0, 180) degrees, got {0}")]
    TrailingEdgeAngle(f64),
    #[error("the circle through ζ = 1 with centre ({0}, {1}) does not enclose ζ = −1")]
    CentreOutside(f64, f64),
}

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self {
            re: (self.re * o.re + self.im * o.im) / d,
            im: (self.im * o.re - self.re * o.im) / d,
        }
    }
    /// Principal-branch power `self^n`
    fn powf(self, n: f64) -> Self {
        let r = self.re.hypot(self.im);
        if r == 0.0 {
            return Self { re: 0.0, im: 0.0 };
        }
        let theta = self.im.atan2(self.re);
        let rn = r.powf(n);
        let (s, c) = (n * theta).sin_cos();
        Self { re: rn * c, im: rn * s }
    }
}

impl KarmanTrefftz {
    pub fn new(x_centre: f64, y_centre: f64, te_angle_deg: f64) -> Result<Self, KarmanTrefftzError> {
        if !(0.0..180.0).contains(&te_angle_deg) {
            return Err(KarmanTrefftzError::TrailingEdgeAngle(te_angle_deg));
        }
        let radius = (1.0 - x_centre).hypot(y_centre);
        if (-1.0 - x_centre).hypot(y_centre) >= radius {
            return Err(KarmanTrefftzError::CentreOutside(x_centre, y_centre));
        }
        Ok(Self {
            x_centre,
            y_centre,
            te_angle_deg,
        })
    }

    /// The map exponent n = 2 − τ/π
    pub fn exponent(&self) -> f64 {
        2.0 - self.te_angle_deg.to_radians() / PI
    }

    fn radius(&self) -> f64 {
        (1.0 - self.x_centre).hypot(self.y_centre)
    }

    /// Circle angle of the trailing edge ζ = 1
    fn theta_te(&self) -> f64 {
        (-self.y_centre).atan2(1.0 - self.x_centre)
    }

    /// The image of the circle point at angle θ, before normalisation
    fn map(&self, theta: f64) -> (f64, f64) {
        let r = self.radius();
        let (s, c) = theta.sin_cos();
        let zeta = Complex {
            re: self.x_centre + r * c,
            im: self.y_centre + r * s,
        };
        let w = Complex {
            re: zeta.re - 1.0,
            im: zeta.im,
        }
        .div(Complex {
            re: zeta.re + 1.0,
            im: zeta.im,
        });
        let n = self.exponent();
        let wn = w.powf(n);
        let z = Complex {
            re: 1.0 + wn.re,
            im: wn.im,
        }
        .div(Complex {
            re: 1.0 - wn.re,
            im: -wn.im,
        });
        (n * z.re, n * z.im)
    }

    /// Circle angle of the leading edge (minimum x of the image), found by sampling and a
    /// golden-section refinement, and the minimum x itself
    fn leading_edge(&self) -> (f64, f64) {
        let t0 = self.theta_te();
        let samples = 720;
        let mut best = (t0 + PI, f64::INFINITY);
        for k in 1..samples {
            let th = t0 + 2.0 * PI * k as f64 / samples as f64;
            let x = self.map(th).0;
            if x < best.1 {
                best = (th, x);
            }
        }
        let step = 2.0 * PI / samples as f64;
        let (mut a, mut b) = (best.0 - step, best.0 + step);
        let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
        let (mut c, mut d) = (b - phi * (b - a), a + phi * (b - a));
        let (mut fc, mut fd) = (self.map(c).0, self.map(d).0);
        for _ in 0..200 {
            if fc < fd {
                b = d;
                d = c;
                fd = fc;
                c = b - phi * (b - a);
                fc = self.map(c).0;
            } else {
                a = c;
                c = d;
                fc = fd;
                d = a + phi * (b - a);
                fd = self.map(d).0;
            }
            if b - a < 1e-15 {
                break;
            }
        }
        let th = 0.5 * (a + b);
        (th, self.map(th).0)
    }

    /// Chord-normalised point at circle angle θ
    fn point(&self, theta: f64, x_le: f64, chord: f64) -> (f64, f64) {
        let (x, y) = self.map(theta);
        ((x - x_le) / chord, y / chord)
    }

    /// Upper and lower surface points at the chord stations `x` (each in [0, 1]):
    /// `(x_upper, y_upper, x_lower, y_lower)` per station
    pub fn surface_at(&self, stations: &[f64]) -> Vec<(f64, f64, f64, f64)> {
        let t_te = self.theta_te();
        let (t_le, x_le) = self.leading_edge();
        let chord = self.exponent() - x_le;
        // x(θ) along a surface is monotone between the trailing and leading edges; bisect
        let solve = |x_target: f64, mut a: f64, mut b: f64| -> f64 {
            // x decreases from a to b on the upper surface (a = TE), increases on the lower
            if x_target >= 1.0 {
                return if self.point(a, x_le, chord).0 >= self.point(b, x_le, chord).0 {
                    a
                } else {
                    b
                };
            }
            let fa = self.point(a, x_le, chord).0 - x_target;
            for _ in 0..200 {
                let m = 0.5 * (a + b);
                let fm = self.point(m, x_le, chord).0 - x_target;
                if fm == 0.0 || (b - a).abs() <= f64::EPSILON * m.abs() {
                    return m;
                }
                if (fm > 0.0) == (fa > 0.0) {
                    a = m;
                } else {
                    b = m;
                }
            }
            0.5 * (a + b)
        };
        stations
            .iter()
            .map(|&x| {
                if x >= 1.0 {
                    // the image of ζ = 1 is z = n exactly, and both surfaces end there
                    return (1.0, 0.0, 1.0, 0.0);
                }
                let tu = solve(x, t_te, t_le);
                let tl = solve(x, t_te + 2.0 * PI, t_le);
                let (xu, yu) = self.point(tu, x_le, chord);
                let (xl, yl) = self.point(tl, x_le, chord);
                (xu, yu, xl, yl)
            })
            .collect()
    }

    /// The panelled section: cosine stations in x, trailing edge → upper → leading edge →
    /// lower → trailing edge, with the provenance record attached
    pub fn geometry(&self, n_nodes: usize) -> Geometry {
        let stations = cosine_stations(n_nodes);
        let pts = self.surface_at(&stations);
        let n = pts.len();
        let mut x = Vec::with_capacity(2 * n);
        let mut y = Vec::with_capacity(2 * n);
        for p in pts.iter().rev() {
            x.push(p.0);
            y.push(p.1);
        }
        for p in &pts {
            x.push(p.2);
            y.push(p.3);
        }
        Geometry {
            cm_ref: [0.25, 0.0],
            x,
            y,
            generator: Some(self.record()),
        }
    }

    /// The section panelled as `config` says: `cosine` is the analytic sampling of
    /// [`Self::geometry`] at `n_nodes` chord stations (no bias applies; the record says so);
    /// `pangen` samples the section at `n_buffer_nodes` (default [`PANGEN_BUFFER_NODES`]) and runs
    /// XFOIL's PANGEN on that buffer. The trailing-edge treatment follows the distribution, and
    /// the record gains `panelling`.
    pub fn panelled(&self, config: &PanelConfig) -> Result<Geometry, RepanelError> {
        config.validate()?;
        match config.method {
            PanelMethod::Cosine(_) => {
                let mut used = *config;
                used.method = PanelMethod::Cosine(CosineConfig { te_bias: None });
                let mut out = apply_te_treatment(self.geometry(config.n_nodes), &used);
                out.validate()?;
                record_panelling(&mut out, &used);
                Ok(out)
            }
            PanelMethod::Pangen(_) => {
                let n_buffer_nodes = config.n_buffer_nodes.unwrap_or(PANGEN_BUFFER_NODES);
                let buffer = self.geometry(n_buffer_nodes);
                let mut on_buffer = *config;
                on_buffer.n_buffer_nodes = None;
                let mut out = repanel(&buffer, &on_buffer)?;
                let mut used = *config;
                used.n_buffer_nodes = Some(n_buffer_nodes);
                record_panelling(&mut out, &used);
                Ok(out)
            }
        }
    }

    /// Designation text
    pub fn designation(&self) -> String {
        format!(
            "Karman-Trefftz xc={} yc={} tau={}",
            self.x_centre, self.y_centre, self.te_angle_deg
        )
    }

    /// The provenance record (`Geometry.generator`)
    pub fn record(&self) -> Value {
        json!({
            "yfoil": env!("CARGO_PKG_VERSION"),
            "series": "karman_trefftz",
            "designation": self.designation(),
            "x_centre": self.x_centre,
            "y_centre": self.y_centre,
            "te_angle_deg": self.te_angle_deg,
            "exponent": self.exponent(),
            "sharp_te": true,
            "references": ["karman1918", "joukowski1910"],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn joukowski_symmetric_section_is_symmetric_and_closed() {
        let kt = KarmanTrefftz::new(-0.1, 0.0, 0.0).unwrap();
        let g = kt.geometry(100);
        assert_eq!(g.x.len(), 100);
        assert_eq!((g.x[0], g.y[0]), (1.0, 0.0));
        assert_eq!((g.x[99], g.y[99]), (1.0, 0.0));
        for i in 0..50 {
            assert!((g.x[i] - g.x[99 - i]).abs() < 1e-14);
            assert!((g.y[i] + g.y[99 - i]).abs() < 1e-14);
        }
        assert!(g.x[49] > 0.0 && g.x[49] < 1e-3, "leading edge straddled: {}", g.x[49]);
    }

    #[test]
    fn cambered_section_has_the_requested_te_angle() {
        let kt = KarmanTrefftz::new(-0.1, 0.05, 20.0).unwrap();
        let pts = kt.surface_at(&[0.999, 0.9999]);
        // the surfaces near the TE meet at the TE angle: slopes of the two surfaces differ by τ
        let slope_u = (0.0 - pts[1].1) / (1.0 - pts[1].0);
        let slope_l = (0.0 - pts[1].3) / (1.0 - pts[1].2);
        let angle = (slope_u.atan() - slope_l.atan()).abs().to_degrees();
        assert!((angle - 20.0).abs() < 1.0, "TE angle {angle}");
        // camber: the upper surface is further from the chord line than the lower
        let mid = kt.surface_at(&[0.4])[0];
        assert!(mid.1 > -mid.3);
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(KarmanTrefftz::new(-0.1, 0.0, 180.0).is_err());
        assert!(KarmanTrefftz::new(0.2, 0.0, 10.0).is_err());
    }
}
