//! Aerofoil series generators: every NACA family with a public, reproducible definition, and the
//! Kármán–Trefftz analytic section.
//!
//! A NACA section is a basic thickness form ([`ThicknessForm`]) laid perpendicular to a mean
//! line ([`MeanLine`]) — the NACA definition (Abbott and von Doenhoff 1959, §6.2), as
//! implemented by the NASA ordinate program (Ladson, Brooks, Hill and Sproles 1996, TM-4741) and
//! its public-domain revision `naca456` (Carmichael 2001), which is the reference the
//! `tests/fixtures/naca456/` fixtures come from. [`Section`] parses a designation into the pair,
//! evaluates the surface at any chord station, and panels it with the cosine-in-x spacing
//! shared by every generator, so that sections of different families are compared on the same
//! node distribution.
//!
//! | Series | Designation | Thickness form | Mean line |
//! |---|---|---|---|
//! | 4-digit | `2412` | 4-digit, t = 0.12 | 2-digit, m = 0.02, p = 0.4 |
//! | 4-digit modified | `0012-34`, `2412-63` | 4-digit modified, I = 3, M = 0.4 | 2-digit |
//! | 5-digit | `23012`, `23112` (reflex) | 4-digit | 3-digit, cl = 0.3, p = 0.15 |
//! | 16-series | `16-212` | 4-digit modified, I = 4, M = 0.5 | 6-series a = 1, cl = 0.2 |
//! | 6-series | `63-415`, `64-010` | 63…67 tables | 6-series a (default 1), cl = 0.4 |
//! | 6A-series | `64A010`, `63A415` | 63A…65A tables | 6A modified line |
//!
//! Every generated [`Geometry`] carries a `generator` record (JSON) stating the series, the
//! designation, the thickness form and mean line with their parameter values, how the thickness
//! was applied, and the references that define the family — the provenance of the coordinates.

mod karman_trefftz;
mod mean_line;
mod six_series;
mod six_series_tables;
mod spline_fmm;
mod thickness;

pub use karman_trefftz::{KarmanTrefftz, KarmanTrefftzError};
pub use mean_line::{three_digit_constants, three_digit_reflex_constants, two_digit_mean_line, MeanLine};
pub use six_series::{basic_thickness_form, scale_factor, SixSeriesFamily, SixSeriesForm};
pub use thickness::{four_digit_half_thickness, four_digit_modified_coefficients, four_digit_slope, ThicknessForm};

use super::airfoil::Geometry;
use super::naca::NacaError;
use super::panel::{
    apply_te_treatment, record_panelling, repanel, CosineConfig, PanelConfig, PanelMethod, RepanelError,
    PANGEN_BUFFER_NODES,
};
use serde_json::{json, Value};

/// The NACA families with a designation grammar
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Series {
    FourDigit,
    FourDigitModified,
    FiveDigit,
    Sixteen,
    Six,
    SixA,
}

impl Series {
    pub const ALL: [Series; 6] = [
        Self::FourDigit,
        Self::FourDigitModified,
        Self::FiveDigit,
        Self::Sixteen,
        Self::Six,
        Self::SixA,
    ];

    /// The `series` value of the provenance record
    pub fn slug(self) -> &'static str {
        match self {
            Self::FourDigit => "naca_4_digit",
            Self::FourDigitModified => "naca_4_digit_modified",
            Self::FiveDigit => "naca_5_digit",
            Self::Sixteen => "naca_16",
            Self::Six => "naca_6",
            Self::SixA => "naca_6a",
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Self::FourDigit => "NACA 4-digit",
            Self::FourDigitModified => "NACA 4-digit modified",
            Self::FiveDigit => "NACA 5-digit",
            Self::Sixteen => "NACA 16-series",
            Self::Six => "NACA 6-series",
            Self::SixA => "NACA 6A-series",
        }
    }

    /// Keys into `docs/references.bib` of the reports that define the family
    pub fn references(self) -> &'static [&'static str] {
        match self {
            Self::FourDigit => &["jacobs1933", "abbott1959", "ladson1996"],
            Self::FourDigitModified => &["stack1934", "abbott1959", "ladson1975", "ladson1996"],
            Self::FiveDigit => &["jacobs1935", "abbott1959", "ladson1996"],
            Self::Sixteen => &["stack1943", "lindsey1948", "ladson1975", "ladson1996"],
            Self::Six => &["abbott1945", "abbott1959", "ladson1974", "ladson1996", "carmichael2001"],
            Self::SixA => &["loftin1948", "ladson1974", "ladson1996", "carmichael2001"],
        }
    }
}

/// A NACA section: a thickness form on a mean line, with its designation
#[derive(Debug, Clone)]
pub struct Section {
    pub designation: String,
    pub series: Series,
    pub thickness_form: ThicknessForm,
    pub mean_line: MeanLine,
}

fn digit(s: &str, i: usize) -> Option<f64> {
    s.as_bytes()
        .get(i)
        .filter(|b| b.is_ascii_digit())
        .map(|b| (b - b'0') as f64)
}

fn all_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

fn invalid(msg: impl Into<String>) -> NacaError {
    NacaError::InvalidDesignation(msg.into())
}

impl Section {
    /// Parse a designation (an optional `NACA ` prefix is ignored; case is ignored):
    /// `0012`, `2412`, `0012-34`, `2412-63`, `23012`, `23112`, `16-212`, `63-415`, `63A415`.
    /// A 6-series section gets the `a = 1.0` mean line; use [`Section::with_a`] for another.
    pub fn from_designation(designation: &str) -> Result<Section, NacaError> {
        let s = designation.trim();
        let s = s
            .strip_prefix("NACA ")
            .or_else(|| s.strip_prefix("naca "))
            .unwrap_or(s)
            .trim();
        let upper = s.to_ascii_uppercase();
        let s = upper.as_str();

        // 4-digit and 4-digit modified: dddd or dddd-IM
        if let Some((front, back)) = s.split_once('-') {
            if front.len() == 4 && all_digits(front) {
                if back.len() != 2 || !all_digits(back) {
                    return Err(invalid(format!(
                        "{designation}: a 4-digit modified designation is dddd-IM (leading-edge index, position of maximum thickness in tenths)"
                    )));
                }
                let (m, p, t) = four_digit_parameters(front);
                let le_radius_index = digit(back, 0).unwrap();
                let x_max_thickness = digit(back, 1).unwrap() / 10.0;
                if x_max_thickness < 0.2 || x_max_thickness > 0.6 {
                    return Err(invalid(format!(
                        "{designation}: the position of maximum thickness must be 2–6 tenths of chord (Report 492 covers 0.2–0.6)"
                    )));
                }
                return Ok(Section {
                    designation: format!("NACA {front}-{back}"),
                    series: Series::FourDigitModified,
                    thickness_form: ThicknessForm::FourDigitModified {
                        t,
                        le_radius_index,
                        x_max_thickness,
                    },
                    mean_line: MeanLine::TwoDigit { m, p },
                });
            }
            // 16-series: 16-ctt
            if front == "16" {
                if back.len() != 3 || !all_digits(back) {
                    return Err(invalid(format!(
                        "{designation}: a 16-series designation is 16-ctt (design lift coefficient in tenths, thickness in percent)"
                    )));
                }
                let cl = digit(back, 0).unwrap() / 10.0;
                let t = back[1..3].parse::<f64>().unwrap() / 100.0;
                return Ok(Section {
                    designation: format!("NACA 16-{back}"),
                    series: Series::Sixteen,
                    thickness_form: ThicknessForm::FourDigitModified {
                        t,
                        le_radius_index: 4.0,
                        x_max_thickness: 0.5,
                    },
                    mean_line: MeanLine::SixSeries { a: 1.0, cl },
                });
            }
            // 6-series: 6f-ctt, with an optional bracketed low-drag-range subscript 6f(d)-ctt
            let family_text = front.split('(').next().unwrap_or(front);
            if let Some(family) = SixSeriesFamily::from_label(family_text) {
                if family.is_a_series() {
                    return Err(invalid(format!(
                        "{designation}: 6A-series designations have no hyphen (64A010)"
                    )));
                }
                if back.len() != 3 || !all_digits(back) {
                    return Err(invalid(format!(
                        "{designation}: a 6-series designation is 6f-ctt (family, design lift coefficient in tenths, thickness in percent)"
                    )));
                }
                let cl = digit(back, 0).unwrap() / 10.0;
                let t = back[1..3].parse::<f64>().unwrap() / 100.0;
                return Ok(Section {
                    designation: format!("NACA {}-{back}", family.label()),
                    series: Series::Six,
                    thickness_form: ThicknessForm::six_series(family, t),
                    mean_line: MeanLine::SixSeries { a: 1.0, cl },
                });
            }
            return Err(invalid(format!("{designation}: unrecognised designation")));
        }

        // 6A-series: 6fActt
        if let Some(pos) = s.find('A') {
            let (front, back) = (&s[..pos + 1], &s[pos + 1..]);
            let family = SixSeriesFamily::from_label(front)
                .filter(|f| f.is_a_series())
                .ok_or_else(|| invalid(format!("{designation}: unrecognised 6A-series family {front}")))?;
            if back.len() != 3 || !all_digits(back) {
                return Err(invalid(format!(
                    "{designation}: a 6A-series designation is 6fActt (family, design lift coefficient in tenths, thickness in percent)"
                )));
            }
            let cl = digit(back, 0).unwrap() / 10.0;
            let t = back[1..3].parse::<f64>().unwrap() / 100.0;
            return Ok(Section {
                designation: format!("NACA {}{back}", family.label()),
                series: Series::SixA,
                thickness_form: ThicknessForm::six_series(family, t),
                mean_line: MeanLine::SixSeriesModified { cl },
            });
        }

        if s.len() == 4 && all_digits(s) {
            let (m, p, t) = four_digit_parameters(s);
            return Ok(Section {
                designation: format!("NACA {s}"),
                series: Series::FourDigit,
                thickness_form: ThicknessForm::FourDigit { t },
                mean_line: MeanLine::TwoDigit { m, p },
            });
        }

        if s.len() == 5 && all_digits(s) {
            let cl = digit(s, 0).unwrap() * 3.0 / 20.0;
            let second = digit(s, 1).unwrap();
            let p = second / 20.0;
            let third = digit(s, 2).unwrap();
            let t = s[3..5].parse::<f64>().unwrap() / 100.0;
            if third > 1.0 {
                return Err(invalid(format!(
                    "{designation}: the third digit of a 5-digit designation is 0 (standard) or 1 (reflex)"
                )));
            }
            let reflex = third == 1.0;
            let mean_line = if cl == 0.0 {
                MeanLine::Symmetric
            } else if reflex {
                if !(2.0..=5.0).contains(&second) {
                    return Err(invalid(format!(
                        "{designation}: reflex mean lines are tabulated for the second digit 2–5"
                    )));
                }
                MeanLine::ThreeDigitReflex { cl, p }
            } else {
                if !(1.0..=5.0).contains(&second) {
                    return Err(invalid(format!(
                        "{designation}: 3-digit mean lines are tabulated for the second digit 1–5"
                    )));
                }
                MeanLine::ThreeDigit { cl, p }
            };
            return Ok(Section {
                designation: format!("NACA {s}"),
                series: Series::FiveDigit,
                thickness_form: ThicknessForm::FourDigit { t },
                mean_line,
            });
        }

        Err(invalid(format!(
            "{designation}: not a NACA 4-digit (2412), 4-digit modified (0012-34), 5-digit (23012), 16-series (16-212), 6-series (63-415) or 6A-series (64A010) designation"
        )))
    }

    /// Set the extent of uniform loading `a` of a 6-series mean line (the `a = 0.5` of
    /// "NACA 63-415, a = 0.5"). Only 6-series and 16-series sections carry one.
    pub fn with_a(mut self, a: f64) -> Result<Section, NacaError> {
        match self.mean_line {
            MeanLine::SixSeries { cl, .. } => {
                if !(0.0..=1.0).contains(&a) {
                    return Err(invalid(format!("a = {a}: the extent of uniform loading is 0 ≤ a ≤ 1")));
                }
                self.mean_line = MeanLine::SixSeries { a, cl };
                if a != 1.0 {
                    self.designation = format!("{}, a = {a}", self.designation);
                }
                Ok(self)
            }
            _ => Err(invalid(format!(
                "{}: only 6-series and 16-series mean lines have a loading extent a",
                self.designation
            ))),
        }
    }

    /// Whether the section's thickness form closes at the trailing edge
    pub fn sharp_te(&self) -> bool {
        self.thickness_form.sharp_te()
    }

    /// Upper and lower surface points at chord station `x`, the thickness laid perpendicular to
    /// the mean line: `(x_upper, y_upper, x_lower, y_lower)`
    pub fn surface_at(&self, x: f64) -> (f64, f64, f64, f64) {
        let (yt, _) = self.thickness_form.at(x);
        let (yc, dyc_dx) = self.mean_line.at(x);
        let theta = dyc_dx.atan();
        (
            x - yt * theta.sin(),
            yc + yt * theta.cos(),
            x + yt * theta.sin(),
            yc - yt * theta.cos(),
        )
    }

    /// The panelled section: `n_panels` nodes (even) at cosine stations in x, trailing edge →
    /// upper → leading edge (straddled) → lower → trailing edge, with the provenance record
    pub fn geometry(&self, n_panels: usize) -> Geometry {
        let stations = cosine_stations(n_panels);
        let n_half = stations.len();
        let mut x_upper = Vec::with_capacity(n_half);
        let mut y_upper = Vec::with_capacity(n_half);
        let mut x_lower = Vec::with_capacity(n_half);
        let mut y_lower = Vec::with_capacity(n_half);
        for &x in &stations {
            let (xu, yu, xl, yl) = self.surface_at(x);
            x_upper.push(xu);
            y_upper.push(yu);
            x_lower.push(xl);
            y_lower.push(yl);
        }
        let mut x = Vec::with_capacity(2 * n_half);
        let mut y = Vec::with_capacity(2 * n_half);
        for i in (0..n_half).rev() {
            x.push(x_upper[i]);
            y.push(y_upper[i]);
        }
        for i in 0..n_half {
            x.push(x_lower[i]);
            y.push(y_lower[i]);
        }
        Geometry {
            cm_ref: [0.25, 0.0], // Quarter chord
            x,
            y,
            generator: Some(self.record()),
        }
    }

    /// The section panelled as `config` says: `cosine` is the analytic sampling of
    /// [`Self::geometry`] at `n_nodes` chord stations (no bias applies; the record says so);
    /// `pangen` samples the section at `buffer_nodes` (default [`PANGEN_BUFFER_NODES`]) and runs
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
                let buffer_nodes = config.buffer_nodes.unwrap_or(PANGEN_BUFFER_NODES);
                let buffer = self.geometry(buffer_nodes);
                let mut on_buffer = *config;
                on_buffer.buffer_nodes = None;
                let mut out = repanel(&buffer, &on_buffer)?;
                let mut used = *config;
                used.buffer_nodes = Some(buffer_nodes);
                record_panelling(&mut out, &used);
                Ok(out)
            }
        }
    }

    /// The provenance record written as `generator` in the geometry JSON
    pub fn record(&self) -> Value {
        json!({
            "yfoil": env!("CARGO_PKG_VERSION"),
            "series": self.series.slug(),
            "designation": self.designation,
            "thickness_form": self.thickness_form.record(),
            "mean_line": self.mean_line.record(),
            "thickness_applied": "perpendicular",
            "sharp_te": self.sharp_te(),
            "references": self.series.references(),
        })
    }
}

/// `(m, p, t)` of a 4-digit designation `mptt`
fn four_digit_parameters(s: &str) -> (f64, f64, f64) {
    let m = digit(s, 0).unwrap() / 100.0;
    let p = digit(s, 1).unwrap() / 10.0;
    let t = s[2..4].parse::<f64>().unwrap() / 100.0;
    (m, p, t)
}

/// The chord stations of one surface for an `n_panels`-node section, leading edge first: the
/// cosine distribution `x = ½(1 − cos β)` with β = π(i + ½)/(n/2 − ½), which puts a node exactly
/// at the trailing edge and straddles the leading edge (no node at x = 0, as XFOIL's `PANE`
/// does; a node there has zero vortex strength and the boundary layer fails to start).
pub fn cosine_stations(n_panels: usize) -> Vec<f64> {
    let n_half = n_panels / 2;
    (0..n_half)
        .map(|i| {
            let beta = std::f64::consts::PI * (i as f64 + 0.5) / (n_half as f64 - 0.5);
            let beta = beta.min(std::f64::consts::PI); // Cap at π for last point
            0.5 * (1.0 - beta.cos())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designations_parse_to_the_documented_parameters() {
        let s = Section::from_designation("2412").unwrap();
        assert_eq!(s.series, Series::FourDigit);
        assert_eq!(s.mean_line, MeanLine::TwoDigit { m: 0.02, p: 0.4 });
        assert_eq!(s.thickness_form.t(), 0.12);

        let s = Section::from_designation("naca 0012-34").unwrap();
        assert_eq!(s.series, Series::FourDigitModified);
        assert!(matches!(
            s.thickness_form,
            ThicknessForm::FourDigitModified { t, le_radius_index, x_max_thickness }
                if t == 0.12 && le_radius_index == 3.0 && x_max_thickness == 0.4
        ));
        assert!(s.mean_line.is_symmetric());

        let s = Section::from_designation("23012").unwrap();
        assert_eq!(s.series, Series::FiveDigit);
        assert_eq!(s.mean_line, MeanLine::ThreeDigit { cl: 0.3, p: 0.15 });
        let s = Section::from_designation("23112").unwrap();
        assert_eq!(s.mean_line, MeanLine::ThreeDigitReflex { cl: 0.3, p: 0.15 });

        let s = Section::from_designation("16-212").unwrap();
        assert_eq!(s.series, Series::Sixteen);
        assert_eq!(s.mean_line, MeanLine::SixSeries { a: 1.0, cl: 0.2 });
        assert!(matches!(
            s.thickness_form,
            ThicknessForm::FourDigitModified { t, le_radius_index, x_max_thickness }
                if t == 0.12 && le_radius_index == 4.0 && x_max_thickness == 0.5
        ));

        let s = Section::from_designation("63-415").unwrap().with_a(0.5).unwrap();
        assert_eq!(s.series, Series::Six);
        assert_eq!(s.designation, "NACA 63-415, a = 0.5");
        assert_eq!(s.mean_line, MeanLine::SixSeries { a: 0.5, cl: 0.4 });
        let s = Section::from_designation("64(1)-212").unwrap();
        assert_eq!(s.designation, "NACA 64-212");

        let s = Section::from_designation("64a010").unwrap();
        assert_eq!(s.series, Series::SixA);
        assert_eq!(s.designation, "NACA 64A010");
        assert!(s.mean_line.is_symmetric());
        assert!(s.sharp_te());
    }

    #[test]
    fn bad_designations_are_rejected() {
        for bad in [
            "001", "0012-3", "0012-17", "26012", "23212", "16-12", "68-415", "63A-415", "xyz",
        ] {
            assert!(Section::from_designation(bad).is_err(), "{bad} should be rejected");
        }
        assert!(Section::from_designation("2412").unwrap().with_a(0.5).is_err());
        assert!(Section::from_designation("63-415").unwrap().with_a(1.5).is_err());
    }

    #[test]
    fn sixteen_series_thickness_is_the_0012_45_form() {
        let a = Section::from_designation("16-012").unwrap();
        let b = Section::from_designation("0012-45").unwrap();
        for &x in &[0.01, 0.1, 0.5, 0.9] {
            assert_eq!(a.thickness_form.at(x), b.thickness_form.at(x));
        }
    }

    #[test]
    fn panels_have_the_requested_count_te_at_one_and_a_straddled_le() {
        for d in ["0012", "4412", "0012-34", "23018", "16-212", "63-415", "64A010"] {
            let g = Section::from_designation(d).unwrap().geometry(160);
            assert_eq!(g.x.len(), 160, "{d}");
            // the TE nodes sit at x = 1 − y_t sin θ: exactly 1 for a symmetric or closed section
            assert!((g.x[0] - 1.0).abs() < 1e-3, "{d}: first node at the TE, x = {}", g.x[0]);
            assert!(
                (g.x[159] - 1.0).abs() < 1e-3,
                "{d}: last node at the TE, x = {}",
                g.x[159]
            );
            // the LE is straddled: the two middle nodes sit within a panel of x = 0 (a cambered
            // section's upper node can be at slightly negative x, offset along the mean-line normal)
            assert!(g.x[79] < 1e-3 && g.x[80] < 1e-3, "{d}: LE straddled");
            g.validate().unwrap_or_else(|e| panic!("{d}: {e}"));
            let rec = g.generator.as_ref().unwrap();
            assert_eq!(
                rec["designation"].as_str().unwrap(),
                format!("NACA {}", d.to_uppercase())
            );
        }
    }

    #[test]
    fn symmetric_sections_are_mirror_images() {
        for d in ["0012", "0012-34", "00012", "16-012", "63-012", "64A010"] {
            let g = Section::from_designation(d).unwrap().geometry(120);
            for i in 0..60 {
                assert_eq!(g.x[i], g.x[119 - i], "{d}: x mirror at {i}");
                assert_eq!(g.y[i], -g.y[119 - i], "{d}: y mirror at {i}");
            }
        }
    }
}
