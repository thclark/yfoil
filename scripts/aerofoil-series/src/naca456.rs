//! The generators against the naca456 reference fixtures: the worst relative difference per
//! case and column, for the report. The gates themselves are `tests/naca456_series_tests.rs`;
//! this tabulates what they assert.

use std::collections::HashMap;
use std::path::Path;
use yfoil::geometry::{MeanLine, Section, ThicknessForm};

pub struct CaseComparison {
    pub slug: String,
    pub designation: String,
    pub family: &'static str,
    pub stations: usize,
    /// Worst |Δ| / max(|a|, |b|, 1) over the stations, per column
    pub y_t: f64,
    pub y_t_slope: f64,
    pub y_c: f64,
    pub y_c_slope: f64,
    pub surface: f64,
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1.0)
}

/// Every case of the fixture manifest, in its order
pub fn compare(root: &Path) -> Vec<CaseComparison> {
    let dir = root.join("tests/fixtures/naca456");
    let manifest: serde_json::Value = match std::fs::read_to_string(dir.join("manifest.json")) {
        Ok(s) => serde_json::from_str(&s).unwrap(),
        Err(_) => return Vec::new(),
    };
    let mut out = Vec::new();
    for slug in manifest["cases"].as_array().unwrap() {
        let slug = slug.as_str().unwrap();
        let text = std::fs::read_to_string(dir.join(format!("{slug}.dat"))).unwrap();
        let mut header: HashMap<String, String> = HashMap::new();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for l in text.lines() {
            if let Some(h) = l.strip_prefix("# ") {
                if let Some((k, v)) = h.split_once(' ') {
                    header.insert(k.to_string(), v.trim().to_string());
                }
            } else if !l.trim().is_empty() {
                rows.push(l.split_whitespace().map(|t| t.parse().unwrap()).collect());
            }
        }
        let a: f64 = header["a"].parse().unwrap();
        let mut s = Section::from_designation(&header["name"]).unwrap();
        if matches!(s.mean_line, MeanLine::SixSeries { .. }) && a != 1.0 {
            s = s.with_a(a).unwrap();
        }
        let family = match s.thickness_form {
            ThicknessForm::FourDigit { .. } => "closed form",
            ThicknessForm::FourDigitModified { .. } => "closed form",
            ThicknessForm::SixSeries(_) => "6-series (tabulated map)",
        };
        let aft_4m = |x: f64| -> bool {
            matches!(s.thickness_form, ThicknessForm::FourDigitModified { x_max_thickness, .. } if x >= x_max_thickness)
        };
        let mut c = CaseComparison {
            slug: slug.to_string(),
            designation: s.designation.clone(),
            family,
            stations: rows.len(),
            y_t: 0.0,
            y_t_slope: 0.0,
            y_c: 0.0,
            y_c_slope: 0.0,
            surface: 0.0,
        };
        for r in &rows {
            let x = r[0];
            let (yt, ytp) = s.thickness_form.at(x);
            let (ym, ymp) = s.mean_line.at(x);
            c.y_t = c.y_t.max(rel(yt, r[1]));
            if x > 0.0 && r[2].abs() < 1e10 {
                // naca456 reports the aft slope of the 4-digit modified form as dy/d(1 − x)
                let reference = if aft_4m(x) { -r[2] } else { r[2] };
                c.y_t_slope = c.y_t_slope.max(rel(ytp, reference));
            }
            c.y_c = c.y_c.max(rel(ym, r[3]));
            c.y_c_slope = c.y_c_slope.max(rel(ymp, r[4]));
            let (xu, yu, xl, yl) = s.surface_at(x);
            for (v, w) in [(xu, r[5]), (yu, r[6]), (xl, r[7]), (yl, r[8])] {
                c.surface = c.surface.max(rel(v, w));
            }
        }
        out.push(c);
    }
    out
}
