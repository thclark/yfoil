//! The numbers the tables and figures show, computed from a run's `summary.json` alone.

use crate::drive::val;
use serde_json::Value;

/// The thresholds of the plateau: half-decades from 1e-14 to 1.
pub fn plateau_thresholds() -> Vec<f64> {
    (0..=28).map(|i| cr_libm::pow(10.0, -14.0 + 0.5 * i as f64)).collect()
}

/// The decades of the histogram, `[10^k, 10^(k+1))` for k = −17 … −1; spreads below 1e-17
/// (zero included) fall in the first, spreads of one or more in the last.
pub const DECADES: std::ops::RangeInclusive<i32> = -17..=-1;

/// A spread is flagged at threshold `t` when it exceeds it or is not a number.
pub fn flagged(s: f64, t: f64) -> bool {
    s.is_nan() || s > t
}

/// A point was measured: valid on single-run evidence, so its twins were compared.
pub fn measured(p: &Value) -> bool {
    p["measured"].as_bool() == Some(true)
}

/// A measured point is ill-conditioned: its largest spread is beyond the error threshold.
pub fn ill_conditioned(p: &Value, summary: &Value) -> bool {
    measured(p) && flagged(val(&p["worst_spread"]), thresholds(summary).error)
}

/// A measured point may be ill-conditioned: its largest spread is beyond the warning threshold
/// and within the error threshold.
pub fn may_be_ill_conditioned(p: &Value, summary: &Value) -> bool {
    let t = thresholds(summary);
    measured(p) && !ill_conditioned(p, summary) && flagged(val(&p["worst_spread"]), t.warning)
}

/// The thresholds of the run.
pub fn thresholds(summary: &Value) -> yfoil::output::Thresholds {
    serde_json::from_value(summary["thresholds"].clone()).expect("summary.json carries the thresholds")
}

/// Every point of every case, in order.
pub fn points(summary: &Value) -> impl Iterator<Item = &Value> {
    points_of(summary, None)
}

/// Every point of the cases of `group` (of every case for `None`), in order.
pub fn points_of<'a>(summary: &'a Value, group: Option<&'a str>) -> impl Iterator<Item = &'a Value> {
    summary["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(move |c| group.is_none_or(|g| c["group"].as_str() == Some(g)))
        .flat_map(|c| c["points"].as_array().unwrap().iter())
}

/// The groups of the run, in order.
pub fn groups(summary: &Value) -> Vec<String> {
    summary["groups"]
        .as_array()
        .map(|a| a.iter().filter_map(|g| g.as_str().map(String::from)).collect())
        .unwrap_or_default()
}

/// The worst spread of every measured point: valid on single-run evidence, so its twins were
/// compared.
pub fn worst_spreads(summary: &Value) -> Vec<f64> {
    worst_spreads_of(summary, None)
}

/// The same over the cases of one group.
pub fn worst_spreads_of(summary: &Value, group: Option<&str>) -> Vec<f64> {
    points_of(summary, group)
        .filter(|p| measured(p))
        .map(|p| val(&p["worst_spread"]))
        .collect()
}

/// The decade index of a spread: `None` for NaN, 0 for the first decade, the last for one or more.
pub fn decade(s: f64) -> Option<usize> {
    if s.is_nan() {
        return None;
    }
    let lo = *DECADES.start();
    let k = if s <= 0.0 {
        lo
    } else {
        (cr_libm::log10(s).floor() as i32).clamp(lo, *DECADES.end() + 1)
    };
    Some((k - lo) as usize)
}

/// Counts per decade (`DECADES` plus the "≥ 1" bin) and of NaN.
pub fn histogram(spreads: &[f64]) -> (Vec<usize>, usize) {
    let n = DECADES.count() + 1;
    let mut bins = vec![0; n];
    let mut nan = 0;
    for &s in spreads {
        match decade(s) {
            Some(i) => bins[i.min(n - 1)] += 1,
            None => nan += 1,
        }
    }
    (bins, nan)
}

/// The gap around threshold `t`: the largest finite spread at or below it, the smallest above
/// it (`None` when there is none), and how many are above (NaN included).
pub struct Gap {
    pub below: Option<f64>,
    pub above: Option<f64>,
    pub n_above: usize,
    pub n_nan: usize,
    pub n: usize,
}

pub fn gap(spreads: impl Iterator<Item = f64>, t: f64) -> Gap {
    let mut g = Gap {
        below: None,
        above: None,
        n_above: 0,
        n_nan: 0,
        n: 0,
    };
    for s in spreads {
        g.n += 1;
        if s.is_nan() {
            g.n_nan += 1;
            g.n_above += 1;
        } else if s > t {
            g.n_above += 1;
            g.above = Some(g.above.map_or(s, |a: f64| a.min(s)));
        } else {
            g.below = Some(g.below.map_or(s, |b: f64| b.max(s)));
        }
    }
    g
}

impl Gap {
    /// The interval of thresholds over which the flagged set is the one at `t`: from the largest
    /// spread below (inclusive) to the smallest above (exclusive).
    pub fn interval(&self) -> (f64, f64) {
        (self.below.unwrap_or(0.0), self.above.unwrap_or(f64::INFINITY))
    }

    /// How many decades the threshold sits above the largest spread below it and below the
    /// smallest spread above it (infinite when there is none on that side).
    pub fn clearance(&self, t: f64) -> (f64, f64) {
        let (lo, hi) = self.interval();
        let decades = |a: f64, b: f64| {
            if a > 0.0 && b.is_finite() {
                cr_libm::log10(b / a)
            } else {
                f64::INFINITY
            }
        };
        (decades(lo, t), decades(t, hi))
    }
}

/// Decades to one decimal place, `∞` when unbounded.
pub fn decades(v: f64) -> String {
    if v.is_infinite() {
        "∞".into()
    } else {
        format!("{v:.1}")
    }
}

/// `1.2e-3`, `0`, `∞`, `NaN`.
pub fn sci(v: f64) -> String {
    if v.is_nan() {
        "NaN".into()
    } else if v == 0.0 {
        "0".into()
    } else if v.is_infinite() {
        "∞".into()
    } else {
        format!("{v:.1e}")
    }
}
