//! The conditioning of a point, measured from its twins (`crate::solver::twins`): Class C of
//! `docs/user-guide/reference/solution-validity.md`.
//!
//! The **spread** of a reported quantity is the furthest any twin moved it, in the project's one
//! error metric: `max over twins of |a − b| / max(|a|, |b|, 1)` — absolute for the coefficients,
//! which are all of order one or less, relative for anything larger. A point is
//! **ill-conditioned** when the spread of any quantity exceeds the error threshold, or is not a
//! number (a twin went through a NaN where the base did not), and **may be ill-conditioned** when
//! the largest spread exceeds only the warning threshold. Only values count: a twin that took
//! more iterations, or did not converge, is compared by the values it ended with, and is counted
//! in `twins_not_valid` as a warning rather than a reason.
//!
//! The default thresholds are set from the conditioning study (`docs/validation/conditioning/`),
//! which measures the spread of every point of the series and branch-coverage cases. Spreads fall
//! in two populations: the noise floor, where every twin reaches the same solution, and twins
//! that reach another solution. The error threshold, [`CONDITIONING_ERROR_THRESHOLD`], sits in
//! the gap between them, where which points are flagged does not change. The warning threshold,
//! [`CONDITIONING_WARNING_THRESHOLD`], sits inside the upper reaches of the noise floor: a point
//! above it reached the same solution, but some twin's convergence test tipped on the way.

use serde::{Deserialize, Serialize};

use crate::solver::analysis::PointResult;

/// A point whose twins move a reported quantity further than this is ill-conditioned, and
/// withheld. Set in the gap between the noise floor and the spread of a twin that reached another
/// solution; the derivation is `docs/validation/conditioning/index.md`.
pub const CONDITIONING_ERROR_THRESHOLD: f64 = 1e-4;

/// A valid point whose twins move a reported quantity further than this may be ill-conditioned,
/// and carries a warning. Set inside the upper reaches of the noise floor
/// (`docs/validation/conditioning/index.md`).
pub const CONDITIONING_WARNING_THRESHOLD: f64 = 1e-8;

/// The two spreads a point's twins are judged against (`--spread-warning-threshold`,
/// `--spread-error-threshold`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Thresholds {
    /// Beyond this, a valid point carries a `Warning::IllConditioned`
    pub warning: f64,
    /// Beyond this, the point is `Invalid` with `Reason::IllConditioned`
    pub error: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            warning: CONDITIONING_WARNING_THRESHOLD,
            error: CONDITIONING_ERROR_THRESHOLD,
        }
    }
}

/// The spread of each reported quantity over the twins. A quantity the point does not report
/// (the viscous ones of an inviscid point) is absent; one a twin made non-finite is `null`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Spread {
    /// Angle of attack (radians): the unknown of a fixed-CL point, zero spread otherwise
    pub alpha: Option<f64>,
    pub cl: Option<f64>,
    pub cm: Option<f64>,
    pub cd_pressure: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cd: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cd_friction: Option<f64>,
    /// Transition x/c on the upper side
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transition_upper: Option<f64>,
    /// Transition x/c on the lower side
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transition_lower: Option<f64>,
}

impl Spread {
    /// (name, spread) of every quantity the point reports; a non-finite spread is `NaN`
    pub fn entries(&self) -> Vec<(&'static str, f64)> {
        [
            ("alpha", self.alpha),
            ("cl", self.cl),
            ("cm", self.cm),
            ("cd_pressure", self.cd_pressure),
            ("cd", self.cd),
            ("cd_friction", self.cd_friction),
            ("transition_upper", self.transition_upper),
            ("transition_lower", self.transition_lower),
        ]
        .into_iter()
        .filter_map(|(name, v)| v.map(|v| (name, v)))
        .collect()
    }
}

/// The conditioning evidence of one point: how far its twins moved what it reports.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Conditioning {
    /// The twins' seeds (`crate::geometry::jitter`)
    pub seeds: Vec<u64>,
    /// The thresholds the spreads were judged against
    pub thresholds: Thresholds,
    /// The furthest any twin moved each reported quantity
    pub spread: Spread,
    /// The quantity with the largest spread (a non-finite one first)
    pub worst: String,
    /// Twins whose own solve is invalid (did not converge, left the Kármán–Tsien domain, …)
    /// although the point itself is valid. A warning, not a reason: only values decide.
    pub twins_not_valid: usize,
}

impl Conditioning {
    /// Measure the point's conditioning from its twins; `None` when it has none.
    pub fn measure(p: &PointResult, viscous: bool, thresholds: Thresholds) -> Option<Self> {
        let twins = p.twins.as_ref()?;
        let spread_of = |f: &dyn Fn(&PointResult) -> f64| -> f64 {
            let a = f(p);
            twins.iter().fold(0.0_f64, |worst, (_, t)| {
                let b = f(t);
                let s = (a - b).abs() / a.abs().max(b.abs()).max(1.0);
                // a NaN wins: it is the largest possible disagreement
                if s.is_nan() || worst.is_nan() {
                    f64::NAN
                } else {
                    worst.max(s)
                }
            })
        };
        let spread = Spread {
            alpha: Some(spread_of(&|r| r.alpha)),
            cl: Some(spread_of(&|r| r.cl)),
            cm: Some(spread_of(&|r| r.cm)),
            cd_pressure: Some(spread_of(&|r| r.cd_pressure)),
            cd: viscous.then(|| spread_of(&|r| r.cd)),
            cd_friction: viscous.then(|| spread_of(&|r| r.cd_friction)),
            transition_upper: viscous.then(|| spread_of(&|r| r.transition_upper[0])),
            transition_lower: viscous.then(|| spread_of(&|r| r.transition_lower[0])),
        };
        let worst = spread
            .entries()
            .into_iter()
            .fold(None::<(&str, f64)>, |w, (name, s)| match w {
                Some((_, ws)) if ws.is_nan() || (!s.is_nan() && s <= ws) => w,
                _ => Some((name, s)),
            })
            .map(|(name, _)| name.to_string())
            .unwrap_or_default();
        let twins_not_valid = twins
            .iter()
            .filter(|(_, t)| !super::results::reasons(t, viscous).is_empty())
            .count();
        Some(Self {
            seeds: twins.iter().map(|(s, _)| *s).collect(),
            thresholds,
            spread,
            worst,
            twins_not_valid,
        })
    }

    /// The spread of the worst quantity
    pub fn worst_spread(&self) -> f64 {
        self.spread
            .entries()
            .into_iter()
            .find(|(name, _)| *name == self.worst)
            .map(|(_, s)| s)
            .unwrap_or(0.0)
    }

    /// Some twin moved a reported quantity beyond the error threshold, or made it non-finite
    pub fn ill_conditioned(&self) -> bool {
        self.spread
            .entries()
            .iter()
            .any(|(_, s)| s.is_nan() || *s > self.thresholds.error)
    }

    /// Not ill-conditioned, but some twin moved a reported quantity beyond the warning threshold
    pub fn may_be_ill_conditioned(&self) -> bool {
        !self.ill_conditioned() && self.worst_spread() > self.thresholds.warning
    }
}
