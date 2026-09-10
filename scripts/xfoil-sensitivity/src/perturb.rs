//! Input perturbations.
//!
//! One fixed pattern `u` of uniform deviates in [−1, 1) — one per coordinate, from a seeded
//! SplitMix64 generator — is shared by every geometry level, so the size of the perturbation is
//! the only thing that changes between levels:
//!
//! - the 1-ULP level moves every coordinate to its neighbouring double in the direction of
//!   sign(u), i.e. by one unit in the last place of that coordinate (relative, ~1.1e-16 for x in
//!   [0.5, 1), ~1e-19 for a y of 1e-3);
//! - the scaled levels add ε·u to every coordinate (absolute, chord units).

use yfoil::geometry::Geometry;

/// Marker value meaning "one ULP" in the level tables
pub const ULP: f64 = 0.0;

/// SplitMix64 (Steele, Lea & Flood 2014): a 64-bit state, one output per step
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    /// Uniform in [−1, 1): the top 53 bits as a double in [0, 1), affinely mapped
    pub fn uniform_pm1(&mut self) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * unit - 1.0
    }
}

/// The perturbation pattern for `n` nodes: `(u_x, u_y)` per node, in node order
pub fn pattern(seed: u64, n: usize) -> Vec<(f64, f64)> {
    let mut g = SplitMix64::new(seed);
    (0..n).map(|_| (g.uniform_pm1(), g.uniform_pm1())).collect()
}

/// The next representable double above `x`
pub fn next_up(x: f64) -> f64 {
    step_ulp(x, 1.0)
}

/// The neighbouring double of `x` in the direction of the sign of `dir` (one ULP of `x`; for
/// `x == 0` the smallest subnormal with that sign)
pub fn step_ulp(x: f64, dir: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        return x;
    }
    if x == 0.0 {
        return if dir >= 0.0 {
            f64::from_bits(1)
        } else {
            -f64::from_bits(1)
        };
    }
    let away = (x > 0.0) == (dir >= 0.0);
    let bits = x.to_bits();
    f64::from_bits(if away { bits + 1 } else { bits - 1 })
}

/// Apply level `eps` of the geometry family with pattern `u`
pub fn perturb_geometry(g: &Geometry, u: &[(f64, f64)], eps: f64) -> Geometry {
    assert_eq!(g.x.len(), u.len(), "pattern length");
    let mut out = g.clone();
    for i in 0..g.x.len() {
        if eps == ULP {
            out.x[i] = step_ulp(g.x[i], u[i].0);
            out.y[i] = step_ulp(g.y[i], u[i].1);
        } else {
            out.x[i] = g.x[i] + eps * u[i].0;
            out.y[i] = g.y[i] + eps * u[i].1;
        }
    }
    out
}

/// Directory slug of a level: `ulp`, `1e-15`, …
pub fn eps_slug(eps: f64) -> String {
    if eps == ULP {
        "ulp".to_string()
    } else {
        format!("{eps:.0e}")
    }
}

/// Legend text of a level: `1 ULP`, `1e-15`, …
pub fn eps_label(eps: f64) -> String {
    if eps == ULP {
        "1 ULP".to_string()
    } else {
        format!("{eps:.0e}")
    }
}
