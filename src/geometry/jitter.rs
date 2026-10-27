//! The one-ULP jitter that makes a *twin* (`docs/development/conventions/terminology.md`): every
//! coordinate moved by −1, 0 or +1 ULP, x and y independently, drawn reproducibly from a seed.
//!
//! This is the one definition. `cargo xtask twins` applies it to the reference's `panels.dat`,
//! the fixture tests to a case's panels, and the solver's own twins (`crate::solver::twins`) to
//! the analysed geometry, so twin `k` of yFoil and twin `k` of the reference see the same panels.
//! Moving every coordinate the same way is not a twin: it translates the aerofoil and measures
//! almost nothing.

use super::Geometry;

/// The seeds of the five twins every validation run and every yFoil analysis carries.
pub const TWIN_SEEDS: [u64; 5] = [1, 2, 3, 4, 5];

/// The jitter of item `i` under `seed`: −1, 0 or +1 ULP for each of a pair of values (a node's x
/// and y), drawn from splitmix64 of (seed, i).
pub fn jitter(seed: u64, i: usize) -> (i8, i8) {
    let mut z = (seed << 32) ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let three = |v: u64| -> i8 { (v % 3) as i8 - 1 };
    (three(z & 0xFFFF_FFFF), three(z >> 32))
}

/// `v` moved by `by` ULP (−1, 0 or +1).
pub fn step_ulp(v: f64, by: i8) -> f64 {
    match by {
        1 => v.next_up(),
        -1 => v.next_down(),
        _ => v,
    }
}

/// The geometry with node `i`'s coordinates moved by [`jitter`]`(seed, i)`; everything else kept.
pub fn jitter_geometry(geometry: &Geometry, seed: u64) -> Geometry {
    let mut g = geometry.clone();
    for i in 0..g.x.len() {
        let (jx, jy) = jitter(seed, i);
        g.x[i] = step_ulp(g.x[i], jx);
        g.y[i] = step_ulp(g.y[i], jy);
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jitter_moves_by_at_most_one_ulp_in_both_directions() {
        let g = Geometry {
            cm_ref: [0.25, 0.0],
            x: (0..200).map(|i| 1.0 - i as f64 / 199.0).collect(),
            y: (0..200).map(|i| 0.01 * ((i % 7) as f64 - 3.0)).collect(),
            generator: None,
        };
        let j = jitter_geometry(&g, 3);
        let mut seen = [0usize; 3];
        for (a, b) in g.x.iter().chain(&g.y).zip(j.x.iter().chain(&j.y)) {
            let moved = if a == b {
                1
            } else if *b == a.next_up() {
                2
            } else {
                assert_eq!(*b, a.next_down());
                0
            };
            seen[moved] += 1;
        }
        // every direction occurs: this is a jitter, not a translation
        assert!(seen.iter().all(|&n| n > 50), "{seen:?}");
    }
}
