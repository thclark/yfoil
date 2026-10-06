//! Content-addressed references for solved points.
//!
//! A point's `id` is a hash of everything the solve depended on: the yFoil version, the panel
//! geometry, the flow conditions, the operating point and its iteration limit, and the id of the
//! point whose converged boundary layer seeded it. The version is in there because the same
//! inputs may not give the same numbers across versions, so results from two of them must not
//! claim to be the same computation. So the id is a fingerprint of the computation, not a
//! position in a run.
//!
//! Three properties follow, and all three are wanted:
//!
//! - **Deterministic.** The same run produces the same ids, so a result file is byte-reproducible
//!   and two of them can be diffed.
//! - **Equal ids mean the same computation.** Results from different runs can be held together —
//!   a sweep restarted from a point of an earlier one, or two sweeps merged — and a point that
//!   appears in both carries one reference, not two colliding ones.
//! - **The chain is a Merkle chain.** Because the predecessor's id is hashed in, an id covers the
//!   whole history behind the point. Two points that agree on their operating point but were
//!   reached along different paths are correctly distinguished — which matters, because a polar is
//!   a state machine and the previous alpha's BL is this one's initial condition.
//!
//! It also removes a special case. The polar re-solves 0° after `INIT` to seed the downward leg,
//! and that re-solve is not itself a polar point. It has the same geometry, conditions, operating
//! point, iteration limit and (absent) predecessor as the first 0° solve, so it hashes to the same
//! id — which is exactly right, the fixture test asserts the two solves are identical — and the
//! downward leg's first point cites the recorded 0° point rather than a phantom.
//!
//! The hash is FNV-1a (128-bit), written out here rather than taken from `DefaultHasher`, whose
//! output is explicitly not stable across Rust releases and would break determinism on a toolchain
//! upgrade.

/// Characters of an id. Base-36 over 12 characters is ~62 bits. A content-addressed id is
/// compared across accumulated results rather than within one run, so the birthday bound is what
/// matters: 12 characters stay collision-free to ~2×10⁹ distinct points, where 8 would be ~41 bits
/// and around a 0.2 % chance of a collision by 10⁵ points. A collision would silently merge two
/// different solves, so the margin is worth four characters.
pub const ID_LENGTH: usize = 12;

const ALPHABET: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyz";
const FNV_OFFSET: u128 = 0x6c62272e07bb014262b821756295c58d;
const FNV_PRIME: u128 = 0x0000000001000000000000000000013b;

/// FNV-1a (128-bit) over a byte stream, accumulated field by field.
#[derive(Debug, Clone)]
pub struct IdHasher {
    state: u128,
}

impl Default for IdHasher {
    fn default() -> Self {
        Self { state: FNV_OFFSET }
    }
}

impl IdHasher {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn bytes(&mut self, b: &[u8]) -> &mut Self {
        for byte in b {
            self.state ^= *byte as u128;
            self.state = self.state.wrapping_mul(FNV_PRIME);
        }
        self
    }

    /// A field separator carrying the field's name, so that two different fields cannot alias by
    /// happening to hold the same bytes.
    pub fn field(&mut self, name: &str) -> &mut Self {
        self.bytes(&[0xff]).bytes(name.as_bytes()).bytes(&[0x00])
    }

    /// Exact bits, never a formatted decimal: the id must distinguish operating points that differ
    /// in the last place, which is the resolution everything else in this crate works at.
    pub fn f64(&mut self, name: &str, v: f64) -> &mut Self {
        self.field(name).bytes(&v.to_bits().to_le_bytes())
    }

    pub fn usize(&mut self, name: &str, v: usize) -> &mut Self {
        self.field(name).bytes(&(v as u64).to_le_bytes())
    }

    pub fn u8(&mut self, name: &str, v: u8) -> &mut Self {
        self.field(name).bytes(&[v])
    }

    pub fn str(&mut self, name: &str, v: &str) -> &mut Self {
        self.field(name).bytes(v.as_bytes())
    }

    /// `None` and `Some` are distinguished by a tag, so an absent predecessor cannot collide with
    /// a present one.
    pub fn opt_f64(&mut self, name: &str, v: Option<f64>) -> &mut Self {
        match v {
            Some(x) => self.field(name).bytes(&[1]).bytes(&x.to_bits().to_le_bytes()),
            None => self.field(name).bytes(&[0]),
        }
    }

    pub fn opt_str(&mut self, name: &str, v: Option<&str>) -> &mut Self {
        match v {
            Some(x) => self.field(name).bytes(&[1]).bytes(x.as_bytes()),
            None => self.field(name).bytes(&[0]),
        }
    }

    pub fn f64_slice(&mut self, name: &str, v: &[f64]) -> &mut Self {
        self.field(name).usize("len", v.len());
        for x in v {
            self.bytes(&x.to_bits().to_le_bytes());
        }
        self
    }

    /// The accumulated state as a base-36 id of [`ID_LENGTH`] characters.
    pub fn finish(&self) -> String {
        let mut v = self.state;
        let mut id = String::with_capacity(ID_LENGTH);
        for _ in 0..ID_LENGTH {
            id.push(ALPHABET[(v % ALPHABET.len() as u128) as usize] as char);
            v /= ALPHABET.len() as u128;
        }
        id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_is_the_declared_length_and_base36() {
        let id = IdHasher::new().f64("alpha", 0.5).finish();
        assert_eq!(id.len(), ID_LENGTH);
        assert!(id.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit()));
    }

    #[test]
    fn the_same_inputs_give_the_same_id() {
        let a = IdHasher::new().f64("alpha", 0.5).usize("n", 3).finish();
        let b = IdHasher::new().f64("alpha", 0.5).usize("n", 3).finish();
        assert_eq!(a, b);
    }

    #[test]
    fn one_ulp_of_the_operating_point_changes_the_id() {
        let a = IdHasher::new().f64("alpha", 0.5).finish();
        let b = IdHasher::new()
            .f64("alpha", f64::from_bits(0.5f64.to_bits() + 1))
            .finish();
        assert_ne!(a, b, "ids must resolve to the precision the solver works at");
    }

    #[test]
    fn fields_do_not_alias_by_holding_the_same_bytes() {
        let a = IdHasher::new().f64("mach", 0.3).f64("ncrit", 9.0).finish();
        let b = IdHasher::new().f64("ncrit", 0.3).f64("mach", 9.0).finish();
        assert_ne!(a, b);
    }

    #[test]
    fn an_absent_predecessor_differs_from_a_present_one() {
        let a = IdHasher::new().opt_str("from", None).finish();
        let b = IdHasher::new().opt_str("from", Some("")).finish();
        assert_ne!(a, b);
    }

    #[test]
    fn the_predecessor_is_part_of_the_identity() {
        // the same operating point reached along two different paths is two different solves
        let a = IdHasher::new()
            .f64("alpha", 2.0)
            .opt_str("from", Some("aaaaaaaaaaaa"))
            .finish();
        let b = IdHasher::new()
            .f64("alpha", 2.0)
            .opt_str("from", Some("bbbbbbbbbbbb"))
            .finish();
        assert_ne!(a, b);
    }
}
