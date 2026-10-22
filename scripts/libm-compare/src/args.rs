//! The argument sets: seeded random sweeps (each argument a pure function of the set's seed and
//! its index, so every target and thread count sees the same arguments), special values, and
//! CORE-MATH's worst cases.

use crate::candidates::Arity;
use std::path::Path;

/// How a set's arguments are produced
pub enum Source {
    /// `count` arguments, the `i`-th from `gen(splitmix(seed, i))`
    Random {
        count: u64,
        seed: u64,
        gen: fn(&mut Rng) -> [f64; 2],
    },
    /// An explicit list
    List(Vec<[f64; 2]>),
}

pub struct Set {
    pub name: String,
    pub source: Source,
}

impl Set {
    pub fn len(&self) -> u64 {
        match &self.source {
            Source::Random { count, .. } => *count,
            Source::List(v) => v.len() as u64,
        }
    }

    #[inline]
    pub fn arg(&self, i: u64) -> [f64; 2] {
        match &self.source {
            Source::Random { seed, gen, .. } => gen(&mut Rng::at(*seed, i)),
            Source::List(v) => v[i as usize],
        }
    }
}

/// SplitMix64, seeded per argument
pub struct Rng(u64);

impl Rng {
    fn at(seed: u64, i: u64) -> Self {
        Rng(seed ^ i.wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(17))
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1)
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    fn value(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Any bit pattern (NaNs and infinities included, 1 in 2048)
    fn bits(&mut self) -> f64 {
        f64::from_bits(self.next())
    }

    /// Any non-negative bit pattern
    fn bits_pos(&mut self) -> f64 {
        f64::from_bits(self.next() >> 1)
    }

    /// Uniform in bits over [0, hi], random sign
    fn bits_upto(&mut self, hi: f64) -> f64 {
        let r = self.next();
        let x = f64::from_bits((r >> 1) % (hi.to_bits() + 1));
        if r & 1 == 1 {
            -x
        } else {
            x
        }
    }

    /// Uniform in bits over |x| in [lo, hi], random sign
    fn bits_between(&mut self, lo: f64, hi: f64) -> f64 {
        let r = self.next();
        let x = f64::from_bits(lo.to_bits() + (r >> 1) % (hi.to_bits() - lo.to_bits() + 1));
        if r & 1 == 1 {
            -x
        } else {
            x
        }
    }

    /// Up to `k` ulps either side of `x`
    fn jog(&mut self, x: f64, k: u64) -> f64 {
        let d = self.next() % (2 * k + 1);
        f64::from_bits((x.to_bits() as i64 + d as i64 - k as i64) as u64)
    }
}

const MAX_SUBNORMAL: f64 = f64::from_bits(0x000f_ffff_ffff_ffff);

fn r(count: u64, seed: u64, gen: fn(&mut Rng) -> [f64; 2]) -> Source {
    Source::Random { count, seed, gen }
}

/// The random sets for `function`, `n` arguments each (some sets take a fraction)
fn random_sets(function: &str, n: u64) -> Vec<(&'static str, Source)> {
    let h = n / 2;
    match function {
        "exp" => vec![
            ("bits", r(n, 1, |g| [g.bits(), 0.0])),
            ("value[-745.2,709.8]", r(n, 2, |g| [g.value(-745.2, 709.8), 0.0])),
            ("value[-1,1]", r(h, 3, |g| [g.value(-1.0, 1.0), 0.0])),
            ("bits|x|<=2^-10", r(h, 4, |g| [g.bits_upto(1.0 / 1024.0), 0.0])),
        ],
        "ln" | "log10" => vec![
            ("bits>=0", r(n, 11, |g| [g.bits_pos(), 0.0])),
            ("value[0.5,2]", r(n, 12, |g| [g.value(0.5, 2.0), 0.0])),
            ("value[0,1e3]", r(h, 13, |g| [g.value(0.0, 1e3), 0.0])),
            ("1+-2^20ulp", r(h, 14, |g| [g.jog(1.0, 1 << 20), 0.0])),
        ],
        "pow" => vec![
            (
                "x:bits>=0,y:value[-64,64]",
                r(n, 21, |g| [g.bits_pos(), g.value(-64.0, 64.0)]),
            ),
            (
                "x:value[0,4],y:value[-8,8]",
                r(n, 22, |g| [g.value(0.0, 4.0), g.value(-8.0, 8.0)]),
            ),
            (
                "x:1+-2^30ulp,y:bits|y|<=2^40",
                r(h, 23, |g| [g.jog(1.0, 1 << 30), g.bits_upto(1099511627776.0)]),
            ),
            (
                "x:bits,y:int[-64,64]",
                r(h, 24, |g| [g.bits(), (g.next() % 129) as f64 - 64.0]),
            ),
            (
                "x:bits,y:half-int[-64,64]",
                r(h, 25, |g| [g.bits(), (g.next() % 257) as f64 * 0.5 - 64.0]),
            ),
            ("x:bits,y:bits", r(h, 26, |g| [g.bits(), g.bits()])),
        ],
        "tanh" => vec![
            ("bits", r(n, 31, |g| [g.bits(), 0.0])),
            ("value[-20,20]", r(n, 32, |g| [g.value(-20.0, 20.0), 0.0])),
            ("value[-1,1]", r(h, 33, |g| [g.value(-1.0, 1.0), 0.0])),
        ],
        "sin" | "cos" | "sin_cos" => vec![
            ("bits", r(n, 41, |g| [g.bits(), 0.0])),
            ("value[-10,10]", r(n, 42, |g| [g.value(-10.0, 10.0), 0.0])),
            ("value[-1e6,1e6]", r(h, 43, |g| [g.value(-1e6, 1e6), 0.0])),
            (
                "bits|x|in[2^20,2^1023]",
                r(h, 44, |g| [g.bits_between(1048576.0, f64::MAX), 0.0]),
            ),
            (
                "k*pi/2+-8ulp",
                r(h, 45, |g| {
                    // k up to 2^e, e uniform in [0, 60]
                    let e = g.next() % 61;
                    let k = (g.next() >> (63 - e)) as f64;
                    [g.jog(k * std::f64::consts::FRAC_PI_2, 8), 0.0]
                }),
            ),
        ],
        "atan" => vec![
            ("bits", r(n, 51, |g| [g.bits(), 0.0])),
            ("value[-10,10]", r(n, 52, |g| [g.value(-10.0, 10.0), 0.0])),
            ("value[-1,1]", r(h, 53, |g| [g.value(-1.0, 1.0), 0.0])),
        ],
        "atan2" => vec![
            ("y:bits,x:bits", r(n, 61, |g| [g.bits(), g.bits()])),
            ("value[-1,1]^2", r(n, 62, |g| [g.value(-1.0, 1.0), g.value(-1.0, 1.0)])),
            (
                "y:bits,x:y*value[-2,2]",
                r(h, 63, |g| {
                    let y = g.bits();
                    [y, y * g.value(-2.0, 2.0)]
                }),
            ),
            (
                "subnormal^2",
                r(h, 64, |g| [g.bits_upto(MAX_SUBNORMAL), g.bits_upto(MAX_SUBNORMAL)]),
            ),
            (
                "|y|,|x|in[2^-1032,2^-1012]",
                r(h, 65, |g| {
                    [
                        g.bits_between(pow2(-1032), pow2(-1012)),
                        g.bits_between(pow2(-1032), pow2(-1012)),
                    ]
                }),
            ),
        ],
        "asin" => vec![
            ("value[-1,1]", r(n, 71, |g| [g.value(-1.0, 1.0), 0.0])),
            ("bits|x|<=1", r(n, 72, |g| [g.bits_upto(1.0), 0.0])),
            (
                "+-1-2^20ulp",
                r(h, 73, |g| {
                    let x = f64::from_bits(1f64.to_bits() - g.next() % (1 << 20));
                    [if g.next() & 1 == 1 { -x } else { x }, 0.0]
                }),
            ),
        ],
        "hypot" => vec![
            ("bits^2", r(n, 81, |g| [g.bits(), g.bits()])),
            (
                "value[-10,10]^2",
                r(n, 82, |g| [g.value(-10.0, 10.0), g.value(-10.0, 10.0)]),
            ),
            (
                "x:bits,y:x*value[-1,1]",
                r(h, 83, |g| {
                    let x = g.bits();
                    [x, x * g.value(-1.0, 1.0)]
                }),
            ),
            (
                "x:bits,y:x*2^-27*value[-1,1]",
                r(h, 84, |g| {
                    let x = g.bits();
                    [x, x * 7.450580596923828e-9 * g.value(-1.0, 1.0)]
                }),
            ),
            (
                "subnormal^2",
                r(h, 85, |g| [g.bits_upto(MAX_SUBNORMAL), g.bits_upto(MAX_SUBNORMAL)]),
            ),
            (
                "|x|,|y|in[2^-1032,2^-1012]",
                r(h, 86, |g| {
                    [
                        g.bits_between(pow2(-1032), pow2(-1012)),
                        g.bits_between(pow2(-1032), pow2(-1012)),
                    ]
                }),
            ),
        ],
        other => panic!("no random sets for {other}"),
    }
}

/// Values every function sees: zeros, subnormals, extremes, infinities, NaNs (three payloads),
/// small integers, constants, every power of two, and each within a few ulps
fn unary_specials(function: &str) -> Vec<f64> {
    let mut base = vec![
        0.0,
        f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::INFINITY,
        1.0,
        0.5,
        2.0,
        10.0,
        std::f64::consts::FRAC_PI_4,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        std::f64::consts::E,
        std::f64::consts::LN_2,
        std::f64::consts::LN_10,
    ];
    base.extend((3..=64).map(f64::from));
    base.extend((-1074..=1023).map(pow2));
    match function {
        // overflow, the subnormal range and underflow to zero
        "exp" => base.extend([
            709.782712893384,
            709.78,
            -708.3964185322641,
            -745.1332191019411,
            -745.1332191019412,
        ]),
        // saturation to +-1 and the small-argument branch
        "tanh" => base.extend([19.061547465398498, 22.0, 0.5493061443340549]),
        "sin" | "cos" | "sin_cos" => base.extend([1e22, 1.0e300, 6381956970095103.0 * pow2(797), 5.0e-325]),
        _ => {}
    }
    let mut out = Vec::new();
    for x in base {
        for d in -4i64..=4 {
            let y = f64::from_bits((x.to_bits() as i64 + d) as u64);
            if y >= 0.0 {
                out.push(y);
                out.push(-y);
            }
        }
    }
    out.extend([
        f64::NAN,
        -f64::NAN,
        f64::from_bits(0x7ff0_0000_0000_0001),
        f64::from_bits(0xfff4_0000_0000_0000),
    ]);
    out
}

fn binary_specials(function: &str) -> Vec<[f64; 2]> {
    let mut v: Vec<f64> = vec![
        0.0,
        f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::MIN_POSITIVE,
        f64::EPSILON,
        0.25,
        0.5,
        1.0 - f64::EPSILON / 2.0,
        1.0,
        1.0 + f64::EPSILON,
        1.5,
        2.0,
        3.0,
        10.0,
        1e10,
        pow2(52),
        pow2(52) + 1.0,
        pow2(53),
        pow2(53) + 2.0,
        1e300,
        f64::MAX,
        f64::INFINITY,
    ];
    v.extend(v.clone().into_iter().map(|x| -x).collect::<Vec<_>>());
    v.extend([f64::NAN, -f64::NAN, f64::from_bits(0x7ff0_0000_0000_0001)]);
    let mut out: Vec<[f64; 2]> = v.iter().flat_map(|&x| v.iter().map(move |&y| [x, y])).collect();
    if function == "pow" {
        let mut ys: Vec<f64> = (-40..=40).map(|k| k as f64 * 0.5).collect();
        ys.extend([
            pow2(52) + 1.0,
            -(pow2(52) + 1.0),
            pow2(53),
            -pow2(53),
            1e300,
            -1e300,
            1075.0,
            -1075.0,
            1024.0,
            -1024.0,
        ]);
        let xs = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            2.0,
            -2.0,
            0.5,
            -0.5,
            -8.0,
            -27.0,
            10.0,
            1e-300,
            -1e-300,
            1e300,
            -1e300,
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::INFINITY,
            -f64::INFINITY,
            f64::NAN,
            std::f64::consts::E,
        ];
        for &x in &xs {
            for &y in &ys {
                out.push([x, y]);
            }
        }
        // exact results: (k, 1/n) and (k^n, 1/n)
        for k in 1..=40 {
            for n in [2.0, 3.0, 4.0, 5.0] {
                let kf = k as f64;
                out.push([kf, 1.0 / n]);
                out.push([kf * kf, 0.5]);
                out.push([kf * kf * kf, 1.0 / 3.0]);
                out.push([kf, n]);
            }
        }
    }
    out
}

fn pow2(e: i32) -> f64 {
    if e >= -1022 {
        f64::from_bits(((e + 1023) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (e + 1074))
    }
}

pub fn sets(function: &str, arity: Arity, n: u64, worst_cases: &[&str], wc_dir: Option<&Path>) -> Vec<Set> {
    let mut out: Vec<Set> = random_sets(function, n)
        .into_iter()
        .map(|(name, source)| Set {
            name: name.to_owned(),
            source,
        })
        .collect();
    let specials = match arity {
        Arity::Unary => unary_specials(function).into_iter().map(|x| [x, 0.0]).collect(),
        Arity::Binary => binary_specials(function),
    };
    out.push(Set {
        name: "specials".into(),
        source: Source::List(specials),
    });
    if let Some(dir) = wc_dir {
        for wc in worst_cases {
            let path = dir.join(format!("{wc}.wc"));
            let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            let list = parse_worst_cases(&text, arity).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            out.push(Set {
                name: format!("worst-cases:{wc}.wc"),
                source: Source::List(list),
            });
        }
    }
    out
}

/// CORE-MATH's `.wc` format: one argument (or `x,y`) per line, C99 hex floats, decimal, `inf`,
/// `nan`, `snan`; `#` starts a comment, and anything after a field's first word is a note
fn parse_worst_cases(text: &str, arity: Arity) -> Result<Vec<[f64; 2]>, String> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        // a field may be followed by a note (`0x1.4f1d73be27a31p+1 44` in sin.wc: the run length)
        let fields: Vec<&str> = line
            .split(',')
            .map(|f| f.split_whitespace().next().unwrap_or(""))
            .collect();
        let want = if arity == Arity::Unary { 1 } else { 2 };
        if fields.len() != want {
            return Err(format!("line {}: {} fields, expected {want}", n + 1, fields.len()));
        }
        let x = parse_float(fields[0]).ok_or_else(|| format!("line {}: {}", n + 1, fields[0]))?;
        let y = if want == 2 {
            parse_float(fields[1]).ok_or_else(|| format!("line {}: {}", n + 1, fields[1]))?
        } else {
            0.0
        };
        out.push([x, y]);
    }
    Ok(out)
}

pub fn parse_float(s: &str) -> Option<f64> {
    let (neg, body) = match s.as_bytes().first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    let x = match body {
        "inf" => f64::INFINITY,
        "nan" => f64::NAN,
        "snan" => f64::from_bits(0x7ff4_0000_0000_0000),
        _ if body.starts_with("0x") || body.starts_with("0X") => parse_hex(&body[2..])?,
        _ => body.parse().ok()?,
    };
    Some(if neg { -x } else { x })
}

/// `h.hhhp±d` (or `h.hhh`), rounded to nearest-even like `strtod` (a few of CORE-MATH's entries
/// carry more digits than the binary64 they stand for, e.g. `0x1.08fa30093e7d8p-1074`)
fn parse_hex(s: &str) -> Option<f64> {
    let (mant, exp) = s.split_once(['p', 'P']).unwrap_or((s, "0"));
    let mut e: i64 = exp.parse().ok()?;
    let mut m: u128 = 0;
    let mut seen_point = false;
    for c in mant.chars() {
        if c == '.' {
            seen_point = true;
            continue;
        }
        m = m.checked_mul(16)? + c.to_digit(16)? as u128;
        if seen_point {
            e -= 4;
        }
    }
    if m == 0 {
        return Some(0.0);
    }
    // value = m * 2^e, in [2^top, 2^(top+1))
    let lead = 127 - m.leading_zeros() as i64;
    let top = lead + e;
    // bits to drop so that 53 remain (fewer for a subnormal, whose last bit is 2^-1074)
    let drop = if top >= -1022 { lead - 52 } else { -1074 - e };
    let (mut q, mut exp2) = if drop <= 0 {
        (m << -drop, e + drop)
    } else if drop >= 127 {
        (0, e + drop)
    } else {
        let q = m >> drop;
        let rem = m & ((1u128 << drop) - 1);
        let half = 1u128 << (drop - 1);
        let up = rem > half || (rem == half && q & 1 == 1);
        (q + up as u128, e + drop)
    };
    if q == 1u128 << 53 {
        q >>= 1;
        exp2 += 1;
    }
    if q >= 1u128 << 52 {
        let biased = exp2 + 52 + 1023;
        if biased >= 2047 {
            return Some(f64::INFINITY);
        }
        Some(f64::from_bits(((biased as u64) << 52) | (q as u64 & ((1 << 52) - 1))))
    } else {
        Some(f64::from_bits(q as u64))
    }
}

#[cfg(test)]
mod tests {
    use super::parse_float;

    #[test]
    fn hex_floats() {
        assert_eq!(parse_float("0x1p+0"), Some(1.0));
        assert_eq!(parse_float("-0x1.8p+1"), Some(-3.0));
        assert_eq!(parse_float("0x1.62e42fefa39efp+9"), Some(709.782712893384));
        assert_eq!(parse_float("0x1p-1074"), Some(f64::from_bits(1)));
        assert_eq!(parse_float("0x0.0000000000001p-1022"), Some(f64::from_bits(1)));
        assert_eq!(parse_float("0x1.fffffffffffffp+1023"), Some(f64::MAX));
        assert_eq!(parse_float("+0"), Some(0.0));
        assert_eq!(parse_float("0x1"), Some(1.0));
        // more digits than a subnormal holds: 1.0351 * 2^-1074 rounds to 2^-1074
        assert_eq!(parse_float("0x1.08fa30093e7d8p-1074"), Some(f64::from_bits(1)));
        assert_eq!(parse_float("0x1.8p-1074"), Some(f64::from_bits(2)));
        assert_eq!(parse_float("0x1.0000000000000cp+0"), Some(1.0 + f64::EPSILON));
        assert_eq!(parse_float("0x1.00000000000008p+0"), Some(1.0));
        assert_eq!(parse_float("0x1.fffffffffffff8p+0"), Some(2.0));
        assert_eq!(parse_float("-inf"), Some(f64::NEG_INFINITY));
    }
}
