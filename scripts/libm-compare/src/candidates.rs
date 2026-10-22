//! The twelve functions yFoil takes from `cr_libm`, and each candidate's implementation of them.

/// One implementation of one function.
#[derive(Clone, Copy)]
pub enum Imp {
    Unary(fn(f64) -> f64),
    Binary(fn(f64, f64) -> f64),
    SinCos(fn(f64) -> (f64, f64)),
}

impl Imp {
    /// The result(s) for `args` (`args[1]` is unused by a unary function; `sin_cos` gives two)
    #[inline(always)]
    pub fn eval(self, args: [f64; 2]) -> [f64; 2] {
        match self {
            Imp::Unary(f) => [f(args[0]), 0.0],
            Imp::Binary(f) => [f(args[0], args[1]), 0.0],
            Imp::SinCos(f) => {
                let (s, c) = f(args[0]);
                [s, c]
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Arity {
    Unary,
    Binary,
}

pub struct Function {
    /// cr_libm's name
    pub name: &'static str,
    pub arity: Arity,
    /// CORE-MATH worst-case files (`<name>.wc` at the pinned commit) that apply to it
    pub worst_cases: &'static [&'static str],
    /// (candidate, implementation); the first is the reference when it is `core-math (C)`
    pub imps: Vec<(&'static str, Imp)>,
}

pub const C_REFERENCE: &str = "core-math-c";

/// metallic 0.3.1's `atan2` misrounds, by up to 6 % relative, when both arguments are subnormal.
/// atan2 depends on y/x only, and scaling two subnormals by 2^54 is exact, so this is metallic
/// with that case moved into the normal range: the workaround `cr-libm` would carry until fixed.
pub fn metallic_atan2_scaled(y: f64, x: f64) -> f64 {
    const TWO_54: f64 = 18014398509481984.0;
    if y.abs() < f64::MIN_POSITIVE && x.abs() < f64::MIN_POSITIVE {
        metallic::atan2(y * TWO_54, x * TWO_54)
    } else {
        metallic::atan2(y, x)
    }
}

pub fn functions() -> Vec<Function> {
    use Imp::*;
    #[allow(unused_mut)]
    let mut fs = vec![
        Function {
            name: "exp",
            arity: Arity::Unary,
            worst_cases: &["exp"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_exp)),
                ("metallic", Unary(metallic::exp)),
                ("libm-0.2.16", Unary(libm::exp)),
                ("libm-main", Unary(libm_main::exp)),
            ],
        },
        Function {
            name: "ln",
            arity: Arity::Unary,
            worst_cases: &["log"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_log)),
                ("metallic", Unary(metallic::log)),
                ("libm-0.2.16", Unary(libm::log)),
                ("libm-main", Unary(libm_main::log)),
            ],
        },
        Function {
            name: "log10",
            arity: Arity::Unary,
            worst_cases: &["log10"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_log10)),
                ("metallic", Unary(metallic::log10)),
                ("libm-0.2.16", Unary(libm::log10)),
                ("libm-main", Unary(libm_main::log10)),
            ],
        },
        Function {
            name: "pow",
            arity: Arity::Binary,
            worst_cases: &["pow"],
            imps: vec![
                ("pxfm", Binary(pxfm::f_pow)),
                ("metallic", Binary(metallic::pow)),
                ("libm-0.2.16", Binary(libm::pow)),
                ("libm-main", Binary(libm_main::pow)),
            ],
        },
        Function {
            name: "tanh",
            arity: Arity::Unary,
            worst_cases: &["tanh"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_tanh)),
                ("metallic", Unary(metallic::tanh)),
                ("libm-0.2.16", Unary(libm::tanh)),
                ("libm-main", Unary(libm_main::tanh)),
            ],
        },
        Function {
            name: "sin",
            arity: Arity::Unary,
            worst_cases: &["sin"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_sin)),
                ("metallic", Unary(metallic::sin)),
                ("libm-0.2.16", Unary(libm::sin)),
                ("libm-main", Unary(libm_main::sin)),
            ],
        },
        Function {
            name: "cos",
            arity: Arity::Unary,
            worst_cases: &["cos"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_cos)),
                ("metallic", Unary(metallic::cos)),
                ("libm-0.2.16", Unary(libm::cos)),
                ("libm-main", Unary(libm_main::cos)),
            ],
        },
        Function {
            name: "sin_cos",
            arity: Arity::Unary,
            worst_cases: &["sincos", "sin", "cos"],
            imps: vec![
                ("pxfm", SinCos(pxfm::f_sincos)),
                ("metallic", SinCos(metallic::sincos)),
                ("libm-0.2.16", SinCos(libm::sincos)),
                ("libm-main", SinCos(libm_main::sincos)),
            ],
        },
        Function {
            name: "atan",
            arity: Arity::Unary,
            worst_cases: &["atan"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_atan)),
                ("metallic", Unary(metallic::atan)),
                ("libm-0.2.16", Unary(libm::atan)),
                ("libm-main", Unary(libm_main::atan)),
            ],
        },
        Function {
            name: "atan2",
            arity: Arity::Binary,
            worst_cases: &["atan2"],
            imps: vec![
                ("pxfm", Binary(pxfm::f_atan2)),
                ("metallic", Binary(metallic::atan2)),
                ("metallic-scaled", Binary(metallic_atan2_scaled)),
                ("libm-0.2.16", Binary(libm::atan2)),
                ("libm-main", Binary(libm_main::atan2)),
            ],
        },
        Function {
            name: "asin",
            arity: Arity::Unary,
            worst_cases: &["asin"],
            imps: vec![
                ("pxfm", Unary(pxfm::f_asin)),
                ("metallic", Unary(metallic::asin)),
                ("libm-0.2.16", Unary(libm::asin)),
                ("libm-main", Unary(libm_main::asin)),
            ],
        },
        Function {
            name: "hypot",
            arity: Arity::Binary,
            worst_cases: &["hypot"],
            imps: vec![
                ("pxfm", Binary(pxfm::f_hypot)),
                ("metallic", Binary(metallic::hypot)),
                ("libm-0.2.16", Binary(libm::hypot)),
                ("libm-main", Binary(libm_main::hypot)),
            ],
        },
    ];
    #[cfg(feature = "c-reference")]
    for f in &mut fs {
        let c = match f.name {
            "exp" => Unary(cr_libm::exp),
            "ln" => Unary(cr_libm::ln),
            "log10" => Unary(cr_libm::log10),
            "pow" => Binary(cr_libm::pow),
            "tanh" => Unary(cr_libm::tanh),
            "sin" => Unary(cr_libm::sin),
            "cos" => Unary(cr_libm::cos),
            "sin_cos" => SinCos(cr_libm::sin_cos),
            "atan" => Unary(cr_libm::atan),
            "atan2" => Binary(cr_libm::atan2),
            "asin" => Unary(cr_libm::asin),
            "hypot" => Binary(cr_libm::hypot),
            other => unreachable!("{other}"),
        };
        f.imps.insert(0, (C_REFERENCE, c));
    }
    fs
}
