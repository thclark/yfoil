//! yFoil's transcendentals, correctly rounded.
//!
//! The standard library's `f64::exp`, `ln`, `powf` and the rest call the host's libm, and the host
//! libms differ: Apple libSystem and glibc disagree by one ULP on 0.1 % (`exp`, `ln`, `pow`) to 18 %
//! (`tanh`) of inputs, which a Newton solve can carry to O(1) past stall. These are CORE-MATH's
//! correctly rounded functions instead: each result is the exact value rounded to nearest, a
//! property of the argument alone, so it is the same on every host. The XFOIL reference is linked
//! against the same sources (`scripts/xfoil-build/core_math_libm.c`), so it and yFoil get the same
//! bits too.
//!
//! The `std` methods are disallowed in the workspace (`clippy.toml`); call these instead. `sqrt` is
//! not here: IEEE 754 requires it correctly rounded, so `f64::sqrt` already is.

#[cfg(not(feature = "metallic"))]
mod ffi {
    extern "C" {
        pub fn cr_exp(x: f64) -> f64;
        pub fn cr_log(x: f64) -> f64;
        pub fn cr_log10(x: f64) -> f64;
        pub fn cr_pow(x: f64, y: f64) -> f64;
        pub fn cr_tanh(x: f64) -> f64;
        pub fn cr_sin(x: f64) -> f64;
        pub fn cr_cos(x: f64) -> f64;
        pub fn cr_sincos(x: f64, sin: *mut f64, cos: *mut f64);
        pub fn cr_atan(x: f64) -> f64;
        pub fn cr_atan2(y: f64, x: f64) -> f64;
        pub fn cr_asin(x: f64) -> f64;
        pub fn cr_hypot(x: f64, y: f64) -> f64;
    }
}

/// The C: CORE-MATH compiled by build.rs
#[cfg(not(feature = "metallic"))]
mod imp {
    use super::ffi;

    // SAFETY (every call): a pure function of its arguments, no pointers
    #[inline]
    pub fn exp(x: f64) -> f64 {
        unsafe { ffi::cr_exp(x) }
    }
    #[inline]
    pub fn ln(x: f64) -> f64 {
        unsafe { ffi::cr_log(x) }
    }
    #[inline]
    pub fn log10(x: f64) -> f64 {
        unsafe { ffi::cr_log10(x) }
    }
    #[inline]
    pub fn pow(x: f64, y: f64) -> f64 {
        unsafe { ffi::cr_pow(x, y) }
    }
    #[inline]
    pub fn tanh(x: f64) -> f64 {
        unsafe { ffi::cr_tanh(x) }
    }
    #[inline]
    pub fn sin(x: f64) -> f64 {
        unsafe { ffi::cr_sin(x) }
    }
    #[inline]
    pub fn cos(x: f64) -> f64 {
        unsafe { ffi::cr_cos(x) }
    }
    #[inline]
    pub fn sin_cos(x: f64) -> (f64, f64) {
        let (mut s, mut c) = (0.0, 0.0);
        // SAFETY: both pointers are to live, exclusively borrowed locals
        unsafe { ffi::cr_sincos(x, &mut s, &mut c) };
        (s, c)
    }
    #[inline]
    pub fn atan(x: f64) -> f64 {
        unsafe { ffi::cr_atan(x) }
    }
    #[inline]
    pub fn atan2(y: f64, x: f64) -> f64 {
        unsafe { ffi::cr_atan2(y, x) }
    }
    #[inline]
    pub fn asin(x: f64) -> f64 {
        unsafe { ffi::cr_asin(x) }
    }
    #[inline]
    pub fn hypot(x: f64, y: f64) -> f64 {
        unsafe { ffi::cr_hypot(x, y) }
    }
}

/// Trial: metallic, a pure-Rust correctly rounded libm, in place of the C. scripts/libm-compare
/// finds it bit-identical to the C on every function and target, except signalling-NaN arguments
/// and metallic 0.3.1's `atan2` with both arguments subnormal, which `atan2` here works around.
#[cfg(feature = "metallic")]
mod imp {
    pub use metallic::{asin, atan, cos, exp, hypot, log as ln, log10, pow, sin, sincos as sin_cos, tanh};

    /// metallic 0.3.1 misrounds `atan2` (by up to 6 % relative) when both arguments are subnormal.
    /// The result depends on y/x alone and scaling two subnormals by 2^54 is exact, so the case is
    /// moved into the normal range, where metallic is correctly rounded.
    #[inline]
    pub fn atan2(y: f64, x: f64) -> f64 {
        const TWO_54: f64 = 18014398509481984.0;
        if y.abs() < f64::MIN_POSITIVE && x.abs() < f64::MIN_POSITIVE {
            metallic::atan2(y * TWO_54, x * TWO_54)
        } else {
            metallic::atan2(y, x)
        }
    }
}

/// e^x (Fortran `EXP`).
#[inline]
pub fn exp(x: f64) -> f64 {
    imp::exp(x)
}

/// Natural logarithm (Fortran `LOG`).
#[inline]
pub fn ln(x: f64) -> f64 {
    imp::ln(x)
}

/// Base-10 logarithm (Fortran `LOG10`).
#[inline]
pub fn log10(x: f64) -> f64 {
    imp::log10(x)
}

/// x^y for a real exponent (Fortran `X**Y` with `Y` real).
#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    imp::pow(x, y)
}

/// Hyperbolic tangent (Fortran `TANH`).
#[inline]
pub fn tanh(x: f64) -> f64 {
    imp::tanh(x)
}

/// Sine (Fortran `SIN`).
#[inline]
pub fn sin(x: f64) -> f64 {
    imp::sin(x)
}

/// Cosine (Fortran `COS`).
#[inline]
pub fn cos(x: f64) -> f64 {
    imp::cos(x)
}

/// `(sin(x), cos(x))`, each correctly rounded, so equal to [`sin`] and [`cos`].
#[inline]
pub fn sin_cos(x: f64) -> (f64, f64) {
    imp::sin_cos(x)
}

/// Arctangent (Fortran `ATAN`).
#[inline]
pub fn atan(x: f64) -> f64 {
    imp::atan(x)
}

/// Four-quadrant arctangent of `y/x` (Fortran `ATAN2(Y, X)`).
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    imp::atan2(y, x)
}

/// Arcsine (Fortran `ASIN`).
#[inline]
pub fn asin(x: f64) -> f64 {
    imp::asin(x)
}

/// `sqrt(x² + y²)` without intermediate overflow.
#[inline]
pub fn hypot(x: f64, y: f64) -> f64 {
    imp::hypot(x, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The binary64 nearest each constant (the well-known bit patterns of e, ln 10, log10 2, √2 and
    // the std consts). They check the binding is wired to the right functions with the right
    // argument order, not CORE-MATH itself.
    #[test]
    fn known_values() {
        assert_eq!(exp(1.0).to_bits(), 0x4005bf0a8b145769); // e
        assert_eq!(ln(10.0).to_bits(), 0x40026bb1bbb55516); // ln 10
        assert_eq!(log10(2.0).to_bits(), 0x3fd34413509f79ff); // log10 2
        assert_eq!(pow(2.0, 0.5).to_bits(), 0x3ff6a09e667f3bcd); // √2
        assert_eq!(4.0 * atan(1.0), std::f64::consts::PI);
        assert_eq!(atan2(1.0, -1.0), 3.0 * std::f64::consts::FRAC_PI_4);
        assert_eq!(atan2(-0.0, 1.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(2.0 * asin(1.0), std::f64::consts::PI);
        assert_eq!(hypot(3.0, 4.0), 5.0);
        assert_eq!(sin(0.0), 0.0);
        assert_eq!(cos(0.0), 1.0);
        assert_eq!(tanh(0.0), 0.0);
        for x in [0.1, 1.0, 2.5, -7.25, 1e5] {
            assert_eq!(sin_cos(x), (sin(x), cos(x)));
        }
    }
}
