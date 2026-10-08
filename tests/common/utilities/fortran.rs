//! Reading the reference's formatted reals.

/// A real the reference wrote with `ES24.16`. With a three-digit exponent the format has no room
/// for the `E` and writes `1.9674408635025233+116`; that form is read too. Anything else panics,
/// naming the token.
pub fn real(s: &str) -> f64 {
    let s = s.trim();
    if let Ok(v) = s.parse() {
        return v;
    }
    s.get(1..)
        .and_then(|t| t.rfind(['+', '-']))
        // the sign at s[j + 1] starts the exponent
        .and_then(|j| format!("{}e{}", &s[..=j], &s[j + 1..]).parse().ok())
        .unwrap_or_else(|| panic!("not a Fortran real: {s:?}"))
}

#[cfg(test)]
mod tests {
    use super::real;

    #[test]
    fn reads_the_three_digit_exponent_form() {
        assert_eq!(real("1.9674408635025233+116"), 1.9674408635025233e116);
        assert_eq!(real("3.0480272680724994-314"), 3.0480272680724994e-314);
        assert_eq!(real(" -1.5000000000000000E+00"), -1.5);
    }
}
