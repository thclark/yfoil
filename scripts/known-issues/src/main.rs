//! Known-issues cases (`docs/validation/known-issues/`): runs that show one of XFOIL's known
//! issues (`docs/xfoil-known-issues.md`) arising, and whether yFoil does the same. The study is the
//! series cases' library run over `group = "known-issues"`; only its descriptor
//! ([`series_cases::KNOWN_ISSUES`]) differs.
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo xtask fixtures --group known-issues
//! cargo run --release -p known-issues -- [--docs]
//! ```

fn main() {
    series_cases::run(&series_cases::KNOWN_ISSUES);
}
