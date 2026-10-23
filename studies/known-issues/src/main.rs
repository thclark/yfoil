//! Known-issues cases (`docs/testing/known-issues/`): runs that show one of XFOIL's known
//! issues (`docs/user-guide/reference/known-xfoil-issues.md`) arising, and whether yFoil does the same. The study is the
//! series cases' library run over `group = "known-issues"`; only its descriptor
//! ([`series_cases::KNOWN_ISSUES`]) differs, and its figures are the series cases' (`series_cases::figures`).
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo xtask fixtures --group known-issues
//! cargo run --release -p known-issues -- [--docs]
//! cargo run --release -p known-issues -- --render studies/known-issues/runs/<UTC> [--docs]
//! ```

fn main() {
    series_cases::run(&series_cases::KNOWN_ISSUES);
}
