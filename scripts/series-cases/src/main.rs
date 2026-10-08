//! Series cases (`docs/validation/series-cases/`): the aerofoil series and flow regimes of interest
//! under normal conditions, yFoil against the reference point by point (the study itself is the
//! library, shared with the known-issues cases).
//!
//! Usage (from the repository root):
//!
//! ```text
//! cargo run --release -p series-cases -- [--docs]
//! ```

fn main() {
    series_cases::run(&series_cases::SERIES);
}
