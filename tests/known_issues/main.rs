//! (d) Known XFOIL weaknesses.
//!
//! yFoil's defined behaviour where XFOIL is known to be incorrect or ill-conditioned; one module
//! per section of `docs/xfoil-known-issues.md`.
//!
//! Categories and rules: `docs/conventions/testing.md`.

#[path = "../common/fixtures/mod.rs"]
mod fixtures;
#[path = "../common/utilities/mod.rs"]
mod utilities;

mod kt_pole_7_7;
mod tgap_6_5;
