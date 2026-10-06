//! (f) Reference integrity.
//!
//! Checks on the reference and its fixtures rather than on yFoil.
//!
//! Categories and rules: `docs/conventions/testing.md`.

#[path = "../common/fixtures/mod.rs"]
mod fixtures;
#[path = "../common/utilities/mod.rs"]
mod utilities;

mod blsolv_conditioning;
mod case_tests;
mod inputs;
mod polar_reseed;
mod tgap_handoff;
