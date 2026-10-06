//! (a) Subroutine equivalence.
//!
//! Each test runs one well-encapsulated yFoil subroutine on inputs XFOIL dumped and compares
//! its outputs with XFOIL's, within a named tolerance (`common/utilities/tolerances.rs`).
//!
//! Categories and rules: `docs/conventions/testing.md`.

#[path = "../common/fixtures/mod.rs"]
mod fixtures;
#[path = "../common/utilities/mod.rs"]
mod utilities;

mod blsolv;
mod closures;
mod ggcalc;
mod mrchdu;
mod mrchue;
mod pane_legacy;
mod pangen;
mod pointers;
mod setbl;
mod speccl;
mod tgap;
mod transition;
mod xywake;
