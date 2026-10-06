//! (b) Execution equivalence.
//!
//! Each test runs an execution of yFoil — one step seeded from XFOIL's dumped state, or a
//! well-conditioned whole run — and requires XFOIL's branch decisions and values.
//!
//! Categories and rules: `docs/conventions/testing.md`.

#[path = "../common/fixtures/mod.rs"]
mod fixtures;
#[path = "../common/utilities/mod.rs"]
mod utilities;

mod analysis;
mod coverage;
mod polar;
mod polar_break;
mod prologue_dij;
mod prologue_pointers;
mod speccl;
mod update;
mod viscal;
