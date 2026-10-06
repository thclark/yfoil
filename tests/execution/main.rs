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
mod branches;
mod events;
mod polar;
mod prologue_dij;
mod prologue_pointers;
mod run;
mod runs;
mod step;
mod steps;
mod update;
mod viscal;
