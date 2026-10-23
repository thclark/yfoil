//! Shared with the tests and the other studies, included by path: the tolerance constants
//! (`tests/common/utilities/`), and the studies' floor comparison
//! (`studies/study-support/records.rs`, which the tests do not use).
#![allow(dead_code)]
#[path = "../../../study-support/records.rs"]
pub mod records;
#[path = "../../../../tests/common/utilities/tolerances.rs"]
pub mod tolerances;
