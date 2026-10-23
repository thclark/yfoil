//! (b) Execution equivalence, branch by branch: every branch of the translated subroutines that
//! a step of a reference run takes is tested by replaying one step that takes it (`step.rs`) and
//! whose route yFoil is observed to follow — the same call count of every translated subroutine
//! in the step as XFOIL's (`cargo xtask route`, `xtask/fixtures-config/route.toml`). The steps are
//! the cover `cargo xtask steps` chose, listed with the branches each takes in
//! `xtask/fixtures-config/step-cover.toml`, which also lists any branch taken only at a step whose
//! route diverges (`docs/development/conventions/testing.md`, 'What "tested" means for a branch').
//!
//! Fixtures: `tests/fixtures/xfoil/<case>/`, the files each step reads — `cargo xtask fixtures
//! --case <case>` (the cover adds the dumps and the files to the case).

use crate::step::{replay, replay_inviscid};

/// `kt_n60_inviscid_a4_cl05_cl8`, inviscid point 3: 29 branches.
#[test]
fn kt_n60_inviscid_a4_cl05_cl8_point_3() {
    replay_inviscid("kt_n60_inviscid_a4_cl05_cl8", 3);
}

/// `naca0012_n60_a2_re1e6_type2`, VISCAL call 1, iteration 1 (SETBL 1): 835 branches.
#[test]
fn naca0012_n60_a2_re1e6_type2_call_1_iteration_1() {
    replay("naca0012_n60_a2_re1e6_type2", 1, 1, 1);
}

/// `naca0012_n60_a2_re1e6_type3`, VISCAL call 1, iteration 1 (SETBL 1): 833 branches.
#[test]
fn naca0012_n60_a2_re1e6_type3_call_1_iteration_1() {
    replay("naca0012_n60_a2_re1e6_type3", 1, 1, 1);
}

/// `naca0012_n60_a2_repeat_re1e6`, VISCAL call 2, iteration 1 (SETBL 6): 560 branches.
#[test]
fn naca0012_n60_a2_repeat_re1e6_call_2_iteration_1() {
    replay("naca0012_n60_a2_repeat_re1e6", 2, 1, 6);
}

/// `naca0012_n60_a8_re1e4_damp`, VISCAL call 1, iteration 3 (SETBL 3): 539 branches.
#[test]
fn naca0012_n60_a8_re1e4_damp_call_1_iteration_3() {
    replay("naca0012_n60_a8_re1e4_damp", 1, 3, 3);
}

/// `naca0012_n60_polar30_re1e6`, VISCAL call 12, every iteration (SETBL 67–91): 701 branches.
#[test]
fn naca0012_n60_polar30_re1e6_call_12() {
    replay("naca0012_n60_polar30_re1e6", 12, 0, 67);
}

/// `naca16-212_n60_a4_re1e6_m07`, VISCAL call 1, iteration 1 (SETBL 1): 916 branches.
#[test]
fn naca16_212_n60_a4_re1e6_m07_call_1_iteration_1() {
    replay("naca16-212_n60_a4_re1e6_m07", 1, 1, 1);
}

/// `naca16-212_n60_a4_re1e6_m07`, VISCAL call 1, iteration 2 (SETBL 2): 465 branches.
#[test]
fn naca16_212_n60_a4_re1e6_m07_call_1_iteration_2() {
    replay("naca16-212_n60_a4_re1e6_m07", 1, 2, 2);
}

/// `naca23012_n100_a1_a2_re1e6_xtr022_0001`, VISCAL call 2, iteration 2 (SETBL 7): 537 branches.
#[test]
fn naca23012_n100_a1_a2_re1e6_xtr022_0001_call_2_iteration_2() {
    replay("naca23012_n100_a1_a2_re1e6_xtr022_0001", 2, 2, 7);
}

/// `naca23012_n60_a4_re1e6_xtr022_0001`, VISCAL call 1, iteration 4 (SETBL 4): 528 branches.
#[test]
fn naca23012_n60_a4_re1e6_xtr022_0001_call_1_iteration_4() {
    replay("naca23012_n60_a4_re1e6_xtr022_0001", 1, 4, 4);
}

/// `naca2412_n60_a14_re1e6_m05_iter60`, VISCAL call 1, iteration 45 (SETBL 45): 542 branches.
#[test]
fn naca2412_n60_a14_re1e6_m05_iter60_call_1_iteration_45() {
    replay("naca2412_n60_a14_re1e6_m05_iter60", 1, 45, 45);
}

/// `naca4412_n160_a16_re1e6_m05`, VISCAL call 1, iteration 8 (SETBL 8): 552 branches.
#[test]
fn naca4412_n160_a16_re1e6_m05_call_1_iteration_8() {
    replay("naca4412_n160_a16_re1e6_m05", 1, 8, 8);
}

/// `naca4412_n60_cl1_clm05_re1e6`, VISCAL call 1, iteration 1 (SETBL 1): 861 branches.
#[test]
fn naca4412_n60_cl1_clm05_re1e6_call_1_iteration_1() {
    replay("naca4412_n60_cl1_clm05_re1e6", 1, 1, 1);
}

/// `naca4412_n60_cl1_clm05_re1e6`, VISCAL call 2, iteration 1 (SETBL 7): 648 branches.
#[test]
fn naca4412_n60_cl1_clm05_re1e6_call_2_iteration_1() {
    replay("naca4412_n60_cl1_clm05_re1e6", 2, 1, 7);
}

/// `naca63-415_n60_a14_re1e6_m05`, VISCAL call 1, iteration 18 (SETBL 18): 565 branches.
#[test]
fn naca63_415_n60_a14_re1e6_m05_call_1_iteration_18() {
    replay("naca63-415_n60_a14_re1e6_m05", 1, 18, 18);
}

/// `naca64a010_n60_inviscid_type2_m03_am23_a0_a19`, inviscid point 1: 187 branches.
#[test]
fn naca64a010_n60_inviscid_type2_m03_am23_a0_a19_point_1() {
    replay_inviscid("naca64a010_n60_inviscid_type2_m03_am23_a0_a19", 1);
}

/// `naca64a010_n60_type2_m03_a0`, VISCAL call 1, iteration 1 (SETBL 1): 911 branches.
#[test]
fn naca64a010_n60_type2_m03_a0_call_1_iteration_1() {
    replay("naca64a010_n60_type2_m03_a0", 1, 1, 1);
}

/// `naca64a010_n60_type2_m03_a0`, VISCAL call 1, iteration 17 (SETBL 17): 520 branches.
#[test]
fn naca64a010_n60_type2_m03_a0_call_1_iteration_17() {
    replay("naca64a010_n60_type2_m03_a0", 1, 17, 17);
}
