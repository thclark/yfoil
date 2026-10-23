---
icon: lucide/list-todo
---

# 8. Open items

- 8.1 `src/bin/generate_subroutine_validation.rs:691-695` still emits an `E24.16` instrumentation
  template (16 significant figures) where the rule is `ES24.16`. Harmless today because nothing
  tests against that generator's output, but it contradicts §1.5.
- 8.2 Two `#[ignore]`d compressible closure unit tests (§2.1, §7.3) await regeneration from the M = 0.3
  coverage case (`tests/IGNORED.txt`).

- 8.3 About twenty unit tests in `src/` (`bl/transition.rs`, `bl/difference.rs`, `bl/station.rs`,
  `bl/params.rs`) compare closures with hard-coded values of unknown provenance, three of them
  `#[ignore]`d (`tests/IGNORED.txt`). The same closures are held to `TOL_PURE` against the reference's
  own log, incompressible and at M = 0.3, by `tests/subroutine/closures.rs`; the unit tests are to be
  moved onto those fixtures (`docs/development/conventions/testing.md`, Layout).
- 8.4 `tests/subroutine/pane_legacy.rs` compares the curvature repanelling with panels XFOIL generated
  itself (`tests/fixtures/naca0012/panels*.json`, contrary to CLAUDE.md Rule 4); PANGEN is tested against
  yFoil-generated panels by `tests/subroutine/pangen.rs`.
- 8.5 One step of the branch candidates is a divergent comparison
  (`xtask/fixtures-config/route.toml`: NACA 63-415 M 0.5, SETBL 2, MRCHDU's rejected stations differ);
  its branches are tested by other steps. The step tests compare a step's outputs, not the route
  through it; observing the route is the next stage (`docs/development/conventions/testing.md`, 'What "tested" means for a branch').
- 8.6 The fixture runner's watchdog used to end a reference run after five minutes of wall-clock, so
  where it cut a run still going depended on the machine's load: regenerated from scratch on
  2026-10-07, the twins study reproduced 7 of its 11 N = 240 cases exactly, and the 4 the cap cut
  differed (the 23012 at 2 recorded points instead of 37). The clock no longer ends a run: a run ends
  on its own, at XFOIL's plot-label hang (120 s without output, a fixed point of the run), or where a
  case asks for a deterministic stop (`stop_on_nonfinite`, `stop_after_setbl`, instrumentation
  patch 20). The twins study was removed on 2026-10-08, superseded by the twins every validation run
  carries at its own panelling; its slow cases live on as the `pathological` copies, which ask for no
  stop and run only on request.
- 8.7 `xtask/fixtures-config/coverage.toml`'s `[[unreachable]]` annotations are code-reading arguments,
  and one was wrong: TRCHEK2's `IF(AX .LE. 0.0) GO TO 101` (`xblsys.f:390`) was annotated
  structurally unreachable, on the argument that AX = AXA + DAX cannot be ≤ 0 while the momentum
  thicknesses stay positive. A run that has gone non-finite takes it: found 2026-10-07 by comparing
  the twins-study cases' branches with the tracked cases' (the N = 240 23012 sweep with trips), and
  reproduced by `naca23012_n100_a1_a2_re1e6_xtr022_0001` at call 2, iteration 2 (SETBL 7), the
  iteration in which the boundary layer turns non-finite; reached from the 1° solution only, not by a
  fresh march at 2°, and at N ≥ 100 only. The annotation is removed and the branch is tested (a
  cover step, `tests/execution/branches.rs`, and every iteration of the call, `tests/execution/steps.rs`).
  Those replays require XFOIL's route and values, but nothing in them names this exit, so the
  reference also records every TRCHEK2 call that takes it (`trchek2_ax_exit.dat`, instrumentation
  patch 23: 28 calls at SETBL 7, all with finite inputs) and
  `tests/subroutine/transition.rs::test_trchek2_no_amplification_exit_matches_xfoil` replays each one,
  requiring the same exit at the same iterate with AX, AMPL2 and XT within `TOL_PURE`.
  A wrong `[[unreachable]]` hides a reachable branch from the coverage search, so every annotation is
  to be proved from the source's syntax tree rather than argued (open).
- 8.8 **Closed-trailing-edge 6-series sections do not converge in a routine viscous sweep** (2026-10-08). The
  NACA 63-415 and 64A010 from yFoil's generators have a trailing-edge gap of exactly zero; at N = 160
  (cosine sampling), ITER 100, Re 1e6, M 0, TYPE 1, neither converges at any point of a polar — the
  first call's RMSBL is about 3.5e3 from its first iteration and MRCHUE fails at the trailing-edge
  station of both sides — where every blunt 4- and 5-digit section converges. Whether this is XFOIL's
  sharp-TE path, the cosine sampling's trailing-edge panels or the generators' closure is open; the
  series has no 6- or 6A-series member until it is settled. Both runs are known-issues cases
  (`naca63-415_n160_polar30_re1e6_iter100`, `naca64a010_n160_polar30_re1e6_iter100`).
- 8.9 The `pathological` cases (§9) are to be studied: why each is slow or hangs, why the 4412 sweeps
  hang with a finite boundary layer, and whether yFoil does the same — the first part of a study
  that replicates every known issue with its yFoil equivalent (open).
- 8.10 The tolerances study, `studies/tolerances/tolerances.sh` (one all-coordinates +1-ULP twin of the reference case, the origin of the
  tolerance constants in `tests/common/utilities/tolerances.rs`) predates the seeded twins and moves
  every coordinate the same way, which `docs/development/conventions/terminology.md` notes measures little. The
  reference case's own five twins, which `cargo xtask fixtures` can now run (`twins = true`), are to
  derive `docs/validation/tolerances/index.md` instead, after which the script can be retired (open).

## Runs that show it

- [NACA 63-415, closed trailing edge](../runs/naca-63-415-closed-trailing-edge.md)
- [NACA 64A010, closed trailing edge](../runs/naca-64a010-closed-trailing-edge.md)
