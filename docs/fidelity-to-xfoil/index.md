---
icon: lucide/microscope
---

# Fidelity to XFOIL

This tab is for developers, and for anyone who wants the evidence that yFoil is a faithful
**translation** of XFOIL 6.99 — not a reimplementation, not "something like XFOIL", but the same
computation, statement by statement, to within rounding. Whether its results are fit for
engineering work is the previous tab, [Engineering validation](../engineering-validation/index.md);
this one asks a narrower and stricter question, against a stricter comparison target: an
instrumented build of XFOIL itself.

## The argument, in order

Each page builds on the ones before it, so they are listed in the order the argument runs.

1. **What is compared against.** [The reference XFOIL build](the-reference-xfoil-build/index.md):
   XFOIL in double precision with pinned flags and correctly rounded maths, instrumented to dump its
   state, and proved sound.
2. **What is captured from it.** [Fixtures](fixtures.md): the reference runs the tests read, on
   panels yFoil generated and XFOIL loaded bitwise.
3. **How close counts as equal.** [Tolerances](tolerances/index.md): derived from the reference's
   own noise floor, with the branch trace compared exactly.
4. **Each routine alone.** [Subroutine equivalence](subroutine-equivalence/index.md): one
   translated subroutine on XFOIL's own inputs.
5. **The routines together.** [Execution equivalence](execution-equivalence/index.md): whole runs
   where the solution is well-conditioned, single replayed steps where it is not.
6. **Over every branch.** [Branch coverage](branch-coverage/index.md): every reachable branch of
   the translated subroutines taken, and each tested by a step on XFOIL's route.

Two reference sections support all of them:

- [XFOIL's known issues](xfoils-known-issues/index.md) — the quirks and bugs a faithful
  translation reproduces, the one deliberate divergence, and runs that show the issues arising;
- [Translation map](translation-map/index.md) — where each XFOIL variable and routine lives in
  yFoil, and how XFOIL is put together.

## The rule the tab enforces

A comparison passes only when the branch trace is identical **and** the values are within
tolerance. Under those conditions, a difference larger than the tolerance is a bug — 1 %, 0.1 %,
1e-6 alike; there is no "acceptable engineering tolerance". Where the two codes take different
routes because a branch's operands are within rounding of its threshold, the comparison says
nothing about the translation, so no test is one; those solutions are tested one step at a time
instead. The words used for all this — *divergent comparison*, *ill-conditioned*, *twin*,
*noise floor* — have exact meanings, defined in the
[terminology](../contributing/conventions/terminology.md).
