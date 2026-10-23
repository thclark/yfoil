---
icon: lucide/infinity
---

# The non-finite set

These are the tracked candidates whose reference run printed a non-finite value. They are not candidates
for the finite set: a non-finite trajectory cannot be compared as a whole run, because NaN comparisons
and `MIN`/`MAX` are where two correct codes can legitimately differ. The branches that only these runs
take are tested one call at a time instead, as events.

--8<-- "validation/branch-coverage/tables/non-finite-set.md"

## Branches only these runs take

The non-finite set takes these branches and the finite set does not. The last column names the test that replays each one: as an event, seeding yFoil with XFOIL's state dumped entering the SETBL call at which `events.dat` records the event and requiring the same event and XFOIL's state after it, or, for TRCHEK2's `AX ≤ 0` exit, each TRCHEK2 call the reference recorded taking it.

--8<-- "validation/branch-coverage/tables/branches-non-finite.md"

## How far yFoil follows each run

This is a measurement, not a claim — the reference cannot reproduce itself on these runs — of how far yFoil follows each one. The comparison above is applied call by call: *parts at* is the first iteration at which the runs were allowed to part (the reference's own 1-ULP spread exceeds `DIVERGENCE_FLOOR` there and the values differ) or, for a mismatch, the iteration after the last one compared; *compared* counts the iterations matched before that. The NaN count is the number of `NaN` tokens in the reference's own output for the run.

--8<-- "validation/branch-coverage/tables/non-finite-comparison.md"
