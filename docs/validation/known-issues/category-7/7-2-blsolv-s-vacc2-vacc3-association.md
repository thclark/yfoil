# 7.2 BLSOLV's VACC2/VACC3 association — Replicated

`VACC2 = VACC3 = (VACCEL*2.0)/(S(N)-S(1))` decides the sparse-elimination skips; a 1-ULP difference in the
threshold flips a branch (`src/bl/blsolv.rs:122-124`). A historic "~1 % BLSOLV error" was a yFoil test
hard-coding `S(N)−S(1) = 2.0` instead of the fixture's 2.0387 (`tests/subroutine/blsolv.rs`);
BLSOLV is bit-identical.
