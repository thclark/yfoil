# 5.1 Subroutines never called — Not translated

`UECALC` and `DSSET` (`xpanel.f`) are never called anywhere in 6.99, and `DIT`'s only call site is
commented out (`xblsys.f:944`). None of the three is translated: they were carried for completeness
until 2026-09-13 and removed once the branch-coverage study confirmed no case can reach them, so the
translated set of `xtask/fixtures-config/coverage.toml` no longer names them.
