# 7.5 EQUIVALENCE aliasing of UNEW/QNEW onto VA/VB — Overcome (structure only)

UPDATE's `UNEW`/`QNEW` alias the factored `VA`/`VB` blocks (`xbl.f`), which are dead after BLSOLV. yFoil
uses separate locals and a consuming `BlsolvInput`, same numerics, so the alias cannot be read by
accident (`src/bl/blsolv.rs:92-94`, `src/solver/update.rs:6-8`).
