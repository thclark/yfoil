# 5.6 Structurally dead branches — Out of scope

Coincident/duplicate consecutive nodes (ABCOPY deletes them on LOAD, PANGEN never produces them):
`xpanel.f:29, 66, 77, 81, 180, 861, 867, 875`, `spline.f:542-543`. SPLIND's `999`/specified-slope end
conditions (only SEGSPL's `-999/-999` is reached on the analysis path): `spline.f:96, 101, 113, 117`.
Array-bound `STOP`s that yFoil, having no fixed dimensions, cannot hit. Interactive prompts (`NITER = 0`,
`NACA` without a designation). Flap hinge moments (`LFLAP`, `LBFLAP`). `LIPAN` cleared with `LBLINI`
kept (`xqdes.f:497` only). The impossible flag combinations `LADIJ = F ∧ LWDIJ = T` and
`LQAIJ = F ∧ LGAMU = T`. Each is an `[[unreachable]]` entry with a class and reason in `coverage.toml`;
84 in total.
