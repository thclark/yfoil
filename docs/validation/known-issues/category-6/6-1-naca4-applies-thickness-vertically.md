# 6.1 NACA4 applies thickness vertically — Overcome (the one live divergence)

`naca.f:62`: `YB(IB) = YC(I) + YT(I)`, i.e. thickness added vertically rather than perpendicular to
the camber line, which is not the NACA definition (irrelevant for symmetric sections). yFoil's
generators lay the thickness perpendicular to the mean line by default, whatever the panelling
method (`yfoil geometry naca …`, `Section`). XFOIL's model is `--thickness vertical`
(`naca_4digit_vertical`, `naca_5digit_vertical`): XFOIL's own 245-point NACA4/NACA5 buffer, always
PANGEN-panelled, and it exists only to replicate the output of XFOIL's `NACA` command bitwise
(`tests/subroutine/pangen.rs`); the geometry file records `thickness_applied: "vertical"`. This
never enters a comparison because yFoil generates the panels and XFOIL consumes them (Rule 4).
Recorded in CLAUDE.md's divergence table, the only row.
