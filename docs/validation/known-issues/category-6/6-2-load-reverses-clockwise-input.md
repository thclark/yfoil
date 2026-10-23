# 6.2 LOAD reverses clockwise input — Constraint

`xfoil.f:1251-1260`: a negative signed area flips the node order. yFoil writes counter-clockwise
(TE → upper → LE → lower → TE) so the bitwise handoff holds; `LNORM` defaults to false
(`xfoil.f:495`) so LOAD does not rescale.
