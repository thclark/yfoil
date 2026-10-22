# Deliberate divergences

yFoil never modifies XFOIL's algorithm: no relaxation factors, iteration limits, convergence
bail-outs, clamps or safeguards that XFOIL does not have, no simplified formulas, no changed signs
or coefficients, no coupling term disabled "temporarily". If XFOIL does X, yFoil does X — including
XFOIL's own bugs and quirks, which is why [the issue register](the-issue-register.md) marks most
of them *Replicated*. A `no-deviations` check in CI greps for the tell-tale vocabulary.

A deliberate divergence from XFOIL is allowed only when it is listed here, and adding one requires
the same standard of evidence as a bug report against XFOIL. There is one.

| Divergence | Why | Where |
|---|---|---|
| NACA 4- and 5-digit thickness applied perpendicular to the camber line | XFOIL's `NACA4` (`naca.f:62`) applies it vertically, which is not the NACA definition ([known issues §6.1](the-issue-register.md)). Irrelevant for the comparison with XFOIL, because yFoil generates the panels and XFOIL consumes them ([the bitwise geometry handoff](../fixtures.md#the-bitwise-geometry-handoff)). XFOIL's variant is available as `--thickness vertical`. | `src/geometry/naca.rs` |

What the divergence means for a user, and how to get XFOIL's sections, is
[thickness perpendicular to the camber line](../../using-yfoil/aerofoil-families/thickness-perpendicular-to-the-camber-line.md).

The other entries marked *Overcome* in the register are not divergences in this sense: they are
on the output side (yFoil writes a quantity XFOIL prints wrongly or imprecisely) or in the
instrumentation, and no solver value differs from XFOIL's because of them.
