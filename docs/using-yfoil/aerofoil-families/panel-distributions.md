# Panel distributions

A section is a curve; an analysis runs on a finite set of nodes on it. How the nodes are placed
matters to the answer — XFOIL's results can depend on the panel count in ways that are not small
([known issues §7.13](../../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md)) — so yFoil
offers XFOIL's own method and a simple analytic one, and records which was used. The flags are on
[making and converting geometry](../commands/making-and-converting-geometry.md#repanel).

## PANGEN, XFOIL's method

**`pangen` (default)** is XFOIL's `PANGEN`, translated line for line and tested
against XFOIL's own output ([XFOIL's own geometry
routines](../../fidelity-to-xfoil/subroutine-equivalence/xfoils-own-geometry-routines.md)).
It splines the input, forms the curvature along the arc length, smooths it, adds
a *fictitious* curvature at the trailing edge and inside optional refinement
windows, and then places the nodes so that
\((1 + 6\,\mathrm{CVPAR}\cdot\kappa)\,\Delta s\) is the same on every
panel. Its parameters are XFOIL's `PPAR` menu, here by descriptive names:

| `PPAR` key | XFOIL variable | Flag | Default | Meaning |
| --- | --- | --- | --- | --- |
| `N` | `NPAN` | `-n`, `--nodes` | 160 | Number of panel nodes, trailing edge round to trailing edge (N − 1 surface panels plus the trailing-edge panel) |
| `P` | `CVPAR` | `--curvature-bunching` | 1.0 | Curvature attraction; 0 gives uniform arc-length spacing |
| `T` | `CTERAT` | `--te-curvature-ratio` | 0.15 | Fictitious trailing-edge curvature as a fraction of the leading-edge curvature ("TE/LE panel density ratio") |
| `R` | `CTRRAT` | `--refined-curvature-ratio` | 0.2 | Fictitious curvature inside the refinement windows as a fraction of the leading-edge curvature |
| `XT` | `XSREF1`, `XSREF2` | `--refine-upper X1,X2` | off | Upper-surface refinement window in \(x/c\) |
| `XB` | `XPREF1`, `XPREF2` | `--refine-lower X1,X2` | off | Lower-surface refinement window in \(x/c\) |

A generated section is sampled first at `n_buffer_nodes` stations (246 by default, XFOIL's own
buffer density) and PANGEN runs on those, as XFOIL's would on a loaded buffer aerofoil.

## Cosine spacing, yFoil's method

**`cosine`** is yFoil's own method and has no XFOIL equivalent: cosine spacing in
arc length on each surface, its parameter warped by a power law set by
`--cosine-te-bias` (1 is a plain cosine; below 1 coarser at the trailing edge
and finer at the leading edge; above 1 finer at the trailing edge; default 0.15,
clamped to 0.05…2). It places a node at the leading edge, so the upper surface
carries one node more than the lower. For a generated section, `--method cosine` uses the
generator's analytic cosine sampling in \(x\) directly, and no bias applies; this is the sampling
the tracked XFOIL comparisons were made with.

Next: what can be done to the [trailing edge](trailing-edge-gap-and-sharpening.md) once the nodes
are placed.
