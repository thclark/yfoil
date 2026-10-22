# Trailing-edge gap and sharpening

## Which sections close

The 4-digit families keep a finite trailing-edge thickness; the 6-series, 6A-series and
Kármán–Trefftz sections close exactly. XFOIL does not blunt a sharp trailing edge anywhere:
`PANGEN` keeps the buffer's end points, `TECALC` detects sharpness geometrically
(`SHARP = DSTE < 10⁻⁴ CHORD`) and the panel method, the wake and the boundary layer switch
branches on it, all of which yFoil translates and tests. yFoil therefore generates each section
as defined.

## Closing or opening a trailing edge

Two adjustments are available, both recorded in the [provenance
record](provenance-of-a-generated-section.md), and they are exclusive:

- `--sharp` moves the two trailing-edge nodes to their midpoint (a closed edge, XFOIL's `SHARP`
  path; a closed 4-digit section, for instance);
- `--te-gap GAP [--te-blend F]` is XFOIL's `TGAP`: the surfaces are moved apart along the
  existing gap's direction (or, for a closed edge, the mean end tangent) by
  \(\tfrac{1}{2}\Delta\,(x/c)\,e^{-(1 - x/c)(1/F - 1)}\) each, blended over the fraction \(F\)
  of the chord from the trailing edge.

Both are applied **after** the nodes are distributed, so the blend profile is evaluated exactly at
every output node rather than splined through a coarse input; XFOIL applies `TGAP` to the buffer
aerofoil and `PANE` follows, and the two orders differ by that interpolation.

On a closed edge XFOIL delivers \(\cos(\tau/2)\) times the requested gap, \(\tau\) the
trailing-edge angle, because the direction it moves the surfaces along is not of unit length
([known issues §6.5](../../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md)); yFoil
reproduces that, so a `--te-gap 0.002` on a NACA 63-415 gives a gap of 2.0012e-3.

``` sh
yfoil geometry naca 4412 -n 160 --sharp -o naca4412_sharp.json
yfoil geometry naca 64A010 --te-gap 0.002 -o naca64a010_gap.json
```
