---
icon: lucide/ruler
---

# 6. Geometry

Where XFOIL's geometry handling differs from what its names, or the NACA definitions, suggest.

| | Issue | Description | Handling |
|---|---|---|---|
| [6.1](6-1-naca4-applies-thickness-vertically.md) | [NACA4 applies thickness vertically](6-1-naca4-applies-thickness-vertically.md) | XFOIL's NACA command adds thickness vertically rather than perpendicular to the camber line. | Overcome (the one live divergence) |
| [6.2](6-2-load-reverses-clockwise-input.md) | [LOAD reverses clockwise input](6-2-load-reverses-clockwise-input.md) | LOAD reverses an aerofoil given clockwise. | Constraint |
| [6.3](6-3-chord-and-apanel-are-not-what-their-names-suggest.md) | [CHORD and APANEL are not what their names suggest](6-3-chord-and-apanel-are-not-what-their-names-suggest.md) | CHORD is not exactly the chord, and APANEL is a panel-normal angle. | Replicated |
| [6.4](6-4-nw-n-12-10int-waklen.md) | [NW = N/12 + 10·INT(WAKLEN)](6-4-nw-n-12-10int-waklen.md) | The number of wake nodes is computed with integer truncation. | Replicated |
| [6.5](6-5-tgap-on-a-sharp-trailing-edge-does-not-deliver-the-requested-gap.md) | [TGAP on a sharp trailing edge does not deliver the requested gap](6-5-tgap-on-a-sharp-trailing-edge-does-not-deliver-the-requested-gap.md) | TGAP on a sharp trailing edge opens a smaller gap than the one asked for. | Replicated |
| [6.6](6-6-naca-generator-constants-the-ordinate-program-s-precision.md) | [NACA generator constants: the ordinate program's precision](6-6-naca-generator-constants-the-ordinate-program-s-precision.md) | The NACA ordinate program yFoil is checked against has limited precision of its own. | Reference note |
