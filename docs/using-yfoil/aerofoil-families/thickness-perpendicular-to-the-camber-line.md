# Thickness perpendicular to the camber line

yFoil lays the thickness form **perpendicular to the camber line**, as the NACA definition does
([how a NACA section is built](index.md#how-a-naca-section-is-built)), whatever the panelling
method. XFOIL's own `NACA` command does not: its `NACA4` and `NACA5` add the thickness
*vertically*, \(y = y_c \pm y_t\). For a symmetric section the two agree, since the camber slope is
zero; for a cambered one they are different sections, and the difference grows with the camber.

## Getting XFOIL's version

`--thickness vertical` reproduces XFOIL's `NACA` command for the 4- and 5-digit families:

``` sh
yfoil geometry naca 4412 -n 160 --thickness vertical -o naca4412_xfoil.json
```

It uses XFOIL's own 245-point buffer and is always panelled by `PANGEN`, because it exists only to
replicate what XFOIL's `NACA` command produces; the geometry file records
`thickness_applied: "vertical"`. It reproduces XFOIL's output to the last bit or two
([XFOIL's own geometry routines](../../fidelity-to-xfoil/subroutine-equivalence/xfoils-own-geometry-routines.md)).

## Why this does not affect the comparison with XFOIL

Every comparison of yFoil with XFOIL runs both codes on the same panels: yFoil generates them and
XFOIL loads them unchanged ([the bitwise geometry
handoff](../../fidelity-to-xfoil/fixtures.md#the-bitwise-geometry-handoff)). So the choice of
thickness application never enters a comparison. It is nonetheless the only place yFoil
knowingly departs from XFOIL, and it is recorded, with the evidence, among the
[deliberate divergences](../../fidelity-to-xfoil/xfoils-known-issues/deliberate-divergences.md).
