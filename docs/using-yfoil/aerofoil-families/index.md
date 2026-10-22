---
icon: lucide/layers
---

# Aerofoil families

What a designation produces. yFoil generates every NACA family that has a public, reproducible
definition, and one analytic section; this section gives each definition, shows each family
across its parameters, and then covers what happens between a definition and the nodes an
analysis runs on — how the thickness is laid, how the nodes are distributed, what is done to the
trailing edge, and the record a generated file keeps of all of it.

Whether the generators reproduce the definitions is checked separately, against NASA's
public-domain ordinate program: [generated sections against
naca456](../../engineering-validation/generated-sections-against-naca456.md).

## How a NACA section is built

A NACA section is a **basic thickness form** \(y_t(x)\) laid **perpendicular** to a
**mean line** \(y_c(x)\) ([abbott1959](../bibliography.md#abbott1959), §6.2):

\[
x_u = x - y_t \sin\theta,\quad y_u = y_c + y_t \cos\theta,\qquad
x_l = x + y_t \sin\theta,\quad y_l = y_c - y_t \cos\theta,\qquad
\theta = \arctan\frac{\mathrm{d}y_c}{\mathrm{d}x}
\]

The formulas on these pages are those of the NASA ordinate program
([ladson1996](../bibliography.md#ladson1996)) and its public-domain revision `naca456`
([carmichael2001](../bibliography.md#carmichael2001)). Every generated file carries a
[provenance record](provenance-of-a-generated-section.md) naming the family, the parameters and
these references.

## The families

| Page | Designation | Thickness form | Mean line | Command |
| --- | --- | --- | --- | --- |
| [NACA 4-digit and 4-digit modified](naca-4-digit-and-4-digit-modified.md) | `2412`, `0012-34` | 4-digit, 4-digit modified | 2-digit | `yfoil geometry naca 2412` |
| [NACA 5-digit, including reflex](naca-5-digit-including-reflex.md) | `23012`, `23112` | 4-digit | 3-digit, reflex | `yfoil geometry naca 23012` |
| [NACA 16-series](naca-16-series.md) | `16-212` | 4-digit modified (I = 4, M = 0.5) | 6-series, a = 1 | `yfoil geometry naca 16-212` |
| [NACA 6- and 6A-series](naca-6-and-6a-series.md) | `63-415`, `64A010` | 63…67, 63A…65A (tabulated) | 6-series a, 6A | `yfoil geometry naca 63-415 [--a 0.5]` |
| [Kármán–Trefftz sections](karman-trefftz-sections.md) | — | conformal map | — | `yfoil geometry karman-trefftz` |

Each family page draws the family across its parameters: one panel per parameter, the family's
default section in black on top, the sections that vary the parameter behind it coloured by value
(YlGnBu, lightest = smallest), 160 nodes each, every panel on the same linked axes (x from −0.2 to
1.2, y from −0.5 to 0.5, at equal scale). The figures are drawn by the geometry validation study
(`cargo run --release -p geometry-validation -- --docs`).

All of the sections are sampled at the same cosine stations in \(x\) (\(n/2\) per surface, a node
exactly at the trailing edge, the leading edge straddled) before they are panelled, so sections
of different families are compared on the same node distribution.

## From a definition to the nodes

| Page | What it covers |
| --- | --- |
| [Thickness perpendicular to the camber line](thickness-perpendicular-to-the-camber-line.md) | Why yFoil lays the thickness as the NACA definition does, unlike XFOIL's `NACA` command, and how to get XFOIL's version |
| [Panel distributions](panel-distributions.md) | XFOIL's `PANGEN` and yFoil's cosine spacing, and the parameters that shape each |
| [Trailing-edge gap and sharpening](trailing-edge-gap-and-sharpening.md) | Which sections close, and opening or closing a trailing edge after panelling, as XFOIL's `TGAP` does |
| [Provenance of a generated section](provenance-of-a-generated-section.md) | The `generator` record that says how a geometry file was made |

The flags themselves are on [making and converting geometry](../commands/making-and-converting-geometry.md).
