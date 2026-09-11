---
icon: lucide/spline
---

# Geometry

``` sh
yfoil geometry <SUBCOMMAND>   # alias: yfoil geom
```

Five subcommands: [`naca`](#naca), [`karman-trefftz`](#karman-trefftz),
[`convert`](#convert), [`repanel`](#repanel) and [`info`](#info).

## Conventions

- Coordinates are normalised by chord.
- Panel ordering is **trailing edge → upper surface → leading edge → lower
  surface → trailing edge**, counter-clockwise.
- Written `.dat` files carry 17 significant figures, so a JSON → `.dat` → JSON
  round trip is bit-exact.
- A generated geometry carries a `generator` record: the series, designation
  and parameters it was made from (see [the data model](data.md#geometry)).

## naca

Generate a NACA section from its designation. Every family with a public,
reproducible definition is covered; the formulas, what each family is for and
the references are in [aerofoil series](aerofoil-series.md).

``` sh
yfoil geometry naca 0012 -n 160 -o naca0012.json
yfoil geometry naca 4412 -n 200 --sharp -o naca4412.json
yfoil geometry naca 23018 -o naca23018.json
yfoil geometry naca 0012-34 -o naca0012-34.json
yfoil geometry naca 16-212 -o naca16-212.json
yfoil geometry naca 63-415 --a 0.5 -o naca63-415-a05.json
yfoil geometry naca 64A010 --te-gap 0.002 -o naca64a010.json
```

| Designation | Series | Parameters |
| --- | --- | --- |
| `2412` | 4-digit | camber `2`%, at `4` tenths, thickness `12`% |
| `0012-34` | 4-digit modified | as above, then leading-edge radius index `3`, maximum thickness at `4` tenths |
| `23012`, `23112` | 5-digit | design \(C_L\) `2` × 0.15, camber position `3` × 0.05, `0` standard or `1` reflex, thickness `12`% |
| `16-212` | 16-series | design \(C_L\) `2` tenths, thickness `12`% |
| `63-415` | 6-series | family `63`, design \(C_L\) `4` tenths, thickness `15`%; `--a` sets the loading extent |
| `64A010` | 6A-series | family `64A`, design \(C_L\) `0`, thickness `10`% |

An optional `NACA ` prefix and letter case are ignored, as is a bracketed
low-drag-range subscript (`64(1)-212`).

| Option | Default | Meaning |
| --- | --- | --- |
| `-n`, `--panels` | 160 | Number of panels |
| `--to` | `json` | Output format: `json` or `dat` |
| `--a` | 1.0 | Extent of uniform loading of a 6-series or 16-series mean line, 0…1 |
| `--sharp` | off | Close the trailing edge (move both TE nodes to their midpoint) |
| `--te-gap` | — | Set the trailing-edge gap (chord units) with XFOIL's `TGAP` blending |
| `--te-blend` | 1.0 | Blending distance/c of `--te-gap`, 0…1 (`TGAP`'s second argument) |
| `--thickness` | `perpendicular` | `perpendicular` or `vertical` — see below |
| `-o`, `--output` | `naca<designation>.<format>` | Output file |

!!! info "Thickness distribution"

    yFoil applies the thickness distribution **perpendicular to the camber
    line**, which is the NACA definition. XFOIL's `NACA4`/`NACA5` apply it
    vertically; that variant is available as `--thickness vertical` for
    4- and 5-digit sections, for comparison.

    This is the only deliberate numerical divergence from XFOIL, and it does not
    affect validation, because yFoil generates the panels and XFOIL consumes
    them.

!!! info "Trailing edges"

    The 4-digit families have a finite trailing-edge thickness; the 6-series,
    6A-series and Kármán–Trefftz sections close exactly. XFOIL does not blunt a
    sharp trailing edge (it detects one and switches branches), so yFoil
    generates each section as defined. `--te-gap` is XFOIL's `TGAP` for those
    who want a gap; on a closed edge XFOIL, and therefore yFoil, delivers
    \(\cos(\tau/2)\) times the requested gap
    ([known issues §6.5](../xfoil-known-issues.md)).

## karman-trefftz

Generate a Kármán–Trefftz section: the conformal map of a circle through
\(\zeta = 1\), an analytic section with an exact potential-flow solution and a
sharp trailing edge of angle \(\tau\) (\(\tau = 0\) is the Joukowski cusp).

``` sh
yfoil geometry karman-trefftz --x-centre -0.1 --y-centre 0.05 --te-angle 10 -o kt.json
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--x-centre` | −0.1 | Circle centre x (negative; sets the thickness) |
| `--y-centre` | 0.05 | Circle centre y (sets the camber) |
| `--te-angle` | 10 | Trailing-edge angle in degrees, 0 ≤ τ < 180 |
| `-n`, `--panels` | 160 | Number of panels |
| `--to` | `json` | Output format: `json` or `dat` |
| `--te-gap`, `--te-blend` | — | As for `naca` |
| `-o`, `--output` | `karman-trefftz.<format>` | Output file |

## convert

``` sh
yfoil geometry convert naca0012.json --to dat -o naca0012.dat
yfoil geometry convert e387.dat --to json -o e387.json
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--to` | `json` | Output format: `json` or `dat` |
| `--name` | `Aerofoil` | Name written into the `.dat` header |
| `-o`, `--output` | stdout | Output file |

## repanel

``` sh
yfoil geometry repanel e387.json -n 160 --method curvature -o e387_repanelled.json
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-n`, `--panels` | 160 | Target panel count |
| `--method` | `curvature` | `curvature` (XFOIL's `PANE`) or `cosine` |
| `--te-le-ratio` | 0.15 | TE/LE panel density ratio (`cosine` only) |

## info

``` sh
yfoil geometry info naca4412.json
yfoil geometry info naca4412.json -o info.json
```

Prints a summary of the section — with `-o`, writes it as JSON instead.
