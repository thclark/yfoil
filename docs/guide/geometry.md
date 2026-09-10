---
icon: lucide/spline
---

# Geometry

``` sh
yfoil geometry <SUBCOMMAND>   # alias: yfoil geom
```

Four subcommands: [`naca`](#naca), [`convert`](#convert), [`repanel`](#repanel)
and [`info`](#info).

## Conventions

- Coordinates are normalised by chord.
- Panel ordering is **trailing edge → upper surface → leading edge → lower
  surface → trailing edge**, counter-clockwise.
- Written `.dat` files carry 17 significant figures, so a JSON → `.dat` → JSON
  round trip is bit-exact.

## naca

Generate a 4- or 5-digit NACA section.

``` sh
yfoil geometry naca 0012 -n 160 -o naca0012.json
yfoil geometry naca 4412 -n 200 --sharp -o naca4412.json
yfoil geometry naca 23015 -o naca23015.json
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-n`, `--panels` | 160 | Number of panels |
| `--to` | `json` | Output format: `json` or `dat` |
| `--sharp` | off | Close the trailing edge (zero TE gap) |
| `--thickness` | `perpendicular` | `perpendicular` or `vertical` — see below |
| `-o`, `--output` | stdout | Output file |

!!! info "Thickness distribution"

    yFoil applies the thickness distribution **perpendicular to the camber
    line**, which is the NACA definition. XFOIL's `NACA4` applies it vertically;
    that variant is available as `--thickness vertical` for comparison.

    This is the only deliberate numerical divergence from XFOIL, and it does not
    affect validation, because yFoil generates the panels and XFOIL consumes
    them.

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
