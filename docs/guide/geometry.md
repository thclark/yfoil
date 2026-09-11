---
icon: lucide/spline
---

# Geometry

``` sh
yfoil geometry <SUBCOMMAND>   # alias: yfoil geom
```

Five subcommands: [`naca`](#naca), [`karman-trefftz`](#karman-trefftz),
[`convert`](#convert), [`repanel`](#repanel) and [`info`](#info). Every one of
them prints its full help when run without arguments; `--help` gives the long
form.

## Conventions

- Coordinates are normalised by chord.
- Panel ordering is **trailing edge → upper surface → leading edge → lower
  surface → trailing edge**, counter-clockwise.
- Written `.dat` files carry 17 significant figures, so a JSON → `.dat` → JSON
  round trip is bit-exact.
- A generated or repanelled geometry carries a `generator` record: the series,
  designation and parameters it was made from, and how it was panelled (see
  [the data model](data.md#geometry)).

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
yfoil geometry naca 63-415 -n 160 --te-curvature-ratio 0.25 -o naca63-415.json   # PANGEN, tuned
yfoil geometry naca 63-415 -n 160 --method cosine -o naca63-415_cosine.json      # analytic sampling
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
| `--to` | `json` | Output format: `json` or `dat` |
| `--a` | 1.0 | Extent of uniform loading of a 6-series or 16-series mean line, 0…1 |
| `--thickness` | `perpendicular` | `perpendicular` or `vertical` — see below |
| `-n`, `--method`, PPAR flags, `--sharp`, `--te-gap`, `--te-blend`, `--panelling` | as [`repanel`](#repanel) | The panelling: `pangen` (default) samples the section at 246 stations and runs XFOIL's PANGEN on them; `cosine` is the analytic cosine sampling in \(x\) (no bias applies) |
| `-o`, `--output` | `naca<designation>.<format>` | Output file |

!!! info "Thickness distribution"

    yFoil applies the thickness distribution **perpendicular to the camber
    line**, which is the NACA definition, and does so whatever the panelling
    method. XFOIL's `NACA4`/`NACA5` apply it vertically
    ([known issues §6.1](../xfoil-known-issues.md)); that model is available as
    `--thickness vertical` for 4- and 5-digit sections, on XFOIL's own 245-point
    buffer and always PANGEN-panelled, only for replicating XFOIL's `NACA`
    command output.

    This is the only deliberate numerical divergence from XFOIL, and it does not
    affect validation, because yFoil generates the panels and XFOIL consumes
    them.

!!! info "Trailing edges"

    The 4-digit families have a finite trailing-edge thickness; the 6-series,
    6A-series and Kármán–Trefftz sections close exactly. XFOIL does not blunt a
    sharp trailing edge (it detects one and switches branches), so yFoil
    generates each section as defined; `--sharp` and `--te-gap` are the
    adjustments, applied after panelling as under [`repanel`](#repanel).

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
| `--to` | `json` | Output format: `json` or `dat` |
| `-n`, `--method`, PPAR flags, `--sharp`, `--te-gap`, `--te-blend`, `--panelling` | as [`repanel`](#repanel) | The panelling, `pangen` by default |
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

Redistribute the nodes of an existing geometry, the operation XFOIL's `PANE` and
`PPAR` commands perform on a loaded aerofoil. The input is a `.json` or `.dat`
geometry; the output is JSON, `<input stem>_repanelled.json` beside the input
unless `-o` says otherwise, and it carries the panelling under `generator`.

``` sh
yfoil geometry repanel                                  # prints the full help
yfoil geometry repanel e387.dat -n 160                  # XFOIL's PANGEN, XFOIL's defaults
yfoil geometry repanel e387.dat -n 200 --curvature-bunching 1.5 --te-curvature-ratio 0.3 \
    --refined-curvature-ratio 0.5 --refine-upper 0.2,0.4 --refine-lower 0.3,0.6 -o e387_200.json
yfoil geometry repanel e387.dat -n 160 --method cosine --cosine-te-bias 0.3 -o e387_cosine.json
yfoil geometry repanel sharp.dat -n 160 --te-gap 0.002 -o sharp_gap.json
yfoil geometry repanel e387.dat --panelling panelling.json   # every option from one file
```

### Methods

**`pangen` (default)** is XFOIL's `PANGEN`, translated line for line and gated
against XFOIL's own output ([validation](../validation/geometry/README.md)).
It splines the input, forms the curvature along the arc length, smooths it, adds
a *fictitious* curvature at the trailing edge and inside optional refinement
windows, and then places the nodes so that
\((1 + 6\,\mathrm{CVPAR}\cdot\kappa)\,\Delta s\) is the same on every
panel. Its parameters are XFOIL's `PPAR` menu, here by descriptive names:

| `PPAR` key | XFOIL variable | Flag | Default | Meaning |
| --- | --- | --- | --- | --- |
| `N` | `NPAN` | `-n`, `--panels` | 160 | Number of panel nodes |
| `P` | `CVPAR` | `--curvature-bunching` | 1.0 | Curvature attraction; 0 gives uniform arc-length spacing |
| `T` | `CTERAT` | `--te-curvature-ratio` | 0.15 | Fictitious trailing-edge curvature as a fraction of the leading-edge curvature ("TE/LE panel density ratio") |
| `R` | `CTRRAT` | `--refined-curvature-ratio` | 0.2 | Fictitious curvature inside the refinement windows as a fraction of the leading-edge curvature |
| `XT` | `XSREF1`, `XSREF2` | `--refine-upper X1,X2` | off | Upper-surface refinement window in \(x/c\) |
| `XB` | `XPREF1`, `XPREF2` | `--refine-lower X1,X2` | off | Lower-surface refinement window in \(x/c\) |

**`cosine`** is yFoil's own method and has no XFOIL equivalent: cosine spacing in
arc length on each surface, its parameter warped by a power law set by
`--cosine-te-bias` (1 is a plain cosine; below 1 coarser at the trailing edge
and finer at the leading edge; above 1 finer at the trailing edge; default 0.15,
clamped to 0.05…2). It writes N + 1 nodes, its historic behaviour, which is frozen
because test fixtures were derived with it.

### Trailing edge

`--sharp` moves the two trailing-edge nodes to their midpoint (a closed edge,
XFOIL's `SHARP` path); `--te-gap GAP [--te-blend F]` sets the gap with XFOIL's
`TGAP`, moving the surfaces apart by
\(\tfrac{1}{2}\Delta\,(x/c)\,e^{-(1 - x/c)(1/F - 1)}\) each. The two are
exclusive. Both are applied **after** the nodes are distributed, so the blend
profile is evaluated exactly at every output node rather than splined through a
coarse input; XFOIL applies `TGAP` to the buffer aerofoil and `PANE` follows,
and the two orders differ by that interpolation. On a closed edge `TGAP`
delivers \(\cos(\tau/2)\) times the requested gap
([known issues §6.5](../xfoil-known-issues.md)).

### One file instead of flags

`--panelling FILE` reads every option from a JSON file in exactly the shape a
geometry file records under `generator.panelling`, so a record can be copied
out of one geometry and applied to another. It cannot be combined with any
other panelling flag; keys of the other method and unknown keys are errors, and
PPAR keys left out take XFOIL's defaults.

``` json
{"method": "pangen", "n_nodes": 160, "sharp_te": false, "te_gap": {"gap": 0.002, "blend": 1.0},
 "curvature_bunching": 1.5, "te_curvature_ratio": 0.3, "refined_curvature_ratio": 0.5,
 "refine_upper": [0.2, 0.4], "refine_lower": [0.3, 0.6]}
```

### Conflicts

Flags of the other method (`--cosine-te-bias` with `pangen`, any PPAR flag with
`cosine`), `--sharp` with `--te-gap`, `--te-blend` without `--te-gap`, a window
that is not two increasing values, and `--panelling` with any other panelling
flag are all errors, from flags and from a file alike.

## info

``` sh
yfoil geometry info naca4412.json
yfoil geometry info naca4412.json -o info.json
```

Prints a summary of the section — with `-o`, writes it as JSON instead.
