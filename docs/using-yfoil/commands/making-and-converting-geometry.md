---
icon: lucide/spline
---

# Making and converting geometry

``` sh
yfoil geometry <SUBCOMMAND>   # alias: yfoil geom
```

Five subcommands: [`naca`](#naca), [`karman-trefftz`](#karman-trefftz),
[`convert`](#convert), [`repanel`](#repanel) and [`info`](#info). Every one of
them prints its full help when run without arguments; `--help` gives the long
form.

## Geometry conventions

- Coordinates are normalised by chord.
- Panel ordering is **trailing edge → upper surface → leading edge → lower
  surface → trailing edge**, counter-clockwise.
- Written `.dat` files carry 17 significant figures, so a JSON → `.dat` → JSON
  round trip is bit-exact.
- A generated or repanelled geometry carries a `generator` record: the series,
  designation and parameters it was made from, and how it was panelled (see
  [provenance of a generated section](../aerofoil-families/provenance-of-a-generated-section.md)).

## The `naca` command { #naca }

Generate a NACA section from its designation. Every family with a public,
reproducible definition is covered; the formulas, what each family is for and
the references are in [aerofoil families](../aerofoil-families/index.md).

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
| `-n`/`--nodes`, `--method`, PPAR flags, `--sharp`, `--te-gap`, `--te-blend`, `--panelling` | as [`repanel`](#repanel) | The panelling: `pangen` (default) samples the section at `n_buffer_nodes` stations (246, XFOIL's own buffer density; settable in a `--panelling` file) and runs XFOIL's PANGEN on them; `cosine` is the analytic cosine sampling in \(x\) (no bias applies) |
| `-o`, `--output` | `naca<designation>.<format>` | Output file |

`--thickness vertical` reproduces XFOIL's own `NACA` command, which lays the thickness
vertically rather than perpendicular to the camber line; why yFoil's default differs is
[thickness perpendicular to the camber line](../aerofoil-families/thickness-perpendicular-to-the-camber-line.md).
Each section is generated as defined, sharp trailing edges included; `--sharp` and `--te-gap`
change the trailing edge after panelling, as under [`repanel`](#repanel)
([trailing-edge gap and sharpening](../aerofoil-families/trailing-edge-gap-and-sharpening.md)).

## The `karman-trefftz` command { #karman-trefftz }

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
| `-n`/`--nodes`, `--method`, PPAR flags, `--sharp`, `--te-gap`, `--te-blend`, `--panelling` | as [`repanel`](#repanel) | The panelling, `pangen` by default |
| `-o`, `--output` | `karman-trefftz.<format>` | Output file |

## The `convert` command { #convert }

``` sh
yfoil geometry convert naca0012.json --to dat -o naca0012.dat
yfoil geometry convert e387.dat --to json -o e387.json
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--to` | `json` | Output format: `json` or `dat` |
| `--name` | `Aerofoil` | Name written into the `.dat` header |
| `-o`, `--output` | stdout | Output file |

## The `repanel` command { #repanel }

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

### Panelling flags

| Flag | Default | Meaning |
| --- | --- | --- |
| `-n`, `--nodes` | 160 | Number of panel nodes, trailing edge round to trailing edge (XFOIL's `PPAR N`) |
| `--method` | `pangen` | `pangen` (XFOIL's `PANGEN`) or `cosine` (yFoil's own) |
| `--curvature-bunching`, `--te-curvature-ratio`, `--refined-curvature-ratio`, `--refine-upper X1,X2`, `--refine-lower X1,X2` | XFOIL's | `pangen` only: XFOIL's `PPAR` parameters `P`, `T`, `R`, `XT`, `XB` |
| `--cosine-te-bias` | 0.15 | `cosine` only: the trailing-edge bias of the cosine spacing |
| `--sharp` | off | Close the trailing edge after panelling |
| `--te-gap GAP`, `--te-blend F` | off | Set the trailing-edge gap after panelling, as XFOIL's `TGAP` does |
| `--panelling FILE` | — | Every option above from one file (below) |

What each method does and what its parameters mean is
[panel distributions](../aerofoil-families/panel-distributions.md); what `--sharp` and
`--te-gap` do to the nodes is
[trailing-edge gap and sharpening](../aerofoil-families/trailing-edge-gap-and-sharpening.md).

### One file instead of flags

`--panelling FILE` reads every option from a JSON file in exactly the shape a
geometry file records under `generator.panelling`, so a configuration can be
kept with a study, applied to any section or loaded geometry, and lifted back
out of a geometry that was made with it. It cannot be combined with any other
panelling flag; keys of the other method and unknown keys are errors, and keys
left out take the defaults.

A complete configuration, `fine-te.json` — 200 nodes, stronger curvature
bunching, a denser trailing edge, refinement windows over the region where
transition is expected on both surfaces, a 0.2% chord trailing-edge gap blended
over 80% of the chord, and the analytic section sampled at 400 stations before
PANGEN:

``` json
{
  "method": "pangen",
  "n_nodes": 200,
  "n_buffer_nodes": 400,
  "sharp_te": false,
  "te_gap": {"gap": 0.002, "blend": 0.8},
  "curvature_bunching": 1.3,
  "te_curvature_ratio": 0.3,
  "refined_curvature_ratio": 0.5,
  "refine_upper": [0.25, 0.45],
  "refine_lower": [0.35, 0.6]
}
```

Generating from it needs no other flag, and the same file repanels a loaded
geometry (drop `n_buffer_nodes` for that: an existing geometry has no section
to sample):

``` sh
yfoil geometry naca 63-415 --panelling fine-te.json -o naca63-415_fine-te.json
yfoil geometry karman-trefftz --te-angle 12 --panelling fine-te.json -o kt_fine-te.json
yfoil geometry repanel e387.dat --panelling fine-te-repanel.json -o e387_fine-te.json
```

The output's `generator.panelling` is the file's contents again (with
`te_bias: null` and `n_buffer_nodes` filled in as applicable), which is what
makes a run reproducible from its result alone.

### Conflicts

Flags of the other method (`--cosine-te-bias` with `pangen`, any PPAR flag with
`cosine`), `--sharp` with `--te-gap`, `--te-blend` without `--te-gap`, a window
that is not two increasing values, and `--panelling` with any other panelling
flag are all errors, from flags and from a file alike.

## The `info` command { #info }

``` sh
yfoil geometry info naca4412.json
yfoil geometry info naca4412.json -o info.json
```

Prints a summary of the section — with `-o`, writes it as JSON instead.
