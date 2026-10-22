# Provenance of a generated section

Every geometry file yFoil generates carries `generator`, the record of how it was made: the
series, the designation, the thickness form and mean line with their parameter values, how the
thickness was applied, whether the trailing edge is sharp, any trailing-edge adjustment, the
yFoil version and the keys of the [references](../bibliography.md) that define the family — and,
under `panelling`, how the nodes were distributed.

``` json
"generator": {
  "yfoil": "0.1.0",
  "series": "naca_6",
  "designation": "NACA 63-415",
  "thickness_form": { "family": "63", "t": 0.15 },
  "mean_line": { "family": "6", "cl": 0.4, "a": 1.0 },
  "thickness_applied": "perpendicular",
  "sharp_te": true,
  "references": ["abbott1945", "abbott1959", "ladson1974", "ladson1996", "carmichael2001"],
  "panelling": {
    "method": "pangen", "n_nodes": 160, "n_buffer_nodes": 246,
    "sharp_te": false, "te_gap": null,
    "curvature_bunching": 1.0, "te_curvature_ratio": 0.15, "refined_curvature_ratio": 0.2,
    "refine_upper": null, "refine_lower": null
  }
}
```

## The section

`series` is one of `naca_4_digit`, `naca_4_digit_modified`, `naca_5_digit`,
`naca_16`, `naca_6`, `naca_6a` or `karman_trefftz` (which carries `x_centre`,
`y_centre`, `te_angle_deg` and `exponent` instead of the NACA entries); the
thickness-form families are `4`, `4M` (with `le_radius_index`,
`x_max_thickness`) and `63`…`67`, `63A`…`65A`; the mean-line families `2`
(`m`, `p`), `3` and `3R` (`cl`, `p`), `6` (`cl`, `a`) and `6A` (`cl`), or
`null` for a symmetric section. `thickness_applied` is `perpendicular` or `vertical`
([thickness perpendicular to the camber line](thickness-perpendicular-to-the-camber-line.md)).
`sharp_te` at the top level is the section *definition's* property (the 6-series closes).

## The panelling

`panelling` is the complete recipe from the section (or the loaded nodes) to
the output nodes: `method` is `pangen` (XFOIL's PANGEN, with its `PPAR`
parameters as the following keys) or `cosine` (yFoil's own, with `te_bias`,
`null` when a generator's analytic sampling was used); `n_nodes` (panel nodes, trailing edge round to trailing edge); `sharp_te` and
`te_gap` (`{"gap", "blend"}` or `null`) are the trailing-edge treatment applied
after the distribution, so `panelling.sharp_te` says what was done to the nodes; `n_buffer_nodes`
appears when a generator sampled its section before PANGEN (default 246; a `--panelling` file may
set it). What each key means is [panel distributions](panel-distributions.md) and
[trailing-edge gap and sharpening](trailing-edge-gap-and-sharpening.md).

The same JSON is what `--panelling FILE` accepts, so a record can be lifted from one geometry and
applied to another ([one file instead of flags](../commands/making-and-converting-geometry.md#one-file-instead-of-flags)).

## What survives a repanel or a conversion

`yfoil geometry repanel` keeps the section record and replaces `panelling`; `--te-gap` changes the
nodes, not the section, so the section record survives it too. A repanelled `.dat` carries only
`yfoil` and `panelling`, since a `.dat` says nothing about its section, and a file converted from
`.dat` has no section record at all.

The rest of a geometry file is in [output data](../output-data.md#geometry-files).
