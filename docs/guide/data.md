---
icon: lucide/braces
---

# The data model

Everything yFoil reads and writes is JSON, and the keys **are** the variable
names — no abbreviations, no serialisation aliases. There is one analysis shape,
whether it comes from `yfoil analyse -o` or is embedded in a polar written with
`--distributions`.

## Geometry

``` json
{
  "x": [1.0, 0.999, ...],
  "y": [0.0013, 0.0014, ...],
  "cm_ref": [0.25, 0.0]
}
```

`x` and `y` are the panel nodes, trailing edge → upper → leading edge → lower →
trailing edge. `cm_ref` is the moment reference point.

A geometry yFoil generated also carries `generator`, the record of how: the
series, the designation, the thickness form and mean line with their parameter
values, how the thickness was applied, whether the trailing edge is sharp, any
trailing-edge adjustment, the yFoil version and the keys of the
[references](../references.md) that define the family — and, under
`panelling`, how the nodes were distributed. `yfoil geometry repanel` keeps the
section record and replaces `panelling`; a repanelled `.dat` carries only
`yfoil` and `panelling`, since a `.dat` says nothing about its section.

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

`panelling` is the complete recipe from the section (or the loaded nodes) to
the output nodes: `method` is `pangen` (XFOIL's PANGEN, with its `PPAR`
parameters as the following keys) or `cosine` (yFoil's own, with `te_bias`,
`null` when a generator's analytic sampling was used); `n_nodes` (panel nodes, trailing edge round to trailing edge); `sharp_te` and
`te_gap` (`{"gap", "blend"}` or `null`) are the trailing-edge treatment applied
after the distribution; `n_buffer_nodes` appears when a generator sampled its
section before PANGEN (default 246; a `--panelling` file may set it). The same JSON is what `--panelling FILE` accepts, so a
record can be lifted from one geometry and applied to another (see
[repanel](geometry.md#repanel)). `sharp_te` at the top level is the section
*definition's* property (the 6-series closes); `panelling.sharp_te` says what was
done to the nodes.

`series` is one of `naca_4_digit`, `naca_4_digit_modified`, `naca_5_digit`,
`naca_16`, `naca_6`, `naca_6a` or `karman_trefftz` (which carries `x_centre`,
`y_centre`, `te_angle_deg` and `exponent` instead of the NACA entries); the
thickness-form families are `4`, `4M` (with `le_radius_index`,
`x_max_thickness`) and `63`…`67`, `63A`…`65A`; the mean-line families `2`
(`m`, `p`), `3` and `3R` (`cl`, `p`), `6` (`cl`, `a`) and `6A` (`cl`), or
`null` for a symmetric section. See [aerofoil series](aerofoil-series.md).

## Analysis record

``` text
foil                 name of the section
conditions           re, mach, ncrit, max_iterations, wake_length,
                     elimination_threshold, x_trip, mach_cl_dependence,
                     re_cl_dependence, amplification_model
results              alpha_deg, cl, cm, cd, cd_friction, cd_pressure, ldratio,
                     transition_upper, transition_lower, converged, iterations,
                     residual
geometry             x, y, s, normal_x, normal_y, chord, x_le, y_le, x_te, y_te,
                     s_le, i_le_node, te_thickness_normal, sharp_te, wake{…}
surface              q, cp                       (one value per panel node)
boundary_layer       upper{…}, lower{…}, wake{…}, plus i_te_station,
                     i_transition_station, n_wake_nodes, qinf, stagnation{…},
                     transition[…], wake_split
```

Each boundary-layer side carries its station and node indices, the station
coordinates (`x`, `y`, `xi`, `cp`), the **primaries** the solver actually marches

``` text
ue  theta  dstar  sqrtctau  mass_defect
```

and the **closures** evaluated on them

``` text
ue_compressible  h  hk  hstar  cf  cdiss  delta  sqrtctaueq  us  retheta
machsqd_edge
```

!!! note "Lagged closures"

    XFOIL also keeps lagged copies of some closure arrays. They are not written
    by default, because they are an artefact of the solver rather than a result;
    pass `--include-lagged-closures` to emit them under each side's
    `lagged_closures`.

## Polar record

``` text
foil          name of the section
conditions    as above
results       one entry per angle, each the same shape as an analysis `results`
summary       cl_max, alpha_at_cl_max, ldratio_max, cl_at_ldratio_max, cd0,
              n_converged, n_failed
completed     whether the sweep ran to both limits
```

With `--distributions`, every entry also carries the full analysis record for
that angle.

## Why JSON, and why these names

XFOIL's output formats are fixed-width Fortran (`G15.7`, `F9.4`, `F11.5`). They
cannot represent a result to better than about 1e-5, which makes them useless as
a validation gate and awkward as an automation interface. JSON with full
round-tripping precision is both.

The names follow [the naming conventions](../conventions/naming.md): symbols
from the equations, not XFOIL's six-character Fortran abbreviations. The XFOIL
name of every variable is recorded in the [mapping
table](../xfoil-reference.md#xfoil-yfoil-mapping).
