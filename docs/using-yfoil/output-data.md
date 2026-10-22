---
icon: lucide/braces
---

# Output data

Everything yFoil reads and writes is JSON, and the keys **are** the variable
names — no abbreviations, no serialisation aliases. There is one analysis shape,
whether it comes from `yfoil analyse -o` or is embedded in a polar written with
`--distributions`. This page lists the fields of each file; what a point's `status` means, and
when its numbers are withheld, is [whether to trust a result](whether-to-trust-a-result.md).

## Geometry files

``` json
{
  "x": [1.0, 0.999, ...],
  "y": [0.0013, 0.0014, ...],
  "cm_ref": [0.25, 0.0]
}
```

`x` and `y` are the panel nodes, trailing edge → upper → leading edge → lower →
trailing edge. `cm_ref` is the moment reference point. A geometry yFoil generated or repanelled
also carries `generator`, the record of the section and its panelling, described field by field
in [provenance of a generated section](aerofoil-families/provenance-of-a-generated-section.md).

## Analysis records

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

## Polar records

``` text
foil          name of the section
conditions    as above
results       one entry per angle, each the same shape as an analysis `results`
summary       cl_max, alpha_at_cl_max, ldratio_max, cl_at_ldratio_max, cd0,
              n_converged, n_failed
completed     whether the sweep ran to both limits
```

With `--distributions`, every entry also carries the full analysis record for
that angle. Every angle that was asked for has an entry, including the ones a sweep that halted
never reached; results are ordered ascending in alpha.

## Point identifiers and the initialisation chain

Every point record carries `id` and `initialised_from`: the `id` of the point whose converged
boundary layer seeded it, absent for a fresh march. A polar is a state machine whose results are
presented ascending in alpha, so this chain is the only record of the order the points were
solved in: two legs, both rooted at the 0° solve ([sweeping a polar](commands/sweeping-a-polar.md#how-the-sweep-runs)).

`id` is **content-addressed**: a 12-character base-36 hash (FNV-1a, 128-bit) of everything the
solve depended on — the yFoil version, the panel geometry, the flow conditions, the operating
point, the iteration limit, and the predecessor's id. So:

- it is deterministic, and a result file is byte-reproducible;
- equal ids mean the same computation, so results from different runs can be merged, or a sweep
  restarted from a point of an earlier one, without colliding references;
- the chain is a Merkle chain: the same operating point reached along different paths has
  different ids.

The re-solved 0° that seeds a polar's downward leg has the same inputs as the first 0° solve and
no predecessor, so it has the same id, and the downward leg cites the recorded point rather than a
phantom one. That is sound because the two solves are identical, which a test asserts
([polars and reinitialisation](../fidelity-to-xfoil/execution-equivalence/polars-and-reinitialisation.md)).

## Why JSON, and why these names

XFOIL's output formats are fixed-width Fortran (`G15.7`, `F9.4`, `F11.5`). They
cannot represent a result to better than about 1e-5, which makes them useless as
a validation reference and awkward as an automation interface. JSON with full
round-tripping precision is both.

The names follow [the naming conventions](../contributing/conventions/naming.md): symbols
from the equations, not XFOIL's six-character Fortran abbreviations. The XFOIL
name of every variable is recorded in the [translation
map](../fidelity-to-xfoil/translation-map/variable-and-routine-names.md).
