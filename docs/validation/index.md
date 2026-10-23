---
icon: lucide/badge-check
---

# Validation

The purpose of this section is to give interested users a detailed insight into the studies undertaken
when translating yFoil from XFOIL, to ensure correctness. 

A direct 1:1 translation is not possible, due to variations between platforms (regardless of translation,
XFOIL itself yields different results on different platforms due to tiny variations in system math libraries),
making extensive validation and testing necessary to build confidence that yFoil matches the output of XFOIL, not 
just in **most**, but in **all**, possible cases.

The process is covered in depth in the written paper (see ['Citing yFoil'](../user-guide/get-started/citing-yfoil.md)).

!!! note "The meaning of 'Validation'"

    yFoil is a translation of XFOIL 6.99, so 'validation', here, means showing that it
    - gives XFOIL's answers where possible, and 
    - fails in a controlled way otherwise.
    
    This is **not** a validation that yFoil agrees with wind-tunnel data: there are plenty of papers
    detailing XFOIL's results against experiments which commutatively apply to yFoil too.


## Pristine and reference builds

To enable deep comparison between internal values XFOIL source code is copied, compiled in double precision, linked to 
an explicit system maths library (to deliver identical results across platforms) and instrumented to:
 - accept input values (like panel nodes) at double precision (17 significant figures),
 - dump intermediate state (like boundary layer variables or solver residuals) at double precision,
 - record events (like entry or exit of loops) for branch coverage studies, and
 - dump output state at double precision.

This is the **reference** build, against which we compare yFoil. However, to create the reference, the XFOIL source code was modified. 
To ensure those modifications did not inadvertently alter XFOIL's results, a complete set of tests was undertaken comparing the reference
XFOIL to an unaltered copy of the source code compiled with the same options (the **pristine** build).

A comparison between yFoil and the reference:
 - passes when both codes take the same branches, and every value is within a tolerance taken from XFOIL's own noise floor (see [Tolerances](tolerances/index.md))
 - is invalid where XFOIL is ill-conditioned (severely post-stall, separation, the Kármán–Tsien pole) 
 - is indeterminate where branches diverge only due to extremely small variations between yFoil and XFOIL (theoretically possible but not seen in the test set)
 - fails otherwise.


## The studies

Validation processes evolved during the build of yFoil, but the main steps have been captured and turned into a compact set of
reproducible studies:

- [Aerofoil series](aerofoil-series/index.md): Selects a set of different aerofoil geometry families ('series') for use in the
  other validation studies, each with a physical rationale for why they've been selected. Correct generation of aerofoil geometry 
  (panelling) is then validated against externally defined references like the NACA456 geometry program.

- [Series cases](series-cases/index.md): For each series, runs simulations a typical aerodynamicist would run (full polars swept 
  through stall at a range of operating points appropriate to the series), then compares those polars against the reference.
 
- [Known issues](known-issues/index.md): Catalogues a wide number of issues with XFOIL's behaviour (from minor 
  issues like dead code or outdated comments to pathological behaviour like program crashes) that emerged during development of yFoil.
  
- [Branch coverage](branch-coverage/index.md): Extracts the smallest set of analysis cases that takes **every reachable branch** of XFOIL's 
  analysis path. This study was used as a basis to generate unit test fixtures that comprehensively demonstrate yFoil's consistency with XFOIL
  in every situation.

- [Tolerances](tolerances/index.md): Determines much XFOIL's internal and output values move under 1-ULP ('Unit of Least Precision') changes to its
  inputs. This study was used as a basis to generate tolerances for comparisons between yFoil and XFOIL (for the series cases and for unit test fixtures).


## Reproducing studies

Each study page states the run it shows. To make a run yourself, generate the study's fixtures (where it has them), then run the
study. Every command runs from the repository root, and each run writes a dated folder under `studies/<study>/runs/<UTC>/`.

```sh
# Aerofoil series (no fixtures needed)
cargo run --release -p aerofoil-series

# Series cases
cargo xtask fixtures --group series
cargo run --release -p series-cases

# Known issues
cargo xtask fixtures --group known-issues
cargo run --release -p known-issues

# Branch coverage
cargo xtask fixtures --group branch-coverage
cargo xtask fixtures --group non-finite
cargo run --release -p branch-coverage

# Tolerances (builds the reference if it is missing)
studies/tolerances/tolerances.sh
```

A run never changes documentation pages. How the pages and the runs fit together is set out in 
[Documentation](../development/conventions/documentation.md). To update a page's tables and figures from 
a run, publish it:

```sh
cargo xtask publish <study> --dry-run       # Shows what would change
cargo xtask publish <study>                 # Copies results across to the docs from the newest run...
cargo xtask publish <study> --run <UTC>     # ...or a chosen historical run
```

`<study>` is `aerofoil-series`, `series-cases`, `known-issues`, `branch-coverage` or `tolerances`.


!!! note "Developer note"
    A run folder made before a change to a study's output can be redrawn without rerunning the study, using
    `cargo run --release -p <study> -- --render studies/<study>/runs/<UTC>` (every study but tolerances).


