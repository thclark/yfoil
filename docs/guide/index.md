---
icon: lucide/book-open
---

# User guide

## Start here

[Install yFoil](install.md), then [run your first analysis](first-analysis.md).
Ten minutes, and you will have a converged viscous solution and a plot of it.

If you already know XFOIL, the mapping is direct — `yfoil analyse` is one OPER
point, `yfoil polar` is a `PACC`/`ASEQ` sweep — and each command below gives the
equivalent XFOIL session.

## The commands

yFoil's CLI is four commands, and they compose through files:

``` mermaid
graph LR
  A[yfoil geometry] -->|geometry.json| B[yfoil analyse]
  A -->|geometry.json| C[yfoil polar]
  B -->|analysis.json| D[yfoil plot]
  C -->|polar.json| D
```

| Command | What it does |
| --- | --- |
| [`yfoil geometry`](geometry.md) | Generate, convert, repanel and inspect aerofoil geometry |
| [`yfoil analyse`](analyse.md) | Solve a single operating point (fixed α or fixed \(C_L\)) |
| [`yfoil polar`](polar.md) | Sweep α as a state machine, exactly as XFOIL's `ASEQ` does |
| [`yfoil plot`](plot.md) | Draw polars, \(C_p\)/\(U_e\) distributions and the foil itself |

Every command that reads a foil takes **JSON geometry**. Use
[`yfoil geometry convert`](geometry.md#convert) to bring in a `.dat` file.

[Bindings](../bindings.md) covers calling yFoil from other languages, and
the [data model](data.md) documents the JSON that comes out.

## Scope

yFoil reproduces XFOIL's OPER *analysis* functionality: fixed-alpha and fixed-CL
points, forced transition, polars and Kármán–Tsien compressibility. Inverse
design, interactive mode, multi-element sections and flap hinge moments are out
of scope.

Everything about how closely it reproduces XFOIL, and how that is measured, is
in [validation](../validation/README.md).

## Open source

!!! note "Stub — to be written"

    yFoil is GPL-3.0-or-later, because XFOIL is vendored in-tree and a
    translation is a derivative work. This section will cover what that means
    for you: what you can do with yFoil, what you owe back, and how to
    contribute.

## FAIR data and software

!!! note "Stub — to be written"

    yFoil's file formats are designed for Findable, Accessible, Interoperable
    and Reusable practice — one documented JSON schema in and out, full
    round-tripping precision, and no interactive state to reproduce. This
    section will set out how yFoil supports FAIR workflows.

## How to cite yFoil

TODO - add citations for our own papers and Mark Drela's.

Years work went into building and validating XFOIL. Authors, maintainers and
users of yFoil owe a greate debt of thanks to Prof. Drela, Harold Youngren and
any other contributors to the original XFOIL.
