---
icon: lucide/book-open
---

# User guide

## Start here

[Install yFoil](get-started/install.md), then [run your first analysis](get-started/your-first-analysis.md).
Ten minutes, and you will have a converged viscous solution and a plot of it, like
this comparison of two polars drawn with [`yfoil plot`](commands/plotting.md) as an
[ironlab](https://github.com/thclark/ironlab) figure (pan, zoom, or click the key to hide a curve):

<ironlab-figure src="figures/intro.fig" format="fig" name="intro" alt="NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha">NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha</ironlab-figure>

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
| [`yfoil geometry`](commands/geometry.md) | Generate, convert, repanel and inspect aerofoil geometry |
| [`yfoil analyse`](commands/analyse-a-single-point.md) | Solve a single operating point (fixed α or fixed \(C_L\)) |
| [`yfoil polar`](commands/polars.md) | Sweep α as a state machine, exactly as XFOIL's `ASEQ` does |
| [`yfoil plot`](commands/plotting.md) | Draw polars, \(C_p\)/\(U_e\) distributions and the foil itself as ironlab figures |

Every command that reads a foil takes **JSON geometry**. Use
[`yfoil geometry convert`](commands/geometry.md#convert) to bring in a `.dat` file.

[Bindings](commands/bindings.md) covers calling yFoil from other languages, and
the [data model](reference/the-data-model.md) documents the JSON that comes out.
[Solution validity](reference/solution-validity.md) explains why a point may come back with a status and no numbers,
and what each reason code means.

## Scope

yFoil reproduces XFOIL's OPER *analysis* functionality: fixed-alpha and fixed-CL
points, forced transition, polars and Kármán–Tsien compressibility. Inverse
design, interactive mode, multi-element sections and flap hinge moments are out
of scope.

Everything about how closely it reproduces XFOIL, and how that is measured, is
in [validation](../validation/index.md).

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

Please cite both yFoil and XFOIL: see [Citing yFoil](get-started/citing-yfoil.md).
