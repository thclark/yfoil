---
icon: lucide/square-terminal
---

# Overview

If you already know XFOIL, using yFoil's CLI will be quite familiar  — `yfoil analyse` is one OPER
point, `yfoil polar` is a `PACC`/`ASEQ` sweep... each command below gives the
equivalent XFOIL session.

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
| [`yfoil analyse`](analyse-a-single-point.md) | Solve a single operating point (fixed α or fixed \(C_L\)) |
| [`yfoil polar`](polars.md) | Sweep α as a state machine, exactly as XFOIL's `ASEQ` does |
| [`yfoil plot`](plotting.md) | Draw polars, \(C_p\)/\(U_e\) distributions and the foil itself as ironlab figures |

Every command that reads a foil takes **JSON geometry**. Use
[`yfoil geometry convert`](geometry.md#convert) to bring in a `.dat` file.

[Bindings](bindings.md) covers calling yFoil from other languages, and
the [data model](../reference/the-data-model.md) documents the JSON that comes out.
[Solution validity](../reference/solution-validity.md) explains why a point may come back with a status and no numbers,
and what each reason code means.
