---
icon: lucide/terminal
---

# Commands

yFoil's command-line tool is four commands, and they compose through files:

``` mermaid
graph LR
  A[yfoil geometry] -->|geometry.json| B[yfoil analyse]
  A -->|geometry.json| C[yfoil polar]
  B -->|analysis.json| D[yfoil plot]
  C -->|polar.json| D
```

| Command | Page | What it does |
| --- | --- | --- |
| `yfoil geometry` | [Making and converting geometry](making-and-converting-geometry.md) | Generate, convert, repanel and inspect aerofoil geometry |
| `yfoil analyse` | [Analysing one operating point](analysing-one-operating-point.md) | Solve a single operating point (fixed α or fixed \(C_L\)) |
| `yfoil polar` | [Sweeping a polar](sweeping-a-polar.md) | Sweep α as a state machine, exactly as XFOIL's `ASEQ` does |
| `yfoil plot` | [Plotting results](plotting-results.md) | Draw polars, \(C_p\)/\(U_e\) distributions and the foil itself as ironlab figures |

Every command that reads a foil takes **JSON geometry**. Use
[`yfoil geometry convert`](making-and-converting-geometry.md#convert) to bring in a `.dat` file.
Every command prints its full help when run without arguments; `--help` gives the long form.

The same analyses can be run from a program instead: [calling yFoil from code](calling-yfoil-from-code.md).

What the commands write is documented field by field in [output data](../output-data.md), and
why a point can come back with a status and no numbers in
[whether to trust a result](../whether-to-trust-a-result.md).
