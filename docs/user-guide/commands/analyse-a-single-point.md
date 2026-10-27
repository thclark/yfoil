---
icon: lucide/activity
---

# Analyse a single point

``` sh
yfoil analyse <GEOMETRY.json> [OPTIONS]
```

## Fixed angle of attack

``` sh
yfoil analyse naca0012.json --alpha 5 --reynolds 1e6 --ncrit 9 -o a5.json
```

## Fixed lift coefficient

``` sh
yfoil analyse naca0012.json --cl 0.6 --reynolds 1e6 -o cl06.json
```

α becomes the unknown; `--cl` overrides `--alpha`.

## Inviscid

``` sh
yfoil analyse naca0012.json --alpha 5 --inviscid
```

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `-a`, `--alpha` | 0 | Angle of attack, degrees |
| `--cl` | — | Fixed-\(C_L\) mode; overrides `--alpha` |
| `-r`, `--reynolds` | 1e6 | Reynolds number |
| `-m`, `--mach` | 0 | Mach number (Kármán–Tsien compressibility) |
| `--ncrit` | 9 | Critical amplification factor for transition |
| `--inviscid` | off | Inviscid solution only |
| `--max-iterations` | 20 | Maximum viscous Newton iterations (XFOIL's `ITER`) |
| `-o`, `--output` | — | Write the full result record as JSON |
| `--include-lagged-closures` | off | Also emit XFOIL's lagged closure arrays |
| `--allow-invalid` | off | Report the numbers of invalid points too; their status is unchanged |
| `--no-twins` | off | Do not run the twins: skip the conditioning check and report every converged, in-domain point as XFOIL would |
| `--spread-warning-threshold` | 1e-8 | Spread over the twins beyond which a valid point carries a warning that it may be ill-conditioned |
| `--spread-error-threshold` | 1e-4 | Spread over the twins beyond which a point is ill-conditioned: `Invalid`, numbers withheld |

Without `-o`, a summary is printed: \(C_L\), \(C_D\) (split into friction and
pressure), \(C_M\), transition locations and the iteration count. With `-o`, the
[full record](../reference/the-data-model.md) is written.

## Ill-conditioned points

By default every point is solved six times: once as asked, and five more on copies of the panels
jittered by one ULP (its [twins](../reference/solution-validity.md#twins)), in parallel. Where the
twins land far from the point, a rounding error would have given a different answer, so the point is
reported `Invalid (IllConditioned)` and its numbers are withheld. A NACA 0012 at −16°:

``` sh
yfoil geometry naca 0012 -n 160 --method cosine -o naca0012.json
yfoil analyse naca0012.json --alpha -16 --max-iterations 100
```

``` text
Status:  Invalid (IllConditioned); warning: solution converged but 2 twins did not

Results withheld: the point is ill-conditioned — a 1-ULP jitter of the panel coordinates moves cl by 7.2e-1 (error threshold 1e-4).
Pass --allow-invalid to report them anyway, or --no-twins to skip the check.

Diagnostics:
  Karman-Tsien margin (forces):   +1.0000e0
  Karman-Tsien margin (pressure): +1.0000e0
  Iterations: 78
  rms:        9.14e-5
  Twins:      largest spread 7.23e-1 (cl); warning threshold 1e-8, error threshold 1e-4
```

The point converged, so XFOIL would report it, and with `--no-twins` yFoil does too:
\(C_L = -0.383\). The twins are not changes to that answer — the point itself is the same solve with
or without them — but five more solves beside it. Three converge to within 2e-5 of it; the other
two run out of iterations, one at \(C_L = -0.163\) and one at \(C_L = -1.381\). It is a value that
decides, so the point is ill-conditioned; that two twins did not converge is the warning.

A point whose twins agree, but less closely than the warning threshold, keeps its status and
numbers and carries a warning instead, such as
`Valid; warning: solution may be ill-conditioned: its twins spread cl by 1.86e-7`. The two
thresholds are [explained with the twins](../reference/solution-validity.md#twins).

`--allow-invalid` shows the numbers and keeps the verdict; `--no-twins` skips the twins altogether,
roughly halving the run time, and reports such points as XFOIL does — numbers with no verdict on
their conditioning. The twins never change the numbers of the point itself.

## The XFOIL equivalent

``` text
LOAD naca0012.dat
OPER
VISC 1000000
ITER 20
ALFA 5
```
