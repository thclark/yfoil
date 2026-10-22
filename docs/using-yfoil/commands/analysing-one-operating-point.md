---
icon: lucide/activity
---

# Analysing one operating point

One operating point, solved from scratch: a fixed angle of attack or a fixed lift coefficient,
viscous or inviscid.

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

## An inviscid point

``` sh
yfoil analyse naca0012.json --alpha 5 --inviscid
```

## Options of `analyse`

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

Without `-o`, a summary is printed: \(C_L\), \(C_D\) (split into friction and
pressure), \(C_M\), transition locations and the iteration count. With `-o`, the
[full record](../output-data.md#analysis-records) is written.

## The XFOIL equivalent of `analyse`

``` text
LOAD naca0012.dat
OPER
VISC 1000000
ITER 20
ALFA 5
```

A sequence of points, each starting from the last, is a polar: [sweeping a polar](sweeping-a-polar.md).
