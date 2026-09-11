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

Without `-o`, a summary is printed: \(C_L\), \(C_D\) (split into friction and
pressure), \(C_M\), transition locations and the iteration count. With `-o`, the
[full record](data.md) is written.

## The XFOIL equivalent

``` text
LOAD naca0012.dat
OPER
VISC 1000000
ITER 20
ALFA 5
```
