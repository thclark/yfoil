---
icon: lucide/play
---

# Your first analysis

## 1. Make a geometry

Every analysis command takes a JSON geometry file. Generate one:

``` sh
yfoil geometry naca 0012 -n 160 -o naca0012.json
```

That is a NACA 0012 panelled with 160 panels, ordered from the trailing edge over
the upper surface, round the leading edge, and back along the lower surface.

## 2. Solve one operating point

``` sh
yfoil analyse naca0012.json --alpha 5 --reynolds 1e6 -o naca0012_a5.json
```

The output JSON carries the operating conditions, the integrated forces, the
panel geometry and wake, the surface \(q\) and \(C_p\) distributions, and every
boundary-layer quantity at every station.

## 3. Sweep a polar

``` sh
yfoil polar naca0012.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --label "NACA 0012" -o naca0012_polar.json
```

A polar is a state machine, not a set of independent solves: each angle is
initialised from the previous one, exactly as XFOIL does it.

## 4. Look at it

``` sh
yfoil plot polar naca0012_polar.json -o naca0012_polar.svg
yfoil plot analysis naca0012_a5.json -o naca0012_a5.svg
```

## Next

The [user guide](index.md) covers each command in full, and the [data
model](data.md) documents what is in the JSON.
