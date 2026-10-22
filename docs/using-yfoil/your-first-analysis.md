---
icon: lucide/play
---

# Your first analysis

One session from start to finish: make a section, solve one point, sweep a polar and plot it.
Each step links to the page that covers its command in full; you need none of them to finish.

## 1. Make a geometry

Every analysis command takes a JSON geometry file. Generate one:

``` sh
yfoil geometry naca 0012 -n 160 -o naca0012.json
```

That is a NACA 0012 panelled with 160 panels, ordered from the trailing edge over
the upper surface, round the leading edge, and back along the lower surface.

## 2. Solve one operating point

``` sh
yfoil analyse naca0012.json --alpha 5 --reynolds 1e9 -o naca0012_a5.json
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
yfoil plot polar naca0012_polar.json -o naca0012_polar.fig
yfoil plot analysis naca0012_a5.json -o naca0012_a5.fig
ironlab-viewer naca0012_polar.fig naca0012_a5.fig
```

The figures open in the [ironlab viewer](commands/plotting-results.md#viewing-and-exporting), which exports
each as a vector PDF. Using [`yfoil plot`](commands/plotting-results.md) and [ironlab](https://github.com/thclark/ironlab)
to visualise them, this analysis should produce the polar and distributions below.

<ironlab-figure src="your-first-analysis/naca0012_polar.fig" format="fig" name="naca0012_polar" alt="NACA 0012 polar at Re 1e6, alpha -5 to 15 degrees">NACA 0012 polar at Re 1e6, alpha -5 to 15 degrees</ironlab-figure>

<ironlab-figure src="your-first-analysis/naca0012_a5.fig" format="fig" name="naca0012_a5" alt="NACA 0012 at alpha 5 degrees, Re 1e9: -Cp and Ue over x/c, upper and lower surfaces, above the aerofoil outline">NACA 0012 at alpha 5 degrees, Re 1e9: -Cp and Ue over x/c, upper and lower surfaces, above the aerofoil outline</ironlab-figure>

## Next

- [Commands](commands/index.md) covers every option of the four commands used above.
- [Aerofoil families](aerofoil-families/index.md) says exactly which section `naca 0012` made, and how
  its 160 nodes were placed.
- [Output data](output-data.md) documents what is in the JSON files.
- [Whether to trust a result](whether-to-trust-a-result.md) explains what a point's `status` says
  before you rely on its numbers.
