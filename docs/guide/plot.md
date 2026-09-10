---
icon: lucide/chart-line
---

# Plotting

!!! note "Feature-gated"

    Plotting is behind a cargo feature. Build with `--features plotting`.

All plots are SVG by default; PNG is produced by rasterising the same SVG, so
there is only ever one drawing path.

## polar

``` sh
yfoil plot polar naca0012_polar.json -o naca0012.svg
yfoil plot polar naca0012_polar.json naca4412_polar.json \
    --title "NACA 0012 vs 4412" -o compare.svg
```

Draws \(C_L\)–α, \(C_L\)–\(C_D\), \(C_M\)–α and \(C_D\)–α. Several files are
overlaid for comparison, each labelled by its `--label`.

## analysis

``` sh
yfoil plot analysis a5.json -o a5.svg
```

\(C_p\) and \(U_e\) distributions from a single-point analysis.

## foil

``` sh
yfoil plot foil naca4412.json --panels notches -o panels.svg

yfoil plot foil naca4412_polar.json --alpha 0,5,10 \
    --quantity dstar,theta --wake -o foil.svg
```

Draws the section itself: panels, wake, boundary-layer quantities plotted normal
to the surface, and the stagnation, transition and separation markers. Input can
be geometry files (panels only), analysis JSONs (one design point each) or one
polar JSON written with `--distributions` (one design point per α).

| Option | Default | Meaning |
| --- | --- | --- |
| `--alpha` | all | Polar input only: which embedded angles to draw |
| `--panels` | `notches` for geometry, else `none` | `none`, `notches` or `dots` |
| `--notch-length` | 0.01 | Notch length as a fraction of chord |
| `--wake` | off | Draw the wake panels, and the wake band with `--quantity` |
| `-q`, `--quantity` | — | `dstar theta delta h hk hs ue cf cdis ctau ctq uslp cp mass` |
| `--scale` | auto | Fixed offset multiplier (1 draws δ\*, θ, δ at true size) |
| `--max-offset` | 0.1 | Auto-scaling: largest value sits this fraction of chord off the surface |
| `--markers` / `--no-markers` | all | `stagnation`, `transition`, `separation` |
| `--label` | file stem | Design-point labels, in input order |

Each quantity gets its own hue and its own scale; each design point gets a tone
and dash pattern, so several angles overlay legibly.
