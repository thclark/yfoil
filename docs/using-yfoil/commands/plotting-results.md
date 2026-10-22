---
icon: lucide/chart-line
---

# Plotting

!!! note "Feature-gated"

    Plotting is behind a cargo feature. Build with `--features plotting`.

`yfoil plot` writes [ironlab](https://github.com/thclark/ironlab) figures: files that hold
the figure's own data structure (every curve's numbers, its axes, key and labels), not a
picture of it. `-o` names the file; a `.fig` extension writes the compact binary form and
`.fig.json` the same figure as JSON. Any other extension is refused.

## Viewing and exporting

Open a figure in the ironlab viewer, a desktop application that pans, zooms, reads values off
the curves and hides series from the key:

``` sh
cargo install ironlab-viewer --version 0.12.0
ironlab-viewer naca0012.fig
ironlab-viewer *.fig           # several: the figure browser lists them
```

**Export PDF…** in the viewer's toolbar writes the figure as a vector PDF, as it is
currently shown (limits, hidden series and all), for a paper or a report. **Save figure…**
saves it back as a `.fig`.

Every figure carries **parameters**, which the viewer's figure browser sorts, filters and
searches by: `study = yfoil-plot`, `figure` (`foil`, `analysis` or `polar`), the flow
conditions (`re`, `mach`, `ncrit`, `x_trip_upper`, …, named as in the `conditions` record
of the JSON output) and, for a single operating point, its results (`alpha_deg`, `cl`,
`cd`, `cm`, `x_transition_upper`, `x_separation_upper`, `dstar_te_upper`, …). A figure of
several operating points or polars carries only the conditions they share. Figures are also
**labelled** by regime: `inviscid`, `compressible`, `forced-transition`, `low-re`. So a
directory of figures from many runs can be opened together and narrowed with, for example,
`re>=1000000 label:compressible`.

## The `plot polar` command { #polar }

``` sh
yfoil plot polar naca0012_polar.json                    # writes naca0012_polar.fig
yfoil plot polar naca0012_polar.json naca4412_polar.json \
    --title "NACA 0012 vs 4412" -o compare.fig
```

The comparison, from the polars of [your first analysis](../get-started/your-first-analysis.md) and [Polars](polars.md):

<ironlab-figure src="../figures/compare.fig" format="fig" name="compare" alt="NACA 0012 and NACA 4412 polars at Re 1e6: CL against alpha, CL against CD, CM against alpha and CD against alpha">NACA 0012 and NACA 4412 polars at Re 1e6: CL against alpha, CL against CD, CM against alpha and CD against alpha</ironlab-figure>

Draws \(C_L\)–α, \(C_L\)–\(C_D\), \(C_M\)–α and \(C_D\)–α, with one key below the
four panels. Several files are overlaid for comparison, each labelled by its `--label`
(or its aerofoil name) and drawn in its own colour and marker. α is linked across the
three α panels and \(C_L\) between the two \(C_L\) panels, so zooming one in the viewer
moves the others with it. Only valid points are drawn. Several files default to
`polar_comparison.fig`.

## The `plot analysis` command { #analysis }

``` sh
yfoil plot analysis a5.json                              # writes a5.fig
```

For the NACA 4412 at α = 5°, Re = 10⁶:

<ironlab-figure src="../figures/a5.fig" format="fig" name="a5" alt="NACA 4412 at alpha 5 degrees, Re 1e6: -Cp and Ue over x/c, upper and lower surfaces, above the aerofoil outline">NACA 4412 at alpha 5 degrees, Re 1e6: -Cp and Ue over x/c, upper and lower surfaces, above the aerofoil outline</ironlab-figure>

\(-C_p\) and \(U_e/U_\infty\) over \(x/c\) from a single-point analysis, upper surface
blue and lower red, above a strip with the aerofoil outline; the three share a linked
\(x\) axis. \(C_p\) is plotted negated, with suction up, because ironlab has no reversed
axis. A point whose numbers were withheld (an invalid point analysed without
`--allow-invalid`) has no distributions to draw, and is refused.

## The `plot foil` command { #foil }

``` sh
yfoil plot foil naca4412.json --panels notches          # writes naca4412_foil.fig

yfoil plot foil naca4412_polar.json --alpha 0,5,10 \
    --quantity dstar,theta --wake -o foil.fig
```

The two commands above, on the NACA 4412 of [Polars](polars.md):

<ironlab-figure src="../figures/naca4412_foil.fig" format="fig" name="naca4412_foil" alt="NACA 4412 panels, drawn as notches">NACA 4412 panels, drawn as notches</ironlab-figure>

<ironlab-figure src="../figures/foil.fig" format="fig" name="foil" alt="NACA 4412 at alpha 0, 5 and 10 degrees: displacement and momentum thickness normal to the surface, with the wake, and the stagnation, transition and separation markers">NACA 4412 at alpha 0, 5 and 10 degrees: displacement and momentum thickness normal to the surface, with the wake, and the stagnation, transition and separation markers</ironlab-figure>

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
| `-q`, `--quantity` | — | `dstar theta delta h hk hstar ue cf cdiss sqrtctau sqrtctaueq us cp mass_defect` |
| `--scale` | auto | Fixed offset multiplier (1 draws δ\*, θ, δ at true size) |
| `--max-offset` | 0.1 | Auto-scaling: largest value sits this fraction of chord off the surface |
| `--markers` / `--no-markers` | all | `stagnation`, `transition`, `separation` |
| `--label` | file stem | Design-point labels, in input order |
| `--title` | the design-point label | Figure title |
| `-o` | `<stem>_foil.fig` | Output figure, `.fig` or `.fig.json` |

Each quantity gets its own hue and its own scale; each design point gets a tone
and dash pattern, so several angles overlay legibly. The markers are open shapes, one per
kind (○ stagnation, □ transition, ◇ forced transition, △ separation, ▽ reattachment),
yellow on the upper surface and green on the lower, and the key below the axes names each
kind and side drawn. The axes are scaled equally in \(x\) and \(y\) at the figure's own
size, so the section is undistorted; zooming in the viewer changes the two independently.
