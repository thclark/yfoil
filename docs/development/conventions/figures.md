---
icon: lucide/image
---

# Figures

Every yFoil figure is an [ironlab](https://github.com/thclark/ironlab) figure: a file holding the figure's own
data structure (`.fig`, Protocol Buffers; `.fig.json`, the same model as JSON), not a picture of it. It is opened
in the ironlab viewer, embedded live in this site, and exported to PDF for publication. The shared conventions
are implemented once, in `src/output/figures.rs` (feature `plotting`), and every figure is built through it.

The exception is the branch-coverage study's two occupancy figures (`docs/validation/branch-coverage/`): the branch map
and the case × subroutine matrix are not charts, and stay SVG drawn by matplotlib (`studies/branch-coverage/plot.py`
through `studies/figures/render.sh`) until ironlab has a figure type for them.

## Who computes what

The Rust that builds a figure is the only place its numbers come from: data, derived quantities, selections and
axis extents. There is no presentation layer that transforms data.

## Size, fonts and colours

- **Size and font** are ironlab's defaults (set up for AMS journals). A figure that needs another aspect ratio, such
  as a tall multi-row sheet, keeps the default width (`figures::WIDTH_MM`) and sets only its height. Fonts are never
  set.
- **The reference and yFoil**: XFOIL is `figures::RED`, dashed, × markers; yFoil is `figures::BLUE`, solid, circles.
  Non-finite behaviour is `figures::ORANGE`.
- **Ordered families** (members of an aerofoil series, perturbation levels) get one explicit colour per line from
  `figures::ylgnbu(t)`, lightest = smallest, sampled exactly as matplotlib's YlGnBu was. A figure never relies on a
  colourmap for line colours. The family's default or base case is black and drawn on top.
- **Unordered families of up to four** use `figures::CATEGORY`.

## Parameters and labels

The viewer sorts, filters and groups figures by their **parameters** (name → value) and **labels** (plain words).
Every figure carries both, so that all the figures every study produces can be loaded into one viewer and sliced.

**Parameters** say *which* figure this is and *what went in and came out*:

| Group | Names | Notes |
|---|---|---|
| Provenance | `study`, `figure`, `case`, `series`, `foil` | `study` is the producing crate (`series-cases`, `known-issues`, `aerofoil-series`, `input-sensitivity`, `xfoil-instrumentation-check`, `yfoil-plot`); `figure` names the figure type (`polar-sheet`, `transition-strip`, `family-sheet`, `sensitivity-sheet`, `instrumentation-check`, `foil`, `analysis`, `polar`); `case` is the fixture case name; `series` is the aerofoil-series slug (`naca-4-digit`, …) |
| Inputs | `n_nodes`, `viscous`, `re`, `mach`, `ncrit`, `max_iterations`, `x_trip_upper`, `x_trip_lower`, `mach_cl_dependence`, `re_cl_dependence`, `amplification_model`, `alpha_deg`, `cl_target` | Named as in yFoil's `conditions` record (`figures::flow_condition_parameters`); `re` is absent for an inviscid figure |
| Outputs | `cl`, `cd`, `cm`, `converged`, `x_transition_upper`, `x_transition_lower`, `x_separation_upper`, `x_separation_lower`, `dstar_te_upper`, `dstar_te_lower`, `theta_te_upper`, `theta_te_lower` | For a figure of one operating point. A value that does not exist (no separation, an inviscid δ*) is left out (`figures::optional_parameter`), never written as a placeholder |
| Summary | `points`, `matched`, `parted`, `cl_max`, `alpha_cl_max_deg` | For a figure of a sweep: counts and extremes the study already computes |

A name that exists as a key in yFoil's own JSON output is spelled exactly as that key. A new value at a place
is named `<quantity>_<place>_<side>`.

**Labels** are for tags a reader filters on, not for values:

- the documentation section a figure appears in: `validation`;
- `xfoil-comparison` for a figure that plots the reference against yFoil;
- regime tags taken from the inputs: `inviscid`, `compressible` (M > 0), `forced-transition`, `fixed-cl`,
  `low-re` (Re < 1e5).

## Linked axes

Panels that share an independent variable (α, x/c) have that axis linked across the whole figure. Panels that plot
the same quantity, such as an upper/lower pair or the same quantity under two perturbation families, have their
dependent axis linked too. Two axes with different quantities are never linked, even if their ranges agree.

## Legends

Where a figure has one key for many panels, the legend goes in a **panel legend** (`fig.panel_legend(row, col)` +
`ax.legend_in(key)`) in a tile of its own. This replaces the hand-placed legend strips of the matplotlib figures.
A legend that belongs to one panel stays in that panel.

## Saving

- `figures::save(fig, path)` validates the figure and refuses to write one with validation errors. The encoding
  follows the extension.
- `figures::save_with_pdf(fig, stem)` also writes the PDF, for the run folders' LaTeX reports.
- Study run folders hold `.fig` and `.pdf`. The documentation holds `.fig` only.
- The user guide's figures (`docs/user-guide/figures/`) are the output of the very commands the guide prints, run by
  `scripts/guide-figures.sh`; a guide page never embeds a figure no command on it produces.
- `yfoil plot` writes what `-o` names (`.fig` or `.fig.json`), `.fig` by default.

## Embedding

A page embeds each of its figures with `figures::embed_html(path, alt)`, `path` relative to the page.
A validation study's figures sit in its section's `figures/` directory (`validation/series-cases/figures/`),
copied there from a study run by `cargo xtask publish`:

```html
<ironlab-figure src="figures/polars-naca0012_n160_polar30_re1e6_iter100.fig" name="polars-naca0012_n160_polar30_re1e6_iter100"
  alt="CL, CD, RMSBL, iterations and transition against alpha, XFOIL and yFoil">…alt text…</ironlab-figure>
```

The site loads the element from the pinned ironlab-web bundle, which `scripts/docs.sh` fetches and checksums
into `docs/javascripts/ironlab/` (gitignored). The alt text is also the fallback content. No PNG or PDF fallback
is committed.
