"""Presentation rules of the branch-coverage study's occupancy figures (branch map, case matrix),
the only figures still drawn by matplotlib; every other figure is an ironlab figure built in Rust
(docs/conventions/figures.md).

The Rust study writes every number (data, derived quantities, extents, selections) into its run
folder; the plot script only presents them:

- figures are drawn at their final physical size and go into the paper at natural size
  (`\\includegraphics{...pdf}` with no `width=` and no `\\resizebox`);
- width is TEXT_WIDTH_PT, the `\\textwidth` of the target template (A4 with 1 in margins), or less;
- Times New Roman throughout; tick values and legends at 8 pt (the journal minimum), axis labels
  at 9 pt; math through mathtext in a Times-compatible face;
- SVG with text kept as text, and PDF with the fonts embedded, from the same figure object.

Run through studies/figures/render.sh, which pins matplotlib.
"""


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEXT_WIDTH_PT = 451.0
FONT = "Times New Roman"
TICK_PT = 8.0
LEGEND_PT = 8.0
AXIS_LABEL_PT = 9.0
PT = 1.0 / 72.0  # inches per point
GRID_COLOUR = "#e4e4e4"

def apply() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [FONT, "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": TICK_PT,
            "axes.labelsize": AXIS_LABEL_PT,
            "axes.titlesize": AXIS_LABEL_PT,
            "xtick.labelsize": TICK_PT,
            "ytick.labelsize": TICK_PT,
            "legend.fontsize": LEGEND_PT,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.75,
            "ytick.major.width": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID_COLOUR,
            "grid.linewidth": 0.5,
            "axes.axisbelow": True,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "figure.dpi": 100,
        }
    )


def figure(width_pt: float, height_pt: float):
    return plt.figure(figsize=(width_pt * PT, height_pt * PT))


def axes(fig, x0_pt: float, y0_pt: float, w_pt: float, h_pt: float):
    """An axes placed in points from the figure's bottom-left corner."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([x0_pt * PT / fw, y0_pt * PT / fh, w_pt * PT / fw, h_pt * PT / fh])
