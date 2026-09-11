"""Publication figure conventions shared by the validation studies' plot scripts.

The Rust studies write every number (data, derived quantities, extents, selections) into their
run folder; the plot scripts only present them. This module carries the presentation rules so
that each script stays short and the figures look alike:

- figures are drawn at their final physical size and go into the paper at natural size
  (`\\includegraphics{...pdf}` with no `width=` and no `\\resizebox`);
- width is TEXT_WIDTH_PT, the `\\textwidth` of the target template (A4 with 1 in margins), or less;
- Times New Roman throughout; tick values and legends at 8 pt (the journal minimum), axis labels
  at 9 pt; math through mathtext in a Times-compatible face;
- SVG with text kept as text, and PDF with the fonts embedded, from the same figure object;
- YlGnBu (ColorBrewer) for ordered families, lightest = smallest, avoiding the near-white end.

Run through scripts/figures/render.sh, which pins matplotlib.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colormaps  # noqa: E402

TEXT_WIDTH_PT = 451.0
FONT = "Times New Roman"
TICK_PT = 8.0
LEGEND_PT = 8.0
AXIS_LABEL_PT = 9.0
PT = 1.0 / 72.0  # inches per point
GRID_COLOUR = "#e4e4e4"

# plotly's D3 category colours used by the instrumentation check
RED = "#d62728"
BLUE = "#1f5ad6"
CATEGORY = ["#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]


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


def ylgnbu(t: float):
    """Colour at t in [0, 1] on YlGnBu, avoiding the near-white end."""
    return colormaps["YlGnBu"](0.2 + 0.8 * max(0.0, min(1.0, t)))


def figure(width_pt: float, height_pt: float):
    return plt.figure(figsize=(width_pt * PT, height_pt * PT))


def axes(fig, x0_pt: float, y0_pt: float, w_pt: float, h_pt: float):
    """An axes placed in points from the figure's bottom-left corner."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([x0_pt * PT / fw, y0_pt * PT / fh, w_pt * PT / fw, h_pt * PT / fh])


def save(fig, svg: Path) -> None:
    """SVG and PDF side by side, then close."""
    fig.savefig(svg)
    fig.savefig(svg.with_suffix(".pdf"))
    plt.close(fig)
    print("wrote", svg)
    print("wrote", svg.with_suffix(".pdf"))


def legend_below(fig, x0_pt, top_pt, w_pt, entries, columns, row_pt=10.0, first_alone=False):
    """A legend strip whose top-left corner is at (x0_pt, top_pt): `entries` are (handle, text),
    laid out row by row in `columns`; with `first_alone` the first entry gets a row of its own.
    Returns the strip's height in points."""
    import math

    head, rest = (entries[:1], entries[1:]) if first_alone else ([], entries)
    rows = len(head) + math.ceil(len(rest) / columns)
    h = row_pt * rows + 4.0
    lax = axes(fig, x0_pt, top_pt - h, w_pt, h)
    lax.axis("off")
    # rows `row_pt` apart: the text is ~1.2 em tall, the rest is label spacing (in em)
    spacing = max(0.0, (row_pt - 1.2 * LEGEND_PT) / LEGEND_PT)
    common = dict(
        frameon=False,
        borderaxespad=0,
        borderpad=0,
        handlelength=1.5,
        columnspacing=1.0,
        handletextpad=0.5,
        labelspacing=spacing,
    )
    y = 1.0
    if head:
        top = lax.legend([head[0][0]], [head[0][1]], loc="upper left", bbox_to_anchor=(0, y), **common)
        lax.add_artist(top)
        y -= row_pt / h
    if rest:
        # matplotlib fills legend columns first; reorder so the entries read row by row
        n = len(rest)
        order = [i for col in range(columns) for i in range(col, n, columns)]
        rest_ordered = [rest[i] for i in order]
        lax.legend(
            [e[0] for e in rest_ordered],
            [e[1] for e in rest_ordered],
            loc="upper left",
            bbox_to_anchor=(0, y),
            ncol=columns,
            **common,
        )
    return h
