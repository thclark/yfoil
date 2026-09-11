#!/usr/bin/env python3
"""Draw the aerofoil-series figures from a run folder: one figure per family, one panel per
parameter, the family default in black on top, the members behind it coloured by parameter value
(YlGnBu, lightest = smallest), every panel on the same equal-scale axes.

Presentation only: `summary.json` names every section and its geometry file (already ordered by
parameter value), and the geometry files hold the coordinates. Run via scripts/figures/render.sh.
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402

X_RANGE = (-0.2, 1.2)
Y_RANGE = (-0.5, 0.5)
COLUMNS = 2
GUTTER_PT = 14.0
Y_LABEL_AREA_PT = 30.0
X_LABEL_AREA_PT = 24.0
TITLE_PT = 13.0
LEGEND_ROW_PT = 10.0
LEGEND_COLS = 3
LEGEND_PAD_PT = 4.0
MARGIN_PT = 3.0
RIGHT_TICK_PT = 6.0


def geometry(path: Path):
    d = json.loads(path.read_text())
    return d["x"], d["y"]


def draw_series(run_dir: Path, series: dict) -> None:
    panels = series["panels"]
    n = len(panels)
    rows = math.ceil(n / COLUMNS)
    col_w = (style.TEXT_WIDTH_PT - GUTTER_PT * (COLUMNS - 1)) / COLUMNS
    plot_w = col_w - Y_LABEL_AREA_PT - 2 * MARGIN_PT - RIGHT_TICK_PT
    plot_h = plot_w * (Y_RANGE[1] - Y_RANGE[0]) / (X_RANGE[1] - X_RANGE[0])
    legend_rows = [1 + math.ceil(len(p["members"]) / LEGEND_COLS) for p in panels]
    cell_h = [TITLE_PT + plot_h + X_LABEL_AREA_PT + 2 * MARGIN_PT + LEGEND_PAD_PT + LEGEND_ROW_PT * r + 4.0 for r in legend_rows]
    row_h = [max(cell_h[r * COLUMNS : (r + 1) * COLUMNS]) for r in range(rows)]
    fig_h = sum(row_h)
    fig = style.figure(style.TEXT_WIDTH_PT, fig_h)

    default_x, default_y = geometry(run_dir / series["slug"] / series["default"]["file"])
    y_top = fig_h
    for k, panel in enumerate(panels):
        r, c = divmod(k, COLUMNS)
        if c == 0 and k > 0:
            y_top -= row_h[r - 1]
        x0 = c * (col_w + GUTTER_PT) + Y_LABEL_AREA_PT + MARGIN_PT
        y0 = y_top - TITLE_PT - plot_h - MARGIN_PT
        ax = style.axes(fig, x0, y0, plot_w, plot_h)
        ax.set_xlim(*X_RANGE)
        ax.set_ylim(*Y_RANGE)
        ax.set_aspect("equal")
        ax.set_xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
        ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])
        ax.set_xlabel(r"$x/c$")
        if c == 0:
            ax.set_ylabel(r"$y/c$")
        ax.set_title(panel["title"], loc="left", pad=3)

        members = panel["members"]  # already ordered by parameter value
        nm = max(len(members), 2) - 1
        entries = []
        for i, m in enumerate(members):
            x, y = geometry(run_dir / series["slug"] / panel["slug"] / m["file"])
            (h,) = ax.plot(x, y, color=style.ylgnbu(i / nm), linewidth=0.75)
            entries.append((h, m["label"]))
        (hd,) = ax.plot(default_x, default_y, color="black", linewidth=1.5)
        style.legend_below(
            fig,
            x0,
            y0 - X_LABEL_AREA_PT - LEGEND_PAD_PT,
            plot_w,
            [(hd, f"{series['default']['label']} (default)")] + entries,
            LEGEND_COLS,
            row_pt=LEGEND_ROW_PT,
            first_alone=True,
        )

    style.save(fig, run_dir / f"{series['slug']}.svg")


def main():
    style.apply()
    run_dir = Path(sys.argv[1])
    summary = json.loads((run_dir / "summary.json").read_text())
    for series in summary["series"]:
        draw_series(run_dir, series)


if __name__ == "__main__":
    main()
