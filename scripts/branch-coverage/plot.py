#!/usr/bin/env python3
"""Draw the branch-coverage study's figures from a run folder's summary.json.

Presentation only: `branch_map` lists every analysis-path subroutine with the status of each of
its branches, `matrix` the fraction of each subroutine's branches every case of the minimal set
takes. Run via scripts/figures/render.sh branch-coverage <run_dir>.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402
from matplotlib import colormaps  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

STATUS_COLOUR = {
    "cover": "#1f5ad6",       # taken by the minimal set
    "finite": "#a6c8ff",      # taken only by finite cases outside the set
    "nonfinite": "#ff7f0e",   # taken only in non-finite reference runs
    "open": "#d62728",        # never taken
    "unreachable": "#9a9a9a", # annotated unreachable
    "loop_entry": "#dddddd",  # DO-line loop-control edge never taken
}
STATUS_LABEL = {
    "cover": "taken by the minimal set",
    "finite": "taken only outside the set",
    "nonfinite": "taken only in non-finite runs",
    "open": "open",
    "unreachable": "annotated unreachable",
    "loop_entry": "DO-line edge never taken (loop left early or never reached)",
}
CELL_PT = 3.3
ROW_PT = 5.4
LABEL_PT = 78.0
LEGEND_ROW_PT = 11.0
LEGEND_ROWS = 2
MARGIN_PT = 4.0
FILE_GAP_PT = 3.0


def branch_map(run_dir: Path, s: dict) -> None:
    rows = [r for r in s["branch_map"] if len(r["statuses"]) > 0]
    files = []
    for r in rows:
        if r["file"] not in files:
            files.append(r["file"])
    legend = LEGEND_ROW_PT * LEGEND_ROWS + 4
    height = legend + 2 * MARGIN_PT + ROW_PT * (len(rows) + len(files)) + FILE_GAP_PT * len(files)
    fig = style.figure(style.TEXT_WIDTH_PT, height)
    ax = style.axes(fig, MARGIN_PT, MARGIN_PT + legend, style.TEXT_WIDTH_PT - 2 * MARGIN_PT,
                    height - legend - 2 * MARGIN_PT)
    ax.set_axis_off()
    ax.set_xlim(0, style.TEXT_WIDTH_PT - 2 * MARGIN_PT)
    ax.set_ylim(0, height - legend - 2 * MARGIN_PT)
    ax.grid(False)
    y = height - legend - 2 * MARGIN_PT
    current_file = None
    for r in rows:
        if r["file"] != current_file:
            current_file = r["file"]
            y -= FILE_GAP_PT
            ax.text(0, y - ROW_PT * 0.15, current_file, fontsize=style.TICK_PT, fontweight="bold",
                    ha="left", va="top", family="monospace")
            y -= ROW_PT
        y -= ROW_PT
        ax.text(LABEL_PT - 4, y + ROW_PT / 2, r["name"], fontsize=style.TICK_PT - 1.5, ha="right",
                va="center", family="monospace")
        for i, st in enumerate(r["statuses"]):
            ax.add_patch(Rectangle((LABEL_PT + i * CELL_PT, y + 0.9), CELL_PT - 0.6, ROW_PT - 1.8,
                                   facecolor=STATUS_COLOUR[st], edgecolor="none"))
    # legend, LEGEND_ROWS rows of three entries
    lax = style.axes(fig, MARGIN_PT, MARGIN_PT, style.TEXT_WIDTH_PT - 2 * MARGIN_PT, legend)
    lax.set_axis_off()
    lax.grid(False)
    lax.set_xlim(0, style.TEXT_WIDTH_PT - 2 * MARGIN_PT)
    lax.set_ylim(0, legend)
    keys = ["cover", "finite", "nonfinite", "open", "unreachable", "loop_entry"]
    per_row = -(-len(keys) // LEGEND_ROWS)
    col_w = (style.TEXT_WIDTH_PT - 2 * MARGIN_PT) / per_row
    for n, key in enumerate(keys):
        x = (n % per_row) * col_w
        yy = legend - 4 - (n // per_row) * LEGEND_ROW_PT - LEGEND_ROW_PT / 2
        lax.add_patch(Rectangle((x, yy - 4), 8, 8, facecolor=STATUS_COLOUR[key], edgecolor="none"))
        lax.text(x + 11, yy, STATUS_LABEL[key], fontsize=style.LEGEND_PT, va="center")
    for ext in ("svg", "pdf"):
        fig.savefig(run_dir / f"branch-map.{ext}")


def case_matrix(run_dir: Path, s: dict) -> None:
    """Subroutines as rows, the cases of the set as columns."""
    m = s["matrix"]
    cases = m["cases"]
    keep = [i for i, row in enumerate(s["branch_map"]) if len(row["statuses"]) > 0]
    subs = [m["subroutines"][i] for i in keep]
    frac = [[m["fraction"][j][i] for i in keep] for j in range(len(cases))]
    unique = {c["name"]: set(c["unique_by_subroutine"].keys()) for c in s["cases"]}
    # three numeric columns to the right of the grid, per subroutine: every branch, the reachable
    # ones (not annotated unreachable, not a never-taken DO-line edge), and those taken on any
    # measured run (the set, the other finite candidates and the non-finite probes)
    rows = {r["key"]: r for r in s["subroutines"]}
    counts = []
    for key in subs:
        r = rows[key]
        total = r["branches"]
        reachable = total - r["unreachable"] - r["loop_entry"]
        taken_all = r["cover"] + r["finite"] + r["nonfinite"]
        counts.append((total, reachable, taken_all))
    col_names = ["branches", "reachable", "taken, all cases"]
    col_x = [16.0, 46.0, 76.0]
    cols_w = 84.0
    cell_h = 7.2
    left, top, bottom = 62.0, 168.0, 50.0
    width = style.TEXT_WIDTH_PT
    cell_w = (width - left - 4 - cols_w) / max(1, len(cases))
    height = top + bottom + cell_h * len(subs)
    fig = style.figure(width, height)
    ax = style.axes(fig, left, bottom, cell_w * len(cases), cell_h * len(subs))
    ax.grid(False)
    ax.set_xlim(0, len(cases))
    ax.set_ylim(0, len(subs))
    cmap = colormaps["YlGnBu"]
    for i, key in enumerate(subs):
        yy = len(subs) - 1 - i
        for j, name in enumerate(cases):
            f = frac[j][i]
            colour = cmap(0.15 + 0.85 * f) if f > 0 else "#f4f4f4"
            ax.add_patch(Rectangle((j, yy), 1, 1, facecolor=colour, edgecolor="white", linewidth=0.4))
            if key.split(":")[1] in unique.get(name, set()):
                ax.plot(j + 0.5, yy + 0.5, marker="o", markersize=2.2, color="black")
    ax.set_yticks([len(subs) - 1 - i + 0.5 for i in range(len(subs))])
    ax.set_yticklabels([k.split(":")[1] for k in subs], fontsize=style.TICK_PT - 2, family="monospace")
    ax.set_xticks([j + 0.5 for j in range(len(cases))])
    ax.set_xticklabels(cases, rotation=90, fontsize=style.TICK_PT - 1.5, family="monospace")
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # the numeric columns
    nax = style.axes(fig, left + cell_w * len(cases) + 2, bottom, cols_w, cell_h * len(subs) + top)
    nax.set_axis_off()
    nax.grid(False)
    nax.set_xlim(0, cols_w)
    nax.set_ylim(0, cell_h * len(subs) + top)
    for x, name in zip(col_x, col_names):
        nax.text(x, cell_h * len(subs) + 4, name, rotation=90, fontsize=style.TICK_PT - 1.5, ha="center",
                 va="bottom", family="monospace")
    for i, (total, reachable, taken_all) in enumerate(counts):
        yy = (len(subs) - 1 - i + 0.5) * cell_h
        for x, v in zip(col_x, (total, reachable, taken_all)):
            nax.text(x + 8, yy, str(v), fontsize=style.TICK_PT - 2, ha="right", va="center", family="monospace")
    # key: the colour scale, the "not entered" cell and the unique-branch dot
    key_y = 16.0
    bar_w = 150.0
    cax = style.axes(fig, left, key_y, bar_w, 7.0)
    cax.grid(False)
    for sp in cax.spines.values():
        sp.set_visible(False)
    steps = 40
    for k in range(steps):
        f = (k + 0.5) / steps
        cax.add_patch(Rectangle((k / steps, 0), 1 / steps, 1, facecolor=cmap(0.15 + 0.85 * f), edgecolor="none"))
    cax.set_xlim(0, 1)
    cax.set_ylim(0, 1)
    cax.set_yticks([])
    cax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    cax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=style.TICK_PT - 1)
    cax.tick_params(length=2, pad=1)
    fig.text((left) * style.PT / fig.get_size_inches()[0], (key_y + 20) * style.PT / fig.get_size_inches()[1],
             "fraction of the subroutine's branches the case took", fontsize=style.LEGEND_PT, va="center")
    kax = style.axes(fig, left + bar_w + 24, key_y - 2, width - left - bar_w - 30, 24)
    kax.set_axis_off()
    kax.grid(False)
    kax.set_xlim(0, width - left - bar_w - 30)
    kax.set_ylim(0, 24)
    kax.add_patch(Rectangle((0, 12), 9, 9, facecolor="#f4f4f4", edgecolor="#cccccc", linewidth=0.4))
    kax.text(13, 16.5, "subroutine not entered by the case", fontsize=style.LEGEND_PT, va="center")
    kax.add_patch(Rectangle((0, 0), 9, 9, facecolor=cmap(0.6), edgecolor="none"))
    kax.plot(4.5, 4.5, marker="o", markersize=2.2, color="black")
    kax.text(13, 4.5, "takes a branch no other finite candidate takes",
             fontsize=style.LEGEND_PT, va="center")
    for ext in ("svg", "pdf"):
        fig.savefig(run_dir / f"case-matrix.{ext}")


def main() -> None:
    run_dir = Path(sys.argv[1])
    style.apply()
    s = json.loads((run_dir / "summary.json").read_text())
    branch_map(run_dir, s)
    case_matrix(run_dir, s)
    print(f"figures written to {run_dir}")


if __name__ == "__main__":
    main()
