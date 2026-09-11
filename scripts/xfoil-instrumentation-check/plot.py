#!/usr/bin/env python3
"""Draw the instrumentation-check figure from a run folder: a 2×2 figure of CL, CD and CM against
α (pristine red with crosses, instrumented blue with open circles) and the upper-surface
trailing-edge boundary layer against α on a logarithmic axis (one colour per variable, pristine
solid, instrumented dashed).

Presentation only: `summary.json` holds each build's plotted series (the unbroken converged run,
selected by the Rust program) and the extents of every metric over those series. Run via
scripts/figures/render.sh.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402

COEFFICIENTS = [("cl", r"$C_L$"), ("cd", r"$C_D$"), ("cm", r"$C_M$")]
# upper-surface trailing-edge variables of the fourth panel (Cf omitted: ≤ 0 once separated)
BL_VARS = [
    ("te_upper.ue", r"$U_e/V_\infty$"),
    ("te_upper.dstar", r"$\delta^*$"),
    ("te_upper.theta", r"$\theta$"),
    ("te_upper.hk", r"$H_k$"),
]
BUILD = {
    "pristine": dict(color=style.RED, marker="x", linestyle="-"),
    "instrumented": dict(color=style.BLUE, marker="o", linestyle="--"),
}

LEGEND_ROW_PT = 11.0
CELL_H_PT = 172.0
Y_LABEL_AREA_PT = 34.0
GUTTER_PT = 16.0
X_LABEL_AREA_PT = 28.0
MARGIN_PT = 3.0
MARKER_PT = 4.0


def value(point: dict, key: str) -> float:
    for k in key.split("."):
        point = point[k]
    return point


def padded(extents: dict, key: str, log: bool = False):
    lo, hi = extents[key]["min"], extents[key]["max"]
    if log:
        return lo / 1.3, hi * 1.3
    pad = 0.06 * max(hi - lo, 1e-9)
    return lo - pad, hi + pad


def main():
    style.apply()
    run_dir = Path(sys.argv[1])
    s = json.loads((run_dir / "summary.json").read_text())
    builds = s["builds"]
    extents = s["extents"]
    legend_h = LEGEND_ROW_PT * len(builds) + 4.0
    width = style.TEXT_WIDTH_PT
    height = legend_h + 2 * CELL_H_PT
    fig = style.figure(width, height)
    alpha_max = extents["alpha_deg"]["max"]
    col_w = (width - GUTTER_PT) / 2
    plot_w = col_w - Y_LABEL_AREA_PT - 2 * MARGIN_PT - 4.0

    def cell(i, ylim, log=False, top_strip=0.0):
        r, c = divmod(i, 2)
        x0 = c * (col_w + GUTTER_PT) + Y_LABEL_AREA_PT + MARGIN_PT
        plot_h = CELL_H_PT - X_LABEL_AREA_PT - 2 * MARGIN_PT - top_strip
        y0 = height - legend_h - CELL_H_PT * (r + 1) + X_LABEL_AREA_PT + MARGIN_PT
        ax = style.axes(fig, x0, y0, plot_w, plot_h)
        ax.set_xlim(-0.5, alpha_max + 0.5)
        ax.set_ylim(*ylim)
        if log:
            ax.set_yscale("log")
        ax.set_xlabel(r"$\alpha$ (°)")
        return ax, x0, y0 + plot_h

    build_handles = []
    for i, (key, label) in enumerate(COEFFICIENTS):
        ax, _, _ = cell(i, padded(extents, key))
        ax.set_ylabel(label)
        for b in builds:
            st = BUILD[b["build"].lower()]
            a = [p["alpha_deg"] for p in b["series"]]
            v = [value(p, key) for p in b["series"]]
            (h,) = ax.plot(a, v, color=st["color"], linewidth=0.75, marker=st["marker"], markersize=MARKER_PT, markerfacecolor="none")
            if i == 0:
                build_handles.append(h)

    lo = min(extents[k]["min"] for k, _ in BL_VARS)
    hi = max(extents[k]["max"] for k, _ in BL_VARS)
    ax, x0, top = cell(3, (lo / 1.3, hi * 1.3), log=True, top_strip=LEGEND_ROW_PT)
    ax.set_ylabel("value, upper-surface TE")
    var_handles = []
    for k, (key, label) in enumerate(BL_VARS):
        colour = style.CATEGORY[k]
        for b in builds:
            st = BUILD[b["build"].lower()]
            a = [p["alpha_deg"] for p in b["series"]]
            v = [value(p, key) for p in b["series"]]
            (h,) = ax.plot(a, v, color=colour, linewidth=1.5, linestyle=st["linestyle"], marker=st["marker"], markersize=MARKER_PT, markerfacecolor="none")
        var_handles.append((h, label))
    style.legend_below(fig, x0, top + LEGEND_ROW_PT, plot_w, var_handles, len(var_handles), row_pt=LEGEND_ROW_PT)

    # legend at the top: one row per build
    entries = []
    for h, b in zip(build_handles, builds):
        end = b["series_end_alpha_deg"]
        end = f"{end:.1f}" if end is not None else "—"
        entries.append((h, f"{b['label']}: {len(b['points'])} converged points, series to {end}°"))
    style.legend_below(fig, 2.0, height, width - 4.0, entries, 1, row_pt=LEGEND_ROW_PT)

    style.save(fig, run_dir / f"naca{s['foil']}_check.svg")


if __name__ == "__main__":
    main()
