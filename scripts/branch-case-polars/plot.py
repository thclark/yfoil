#!/usr/bin/env python3
"""Draw the branch-case-polars figures from a run folder's summary.json, one per case: CL and CD
(top), the last iteration's RMSBL and the iterations taken (bottom) against alpha, XFOIL in red
(dashed, crosses), yFoil in blue (solid, circles), non-finite behaviour in orange, departures as
diamonds.

Presentation only: every number comes from `summary.json`. Run via
scripts/figures/render.sh branch-case-polars <run_dir>.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402

XFOIL_COLOUR = "#d62728"
YFOIL_COLOUR = "#1f5ad6"
NONFINITE_COLOUR = "#ff7f0e"


def polars(run_dir: Path, s: dict) -> None:
    """One figure per case: CL, CD (top) and the last RMSBL, iterations (bottom) against alpha,
    XFOIL in red (dashed, crosses), yFoil in blue (solid, circles filled where yFoil followed the
    reference to the end of the point); an orange outer square where that code went through a
    non-finite value within the point, an orange tick at the top of the panel where the plotted
    value is itself NaN; a black outer diamond where the runs part and the reference is unstable
    against its 1-ULP twins, a filled grey diamond where they part and it is stable."""
    for run in s.get("polars", []):
        polar_figure(run_dir / f"polars-{run['name']}", run)


def polar_figure(stem: Path, run: dict) -> None:
    from matplotlib.lines import Line2D
    width = style.TEXT_WIDTH_PT
    gutter, left, bottom, top = 34.0, 44.0, 28.0, 14.0
    legend_h = 11.0 * 8 + 8.0
    panel_w = (width - left - gutter - 8) / 2
    panel_h = 118.0
    # the state row only for a case that kept the reference's per-call state records
    has_state = any(p.get("state") for p in run["points"])
    panels = [("cl", "$C_L$", False), ("cd", "$C_D$", False),
              ("rmsbl", "RMSBL at the last iteration", True),
              ("iterations", "iterations", False),
              ("xtr_upper", "$x_{tr}/c$, upper side", False),
              ("xtr_lower", "$x_{tr}/c$, lower side", False)]
    panels += [("trpos_upper", "transition position in its interval, upper", False),
               ("trpos_lower", "transition position in its interval, lower", False)]
    if has_state:
        panels += [("state_ratio", "worst state array: $|\\Delta|$ / gate", True),
                   ("state_diff", "worst state array: $|\\Delta|$, floor (red)", True)]
    rows = (len(panels) + 1) // 2
    height = top + rows * panel_h + bottom + (rows - 1) * 30.0 + legend_h
    fig = style.figure(width, height)
    axes = {}
    for r, (key, ylabel, log) in enumerate(panels):
        col = r % 2
        row = r // 2
        x0 = left + col * (panel_w + gutter)
        y0 = height - top - (row + 1) * panel_h - row * 30.0
        ax = style.axes(fig, x0, y0, panel_w, panel_h)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(r"$\alpha$ / deg")
        ax.set_xlim(-30, 30)
        if log:
            ax.set_yscale("log")
        axes[key] = ax
    fig.text(0.5, 1 - 4.0 * style.PT / fig.get_size_inches()[1], f"{run['section']}: {run['name']}",
             ha="center", va="top", fontsize=style.AXIS_LABEL_PT, family="monospace")
    nan_ticks = {key: [] for key in axes}

    def value(p, code, key):
        if key in ("trpos_upper", "trpos_lower"):
            # yFoil's state only (the reference agrees in ITRAN wherever the point was followed)
            if code == "xfoil":
                return None
            t = p["transition_position"][0 if key == "trpos_upper" else 1]
            return None if t is None or t["at_te"] else t["fraction"]
        if key == "state_ratio":
            # the reference has no ratio: its own twin's spread is the floor the gate is built on
            return (p["state"][0]["ratio"] if p.get("state") else None) if code == "yfoil" else None
        if key == "state_diff":
            if not p.get("state"):
                return None
            return p["state"][0]["diff"] if code == "yfoil" else p["state"][0]["floor"]
        if key == "xtr_upper":
            return (p["xfoil_transition"][0] if p.get("xfoil_transition") else None) if code == "xfoil" else p["yfoil_transition"][0]
        if key == "xtr_lower":
            return (p["xfoil_transition"][1] if p.get("xfoil_transition") else None) if code == "xfoil" else p["yfoil_transition"][1]
        return p[code][key]

    for leg in (1, 2):
        pts = [p for p in run["points"] if p["leg"] == leg]
        pts.sort(key=lambda p: p["alpha_deg"])
        a = [p["alpha_deg"] for p in pts]
        for key, ax in axes.items():
            # the reference first, so that yFoil's markers sit on top where they coincide
            # the reference has no series in the ratio panel (its twin's spread is the gate)
            xp = [] if key in ("state_ratio", "trpos_upper", "trpos_lower") else [p for p in pts if p["xfoil"] is not None]
            xa = [p["alpha_deg"] for p in xp]
            xy = [value(p, "xfoil", key) for p in xp]
            ax.plot(xa, xy, "--", color=XFOIL_COLOUR, linewidth=0.9, marker="x", markersize=3.2, markeredgewidth=0.7)
            for p, yy in zip(xp, xy):
                if yy is None or yy != yy:
                    nan_ticks[key].append((p["alpha_deg"], XFOIL_COLOUR))
                    continue
                if p.get("xfoil_nonfinite"):
                    ax.plot(p["alpha_deg"], yy, marker="s", markersize=5.2, color=NONFINITE_COLOUR,
                            markerfacecolor="none", markeredgewidth=0.7, linestyle="none")
            y = [value(p, "yfoil", key) for p in pts]
            ax.plot(a, y, "-", color=YFOIL_COLOUR, linewidth=0.8)
            for p, yy in zip(pts, y):
                filled = p.get("followed", False)
                if yy is None or yy != yy:
                    nan_ticks[key].append((p["alpha_deg"], YFOIL_COLOUR))
                    continue
                ax.plot(p["alpha_deg"], yy, marker="o", markersize=2.4, color=YFOIL_COLOUR,
                        markerfacecolor=YFOIL_COLOUR if filled else "none", markeredgewidth=0.6, linestyle="none")
                if p.get("yfoil_nonfinite"):
                    ax.plot(p["alpha_deg"], yy, marker="s", markersize=5.2, color=NONFINITE_COLOUR,
                            markerfacecolor="none", markeredgewidth=0.7, linestyle="none")
                if p.get("outcome") == "divergent" and not filled:
                    ax.plot(p["alpha_deg"], yy, marker="D", markersize=5.8, color="black",
                            markerfacecolor="none", markeredgewidth=0.5, linestyle="none")
                elif p.get("outcome") == "mismatch":
                    ax.plot(p["alpha_deg"], yy, marker="D", markersize=5.8, color="black",
                            markerfacecolor="black", markeredgewidth=0.5, linestyle="none", alpha=0.35)
                if key in ("xtr_upper", "xtr_lower") and p.get("stations_differ"):
                    # the final stagnation or transition station differs between the codes
                    ax.plot(p["alpha_deg"], yy, marker="v", markersize=6.5, color="black",
                            markerfacecolor="none", markeredgewidth=0.5, linestyle="none")
    for key, ax in axes.items():
        if key == "rmsbl":
            ax.set_ylim(bottom=1e-8)
        if key == "state_ratio":
            # the gate itself: above it the array is outside what the twin allows
            ax.axhline(1.0, color="black", linewidth=0.5, linestyle=":")
        if key in ("trpos_upper", "trpos_lower"):
            ax.set_ylim(-0.02, 1.02)
            for y in (0.0, 0.5, 1.0):
                ax.axhline(y, color="black", linewidth=0.4, linestyle=":")
        for alpha, colour in nan_ticks[key]:
            ax.plot(alpha, 0.985, marker="|", markersize=4, color=NONFINITE_COLOUR, markeredgewidth=0.8,
                    transform=ax.get_xaxis_transform(), linestyle="none", clip_on=False)
    handles = [
        Line2D([0], [0], color=XFOIL_COLOUR, linestyle="--", marker="x", markersize=3.2, linewidth=0.9),
        Line2D([0], [0], color=YFOIL_COLOUR, linestyle="-", marker="o", markersize=2.4, markerfacecolor=YFOIL_COLOUR, linewidth=0.8),
        Line2D([0], [0], color=YFOIL_COLOUR, linestyle="-", marker="o", markersize=2.4, markerfacecolor="none", linewidth=0.8),
        Line2D([0], [0], color=NONFINITE_COLOUR, linestyle="none", marker="s", markersize=5.2, markerfacecolor="none"),
        Line2D([0], [0], color=NONFINITE_COLOUR, linestyle="none", marker="|", markersize=4),
        Line2D([0], [0], color="black", linestyle="none", marker="D", markersize=5.8, markerfacecolor="none"),
        Line2D([0], [0], color="black", linestyle="none", marker="D", markersize=5.8, markerfacecolor="black", alpha=0.35),
        Line2D([0], [0], color="black", linestyle="none", marker="v", markersize=6.5, markerfacecolor="none"),
    ]
    labels = [
        "XFOIL (reference)",
        "yFoil, followed to the end of the point",
        "yFoil, not followed, or not gated after an earlier departure",
        "went through a non-finite value within the point",
        "plotted value is NaN (tick at its α)",
        "threshold: run departs where the reference is unstable against its 1-ULP twins",
        "mismatch: run departs where the reference is stable against its 1-ULP twins",
        "final stagnation or transition station (IST, ITRAN) differs between the codes",
    ]
    fig.legend(handles, labels, loc="lower left", ncol=1, fontsize=style.LEGEND_PT, frameon=False,
               bbox_to_anchor=(left * style.PT / fig.get_size_inches()[0], 0.0))
    for ext in ("svg", "pdf"):
        fig.savefig(f"{stem}.{ext}")


def transition_positions(stem: Path, s: dict) -> None:
    """Every recorded point of every case: where transition sits within its station interval
    (0 = at the station before, 1 = at the transition station), one row per case, one panel per
    side; the points at which the run parts drawn as diamonds (black where the reference is
    unstable against its 1-ULP twins, grey where it is stable)."""
    from matplotlib.lines import Line2D
    width = style.TEXT_WIDTH_PT
    runs = s.get("polars", [])
    left, right, top, bottom, gutter = 206.0, 8.0, 14.0, 26.0, 16.0
    row_h = 11.0
    panel_h = row_h * max(len(runs), 1) + 8.0
    panel_w = (width - left - right - gutter) / 2
    legend_h = 11.0 * 4 + 6.0
    height = top + panel_h + bottom + legend_h
    fig = style.figure(width, height)
    for col, (side, label) in enumerate([(0, "upper side"), (1, "lower side")]):
        ax = style.axes(fig, left + col * (panel_w + gutter), height - top - panel_h, panel_w, panel_h)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.5, len(runs) - 0.5)
        ax.set_xlabel(f"position of transition in its station interval, {label}" if col == 0 else label)
        ax.set_yticks(range(len(runs)))
        ax.set_yticklabels([r["name"] for r in runs] if col == 0 else [""] * len(runs), family="monospace")
        ax.invert_yaxis()
        for x in (0.0, 0.5, 1.0):
            ax.axvline(x, color="black", linewidth=0.4, linestyle=":")
        for row, run in enumerate(runs):
            for p in run["points"]:
                t = p["transition_position"][side]
                if t is None or t["at_te"]:
                    continue
                x = t["fraction"]
                if p.get("followed"):
                    ax.plot(x, row, marker="o", markersize=2.0, color=YFOIL_COLOUR, linestyle="none", alpha=0.6)
                elif p.get("outcome") == "mismatch":
                    ax.plot(x, row, marker="D", markersize=5.0, color="black", markerfacecolor="black",
                            markeredgewidth=0.5, linestyle="none", alpha=0.35)
                elif p.get("outcome") == "divergent" and p.get("parted_iteration"):
                    ax.plot(x, row, marker="D", markersize=5.0, color="black", markerfacecolor="none",
                            markeredgewidth=0.5, linestyle="none")
                else:
                    ax.plot(x, row, marker="o", markersize=2.0, color=YFOIL_COLOUR, markerfacecolor="none",
                            markeredgewidth=0.5, linestyle="none", alpha=0.6)
    handles = [
        Line2D([0], [0], color=YFOIL_COLOUR, linestyle="none", marker="o", markersize=2.0, alpha=0.6),
        Line2D([0], [0], color=YFOIL_COLOUR, linestyle="none", marker="o", markersize=2.0, markerfacecolor="none", alpha=0.6),
        Line2D([0], [0], color="black", linestyle="none", marker="D", markersize=5.0, markerfacecolor="none"),
        Line2D([0], [0], color="black", linestyle="none", marker="D", markersize=5.0, markerfacecolor="black", alpha=0.35),
    ]
    labels = [
        "point followed to the end",
        "point not followed, or not gated after an earlier departure",
        "threshold: run departs where the reference is unstable against its 1-ULP twins",
        "mismatch: run departs where the reference is stable against its 1-ULP twins",
    ]
    fig.legend(handles, labels, loc="lower left", ncol=1, fontsize=style.LEGEND_PT, frameon=False,
               bbox_to_anchor=(0.0, 0.0))
    for ext in ("svg", "pdf"):
        fig.savefig(f"{stem}.{ext}")


def main() -> None:
    run_dir = Path(sys.argv[1])
    style.apply()
    s = json.loads((run_dir / "summary.json").read_text())
    polars(run_dir, s)
    transition_positions(run_dir / "transition-position", s)
    print(f"figures written to {run_dir}")


if __name__ == "__main__":
    main()
