#!/usr/bin/env python3
"""Draw the XFOIL input-sensitivity figures from a run folder: seven rows of per-alpha
quantities, one column per input family, the base case in front (black, filled markers) and the
perturbation levels behind it coloured by perturbation size (YlGnBu, lightest = smallest), open
circles on unconverged points. The node and alpha-step families share a two-column sheet; every
family also gets a one-column figure of the same column width.

Presentation only: `metrics.json` holds every level with its points and the converged extents of
every quantity (the axis ranges come from those), `metadata.json` the base case. Run via
scripts/figures/render.sh.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import style  # noqa: E402

# rows, top to bottom: metrics.json key, axis label, log axis
ROWS = [
    ("cl", r"$C_L$", False),
    ("cd", r"$C_D$", False),
    ("dstar_te_upper", r"$\delta^*$ at TE", False),
    ("h_te_upper", r"$H$ at TE", False),
    ("xtr_upper", r"$x/c$ transition", False),
    ("x_sep_upper", r"$x/c$ separation", False),
    ("rmsbl", "RMSBL", True),
]
FAMILIES = {
    "geometry": "Node coordinates perturbed",
    "panels": "Panel count reduced",
    "alpha-step": "Alpha step perturbed",
}
SHEET = ["geometry", "alpha-step"]

GUTTER_PT = 16.0
ROW_H_PT = 82.0
X_LABEL_AREA_PT = 28.0
X_STUB_PT = 4.0
Y_LABEL_AREA_PT = 34.0
MARGIN_PT = 3.0
LEGEND_ROW_PT = 10.0
LEGEND_TITLE_PT = 12.0
LEGEND_COLS = 2
MARKER_PT = 2.4
OPEN_MARKER_PT = 4.4


def exponent(v: float) -> str:
    return str(-round(math.log10(v)))


def level_label(level: dict, meta: dict) -> str:
    """Legend text from the level's fields."""
    base = meta["base"]
    if level["magnitude"] == 0.0:
        return rf"base: $N$ = {base['panels']}, $\Delta\alpha$ = {base['alpha_step_deg']}°"
    ulp = meta["ulp_marker"]
    fam = level["family"]
    if fam == "geometry":
        eps = level["geometry_eps"]
        text = "nodes ± 1 ULP" if eps == ulp else rf"nodes ± $10^{{-{exponent(eps)}}}$"
    elif fam == "panels":
        n, req = level["actual_panels"], level["panels"]
        text = rf"$N$ = {n}" + (f" (req. {req})" if n != req else "")
    else:
        d = level["raw"]
        text = r"$\Delta\alpha$ + 1 ULP" if d == ulp else rf"$\Delta\alpha + 10^{{-{exponent(d)}}}{{}}^\circ$"
    if level["ending"] != "Completed":
        text += " (hung)"
    return text


def ordered(levels):
    """Draw order: the base first (drawn last, on top), then levels by increasing magnitude."""
    base = [lv for lv in levels if lv["magnitude"] == 0.0]
    rest = sorted((lv for lv in levels if lv["magnitude"] > 0.0), key=lambda lv: lv["magnitude"])
    return base + rest


def y_limits(extents: dict, key: str, log: bool):
    lo, hi = extents[key]["min"], extents[key]["max"]
    if log:
        return lo / 10.0, hi * 1e4
    pad = 0.5 * max(hi - lo, 1e-12)
    return lo - pad, hi + pad


def column_width_pt() -> float:
    return (style.TEXT_WIDTH_PT - GUTTER_PT * (len(SHEET) - 1)) / len(SHEET)


def legend_rows(all_levels) -> int:
    n = max(len(v) for v in all_levels.values())
    return 1 + math.ceil((n - 1) / LEGEND_COLS)


def draw(run_dir: Path, foil: str, meta: dict, metrics: dict, families: list, name: str) -> None:
    levels_by_family = metrics["levels_by_family"]
    extents = metrics["extents"]
    alpha_max = extents["alpha_deg"]["max"]
    ncol = len(families)
    col_w = column_width_pt()
    legend_h = LEGEND_TITLE_PT + LEGEND_ROW_PT * legend_rows(levels_by_family) + 4.0
    width = col_w * ncol + GUTTER_PT * (ncol - 1)
    height = legend_h + ROW_H_PT * len(ROWS) + X_LABEL_AREA_PT - X_STUB_PT
    fig = style.figure(width, height)

    for c, family in enumerate(families):
        levels = ordered(levels_by_family[family])
        x_left = c * (col_w + GUTTER_PT)
        n = max(len(levels) - 1, 1) - 1
        handles = []
        for r, (key, label, log) in enumerate(ROWS):
            bottom = r == len(ROWS) - 1
            x_area = X_LABEL_AREA_PT if bottom else X_STUB_PT
            row_top = height - legend_h - ROW_H_PT * r
            plot_h = ROW_H_PT - x_area - MARGIN_PT * 2 + (0 if bottom else 0)
            ax = style.axes(
                fig,
                x_left + Y_LABEL_AREA_PT + MARGIN_PT,
                row_top - MARGIN_PT - plot_h,
                col_w - Y_LABEL_AREA_PT - 2 * MARGIN_PT - 4.0,
                plot_h,
            )
            lo, hi = y_limits(extents, key, log)
            ax.set_xlim(-0.5, alpha_max + 0.5)
            ax.set_ylim(lo, hi)
            if log:
                ax.set_yscale("log")
            ax.set_ylabel(label)
            if bottom:
                ax.set_xlabel(r"$\alpha$ (°)")
            else:
                ax.tick_params(labelbottom=False)
            # stack, largest perturbation first so the base lands on top
            row_handles = []
            for i, lv in reversed(list(enumerate(levels))):
                is_base = lv["magnitude"] == 0.0
                colour = "black" if is_base else style.ylgnbu((i - 1) / n)
                pts = lv["points"]
                a = np.array([p["alpha_deg"] for p in pts])
                v = np.array([p[key] for p in pts], dtype=float)
                conv = np.array([p["converged"] for p in pts])
                ok = np.isfinite(v) & ((v > 0) if log else True)
                # values beyond the converged range are drawn on the axis edge (presentation)
                vc = np.clip(v, lo, hi)
                (h,) = ax.plot(a[ok], vc[ok], color=colour, linewidth=1.5 if is_base else 0.75)
                if is_base:
                    ax.plot(a[ok], vc[ok], "o", color=colour, markersize=MARKER_PT, linestyle="none")
                unc = ok & ~conv
                ax.plot(a[unc], vc[unc], "o", markerfacecolor="none", markeredgecolor=colour, markersize=OPEN_MARKER_PT, linestyle="none")
                row_handles.append((i, h))
            if r == 0:
                handles = [h for _, h in sorted(row_handles)]

        # legend above the column: the title, the base on its own row, the levels in columns
        entries = [(h, level_label(lv, meta)) for h, lv in zip(handles, levels)]
        strip_top = height - LEGEND_TITLE_PT
        fig.text(
            (x_left + Y_LABEL_AREA_PT + MARGIN_PT) * style.PT / fig.get_size_inches()[0],
            (height - LEGEND_TITLE_PT + 3.0) * style.PT / fig.get_size_inches()[1],
            FAMILIES[family],
            fontsize=style.AXIS_LABEL_PT,
        )
        style.legend_below(
            fig,
            x_left + Y_LABEL_AREA_PT + MARGIN_PT,
            strip_top,
            col_w - Y_LABEL_AREA_PT - 2 * MARGIN_PT,
            entries,
            LEGEND_COLS,
            row_pt=LEGEND_ROW_PT,
            first_alone=True,
        )

    style.save(fig, run_dir / f"naca{foil}_{name}.svg")


def main():
    style.apply()
    run_dir = Path(sys.argv[1])
    meta = json.loads((run_dir / "metadata.json").read_text())
    metrics = json.loads((run_dir / "metrics.json").read_text())
    foil = meta["foil"].removeprefix("naca")
    families = list(metrics["levels_by_family"])
    for family in families:
        draw(run_dir, foil, meta, metrics, [family], family)
    if all(f in families for f in SHEET):
        draw(run_dir, foil, meta, metrics, SHEET, "sheet")


if __name__ == "__main__":
    main()
