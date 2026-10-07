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


# ---------------------------------------------------------------------------------------------
# The twins mode: one sheet per case from twins.json
# ---------------------------------------------------------------------------------------------

TWIN_ROWS = [
    ("cl", r"$C_L$"),
    ("cd", r"$C_D$"),
    ("xtr_upper", r"$x/c$ transition, upper"),
    ("dstar_te_upper", r"$\delta^*$ at TE, upper"),
]
ENVELOPE_COLOUR = "#d62728"
WARNING_COLOUR = "#ff7f0e"
TWIN_ROW_H_PT = 96.0
TWIN_GUTTER_PT = 30.0


def stitched(points):
    """The points of both legs as one curve ascending in alpha (leg 2 reversed, then leg 1)."""
    leg2 = sorted((p for p in points if p["leg"] == 2), key=lambda p: p["alpha_deg"])
    leg1 = sorted((p for p in points if p["leg"] == 1), key=lambda p: p["alpha_deg"])
    # the seed alpha of leg 2 repeats leg 1's first point: keep leg 1's
    if leg2 and leg1 and abs(leg2[-1]["alpha_deg"] - leg1[0]["alpha_deg"]) < 1e-6:
        leg2 = leg2[:-1]
    return leg2 + leg1


def twins_sheet(run_dir: Path, case: dict) -> None:
    from matplotlib.lines import Line2D
    name = case["name"]
    ref = case["reference"]["points"]
    twins = case["twins"]
    env = {(e["leg"], round(e["alpha_deg"] * 1000)): e for e in case["envelope"]}
    ext = case["extents"]
    nseeds = max(len(twins), 1)
    colour_of = {t["seed"]: style.ylgnbu(i / max(nseeds - 1, 1)) for i, t in enumerate(twins)}

    rows = [tuple(r) for r in case.get("plot_rows", TWIN_ROWS)]
    col_w = (style.TEXT_WIDTH_PT - TWIN_GUTTER_PT) / 2
    legend_h = LEGEND_TITLE_PT + LEGEND_ROW_PT * (1 + math.ceil((nseeds + 4) / 2)) + 4.0
    height = legend_h + TWIN_ROW_H_PT * len(rows) + X_LABEL_AREA_PT - X_STUB_PT
    fig = style.figure(style.TEXT_WIDTH_PT, height)
    a_lo, a_hi = ext["alpha_deg"]["min"], ext["alpha_deg"]["max"]
    handles = {}
    for r, (key, label) in enumerate(rows):
        bottom = r == len(rows) - 1
        x_area = X_LABEL_AREA_PT if bottom else X_STUB_PT
        row_top = height - legend_h - TWIN_ROW_H_PT * r
        plot_h = TWIN_ROW_H_PT - x_area - MARGIN_PT * 2
        for c in range(2):
            x_left = c * (col_w + TWIN_GUTTER_PT)
            ax = style.axes(fig, x_left + Y_LABEL_AREA_PT + MARGIN_PT, row_top - MARGIN_PT - plot_h,
                            col_w - Y_LABEL_AREA_PT - 2 * MARGIN_PT - 4.0, plot_h)
            ax.set_xlim(a_lo - 0.5, a_hi + 0.5)
            if bottom:
                ax.set_xlabel(r"$\alpha$ (°)")
            else:
                ax.tick_params(labelbottom=False)
            if c == 0:
                # the baseline: twins beneath (higher seed drawn first), the reference on top
                lo, hi = ext[key]["min"], ext[key]["max"]
                if lo is None or hi is None:
                    lo, hi = 0.0, 1.0
                pad = 0.5 * max(hi - lo, 1e-12)
                lo, hi = lo - pad, hi + pad
                ax.set_ylim(lo, hi)
                ax.set_ylabel(label)
                for t in reversed(twins):
                    pts = stitched(t["points"])
                    a = np.array([p["alpha_deg"] for p in pts])
                    v = np.array([p[key] for p in pts], dtype=float)
                    ok = np.isfinite(v)
                    conv = np.array([p["converged"] for p in pts])
                    (h,) = ax.plot(a[ok], np.clip(v[ok], lo, hi), color=colour_of[t["seed"]], linewidth=0.75)
                    unc = ok & ~conv
                    ax.plot(a[unc], np.clip(v[unc], lo, hi), "o", markerfacecolor="none",
                            markeredgecolor=colour_of[t["seed"]], markersize=OPEN_MARKER_PT, linestyle="none")
                    handles.setdefault(("twin", t["seed"]), h)
                pts = stitched(ref)
                a = np.array([p["alpha_deg"] for p in pts])
                v = np.array([p[key] for p in pts], dtype=float)
                ok = np.isfinite(v)
                conv = np.array([p["converged"] for p in pts])
                capped = np.array([p["capped"] for p in pts])
                (h,) = ax.plot(a[ok], np.clip(v[ok], lo, hi), color="black", linewidth=1.5)
                ax.plot(a[ok], np.clip(v[ok], lo, hi), "o", color="black", markersize=MARKER_PT, linestyle="none")
                unc = ok & ~conv
                ax.plot(a[unc], np.clip(v[unc], lo, hi), "o", markerfacecolor="none", markeredgecolor="black",
                        markersize=OPEN_MARKER_PT, linestyle="none")
                handles.setdefault("reference", h)
                # the reference's completion: an orange tick at the top where it stopped at its limit
                for aa in a[capped]:
                    ax.plot(aa, 0.985, marker="|", markersize=4, color=WARNING_COLOUR, markeredgewidth=0.9,
                            transform=ax.get_xaxis_transform(), linestyle="none", clip_on=False)
            else:
                # the error column: |twin − reference| per twin, the envelope in red
                dlo, dhi = ext[key]["diff_min"], ext[key]["diff_max"]
                if dlo is None or dhi is None or not (np.isfinite(dlo) and np.isfinite(dhi)):
                    dlo, dhi = 1e-16, 1.0
                ax.set_yscale("log")
                ax.set_ylim(dlo / 10.0, dhi * 10.0)
                ax.set_ylabel(r"$|\Delta|$ " + label)
                for t in twins:
                    pts = stitched(t["points"])
                    a = np.array([p["alpha_deg"] for p in pts])
                    d = np.array([p["diff"][key] for p in pts], dtype=float)
                    ok = np.isfinite(d) & (d > 0)
                    (h,) = ax.plot(a[ok], d[ok], color=colour_of[t["seed"]], linewidth=0.75)
                    differs = np.array([p["completion_differs"] for p in pts])
                    m = ok & differs
                    ax.plot(a[m], d[m], marker="D", markersize=5.0, color=WARNING_COLOUR, markerfacecolor="none",
                            markeredgewidth=0.9, linestyle="none")
                    # a zero difference has no place on a log axis: a tick at the bottom
                    z = np.isfinite(d) & (d == 0)
                    ax.plot(a[z], np.full(z.sum(), 0.015), marker="|", markersize=3, color=colour_of[t["seed"]],
                            transform=ax.get_xaxis_transform(), linestyle="none", clip_on=False)
                    handles.setdefault(("twin", t["seed"]), h)
                pts = stitched(case["envelope"])
                a = np.array([p["alpha_deg"] for p in pts])
                d = np.array([p["scalars"][key] for p in pts], dtype=float)
                ok = np.isfinite(d) & (d > 0)
                (h,) = ax.plot(a[ok], d[ok], color=ENVELOPE_COLOUR, linewidth=1.3)
                handles.setdefault("envelope", h)
                # a non-finite envelope (a twin non-finite where the reference was not): a red tick
                nan = ~np.isfinite(d)
                for aa in a[nan]:
                    ax.plot(aa, 0.985, marker="|", markersize=4, color=ENVELOPE_COLOUR, markeredgewidth=0.9,
                            transform=ax.get_xaxis_transform(), linestyle="none", clip_on=False)
    entries = [(handles["reference"], "reference")]
    for t in twins:
        entries.append((handles[("twin", t["seed"])], f"twin, seed {t['seed']}"))
    entries.append((handles["envelope"], "envelope (max over the twins)"))
    entries.append((Line2D([0], [0], color=WARNING_COLOUR, marker="D", markerfacecolor="none", linestyle="none", markersize=5.0),
                    "twin finished differently from the reference"))
    entries.append((Line2D([0], [0], color=WARNING_COLOUR, marker="|", linestyle="none", markersize=4),
                    "reference stopped at its iteration limit"))
    entries.append((Line2D([0], [0], color="black", marker="o", markerfacecolor="none", linestyle="none", markersize=OPEN_MARKER_PT),
                    "unconverged point"))
    fig.text((Y_LABEL_AREA_PT + MARGIN_PT) * style.PT / fig.get_size_inches()[0],
             (height - LEGEND_TITLE_PT + 3.0) * style.PT / fig.get_size_inches()[1],
             f"{case['section']}: {name}  (N = {case.get('n_nodes', '?')}, "
             + ("inviscid)" if case.get("inviscid") else f"ITER {case.get('max_iterations', '?')})"),
             fontsize=style.AXIS_LABEL_PT, family="monospace")
    style.legend_below(fig, Y_LABEL_AREA_PT + MARGIN_PT, height - LEGEND_TITLE_PT,
                       style.TEXT_WIDTH_PT - Y_LABEL_AREA_PT - 2 * MARGIN_PT, entries, 2,
                       row_pt=LEGEND_ROW_PT, first_alone=True)
    style.save(fig, run_dir / f"twins-{name}.svg")


def main():
    style.apply()
    run_dir = Path(sys.argv[1])
    if (run_dir / "twins.json").exists():
        j = json.loads((run_dir / "twins.json").read_text())
        for case in j["cases"]:
            twins_sheet(run_dir, case)
        return
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
