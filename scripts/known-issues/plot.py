#!/usr/bin/env python3
"""Draw the known-issues figures from a run folder's summary.json. The run folder has the series
cases' shape (the two studies share their Rust library), so their figures are drawn by the series
cases' plot.py. Run via scripts/figures/render.sh known-issues <run_dir>.
"""
import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parents[1] / "series-cases" / "plot.py"), run_name="__main__")
