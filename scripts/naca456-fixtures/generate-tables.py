#!/usr/bin/env python3
"""Write src/geometry/series/six_series_tables.rs from naca456's epspsi.f90.

Usage: generate-tables.py <path to epspsi.f90> > src/geometry/series/six_series_tables.rs
"""
import re
import sys

src = open(sys.argv[1]).read()


def table(name):
    m = re.search(r"DIMENSION\(201\)::\s*%s\s*=\s*\(/(.*?)/\)" % name, src, re.S)
    nums = [float(v) for v in re.findall(r"-?\d+\.\d+", m.group(1))]
    assert len(nums) == 201, (name, len(nums))
    return nums


out = [
    "//! The ε and ψ functions of the eight NACA 6- and 6A-series thickness families, 201 values each at",
    "//! φ = kπ/200, k = 0…200, as tabulated in the NASA program of TM-4741 (Ladson, Brooks, Hill and",
    "//! Sproles, 1996) and carried by PDAS `naca456` (`epspsi.f90`, Carmichael 2001). The conformal",
    "//! mapping in `six_series.rs` turns each pair into the basic thickness form; the scale factor that",
    "//! sets t/c is a polynomial fit (`scale_factor`). These numbers are the definition of the families:",
    "//! TM-4741 records that the original functions were plotted on large sheets of graph paper and read",
    "//! off for interpolation, and that the graphs are lost.",
    "//!",
    "//! Generated from `epspsi.f90` (naca456.zip, sha256 327f9dae…) by",
    "//! `scripts/naca456-fixtures/generate-tables.py`; do not edit.",
    "",
    "/// Number of φ stations of every table",
    "pub const N_PHI: usize = 201;",
    "",
]
families = [("63", "1"), ("64", "2"), ("65", "3"), ("66", "4"), ("67", "5"), ("63A", "6"), ("64A", "7"), ("65A", "8")]
for fam, idx in families:
    for kind in ("EPS", "PSI"):
        vals = table(kind + idx)
        sym = "EPSILON" if kind == "EPS" else "PSI"
        out.append(f"/// {'ε' if kind == 'EPS' else 'ψ'}(φ) of the {fam} family (`{kind}{idx}` in `epspsi.f90`)")
        out.append(f"pub const {sym}_{fam}: [f64; N_PHI] = [")
        for i in range(0, 201, 8):
            out.append("    " + ", ".join(f"{v:.5f}" for v in vals[i : i + 8]) + ",")
        out.append("];")
        out.append("")
sys.stdout.write("\n".join(out))
