#!/usr/bin/env python3
"""One-off: translate tracked fixture *inputs* to the names of docs/conventions/naming.md.
Textual key renames only, so formatting (and hence `cargo xtask fixtures --verify`) is preserved.
XFOIL's own dumps under tests/fixtures/xfoil/ carry Fortran names by design and are untouched."""
import glob, re
def sub_keys(path, pairs):
    s = open(path).read(); t = s
    for a, b in pairs:
        t = re.sub(r'"%s"(\s*[:=])' % re.escape(a), r'"%s"\1' % b, t)
    if t != s: open(path, "w").write(t); print("translated", path)
geom = [("reference", "cm_ref"), ("x_c", "x"), ("y_c", "y")]
for p in glob.glob("tests/fixtures/xfoil/*/panels.json") + ["tests/fixtures/aerofoil_with_invalid_last_point.json"]:
    sub_keys(p, geom)
for p in glob.glob("tests/fixtures/*/*/input.json"):
    sub_keys(p, [("n_crit", "ncrit")])
for p in glob.glob("tests/fixtures/*/*/bl_stations.json"):
    sub_keys(p, [("delta_star", "dstar")])
wake = [("THET_TE1", "theta_te_station1"), ("THET_TE2", "theta_te_station2"), ("DSTR_TE1", "dstar_te_station1"),
        ("DSTR_TE2", "dstar_te_station2"), ("CTAU_TE1", "sqrtctau_te_station1"), ("CTAU_TE2", "sqrtctau_te_station2"),
        ("UEDG_TE1", "ue_te_station1"), ("UEDG_TE2", "ue_te_station2"), ("ANTE", "te_thickness_normal"),
        ("TTE", "theta_te"), ("DTE", "dstar_te"), ("CTE", "sqrtctau_te"), ("IW", "i_wake"), ("XSSI", "xi"),
        ("UEDG", "ue"), ("THET", "theta"), ("DSTR", "dstar"), ("MASS", "mass_defect"), ("CTAU", "sqrtctau"),
        ("T2", "theta_station2"), ("D2", "dstar_station2"), ("U2", "ue_station2"), ("HK2", "hk_station2"),
        ("HS2", "hstar_station2"), ("CF2", "cf_station2"), ("DI2", "cdiss_station2"), ("airfoil", "foil"), ("reynolds", "re")]
for p in glob.glob("tests/fixtures/*/wake_*.json"):
    sub_keys(p, wake)
for p in glob.glob("tests/fixtures/xfoil/*/manifest.json"):
    sub_keys(p, [("airfoil", "foil"), ("iter", "max_iterations")])
s = open("xtask/fixtures-config/cases.toml").read()
s = re.sub(r"^airfoil = ", "foil = ", s, flags=re.M); s = re.sub(r"^iter = ", "max_iterations = ", s, flags=re.M)
open("xtask/fixtures-config/cases.toml", "w").write(s); print("translated cases.toml")
