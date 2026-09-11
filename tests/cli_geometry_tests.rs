//! `yfoil geometry` end to end through the binary: the help every subcommand prints when run
//! bare, `repanel` with XFOIL's PANGEN (PPAR parameters by flag and by file) and yFoil's cosine
//! method, the generators' panelling options, the provenance record, and every conflict the
//! panelling flags and files must reject.

use std::path::{Path, PathBuf};
use std::process::Command;

fn yfoil() -> Command {
    Command::new(env!("CARGO_BIN_EXE_yfoil"))
}

fn run(args: &[&str]) -> (bool, String, String) {
    let out = yfoil().args(args).output().expect("run yfoil");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

fn ok(args: &[&str]) -> String {
    let (success, stdout, stderr) = run(args);
    assert!(success, "yfoil {:?} failed:\n{stdout}\n{stderr}", args);
    stdout
}

fn fails(args: &[&str], needle: &str) {
    let (success, stdout, stderr) = run(args);
    assert!(!success, "yfoil {:?} unexpectedly succeeded:\n{stdout}", args);
    assert!(
        stderr.contains(needle),
        "yfoil {:?}: stderr lacks '{needle}':\n{stderr}",
        args
    );
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

struct Case {
    dir: tempfile::TempDir,
}

impl Case {
    fn new() -> Self {
        Self {
            dir: tempfile::tempdir().unwrap(),
        }
    }
    fn path(&self, name: &str) -> PathBuf {
        self.dir.path().join(name)
    }
    fn p(&self, name: &str) -> String {
        self.path(name).to_string_lossy().to_string()
    }
    /// A 60-node NACA 0012 `.dat`, the "arbitrary loaded geometry" of the tests
    fn source_dat(&self) -> String {
        ok(&[
            "geometry",
            "naca",
            "0012",
            "-n",
            "60",
            "--method",
            "cosine",
            "--to",
            "dat",
            "-o",
            &self.p("src.dat"),
        ]);
        self.p("src.dat")
    }
}

#[test]
fn every_geometry_subcommand_prints_its_help_when_run_bare() {
    for sub in ["convert", "naca", "karman-trefftz", "repanel", "info"] {
        let (success, stdout, stderr) = run(&["geometry", sub]);
        let text = format!("{stdout}{stderr}");
        assert!(!success, "{sub}: a bare call is not a run");
        assert!(
            text.contains("Usage:") && text.contains("Options:"),
            "{sub}: no help printed:\n{text}"
        );
    }
    let (_, stdout, stderr) = run(&["geometry", "repanel"]);
    let text = format!("{stdout}{stderr}");
    for flag in [
        "--method",
        "--panelling",
        "--curvature-bunching",
        "--te-curvature-ratio",
        "--refined-curvature-ratio",
        "--refine-upper",
        "--refine-lower",
        "--cosine-te-bias",
        "--sharp",
        "--te-gap",
        "--te-blend",
    ] {
        assert!(text.contains(flag), "repanel help lacks {flag}");
    }
    for heading in ["PANGEN parameters", "Cosine parameters", "Trailing edge"] {
        assert!(text.contains(heading), "repanel help lacks the heading {heading}");
    }
    // clap prints the bare-command help on stderr with exit code 2
    let (_, stdout, stderr) = run(&["geometry"]);
    let text = format!("{stdout}{stderr}");
    assert!(text.contains("repanel") && text.contains("karman-trefftz"), "{text}");
}

#[test]
fn repanel_defaults_to_pangen_with_xfoil_defaults_and_records_it() {
    let c = Case::new();
    let src = c.source_dat();
    let out = ok(&["geometry", "repanel", &src]);
    assert!(out.contains("Repanelled from 60 to 160 nodes with pangen"));
    let g = json(&c.path("src_repanelled.json"));
    assert_eq!(g["x"].as_array().unwrap().len(), 160);
    let p = &g["generator"]["panelling"];
    assert_eq!(p["method"], "pangen");
    assert_eq!(p["n_nodes"], 160);
    assert_eq!(p["curvature_bunching"], 1.0);
    assert_eq!(p["te_curvature_ratio"], 0.15);
    assert_eq!(p["refined_curvature_ratio"], 0.2);
    assert!(p["refine_upper"].is_null() && p["refine_lower"].is_null());
    assert_eq!(p["sharp_te"], false);
    assert!(p["te_gap"].is_null());
    // a .dat carries no section record: only yfoil and the panelling
    assert!(g["generator"].get("designation").is_none());
    assert!(g["generator"]["yfoil"].is_string());
}

#[test]
fn ppar_flags_are_recorded_and_a_file_reproduces_them_bitwise() {
    let c = Case::new();
    let src = c.source_dat();
    ok(&[
        "geometry",
        "repanel",
        &src,
        "-n",
        "200",
        "--curvature-bunching",
        "1.5",
        "--te-curvature-ratio",
        "0.3",
        "--refined-curvature-ratio",
        "0.5",
        "--refine-upper",
        "0.2,0.4",
        "--refine-lower",
        "0.3,0.6",
        "-o",
        &c.p("flags.json"),
    ]);
    let g = json(&c.path("flags.json"));
    let p = &g["generator"]["panelling"];
    assert_eq!(p["n_nodes"], 200);
    assert_eq!(p["curvature_bunching"], 1.5);
    assert_eq!(p["te_curvature_ratio"], 0.3);
    assert_eq!(p["refined_curvature_ratio"], 0.5);
    assert_eq!(p["refine_upper"], serde_json::json!([0.2, 0.4]));
    assert_eq!(p["refine_lower"], serde_json::json!([0.3, 0.6]));
    // the record is itself a valid panelling file
    std::fs::write(c.path("panelling.json"), serde_json::to_string(p).unwrap()).unwrap();
    ok(&[
        "geometry",
        "repanel",
        &src,
        "--panelling",
        &c.p("panelling.json"),
        "-o",
        &c.p("file.json"),
    ]);
    let f = json(&c.path("file.json"));
    assert_eq!(f["x"], g["x"]);
    assert_eq!(f["y"], g["y"]);
}

#[test]
fn cosine_method_writes_n_plus_one_nodes_and_records_its_bias() {
    let c = Case::new();
    let src = c.source_dat();
    let out = ok(&[
        "geometry",
        "repanel",
        &src,
        "-n",
        "160",
        "--method",
        "cosine",
        "--cosine-te-bias",
        "0.3",
        "-o",
        &c.p("c.json"),
    ]);
    assert!(out.contains("Repanelled from 60 to 161 nodes with cosine"));
    let p = json(&c.path("c.json"))["generator"]["panelling"].clone();
    assert_eq!(p["method"], "cosine");
    assert_eq!(p["te_bias"], 0.3);
    assert!(p.get("curvature_bunching").is_none());
    // the default bias is recorded when none is given
    ok(&["geometry", "repanel", &src, "--method", "cosine", "-o", &c.p("d.json")]);
    assert_eq!(json(&c.path("d.json"))["generator"]["panelling"]["te_bias"], 0.15);
}

#[test]
fn trailing_edge_treatment_is_applied_after_panelling() {
    let c = Case::new();
    let src = c.source_dat();
    ok(&[
        "geometry",
        "repanel",
        &src,
        "-n",
        "100",
        "--te-gap",
        "0.004",
        "--te-blend",
        "0.8",
        "-o",
        &c.p("g.json"),
    ]);
    let g = json(&c.path("g.json"));
    let (x, y) = (g["x"].as_array().unwrap(), g["y"].as_array().unwrap());
    let d = (x[0].as_f64().unwrap() - x[99].as_f64().unwrap()).hypot(y[0].as_f64().unwrap() - y[99].as_f64().unwrap());
    assert!((d - 0.004).abs() < 1e-12, "gap {d}");
    assert_eq!(
        g["generator"]["panelling"]["te_gap"],
        serde_json::json!({ "gap": 0.004, "blend": 0.8 })
    );
    ok(&[
        "geometry",
        "repanel",
        &src,
        "-n",
        "100",
        "--sharp",
        "-o",
        &c.p("s.json"),
    ]);
    let s = json(&c.path("s.json"));
    assert_eq!(s["x"][0], s["x"][99]);
    assert_eq!(s["y"][0], s["y"][99]);
    assert_eq!(s["generator"]["panelling"]["sharp_te"], true);
}

#[test]
fn panelling_conflicts_are_rejected_by_flags_and_by_file() {
    let c = Case::new();
    let src = c.source_dat();
    let x = c.p("x.json");
    fails(
        &["geometry", "repanel", &src, "--sharp", "--te-gap", "0.002", "-o", &x],
        "cannot be used with",
    );
    fails(
        &["geometry", "repanel", &src, "--te-blend", "0.5", "-o", &x],
        "--te-gap <GAP>",
    );
    fails(
        &["geometry", "repanel", &src, "--cosine-te-bias", "0.3", "-o", &x],
        "belongs to --method cosine",
    );
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--method",
            "cosine",
            "--curvature-bunching",
            "2",
            "-o",
            &x,
        ],
        "belongs to --method pangen",
    );
    fails(
        &["geometry", "repanel", &src, "--refine-upper", "0.4,0.2", "-o", &x],
        "refinement window",
    );
    fails(
        &["geometry", "repanel", &src, "--refine-upper", "0.2", "-o", &x],
        "two x/c values",
    );
    fails(&["geometry", "repanel", &src, "-n", "300", "-o", &x], "too many");
    std::fs::write(c.path("ok.json"), r#"{"method": "pangen", "n_nodes": 100}"#).unwrap();
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--panelling",
            &c.p("ok.json"),
            "-n",
            "100",
            "-o",
            &x,
        ],
        "cannot be used with",
    );
    let file = |name: &str, body: &str| {
        std::fs::write(c.path(name), body).unwrap();
        c.p(name)
    };
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--panelling",
            &file(
                "a.json",
                r#"{"method": "cosine", "n_nodes": 100, "curvature_bunching": 2}"#,
            ),
            "-o",
            &x,
        ],
        "belongs to the pangen method",
    );
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--panelling",
            &file(
                "b.json",
                r#"{"method": "pangen", "n_nodes": 100, "sharp_te": true, "te_gap": {"gap": 0.002}}"#,
            ),
            "-o",
            &x,
        ],
        "contradictory",
    );
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--panelling",
            &file("c.json", r#"{"method": "pangen", "n_nodes": 100, "typo": 1}"#),
            "-o",
            &x,
        ],
        "unknown panelling key typo",
    );
    fails(
        &[
            "geometry",
            "repanel",
            &src,
            "--panelling",
            &file(
                "d.json",
                r#"{"method": "pangen", "n_nodes": 100, "n_buffer_nodes": 246}"#,
            ),
            "-o",
            &x,
        ],
        "n_buffer_nodes applies to generated sections only",
    );
}

#[test]
fn generators_default_to_pangen_and_cosine_is_the_analytic_sampling() {
    let c = Case::new();
    ok(&[
        "geometry",
        "naca",
        "63-415",
        "-n",
        "120",
        "--te-curvature-ratio",
        "0.25",
        "-o",
        &c.p("p.json"),
    ]);
    let g = json(&c.path("p.json"));
    assert_eq!(g["x"].as_array().unwrap().len(), 120);
    let p = &g["generator"]["panelling"];
    assert_eq!(p["method"], "pangen");
    assert_eq!(p["n_buffer_nodes"], 246);
    assert_eq!(p["te_curvature_ratio"], 0.25);
    assert_eq!(g["generator"]["designation"], "NACA 63-415");

    ok(&[
        "geometry",
        "naca",
        "63-415",
        "-n",
        "120",
        "--method",
        "cosine",
        "-o",
        &c.p("c.json"),
    ]);
    let g = json(&c.path("c.json"));
    let lib = yfoil::geometry::Section::from_designation("63-415")
        .unwrap()
        .geometry(120);
    let x: Vec<f64> = g["x"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    assert_eq!(x, lib.x, "cosine is the library's analytic sampling, bitwise");
    assert_eq!(g["generator"]["panelling"]["method"], "cosine");
    assert!(g["generator"]["panelling"]["te_bias"].is_null());

    ok(&[
        "geometry",
        "karman-trefftz",
        "-n",
        "100",
        "--te-angle",
        "10",
        "-o",
        &c.p("kt.json"),
    ]);
    let k = json(&c.path("kt.json"));
    assert_eq!(k["x"].as_array().unwrap().len(), 100);
    assert_eq!(k["generator"]["panelling"]["method"], "pangen");

    fails(
        &[
            "geometry",
            "naca",
            "0012",
            "--method",
            "cosine",
            "--cosine-te-bias",
            "0.3",
            "-o",
            &c.p("x.json"),
        ],
        "has no bias",
    );
}

#[test]
fn vertical_thickness_is_xfoil_model_always_pangen() {
    let c = Case::new();
    ok(&[
        "geometry",
        "naca",
        "4412",
        "-n",
        "160",
        "--thickness",
        "vertical",
        "--te-curvature-ratio",
        "0.3",
        "-o",
        &c.p("v.json"),
    ]);
    let g = json(&c.path("v.json"));
    assert_eq!(g["generator"]["thickness_applied"], "vertical");
    assert_eq!(g["generator"]["buffer"], "NACA4");
    assert_eq!(g["generator"]["panelling"]["n_buffer_nodes"], 245);
    assert_eq!(g["generator"]["panelling"]["te_curvature_ratio"], 0.3);
    fails(
        &[
            "geometry",
            "naca",
            "4412",
            "--method",
            "cosine",
            "--thickness",
            "vertical",
            "-o",
            &c.p("x.json"),
        ],
        "always panelled with PANGEN",
    );
    fails(
        &[
            "geometry",
            "naca",
            "63-415",
            "--thickness",
            "vertical",
            "-o",
            &c.p("x.json"),
        ],
        "4- and 5-digit",
    );
}
