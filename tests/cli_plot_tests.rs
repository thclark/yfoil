//! `yfoil plot foil` end to end through the binary: geometry, analysis and polar inputs, the
//! SVG and PNG renderers, and the input-dispatch errors. Needs the `plotting` feature.
#![cfg(feature = "plotting")]

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

struct Case {
    dir: tempfile::TempDir,
}

impl Case {
    fn new() -> Self {
        let dir = tempfile::tempdir().unwrap();
        let c = Self { dir };
        ok(&["geometry", "naca", "0012", "-n", "60", "-o", c.p("n0012.json")]);
        ok(&["geometry", "naca", "0012", "-n", "40", "-o", c.p("n0012_40.json")]);
        ok(&[
            "analyze",
            c.p("n0012.json"),
            "--alpha",
            "4",
            "-r",
            "1e6",
            "-o",
            c.p("a4.json"),
        ]);
        ok(&[
            "analyze",
            c.p("n0012.json"),
            "--alpha",
            "4",
            "-r",
            "1e6",
            "-m",
            "0.3",
            "-o",
            c.p("a4m3.json"),
        ]);
        c
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.path().join(name)
    }

    fn p(&self, name: &str) -> &'static str {
        Box::leak(self.path(name).to_string_lossy().into_owned().into_boxed_str())
    }

    fn svg(&self, name: &str) -> String {
        std::fs::read_to_string(self.path(name)).unwrap()
    }
}

fn count(hay: &str, needle: &str) -> usize {
    hay.matches(needle).count()
}

#[test]
fn geometry_inputs_overlay_panels_only() {
    let c = Case::new();
    ok(&["plot", "foil", c.p("n0012.json"), "-o", c.p("one.svg")]);
    let svg = c.svg("one.svg");
    // outline + 60 notches (default for geometry input), no legend
    assert_eq!(count(&svg, "<polyline"), 1);
    let notch_group = svg
        .split("<g stroke=\"rgb(64,64,64)\" stroke-width=\"1\">")
        .nth(1)
        .expect("notch group");
    assert_eq!(count(notch_group.split("</g>").next().unwrap(), "<line"), 60);
    assert!(!svg.contains("<rect x=\"70"), "no legend for a single geometry");

    ok(&[
        "plot",
        "foil",
        c.p("n0012.json"),
        c.p("n0012_40.json"),
        "--panels",
        "dots",
        "-o",
        c.p("two.svg"),
    ]);
    let svg = c.svg("two.svg");
    assert_eq!(count(&svg, "<polyline"), 2);
    assert_eq!(count(&svg, "<circle"), 100);
    assert!(svg.contains("n0012_40"), "legend names each geometry");

    fails(
        &["plot", "foil", c.p("n0012.json"), "--quantity", "dstar"],
        "geometry input has no boundary layer",
    );
    fails(
        &["plot", "foil", c.p("n0012.json"), "--wake"],
        "geometry input has no wake",
    );
    fails(&["plot", "foil", c.p("n0012.json"), "--alpha", "4"], "--alpha");
}

#[test]
fn analysis_inputs_draw_quantities_wake_and_markers() {
    let c = Case::new();
    ok(&[
        "plot",
        "foil",
        c.p("a4.json"),
        "--panels",
        "notches",
        "--wake",
        "--quantity",
        "dstar,theta",
        "-o",
        c.p("one.svg"),
    ]);
    let svg = c.svg("one.svg");
    // surface + wake line, then per quantity: wake upper, wake lower, upper, lower
    assert_eq!(count(&svg, "<polyline"), 2 + 2 * 4);
    assert!(svg.contains("δ*, scaled ×"));
    assert!(svg.contains("θ, scaled ×"));
    // stagnation circle, and a marker legend naming what is drawn (this case has a laminar
    // separation bubble on the lower surface)
    assert_eq!(count(&svg, "fill=\"none\" stroke=\"rgb(40,40,40)\""), 1);
    for text in [
        "stagnation",
        "transition",
        "separation",
        "reattachment",
        "upper surface",
        "lower surface",
    ] {
        assert!(svg.contains(&format!(">{text}<")), "{text} in the marker legend");
    }

    // two design points × two quantities: four legend entries, four distinct styles
    ok(&[
        "plot",
        "foil",
        c.p("a4.json"),
        c.p("a4m3.json"),
        "--quantity",
        "dstar,theta",
        "--label",
        "M=0",
        "--label",
        "M=0.3",
        "--no-markers",
        "-o",
        c.p("grid.svg"),
    ]);
    let svg = c.svg("grid.svg");
    assert!(svg.contains("δ* at M=0, scaled ×"));
    assert!(svg.contains("θ at M=0.3, scaled ×"));
    assert_eq!(
        count(&svg, "stroke-dasharray=\"2,3\""),
        4 + 2,
        "second point dotted: 4 curves + 2 legend keys"
    );
    assert_eq!(count(&svg, "fill=\"none\" stroke=\"rgb(40,40,40)\""), 0, "no markers");
    assert!(!svg.contains(">stagnation<"), "no marker legend");

    // PNG through the same SVG
    ok(&[
        "plot",
        "foil",
        c.p("a4.json"),
        "--quantity",
        "dstar",
        "-o",
        c.p("one.png"),
    ]);
    let png = std::fs::read(c.path("one.png")).unwrap();
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    ok(&[
        "plot",
        "foil",
        c.p("a4.json"),
        "--format",
        "png",
        "-o",
        c.p("forced.out"),
    ]);
    assert_eq!(&std::fs::read(c.path("forced.out")).unwrap()[..4], b"\x89PNG");

    // mismatched panels cannot be overlaid with quantities
    ok(&[
        "analyze",
        c.p("n0012_40.json"),
        "--alpha",
        "4",
        "-r",
        "1e6",
        "-o",
        c.p("a4_40.json"),
    ]);
    fails(
        &["plot", "foil", c.p("a4.json"), c.p("a4_40.json"), "--quantity", "dstar"],
        "different panel geometry",
    );
    fails(&["plot", "foil", c.p("a4.json"), "--alpha", "4"], "needs a polar JSON");
    fails(&["plot", "foil", c.p("a4.json"), c.p("n0012.json")], "one kind");
    fails(
        &["plot", "foil", c.p("a4.json"), "--quantity", "nope"],
        "unknown BL quantity",
    );
}

#[test]
fn polar_distributions_select_design_points() {
    let c = Case::new();
    ok(&[
        "polar",
        c.p("n0012.json"),
        "--alpha-min",
        "-2",
        "--alpha-max",
        "4",
        "--alpha-step",
        "2",
        "-r",
        "1e6",
        "--distributions",
        "-o",
        c.p("polar.json"),
    ]);
    let text = std::fs::read_to_string(c.path("polar.json")).unwrap();
    assert!(text.contains("\"distributions\""));

    ok(&[
        "plot",
        "foil",
        c.p("polar.json"),
        "--alpha",
        "0,4",
        "--quantity",
        "dstar",
        "--wake",
        "-o",
        c.p("sweep.svg"),
    ]);
    let svg = c.svg("sweep.svg");
    assert!(svg.contains("δ* at α=0°, scaled ×"));
    assert!(svg.contains("δ* at α=4°, scaled ×"));
    assert!(svg.contains(">polar<"), "the polar stem is the title");
    assert!(!svg.contains("α=2°"));
    // one scale factor for both points
    let ks: Vec<&str> = svg
        .split(", scaled ×")
        .skip(1)
        .map(|rest| rest.split('<').next().unwrap().trim())
        .collect();
    assert_eq!(ks.len(), 2);
    assert_eq!(ks[0], ks[1]);

    // all embedded points by default (4 alphas: -2, 0, 2, 4)
    ok(&[
        "plot",
        "foil",
        c.p("polar.json"),
        "--quantity",
        "hk",
        "-o",
        c.p("all.svg"),
    ]);
    assert_eq!(count(&c.svg("all.svg"), "Hk at α="), 4);

    fails(
        &["plot", "foil", c.p("polar.json"), "--alpha", "3"],
        "no distribution at alpha = 3",
    );
    ok(&[
        "polar",
        c.p("n0012.json"),
        "--alpha-min",
        "0",
        "--alpha-max",
        "1",
        "--alpha-step",
        "1",
        "-r",
        "1e6",
        "-o",
        c.p("plain.json"),
    ]);
    fails(&["plot", "foil", c.p("plain.json")], "carries no distributions");
}

#[test]
fn help_lists_plot_foil() {
    let stdout = ok(&["plot", "foil", "--help"]);
    assert!(stdout.contains("--quantity"));
    assert!(stdout.contains("--max-offset"));
    let stdout = ok(&["polar", "--help"]);
    assert!(stdout.contains("--distributions"));
    assert!(Path::new(env!("CARGO_BIN_EXE_yfoil")).exists());
}
