//! `cargo xtask cross-host`: the reference is the same on another host, bit for bit.
//!
//! Both codes take their transcendentals from CORE-MATH's correctly rounded functions, not from the
//! host's libm (CLAUDE.md Rule 1, `scripts/xfoil-build/core_math_libm.c`), so a fixture generated
//! on one host is reproduced exactly on any other and the tests hold every value to the same
//! tolerance everywhere. This command checks that claim. Run on a host other than the fixtures'
//! (`scripts/cross-host.sh` does so in a Linux container), it reruns the instrumented reference on
//! each tracked run case's own `panels.dat` and `xfoil.inp` and requires every file the case tracks
//! that the run writes again — the records the run tests compare, and every dump kept — to be
//! byte-identical to the tracked one. It exits non-zero, naming the files, if any differs.

use crate::{run_xfoil, Cases, RunEnd};
use std::fs;
use std::process::Command;

fn host() -> String {
    let u = |a: &str| {
        Command::new("uname")
            .arg(a)
            .output()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .unwrap_or_default()
    };
    format!("{}-{}", u("-m"), u("-s"))
}

pub(crate) fn run(_flags: &[String]) {
    let root = crate::root();
    let cases: Cases = toml::from_str(
        &fs::read_to_string(root.join("xtask/fixtures-config/cases.toml")).expect("xtask/fixtures-config/cases.toml"),
    )
    .expect("parse cases.toml");
    let xfoil = crate::ensure_reference("instrumented");
    let here = host();
    let mut differing: Vec<String> = vec![];
    let mut compared = 0;
    for case in cases.cases.iter().filter(|c| c.track && c.run) {
        let fixture = root.join("tests/fixtures/xfoil").join(&case.name);
        let work = root.join("target/cross-host").join(&case.name);
        let _ = fs::remove_dir_all(&work);
        fs::create_dir_all(&work).unwrap();
        for f in ["panels.dat", "xfoil.inp"] {
            fs::copy(fixture.join(f), work.join(f)).unwrap_or_else(|e| panic!("{}: {f}: {e}", case.name));
        }
        print!("{} ... ", case.name);
        match run_xfoil(&xfoil, &work) {
            Ok(RunEnd::Finished) | Ok(RunEnd::Hung) => {}
            other => panic!("{}: the reference did not run to its end: {other:?}", case.name),
        }
        // every tracked file the run wrote again; the inputs and the manifest (which names the host)
        // are not outputs of the run
        let mut files: Vec<String> = fs::read_dir(&fixture)
            .unwrap()
            .flatten()
            .map(|e| e.file_name().to_string_lossy().to_string())
            .filter(|n| !matches!(n.as_str(), "manifest.json" | "panels.dat" | "panels.json" | "xfoil.inp"))
            .filter(|n| work.join(n).exists())
            .collect();
        files.sort();
        let differ: Vec<&String> = files
            .iter()
            .filter(|f| fs::read(fixture.join(f)).ok() != fs::read(work.join(f)).ok())
            .collect();
        compared += files.len();
        if differ.is_empty() {
            println!("{} files byte-identical", files.len());
        } else {
            println!("DIFFERENT: {differ:?}");
            differing.extend(differ.iter().map(|f| format!("{}/{f}", case.name)));
        }
    }
    if !differing.is_empty() {
        eprintln!(
            "on {here} the reference differs from the tracked fixtures in {} of {compared} files: {differing:?}",
            differing.len()
        );
        std::process::exit(1);
    }
    println!("on {here} the reference reproduces all {compared} tracked records byte for byte");
}
