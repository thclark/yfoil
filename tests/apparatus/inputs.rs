//! (f) Reference integrity: every tracked reference case's inputs agree with each other.
//!
//! A case directory under `tests/fixtures/xfoil/` carries the inputs the reference was run from
//! beside the files its tests read: `manifest.json` (the case definition and the reference build),
//! `xfoil.inp` (the OPER script) and, for a case run on yFoil's panels, `panels.json` and the
//! `panels.dat` XFOIL `LOAD`s. For every case: the manifest names its own directory, the script
//! disables graphics before anything else, a `LOAD`ing script loads `panels.dat`, and that file is
//! byte for byte what yFoil's `.dat` writer makes of `panels.json` (CLAUDE.md Rule 4: yFoil
//! generates, XFOIL consumes; the 17-significant-digit rendering round-trips every coordinate).
//! The handoff into XFOIL's own arrays is asserted when the fixture is generated
//! (`cargo xtask fixtures`, against XFOIL's dump of X/Y after `LOAD`).
//!
//! Fixtures: every `tests/fixtures/xfoil/<case>/` — `cargo xtask fixtures`.

use std::path::Path;
use yfoil::geometry::{read_geometry_from_file, write_dat_file};

fn check(dir: &Path) {
    let case = dir.file_name().unwrap().to_string_lossy().to_string();
    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).expect("manifest.json")).unwrap();
    assert_eq!(
        manifest["case"]["name"],
        case.as_str(),
        "{case}: manifest names its directory"
    );
    let script = std::fs::read_to_string(dir.join("xfoil.inp")).expect("xfoil.inp");
    let lines: Vec<&str> = script.lines().map(str::trim).filter(|l| !l.is_empty()).collect();
    assert_eq!(&lines[..2], ["PLOP", "G F"], "{case}: graphics off first");
    if !lines.contains(&"LOAD panels.dat") {
        assert!(
            !dir.join("panels.dat").exists(),
            "{case}: panels.dat is not loaded by the script"
        );
        return;
    }
    let geometry = read_geometry_from_file(dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let written = tempfile::NamedTempFile::new().unwrap();
    write_dat_file(&geometry, &case, written.path()).unwrap();
    let tracked = std::fs::read_to_string(dir.join("panels.dat")).expect("panels.dat");
    let ours = std::fs::read_to_string(written.path()).unwrap();
    // the first line is the name XFOIL prints; the coordinates are the handoff
    assert_eq!(
        tracked.lines().skip(1).collect::<Vec<_>>(),
        ours.lines().skip(1).collect::<Vec<_>>(),
        "{case}: panels.dat is yFoil's rendering of panels.json"
    );
}

#[test]
fn every_tracked_case_has_consistent_inputs() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/xfoil");
    let mut n = 0;
    for entry in std::fs::read_dir(&root).unwrap().flatten() {
        if entry.path().is_dir() {
            check(&entry.path());
            n += 1;
        }
    }
    assert!(n > 0, "no tracked cases");
}
