//! Which host a fixture was generated on, and whether that is this one.
//!
//! Both codes take `EXP`/`LOG`/`**`/`ATAN2`/`SIN`/`COS`/`TANH` from the host C library (gfortran's
//! XFOIL through libgfortran, yFoil through `std`'s `f64` methods), and the libraries differ:
//! measured 2026-09-11 on 20 000 inputs per function, Apple libSystem 1345.120.2 and glibc 2.39
//! disagree by exactly 1 ULP on 0.1 % (`exp`, `ln`, `powf`) to 18 % (`tanh`) of inputs, on
//! `x86_64` and `aarch64` alike (glibc's results are identical on both). So bit-identity with a
//! fixture is a property of the host it was generated on (CLAUDE.md Rule 1); on any other host
//! the comparison is an ULP budget, and the pins that only a bit-identical trajectory can hold
//! (the iteration at which a threshold-straddling run parts, a one-step replay at `TOL_SOLVER`
//! in a hypersensitive state) are reported instead of asserted.

#![allow(dead_code)]

use std::path::Path;

/// A host as the fixture manifests record it (`uname -m`-`uname -s`-…): the pair that decides
/// which libm the transcendentals came from, normalised so that the manifest's `arm64-Darwin`
/// and Rust's `aarch64`/`macos` compare equal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Host {
    pub arch: String,
    pub os: String,
}

impl Host {
    fn normalise(arch: &str, os: &str) -> Host {
        let arch = match arch {
            "arm64" => "aarch64",
            "amd64" => "x86_64",
            a => a,
        }
        .to_string();
        let os = match os.to_ascii_lowercase().as_str() {
            "darwin" => "macos".to_string(),
            o => o.to_string(),
        };
        Host { arch, os }
    }

    /// The host these tests are running on.
    pub fn running() -> Host {
        Host::normalise(std::env::consts::ARCH, std::env::consts::OS)
    }

    /// The host that generated the fixture in `dir`, from its `manifest.json` (`xfoil_ref` →
    /// `host: <arch>-<os>-<release>`).
    pub fn of_fixture(dir: &Path) -> Host {
        let m: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).expect("manifest.json"))
                .expect("manifest.json");
        let line = m["xfoil_ref"]
            .as_array()
            .and_then(|a| a.iter().filter_map(|v| v.as_str()).find(|s| s.starts_with("host:")))
            .unwrap_or_else(|| panic!("{}: manifest.json has no `host:` line", dir.display()));
        let triple = line["host:".len()..].trim();
        let mut parts = triple.splitn(3, '-');
        let (arch, os) = (parts.next().unwrap_or(""), parts.next().unwrap_or(""));
        Host::normalise(arch, os)
    }

    /// The `libm:` line of the fixture's manifest, for reports.
    pub fn libm_of_fixture(dir: &Path) -> String {
        let m: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).expect("manifest.json"))
                .expect("manifest.json");
        m["xfoil_ref"]
            .as_array()
            .and_then(|a| a.iter().filter_map(|v| v.as_str()).find(|s| s.starts_with("libm:")))
            .map(|s| s["libm:".len()..].trim().to_string())
            .unwrap_or_default()
    }
}

impl std::fmt::Display for Host {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.arch, self.os)
    }
}

/// Whether the fixture in `dir` was generated on this host. When it was not, prints why the
/// same-host pins are being reported rather than asserted, once per call.
pub fn same_host(dir: &Path) -> bool {
    let (fixture, running) = (Host::of_fixture(dir), Host::running());
    if fixture != running {
        println!(
            "CROSS-HOST: fixture generated on {fixture} ({}), running on {running}: transcendentals come from a different libm (1-ULP differences), so the same-host pins are reported, not asserted",
            Host::libm_of_fixture(dir)
        );
    }
    fixture == running
}
