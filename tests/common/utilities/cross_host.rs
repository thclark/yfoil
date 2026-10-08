//! The cross-host gate of a whole run: on a host other than the fixture's, a value cannot be met
//! more closely than the reference itself moves between the two hosts' maths libraries, so each
//! run-test value is gated at `max(tolerance, TOL_CROSS_HOST, CROSS_HOST_FACTOR × spread)`, the
//! spread measured per case and value by `cargo xtask cross-host` (`scripts/cross-host.sh`) into
//! `tests/fixtures/cross-host/spread.json` (CLAUDE.md Rule 1).
//!
//! On every host the measurement must exist for the case and have been taken against the
//! fixture's current records (their length and FNV-1a hash are compared), so a regenerated run case
//! fails until it is remeasured. On the fixture's own host the tolerance then stands. On another
//! host the measurement must also cover it — the same libm as the measurement's other host — or
//! the test fails and says what to measure. A case whose reference takes a different route on the
//! two hosts cannot be a cross-host run test, and says so.

use super::host::Host;
use super::tolerances::{CROSS_HOST_FACTOR, TOL_CROSS_HOST};
use std::path::Path;
use std::sync::OnceLock;

fn spread_file() -> &'static serde_json::Value {
    static DOC: OnceLock<serde_json::Value> = OnceLock::new();
    DOC.get_or_init(|| {
        let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cross-host/spread.json");
        let text =
            std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e} — run scripts/cross-host.sh", p.display()));
        serde_json::from_str(&text).expect("spread.json")
    })
}

/// FNV-1a (64-bit), as `xtask/src/cross_host.rs` computes it.
fn fnv64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The gate for one value of a run test and, when the measured spread sets it, a note saying so
/// for the failure message.
pub struct Gate {
    pub tol: f64,
    pub note: String,
}

/// The gate of `value` (`level` = `iteration`, `point` or `node`) of case `case` in `dir`, whose
/// tolerance on the fixture's own host is `tol`.
pub fn gate(case: &str, dir: &Path, level: &str, value: &str, tol: f64) -> Gate {
    let doc = spread_file();
    let measured = doc["measured_host"]["host"].as_str().unwrap_or("");
    // on every host: the measurement exists for the case and was taken against its current records
    let c = &doc["cases"][case];
    assert!(
        c.is_object(),
        "{case}: no cross-host measurement — run scripts/cross-host.sh"
    );
    for (file, rec) in c["records"].as_object().into_iter().flatten() {
        let bytes = std::fs::read(dir.join(file)).unwrap_or_else(|e| panic!("{case}: {file}: {e}"));
        let now = format!("{:016x}", fnv64(&bytes));
        assert!(
            rec["bytes"].as_u64() == Some(bytes.len() as u64) && rec["fnv64"].as_str() == Some(now.as_str()),
            "{case}: {file} has changed since its cross-host spread was measured — rerun scripts/cross-host.sh"
        );
    }
    if Host::of_fixture(dir) == Host::running() {
        return Gate {
            tol,
            note: String::new(),
        };
    }
    let running = Host::running();
    let os = measured.rsplit('-').next().unwrap_or("").to_ascii_lowercase();
    assert!(
        running.os == os,
        "{case}: running on {running}, but the cross-host spread was measured on {measured}; \
         measure this host's libm with `cargo xtask cross-host` on it"
    );
    if let Some(r) = c["route_differs"].as_str() {
        panic!("{case}: the reference takes a different route on {measured} ({r}); it cannot be a cross-host run test");
    }
    let s = &c["spread"][level][value];
    let spread = s
        .as_f64()
        .unwrap_or_else(|| panic!("{case}: {level} {value}: cross-host spread {s} (non-finite on one host only)"));
    let floor = CROSS_HOST_FACTOR * spread;
    let tol_x = tol.max(TOL_CROSS_HOST);
    if floor > tol_x {
        Gate {
            tol: floor,
            note: format!(
                " [cross-host gate {CROSS_HOST_FACTOR} × the reference's own {} ↔ {measured} spread {spread:.2e}]",
                doc["fixture_host"]["host"].as_str().unwrap_or("")
            ),
        }
    } else {
        Gate {
            tol: tol_x,
            note: String::new(),
        }
    }
}
