//! A tracked reference case as the (a) and (b) tests drive it: its flow conditions and OPER
//! script from `manifest.json` and `xfoil.inp`, its VISCAL calls' records (`viscal_points.dat`),
//! and the seeding of a yFoil session with XFOIL's exact state on entry to one SETBL call — the
//! inviscid side built by yFoil from the panels and the call's operating point (SPECAL or SPECCL,
//! then VISCAL's set-up with no iterations: wake, QINV, pointers, UINV, DIJ), the BL side and the
//! pointers from `mrchdu_input_<k>.dat`.

use super::mrchdu_fixtures::{parse_bl_dump, BlDump};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use yfoil::bl::system::{AmplificationModel, MachClDependence, ReClDependence};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{FlowConditions, Session};
use yfoil::solver::clcalc::set_compressibility;
use yfoil::solver::pointers::{map_stations_to_nodes, map_stations_to_rows, set_station_xi};
use yfoil::solver::setbl::set_mach_re_from_cl;
use yfoil::solver::specal::{alpha_command, cl_command, sequence_command};
use yfoil::solver::velocity::set_ue_inviscid;
use yfoil::solver::viscal::solve_viscous;

/// `|a − b| ≤ tol · max(|a|, |b|, scale)` — the Rule 1 metric (`utilities::tolerances`).
fn assert_within(a: f64, b: f64, tol: f64, scale: f64, what: &str) {
    let err = (a - b).abs() / a.abs().max(b.abs()).max(scale);
    assert!(
        err <= tol,
        "{what}: yfoil={a:.17e} xfoil={b:.17e} err={err:.3e} > tol={tol:.1e}"
    );
}

/// One operating point of the reference script, as yFoil drives it.
#[derive(Clone, Copy, Debug)]
pub enum Op {
    Alfa(f64),
    Cl(f64),
    /// a point of an `ASEQ`, alpha in degrees as `compute_polar` forms it
    Aseq(f64),
}

pub struct Case {
    pub name: String,
    pub dir: PathBuf,
    pub spec: FlowConditions,
    /// the operating points XFOIL ran, one per VISCAL (or SPECAL/SPECCL) call, in order — the
    /// script's points less those an `ASEQ` never reached because it halted (`executed`)
    pub ops: Vec<Op>,
    /// the number of calls before each `INIT` of the script
    pub inits: Vec<usize>,
    /// per VISCAL call: NITDONE, LVCONV and the converged point's values
    pub points: Vec<HashMap<String, String>>,
}

pub fn ops(script: &str) -> Vec<Op> {
    let mut out = vec![];
    for l in script.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        match t.first().copied() {
            Some("ALFA") => out.push(Op::Alfa(t[1].parse().unwrap())),
            Some("CL") => out.push(Op::Cl(t[1].parse().unwrap())),
            Some("ASEQ") => {
                let v: Vec<f64> = t[1..4].iter().map(|x| x.parse().unwrap()).collect();
                let np = ((v[1] - v[0]) / v[2] + 0.5).floor() as usize + 1;
                out.extend((0..np).map(|i| Op::Aseq(v[0] + v[2] * i as f64)));
            }
            _ => {}
        }
    }
    out
}

/// The points XFOIL ran and the number of calls before each `INIT`. An `ASEQ` halts after NSEQEX
/// = 4 consecutive points that do not converge (`xoper.f`, 'Sequence halted'), skipping the rest of
/// that sequence; which points converged is read from the calls' `LVCONV` in `points`. Without
/// VISCAL records (an inviscid run) every point runs.
fn executed(script: &str, points: &[HashMap<String, String>]) -> (Vec<Op>, Vec<usize>) {
    const NSEQEX: usize = 4;
    let mut out = vec![];
    let mut inits = vec![];
    for l in script.lines() {
        let t: Vec<&str> = l.split_whitespace().collect();
        if l.trim() == "INIT" {
            inits.push(out.len());
            continue;
        }
        let seq: Vec<Op> = match t.first().copied() {
            Some("ASEQ") => ops(l),
            Some("ALFA") | Some("CL") => {
                out.extend(ops(l));
                continue;
            }
            _ => continue,
        };
        let mut failures = 0;
        for op in seq {
            let n = out.len();
            out.push(op);
            let converged = points.get(n).is_none_or(|p| p.get("LVCONV").is_none_or(|v| v == "T"));
            failures = if converged { 0 } else { failures + 1 };
            if !points.is_empty() && failures >= NSEQEX {
                break;
            }
        }
    }
    (out, inits)
}

pub fn blocks(path: &Path) -> Vec<HashMap<String, String>> {
    let mut out: Vec<HashMap<String, String>> = vec![];
    for l in std::fs::read_to_string(path).unwrap().lines() {
        let Some((k, v)) = l.split_once('=') else { continue };
        if k.trim() == "CALL" {
            out.push(HashMap::new());
        }
        if let Some(cur) = out.last_mut() {
            cur.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

pub fn load(name: &str) -> Case {
    // the route measurement (`cargo xtask route`) reads its own reference runs, every SETBL call
    // dumped, from the directory it names
    #[cfg(yfoil_route)]
    let dir = match std::env::var("YFOIL_ROUTE_FIXTURES") {
        Ok(root) => PathBuf::from(root).join(name),
        Err(_) => super::require_fixture(&format!("tests/fixtures/xfoil/{name}")),
    };
    #[cfg(not(yfoil_route))]
    let dir = super::require_fixture(&format!("tests/fixtures/xfoil/{name}"));
    let m: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("manifest.json")).unwrap()).unwrap();
    let c = &m["case"];
    let matyp = c["matyp"].as_u64().unwrap_or(0);
    let xtr: Vec<f64> = c["xtr"]
        .as_array()
        .map(|a| a.iter().map(|v| v.as_f64().unwrap()).collect())
        .unwrap_or_default();
    let spec = FlowConditions {
        re: (!c["inviscid"].as_bool().unwrap_or(false)).then(|| c["re"].as_f64().unwrap()),
        mach: c["mach"].as_f64().unwrap(),
        ncrit: c["ncrit"].as_f64().unwrap(),
        max_iterations: c["max_iterations"].as_u64().unwrap() as usize,
        x_trip: if xtr.is_empty() { [1.0, 1.0] } else { [xtr[0], xtr[1]] },
        amplification_model: if c["damp"].as_bool().unwrap_or(false) {
            AmplificationModel::ModifiedEnvelope
        } else {
            AmplificationModel::Envelope
        },
        mach_cl_dependence: if matyp == 2 {
            MachClDependence::InverseSqrtCl
        } else {
            MachClDependence::Fixed
        },
        re_cl_dependence: match matyp {
            2 => ReClDependence::InverseSqrtCl,
            3 => ReClDependence::InverseCl,
            _ => ReClDependence::Fixed,
        },
        ..FlowConditions::default()
    };
    let script = std::fs::read_to_string(dir.join("xfoil.inp")).unwrap();
    let points = if dir.join("viscal_points.dat").exists() {
        blocks(&dir.join("viscal_points.dat"))
    } else {
        vec![]
    };
    let (ops, inits) = executed(&script, &points);
    Case {
        name: name.to_string(),
        spec,
        ops,
        inits,
        points,
        dir,
    }
}

pub fn dump(c: &Case, file: &str) -> BlDump {
    let path = c.dir.join(file);
    let rel = path
        .strip_prefix(env!("CARGO_MANIFEST_DIR"))
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|_| path.to_string_lossy().to_string());
    parse_bl_dump(&super::require_fixture(&rel))
}

/// A session at the prologue of VISCAL call `call`: the call's operating point (SPECAL or SPECCL
/// at the alpha or CL it was given), then VISCAL's set-up with no iterations (wake, QINV,
/// pointers, UINV, DIJ). `alpha` overrides an `ALFA`/`ASEQ` point's angle with XFOIL's own.
pub fn prologue(c: &Case, call: usize, alpha: Option<f64>) -> Session {
    prologue_jogged(c, call, alpha, None)
}

/// The one-ULP jog of item `i` under `seed`: −1, 0 or +1 ULP for each of a pair of values, drawn
/// from splitmix64 of (seed, i) — the jog `cargo xtask twins` applies to panel coordinates.
pub fn jog(seed: u64, i: usize) -> (i8, i8) {
    let mut z = (seed << 32) ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    let three = |v: u64| -> i8 { (v % 3) as i8 - 1 };
    (three(z & 0xFFFF_FFFF), three(z >> 32))
}

/// `v` moved by `by` ULP (−1, 0 or +1).
pub fn ulp(v: f64, by: i8) -> f64 {
    match by {
        1 => v.next_up(),
        -1 => v.next_down(),
        _ => v,
    }
}

/// A dump with every BL entry but XSSI (which the pointers rebuild) moved by −1, 0 or +1 ULP.
pub fn jog_dump(d: &BlDump, seed: u64) -> BlDump {
    let mut o = d.clone();
    let mut n = 0;
    for is in 1..=2 {
        for row in o.bl[is].iter_mut() {
            for m in (1..6).step_by(2) {
                let (a, b) = jog(seed, n);
                n += 1;
                row[m] = ulp(row[m], a);
                if m + 1 < 6 {
                    row[m + 1] = ulp(row[m + 1], b);
                }
            }
        }
    }
    o
}

/// `prologue` on the case's panels, with every coordinate jogged by −1, 0 or +1 ULP under `seed`
/// when one is given.
pub fn prologue_jogged(c: &Case, call: usize, alpha: Option<f64>, seed: Option<u64>) -> Session {
    let mut g = read_geometry_from_file(c.dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    if let Some(seed) = seed {
        for i in 0..g.x.len() {
            let (a, b) = jog(seed, i);
            g.x[i] = ulp(g.x[i], a);
            g.y[i] = ulp(g.y[i], b);
        }
    }
    let mut session = Session::new(&panel_foil(&g), c.spec.clone());
    let wake_length = c.spec.wake_length;
    let (st, sys) = session.parts_mut();
    match c.ops[call - 1] {
        Op::Alfa(a) => alpha_command(st, sys, alpha.unwrap_or(a.to_radians())),
        Op::Aseq(a) => sequence_command(st, sys, alpha.unwrap_or(a.to_radians())),
        Op::Cl(cl) => {
            cl_command(st, sys, cl);
        }
    }
    solve_viscous(st, sys.as_mut(), 0, wake_length, None);
    session
}

/// Overwrite the session's BL side with XFOIL's state entering SETBL call `k`.
pub fn seed(session: &mut Session, d: &BlDump, tol: f64) {
    let st = session.state_mut();
    st.bl_initialised = true;
    st.converged = false;
    st.alpha = d.real("ALFA");
    st.alpha_specified = d.logical("LALFA");
    st.i_stagnation_node = d.int("IST");
    st.s_stagnation = d.real("SST");
    st.s_stagnation_d_gamma_node0 = d.real("SST_GO");
    st.s_stagnation_d_gamma_node1 = d.real("SST_GP");
    map_stations_to_nodes(st);
    set_station_xi(st);
    map_stations_to_rows(st);
    assert_eq!(
        [d.int("NBL1"), d.int("NBL2")],
        [st.n_stations[1], st.n_stations[2]],
        "NBL from the dumped IST"
    );
    st.i_transition_station = [0, d.int("ITRAN1"), d.int("ITRAN2")];
    st.cl = d.real("CLMR");
    let (m_cl, re_cl) = set_mach_re_from_cl(st, st.cl);
    st.mach_d_cl = m_cl;
    st.re_d_cl = re_cl;
    set_compressibility(st);
    set_ue_inviscid(st);
    for is in 1..=2 {
        for ibl in 1..=st.n_stations[is] {
            let r = d.bl[is][ibl];
            assert_within(
                st.xi[is][ibl],
                r[0],
                tol,
                1.0,
                &format!("XSSI({ibl},{is}) rebuilt from the dumped IST/SST"),
            );
            st.xi[is][ibl] = r[0];
            st.ue[is][ibl] = r[1];
            st.theta[is][ibl] = r[2];
            st.dstar[is][ibl] = r[3];
            st.sqrtctau[is][ibl] = r[4];
            st.mass_defect[is][ibl] = r[5];
        }
    }
}
