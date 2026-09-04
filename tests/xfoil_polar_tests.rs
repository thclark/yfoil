//! S10 gate: the polar state machine against the reference.
//!
//! The tracked polar case runs XFOIL's OPER script `ALFA 0 / ASEQ 1 5 1 / INIT / ALFA -1 /
//! ASEQ -2 -5 -1` on YFoil's panels; `viscal_points.dat` records every VISCAL call (ITMAX
//! passed, LBLINI on entry, iterations, LVCONV, CL/CM/CD/CDF/CDP/XOCTR/IST/ITRAN). The same
//! sequence is replayed through one persistent `Session`, and every record must match:
//! identical iteration counts, converged flags, IST and ITRAN; converged forces within TOL_SOLVER;
//! per-iteration transients (RLX/RMSBL/CL/CD) within TOL_TRANSIENT (see tolerances.rs).

mod fixtures;
mod utilities;

use std::collections::HashMap;
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_SOLVER, TOL_TRANSIENT};
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::solver::analysis::{compute_polar, FlowSpec, OperatingPoint, PolarConfig, Session};

const CASE: &str = "tests/fixtures/xfoil/naca0012_n60_polar_re1e6";

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{CASE}/{name}"))
}

fn parse_points(path: &PathBuf) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out: Vec<HashMap<String, String>> = Vec::new();
    for l in text.lines() {
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

fn spec() -> FlowSpec {
    FlowSpec {
        re: 1.0e6,
        mach: 0.0,
        ncrit: 9.0,
        itmax: 20,
        ..FlowSpec::default()
    }
}

/// `IT k i RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2` lines of viscal_iters_all.dat, by call.
fn parse_iters_all(path: &PathBuf) -> HashMap<usize, Vec<Vec<f64>>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out: HashMap<usize, Vec<Vec<f64>>> = HashMap::new();
    for l in text.lines() {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(|t| t.parse().unwrap()).collect();
        out.entry(v[0] as usize).or_default().push(v[1..].to_vec());
    }
    out
}

fn check_point(k: usize, p: &OperatingPoint, x: &HashMap<String, String>, entry_lblini: bool, niter1: usize) {
    let ctx = format!("call {k} (alpha {:.1}°)", p.alpha.to_degrees());
    // per-iteration first: the first diverging iteration is the diagnostic
    let all = parse_iters_all(&fixture_path("viscal_iters_all.dat"));
    let xi = &all[&k];
    for (y, xv) in p.trace.iter().zip(xi) {
        let ictx = format!("{ctx} iteration {}", y.iter);
        println!(
            "  {ictx}: yfoil rms {:.6e} rlx {:.4} CL {:.10} IST {} | xfoil rms {:.6e} rlx {:.4} CL {:.10} ISTB {} IST {}",
            y.rmsbl, y.rlx, y.cl, y.ist, xv[1], xv[2], xv[3], xv[6] as usize, xv[7] as usize
        );
        assert_within(y.rmsbl, xv[1], TOL_TRANSIENT, 1.0, &format!("{ictx}: RMSBL"));
        assert_within(y.rlx, xv[2], TOL_TRANSIENT, 1.0, &format!("{ictx}: RLX"));
        assert_within(y.cl, xv[3], TOL_TRANSIENT, 1.0, &format!("{ictx}: CL"));
        assert_within(y.cd, xv[4], TOL_TRANSIENT, 1.0, &format!("{ictx}: CD"));
        assert_eq!(y.ist, xv[7] as usize, "{ictx}: IST after STMOVE");
        assert_eq!(y.itran[1..], [xv[8] as usize, xv[9] as usize], "{ictx}: ITRAN");
    }
    assert_eq!(
        niter1,
        x["NITER1"].parse::<usize>().unwrap(),
        "{ctx}: ITMAX passed to VISCAL"
    );
    assert_eq!(entry_lblini, x["LBLINI0"] == "T", "{ctx}: LBLINI on entry");
    assert_within(
        p.alpha,
        x["ALFA"].parse().unwrap(),
        TOL_SOLVER,
        1.0,
        &format!("{ctx}: ALFA"),
    );
    assert_eq!(
        p.iterations,
        x["NITDONE"].parse::<usize>().unwrap(),
        "{ctx}: iteration count"
    );
    assert_eq!(p.converged, x["LVCONV"] == "T", "{ctx}: LVCONV");
    for (name, ours) in [
        ("RMSBL", p.rmsbl),
        ("CL", p.cl),
        ("CM", p.cm),
        ("CD", p.cd),
        ("CDF", p.cdf),
        ("CDP", p.cdp),
        ("CL_ALF", p.cl_alf),
        ("XOCTR1", p.xtr_upper),
        ("XOCTR2", p.xtr_lower),
    ] {
        assert_within(
            ours,
            x[name].parse().unwrap(),
            TOL_SOLVER,
            1.0,
            &format!("{ctx}: {name}"),
        );
    }
    assert_eq!(
        p.itran[1..],
        [
            x["ITRAN1"].parse::<usize>().unwrap(),
            x["ITRAN2"].parse::<usize>().unwrap()
        ],
        "{ctx}: ITRAN"
    );
    println!(
        "{ctx}: {:2} iterations, converged={} CL {:.8} CD {:.8} CM {:+.8} XTR {:.5}/{:.5} — match",
        p.iterations, p.converged, p.cl, p.cd, p.cm, p.xtr_upper, p.xtr_lower
    );
}

#[test]
fn test_polar_sequence_matches_xfoil_point_by_point() {
    let geometry = read_geometry_from_file(fixture_path("panels.json").to_str().unwrap()).expect("panels.json");
    let airfoil = create_paneled_airfoil(&geometry);
    let xf = parse_points(&fixture_path("viscal_points.dat"));
    assert_eq!(xf.len(), 11, "expected 11 VISCAL calls (0..5, INIT, -1..-5)");

    let mut session = Session::new(&airfoil, spec());
    let mut k = 0;
    let run = |session: &mut Session, alpha_deg: f64, aseq: bool, k: &mut usize| {
        let entry_lblini = session.st.lblini;
        let p = if aseq {
            session.aseq(alpha_deg.to_radians())
        } else {
            session.alfa(alpha_deg.to_radians())
        };
        let niter1 = if aseq { 25 } else { 20 };
        check_point(*k + 1, &p, &xf[*k], entry_lblini, niter1);
        *k += 1;
        p
    };

    // ALFA 0 / ASEQ 1 5 1
    let mut seq_points = vec![run(&mut session, 0.0, false, &mut k)];
    for a in 1..=5 {
        seq_points.push(run(&mut session, a as f64, true, &mut k));
    }
    // INIT / ALFA -1 / ASEQ -2 -5 -1
    session.init();
    assert!(!session.st.lblini && !session.st.lipan, "INIT clears LBLINI and LIPAN");
    seq_points.push(run(&mut session, -1.0, false, &mut k));
    for a in 2..=5 {
        seq_points.push(run(&mut session, -(a as f64), true, &mut k));
    }
    assert_eq!(k, 11);

    // compute_polar drives exactly this sequence and stitches ascending
    let polar = compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 5.0,
            alpha_min: -5.0,
            alpha_step: 1.0,
            spec: spec(),
            ..PolarConfig::default()
        },
    );
    assert!(polar.completed);
    let converged: Vec<&OperatingPoint> = seq_points.iter().filter(|p| p.converged).collect();
    assert_eq!(polar.points.len(), converged.len(), "compute_polar point count");
    for p in &polar.points {
        let q = converged
            .iter()
            .find(|q| (q.alpha - p.alpha).abs() < 1e-12)
            .expect("compute_polar alpha present in the manual sequence");
        assert_eq!(
            p.cl.to_bits(),
            q.cl.to_bits(),
            "compute_polar CL at {:.1}°",
            p.alpha.to_degrees()
        );
        assert_eq!(
            p.cd.to_bits(),
            q.cd.to_bits(),
            "compute_polar CD at {:.1}°",
            p.alpha.to_degrees()
        );
        assert_eq!(p.iterations, q.iterations);
    }
    assert!(
        polar.points.windows(2).all(|w| w[0].alpha < w[1].alpha),
        "points ascend in alpha"
    );
}
