//! The whole-run comparison used by the (b) run tests: yFoil driven through the case's own OPER
//! script — every `ALFA`, `CL`, `ASEQ` point and `INIT`, in one session — and required to take
//! XFOIL's route through the whole run.
//!
//! Per VISCAL call (`viscal_points.dat`, `viscal_iters_all.dat`): the iteration count and
//! convergence exactly; per iteration, IST and ITRAN exactly and RMSBL, RLX, CL, CD, CM, ALFA,
//! MINF and REINF within the tolerance; the converged point (ALFA, CL, CM, CD, CDF, CDP, the
//! transition points, MINF, REINF; IST and ITRAN exactly) within `TOL_SOLVER`. The per-iteration
//! tolerance is `TOL_SOLVER` for a run's first call and `TOL_TRANSIENT` for the calls after it,
//! which carry the previous point's state. On a host other than the fixture's both become
//! `TOL_CROSS_HOST` at least (`common/utilities/host.rs`). There is no third outcome: a run whose
//! route differs from XFOIL's anywhere is a divergent comparison and not a run test
//! (`docs/conventions/testing.md`, rule 3) — such cases are tested by single steps. An inviscid
//! case's points are compared one by one (`step.rs`, `replay_inviscid`).

use crate::fixtures::cases::{load, Op};
use crate::step::replay_inviscid;
use crate::utilities::fortran::real;
use crate::utilities::host::same_host;
use crate::utilities::tolerances::{assert_within, TOL_CROSS_HOST, TOL_SOLVER, TOL_TRANSIENT};
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{PointResult, Session};

/// `IT k i RMSBL RLX CL CD CM ISTB IST ITRAN1 ITRAN2 ALFA MINF REINF` rows, by call.
fn iterations(dir: &std::path::Path) -> Vec<Vec<Vec<f64>>> {
    let mut out: Vec<Vec<Vec<f64>>> = vec![];
    for l in std::fs::read_to_string(dir.join("viscal_iters_all.dat"))
        .unwrap()
        .lines()
    {
        let Some(r) = l.strip_prefix("IT ") else { continue };
        let v: Vec<f64> = r.split_whitespace().map(real).collect();
        let k = v[0] as usize;
        while out.len() < k {
            out.push(vec![]);
        }
        out[k - 1].push(v[1..].to_vec());
    }
    out
}

/// Drive `case` through its script and require XFOIL's route and values through VISCAL call
/// `through` (0: every call); returns the points compared.
pub fn compare_run(case: &str, through: usize) -> Vec<PointResult> {
    let c = load(case);
    if c.spec.re.is_none() {
        // an inviscid run: every point is one SPECAL or SPECCL, compared with XFOIL's record
        let last = if through == 0 { c.ops.len() } else { through };
        for call in 1..=last {
            replay_inviscid(case, call);
        }
        return vec![];
    }
    let host_tol = |t: f64| if same_host(&c.dir) { t } else { t.max(TOL_CROSS_HOST) };
    let its = iterations(&c.dir);
    assert_eq!(
        c.points.len(),
        c.ops.len(),
        "{case}: one VISCAL call per operating point"
    );
    let g = read_geometry_from_file(c.dir.join("panels.json").to_str().unwrap()).expect("panels.json");
    let mut session = Session::new(&panel_foil(&g), c.spec.clone());
    let mut out = vec![];
    let last = if through == 0 { c.ops.len() } else { through };
    for (n, op) in c.ops.iter().enumerate().take(last) {
        let k = n + 1;
        if c.inits.contains(&n) {
            session.init();
        }
        let p = match *op {
            Op::Alfa(a) => session.alpha(a.to_radians()),
            Op::Cl(cl) => session.cl(cl),
            Op::Aseq(a) => session.sequence_point(a.to_radians()),
        };
        let x = &c.points[n];
        let ctx = format!("{case} call {k}");
        assert_eq!(
            p.iterations,
            x["NITDONE"].parse::<usize>().unwrap(),
            "{ctx}: iteration count"
        );
        assert_eq!(p.converged, x["LVCONV"] == "T", "{ctx}: convergence");
        let tol = host_tol(if k == 1 { TOL_SOLVER } else { TOL_TRANSIENT });
        for (y, r) in p.iteration_records.iter().zip(&its[n]) {
            let ictx = format!("{ctx} iteration {}", y.iteration);
            assert_eq!(y.i_stagnation_node, r[7] as usize, "{ictx}: IST");
            assert_eq!(
                y.i_transition_station[1..],
                [r[8] as usize, r[9] as usize],
                "{ictx}: ITRAN"
            );
            for (name, ours, theirs, scale) in [
                ("RMSBL", y.residual, r[1], 1.0),
                ("RLX", y.relaxation, r[2], 1.0),
                ("CL", y.cl, r[3], 1.0),
                ("CD", y.cd, r[4], 1.0),
                ("CM", y.cm, r[5], 1.0),
                ("ALFA", y.alpha, r[10], 1.0),
                ("MINF", y.mach, r[11], 1.0),
                ("REINF", y.re, r[12], r[12]),
            ] {
                assert_within(ours, theirs, tol, scale, &format!("{ictx}: {name}"));
            }
        }
        let st = session.state();
        let pt = host_tol(TOL_SOLVER);
        for (name, ours, scale) in [
            ("ALFA", p.alpha, 1.0),
            ("CL", p.cl, 1.0),
            ("CM", p.cm, 1.0),
            ("CD", p.cd, 1.0),
            ("CDF", p.cd_friction, 1.0),
            ("CDP", p.cd_pressure, 1.0),
            ("XOCTR1", p.transition_upper[0], 1.0),
            ("XOCTR2", p.transition_lower[0], 1.0),
            ("MINF", st.mach, 1.0),
            ("REINF", st.re, st.re),
        ] {
            assert_within(ours, real(&x[name]), pt, scale, &format!("{ctx}: {name}"));
        }
        assert_eq!(st.i_stagnation_node, x["IST"].parse::<usize>().unwrap(), "{ctx}: IST");
        assert_eq!(
            p.i_transition_station[1..],
            [
                x["ITRAN1"].parse::<usize>().unwrap(),
                x["ITRAN2"].parse::<usize>().unwrap()
            ],
            "{ctx}: ITRAN"
        );
        out.push(p);
    }
    out
}
