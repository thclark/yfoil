//! (c) yFoil functionality: the polar driver: `compute_polar` runs XFOIL's polar procedure (0 → max, `INIT`, re-solve
//! 0°, 0 → min) through one session, keeps the re-solved 0° out of the points, and stitches the
//! legs ascending; the re-solved 0° is bit-identical to the first.
//!
//! Fixtures: the panels of `tests/fixtures/xfoil/naca0012_n60_polar_re1e6/` — `cargo xtask fixtures --case naca0012_n60_polar_re1e6`.

use crate::fixtures;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::{compute_polar, FlowConditions, PointResult, PolarConfig, Session};

fn spec() -> FlowConditions {
    FlowConditions {
        re: Some(1.0e6),
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    }
}

#[test]
fn test_polar_driver_reproduces_the_manual_sequence() {
    let path = fixtures::require_fixture("tests/fixtures/xfoil/naca0012_n60_polar_re1e6/panels.json");
    let geometry = read_geometry_from_file(path.to_str().unwrap()).expect("panels.json");
    let airfoil = panel_foil(&geometry);

    // the manual sequence: ALFA 0 / ASEQ 1 5 1 / INIT / ALFA 0 / ASEQ -1 -5 -1
    let mut session = Session::new(&airfoil, spec());
    let mut seq_points = vec![session.alpha(0.0)];
    for a in 1..=5 {
        seq_points.push(session.sequence_point((a as f64).to_radians()));
    }
    session.init();
    assert!(
        !session.state().bl_initialised && !session.state().pointers_built,
        "INIT clears LBLINI and LIPAN"
    );
    let reseed = session.alpha(0.0);
    for a in 1..=5 {
        seq_points.push(session.sequence_point(-(a as f64).to_radians()));
    }

    // the re-solved 0° is the first 0° again, bit for bit
    let first = &seq_points[0];
    assert_eq!(
        reseed.iterations, first.iterations,
        "yFoil re-solved 0°: iteration count"
    );
    assert_eq!(reseed.converged, first.converged);
    for (name, a, b) in [
        ("CL", reseed.cl, first.cl),
        ("CD", reseed.cd, first.cd),
        ("CM", reseed.cm, first.cm),
        ("CDF", reseed.cd_friction, first.cd_friction),
        ("CDP", reseed.cd_pressure, first.cd_pressure),
        ("CL_ALF", reseed.cl_d_alpha, first.cl_d_alpha),
        ("RMSBL", reseed.residual, first.residual),
        ("XOCTR1", reseed.transition_upper[0], first.transition_upper[0]),
        ("XOCTR2", reseed.transition_lower[0], first.transition_lower[0]),
    ] {
        assert_eq!(a.to_bits(), b.to_bits(), "yFoil re-solved 0° after INIT: {name}");
    }

    // compute_polar drives exactly this sequence, keeps the 0° seed out of the points, and
    // stitches ascending
    let polar = compute_polar(
        &airfoil,
        &PolarConfig {
            alpha_max: 5.0,
            alpha_min: -5.0,
            alpha_step: 1.0,
            conditions: spec(),
            ..PolarConfig::default()
        },
    );
    assert!(polar.completed);
    let converged: Vec<&PointResult> = seq_points.iter().filter(|p| p.converged).collect();
    assert_eq!(polar.results.len(), converged.len(), "compute_polar point count");
    for p in &polar.results {
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
        polar.results.windows(2).all(|w| w[0].alpha < w[1].alpha),
        "points ascend in alpha"
    );
}
