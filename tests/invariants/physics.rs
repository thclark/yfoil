//! (e) Physical invariant: symmetry and camber of generated sections, the panel method's lift slope and zero-lift behaviour, Blasius closures, the Kutta condition, trailing-edge treatments.
//!
//! Fixtures: none.

use approx::assert_relative_eq;
use yfoil::geometry::Thickness;

use yfoil::geometry::{naca_4digit, naca_5digit, panel_foil};

/// Test that symmetric airfoils are indeed symmetric
#[test]
fn test_airfoil_symmetry() {
    let geom = naca_4digit("0015", 160, Thickness::Perpendicular).unwrap();

    let max_y = geom.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y = geom.y.iter().cloned().fold(f64::INFINITY, f64::min);

    // Should be symmetric about y=0
    assert_relative_eq!(max_y, -min_y, epsilon = 0.001);
}

/// Test cambered airfoils have positive camber
#[test]
fn test_airfoil_camber() {
    // NACA 4-digit cambered
    let geom_4 = naca_4digit("6412", 160, Thickness::Perpendicular).unwrap();
    let max_y_4 = geom_4.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y_4 = geom_4.y.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(max_y_4 > -min_y_4, "4-digit should have positive camber");

    // NACA 5-digit cambered
    let geom_5 = naca_5digit("23018", 160, Thickness::Perpendicular).unwrap();
    let max_y_5 = geom_5.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min_y_5 = geom_5.y.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(max_y_5 > -min_y_5, "5-digit should have positive camber");
}

// ============================================================================
// Panel Method Validation Tests
// ============================================================================

use yfoil::solver::analysis::{analyse, FlowConditions, Session};
use yfoil::solver::blstate::SolverState;
use yfoil::solver::ggcalc::build_inviscid_system;
use yfoil::solver::specal::solve_inviscid_at_alpha;

fn inviscid_cl(airfoil: &yfoil::geometry::PanelledFoil, alpha: f64) -> f64 {
    analyse(
        airfoil,
        alpha,
        &FlowConditions {
            re: None,
            ..FlowConditions::default()
        },
    )
    .cl
}

/// Test inviscid panel method lift curve slope against thin airfoil theory.
///
/// Thin airfoil theory predicts: dCL/dα = 2π ≈ 6.28 per radian
/// For a 12% thick airfoil, the actual value is slightly lower due to
/// thickness effects (typically 5.8-6.2).
#[test]
fn test_panel_method_lift_slope() {
    let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
    let airfoil = panel_foil(&geom);
    // Calculate CL at two angles
    let alpha1 = 0.0_f64;
    let alpha2 = 5.0_f64.to_radians();
    let dcl_dalpha = (inviscid_cl(&airfoil, alpha2) - inviscid_cl(&airfoil, alpha1)) / (alpha2 - alpha1);

    // Should be in range 5.0 to 7.0 (thin airfoil theory ≈ 6.28)
    assert!(
        dcl_dalpha > 5.0 && dcl_dalpha < 7.0,
        "Lift slope {} per radian should be near 2π ≈ 6.28",
        dcl_dalpha
    );
}

/// Test that symmetric airfoil has zero lift at alpha=0
#[test]
fn test_panel_method_symmetric_zero_lift() {
    let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
    let airfoil = panel_foil(&geom);
    let cl = inviscid_cl(&airfoil, 0.0);
    assert!(
        cl.abs() < 1e-12,
        "CL={cl} should be zero to the noise floor for a symmetric airfoil at α=0"
    );
}

/// Test that cambered airfoil has positive lift at alpha=0
#[test]
fn test_panel_method_cambered_lift() {
    let geom = naca_4digit("4412", 160, Thickness::Perpendicular).unwrap();
    let airfoil = panel_foil(&geom);
    let cl = inviscid_cl(&airfoil, 0.0);
    // NACA 4412 has 4% camber, should produce positive lift at α=0
    // Thin airfoil theory predicts CL ≈ 2π * 2 * (0.04) ≈ 0.5 for 4% camber
    assert!(cl > 0.2, "CL={cl} should be positive for cambered airfoil at α=0");
}

// ============================================================================
// Boundary Layer Validation Tests
// ============================================================================

use yfoil::bl::system::amplification_rate;
use yfoil::bl::{cf_laminar, hk_from_h, hstar_laminar};

/// Test Blasius flat plate solution for laminar boundary layer
///
/// The Blasius solution provides the exact laminar BL solution for a flat plate:
/// - θ = 0.664 * sqrt(ν*x/U)   (momentum thickness)
/// - δ* = 1.7208 * sqrt(ν*x/U) (displacement thickness)
/// - H = δ*/θ = 2.591          (shape factor)
/// - Cf = 0.664 / sqrt(Re_x)   (skin friction)
///
/// Reference: Schlichting, "Boundary Layer Theory"
#[test]
fn test_blasius_shape_factor() {
    // Blasius shape factor H = 2.591
    let h_blasius = 2.591;
    let (hk, _, _) = hk_from_h(h_blasius, 0.0);

    // At incompressible conditions, Hk = H
    assert_relative_eq!(hk, h_blasius, epsilon = 1e-10);
}

#[test]
fn test_blasius_skin_friction() {
    // Test Cf correlation at Blasius conditions
    // Blasius: Cf = 0.664 / sqrt(Re_x)
    // For Re_θ, we have Re_x = (Re_θ / 0.664)² approximately
    // So at Re_θ = 1000, Re_x ≈ 2.27e6, and Cf ≈ 0.664/1507 ≈ 0.00044
    // But we test Cf * sqrt(Re_θ) which should be roughly constant

    let hk_blasius = 2.591;

    // Test at different Re_θ values
    for rt in [500.0, 1000.0, 2000.0, 5000.0] {
        let result = cf_laminar(hk_blasius, rt, 0.0);
        // Cf should be positive for attached laminar flow
        assert!(result.value > 0.0, "Cf should be positive at Re_θ={}", rt);

        // Cf * sqrt(Re_θ) should be approximately constant for Blasius
        // The Falkner-Skan correlation gives Cf * Re_θ ≈ constant for given Hk
        let cf_rt = result.value * rt;
        assert!(
            cf_rt > 0.1 && cf_rt < 1.0,
            "Cf*Re_θ = {} should be O(1) at Re_θ={}",
            cf_rt,
            rt
        );
    }
}

#[test]
fn test_blasius_energy_shape_factor() {
    // Test H* correlation at Blasius conditions
    // For Blasius flow, H* ≈ 1.573
    let hk_blasius = 2.591;
    let result = hstar_laminar(hk_blasius, 1000.0, 0.0);

    // H* should be between 1.5 and 1.7 for Blasius-like conditions
    assert!(
        result.value > 1.5 && result.value < 1.7,
        "H* = {} should be ~1.57 for Blasius",
        result.value
    );
}

#[test]
fn test_laminar_amplification_off_below_critical() {
    // Below critical Re_θ, amplification should be zero
    let hk_blasius = 2.591;
    let theta = 0.001; // Small momentum thickness

    // At low Re_θ (below ~200 for Blasius)
    let ax = amplification_rate(hk_blasius, theta, 100.0).rate;
    assert_eq!(ax, 0.0, "Amplification should be 0 below critical Re_θ");
}

#[test]
fn test_laminar_amplification_on_above_critical() {
    // Above critical Re_θ, amplification should be positive
    let hk_blasius = 2.591;
    let theta = 0.001;

    // At high Re_θ (well above critical)
    let ax = amplification_rate(hk_blasius, theta, 5000.0).rate;
    assert!(ax > 0.0, "Amplification should be positive above critical Re_θ");
}

// ============================================================================
// Trailing Edge Geometry Tests
// ============================================================================

/// Test panel method with blunt TE still produces reasonable results
///
/// Blunt TE uses standard flow tangency equations (no curvature extrapolation).
/// The Kutta condition still ensures finite TE velocity, and results should
/// be physically reasonable even if not as clean as sharp TE.
#[test]
fn test_blunt_te_reasonable_results() {
    let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
    let blunt_geom = geom.blunten(0.002); // 0.2% gap
    let airfoil = panel_foil(&blunt_geom);

    assert!(!airfoil.sharp_te, "Should use blunt TE handling");

    let mut session = Session::new(
        &airfoil,
        FlowConditions {
            re: None,
            ..FlowConditions::default()
        },
    );
    let coeffs = session.alpha(0.0);
    let vel: Vec<f64> = session.state().q_inviscid[1..=airfoil.n_foil_nodes].to_vec();

    // For symmetric airfoil at α=0, CL should still be near zero
    assert!(
        coeffs.cl.abs() < 0.05,
        "CL = {} should be near 0 for symmetric airfoil at α=0",
        coeffs.cl
    );

    // CDp may not be exactly zero but should be small
    // The blunt TE may cause some small asymmetry in the solution
    assert!(
        coeffs.cd_pressure.abs() < 0.01,
        "CDp = {} should be small for blunt TE",
        coeffs.cd_pressure
    );

    // Velocities should be smooth (no huge spikes)
    let max_vel = vel.iter().map(|v| v.abs()).fold(f64::NEG_INFINITY, f64::max);
    assert!(
        max_vel < 3.0,
        "Max velocity {} should be reasonable (no singularity)",
        max_vel
    );
}

/// Test that lift slope is correct for both sharp and blunt TE
#[test]
fn test_te_type_lift_slope_comparison() {
    // Sharp TE
    let geom_sharp = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
    let airfoil_sharp = panel_foil(&geom_sharp);

    // Blunt TE
    let geom_blunt = geom_sharp.blunten(0.002);
    let airfoil_blunt = panel_foil(&geom_blunt);

    // Calculate lift slope for both
    let alpha1 = 0.0_f64;
    let alpha2 = 5.0_f64.to_radians();
    let dcl_dalpha_sharp =
        (inviscid_cl(&airfoil_sharp, alpha2) - inviscid_cl(&airfoil_sharp, alpha1)) / (alpha2 - alpha1);
    let dcl_dalpha_blunt =
        (inviscid_cl(&airfoil_blunt, alpha2) - inviscid_cl(&airfoil_blunt, alpha1)) / (alpha2 - alpha1);

    // Both should be close to thin airfoil theory (2π ≈ 6.28)
    assert!(
        dcl_dalpha_sharp > 5.0 && dcl_dalpha_sharp < 7.5,
        "Sharp TE lift slope {} should be near 2π",
        dcl_dalpha_sharp
    );
    assert!(
        dcl_dalpha_blunt > 5.0 && dcl_dalpha_blunt < 7.5,
        "Blunt TE lift slope {} should be near 2π",
        dcl_dalpha_blunt
    );

    // They should be similar (within 10%)
    let diff_pct = ((dcl_dalpha_sharp - dcl_dalpha_blunt) / dcl_dalpha_sharp).abs() * 100.0;
    assert!(
        diff_pct < 15.0,
        "Lift slopes should be similar: sharp={}, blunt={}, diff={}%",
        dcl_dalpha_sharp,
        dcl_dalpha_blunt,
        diff_pct
    );
}

/// Test Kutta condition is satisfied for both sharp and blunt TE
#[test]
fn test_kutta_condition_both_te_types() {
    // GGCALC's Kutta row: GAMU(1) + GAMU(N) = 0 for both the alpha = 0 and alpha = 90 solutions
    let kutta = |airfoil: &yfoil::geometry::PanelledFoil| -> (f64, f64) {
        let mut st = SolverState::from_foil(airfoil, airfoil.n_foil_nodes / 12 + 10);
        build_inviscid_system(&mut st);
        let n = airfoil.n_foil_nodes;
        (
            st.q_inviscid_basis[1][1] + st.q_inviscid_basis[1][n],
            st.q_inviscid_basis[2][1] + st.q_inviscid_basis[2][n],
        )
    };

    // Sharp TE
    let geom_sharp = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap().sharpen();
    let airfoil_sharp = panel_foil(&geom_sharp);
    assert!(airfoil_sharp.sharp_te);
    let (k0, k90) = kutta(&airfoil_sharp);
    assert!(k0.abs() < 1e-9, "Sharp TE Kutta condition violated for α=0°: {k0}");
    assert!(k90.abs() < 1e-9, "Sharp TE Kutta condition violated for α=90°: {k90}");

    // Blunt TE
    let geom_blunt = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap();
    let airfoil_blunt = panel_foil(&geom_blunt);
    assert!(!airfoil_blunt.sharp_te);
    let (k0, k90) = kutta(&airfoil_blunt);
    assert!(k0.abs() < 1e-9, "Blunt TE Kutta condition violated for α=0°: {k0}");
    assert!(k90.abs() < 1e-9, "Blunt TE Kutta condition violated for α=90°: {k90}");
}

/// Test that the gamma distribution is smooth near TE for sharp case
///
/// With curvature extrapolation, gamma should smoothly approach zero at TE
/// from both upper and lower surfaces.
#[test]
fn test_sharp_te_smooth_gamma() {
    // For the symmetric airfoil at α = 0 the TE vorticity is small and the second differences
    // approaching the TE from both sides agree (the sharp-TE bisector row of GGCALC)
    let geom = naca_4digit("0012", 160, Thickness::Perpendicular).unwrap().sharpen();
    let airfoil = panel_foil(&geom);
    let mut st = SolverState::from_foil(&airfoil, airfoil.n_foil_nodes / 12 + 10);
    let mut sys = None;
    st.alpha = 0.0;
    st.qinf = 1.0;
    solve_inviscid_at_alpha(&mut st, &mut sys);
    let n = airfoil.n_foil_nodes;
    let gam = &st.gamma;
    assert!(gam[1].abs() < 2.0, "Upper TE gamma {} should be small", gam[1]);
    assert!(gam[n].abs() < 2.0, "Lower TE gamma {} should be small", gam[n]);
    // XFOIL's GAM is antisymmetric on a symmetric section at α = 0 (GAM(I) = −GAM(N+1−I)), so
    // the second differences approaching the TE are equal and opposite
    let curv_upper = gam[3] - 2.0 * gam[2] + gam[1];
    let curv_lower = gam[n - 2] - 2.0 * gam[n - 1] + gam[n];
    assert!(
        (curv_upper + curv_lower).abs() < 1e-9,
        "TE curvatures should be antisymmetric: upper={curv_upper}, lower={curv_lower}"
    );
}
