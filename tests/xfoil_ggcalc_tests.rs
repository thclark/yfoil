//! S4 gate: GGCALC (unit vorticity distributions through LUDCMP/BAKSUB) and QDCALC (the full
//! (N+NW)×(N+NW) source influence matrix) from YFoil's own geometry, against the reference
//! run's QINVU and DIJ dumps. The whole VISCAL prologue is exercised: GGCALC → QISET →
//! SPECAL's GAM → XYWAKE → QWCALC → QDCALC.
//!
//! DIJ is gated with the row-scaled metric (plan addendum R3): `|a−b| ≤ tol·max(|a|,|b|,rowmax)`.
//! Element-wise relative error is meaningless on its near-zero entries — the noise-floor
//! measurement shows a 1-ULP geometry perturbation moves some of them by 5e-6 relative.

mod fixtures;
mod utilities;

use fixtures::pointers_fixtures::{parse_dij, parse_pointers, parse_uinv};
use std::path::PathBuf;
use utilities::tolerances::{assert_within, TOL_LINALG, TOL_PURE, TOL_SOLVER};
use yfoil::geometry::{create_paneled_airfoil, read_geometry_from_file};
use yfoil::solver::blstate::BlState;
use yfoil::solver::ggcalc::{ggcalc, InviscidSystem};
use yfoil::solver::qdcalc::qdcalc;
use yfoil::solver::velocity::qiset;
use yfoil::solver::xywake::{qwcalc, xywake};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// GGCALC on YFoil's geometry, then the alpha superposition of SPECAL.
fn prologue_to_specal() -> (BlState, InviscidSystem, fixtures::pointers_fixtures::UinvFixture) {
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    let u = parse_uinv(&fixture_path("xfoil_uinv.dat"), 1);
    let geom = read_geometry_from_file(fixture_path("panels.json")).unwrap();
    let af = create_paneled_airfoil(&geom);
    let mut st = BlState::from_airfoil(&af, f.nw);
    let sys = ggcalc(&mut st);
    st.alfa = u.alfa;
    qiset(&mut st, u.alfa);
    // SPECAL: GAM(I) = COSA*GAMU(I,1) + SINA*GAMU(I,2)  (= QINV on the airfoil)
    for i in 1..=st.n {
        st.gam[i] = st.qinv[i];
    }
    (st, sys, u)
}

#[test]
fn test_ggcalc_matches_xfoil_qinvu() {
    let (st, _, u) = prologue_to_specal();
    let mut worst = 0.0_f64;
    for i in 1..=st.n {
        assert_within(st.qinvu[1][i], u.qinvu1[i], TOL_LINALG, 1.0, &format!("QINVU({i},1)"));
        assert_within(st.qinvu[2][i], u.qinvu2[i], TOL_LINALG, 1.0, &format!("QINVU({i},2)"));
        worst = worst
            .max((st.qinvu[1][i] - u.qinvu1[i]).abs())
            .max((st.qinvu[2][i] - u.qinvu2[i]).abs());
    }
    println!("ggcalc: airfoil QINVU within {TOL_LINALG:.0e} (worst abs diff {worst:.2e})");
}

/// Compare a DIJ against the reference with the row-scaled metric; returns the worst scaled error.
fn check_dij(st: &BlState, tol: f64, what: &str) -> (f64, usize, usize) {
    let (n, nw, xd) = parse_dij(&fixture_path("xfoil_dij.dat"));
    assert_eq!((n, nw), (st.n, st.nw));
    let np = n + nw;
    let mut worst = (0.0_f64, 0usize, 0usize);
    for i in 1..=np {
        let rowmax = (1..=np).map(|j| xd[i][j].abs()).fold(0.0_f64, f64::max);
        for j in 1..=np {
            let (a, b) = (st.dij[i][j], xd[i][j]);
            let scaled = (a - b).abs() / a.abs().max(b.abs()).max(rowmax);
            if scaled > worst.0 {
                worst = (scaled, i, j);
            }
            assert_within(
                a,
                b,
                tol,
                rowmax,
                &format!("{what}: DIJ({i},{j}) [rowmax {rowmax:.3e}]"),
            );
        }
    }
    // structural: first wake row equals the TE row, bitwise
    for j in 1..=np {
        assert_eq!(
            st.dij[n + 1][j].to_bits(),
            st.dij[n][j].to_bits(),
            "DIJ(N+1,{j}) != DIJ(N,{j})"
        );
    }
    println!(
        "{what}: {np}x{np} DIJ within {tol:.0e} row-scaled (worst {:.2e} at ({},{}))",
        worst.0, worst.1, worst.2
    );
    worst
}

/// QDCALC in isolation: XFOIL's own wake nodes (position, normals, angles, QINVU) are
/// injected, so only PSWLIN/PSILIN/BAKSUB and the DIJ assembly are under test.
#[test]
fn test_qdcalc_matches_xfoil_dij_given_xfoil_wake() {
    let (mut st, mut sys, u) = prologue_to_specal();
    let f = parse_pointers(&fixture_path("xfoil_pointers.dat"), 1);
    for i in (st.n + 1)..=(st.n + st.nw) {
        st.x[i] = f.x[i];
        st.y[i] = f.y[i];
        st.s[i] = f.s[i];
        st.nx[i] = f.nx[i];
        st.ny[i] = f.ny[i];
        st.apanel[i] = f.apanel[i];
        st.qinvu[1][i] = u.qinvu1[i];
        st.qinvu[2][i] = u.qinvu2[i];
    }
    qdcalc(&mut st, &mut sys);
    check_dij(&st, TOL_LINALG, "qdcalc (XFOIL wake)");
}

/// The whole prologue from YFoil's geometry. DIJ's sensitivity to a wake node position is
/// O(1/distance) ≈ 3e2 here, and XYWAKE places the nodes within ~1e-12 of XFOIL's (its input
/// GAM is 4e-14 off after the LU solve), so the DIJ floor for this chained test is ~3e-10:
/// gated at TOL_SOLVER, not TOL_LINALG. The isolated test above is the translation gate.
#[test]
fn test_prologue_dij_from_yfoil_geometry() {
    let (mut st, mut sys, u) = prologue_to_specal();
    xywake(&mut st, 1.0);
    qwcalc(&mut st);
    for i in (st.n + 1)..=(st.n + st.nw) {
        assert_within(
            st.qinvu[1][i],
            u.qinvu1[i],
            TOL_PURE,
            1.0,
            &format!("wake QINVU({i},1)"),
        );
    }
    qdcalc(&mut st, &mut sys);
    check_dij(&st, TOL_SOLVER, "qdcalc (YFoil wake)");
}
