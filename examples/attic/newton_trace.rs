//! Newton trace - outputs same format as instrumented XFOIL for comparison
//!
//! Run: cargo run --release --example newton_trace > yfoil_newton_trace.dat

use std::fs::File;
use std::io::BufReader;
use yfoil::bl::{cf_laminar, cf_turbulent, cdiss_laminar, hkin, hstar_laminar, hstar_turbulent, Closure, FlowConditions, FlowRegime};
use yfoil::geometry::{panel_foil, Geometry};
use yfoil::panel::solve_inviscid;
use yfoil::solver::{extract_upper_surface, find_stagnation_point};

/// Load geometry from JSON file
fn load_geometry(path: &str) -> Geometry {
    let file = File::open(path).expect("Failed to open geometry file");
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).expect("Failed to parse geometry JSON")
}

/// Instrumented Newton solve for a single station
fn solve_station_instrumented(
    ibl: usize,
    theta1: f64,
    h1: f64,
    ue1: f64,
    n1: f64,
    ue2: f64,
    s1: f64,
    s2: f64,
    regime: FlowRegime,
    cond: &FlowConditions,
) -> (f64, f64, f64, f64, f64, f64, bool) {
    // Configuration matching XFOIL
    let max_iter = 25;
    let tol = 1e-5;
    let hk_min = 1.00005;
    let hk_max = 12.0;
    let theta_min = 1e-10;

    let ds = s2 - s1;
    let due_ds = (ue2 - ue1) / ds;
    let ue_avg = 0.5 * (ue1 + ue2);

    // Upstream quantities
    let (hk1, _, _) = hkin(h1, cond.msq);
    let hk1 = hk1.clamp(hk_min, hk_max);
    let rt1 = (ue1 * theta1 / cond.nu).max(1.0);

    // Initial guess
    let cf1_guess = match regime {
        FlowRegime::Laminar => cf_laminar(hk1, rt1, cond.msq).val,
        _ => cf_turbulent(hk1, rt1.max(200.0), cond.msq, 1.0).val,
    };
    let shape_term = (h1 + 2.0 - cond.msq) * theta1 / ue_avg * due_ds;
    let dtheta_ds_pred = cf1_guess / 2.0 - shape_term;
    let mut theta2 = (theta1 + ds * dtheta_ds_pred).clamp(theta1 * 0.5, theta1 * 2.0);

    let h_init = match regime {
        FlowRegime::Turbulent | FlowRegime::Wake if h1 > 2.0 => 1.5,
        FlowRegime::Turbulent | FlowRegime::Wake => h1.clamp(1.2, 2.5),
        FlowRegime::Laminar => h1,
    };
    let mut dstar2 = h_init * theta2;

    // Initial Ctau for turbulent (S2 in XFOIL)
    let ctau2 = 0.03;

    // Upstream closures
    let (hs1, cf1, di1) = match regime {
        FlowRegime::Laminar => (hstar_laminar(hk1, rt1, cond.msq), cf_laminar(hk1, rt1, cond.msq), cdiss_laminar(hk1, rt1)),
        _ => {
            let cf = cf_turbulent(hk1, rt1.max(200.0), cond.msq, 1.0);
            let hs = hstar_turbulent(hk1, rt1.max(200.0), cond.msq);
            let di = Closure {
                val: 0.5 * cf.val * ue1 / hs.val,
                val_hk: 0.0,
                val_rt: 0.0,
                val_msq: 0.0,
            };
            (hs, cf, di)
        }
    };

    // Newton iteration
    let mut converged = false;
    let mut _iterations = 0;

    for iter in 0..max_iter {
        _iterations = iter + 1;

        let h2 = dstar2 / theta2;
        let (hk2_raw, hk2_h, _) = hkin(h2, cond.msq);
        let hk2 = hk2_raw.clamp(hk_min, hk_max);
        let rt2 = (ue2 * theta2 / cond.nu).max(1.0);

        // XFOIL kinematic viscosity - for incompressible: V2 = 1/REYBL
        let v2 = cond.nu;

        // Get closure relations
        let (cf_res, hs_res, di_res) = match regime {
            FlowRegime::Laminar => {
                let cf = cf_laminar(hk2, rt2, cond.msq);
                let hs = hstar_laminar(hk2, rt2, cond.msq);
                let di = cdiss_laminar(hk2, rt2);
                (cf, hs, di)
            }
            FlowRegime::Turbulent | FlowRegime::Wake => {
                let cf = cf_turbulent(hk2, rt2.max(200.0), cond.msq, 1.0);
                let hs = hstar_turbulent(hk2, rt2.max(200.0), cond.msq);
                let di = Closure {
                    val: 0.5 * cf.val * ue2 / hs.val,
                    val_hk: 0.0,
                    val_rt: 0.0,
                    val_msq: 0.0,
                };
                (cf, hs, di)
            }
        };

        // XFOIL Us and Cq (slip velocity and equilibrium Ctau)
        // US2 = 0.5*HS2*(1 - (HK2-1)/(GBCON*H2))
        let gbcon = 6.0;
        let us2 = 0.5 * hs_res.val * (1.0 - (hk2 - 1.0) / (gbcon * h2)).min(0.98);

        // CQ2 = sqrt(CTCON*HS2*HKB*HKC^2 / (USB*H2*HK2^2))
        let ctcon = 0.5 / 6.7 / 6.7 / 0.75;
        let hkb = hk2 - 1.0;
        let hkc = (hk2 - 1.0).max(0.01);
        let usb = (1.0 - us2).max(0.02);
        let cq2 = (ctcon * hs_res.val * hkb * hkc * hkc / (usb * h2 * hk2 * hk2))
            .max(0.0)
            .sqrt();

        // Output in XFOIL format
        // PRIMARY: X2, U2, T2, D2, S2
        println!("STATION {:>3}   1 {:>3}", ibl, iter + 1);
        println!(
            "PRIMARY:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            s2, ue2, theta2, dstar2, ctau2
        );
        // KINEMATIC: M2, H2, HK2, RT2, V2
        println!(
            "KINEMATIC:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            cond.msq.sqrt(),
            h2,
            hk2,
            rt2,
            v2
        );
        // CLOSURE: HS2, US2, CQ2, CF2, DI2
        println!(
            "CLOSURE:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            hs_res.val, us2, cq2, cf_res.val, di_res.val
        );

        // Closure derivatives in (θ, δ*) coordinates
        let h2_t2 = -h2 / theta2;
        let h2_d2 = 1.0 / theta2;
        let hk2_t2 = hk2_h * h2_t2;
        let hk2_d2 = hk2_h * h2_d2;
        let rt2_t2 = rt2 / theta2;

        let cf2_t2 = cf_res.val_hk * hk2_t2 + cf_res.val_rt * rt2_t2;
        let cf2_d2 = cf_res.val_hk * hk2_d2;
        let hs2_t2 = hs_res.val_hk * hk2_t2 + hs_res.val_rt * rt2_t2;
        let hs2_d2 = hs_res.val_hk * hk2_d2;
        let di2_t2 = di_res.val_hk * hk2_t2 + di_res.val_rt * rt2_t2;
        let di2_d2 = di_res.val_hk * hk2_d2;

        // Average quantities
        let theta_avg = 0.5 * (theta1 + theta2);
        let h_avg = 0.5 * (h1 + h2);
        let s_avg = 0.5 * (s1 + s2);
        let hs_avg = 0.5 * (hs1.val + hs_res.val);

        // Midpoint Cf
        let hk_avg = 0.5 * (hk1 + hk2);
        let rt_avg = 0.5 * (rt1 + rt2);
        let cf_mid = match regime {
            FlowRegime::Laminar => cf_laminar(hk_avg, rt_avg, cond.msq),
            _ => cf_turbulent(hk_avg, rt_avg.max(200.0), cond.msq, 1.0),
        };

        let use_log_form = s1 > 1e-6;

        let (r1, dr1_dt2, dr1_dd2) = if use_log_form {
            let tlog = (theta2 / theta1).ln();
            let ulog = (ue2 / ue1).ln();
            let xlog = (s2 / s1).ln();
            let btmp = h_avg + 2.0 - cond.msq;

            let cfx = 0.5 * cf_mid.val * s_avg / theta_avg + 0.25 * (cf1.val * s1 / theta1 + cf_res.val * s2 / theta2);

            let r1_val = tlog + btmp * ulog - xlog * 0.5 * cfx;

            let cfx_t2 = -0.5 * cf_mid.val * s_avg / theta_avg.powi(2) * 0.5 - 0.25 * cf_res.val * s2 / theta2.powi(2)
                + 0.25 * cf2_t2 * s2 / theta2;
            let cfx_d2 = 0.25 * cf2_d2 * s2 / theta2;

            let dr1_dt2_val = 1.0 / theta2 + 0.5 * ulog * h2_t2 - xlog * 0.5 * cfx_t2;
            let dr1_dd2_val = 0.5 * ulog * h2_d2 - xlog * 0.5 * cfx_d2;

            (r1_val, dr1_dt2_val, dr1_dd2_val)
        } else {
            let btmp = h_avg + 2.0 - cond.msq;
            let mom_coef = btmp * theta_avg / ue_avg;
            let r1_val = (theta2 - theta1) / ds + mom_coef * due_ds - cf_mid.val / 2.0;

            let dr1_dt2_val =
                1.0 / ds + 0.5 * btmp / ue_avg * due_ds + 0.5 * h2_t2 * theta_avg / ue_avg * due_ds - cf2_t2 / 2.0;
            let dr1_dd2_val = 0.5 * h2_d2 * theta_avg / ue_avg * due_ds - cf2_d2 / 2.0;

            (r1_val, dr1_dt2_val, dr1_dd2_val)
        };

        // Shape equation with upwinding
        let upw = {
            let hdcon = 5.0 / hk2.powi(2);
            let arg = ((hk2 - 1.0) / (hk1 - 1.0).max(0.01)).abs();
            let hl = arg.ln();
            let hlsq = hl.powi(2).min(15.0);
            1.0 - 0.5 * (-hlsq * hdcon).exp()
        };

        let (r2, dr2_dt2, dr2_dd2) = if use_log_form {
            let hlog = (hs_res.val / hs1.val).ln();
            let ulog = (ue2 / ue1).ln();
            let xlog = (s2 / s1).ln();
            let btmp_h = 1.0 - h_avg;

            let xot1 = s1 / theta1;
            let xot2 = s2 / theta2;
            let dix = (1.0 - upw) * di1.val * xot1 + upw * di_res.val * xot2;
            let cfx = (1.0 - upw) * cf1.val * xot1 + upw * cf_res.val * xot2;

            let r2_val = hlog + btmp_h * ulog + xlog * (0.5 * cfx - dix);

            let dhlog_dt2 = hs2_t2 / hs_res.val;
            let dhlog_dd2 = hs2_d2 / hs_res.val;

            let dix_t2 = upw * (-di_res.val * xot2 / theta2 + di2_t2 * xot2);
            let dix_d2 = upw * di2_d2 * xot2;
            let cfx_t2 = upw * (-cf_res.val * xot2 / theta2 + cf2_t2 * xot2);
            let cfx_d2 = upw * cf2_d2 * xot2;

            let dr2_dt2_val = dhlog_dt2 - 0.5 * h2_t2 * ulog + xlog * (0.5 * cfx_t2 - dix_t2);
            let dr2_dd2_val = dhlog_dd2 - 0.5 * h2_d2 * ulog + xlog * (0.5 * cfx_d2 - dix_d2);

            (r2_val, dr2_dt2_val, dr2_dd2_val)
        } else {
            let dhs_ds = (hs_res.val - hs1.val) / ds;
            let shape_coef = hs_avg * (1.0 - h_avg) * theta_avg / ue_avg;
            let r2_val = theta_avg * dhs_ds + shape_coef * due_ds - 2.0 * di_res.val + hs_avg * cf_res.val / 2.0;

            let dr2_dt2_val = 0.5 * dhs_ds
                + theta_avg / ds * hs2_t2
                + (0.5 * hs2_t2 * (1.0 - h_avg) - 0.5 * hs_avg * h2_t2) * theta_avg / ue_avg * due_ds
                + hs_avg * cf2_t2 / 2.0
                - 2.0 * di2_t2;

            let dr2_dd2_val = theta_avg / ds * hs2_d2
                + (0.5 * hs2_d2 * (1.0 - h_avg) - 0.5 * hs_avg * h2_d2) * theta_avg / ue_avg * due_ds
                + hs_avg * cf2_d2 / 2.0
                - 2.0 * di2_d2;

            (r2_val, dr2_dt2_val, dr2_dd2_val)
        };

        // Output residual in XFOIL format (VS2 reduced to 2x2 for direct mode)
        // XFOIL outputs 4 residuals, but for direct mode (dUe=0), 4th is always 0
        println!("RESIDUAL:  {:.10E}  {:.10E}  {:.10E}  {:.10E}", 0.0, r1, r2, 0.0);
        // VS2 matrix - in direct mode, this is effectively 2x2 embedded in 3x5
        println!(
            "VS2_1:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            1.0, 0.0, 0.0, 0.0, 0.0
        );
        println!(
            "VS2_2:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            0.0, dr1_dt2, dr1_dd2, 0.0, 0.0
        );
        println!(
            "VS2_3:  {:.10E}  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            0.0, dr2_dt2, dr2_dd2, 0.0, 0.0
        );

        // Solve 2x2 Newton system
        let det = dr1_dt2 * dr2_dd2 - dr1_dd2 * dr2_dt2;

        let h_min = match regime {
            FlowRegime::Turbulent | FlowRegime::Wake => 1.2,
            FlowRegime::Laminar => 1.0,
        };
        let theta_max = theta1 * 3.0;

        let (d_theta, d_dstar) = if det.abs() < 1e-20 {
            (-r1 * 0.1, -r2 * 0.1)
        } else {
            (
                (-r1 * dr2_dd2 + r2 * dr1_dd2) / det,
                (-r2 * dr1_dt2 + r1 * dr2_dt2) / det,
            )
        };

        println!(
            "SOLUTION:  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            0.0, d_theta, d_dstar, 0.0
        );

        // XFOIL-style unified relaxation
        let dmax = (d_theta / theta2).abs().max((d_dstar / dstar2).abs());
        let rlx = if dmax > 0.3 { 0.3 / dmax } else { 1.0 };

        println!("RELAX:  {:.10E}  {:.10E}", dmax, rlx);

        theta2 = (theta2 + rlx * d_theta).clamp(theta_min, theta_max);
        dstar2 = (dstar2 + rlx * d_dstar).clamp(h_min * theta2, hk_max * theta2);

        println!(
            "UPDATED:  {:.10E}  {:.10E}  {:.10E}  {:.10E}",
            ctau2, theta2, dstar2, ue2
        );

        let residual = (r1.powi(2) + r2.powi(2)).sqrt();
        if residual < tol {
            converged = true;
            println!("CONVERGED: {:>3}  {:.10E}", iter + 1, dmax);
            break;
        }
    }

    let h2 = dstar2 / theta2;
    let (hk2, _, _) = hkin(h2, cond.msq);
    let rt2 = (ue2 * theta2 / cond.nu).max(1.0);
    let cf2 = match regime {
        FlowRegime::Laminar => cf_laminar(hk2, rt2, cond.msq).val,
        _ => cf_turbulent(hk2, rt2.max(200.0), cond.msq, 1.0).val,
    };

    (theta2, dstar2, h2, hk2, cf2, n1, converged)
}

fn main() {
    // Load same geometry as XFOIL
    let geom = load_geometry(".tmp/naca0012_xfoil_paneled.json");
    let airfoil = panel_foil(&geom);

    // Solve inviscid
    let inviscid = solve_inviscid(&airfoil);
    let qinv = inviscid.velocity_at_nodes(0.0);
    let stag_idx = find_stagnation_point(&airfoil, &qinv);

    let (_x_upper, _, s_upper, ue_upper) = extract_upper_surface(&airfoil, &qinv, stag_idx);

    let cond = FlowConditions::new(1_000_000.0, 0.0, 9.0, 1.0);

    // Header
    println!("# Newton iteration trace");
    println!("# FORMAT: STATION IBL IS ITBL");
    println!("# PRIMARY: X2 U2 T2 D2 S2");
    println!("# KINEMATIC: M2 H2 HK2 RT2 V2");
    println!("# CLOSURE: HS2 US2 CQ2 CF2 DI2");
    println!("# RESIDUAL: RES1 RES2 RES3");
    println!("# UPDATE: DZ1 DZ2 DZ3 DMAX RLX");

    // Initialize at stagnation point using XFOIL's Thwaites formula
    let s_init = s_upper[1].max(1e-10);
    let ue_init = ue_upper[1].max(0.01);
    let theta_init = (0.45 * cond.nu * s_init / (6.0 * ue_init)).sqrt();
    let h_init = 2.2;

    let mut theta = theta_init;
    let mut h = h_init;
    let mut n_amp = 0.0;
    let regime = FlowRegime::Laminar;

    // March stations
    for i in 1..ue_upper.len().min(20) {
        // Only first 20 stations for comparison
        let s1 = s_upper[i - 1];
        let s2 = s_upper[i];
        let ue1 = ue_upper[i - 1];
        let ue2 = ue_upper[i];

        let (theta2, _dstar2, h2, _hk2, _cf2, n2, _converged) = solve_station_instrumented(
            i + 1, // IBL = i+1 (1-indexed, starting at 2)
            theta,
            h,
            ue1,
            n_amp,
            ue2,
            s1,
            s2,
            regime,
            &cond,
        );

        theta = theta2;
        h = h2;
        n_amp = n2;
    }
}
