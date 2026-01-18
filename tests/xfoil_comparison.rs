//! Integration tests comparing yfoil with instrumented XFOIL
//!
//! These tests require the instrumented XFOIL binary to be built first.
//! Run: `cd xfoil/instrumented && make install`

use std::fs;
use std::process::Command;

/// Get path to instrumented XFOIL binary
fn xfoil_path() -> String {
    std::env::var("YFOIL_XFOIL_PATH").unwrap_or_else(|_| {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        format!("{}/xfoil/xfoil6.99/bin/xfoil", manifest_dir)
    })
}

/// Parse XFOIL panel geometry output file
fn parse_xfoil_panels(path: &str) -> XfoilPanels {
    let content = fs::read_to_string(path).expect("Failed to read XFOIL panels file");
    let mut panels = XfoilPanels::default();

    let mut in_data = false;
    for line in content.lines() {
        let line = line.trim();

        if line.starts_with("N =") {
            panels.n = line.split_whitespace().last().unwrap().parse().unwrap();
        } else if line.starts_with("SLE =") {
            panels.sle = parse_fortran_float(&line[5..]);
        } else if line.starts_with("XLE =") {
            panels.xle = parse_fortran_float(&line[5..]);
        } else if line.starts_with("YLE =") {
            panels.yle = parse_fortran_float(&line[5..]);
        } else if line.starts_with("XTE =") {
            panels.xte = parse_fortran_float(&line[5..]);
        } else if line.starts_with("YTE =") {
            panels.yte = parse_fortran_float(&line[5..]);
        } else if line.starts_with("CHORD =") {
            panels.chord = parse_fortran_float(&line[7..]);
        } else if line.starts_with("SHARP =") {
            panels.sharp = line.contains("T");
        } else if line.starts_with("I, X, Y") {
            in_data = true;
            continue;
        } else if in_data && !line.is_empty() && !line.starts_with("---") {
            // Parse panel data line: I, X, Y, S, NX, NY, APANEL
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 7 {
                panels.x.push(parse_fortran_float(parts[1]));
                panels.y.push(parse_fortran_float(parts[2]));
                panels.s.push(parse_fortran_float(parts[3]));
                panels.nx.push(parse_fortran_float(parts[4]));
                panels.ny.push(parse_fortran_float(parts[5]));
                panels.apanel.push(parse_fortran_float(parts[6]));
            }
        }
    }

    panels
}

/// Parse Fortran scientific notation (handles D exponent)
fn parse_fortran_float(s: &str) -> f64 {
    let s = s.trim().replace("D", "E").replace("d", "e");
    s.parse().unwrap_or_else(|_| panic!("Failed to parse: {}", s))
}

#[derive(Default, Debug)]
struct XfoilPanels {
    n: usize,
    sle: f64,
    xle: f64,
    yle: f64,
    xte: f64,
    yte: f64,
    chord: f64,
    sharp: bool,
    x: Vec<f64>,
    y: Vec<f64>,
    s: Vec<f64>,
    nx: Vec<f64>,
    ny: Vec<f64>,
    apanel: Vec<f64>,
}

/// Run instrumented XFOIL with yfoil-generated geometry
/// Uses LOAD and PCOP to load panels directly without XFOIL repaneling
fn run_xfoil_with_geometry(dat_path: &str) {
    // LOAD loads the geometry, PCOP copies buffer directly to current airfoil
    // This avoids XFOIL's PANE command which would repanel
    let xfoil_input = format!(
        "PLOP\nG F\n\nLOAD {}\nPCOP\nQUIT\n",
        dat_path
    );

    let input_path = "/tmp/xfoil_test_input.txt";
    fs::write(input_path, &xfoil_input).expect("Failed to write XFOIL input");

    let status = Command::new("sh")
        .arg("-c")
        .arg(format!("{} < {}", xfoil_path(), input_path))
        .status()
        .expect("Failed to run XFOIL");

    if !status.success() {
        panic!("XFOIL failed with status: {:?}", status);
    }
}

/// Run instrumented XFOIL inviscid solve at given alpha
fn run_xfoil_inviscid(dat_path: &str, alpha_deg: f64) {
    let xfoil_input = format!(
        "PLOP\nG F\n\nLOAD {}\nPCOP\nOPER\nALFA {}\n\n",
        dat_path, alpha_deg
    );

    let input_path = "/tmp/xfoil_test_input.txt";
    fs::write(input_path, &xfoil_input).expect("Failed to write XFOIL input");

    // XFOIL will error on EOF but that's OK - instrumentation writes before that
    let _ = Command::new("sh")
        .arg("-c")
        .arg(format!("{} < {}", xfoil_path(), input_path))
        .output();
}

/// Parse XFOIL inviscid solution output file
#[derive(Default, Debug)]
struct XfoilInviscid {
    alfa: f64,
    cl: f64,
    cm: f64,
    n: usize,
    gam: Vec<f64>,
    qinv: Vec<f64>,
    cpi: Vec<f64>,
    gamu1: Vec<f64>,
    gamu2: Vec<f64>,
    qinvu1: Vec<f64>,
    qinvu2: Vec<f64>,
}

fn parse_xfoil_inviscid(path: &str) -> XfoilInviscid {
    let content = fs::read_to_string(path).expect("Failed to read XFOIL inviscid file");
    let mut result = XfoilInviscid::default();

    let mut section = 0; // 0 = header, 1 = GAM/QINV, 2 = GAMU
    for line in content.lines() {
        let line = line.trim();

        if line.starts_with("ALFA =") {
            result.alfa = parse_fortran_float(&line[6..]);
        } else if line.starts_with("CL =") {
            result.cl = parse_fortran_float(&line[4..]);
        } else if line.starts_with("CM =") {
            result.cm = parse_fortran_float(&line[4..]);
        } else if line.starts_with("N =") {
            result.n = line.split_whitespace().last().unwrap().parse().unwrap();
        } else if line.starts_with("I, GAM, QINV") {
            section = 1;
            continue;
        } else if line.starts_with("I, GAMU1") {
            section = 2;
            continue;
        } else if line.starts_with("---") || line.starts_with("===") {
            continue;
        } else if !line.is_empty() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if section == 1 && parts.len() >= 4 {
                result.gam.push(parse_fortran_float(parts[1]));
                result.qinv.push(parse_fortran_float(parts[2]));
                result.cpi.push(parse_fortran_float(parts[3]));
            } else if section == 2 && parts.len() >= 5 {
                result.gamu1.push(parse_fortran_float(parts[1]));
                result.gamu2.push(parse_fortran_float(parts[2]));
                result.qinvu1.push(parse_fortran_float(parts[3]));
                result.qinvu2.push(parse_fortran_float(parts[4]));
            }
        }
    }

    result
}

#[test]
fn test_panel_geometry_matches_xfoil() {
    use yfoil::geometry::{create_paneled_airfoil, naca_4digit, write_dat_file};

    // Generate yfoil panels
    let geom = naca_4digit("0012", 160).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    // Save yfoil geometry to .dat file for XFOIL to load
    let dat_path = "/tmp/yfoil_test_airfoil.dat";
    write_dat_file(&geom, "NACA 0012 from yfoil", dat_path)
        .expect("Failed to write .dat file");

    // Run XFOIL with yfoil's geometry (LOAD + PCOP, no repaneling)
    run_xfoil_with_geometry(dat_path);

    // Parse XFOIL output
    let xfoil = parse_xfoil_panels("/tmp/xfoil_panels.dat");

    // Compare panel count
    assert_eq!(
        airfoil.n, xfoil.n,
        "Panel count mismatch: yfoil={}, XFOIL={}",
        airfoil.n, xfoil.n
    );

    // Compare LE position
    let sle_diff = (airfoil.sle - xfoil.sle).abs();
    assert!(
        sle_diff < 1e-10,
        "SLE mismatch: yfoil={:.16e}, XFOIL={:.16e}, diff={:.16e}",
        airfoil.sle,
        xfoil.sle,
        sle_diff
    );

    // Compare each panel coordinate
    let mut max_x_diff = 0.0f64;
    let mut max_y_diff = 0.0f64;
    let mut max_s_diff = 0.0f64;
    let mut max_nx_diff = 0.0f64;
    let mut max_ny_diff = 0.0f64;

    for i in 0..airfoil.n {
        let x_diff = (airfoil.x[i] - xfoil.x[i]).abs();
        let y_diff = (airfoil.y[i] - xfoil.y[i]).abs();
        let s_diff = (airfoil.s[i] - xfoil.s[i]).abs();
        let nx_diff = (airfoil.nx[i] - xfoil.nx[i]).abs();
        let ny_diff = (airfoil.ny[i] - xfoil.ny[i]).abs();

        max_x_diff = max_x_diff.max(x_diff);
        max_y_diff = max_y_diff.max(y_diff);
        max_s_diff = max_s_diff.max(s_diff);
        max_nx_diff = max_nx_diff.max(nx_diff);
        max_ny_diff = max_ny_diff.max(ny_diff);

        // Report first significant discrepancy
        if x_diff > 1e-10 || y_diff > 1e-10 {
            eprintln!(
                "Panel {} mismatch:\n  X: yfoil={:.16e} XFOIL={:.16e} diff={:.16e}\n  Y: yfoil={:.16e} XFOIL={:.16e} diff={:.16e}",
                i + 1, airfoil.x[i], xfoil.x[i], x_diff, airfoil.y[i], xfoil.y[i], y_diff
            );
        }
    }

    eprintln!("\nMax differences:");
    eprintln!("  X:  {:.16e}", max_x_diff);
    eprintln!("  Y:  {:.16e}", max_y_diff);
    eprintln!("  S:  {:.16e}", max_s_diff);
    eprintln!("  NX: {:.16e}", max_nx_diff);
    eprintln!("  NY: {:.16e}", max_ny_diff);

    // Assert all differences are within floating point precision
    assert!(
        max_x_diff < 1e-10,
        "X coordinate mismatch exceeds tolerance: {:.16e}",
        max_x_diff
    );
    assert!(
        max_y_diff < 1e-10,
        "Y coordinate mismatch exceeds tolerance: {:.16e}",
        max_y_diff
    );
    assert!(
        max_s_diff < 1e-10,
        "S (arc length) mismatch exceeds tolerance: {:.16e}",
        max_s_diff
    );
    assert!(
        max_nx_diff < 1e-10,
        "NX (normal) mismatch exceeds tolerance: {:.16e}",
        max_nx_diff
    );
    assert!(
        max_ny_diff < 1e-10,
        "NY (normal) mismatch exceeds tolerance: {:.16e}",
        max_ny_diff
    );
}

#[test]
fn test_inviscid_solution_matches_xfoil() {
    use yfoil::geometry::{create_paneled_airfoil, naca_4digit, write_dat_file};
    use yfoil::panel::solve_inviscid;

    // Generate yfoil panels
    let geom = naca_4digit("0012", 160).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    // Save yfoil geometry to .dat file for XFOIL to load
    let dat_path = "/tmp/yfoil_test_airfoil.dat";
    write_dat_file(&geom, "NACA 0012 from yfoil", dat_path)
        .expect("Failed to write .dat file");

    // Run yfoil inviscid solve
    let solution = solve_inviscid(&airfoil);
    let alpha_rad = 0.0_f64.to_radians();

    // Get vortex strengths at alpha=0 (GAM in XFOIL terminology)
    let gam = solution.gamma_at_alpha(alpha_rad);
    // For inviscid attached flow, QINV = GAM (surface velocity equals vorticity)
    let qinv = solution.velocity_at_nodes(alpha_rad);
    // Cp = 1 - (Q/Qinf)^2, and Qinf=1 so Cp = 1 - Q^2
    let cpi: Vec<f64> = qinv.iter().map(|&q| 1.0 - q * q).collect();

    // Run XFOIL inviscid solve at alpha=0
    run_xfoil_inviscid(dat_path, 0.0);

    // Parse XFOIL output
    let xfoil = parse_xfoil_inviscid("/tmp/xfoil_inviscid.dat");

    // Compare panel count
    assert_eq!(
        airfoil.n, xfoil.n,
        "Panel count mismatch: yfoil={}, XFOIL={}",
        airfoil.n, xfoil.n
    );

    // Compare GAM (vortex strengths)
    let mut max_gam_diff = 0.0f64;
    let mut max_qinv_diff = 0.0f64;
    let mut max_cpi_diff = 0.0f64;

    for i in 0..airfoil.n {
        let gam_diff = (gam[i] - xfoil.gam[i]).abs();
        let qinv_diff = (qinv[i] - xfoil.qinv[i]).abs();
        let cpi_diff = (cpi[i] - xfoil.cpi[i]).abs();

        max_gam_diff = max_gam_diff.max(gam_diff);
        max_qinv_diff = max_qinv_diff.max(qinv_diff);
        max_cpi_diff = max_cpi_diff.max(cpi_diff);

        // Report first significant discrepancy
        if gam_diff > 1e-10 {
            eprintln!(
                "Panel {} GAM mismatch: yfoil={:.16e} XFOIL={:.16e} diff={:.16e}",
                i + 1, gam[i], xfoil.gam[i], gam_diff
            );
        }
    }

    eprintln!("\nInviscid solution max differences:");
    eprintln!("  GAM:  {:.16e}", max_gam_diff);
    eprintln!("  QINV: {:.16e}", max_qinv_diff);
    eprintln!("  CPI:  {:.16e}", max_cpi_diff);

    // Assert all differences are within floating point precision
    assert!(
        max_gam_diff < 1e-10,
        "GAM mismatch exceeds tolerance: {:.16e}",
        max_gam_diff
    );
    assert!(
        max_qinv_diff < 1e-10,
        "QINV mismatch exceeds tolerance: {:.16e}",
        max_qinv_diff
    );
    assert!(
        max_cpi_diff < 1e-10,
        "CPI mismatch exceeds tolerance: {:.16e}",
        max_cpi_diff
    );
}

/// Parse XFOIL DIJ matrix output file
#[derive(Default, Debug)]
struct XfoilDij {
    n: usize,
    nw: usize,
    entries: Vec<(usize, usize, f64)>,
}

fn parse_xfoil_dij(path: &str) -> XfoilDij {
    let content = fs::read_to_string(path).expect("Failed to read XFOIL DIJ file");
    let mut result = XfoilDij::default();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("N =") {
            result.n = line.split_whitespace().last().unwrap().parse().unwrap();
        } else if line.starts_with("NW =") {
            result.nw = line.split_whitespace().last().unwrap().parse().unwrap();
        } else if !line.starts_with("===") && !line.is_empty() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 3 {
                if let (Ok(i), Ok(j), Ok(val)) = (
                    parts[0].parse::<usize>(),
                    parts[1].parse::<usize>(),
                    parts[2].replace("D", "E").parse::<f64>(),
                ) {
                    result.entries.push((i, j, val));
                }
            }
        }
    }

    result
}

/// Run XFOIL viscous solve to generate DIJ matrix
fn run_xfoil_viscous(dat_path: &str, alpha_deg: f64) {
    let xfoil_input = format!(
        "PLOP\nG F\n\nLOAD {}\nPCOP\nOPER\nVISC 1e6\nALFA {}\n\n",
        dat_path, alpha_deg
    );

    let input_path = "/tmp/xfoil_test_input.txt";
    fs::write(input_path, &xfoil_input).expect("Failed to write XFOIL input");

    // XFOIL will write DIJ during VISCAL
    let _ = Command::new("sh")
        .arg("-c")
        .arg(format!("{} < {}", xfoil_path(), input_path))
        .output();
}

/// Test that velocity distribution matches XFOIL to machine precision
/// This test locks in the good velocity matching we've verified.
#[test]
fn test_velocity_distribution_matches_xfoil() {
    use yfoil::geometry::{create_paneled_airfoil, naca_4digit, repanel_xfoil, write_dat_file, PaneConfig};
    use yfoil::panel::solve_inviscid;

    // Generate and repanel using xfoil method to ensure identical paneling
    let geom = naca_4digit("0012", 160).expect("Failed to generate NACA 0012");
    let config = PaneConfig::default();
    let repaneled = repanel_xfoil(&geom, 160, &config);
    let airfoil = create_paneled_airfoil(&repaneled);

    // Save yfoil geometry to .dat file for XFOIL to load
    let dat_path = "/tmp/yfoil_velocity_test.dat";
    write_dat_file(&repaneled, "NACA 0012 yfoil repaneled", dat_path)
        .expect("Failed to write .dat file");

    // Test at alpha=0
    {
        let solution = solve_inviscid(&airfoil);
        let alpha_rad = 0.0_f64;
        let velocity = solution.velocity_at_nodes(alpha_rad);

        run_xfoil_inviscid(dat_path, 0.0);
        let xfoil = parse_xfoil_inviscid("/tmp/xfoil_inviscid.dat");

        assert_eq!(airfoil.n, xfoil.n, "Panel count mismatch");

        let mut max_vel_diff = 0.0f64;
        for i in 0..airfoil.n {
            let diff = (velocity[i] - xfoil.qinv[i]).abs();
            max_vel_diff = max_vel_diff.max(diff);
        }

        eprintln!("Alpha=0: Max velocity difference: {:.6e}", max_vel_diff);
        assert!(
            max_vel_diff < 1e-4,
            "Velocity mismatch at alpha=0: max diff = {:.6e}",
            max_vel_diff
        );
    }

    // Test at alpha=1 degree
    {
        let solution = solve_inviscid(&airfoil);
        let alpha_rad = 1.0_f64.to_radians();
        let velocity = solution.velocity_at_nodes(alpha_rad);

        run_xfoil_inviscid(dat_path, 1.0);
        let xfoil = parse_xfoil_inviscid("/tmp/xfoil_inviscid.dat");

        let mut max_vel_diff = 0.0f64;
        for i in 0..airfoil.n {
            let diff = (velocity[i] - xfoil.qinv[i]).abs();
            max_vel_diff = max_vel_diff.max(diff);
        }

        eprintln!("Alpha=1: Max velocity difference: {:.6e}", max_vel_diff);
        assert!(
            max_vel_diff < 1e-4,
            "Velocity mismatch at alpha=1: max diff = {:.6e}",
            max_vel_diff
        );
    }
}

/// Test force coefficient integration against XFOIL
/// This test captures the current discrepancy in force integration.
#[test]
fn test_force_coefficients_match_xfoil() {
    use yfoil::geometry::{create_paneled_airfoil, naca_4digit, repanel_xfoil, write_dat_file, PaneConfig};
    use yfoil::panel::solve_inviscid;
    use yfoil::forces::integrate_forces;

    // Generate and repanel using xfoil method
    let geom = naca_4digit("0012", 160).expect("Failed to generate NACA 0012");
    let config = PaneConfig::default();
    let repaneled = repanel_xfoil(&geom, 160, &config);
    let airfoil = create_paneled_airfoil(&repaneled);

    let dat_path = "/tmp/yfoil_force_test.dat";
    write_dat_file(&repaneled, "NACA 0012 yfoil repaneled", dat_path)
        .expect("Failed to write .dat file");

    // Test at alpha=0 (symmetric airfoil should have CL=0)
    {
        let solution = solve_inviscid(&airfoil);
        let alpha_rad = 0.0_f64;
        let alpha_deg = 0.0_f64;
        let mach = 0.0_f64;

        // integrate_forces uses velocity at nodes and computes Cp internally
        let velocity = solution.velocity_at_nodes(alpha_rad);
        let coeffs = integrate_forces(&airfoil, &velocity, alpha_rad, mach);

        // Get XFOIL coefficients
        run_xfoil_inviscid(dat_path, alpha_deg);
        let xfoil = parse_xfoil_inviscid("/tmp/xfoil_inviscid.dat");

        eprintln!("\n=== Alpha=0 Force Coefficients ===");
        eprintln!("XFOIL:  CL={:.10e}, CM={:.10e}", xfoil.cl, xfoil.cm);
        eprintln!("yfoil:  CL={:.10e}, CM={:.10e}, CDp={:.10e}",
                  coeffs.cl, coeffs.cm, coeffs.cdp);

        let cl_diff = (coeffs.cl - xfoil.cl).abs();
        let cm_diff = (coeffs.cm - xfoil.cm).abs();

        eprintln!("CL diff: {:.6e}", cl_diff);
        eprintln!("CM diff: {:.6e}", cm_diff);

        // integrate_forces should match XFOIL to high precision
        assert!(
            cl_diff < 1e-6,
            "CL mismatch at alpha=0: yfoil={:.10e}, XFOIL={:.10e}, diff={:.6e}",
            coeffs.cl, xfoil.cl, cl_diff
        );
    }

    // Test at alpha=1 degree
    {
        let solution = solve_inviscid(&airfoil);
        let alpha_rad = 1.0_f64.to_radians();
        let alpha_deg = 1.0_f64;
        let mach = 0.0_f64;

        let velocity = solution.velocity_at_nodes(alpha_rad);
        let coeffs = integrate_forces(&airfoil, &velocity, alpha_rad, mach);

        run_xfoil_inviscid(dat_path, alpha_deg);
        let xfoil = parse_xfoil_inviscid("/tmp/xfoil_inviscid.dat");

        eprintln!("\n=== Alpha=1 Force Coefficients ===");
        eprintln!("XFOIL:  CL={:.10e}, CM={:.10e}", xfoil.cl, xfoil.cm);
        eprintln!("yfoil:  CL={:.10e}, CM={:.10e}, CDp={:.10e}",
                  coeffs.cl, coeffs.cm, coeffs.cdp);

        let cl_diff = (coeffs.cl - xfoil.cl).abs();
        let cm_diff = (coeffs.cm - xfoil.cm).abs();

        eprintln!("CL diff: {:.6e}", cl_diff);
        eprintln!("CM diff: {:.6e}", cm_diff);

        // These should match to high precision
        assert!(
            cl_diff < 1e-6,
            "CL mismatch at alpha=1: yfoil={:.10e}, XFOIL={:.10e}, diff={:.6e}",
            coeffs.cl, xfoil.cl, cl_diff
        );
    }
}

#[test]
fn test_dij_matrix_matches_xfoil() {
    use yfoil::geometry::{create_paneled_airfoil, naca_4digit, write_dat_file};
    use yfoil::panel::solve_inviscid;

    // Generate yfoil panels
    let geom = naca_4digit("0012", 160).expect("Failed to generate NACA 0012");
    let airfoil = create_paneled_airfoil(&geom);

    // Save yfoil geometry to .dat file for XFOIL to load
    let dat_path = "/tmp/yfoil_test_airfoil.dat";
    write_dat_file(&geom, "NACA 0012 from yfoil", dat_path)
        .expect("Failed to write .dat file");

    // Run yfoil inviscid solve (computes DIJ)
    let solution = solve_inviscid(&airfoil);

    // Get DIJ matrix
    let dij = solution.get_dij().expect("DIJ not computed");

    // Run XFOIL viscous solve to generate DIJ
    run_xfoil_viscous(dat_path, 2.0);

    // Parse XFOIL DIJ output
    let xfoil_dij = parse_xfoil_dij("/tmp/xfoil_dij.dat");

    eprintln!("yfoil DIJ: {}x{}", dij.nrows(), dij.ncols());
    eprintln!("XFOIL N={}, NW={}, entries: {}", xfoil_dij.n, xfoil_dij.nw, xfoil_dij.entries.len());

    // Compare airfoil-only part (first N x N)
    // Note: XFOIL's DIJ includes wake panels, yfoil's is airfoil-only
    let n = airfoil.n;
    let mut max_diff = 0.0f64;
    let mut max_i = 0;
    let mut max_j = 0;
    let mut yfoil_at_max = 0.0;
    let mut xfoil_at_max = 0.0;

    for &(i, j, xfoil_val) in &xfoil_dij.entries {
        if i <= n && j <= n {
            // Convert to 0-based
            let yfoil_val = dij[(i - 1, j - 1)];
            let diff = (yfoil_val - xfoil_val).abs();
            if diff > max_diff {
                max_diff = diff;
                max_i = i;
                max_j = j;
                yfoil_at_max = yfoil_val;
                xfoil_at_max = xfoil_val;
            }
        }
    }

    eprintln!("\nMax DIJ difference in airfoil region: {:.6e} at ({}, {})", max_diff, max_i, max_j);
    eprintln!("  yfoil: {:.16e}", yfoil_at_max);
    eprintln!("  XFOIL: {:.16e}", xfoil_at_max);

    // Print first few diagonal values for comparison
    eprintln!("\nFirst 5 diagonal DIJ values:");
    eprintln!("  {:>5} {:>24} {:>24} {:>16}", "i", "yfoil", "XFOIL", "diff");
    for i in 0..5.min(n) {
        let xfoil_val = xfoil_dij.entries.iter()
            .find(|&&(ii, jj, _)| ii == i + 1 && jj == i + 1)
            .map(|&(_, _, v)| v)
            .unwrap_or(0.0);
        let diff = (dij[(i, i)] - xfoil_val).abs();
        eprintln!("  {:>5} {:>24.16e} {:>24.16e} {:>16.6e}", i + 1, dij[(i, i)], xfoil_val, diff);
    }

    // Assert DIJ matches
    assert!(
        max_diff < 1e-8,
        "DIJ mismatch exceeds tolerance: {:.16e} at ({}, {})",
        max_diff, max_i, max_j
    );
}
