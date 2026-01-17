//! Generate test fixtures from instrumented XFOIL
//!
//! This example runs XFOIL with instrumentation and parses the output
//! into JSON fixtures for validation testing.
//!
//! Usage:
//!   cargo run --example generate_fixtures
//!
//! The fixtures are written to tests/fixtures/

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

const XFOIL_BINARY: &str = "xfoil/xfoil6.99/bin/xfoil";
const FIXTURES_DIR: &str = "tests/fixtures";

#[derive(Debug, Serialize, Deserialize)]
struct TestInput {
    airfoil: String,
    alpha_deg: f64,
    reynolds: f64,
    mach: f64,
    n_panels: usize,
    n_crit: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct GeometryFixture {
    n: usize,
    x: Vec<f64>,
    y: Vec<f64>,
    s: Vec<f64>,
    nx: Vec<f64>,
    ny: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
struct InviscidFixture {
    alpha_rad: f64,
    gamma: Vec<f64>,
    qinv: Vec<f64>,
    cpi: Vec<f64>,
    cl_inv: f64,
    cm_inv: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct BLStation {
    ibl: usize,
    x: f64,
    s: f64,
    ue: f64,
    delta_star: f64,
    theta: f64,
    hk: f64,
    cf: f64,
    ctau: f64,
    mass: f64,
    regime: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct BLSide {
    n_stations: usize,
    stations: Vec<BLStation>,
    itran: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct BLFixture {
    upper: BLSide,
    lower: BLSide,
}

#[derive(Debug, Serialize, Deserialize)]
struct ViscalIteration {
    iter: usize,
    alpha_deg: f64,
    cl: f64,
    cd: f64,
    cdf: f64,
    cdp: f64,
    cm: f64,
    rmsbl: f64,
    rmxbl: f64,
    rlx: f64,
}

#[derive(Debug, Serialize, Deserialize)]
struct ViscalFixture {
    n_iterations: usize,
    converged: bool,
    iterations: Vec<ViscalIteration>,
}

#[derive(Debug, Serialize, Deserialize)]
struct FinalFixture {
    alpha_deg: f64,
    cl: f64,
    cd: f64,
    cdf: f64,
    cdp: f64,
    cm: f64,
    xtr_upper: f64,
    xtr_lower: f64,
}

/// Test case configuration
struct TestCase {
    airfoil: &'static str,
    naca_code: &'static str,
    alpha_deg: f64,
    reynolds: f64,
    dir_name: &'static str,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define test cases
    let test_cases = vec![
        TestCase {
            airfoil: "NACA 0012",
            naca_code: "0012",
            alpha_deg: 0.0,
            reynolds: 1e6,
            dir_name: "naca0012/alpha_0_re_1e6",
        },
        TestCase {
            airfoil: "NACA 0012",
            naca_code: "0012",
            alpha_deg: 2.0,
            reynolds: 1e6,
            dir_name: "naca0012/alpha_2_re_1e6",
        },
        TestCase {
            airfoil: "NACA 0012",
            naca_code: "0012",
            alpha_deg: 5.0,
            reynolds: 1e6,
            dir_name: "naca0012/alpha_5_re_1e6",
        },
        TestCase {
            airfoil: "NACA 4412",
            naca_code: "4412",
            alpha_deg: 0.0,
            reynolds: 1e6,
            dir_name: "naca4412/alpha_0_re_1e6",
        },
        TestCase {
            airfoil: "NACA 4412",
            naca_code: "4412",
            alpha_deg: 4.0,
            reynolds: 1e6,
            dir_name: "naca4412/alpha_4_re_1e6",
        },
    ];

    println!("Generating test fixtures from instrumented XFOIL...\n");

    for case in &test_cases {
        println!("Running: {} at alpha = {}°", case.airfoil, case.alpha_deg);
        generate_fixture(case)?;
        println!("  Done.\n");
    }

    println!("All fixtures generated successfully.");
    Ok(())
}

fn generate_fixture(case: &TestCase) -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = Path::new(FIXTURES_DIR).join(case.dir_name);
    fs::create_dir_all(&output_dir)?;

    // Build XFOIL script
    let script = format!(
        r#"PLOP
G F

NACA {}
PANE
OPER
VISC {}
ITER 50
ALFA {}

QUIT
"#,
        case.naca_code, case.reynolds, case.alpha_deg
    );

    // Run XFOIL
    let xfoil_path = Path::new(XFOIL_BINARY);
    if !xfoil_path.exists() {
        return Err(format!("XFOIL binary not found at: {}", XFOIL_BINARY).into());
    }

    let mut child = Command::new(xfoil_path)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    // Write script to XFOIL stdin
    if let Some(stdin) = child.stdin.as_mut() {
        stdin.write_all(script.as_bytes())?;
    }

    let output = child.wait_with_output()?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("XFOIL stderr: {}", stderr);
    }

    // Parse stdout for final results
    let stdout = String::from_utf8_lossy(&output.stdout);

    // Write input.json
    let input = TestInput {
        airfoil: case.airfoil.to_string(),
        alpha_deg: case.alpha_deg,
        reynolds: case.reynolds,
        mach: 0.0,
        n_panels: 160,
        n_crit: 9.0,
    };
    write_json(&output_dir.join("input.json"), &input)?;

    // Parse instrumentation files and write fixtures
    parse_and_write_viscal_fixture(&output_dir)?;
    parse_and_write_final_fixture(&output_dir, &stdout, case)?;

    Ok(())
}

fn parse_and_write_viscal_fixture(output_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let viscal_path = Path::new("/tmp/xfoil_viscal_iter.dat");
    if !viscal_path.exists() {
        eprintln!("  Warning: viscal iteration log not found");
        return Ok(());
    }

    let file = File::open(viscal_path)?;
    let reader = BufReader::new(file);

    let mut iterations = Vec::new();
    let mut current_iter: Option<ViscalIteration> = None;
    let mut converged = false;

    for line in reader.lines() {
        let line = line?;

        if line.contains("--- ITERATION") {
            // Save previous iteration
            if let Some(iter) = current_iter.take() {
                iterations.push(iter);
            }

            // Extract iteration number
            let iter_num: usize = line
                .split_whitespace()
                .find(|s| s.parse::<usize>().is_ok())
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);

            current_iter = Some(ViscalIteration {
                iter: iter_num,
                alpha_deg: 0.0,
                cl: 0.0,
                cd: 0.0,
                cdf: 0.0,
                cdp: 0.0,
                cm: 0.0,
                rmsbl: 0.0,
                rmxbl: 0.0,
                rlx: 1.0,
            });
        }

        if let Some(ref mut iter) = current_iter {
            let trimmed = line.trim();
            if trimmed.starts_with("ALFA =") {
                iter.alpha_deg = parse_fortran_value(trimmed).to_degrees();
            } else if trimmed.starts_with("CL =") {
                iter.cl = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("CD =") && !trimmed.starts_with("CDF") && !trimmed.starts_with("CDP") {
                iter.cd = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("CDF =") {
                iter.cdf = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("CDP =") {
                iter.cdp = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("CM =") {
                iter.cm = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("RMSBL =") {
                iter.rmsbl = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("RMXBL =") {
                iter.rmxbl = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("RLX =") {
                iter.rlx = parse_fortran_value(trimmed);
            }
        }

        let trimmed = line.trim();
        if trimmed.contains("CONVERGED") && trimmed.contains("T") {
            converged = true;
        }
    }

    // Don't forget last iteration
    if let Some(iter) = current_iter {
        iterations.push(iter);
    }

    let fixture = ViscalFixture {
        n_iterations: iterations.len(),
        converged,
        iterations,
    };

    write_json(&output_dir.join("viscal_iters.json"), &fixture)?;
    Ok(())
}

fn parse_and_write_final_fixture(
    output_dir: &Path,
    stdout: &str,
    case: &TestCase,
) -> Result<(), Box<dyn std::error::Error>> {
    // Parse final values from XFOIL output
    let mut final_result = FinalFixture {
        alpha_deg: case.alpha_deg,
        cl: 0.0,
        cd: 0.0,
        cdf: 0.0,
        cdp: 0.0,
        cm: 0.0,
        xtr_upper: 0.0,
        xtr_lower: 0.0,
    };

    // Parse from the viscal log for final values
    let viscal_path = Path::new("/tmp/xfoil_viscal_iter.dat");
    if viscal_path.exists() {
        let content = fs::read_to_string(viscal_path)?;

        // Variables to track last iteration values for CM
        let mut last_cm = 0.0;

        // Look for FINAL CL, CD etc and track CM from iterations
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with("FINAL CL =") {
                final_result.cl = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("FINAL CD =") && !trimmed.contains("CDF") && !trimmed.contains("CDP") {
                final_result.cd = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("FINAL CDF =") {
                final_result.cdf = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("FINAL CDP =") {
                final_result.cdp = parse_fortran_value(trimmed);
            } else if trimmed.starts_with("CM =") {
                last_cm = parse_fortran_value(trimmed);
            }
        }

        // Use last iteration CM value
        final_result.cm = last_cm;
    }

    // Parse transition locations from XFOIL stdout
    // Look for lines like "   Side 1  free  transition at x/c =  0.4123"
    for line in stdout.lines() {
        if line.contains("transition at x/c") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            // Find the value after "="
            for (i, part) in parts.iter().enumerate() {
                if *part == "=" && i + 1 < parts.len() {
                    if let Ok(val) = parts[i + 1].parse::<f64>() {
                        if line.contains("Side 1") {
                            final_result.xtr_upper = val;
                        } else if line.contains("Side 2") {
                            final_result.xtr_lower = val;
                        }
                    }
                }
            }
        }
    }

    write_json(&output_dir.join("final.json"), &final_result)?;
    Ok(())
}

fn parse_fortran_value(line: &str) -> f64 {
    // Parse Fortran E-notation values
    // Example: "CL =  1.234567890123456E-01"
    line.split('=')
        .nth(1)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0)
}

fn write_json<T: Serialize>(path: &Path, data: &T) -> Result<(), Box<dyn std::error::Error>> {
    let json = serde_json::to_string_pretty(data)?;
    fs::write(path, json)?;
    println!("  Wrote: {}", path.display());
    Ok(())
}
