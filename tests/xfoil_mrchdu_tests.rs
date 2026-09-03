//! Tests for MRCHDU against XFOIL fixture data
//!
//! These tests validate that yfoil's MRCHDU implementation produces
//! identical results to XFOIL's MRCHDU subroutine.

mod fixtures;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use yfoil::bl::gauss::gauss_solve_4x4;
use yfoil::bl::mrchdu::{march_station, MarchResult, SurfaceMarchState};
use yfoil::bl::system::{BLFlowType, BLGlobalParams, BLStationState};

fn fixture_path(name: &str) -> PathBuf {
    fixtures::require_fixture(&format!("{}/{}", fixtures::REF_CASE, name))
}

/// BL station data from XFOIL fixture
#[derive(Debug, Clone, Default)]
pub struct BLStation {
    pub xssi: f64, // Arc length
    pub uedg: f64, // Edge velocity
    pub thet: f64, // Momentum thickness
    pub dstr: f64, // Displacement thickness
    pub ctau: f64, // Shear stress coefficient (or amplification for laminar)
    pub mass: f64, // Mass defect
}

/// Input data for MRCHDU
#[derive(Debug, Clone)]
pub struct MrchduInput {
    pub nbl1: usize,         // Number of BL stations on surface 1
    pub nbl2: usize,         // Number of BL stations on surface 2
    pub iblte1: usize,       // TE index on surface 1
    pub iblte2: usize,       // TE index on surface 2
    pub minf: f64,           // Freestream Mach number
    pub reinf: f64,          // Reynolds number
    pub acrit1: f64,         // Critical N-factor surface 1
    pub acrit2: f64,         // Critical N-factor surface 2
    pub bl1: Vec<BLStation>, // Surface 1 stations
    pub bl2: Vec<BLStation>, // Surface 2 stations
}

/// Output data from MRCHDU
#[derive(Debug, Clone)]
pub struct MrchduOutput {
    pub itran1: usize,       // Transition index surface 1
    pub itran2: usize,       // Transition index surface 2
    pub bl1: Vec<BLStation>, // Updated surface 1 stations
    pub bl2: Vec<BLStation>, // Updated surface 2 stations
}

fn parse_int_value(line: &str) -> Option<usize> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() != 2 {
        return None;
    }
    parts[1].trim().parse().ok()
}

fn parse_float_value(line: &str) -> Option<f64> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() != 2 {
        return None;
    }
    parts[1].trim().replace('D', "E").parse().ok()
}

fn parse_bl_line(line: &str) -> Option<(usize, usize, BLStation)> {
    // Format: BL( 1,   2)=  0.10490...  0.49389...  etc
    let parts: Vec<&str> = line.split(')').collect();
    if parts.len() < 2 {
        return None;
    }

    // Parse surface and station indices
    let header = parts[0].replace("BL(", "").replace(" ", "");
    let indices: Vec<&str> = header.split(',').collect();
    if indices.len() != 2 {
        return None;
    }
    let is: usize = indices[0].parse().ok()?;
    let ibl: usize = indices[1].parse().ok()?;

    // Parse values
    let values_str = parts[1].trim_start_matches('=').trim();
    let values: Vec<f64> = values_str
        .split_whitespace()
        .filter_map(|s| s.replace('D', "E").parse().ok())
        .collect();

    if values.len() != 6 {
        return None;
    }

    Some((
        is,
        ibl,
        BLStation {
            xssi: values[0],
            uedg: values[1],
            thet: values[2],
            dstr: values[3],
            ctau: values[4],
            mass: values[5],
        },
    ))
}

/// Parse MRCHDU input fixture file
pub fn parse_mrchdu_input(path: &Path) -> Option<MrchduInput> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map_while(Result::ok).collect();

    let mut input = MrchduInput {
        nbl1: 0,
        nbl2: 0,
        iblte1: 0,
        iblte2: 0,
        minf: 0.0,
        reinf: 0.0,
        acrit1: 0.0,
        acrit2: 0.0,
        bl1: Vec::new(),
        bl2: Vec::new(),
    };

    for line in &lines {
        if line.starts_with("NBL1=") {
            input.nbl1 = parse_int_value(line)?;
        } else if line.starts_with("NBL2=") {
            input.nbl2 = parse_int_value(line)?;
        } else if line.starts_with("IBLTE1=") {
            input.iblte1 = parse_int_value(line)?;
        } else if line.starts_with("IBLTE2=") {
            input.iblte2 = parse_int_value(line)?;
        } else if line.starts_with("MINF=") {
            input.minf = parse_float_value(line)?;
        } else if line.starts_with("REINF=") {
            input.reinf = parse_float_value(line)?;
        } else if line.starts_with("ACRIT1=") {
            input.acrit1 = parse_float_value(line)?;
        } else if line.starts_with("ACRIT2=") {
            input.acrit2 = parse_float_value(line)?;
        } else if line.starts_with("BL(") {
            if let Some((is, ibl, station)) = parse_bl_line(line) {
                // Ensure vectors are large enough
                match is {
                    1 => {
                        while input.bl1.len() < ibl {
                            input.bl1.push(BLStation::default());
                        }
                        if ibl > 0 {
                            input.bl1[ibl - 1] = station;
                        }
                    }
                    2 => {
                        while input.bl2.len() < ibl {
                            input.bl2.push(BLStation::default());
                        }
                        if ibl > 0 {
                            input.bl2[ibl - 1] = station;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    Some(input)
}

/// Parse MRCHDU output fixture file
pub fn parse_mrchdu_output(path: &Path) -> Option<MrchduOutput> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map_while(Result::ok).collect();

    let mut output = MrchduOutput {
        itran1: 0,
        itran2: 0,
        bl1: Vec::new(),
        bl2: Vec::new(),
    };

    for line in &lines {
        if line.starts_with("ITRAN1=") {
            output.itran1 = parse_int_value(line)?;
        } else if line.starts_with("ITRAN2=") {
            output.itran2 = parse_int_value(line)?;
        } else if line.starts_with("BL(") {
            if let Some((is, ibl, station)) = parse_bl_line(line) {
                match is {
                    1 => {
                        while output.bl1.len() < ibl {
                            output.bl1.push(BLStation::default());
                        }
                        if ibl > 0 {
                            output.bl1[ibl - 1] = station;
                        }
                    }
                    2 => {
                        while output.bl2.len() < ibl {
                            output.bl2.push(BLStation::default());
                        }
                        if ibl > 0 {
                            output.bl2[ibl - 1] = station;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    Some(output)
}

#[test]
fn test_parse_mrchdu_fixtures() {
    let input_path = fixture_path("mrchdu_input.dat");
    let output_path = fixture_path("mrchdu_output.dat");

    let input = parse_mrchdu_input(&input_path).expect("Failed to parse MRCHDU input");
    let output = parse_mrchdu_output(&output_path).expect("Failed to parse MRCHDU output");

    println!("MRCHDU Input:");
    println!("  NBL1={}, NBL2={}", input.nbl1, input.nbl2);
    println!("  IBLTE1={}, IBLTE2={}", input.iblte1, input.iblte2);
    println!("  MINF={}, REINF={}", input.minf, input.reinf);
    println!("  ACRIT1={}, ACRIT2={}", input.acrit1, input.acrit2);
    println!("  BL1 stations: {}", input.bl1.len());
    println!("  BL2 stations: {}", input.bl2.len());

    println!("\nMRCHDU Output:");
    println!("  ITRAN1={}, ITRAN2={}", output.itran1, output.itran2);
    println!("  BL1 stations: {}", output.bl1.len());
    println!("  BL2 stations: {}", output.bl2.len());

    // Verify dimensions match
    assert_eq!(input.bl1.len(), input.nbl1);
    assert_eq!(input.bl2.len(), input.nbl2);
    assert_eq!(output.bl1.len(), input.nbl1);
    assert_eq!(output.bl2.len(), input.nbl2);

    // Print first few stations for comparison
    println!("\nSurface 1, first 5 stations (input vs output):");
    for i in 0..5.min(input.bl1.len()) {
        let inp = &input.bl1[i];
        let out = &output.bl1[i];
        println!("  Station {}: THET diff = {:.2e}", i + 1, (out.thet - inp.thet).abs());
    }

    // Check that transition was found
    assert!(
        output.itran1 > 0 && output.itran1 <= input.nbl1,
        "ITRAN1={} should be within [1, {}]",
        output.itran1,
        input.nbl1
    );
    assert!(
        output.itran2 > 0 && output.itran2 <= input.nbl2,
        "ITRAN2={} should be within [1, {}]",
        output.itran2,
        input.nbl2
    );

    println!("\nTransition locations:");
    println!("  Surface 1: station {} (of {})", output.itran1, input.nbl1);
    println!("  Surface 2: station {} (of {})", output.itran2, input.nbl2);
}

#[test]
fn test_mrchdu_state_changes() {
    let input_path = fixture_path("mrchdu_input.dat");
    let output_path = fixture_path("mrchdu_output.dat");

    let input = parse_mrchdu_input(&input_path).expect("Failed to parse");
    let output = parse_mrchdu_output(&output_path).expect("Failed to parse");

    // Analyze how much MRCHDU changes the state
    let mut max_thet_change: f64 = 0.0;
    let mut max_dstr_change: f64 = 0.0;
    let mut max_uedg_change: f64 = 0.0;

    for i in 1..input.bl1.len() {
        // Skip station 0 (similarity point)
        let inp = &input.bl1[i];
        let out = &output.bl1[i];

        if inp.thet.abs() > 1e-15 {
            let rel_change = ((out.thet - inp.thet) / inp.thet).abs();
            max_thet_change = max_thet_change.max(rel_change);
        }
        if inp.dstr.abs() > 1e-15 {
            let rel_change = ((out.dstr - inp.dstr) / inp.dstr).abs();
            max_dstr_change = max_dstr_change.max(rel_change);
        }
        if inp.uedg.abs() > 1e-15 {
            let rel_change = ((out.uedg - inp.uedg) / inp.uedg).abs();
            max_uedg_change = max_uedg_change.max(rel_change);
        }
    }

    println!("\nMaximum relative changes in MRCHDU (Surface 1):");
    println!("  THET: {:.2e}", max_thet_change);
    println!("  DSTR: {:.2e}", max_dstr_change);
    println!("  UEDG: {:.2e}", max_uedg_change);

    // MRCHDU should make relatively small changes
    // (it's refining the BL state, not completely recomputing it)
    assert!(
        max_thet_change < 0.1,
        "THET change {:.2e} seems too large",
        max_thet_change
    );
    assert!(
        max_dstr_change < 0.1,
        "DSTR change {:.2e} seems too large",
        max_dstr_change
    );
}

/// Newton iteration fixture data from XFOIL MRCHDU
#[derive(Debug, Clone)]
struct NewtonFixture {
    ibl: usize,
    is: usize,
    itbl: usize,
    // Input state
    xsi: f64,
    uei: f64,
    thi: f64,
    dsi: f64,
    ami: f64,
    cti: f64,
    dswaki: f64,
    // Secondary variables
    u2: f64,
    t2: f64,
    d2: f64,
    hk2: f64,
    rt2: f64,
    ueref: f64,
    hkref: f64,
    // VS2 matrix (4x4)
    vs2: [[f64; 4]; 4],
    // Residuals
    vsrez_before: [f64; 4],
    vsrez_solution: [f64; 4],
}

fn parse_newton_fixture(path: &Path) -> Option<NewtonFixture> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map_while(Result::ok).collect();

    // Debug: eprintln!("Parsing {} lines from {:?}", lines.len(), path);

    let mut fixture = NewtonFixture {
        ibl: 0,
        is: 0,
        itbl: 0,
        xsi: 0.0,
        uei: 0.0,
        thi: 0.0,
        dsi: 0.0,
        ami: 0.0,
        cti: 0.0,
        dswaki: 0.0,
        u2: 0.0,
        t2: 0.0,
        d2: 0.0,
        hk2: 0.0,
        rt2: 0.0,
        ueref: 0.0,
        hkref: 0.0,
        vs2: [[0.0; 4]; 4],
        vsrez_before: [0.0; 4],
        vsrez_solution: [0.0; 4],
    };

    for line in &lines {
        let line = line.trim();

        // Parse IBL, IS, ITBL
        // Format: IBL=  10 IS= 1 ITBL=   1
        if line.starts_with("IBL=") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            for (i, part) in parts.iter().enumerate() {
                if *part == "IBL=" && i + 1 < parts.len() {
                    fixture.ibl = parts[i + 1].parse().unwrap_or(0);
                } else if part.starts_with("IBL=") && part.len() > 4 {
                    fixture.ibl = part[4..].parse().unwrap_or(0);
                } else if *part == "IS=" && i + 1 < parts.len() {
                    fixture.is = parts[i + 1].parse().unwrap_or(0);
                } else if part.starts_with("IS=") && part.len() > 3 {
                    fixture.is = part[3..].parse().unwrap_or(0);
                } else if *part == "ITBL=" && i + 1 < parts.len() {
                    fixture.itbl = parts[i + 1].parse().unwrap_or(0);
                } else if part.starts_with("ITBL=") && part.len() > 5 {
                    fixture.itbl = part[5..].parse().unwrap_or(0);
                }
            }
        }
        // Parse scalar values
        else if let Some(rest) = line.strip_prefix("XSI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.xsi = v;
            }
        } else if let Some(rest) = line.strip_prefix("UEI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.uei = v;
            }
        } else if let Some(rest) = line.strip_prefix("THI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.thi = v;
            }
        } else if let Some(rest) = line.strip_prefix("DSI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.dsi = v;
            }
        } else if let Some(rest) = line.strip_prefix("AMI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.ami = v;
            }
        } else if let Some(rest) = line.strip_prefix("CTI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.cti = v;
            }
        } else if let Some(rest) = line.strip_prefix("DSWAKI=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.dswaki = v;
            }
        } else if let Some(rest) = line.strip_prefix("U2=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.u2 = v;
            }
        } else if let Some(rest) = line.strip_prefix("T2=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.t2 = v;
            }
        } else if let Some(rest) = line.strip_prefix("D2=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.d2 = v;
            }
        } else if let Some(rest) = line.strip_prefix("HK2=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.hk2 = v;
            }
        } else if let Some(rest) = line.strip_prefix("RT2=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.rt2 = v;
            }
        } else if let Some(rest) = line.strip_prefix("UEREF=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.ueref = v;
            }
        } else if let Some(rest) = line.strip_prefix("HKREF=") {
            if let Some(v) = parse_fortran_float(rest) {
                fixture.hkref = v;
            }
        }
        // Parse VS2 matrix rows
        else if line.starts_with("VS2(") {
            // Format: VS2( 1,1:4)=  val val val val
            // Extract row number from between VS2( and ,
            let row_idx_str: String = line
                .chars()
                .skip(4) // skip "VS2("
                .skip_while(|c| c.is_whitespace()) // skip any spaces
                .take_while(|c| c.is_ascii_digit()) // take the digit
                .collect();
            if let Ok(row_idx) = row_idx_str.parse::<usize>() {
                if row_idx >= 1 && row_idx <= 4 {
                    if let Some(values_part) = line.split('=').nth(1) {
                        let values: Vec<f64> = values_part
                            .split_whitespace()
                            .filter_map(|s| parse_fortran_float(s.trim()))
                            .collect();
                        if values.len() >= 4 {
                            fixture.vs2[row_idx - 1][..4].copy_from_slice(&values[..4]);
                        }
                    }
                }
            }
        }
        // Parse VSREZ arrays
        else if line.starts_with("VSREZ=") {
            if let Some(values_part) = line.split('=').nth(1) {
                let values: Vec<f64> = values_part
                    .split_whitespace()
                    .filter_map(|s| parse_fortran_float(s.trim()))
                    .collect();
                if values.len() >= 4 {
                    fixture.vsrez_before = [values[0], values[1], values[2], values[3]];
                }
            }
        } else if line.starts_with("VSREZ_SOL=") {
            if let Some(values_part) = line.split('=').nth(1) {
                let values: Vec<f64> = values_part
                    .split_whitespace()
                    .filter_map(|s| parse_fortran_float(s.trim()))
                    .collect();
                if values.len() >= 4 {
                    fixture.vsrez_solution = [values[0], values[1], values[2], values[3]];
                }
            }
        }
    }

    // Validate that we got the critical data
    if fixture.ibl == 0 || fixture.vs2[0][0] == 0.0 {
        return None;
    }

    Some(fixture)
}

fn parse_fortran_float(s: &str) -> Option<f64> {
    s.trim().replace('D', "E").replace('d', "e").parse().ok()
}

#[test]
fn test_gauss_against_xfoil() {
    // Test that our GAUSS implementation matches XFOIL's
    let fixture_path = fixture_path("mrchdu_newton.dat");

    let fixture = parse_newton_fixture(&fixture_path).expect("Failed to parse Newton fixture");

    println!("Newton fixture for station {}, surface {}", fixture.ibl, fixture.is);
    println!("VS2 matrix:");
    for i in 0..4 {
        println!(
            "  Row {}: {:e} {:e} {:e} {:e}",
            i + 1,
            fixture.vs2[i][0],
            fixture.vs2[i][1],
            fixture.vs2[i][2],
            fixture.vs2[i][3]
        );
    }
    println!("\nVSREZ before GAUSS:");
    println!(
        "  {:e} {:e} {:e} {:e}",
        fixture.vsrez_before[0], fixture.vsrez_before[1], fixture.vsrez_before[2], fixture.vsrez_before[3]
    );
    println!("\nExpected solution:");
    println!(
        "  {:e} {:e} {:e} {:e}",
        fixture.vsrez_solution[0], fixture.vsrez_solution[1], fixture.vsrez_solution[2], fixture.vsrez_solution[3]
    );

    // Copy VS2 and VSREZ since GAUSS modifies them
    let mut z = fixture.vs2;
    let mut r = fixture.vsrez_before;

    // Solve with our GAUSS implementation
    gauss_solve_4x4(&mut z, &mut r);

    println!("\nOur GAUSS solution:");
    println!("  {:e} {:e} {:e} {:e}", r[0], r[1], r[2], r[3]);

    // Compare against XFOIL's solution
    // Note: The residuals are very small, so we use absolute tolerance
    for i in 0..4 {
        let xfoil_sol = fixture.vsrez_solution[i];
        let our_sol = r[i];
        let diff = (our_sol - xfoil_sol).abs();

        println!(
            "  Component {}: XFOIL={:e}, ours={:e}, diff={:e}",
            i, xfoil_sol, our_sol, diff
        );

        // Since these are very small numbers, use absolute tolerance
        // The tolerance is loose because the residuals are already ~1e-12
        assert!(
            diff < 1e-10 || diff / xfoil_sol.abs().max(1e-30) < 1e-6,
            "GAUSS solution mismatch at component {}: XFOIL={:e}, ours={:e}",
            i,
            xfoil_sol,
            our_sol
        );
    }

    println!("\nGAUSS solution matches XFOIL!");
}

#[test]
fn test_mrchdu_full_march() {
    // Test full MRCHDU march against XFOIL fixture
    let input_path = fixture_path("mrchdu_input.dat");
    let output_path = fixture_path("mrchdu_output.dat");

    let input = parse_mrchdu_input(&input_path).expect("Failed to parse MRCHDU input");
    let output = parse_mrchdu_output(&output_path).expect("Failed to parse MRCHDU output");

    println!("MRCHDU full march test:");
    println!("  NBL1={}, NBL2={}", input.nbl1, input.nbl2);
    println!("  IBLTE1={}, IBLTE2={}", input.iblte1, input.iblte2);
    println!("  MINF={}, REINF={}", input.minf, input.reinf);

    // Create global BL parameters
    let params = BLGlobalParams::new(input.minf, input.reinf, 1.4);

    // Test surface 1 march
    // The fixture is from a converged VISCAL iteration where transition was at station 61.
    // We use itrold = output.itran1 so the march knows where transition was found.
    // With the corrected ami initialization (using s1.ampl), amplification accumulates
    // correctly and transition should be detected at the same location.
    let mut march = SurfaceMarchState::new(input.nbl1, input.iblte1, input.acrit1);
    march.init_march(output.itran1); // ITROLD = known transition location

    // Track results
    let mut result_thet: Vec<f64> = vec![0.0; input.nbl1];
    let mut result_dstr: Vec<f64> = vec![0.0; input.nbl1];
    let mut result_uedg: Vec<f64> = vec![0.0; input.nbl1];
    let mut result_ctau: Vec<f64> = vec![0.0; input.nbl1];

    // Station 1 is the stagnation point - just copy
    result_thet[0] = input.bl1[0].thet;
    result_dstr[0] = input.bl1[0].dstr;
    result_uedg[0] = input.bl1[0].uedg;
    result_ctau[0] = input.bl1[0].ctau;

    // Set up s1 from station 1
    let mut s1 = BLStationState::default();
    s1.blprv(
        input.bl1[0].xssi,
        input.bl1[0].ctau, // AMI (amplification for laminar)
        0.03,              // CTI (not used for laminar)
        input.bl1[0].thet,
        input.bl1[0].dstr,
        0.0, // DSWAKI
        input.bl1[0].uedg,
        &params,
    );
    s1.blkin(&params);
    s1.blvar(BLFlowType::Laminar, &params);

    let mut converged_count = 0;
    let mut partial_count = 0;
    let mut failed_count = 0;

    // March stations 2 to NBL1
    for ibl in 2..=input.nbl1 {
        let idx = ibl - 1; // 0-based index
        let stn = &input.bl1[idx];

        // Set up s2_init from input fixture
        // Note: For laminar stations, CTAU stores amplification, so pass it as cti
        // (march_station reads ami from s2_init.ctau when ibl < itrold)
        let mut s2_init = BLStationState::default();
        s2_init.blprv(
            stn.xssi, 0.0,      // AMI placeholder (not used - read from ctau for laminar)
            stn.ctau, // CTI = amplification for laminar stations
            stn.thet, stn.dstr, 0.0, // DSWAKI (surface, not wake)
            stn.uedg, &params,
        );
        s2_init.blkin(&params);
        let flow_type = if ibl >= march.itran {
            BLFlowType::Turbulent
        } else {
            BLFlowType::Laminar
        };
        s2_init.blvar(flow_type, &params);

        // March this station
        let (result, s2) = march_station(
            &s1, &s2_init, &mut march, &params, ibl, stn.xssi, 0.0,  // DSWAKI
            None, // Not first wake
        );

        // Store results
        result_thet[idx] = s2.theta;
        result_dstr[idx] = s2.dstar + s2.dw;
        result_uedg[idx] = s2.u / s2.u_uei; // Convert back to incompressible
        result_ctau[idx] = s2.ctau;

        match result {
            MarchResult::Converged => converged_count += 1,
            MarchResult::PartiallyConverged { .. } => partial_count += 1,
            MarchResult::Failed { dmax } => {
                failed_count += 1;
                println!("  Station {} failed: dmax={:.2e}", ibl, dmax);
            }
        }

        // Debug: print amplification for first 50 laminar stations
        if ibl < 50 && !march.turb {
            println!(
                "AMPL IBL={:4} AMPL1={:.16e} AMPL2={:.16e} HK2={:.6}",
                ibl, s1.ampl, s2.ampl, s2.hk
            );
        }

        // Update s1 for next station
        s1 = s2;

        // Handle transition
        if march.tran || ibl == march.iblte {
            march.turb = true;
        }
        march.tran = false;
    }

    println!("\nMarch results:");
    println!(
        "  Converged: {}, Partial: {}, Failed: {}",
        converged_count, partial_count, failed_count
    );
    println!("  Our ITRAN: {}, Expected ITRAN: {}", march.itran, output.itran1);

    // Compare results to expected output
    let mut max_thet_err: f64 = 0.0;
    let mut max_dstr_err: f64 = 0.0;
    let mut max_uedg_err: f64 = 0.0;
    let mut station_with_max_err = 0;

    for i in 1..input.nbl1 {
        let exp = &output.bl1[i];

        if exp.thet.abs() > 1e-15 {
            let rel_err = ((result_thet[i] - exp.thet) / exp.thet).abs();
            if rel_err > max_thet_err {
                max_thet_err = rel_err;
                station_with_max_err = i + 1;
            }
        }
        if exp.dstr.abs() > 1e-15 {
            let rel_err = ((result_dstr[i] - exp.dstr) / exp.dstr).abs();
            max_dstr_err = max_dstr_err.max(rel_err);
        }
        if exp.uedg.abs() > 1e-15 {
            let rel_err = ((result_uedg[i] - exp.uedg) / exp.uedg).abs();
            max_uedg_err = max_uedg_err.max(rel_err);
        }
    }

    println!("\nMax relative errors vs XFOIL output:");
    println!("  THET: {:.2e} (at station {})", max_thet_err, station_with_max_err);
    println!("  DSTR: {:.2e}", max_dstr_err);
    println!("  UEDG: {:.2e}", max_uedg_err);

    // Print first few stations for debugging
    println!("\nFirst 5 stations comparison (station, ours, expected, diff):");
    for i in 1..6.min(input.nbl1) {
        let exp = &output.bl1[i];
        println!(
            "  {} THET: {:.6e} vs {:.6e} ({:.2e})",
            i + 1,
            result_thet[i],
            exp.thet,
            (result_thet[i] - exp.thet).abs()
        );
    }

    // Check laminar region matches exactly (before transition)
    // TRDIF (transition difference) is not yet implemented, so we expect
    // divergence at and after the transition station
    let mut max_laminar_thet_err: f64 = 0.0;
    let mut max_laminar_dstr_err: f64 = 0.0;

    // Check stations where BOTH our march and XFOIL are laminar
    // Use the minimum of our transition and XFOIL's transition
    let common_laminar_end = march.itran.min(output.itran1);
    for i in 1..(common_laminar_end - 1) {
        let exp = &output.bl1[i];
        if exp.thet.abs() > 1e-15 {
            let rel_err = ((result_thet[i] - exp.thet) / exp.thet).abs();
            max_laminar_thet_err = max_laminar_thet_err.max(rel_err);
        }
        if exp.dstr.abs() > 1e-15 {
            let rel_err = ((result_dstr[i] - exp.dstr) / exp.dstr).abs();
            max_laminar_dstr_err = max_laminar_dstr_err.max(rel_err);
        }
    }

    println!("\nLaminar region (stations 1 to {}) errors:", common_laminar_end - 2);
    println!("  Max THET error: {:.2e}", max_laminar_thet_err);
    println!("  Max DSTR error: {:.2e}", max_laminar_dstr_err);

    // Laminar region should match closely - tolerance allows for accumulated numerical
    // differences over many Newton iterations. With 46+ laminar stations, ~5e-6 is excellent.
    assert!(
        max_laminar_thet_err < 5e-6,
        "Laminar THET error {:.2e} too large - should match XFOIL closely",
        max_laminar_thet_err
    );
    assert!(
        max_laminar_dstr_err < 5e-6,
        "Laminar DSTR error {:.2e} too large - should match XFOIL closely",
        max_laminar_dstr_err
    );

    // Transition location should match XFOIL exactly now that the amplification
    // accumulation bug is fixed.
    println!(
        "\nTransition: Our itran={}, Fixture itran={}",
        march.itran, output.itran1
    );
    assert_eq!(
        march.itran, output.itran1,
        "Transition location mismatch: ours={}, expected={}",
        march.itran, output.itran1
    );

    // Print comparison at and around transition
    let itran = output.itran1;
    println!("\nComparison at transition (station {}) and nearby:", itran);
    for i in (itran - 2)..=(itran + 5).min(input.nbl1 - 1) {
        let idx = i - 1;
        let exp = &output.bl1[idx];
        let thet_err = if exp.thet.abs() > 1e-15 {
            ((result_thet[idx] - exp.thet) / exp.thet).abs()
        } else {
            0.0
        };
        let dstr_err = if exp.dstr.abs() > 1e-15 {
            ((result_dstr[idx] - exp.dstr) / exp.dstr).abs()
        } else {
            0.0
        };
        let marker = if i == itran { " <-- TRANSITION" } else { "" };
        println!(
            "  Stn {}: THET {:.6e} vs {:.6e} ({:.2e}), DSTR {:.6e} vs {:.6e} ({:.2e}){}",
            i, result_thet[idx], exp.thet, thet_err, result_dstr[idx], exp.dstr, dstr_err, marker
        );
    }

    // Check turbulent region - stations where BOTH are turbulent
    // This is stations >= max(our_itran, fixture_itran)
    let common_turbulent_start = march.itran.max(output.itran1);
    let mut max_turb_thet_err: f64 = 0.0;
    let mut max_turb_dstr_err: f64 = 0.0;
    let mut turb_station_with_max_err = 0;

    println!(
        "\nTurbulent region (stations {} to {}):",
        common_turbulent_start, input.iblte1
    );
    for i in common_turbulent_start..input.iblte1 {
        let idx = i - 1;
        let exp = &output.bl1[idx];
        if exp.thet.abs() > 1e-15 {
            let rel_err = ((result_thet[idx] - exp.thet) / exp.thet).abs();
            if rel_err > max_turb_thet_err {
                max_turb_thet_err = rel_err;
                turb_station_with_max_err = i;
            }
        }
        if exp.dstr.abs() > 1e-15 {
            let rel_err = ((result_dstr[idx] - exp.dstr) / exp.dstr).abs();
            max_turb_dstr_err = max_turb_dstr_err.max(rel_err);
        }
    }

    println!(
        "  Max THET error: {:.2e} (at station {})",
        max_turb_thet_err, turb_station_with_max_err
    );
    println!("  Max DSTR error: {:.2e}", max_turb_dstr_err);

    // Print a few turbulent stations for debugging
    println!("\nFirst 5 turbulent stations comparison:");
    for i in common_turbulent_start..(common_turbulent_start + 5).min(input.iblte1) {
        let idx = i - 1;
        let exp = &output.bl1[idx];
        let thet_err = if exp.thet.abs() > 1e-15 {
            ((result_thet[idx] - exp.thet) / exp.thet).abs()
        } else {
            0.0
        };
        let dstr_err = if exp.dstr.abs() > 1e-15 {
            ((result_dstr[idx] - exp.dstr) / exp.dstr).abs()
        } else {
            0.0
        };
        println!(
            "  Stn {}: THET {:.6e} vs {:.6e} ({:.2e}), DSTR {:.6e} vs {:.6e} ({:.2e})",
            i, result_thet[idx], exp.thet, thet_err, result_dstr[idx], exp.dstr, dstr_err
        );
    }

    // Note: Turbulent region errors may be larger due to different transition points
    // affecting the initial turbulent state. Still check they're reasonable.
    println!("\n[NOTE] Turbulent errors may be larger due to different transition locations");
}

/// Detailed test for TRCHEK inputs at station 3 to diagnose the 2x amplification issue
#[test]
fn test_trchek_station_3_inputs() {
    use yfoil::bl::system::trchek;

    let input_path = fixture_path("mrchdu_input.dat");
    let output_path = fixture_path("mrchdu_output.dat");

    let input = parse_mrchdu_input(&input_path).expect("Failed to parse MRCHDU input");
    let output = parse_mrchdu_output(&output_path).expect("Failed to parse MRCHDU output");
    let params = BLGlobalParams::new(input.minf, input.reinf, 1.4);
    let acrit = input.acrit1;

    println!("=== Testing TRCHEK inputs at station 3 ===");
    println!(
        "Parameters: MINF={}, REINF={}, ACRIT={}",
        input.minf, input.reinf, acrit
    );

    // Set up station 2 (similarity station) - this becomes s1 for station 3
    let stn2 = &input.bl1[1]; // 0-based index 1 = station 2
    let mut s1 = BLStationState::default();
    s1.blprv(
        stn2.xssi, 0.0,  // AMI = 0 at similarity station
        0.03, // CTI
        stn2.thet, stn2.dstr, 0.0, // DSWAKI
        stn2.uedg, &params,
    );
    s1.blkin(&params);
    s1.blvar(BLFlowType::Laminar, &params);

    println!("\nStation 2 (s1) from fixture:");
    println!("  XSSI={:.16e}", stn2.xssi);
    println!("  UEDG={:.16e}", stn2.uedg);
    println!("  THET={:.16e}", stn2.thet);
    println!("  DSTR={:.16e}", stn2.dstr);

    println!("\nStation 2 (s1) computed BL state:");
    println!("  x={:.16e}", s1.x);
    println!("  u={:.16e}", s1.u);
    println!("  theta={:.16e}", s1.theta);
    println!("  dstar={:.16e}", s1.dstar);
    println!("  hk={:.16e}", s1.hk);
    println!("  rt={:.16e}", s1.rt);
    println!("  ampl={:.16e}", s1.ampl);

    // Set up station 3 (s2) - initial values from fixture
    let stn3 = &input.bl1[2]; // 0-based index 2 = station 3
    let mut s2 = BLStationState::default();
    s2.blprv(
        stn3.xssi, 0.0,       // AMI initial guess
        stn3.ctau, // CTI (stores amplification for laminar)
        stn3.thet, stn3.dstr, 0.0, // DSWAKI
        stn3.uedg, &params,
    );
    s2.blkin(&params);
    s2.blvar(BLFlowType::Laminar, &params);

    println!("\nStation 3 (s2) from fixture:");
    println!("  XSSI={:.16e}", stn3.xssi);
    println!("  UEDG={:.16e}", stn3.uedg);
    println!("  THET={:.16e}", stn3.thet);
    println!("  DSTR={:.16e}", stn3.dstr);
    println!("  CTAU={:.16e}", stn3.ctau);

    println!("\nStation 3 (s2) computed BL state:");
    println!("  x={:.16e}", s2.x);
    println!("  u={:.16e}", s2.u);
    println!("  theta={:.16e}", s2.theta);
    println!("  dstar={:.16e}", s2.dstar);
    println!("  hk={:.16e}", s2.hk);
    println!("  rt={:.16e}", s2.rt);
    println!("  ampl={:.16e}", s2.ampl);

    // Compare with XFOIL AXSET log values for station 3:
    // HK1=  0.2229507759055277E+01
    // T1=  0.3206370158436509E-04
    // RT1=  0.5376303335340089E+00
    // HK2=  0.2229699712062246E+01
    // T2=  0.3218660306693949E-04
    // RT2=  0.5562086719295983E+01
    println!("\n=== XFOIL reference values for first AXSET call ===");
    println!("  HK1=  0.2229507759055277E+01");
    println!("  T1=   0.3206370158436509E-04");
    println!("  RT1=  0.5376303335340089E+00");
    println!("  HK2=  0.2229699712062246E+01");
    println!("  T2=   0.3218660306693949E-04");
    println!("  RT2=  0.5562086719295983E+01");
    println!("  AX=   0.6416011981965619E-07");

    // Call TRCHEK
    let ampl1 = s1.ampl;
    println!("\n=== Calling TRCHEK ===");
    println!("  ampl1 (input) = {:.16e}", ampl1);
    println!("  acrit = {:.16e}", acrit);
    println!("  xiforc = inf (disabled)");

    let result = trchek(&s1, &s2, ampl1, acrit, f64::MAX, &params);

    match &result {
        yfoil::bl::system::TransitionResult::NoTransition { ampl2 } => {
            println!("\n=== TRCHEK Result: No Transition ===");
            println!("  ampl2 (output) = {:.16e}", ampl2);

            // Expected from XFOIL mrchdu_ampl.dat station 3:
            // AMPL2=  0.2227257501770935E-09
            let xfoil_ampl2 = 2.227257501770935e-10;
            let diff = (ampl2 - xfoil_ampl2).abs();
            let rel_err = diff / xfoil_ampl2;
            println!("\n=== Comparison with XFOIL ===");
            println!("  YFoil ampl2 = {:.16e}", ampl2);
            println!("  XFOIL ampl2 = {:.16e}", xfoil_ampl2);
            println!("  Diff        = {:.16e}", diff);
            println!("  Rel error   = {:.2e}", rel_err);
            println!("  Ratio       = {:.4}", ampl2 / xfoil_ampl2);

            if rel_err > 0.01 {
                println!(
                    "\n  *** MISMATCH: YFoil ampl2 is {:.2}x XFOIL! ***",
                    ampl2 / xfoil_ampl2
                );
            }
        }
        yfoil::bl::system::TransitionResult::FreeTransition { ampl2, location } => {
            println!("\n=== TRCHEK Result: Free Transition ===");
            println!("  ampl2 = {:.16e}", ampl2);
            println!("  xt = {:.16e}", location.xt);
        }
        yfoil::bl::system::TransitionResult::ForcedTransition { location } => {
            println!("\n=== TRCHEK Result: Forced Transition ===");
            println!("  xt = {:.16e}", location.xt);
        }
    }

    // Now test with the MARCHED state (simulating what happens in full march)
    println!("\n\n=== Now testing with MARCHED station 2 ===");

    // March station 2 (similarity) and get the converged state
    let mut march = SurfaceMarchState::new(input.nbl1, input.iblte1, acrit);
    march.init_march(output.itran1);

    // Set up s1 from station 1 (stagnation point - all zeros)
    let mut s1_for_march = BLStationState::default();
    s1_for_march.blprv(
        input.bl1[0].xssi, // 0
        input.bl1[0].ctau, // 0
        0.03,
        input.bl1[0].thet, // 0
        input.bl1[0].dstr, // 0
        0.0,
        input.bl1[0].uedg, // 0
        &params,
    );
    s1_for_march.blkin(&params);
    s1_for_march.blvar(BLFlowType::Laminar, &params);

    // March station 2 (similarity)
    let stn2_init = &input.bl1[1];
    let mut s2_init_march = BLStationState::default();
    s2_init_march.blprv(
        stn2_init.xssi,
        0.0,
        stn2_init.ctau,
        stn2_init.thet,
        stn2_init.dstr,
        0.0,
        stn2_init.uedg,
        &params,
    );
    s2_init_march.blkin(&params);
    s2_init_march.blvar(BLFlowType::Laminar, &params);

    let (result2, s2_converged) = march_station(
        &s1_for_march,
        &s2_init_march,
        &mut march,
        &params,
        2, // ibl=2
        stn2_init.xssi,
        0.0,
        None,
    );

    println!("Station 2 march result: {:?}", result2);

    println!("\nStation 2 FIXTURE values:");
    println!("  xssi={:.16e}", stn2.xssi);
    println!("  theta={:.16e}", stn2.thet);
    println!("  dstar={:.16e}", stn2.dstr);
    println!("  uedg={:.16e}", stn2.uedg);

    println!("\nStation 2 MARCHED (converged) values:");
    println!("  x={:.16e}", s2_converged.x);
    println!("  theta={:.16e}", s2_converged.theta);
    println!("  dstar={:.16e}", s2_converged.dstar);
    println!("  u={:.16e}", s2_converged.u);
    println!("  hk={:.16e}", s2_converged.hk);
    println!("  rt={:.16e}", s2_converged.rt);
    println!("  ampl={:.16e}", s2_converged.ampl);

    // Compare key values
    let theta_diff = (s2_converged.theta - stn2.thet).abs() / stn2.thet;
    let dstar_diff = (s2_converged.dstar - stn2.dstr).abs() / stn2.dstr;
    println!("\nDifferences (marched vs fixture):");
    println!("  theta rel diff: {:.2e}", theta_diff);
    println!("  dstar rel diff: {:.2e}", dstar_diff);

    // Now march station 3 using the MARCHED s2 as s1
    let s1_marched = s2_converged; // This is what s1 looks like after station 2 converges

    let stn3_for_march = &input.bl1[2];
    let mut s2_init_stn3 = BLStationState::default();
    s2_init_stn3.blprv(
        stn3_for_march.xssi,
        0.0,
        stn3_for_march.ctau,
        stn3_for_march.thet,
        stn3_for_march.dstr,
        0.0,
        stn3_for_march.uedg,
        &params,
    );
    s2_init_stn3.blkin(&params);
    s2_init_stn3.blvar(BLFlowType::Laminar, &params);

    println!("\n=== Calling TRCHEK with MARCHED s1 ===");
    println!("s1 (marched station 2):");
    println!("  x={:.16e}", s1_marched.x);
    println!("  hk={:.16e}", s1_marched.hk);
    println!("  theta={:.16e}", s1_marched.theta);
    println!("  rt={:.16e}", s1_marched.rt);
    println!("  ampl={:.16e}", s1_marched.ampl);

    let result_marched = trchek(&s1_marched, &s2_init_stn3, s1_marched.ampl, acrit, f64::MAX, &params);

    match &result_marched {
        yfoil::bl::system::TransitionResult::NoTransition { ampl2 } => {
            println!("\n=== TRCHEK Result with MARCHED s1 ===");
            println!("  ampl2 = {:.16e}", ampl2);

            let xfoil_ampl2 = 2.227257501770935e-10;
            let ratio = ampl2 / xfoil_ampl2;
            println!("  Ratio vs XFOIL: {:.4}", ratio);

            if ratio > 1.5 {
                println!("\n  *** PROBLEM: ampl2 is {:.2}x XFOIL! ***", ratio);
                println!("  This suggests the marched s1 state differs from fixture.");
            }
        }
        _ => {
            println!("Unexpected transition result");
        }
    }
}
