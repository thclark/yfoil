//! Generate subroutine-level validation fixtures from instrumented XFOIL
//!
//! This example:
//! 1. Runs the instrumented XFOIL binary on test cases
//! 2. Parses the xfoil_subroutine_log.dat output
//! 3. Extracts unique input/output pairs for each subroutine
//! 4. Writes JSON fixtures for testing
//!
//! Run with: cargo run --example generate_subroutine_fixtures

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

/// HKIN fixture: Kinematic shape factor
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HkinFixture {
    pub input: HkinInput,
    pub output: HkinOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HkinInput {
    pub h: f64,
    pub msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HkinOutput {
    pub hk: f64,
    pub hk_h: f64,
    pub hk_msq: f64,
}

/// CFL fixture: Laminar skin friction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CflFixture {
    pub input: CflInput,
    pub output: CflOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CflInput {
    pub hk: f64,
    pub rt: f64,
    pub msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CflOutput {
    pub cf: f64,
    pub cf_hk: f64,
    pub cf_rt: f64,
    pub cf_msq: f64,
}

/// HSL fixture: Laminar energy shape factor
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HslFixture {
    pub input: HslInput,
    pub output: HslOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HslInput {
    pub hk: f64,
    pub rt: f64,
    pub msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HslOutput {
    pub hs: f64,
    pub hs_hk: f64,
    pub hs_rt: f64,
    pub hs_msq: f64,
}

/// DIL fixture: Laminar dissipation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DilFixture {
    pub input: DilInput,
    pub output: DilOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DilInput {
    pub hk: f64,
    pub rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DilOutput {
    pub di: f64,
    pub di_hk: f64,
    pub di_rt: f64,
}

/// HST fixture: Turbulent energy shape factor
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HstFixture {
    pub input: HstInput,
    pub output: HstOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HstInput {
    pub hk: f64,
    pub rt: f64,
    pub msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HstOutput {
    pub hs: f64,
    pub hs_hk: f64,
    pub hs_rt: f64,
    pub hs_msq: f64,
}

/// CFT fixture: Turbulent skin friction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CftFixture {
    pub input: CftInput,
    pub output: CftOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CftInput {
    pub hk: f64,
    pub rt: f64,
    pub msq: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CftOutput {
    pub cf: f64,
    pub cf_hk: f64,
    pub cf_rt: f64,
    pub cf_msq: f64,
}

/// DAMPL fixture: Amplification rate
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DamplFixture {
    pub input: DamplInput,
    pub output: DamplOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DamplInput {
    pub hk: f64,
    pub th: f64,
    pub rt: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DamplOutput {
    pub ax: f64,
    pub ax_hk: f64,
    pub ax_th: f64,
    pub ax_rt: f64,
}

/// BLKIN fixture: Secondary BL variables
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlkinFixture {
    pub input: BlkinInput,
    pub output: BlkinOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BlkinInput {
    pub t2: f64,
    pub d2: f64,
    pub u2: f64,
    pub hstinv: f64,
    pub gm1bl: f64,
    pub rstbl: f64,
    pub hvrat: f64,
    pub reybl: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlkinOutput {
    pub m2: f64,
    pub h2: f64,
    pub hk2: f64,
    pub rt2: f64,
    pub hk2_t2: f64,
    pub hk2_d2: f64,
    pub hk2_u2: f64,
    pub rt2_t2: f64,
    pub rt2_u2: f64,
}

/// BLVAR fixture: Closure calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlvarFixture {
    pub input: BlvarInput,
    pub output: BlvarOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BlvarInput {
    pub ityp: i32,
    pub hk2: f64,
    pub rt2: f64,
    pub m2: f64,
    pub t2: f64,
    pub d2: f64,
    pub s2: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlvarOutput {
    pub hs2: f64,
    pub cf2: f64,
    pub di2: f64,
    pub us2: f64,
    pub cq2: f64,
    pub de2: f64,
}

/// TRCHEK2 fixture: Transition check
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trchek2Fixture {
    pub input: Trchek2Input,
    pub output: Trchek2Output,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Trchek2Input {
    pub x1: f64,
    pub x2: f64,
    pub ampl1: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trchek2Output {
    pub ampl2: f64,
    pub xt: f64,
    pub tran: bool,
    pub turb: bool,
}

/// BLSYS fixture: Newton system assembly
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlsysFixture {
    pub output: BlsysOutput,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlsysOutput {
    pub vsrez: [f64; 4],
    pub vs1: [[f64; 5]; 3],
    pub vs2: [[f64; 5]; 3],
}

/// Parse a floating point value from a line like "HK=  0.2200000000000000E+01"
fn parse_value(line: &str) -> Option<f64> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() == 2 {
        parts[1].trim().parse::<f64>().ok()
    } else {
        None
    }
}

/// Parse a boolean from format like "TRAN=T" or "TRAN=F"
fn parse_bool(line: &str) -> Option<bool> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() == 2 {
        let val = parts[1].trim();
        if val == "T" {
            Some(true)
        } else if val == "F" {
            Some(false)
        } else {
            None
        }
    } else {
        None
    }
}

/// Represents which subroutine we're currently parsing
#[derive(Debug, Clone, Copy, PartialEq)]
enum CurrentSub {
    None,
    Hkin,
    Dil,
    Hsl,
    Cfl,
    Hst,
    Cft,
    Dampl,
    Blkin,
    Blvar,
    Trchek2,
    Blsys,
}

#[derive(Default)]
struct ParsedFixtures {
    hkin: Vec<HkinFixture>,
    cfl: Vec<CflFixture>,
    hsl: Vec<HslFixture>,
    dil: Vec<DilFixture>,
    hst: Vec<HstFixture>,
    cft: Vec<CftFixture>,
    dampl: Vec<DamplFixture>,
    blkin: Vec<BlkinFixture>,
    blvar: Vec<BlvarFixture>,
    trchek2: Vec<Trchek2Fixture>,
    #[allow(dead_code)]
    blsys: Vec<BlsysFixture>,
}

/// Parse the XFOIL subroutine log and extract fixtures
fn parse_log(log_path: &Path) -> Result<ParsedFixtures, Box<dyn std::error::Error>> {
    let file = File::open(log_path)?;
    let reader = BufReader::new(file);

    let mut fixtures = ParsedFixtures::default();
    let mut current_sub = CurrentSub::None;
    let mut current_values: HashMap<String, f64> = HashMap::new();
    let mut current_bools: HashMap<String, bool> = HashMap::new();

    // Track unique inputs to avoid duplicates
    let mut seen_hkin: HashSet<String> = HashSet::new();
    let mut seen_cfl: HashSet<String> = HashSet::new();
    let mut seen_hsl: HashSet<String> = HashSet::new();
    let mut seen_dil: HashSet<String> = HashSet::new();
    let mut seen_hst: HashSet<String> = HashSet::new();
    let mut seen_cft: HashSet<String> = HashSet::new();
    let mut seen_dampl: HashSet<String> = HashSet::new();
    let mut seen_blkin: HashSet<String> = HashSet::new();
    let mut seen_blvar: HashSet<String> = HashSet::new();
    let mut seen_trchek2: HashSet<String> = HashSet::new();
    let mut seen_blsys: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        // Check for subroutine markers
        if line.starts_with("=== ") && line.ends_with(" ===") {
            // Save previous subroutine data before switching
            save_fixture(
                &current_sub,
                &current_values,
                &current_bools,
                &mut fixtures,
                &mut seen_hkin,
                &mut seen_cfl,
                &mut seen_hsl,
                &mut seen_dil,
                &mut seen_hst,
                &mut seen_cft,
                &mut seen_dampl,
                &mut seen_blkin,
                &mut seen_blvar,
                &mut seen_trchek2,
                &mut seen_blsys,
            );

            current_values.clear();
            current_bools.clear();

            current_sub = match line {
                "=== HKIN ===" => CurrentSub::Hkin,
                "=== DIL ===" => CurrentSub::Dil,
                "=== HSL ===" => CurrentSub::Hsl,
                "=== CFL ===" => CurrentSub::Cfl,
                "=== HST ===" => CurrentSub::Hst,
                "=== CFT ===" => CurrentSub::Cft,
                "=== DAMPL ===" => CurrentSub::Dampl,
                "=== BLKIN ===" => CurrentSub::Blkin,
                "=== BLVAR ===" => CurrentSub::Blvar,
                "=== TRCHEK2 ===" => CurrentSub::Trchek2,
                "=== BLSYS ===" => CurrentSub::Blsys,
                _ => CurrentSub::None,
            };
            continue;
        }

        // Parse boolean values first (for TRCHEK2)
        if let Some(val) = parse_bool(line) {
            let key = line.split('=').next().unwrap().trim().to_uppercase();
            current_bools.insert(key, val);
        }
        // Parse value lines
        else if let Some(val) = parse_value(line) {
            let key = line.split('=').next().unwrap().trim().to_uppercase();
            current_values.insert(key, val);
        }
    }

    // Save final subroutine
    save_fixture(
        &current_sub,
        &current_values,
        &current_bools,
        &mut fixtures,
        &mut seen_hkin,
        &mut seen_cfl,
        &mut seen_hsl,
        &mut seen_dil,
        &mut seen_hst,
        &mut seen_cft,
        &mut seen_dampl,
        &mut seen_blkin,
        &mut seen_blvar,
        &mut seen_trchek2,
        &mut seen_blsys,
    );

    Ok(fixtures)
}

fn save_fixture(
    current_sub: &CurrentSub,
    values: &HashMap<String, f64>,
    bools: &HashMap<String, bool>,
    fixtures: &mut ParsedFixtures,
    seen_hkin: &mut HashSet<String>,
    seen_cfl: &mut HashSet<String>,
    seen_hsl: &mut HashSet<String>,
    seen_dil: &mut HashSet<String>,
    seen_hst: &mut HashSet<String>,
    seen_cft: &mut HashSet<String>,
    seen_dampl: &mut HashSet<String>,
    seen_blkin: &mut HashSet<String>,
    seen_blvar: &mut HashSet<String>,
    seen_trchek2: &mut HashSet<String>,
    seen_blsys: &mut HashSet<String>,
) {
    match current_sub {
        CurrentSub::Hkin => {
            if let (Some(&h), Some(&msq), Some(&hk), Some(&hk_h), Some(&hk_msq)) = (
                values.get("H"),
                values.get("MSQ"),
                values.get("HK"),
                values.get("HK_H"),
                values.get("HK_MSQ"),
            ) {
                let key = format!("{:.16e}_{:.16e}", h, msq);
                if !seen_hkin.contains(&key) && fixtures.hkin.len() < 200 {
                    seen_hkin.insert(key);
                    fixtures.hkin.push(HkinFixture {
                        input: HkinInput { h, msq },
                        output: HkinOutput { hk, hk_h, hk_msq },
                    });
                }
            }
        }
        CurrentSub::Cfl => {
            if let (Some(&hk), Some(&rt), Some(&msq), Some(&cf), Some(&cf_hk), Some(&cf_rt), Some(&cf_msq)) = (
                values.get("HK"),
                values.get("RT"),
                values.get("MSQ"),
                values.get("CF"),
                values.get("CF_HK"),
                values.get("CF_RT"),
                values.get("CF_MSQ"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", hk, rt, msq);
                if !seen_cfl.contains(&key) && fixtures.cfl.len() < 200 {
                    seen_cfl.insert(key);
                    fixtures.cfl.push(CflFixture {
                        input: CflInput { hk, rt, msq },
                        output: CflOutput {
                            cf,
                            cf_hk,
                            cf_rt,
                            cf_msq,
                        },
                    });
                }
            }
        }
        CurrentSub::Hsl => {
            if let (Some(&hk), Some(&rt), Some(&msq), Some(&hs), Some(&hs_hk), Some(&hs_rt), Some(&hs_msq)) = (
                values.get("HK"),
                values.get("RT"),
                values.get("MSQ"),
                values.get("HS"),
                values.get("HS_HK"),
                values.get("HS_RT"),
                values.get("HS_MSQ"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", hk, rt, msq);
                if !seen_hsl.contains(&key) && fixtures.hsl.len() < 200 {
                    seen_hsl.insert(key);
                    fixtures.hsl.push(HslFixture {
                        input: HslInput { hk, rt, msq },
                        output: HslOutput {
                            hs,
                            hs_hk,
                            hs_rt,
                            hs_msq,
                        },
                    });
                }
            }
        }
        CurrentSub::Dil => {
            if let (Some(&hk), Some(&rt), Some(&di), Some(&di_hk), Some(&di_rt)) = (
                values.get("HK"),
                values.get("RT"),
                values.get("DI"),
                values.get("DI_HK"),
                values.get("DI_RT"),
            ) {
                let key = format!("{:.16e}_{:.16e}", hk, rt);
                if !seen_dil.contains(&key) && fixtures.dil.len() < 200 {
                    seen_dil.insert(key);
                    fixtures.dil.push(DilFixture {
                        input: DilInput { hk, rt },
                        output: DilOutput { di, di_hk, di_rt },
                    });
                }
            }
        }
        CurrentSub::Hst => {
            if let (Some(&hk), Some(&rt), Some(&msq), Some(&hs), Some(&hs_hk), Some(&hs_rt), Some(&hs_msq)) = (
                values.get("HK"),
                values.get("RT"),
                values.get("MSQ"),
                values.get("HS"),
                values.get("HS_HK"),
                values.get("HS_RT"),
                values.get("HS_MSQ"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", hk, rt, msq);
                if !seen_hst.contains(&key) && fixtures.hst.len() < 200 {
                    seen_hst.insert(key);
                    fixtures.hst.push(HstFixture {
                        input: HstInput { hk, rt, msq },
                        output: HstOutput {
                            hs,
                            hs_hk,
                            hs_rt,
                            hs_msq,
                        },
                    });
                }
            }
        }
        CurrentSub::Cft => {
            if let (Some(&hk), Some(&rt), Some(&msq), Some(&cf), Some(&cf_hk), Some(&cf_rt), Some(&cf_msq)) = (
                values.get("HK"),
                values.get("RT"),
                values.get("MSQ"),
                values.get("CF"),
                values.get("CF_HK"),
                values.get("CF_RT"),
                values.get("CF_MSQ"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", hk, rt, msq);
                if !seen_cft.contains(&key) && fixtures.cft.len() < 200 {
                    seen_cft.insert(key);
                    fixtures.cft.push(CftFixture {
                        input: CftInput { hk, rt, msq },
                        output: CftOutput {
                            cf,
                            cf_hk,
                            cf_rt,
                            cf_msq,
                        },
                    });
                }
            }
        }
        CurrentSub::Dampl => {
            if let (Some(&hk), Some(&th), Some(&rt), Some(&ax), Some(&ax_hk), Some(&ax_th), Some(&ax_rt)) = (
                values.get("HK"),
                values.get("TH"),
                values.get("RT"),
                values.get("AX"),
                values.get("AX_HK"),
                values.get("AX_TH"),
                values.get("AX_RT"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", hk, th, rt);
                if !seen_dampl.contains(&key) && fixtures.dampl.len() < 200 {
                    seen_dampl.insert(key);
                    fixtures.dampl.push(DamplFixture {
                        input: DamplInput { hk, th, rt },
                        output: DamplOutput {
                            ax,
                            ax_hk,
                            ax_th,
                            ax_rt,
                        },
                    });
                }
            }
        }
        CurrentSub::Blkin => {
            if let (
                Some(&t2),
                Some(&d2),
                Some(&u2),
                Some(&hstinv),
                Some(&gm1bl),
                Some(&rstbl),
                Some(&hvrat),
                Some(&reybl),
                Some(&m2),
                Some(&h2),
                Some(&hk2),
                Some(&rt2),
                Some(&hk2_t2),
                Some(&hk2_d2),
                Some(&hk2_u2),
                Some(&rt2_t2),
                Some(&rt2_u2),
            ) = (
                values.get("T2"),
                values.get("D2"),
                values.get("U2"),
                values.get("HSTINV"),
                values.get("GM1BL"),
                values.get("RSTBL"),
                values.get("HVRAT"),
                values.get("REYBL"),
                values.get("M2"),
                values.get("H2"),
                values.get("HK2"),
                values.get("RT2"),
                values.get("HK2_T2"),
                values.get("HK2_D2"),
                values.get("HK2_U2"),
                values.get("RT2_T2"),
                values.get("RT2_U2"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", t2, d2, u2);
                if !seen_blkin.contains(&key) && fixtures.blkin.len() < 200 {
                    seen_blkin.insert(key);
                    fixtures.blkin.push(BlkinFixture {
                        input: BlkinInput {
                            t2,
                            d2,
                            u2,
                            hstinv,
                            gm1bl,
                            rstbl,
                            hvrat,
                            reybl,
                        },
                        output: BlkinOutput {
                            m2,
                            h2,
                            hk2,
                            rt2,
                            hk2_t2,
                            hk2_d2,
                            hk2_u2,
                            rt2_t2,
                            rt2_u2,
                        },
                    });
                }
            }
        }
        CurrentSub::Blvar => {
            // ITYP is stored as float (parsed from "ITYP=    1")
            if let (
                Some(&ityp_f),
                Some(&hk2),
                Some(&rt2),
                Some(&m2),
                Some(&t2),
                Some(&d2),
                Some(&s2),
                Some(&hs2),
                Some(&cf2),
                Some(&di2),
                Some(&us2),
                Some(&cq2),
                Some(&de2),
            ) = (
                values.get("ITYP"),
                values.get("HK2"),
                values.get("RT2"),
                values.get("M2"),
                values.get("T2"),
                values.get("D2"),
                values.get("S2"),
                values.get("HS2"),
                values.get("CF2"),
                values.get("DI2"),
                values.get("US2"),
                values.get("CQ2"),
                values.get("DE2"),
            ) {
                let ityp = ityp_f as i32;
                let key = format!("{}_{:.16e}_{:.16e}_{:.16e}", ityp, hk2, rt2, m2);
                if !seen_blvar.contains(&key) && fixtures.blvar.len() < 200 {
                    seen_blvar.insert(key);
                    fixtures.blvar.push(BlvarFixture {
                        input: BlvarInput {
                            ityp,
                            hk2,
                            rt2,
                            m2,
                            t2,
                            d2,
                            s2,
                        },
                        output: BlvarOutput {
                            hs2,
                            cf2,
                            di2,
                            us2,
                            cq2,
                            de2,
                        },
                    });
                }
            }
        }
        CurrentSub::Trchek2 => {
            if let (Some(&x1), Some(&x2), Some(&ampl1), Some(&ampl2), Some(&xt), Some(&tran), Some(&turb)) = (
                values.get("X1"),
                values.get("X2"),
                values.get("AMPL1"),
                values.get("AMPL2"),
                values.get("XT"),
                bools.get("TRAN"),
                bools.get("TURB"),
            ) {
                let key = format!("{:.16e}_{:.16e}_{:.16e}", x1, x2, ampl1);
                if !seen_trchek2.contains(&key) && fixtures.trchek2.len() < 100 {
                    seen_trchek2.insert(key);
                    fixtures.trchek2.push(Trchek2Fixture {
                        input: Trchek2Input { x1, x2, ampl1 },
                        output: Trchek2Output { ampl2, xt, tran, turb },
                    });
                }
            }
        }
        CurrentSub::Blsys => {
            // BLSYS outputs arrays - parse specially
            if let (Some(&vsrez1), Some(&_vsrez2), Some(&_vsrez3), Some(&_vsrez4)) =
                (values.get("VSREZ"), None::<&f64>, None::<&f64>, None::<&f64>)
            // Skip for now
            {
                let _ = (vsrez1, seen_blsys); // Suppress warnings
                                              // BLSYS has complex array output - skip for initial implementation
            }
        }
        _ => {}
    }
}

/// Run XFOIL with instrumentation and generate log file
fn run_xfoil(xfoil_path: &Path, work_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    // Create XFOIL script for NACA 0012 at multiple angles
    let script = r#"PLOP
G F

NACA 0012
OPER
VISC 1e6
ITER 50
ALFA 0
ALFA 2
ALFA 5
ALFA 8
ALFA 10

QUIT
"#;

    let script_path = work_dir.join("xfoil_script.txt");
    let mut file = File::create(&script_path)?;
    file.write_all(script.as_bytes())?;

    // Run XFOIL
    let output = Command::new(xfoil_path)
        .current_dir(work_dir)
        .stdin(File::open(&script_path)?)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()?;

    if !output.success() {
        eprintln!("Warning: XFOIL returned non-zero exit code");
    }

    Ok(())
}

/// Write fixtures to JSON files
fn write_fixtures(fixtures: &ParsedFixtures, output_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    write_fixture_files(&fixtures.hkin, &output_dir.join("hkin"), "hkin")?;
    write_fixture_files(&fixtures.cfl, &output_dir.join("cfl"), "cfl")?;
    write_fixture_files(&fixtures.hsl, &output_dir.join("hsl"), "hsl")?;
    write_fixture_files(&fixtures.dil, &output_dir.join("dil"), "dil")?;
    write_fixture_files(&fixtures.hst, &output_dir.join("hst"), "hst")?;
    write_fixture_files(&fixtures.cft, &output_dir.join("cft"), "cft")?;
    write_fixture_files(&fixtures.dampl, &output_dir.join("dampl"), "dampl")?;
    write_fixture_files(&fixtures.blkin, &output_dir.join("blkin"), "blkin")?;
    write_fixture_files(&fixtures.blvar, &output_dir.join("blvar"), "blvar")?;
    write_fixture_files(&fixtures.trchek2, &output_dir.join("trchek"), "trchek2")?;

    Ok(())
}

fn write_fixture_files<T: Serialize>(fixtures: &[T], dir: &Path, name: &str) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(dir)?;

    // Clear existing files
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.filter_map(Result::ok) {
            if entry.path().extension().is_some_and(|e| e == "json") {
                let _ = fs::remove_file(entry.path());
            }
        }
    }

    for (i, fixture) in fixtures.iter().enumerate() {
        let filename = format!("case_{:03}.json", i + 1);
        let path = dir.join(&filename);
        let file = File::create(&path)?;
        serde_json::to_writer_pretty(file, fixture)?;
    }

    println!("Wrote {} {} fixtures to {}", fixtures.len(), name, dir.display());
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Paths
    let project_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let xfoil_path = project_root.join("xfoil/xfoil6.99/bin/xfoil");
    let work_dir = project_root.join(".tmp/fixture_gen");
    let fixture_dir = project_root.join("tests/fixtures/subroutines");

    // Create work directory
    fs::create_dir_all(&work_dir)?;

    println!("Running instrumented XFOIL to generate log file...");
    run_xfoil(&xfoil_path, &work_dir)?;

    // Check for log file
    let log_path = work_dir.join("xfoil_subroutine_log.dat");
    if !log_path.exists() {
        return Err(format!("Log file not found at {}", log_path.display()).into());
    }

    println!("Parsing log file: {}", log_path.display());
    let fixtures = parse_log(&log_path)?;

    println!("\nParsed fixtures:");
    println!("  HKIN:    {} cases", fixtures.hkin.len());
    println!("  CFL:     {} cases", fixtures.cfl.len());
    println!("  HSL:     {} cases", fixtures.hsl.len());
    println!("  DIL:     {} cases", fixtures.dil.len());
    println!("  HST:     {} cases", fixtures.hst.len());
    println!("  CFT:     {} cases", fixtures.cft.len());
    println!("  DAMPL:   {} cases", fixtures.dampl.len());
    println!("  BLKIN:   {} cases", fixtures.blkin.len());
    println!("  BLVAR:   {} cases", fixtures.blvar.len());
    println!("  TRCHEK2: {} cases", fixtures.trchek2.len());

    println!("\nWriting fixtures to {}", fixture_dir.display());
    write_fixtures(&fixtures, &fixture_dir)?;

    println!("\nDone!");
    Ok(())
}
