//! BLSOLV test fixtures from instrumented XFOIL runs
//!
//! These fixtures capture the exact inputs and outputs of XFOIL's BLSOLV
//! subroutine for validating the yfoil implementation.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Input data for one BLSOLV call
#[derive(Debug, Clone)]
pub struct BlsolvInput {
    /// Call number (1, 2, 3, ...)
    pub call_number: usize,
    /// Total number of BL stations (both sides + wake - overlaps)
    pub nsys: usize,
    /// Index where upper surface ends at TE (1-based)
    pub iblte1: usize,
    /// Index where lower surface ends at TE (1-based)
    pub iblte2: usize,
    /// Number of BL stations on upper surface
    pub nbl1: usize,
    /// Number of BL stations on lower surface
    pub nbl2: usize,
    /// Acceleration parameter for tridiagonal solve
    pub vaccel: f64,
    /// VA matrix: diagonal blocks, [nsys][3][2]
    /// VA[iv][k][l] is the coefficient for equation k, column l at station iv
    pub va: Vec<[[f64; 2]; 3]>,
    /// VB matrix: sub-diagonal blocks, [nsys][3][2]
    /// VB[iv][k][l] is the coupling from station iv-1 to station iv, equation k, column l
    pub vb: Vec<[[f64; 2]; 3]>,
    /// VDEL before solve: residuals and Reynolds sensitivities, [nsys][3][2]
    /// VDEL[iv][k][0] = residual for equation k at station iv
    /// VDEL[iv][k][1] = Reynolds sensitivity for equation k at station iv
    pub vdel_in: Vec<[[f64; 2]; 3]>,
    /// VM matrix: mass defect coupling, [nsys][nsys][3]
    /// VM[iv][j][k] = coupling from station j's mass to station iv's equation k
    /// Note: XFOIL format has VM(3, NSYS, NSYS), so VM(k, j, iv) in Fortran
    pub vm: Vec<Vec<[f64; 3]>>,
    /// VZ matrix: trailing edge coupling block [3][2]
    /// Couples upper surface TE to lower surface/wake
    pub vz: [[f64; 2]; 3],
    /// Total arc length S(N) - S(1) for VACC scaling
    pub arc_length: Option<f64>,
}

/// Output data for one BLSOLV call
#[derive(Debug, Clone)]
pub struct BlsolvOutput {
    /// Call number (1, 2, 3, ...)
    pub call_number: usize,
    /// Total number of BL stations
    pub nsys: usize,
    /// VDEL after solve: solution deltas, [nsys][3][2]
    /// VDEL[iv][k][0] = solution delta for equation k at station iv
    /// VDEL[iv][k][1] = Reynolds sensitivity solution at station iv
    pub vdel_out: Vec<[[f64; 2]; 3]>,
}

impl BlsolvInput {
    /// Get a small subset for testing (first n_small stations)
    pub fn small_subset(&self, n_small: usize) -> BlsolvInput {
        let n = n_small.min(self.nsys);
        BlsolvInput {
            call_number: self.call_number,
            nsys: n,
            iblte1: self.iblte1.min(n),
            iblte2: self.iblte2.min(n),
            nbl1: self.nbl1.min(n),
            nbl2: self.nbl2.min(n),
            vaccel: self.vaccel,
            va: self.va[..n].to_vec(),
            vb: self.vb[..n].to_vec(),
            vdel_in: self.vdel_in[..n].to_vec(),
            vm: self.vm[..n]
                .iter()
                .map(|row| row[..n].to_vec())
                .collect(),
            vz: self.vz,
            arc_length: self.arc_length,
        }
    }

    /// Compute the system index for upper surface trailing edge (0-based)
    /// In XFOIL: IVTE1 = ISYS(IBLTE(1),1)
    ///
    /// XFOIL ISYS mapping is built by:
    ///   DO IBL=2, NBL(IS)
    ///     IV = IV+1
    ///     ISYS(IBL,IS) = IV
    ///
    /// For surface 1, this gives: ISYS(IBL,1) = IBL - 1 (1-based)
    /// So IVTE1 = ISYS(IBLTE1, 1) = IBLTE1 - 1 (1-based)
    /// For 0-based: ivte1 = IVTE1 - 1 = IBLTE1 - 2
    ///
    /// Example: IBLTE1=88 → IVTE1=87 (1-based) → ivte1=86 (0-based)
    pub fn ivte1_0based(&self) -> usize {
        self.iblte1 - 2
    }

    /// Compute the system index for first wake station (0-based)
    /// In XFOIL: IVZ = ISYS(IBLTE(2)+1,2)
    ///
    /// XFOIL ISYS mapping for surface 2 continues from surface 1:
    ///   For IBL=2 to NBL(2): ISYS(IBL,2) = NBL1 + IBL - 2 (1-based)
    ///
    /// First wake station is BL index IBLTE2+1 on surface 2
    /// IVZ = ISYS(IBLTE2+1, 2) = NBL1 + (IBLTE2+1) - 2 = NBL1 + IBLTE2 - 1 (1-based)
    /// For 0-based: ivz = IVZ - 1 = NBL1 + IBLTE2 - 2
    ///
    /// Example: NBL1=88, IBLTE2=74 → IVZ=161 (1-based) → ivz=160 (0-based)
    pub fn ivz_0based(&self) -> usize {
        self.nbl1 + self.iblte2 - 2
    }
}

/// Parse BLSOLV input fixture file
pub fn parse_blsolv_input(path: &Path, call_number: usize) -> Option<BlsolvInput> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    let mut idx = 0;

    // Find the requested call - format is "=== BLSOLV_CALL     N"
    // Use regex-like matching that handles variable whitespace
    while idx < lines.len() {
        let line = &lines[idx];
        if line.contains("BLSOLV_CALL") {
            // Extract the call number from the line
            if let Some(num_str) = line.split("BLSOLV_CALL").nth(1) {
                if let Ok(found_num) = num_str.trim().parse::<usize>() {
                    if found_num == call_number {
                        break;
                    }
                }
            }
        }
        idx += 1;
    }
    if idx >= lines.len() {
        return None;
    }
    idx += 1; // Move past header

    // Parse NSYS
    let nsys = parse_int_value(&lines[idx])?;
    idx += 1;

    // Parse IBLTE1
    let iblte1 = parse_int_value(&lines[idx])?;
    idx += 1;

    // Parse IBLTE2
    let iblte2 = parse_int_value(&lines[idx])?;
    idx += 1;

    // Parse NBL1
    let nbl1 = parse_int_value(&lines[idx])?;
    idx += 1;

    // Parse NBL2
    let nbl2 = parse_int_value(&lines[idx])?;
    idx += 1;

    // Parse VACCEL
    let vaccel = parse_float_value(&lines[idx])?;
    idx += 1;

    // Parse ARC_LENGTH (optional - may not be present in older fixtures)
    let arc_length = if lines[idx].contains("ARC_LENGTH") {
        let val = parse_float_value(&lines[idx]);
        idx += 1;
        val
    } else {
        None
    };

    let mut va = vec![[[0.0; 2]; 3]; nsys];
    let mut vb = vec![[[0.0; 2]; 3]; nsys];
    let mut vdel_in = vec![[[0.0; 2]; 3]; nsys];
    let mut vm = vec![vec![[0.0; 3]; nsys]; nsys];

    // Parse each station
    for iv in 0..nsys {
        // Skip "--- IV= N" line
        if !lines[idx].contains("IV=") {
            eprintln!("Expected IV= line at index {}, got: {}", idx, &lines[idx]);
            return None;
        }
        idx += 1;

        // Parse VA (6 values as 3 rows × 2 cols)
        let va_vals = parse_labeled_floats(&lines[idx], "VA=")?;
        if va_vals.len() < 6 {
            eprintln!("VA has {} values, expected 6", va_vals.len());
            return None;
        }
        va[iv] = [
            [va_vals[0], va_vals[1]],
            [va_vals[2], va_vals[3]],
            [va_vals[4], va_vals[5]],
        ];
        idx += 1;

        // Parse VB (6 values)
        let vb_vals = parse_labeled_floats(&lines[idx], "VB=")?;
        if vb_vals.len() < 6 {
            eprintln!("VB has {} values, expected 6", vb_vals.len());
            return None;
        }
        vb[iv] = [
            [vb_vals[0], vb_vals[1]],
            [vb_vals[2], vb_vals[3]],
            [vb_vals[4], vb_vals[5]],
        ];
        idx += 1;

        // Parse VDEL (6 values)
        let vdel_vals = parse_labeled_floats(&lines[idx], "VDEL=")?;
        if vdel_vals.len() < 6 {
            eprintln!("VDEL has {} values, expected 6", vdel_vals.len());
            return None;
        }
        vdel_in[iv] = [
            [vdel_vals[0], vdel_vals[1]],
            [vdel_vals[2], vdel_vals[3]],
            [vdel_vals[4], vdel_vals[5]],
        ];
        idx += 1;

        // Parse VM entries for each station j (from station j to equations at station iv)
        for j in 0..nsys {
            // VM( j+1) = val1 val2 val3
            let vm_vals = parse_vm_line(&lines[idx])?;
            if vm_vals.len() < 3 {
                eprintln!("VM has {} values, expected 3 at iv={}, j={}", vm_vals.len(), iv, j);
                return None;
            }
            vm[iv][j] = [vm_vals[0], vm_vals[1], vm_vals[2]];
            idx += 1;
        }
    }

    // Parse VZ block (trailing edge coupling)
    // Format: VZ= val1 val2 val3 val4 val5 val6 (3 rows × 2 cols)
    let vz_vals = parse_labeled_floats(&lines[idx], "VZ=")?;
    if vz_vals.len() < 6 {
        eprintln!("VZ has {} values, expected 6", vz_vals.len());
        return None;
    }
    let vz = [
        [vz_vals[0], vz_vals[1]],
        [vz_vals[2], vz_vals[3]],
        [vz_vals[4], vz_vals[5]],
    ];

    Some(BlsolvInput {
        call_number,
        nsys,
        iblte1,
        iblte2,
        nbl1,
        nbl2,
        vaccel,
        va,
        vb,
        vdel_in,
        vm,
        vz,
        arc_length,
    })
}

/// Parse BLSOLV output fixture file
pub fn parse_blsolv_output(path: &Path, call_number: usize) -> Option<BlsolvOutput> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    let mut idx = 0;

    // Find the requested call - format is "=== BLSOLV_CALL     N"
    while idx < lines.len() {
        let line = &lines[idx];
        if line.contains("BLSOLV_CALL") {
            if let Some(num_str) = line.split("BLSOLV_CALL").nth(1) {
                if let Ok(found_num) = num_str.trim().parse::<usize>() {
                    if found_num == call_number {
                        break;
                    }
                }
            }
        }
        idx += 1;
    }
    if idx >= lines.len() {
        return None;
    }
    idx += 1;

    // Parse NSYS
    let nsys = parse_int_value(&lines[idx])?;
    idx += 1;

    let mut vdel_out = vec![[[0.0; 2]; 3]; nsys];

    // Parse VDEL for each station
    for iv in 0..nsys {
        // VDEL( iv+1)= val1 val2 val3 val4 val5 val6
        let vdel_vals = parse_vdel_out_line(&lines[idx])?;
        if vdel_vals.len() < 6 {
            eprintln!("VDEL output has {} values, expected 6 at iv={}", vdel_vals.len(), iv);
            return None;
        }
        vdel_out[iv] = [
            [vdel_vals[0], vdel_vals[1]],
            [vdel_vals[2], vdel_vals[3]],
            [vdel_vals[4], vdel_vals[5]],
        ];
        idx += 1;
    }

    Some(BlsolvOutput {
        call_number,
        nsys,
        vdel_out,
    })
}

fn parse_int_value(line: &str) -> Option<usize> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() >= 2 {
        parts[1].trim().parse().ok()
    } else {
        None
    }
}

fn parse_float_value(line: &str) -> Option<f64> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() >= 2 {
        parts[1].trim().parse().ok()
    } else {
        None
    }
}

fn parse_labeled_floats(line: &str, label: &str) -> Option<Vec<f64>> {
    let start = line.find(label)? + label.len();
    let remainder = &line[start..];
    let values: Vec<f64> = remainder
        .split_whitespace()
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    if values.is_empty() {
        None
    } else {
        Some(values)
    }
}

fn parse_vm_line(line: &str) -> Option<Vec<f64>> {
    // Format: VM(   1)=  0.0000000000E+00  0.1188548718E+08 -0.8789180762E+07
    let eq_pos = line.find('=')?;
    let remainder = &line[eq_pos + 1..];
    let values: Vec<f64> = remainder
        .split_whitespace()
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    Some(values)
}

fn parse_vdel_out_line(line: &str) -> Option<Vec<f64>> {
    // Format: VDEL(   1)=  0.0000000000E+00  0.0000000000E+00  ...
    let eq_pos = line.find('=')?;
    let remainder = &line[eq_pos + 1..];
    let values: Vec<f64> = remainder
        .split_whitespace()
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    Some(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(".tmp")
            .join(name)
    }

    #[test]
    fn test_parse_blsolv_input_call_1() {
        let path = fixture_path("blsolv_input.dat");
        if !path.exists() {
            eprintln!("Skipping test: fixture file not found at {:?}", path);
            return;
        }

        let input = parse_blsolv_input(&path, 1);
        assert!(input.is_some(), "Failed to parse BLSOLV input");
        let input = input.unwrap();

        assert_eq!(input.call_number, 1);
        assert_eq!(input.nsys, 183);
        // Fixture values from NACA 0012 at Re=1e6, alpha=5 deg
        assert_eq!(input.iblte1, 88);
        assert_eq!(input.iblte2, 74);
        assert_eq!(input.nbl1, 88);
        assert_eq!(input.nbl2, 97);
        assert!((input.vaccel - 0.01).abs() < 1e-10);

        // Check first station VA - value should be ~1.0
        assert!((input.va[0][0][0] - 1.0).abs() < 1e-10, "VA[0][0][0] should be 1.0");

        // Check that we have valid data (not all zeros)
        assert!(input.vdel_in[0][1][0] != 0.0 || input.vdel_in[0][2][0] != 0.0,
            "VDEL should have non-zero values");

        // Check VM structure has non-zero entries
        assert!(input.vm[0][0][2].abs() > 0.0,
            "VM[0][0][2] should be non-zero for mass defect coupling");
    }

    #[test]
    fn test_parse_blsolv_output_call_1() {
        let path = fixture_path("blsolv_output.dat");
        if !path.exists() {
            eprintln!("Skipping test: fixture file not found at {:?}", path);
            return;
        }

        let output = parse_blsolv_output(&path, 1);
        assert!(output.is_some(), "Failed to parse BLSOLV output");
        let output = output.unwrap();

        assert_eq!(output.call_number, 1);
        assert_eq!(output.nsys, 183);
        assert_eq!(output.vdel_out.len(), 183);

        // Check first station output has valid data (not all zeros)
        // The exact values depend on the fixture which may change
        assert!(output.vdel_out[0][1][0] != 0.0 || output.vdel_out[0][2][0] != 0.0,
            "VDEL_out should have non-zero values at first station");
    }

    #[test]
    fn test_small_subset() {
        let path = fixture_path("blsolv_input.dat");
        if !path.exists() {
            return;
        }

        let input = parse_blsolv_input(&path, 1).unwrap();
        let small = input.small_subset(10);

        assert_eq!(small.nsys, 10);
        assert_eq!(small.va.len(), 10);
        assert_eq!(small.vm.len(), 10);
        assert_eq!(small.vm[0].len(), 10);
    }
}
