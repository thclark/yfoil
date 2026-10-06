//! Parsers for the MRCHDU instrumentation: `mrchdu_input_<k>.dat` / `mrchdu_output_<k>.dat`
//! (BL state before/after MRCHDU on SETBL call k) and `xfoil_mrchdu_trace.dat` (per-iteration
//! Newton trace on call 1).

use std::collections::HashMap;
use std::path::Path;

/// BL state dump: a `KEY=value` header followed by `BL(is,ibl)=` (XSSI UEDG THET DSTR CTAU MASS)
/// and, on output dumps, `BLX(is,ibl)=` (TAU DIS CTQ DELT TSTR) lines. 1-based, index 0 unused.
#[derive(Debug, Clone, Default)]
pub struct BlDump {
    pub header: HashMap<String, String>,
    pub bl: [Vec<[f64; 6]>; 3],
    pub blx: [Vec<[f64; 5]>; 3],
}

impl BlDump {
    pub fn int(&self, key: &str) -> usize {
        self.header
            .get(key)
            .unwrap_or_else(|| panic!("missing header `{key}`"))
            .parse()
            .unwrap()
    }
    pub fn real(&self, key: &str) -> f64 {
        self.header
            .get(key)
            .unwrap_or_else(|| panic!("missing header `{key}`"))
            .parse()
            .unwrap()
    }
    pub fn logical(&self, key: &str) -> bool {
        self.header.get(key).unwrap_or_else(|| panic!("missing header `{key}`")) == "T"
    }
    /// NBL(IS) as the number of `BL(` rows seen for that side
    pub fn nbl(&self, is: usize) -> usize {
        self.bl[is].len() - 1
    }
}

pub fn parse_bl_dump(path: &Path) -> BlDump {
    let text = std::fs::read_to_string(path).unwrap();
    let mut d = BlDump::default();
    let mut rows: [Vec<(usize, [f64; 6])>; 3] = Default::default();
    let mut rowsx: [Vec<(usize, [f64; 5])>; 3] = Default::default();
    for l in text.lines() {
        if let Some(r) = l.strip_prefix("BLX(") {
            let (idx, vals) = r.split_once(")=").unwrap();
            let (is, ibl) = idx.split_once(',').unwrap();
            let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
            rowsx[is.trim().parse::<usize>().unwrap()]
                .push((ibl.trim().parse().unwrap(), [v[0], v[1], v[2], v[3], v[4]]));
        } else if let Some(r) = l.strip_prefix("BL(") {
            let (idx, vals) = r.split_once(")=").unwrap();
            let (is, ibl) = idx.split_once(',').unwrap();
            let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
            rows[is.trim().parse::<usize>().unwrap()]
                .push((ibl.trim().parse().unwrap(), [v[0], v[1], v[2], v[3], v[4], v[5]]));
        } else if let Some((k, v)) = l.split_once('=') {
            d.header.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    for is in 1..=2 {
        let n = rows[is].iter().map(|r| r.0).max().unwrap_or(0);
        d.bl[is] = vec![[0.0; 6]; n + 1];
        for (ibl, v) in &rows[is] {
            d.bl[is][*ibl] = *v;
        }
        let nx = rowsx[is].iter().map(|r| r.0).max().unwrap_or(0);
        d.blx[is] = vec![[0.0; 5]; nx + 1];
        for (ibl, v) in &rowsx[is] {
            d.blx[is][*ibl] = *v;
        }
    }
    d
}

/// One MRCHDU Newton iteration from the reference trace.
#[derive(Debug, Clone, Default)]
pub struct MrchduTraceIter {
    pub ibl: usize,
    pub is: usize,
    pub itbl: usize,
    pub ampl: [f64; 4],
    pub tran: bool,
    pub itran: usize,
    pub primary: [f64; 5],
    pub kinematic: [f64; 5],
    pub closure: [f64; 5],
    pub ueref: f64,
    pub hkref: f64,
    pub sens: Option<[f64; 2]>,
    pub residual: [f64; 4],
    pub vs2: [[f64; 5]; 4],
    pub solution: [f64; 4],
    pub dmax: f64,
    pub rlx: f64,
    pub updated: [f64; 5],
    pub converged: bool,
}

#[derive(Debug, Clone, Default)]
pub struct MrchduTraceDump {
    pub iters: Vec<MrchduTraceIter>,
    /// (IBL, IS, DMAX) FAILED records
    pub failed: Vec<(usize, usize, f64)>,
}

fn nums(l: &str) -> Vec<f64> {
    l.split_once(':')
        .unwrap()
        .1
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

pub fn parse_mrchdu_trace(path: &Path) -> MrchduTraceDump {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = MrchduTraceDump::default();
    for l in text.lines() {
        if l.starts_with('#') {
            continue;
        }
        if let Some(r) = l.strip_prefix("STATION") {
            let v: Vec<usize> = r.split_whitespace().map(|t| t.parse().unwrap()).collect();
            out.iters.push(MrchduTraceIter {
                ibl: v[0],
                is: v[1],
                itbl: v[2],
                ..Default::default()
            });
            continue;
        }
        if let Some(r) = l.strip_prefix("FAILED:") {
            let t: Vec<&str> = r.split_whitespace().collect();
            out.failed
                .push((t[0].parse().unwrap(), t[1].parse().unwrap(), t[2].parse().unwrap()));
            continue;
        }
        let Some(cur) = out.iters.last_mut() else { continue };
        if let Some(r) = l.strip_prefix("AMPL:") {
            let t: Vec<&str> = r.split_whitespace().collect();
            for i in 0..4 {
                cur.ampl[i] = t[i].parse().unwrap();
            }
            cur.tran = t[4] == "T";
            cur.itran = t[5].parse().unwrap();
        } else if l.starts_with("PRIMARY:") {
            cur.primary.copy_from_slice(&nums(l)[..5]);
        } else if l.starts_with("KINEMATIC:") {
            cur.kinematic.copy_from_slice(&nums(l)[..5]);
        } else if l.starts_with("CLOSURE:") {
            cur.closure.copy_from_slice(&nums(l)[..5]);
        } else if l.starts_with("REF:") {
            let v = nums(l);
            cur.ueref = v[0];
            cur.hkref = v[1];
        } else if l.starts_with("SENS:") {
            let v = nums(l);
            cur.sens = Some([v[0], v[1]]);
        } else if l.starts_with("RESIDUAL:") {
            cur.residual.copy_from_slice(&nums(l)[..4]);
        } else if let Some(r) = l.strip_prefix("VS2_") {
            let (k, rest) = r.split_once(':').unwrap();
            let k: usize = k.parse().unwrap();
            let v: Vec<f64> = rest.split_whitespace().map(|t| t.parse().unwrap()).collect();
            cur.vs2[k - 1].copy_from_slice(&v[..5]);
        } else if l.starts_with("SOLUTION:") {
            cur.solution.copy_from_slice(&nums(l)[..4]);
        } else if l.starts_with("RELAX:") {
            let v = nums(l);
            cur.dmax = v[0];
            cur.rlx = v[1];
        } else if l.starts_with("UPDATED:") {
            cur.updated.copy_from_slice(&nums(l)[..5]);
        } else if l.starts_with("CONVERGED:") {
            cur.converged = true;
        }
    }
    out
}
