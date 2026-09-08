//! Parsers for the MRCHUE/MRCHDU instrumentation: `mrchdu_input.dat` (state after MRCHUE,
//! before the first MRCHDU) and `xfoil_newton_trace.dat` (MRCHUE per-iteration trace).

use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct BlStateDump {
    pub nbl: [usize; 3],
    pub iblte: [usize; 3],
    pub minf: f64,
    pub reinf: f64,
    pub acrit: [f64; 3],
    pub itran: [usize; 3],
    pub reybl: f64,
    pub hstinv: f64,
    pub gm1bl: f64,
    pub ante: f64,
    /// [is][ibl] -> (XSSI, UEDG, THET, DSTR, CTAU, MASS), 1-based
    pub bl: [Vec<[f64; 6]>; 3],
}

fn val<'a>(l: &'a str, key: &str) -> &'a str {
    l.strip_prefix(key)
        .unwrap_or_else(|| panic!("expected `{key}` in `{l}`"))
        .trim()
}

pub fn parse_bl_state(path: &Path) -> BlStateDump {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let mut d = BlStateDump::default();
    d.nbl[1] = val(lines.next().unwrap(), "NBL1=").parse().unwrap();
    d.nbl[2] = val(lines.next().unwrap(), "NBL2=").parse().unwrap();
    d.iblte[1] = val(lines.next().unwrap(), "IBLTE1=").parse().unwrap();
    d.iblte[2] = val(lines.next().unwrap(), "IBLTE2=").parse().unwrap();
    d.minf = val(lines.next().unwrap(), "MINF=").parse().unwrap();
    d.reinf = val(lines.next().unwrap(), "REINF=").parse().unwrap();
    d.acrit[1] = val(lines.next().unwrap(), "ACRIT1=").parse().unwrap();
    d.acrit[2] = val(lines.next().unwrap(), "ACRIT2=").parse().unwrap();
    // optional extended header (added in S5)
    let mut rest: Vec<&str> = lines.collect();
    while let Some(l) = rest.first() {
        if let Some(v) = l.strip_prefix("ITRAN1=") {
            d.itran[1] = v.trim().parse().unwrap();
        } else if let Some(v) = l.strip_prefix("ITRAN2=") {
            d.itran[2] = v.trim().parse().unwrap();
        } else if let Some(v) = l.strip_prefix("REYBL=") {
            d.reybl = v.trim().parse().unwrap();
        } else if let Some(v) = l.strip_prefix("HSTINV=") {
            d.hstinv = v.trim().parse().unwrap();
        } else if let Some(v) = l.strip_prefix("GM1BL=") {
            d.gm1bl = v.trim().parse().unwrap();
        } else if let Some(v) = l.strip_prefix("ANTE=") {
            d.ante = v.trim().parse().unwrap();
        } else {
            break;
        }
        rest.remove(0);
    }
    let n = d.nbl[1].max(d.nbl[2]) + 1;
    d.bl = [Vec::new(), vec![[0.0; 6]; n], vec![[0.0; 6]; n]];
    for l in rest {
        // BL( 1,   2)=  xssi uedg thet dstr ctau mass
        let Some(r) = l.strip_prefix("BL(") else { continue };
        let (idx, vals) = r.split_once(")=").unwrap();
        let (is, ibl) = idx.split_once(',').unwrap();
        let is: usize = is.trim().parse().unwrap();
        let ibl: usize = ibl.trim().parse().unwrap();
        let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
        d.bl[is][ibl] = [v[0], v[1], v[2], v[3], v[4], v[5]];
    }
    d
}

/// One MRCHUE Newton iteration from the reference trace.
#[derive(Debug, Clone, Default)]
pub struct TraceIter {
    pub ibl: usize,
    pub is: usize,
    pub itbl: usize,
    pub ampl: [f64; 4],
    pub tran: bool,
    pub itran: usize,
    pub primary: [f64; 5],
    pub kinematic: [f64; 5],
    pub closure: [f64; 5],
    pub residual: [f64; 4],
    pub vs2: [[f64; 5]; 3],
    pub solution: [f64; 4],
    pub dmax: f64,
    pub rlx: f64,
    pub updated: [f64; 4],
    pub has_update: bool,
    pub converged: bool,
}

fn nums(l: &str) -> Vec<f64> {
    l.split_once(':')
        .unwrap()
        .1
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect()
}

pub fn parse_newton_trace(path: &Path) -> Vec<TraceIter> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out: Vec<TraceIter> = Vec::new();
    for l in text.lines() {
        if let Some(r) = l.strip_prefix("STATION") {
            let v: Vec<usize> = r.split_whitespace().map(|t| t.parse().unwrap()).collect();
            out.push(TraceIter {
                ibl: v[0],
                is: v[1],
                itbl: v[2],
                ..Default::default()
            });
            continue;
        }
        let Some(cur) = out.last_mut() else { continue };
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
            cur.updated.copy_from_slice(&nums(l)[..4]);
            cur.has_update = true;
        } else if l.starts_with("CONVERGED:") {
            cur.converged = true;
        }
    }
    out
}
