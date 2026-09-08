//! Parsers for the S2 instrumentation dumps: xfoil_pointers.dat (IBLSYS), xfoil_uinv.dat
//! (UICALC) and the per-node GAM block of xfoil_inviscid.dat.

use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct PointersFixture {
    pub call: usize,
    pub ist: usize,
    pub n: usize,
    pub nw: usize,
    pub nsys: usize,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,
    pub ante: f64,
    pub aste: f64,
    pub dste: f64,
    pub chord: f64,
    pub sle: f64,
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    pub sharp: bool,
    pub nbl: [usize; 3],
    pub iblte: [usize; 3],
    pub itran: [usize; 3],
    /// 1-based: x[i], y[i], s[i] for i in 1..=n+nw
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub s: Vec<f64>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub apanel: Vec<f64>,
    /// [is][ibl], 1-based, ibl in 1..=iblte[is]+nw
    pub ipan: [Vec<usize>; 3],
    pub vti: [Vec<f64>; 3],
    pub isys: [Vec<usize>; 3],
    pub xssi: [Vec<f64>; 3],
    /// 1-based, 1..=nw
    pub wgap: Vec<f64>,
}

fn val<'a>(line: &'a str, key: &str) -> &'a str {
    line.strip_prefix(key)
        .unwrap_or_else(|| panic!("expected `{key}` in `{line}`"))
        .trim()
}
fn ints(s: &str) -> Vec<usize> {
    s.split_whitespace().map(|t| t.parse().unwrap()).collect()
}

pub fn parse_pointers(path: &Path, call: usize) -> PointersFixture {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines().peekable();
    // locate call
    loop {
        let l = lines.next().unwrap_or_else(|| panic!("IBLSYS_CALL {call} not found"));
        if let Some(rest) = l.strip_prefix("=== IBLSYS_CALL") {
            if rest.trim().parse::<usize>().unwrap() == call {
                break;
            }
        }
    }
    let mut f = PointersFixture {
        call,
        ..Default::default()
    };
    f.ist = val(lines.next().unwrap(), "IST=").parse().unwrap();
    f.n = val(lines.next().unwrap(), "N=").parse().unwrap();
    f.nw = val(lines.next().unwrap(), "NW=").parse().unwrap();
    f.nsys = val(lines.next().unwrap(), "NSYS=").parse().unwrap();
    f.sst = val(lines.next().unwrap(), "SST=").parse().unwrap();
    f.sst_go = val(lines.next().unwrap(), "SST_GO=").parse().unwrap();
    f.sst_gp = val(lines.next().unwrap(), "SST_GP=").parse().unwrap();
    f.ante = val(lines.next().unwrap(), "ANTE=").parse().unwrap();
    f.aste = val(lines.next().unwrap(), "ASTE=").parse().unwrap();
    f.dste = val(lines.next().unwrap(), "DSTE=").parse().unwrap();
    f.chord = val(lines.next().unwrap(), "CHORD=").parse().unwrap();
    f.sle = val(lines.next().unwrap(), "SLE=").parse().unwrap();
    f.xle = val(lines.next().unwrap(), "XLE=").parse().unwrap();
    f.yle = val(lines.next().unwrap(), "YLE=").parse().unwrap();
    f.xte = val(lines.next().unwrap(), "XTE=").parse().unwrap();
    f.yte = val(lines.next().unwrap(), "YTE=").parse().unwrap();
    f.sharp = val(lines.next().unwrap(), "SHARP=") == "T";
    let v = ints(val(lines.next().unwrap(), "NBL="));
    f.nbl = [0, v[0], v[1]];
    let v = ints(val(lines.next().unwrap(), "IBLTE="));
    f.iblte = [0, v[0], v[1]];
    let v = ints(val(lines.next().unwrap(), "ITRAN="));
    f.itran = [0, v[0], v[1]];
    assert!(lines.next().unwrap().starts_with("--- NODES"));
    let np = f.n + f.nw + 1;
    f.x = vec![0.0; np];
    f.y = vec![0.0; np];
    f.s = vec![0.0; np];
    f.nx = vec![0.0; np];
    f.ny = vec![0.0; np];
    f.apanel = vec![0.0; np];
    for _ in 1..np {
        let p: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let i: usize = p[0].parse().unwrap();
        f.x[i] = p[1].parse().unwrap();
        f.y[i] = p[2].parse().unwrap();
        f.s[i] = p[3].parse().unwrap();
        f.nx[i] = p[4].parse().unwrap();
        f.ny[i] = p[5].parse().unwrap();
        f.apanel[i] = p[6].parse().unwrap();
    }
    let ivx = f.n + f.nw + 2;
    f.ipan = [Vec::new(), vec![0; ivx], vec![0; ivx]];
    f.vti = [Vec::new(), vec![0.0; ivx], vec![0.0; ivx]];
    f.isys = [Vec::new(), vec![0; ivx], vec![0; ivx]];
    f.xssi = [Vec::new(), vec![0.0; ivx], vec![0.0; ivx]];
    for is in 1..=2 {
        let hdr = lines.next().unwrap();
        assert!(hdr.starts_with("--- SIDE"), "got `{hdr}`");
        for _ in 1..=(f.iblte[is] + f.nw) {
            let p: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
            let ibl: usize = p[0].parse().unwrap();
            f.ipan[is][ibl] = p[1].parse().unwrap();
            f.vti[is][ibl] = p[2].parse().unwrap();
            f.isys[is][ibl] = p[3].parse().unwrap();
            f.xssi[is][ibl] = p[4].parse().unwrap();
        }
    }
    assert!(lines.next().unwrap().starts_with("--- WGAP"));
    f.wgap = vec![0.0; f.nw + 1];
    for _ in 1..=f.nw {
        let p: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let iw: usize = p[0].parse().unwrap();
        f.wgap[iw] = p[1].parse().unwrap();
    }
    f
}

#[derive(Debug, Clone, Default)]
pub struct UinvFixture {
    pub call: usize,
    pub alfa: f64,
    pub n: usize,
    pub nw: usize,
    /// 1-based, 1..=n+nw
    pub qinv: Vec<f64>,
    pub qinv_a: Vec<f64>,
    pub qinvu1: Vec<f64>,
    pub qinvu2: Vec<f64>,
    /// [is][ibl], 1..=nbl[is]
    pub uinv: [Vec<f64>; 3],
    pub uinv_a: [Vec<f64>; 3],
}

pub fn parse_uinv(path: &Path, call: usize) -> UinvFixture {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    loop {
        let l = lines.next().unwrap_or_else(|| panic!("UICALC_CALL {call} not found"));
        if let Some(rest) = l.strip_prefix("=== UICALC_CALL") {
            if rest.trim().parse::<usize>().unwrap() == call {
                break;
            }
        }
    }
    let mut f = UinvFixture {
        call,
        ..Default::default()
    };
    f.alfa = val(lines.next().unwrap(), "ALFA=").parse().unwrap();
    let v = ints(val(lines.next().unwrap(), "N NW="));
    f.n = v[0];
    f.nw = v[1];
    assert!(lines.next().unwrap().starts_with("--- NODES"));
    let np = f.n + f.nw + 1;
    f.qinv = vec![0.0; np];
    f.qinv_a = vec![0.0; np];
    f.qinvu1 = vec![0.0; np];
    f.qinvu2 = vec![0.0; np];
    for _ in 1..np {
        let p: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let i: usize = p[0].parse().unwrap();
        f.qinv[i] = p[1].parse().unwrap();
        f.qinv_a[i] = p[2].parse().unwrap();
        f.qinvu1[i] = p[3].parse().unwrap();
        f.qinvu2[i] = p[4].parse().unwrap();
    }
    f.uinv = [Vec::new(), Vec::new(), Vec::new()];
    f.uinv_a = [Vec::new(), Vec::new(), Vec::new()];
    for is in 1..=2 {
        assert!(lines.next().unwrap().starts_with("--- SIDE"));
        let mut u = vec![0.0];
        let mut ua = vec![0.0];
        // read until the next header or EOF; rows are `IBL UINV UINV_A`
        while let Some(l) = lines.clone().next() {
            if l.starts_with("---") || l.starts_with("===") {
                break;
            }
            let l = lines.next().unwrap();
            let p: Vec<&str> = l.split_whitespace().collect();
            u.push(p[1].parse().unwrap());
            ua.push(p[2].parse().unwrap());
        }
        f.uinv[is] = u;
        f.uinv_a[is] = ua;
    }
    f
}

/// GAM(I), 1-based, from the `I, GAM, QINV, CPI` block of xfoil_inviscid.dat.
pub fn parse_inviscid_gam(path: &Path) -> Vec<f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut gam = vec![0.0];
    let mut on = false;
    for l in text.lines() {
        if l.contains("I, GAM, QINV") {
            on = true;
            continue;
        }
        if !on {
            continue;
        }
        let p: Vec<&str> = l.split_whitespace().collect();
        if p.len() < 4 || p[0].parse::<usize>().is_err() {
            break;
        }
        gam.push(p[1].parse().unwrap());
    }
    gam
}

/// DIJ(I,J) from xfoil_dij.dat, 1-based (N+NW)×(N+NW); returns (n, nw, dij).
pub fn parse_dij(path: &Path) -> (usize, usize, Vec<Vec<f64>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut n = 0;
    let mut nw = 0;
    let mut dij: Vec<Vec<f64>> = Vec::new();
    for l in text.lines() {
        let t = l.trim();
        if let Some(v) = t.strip_prefix("N =") {
            n = v.trim().parse().unwrap();
        } else if let Some(v) = t.strip_prefix("NW =") {
            nw = v.trim().parse().unwrap();
            dij = vec![vec![0.0; n + nw + 1]; n + nw + 1];
        } else {
            let p: Vec<&str> = t.split_whitespace().collect();
            if p.len() == 3 {
                if let (Ok(i), Ok(j)) = (p[0].parse::<usize>(), p[1].parse::<usize>()) {
                    dij[i][j] = p[2].parse().unwrap();
                }
            }
        }
    }
    (n, nw, dij)
}

/// The inviscid arrays entering VISCAL's Newton loop (`viscal_inviscid.dat`): per node
/// QINVU(I,1), QINVU(I,2), QINV(I), QINV_A(I), GAM(I), GAM_A(I), 1-based; plus the header.
#[derive(Debug, Clone, Default)]
pub struct ViscalInviscid {
    pub header: std::collections::HashMap<String, String>,
    pub qinvu1: Vec<f64>,
    pub qinvu2: Vec<f64>,
    pub qinv: Vec<f64>,
    pub qinv_a: Vec<f64>,
    pub gam: Vec<f64>,
    pub gam_a: Vec<f64>,
}

pub fn parse_viscal_inviscid(path: &Path) -> ViscalInviscid {
    let text = std::fs::read_to_string(path).unwrap();
    let mut d = ViscalInviscid::default();
    for v in [
        &mut d.qinvu1,
        &mut d.qinvu2,
        &mut d.qinv,
        &mut d.qinv_a,
        &mut d.gam,
        &mut d.gam_a,
    ] {
        v.push(0.0);
    }
    for l in text.lines() {
        if let Some(r) = l.strip_prefix("NODE(") {
            let (_, vals) = r.split_once(")=").unwrap();
            let v: Vec<f64> = vals.split_whitespace().map(|t| t.parse().unwrap()).collect();
            d.qinvu1.push(v[0]);
            d.qinvu2.push(v[1]);
            d.qinv.push(v[2]);
            d.qinv_a.push(v[3]);
            d.gam.push(v[4]);
            d.gam_a.push(v[5]);
        } else if let Some((k, v)) = l.split_once('=') {
            d.header.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    d
}
