//! Bit-for-bit comparison of pure-Rust maths crates against cr-libm (CORE-MATH, the C that yFoil
//! and the XFOIL reference are linked against), on the twelve functions yFoil uses.
//!
//! A correctly rounded function has one right answer per argument, so a candidate that is truly
//! correctly rounded agrees with CORE-MATH on every argument; NaN payloads aside, which IEEE 754
//! leaves open and which are counted apart. This is the gate for replacing the C.
//!
//! ```text
//! libm-compare compare [--millions N] [--worst-cases DIR] [--only exp,pow] [--candidates pxfm,metallic]
//!                      [--digests-out FILE] [--digests-in FILE] [--threads N]
//! libm-compare time [--only ...] [--candidates ...]
//! libm-compare dump FUNCTION SET BLOCK CANDIDATE [--millions N] [--worst-cases DIR]
//! ```
//!
//! `compare` against the C needs the `c-reference` feature (the default). `--digests-out` records
//! every implementation's per-block digests; built without the C (for a target with no C
//! compiler: `--no-default-features`), `--digests-in` checks each candidate against the C's
//! recorded digests and against its own (so a target-dependent result shows even where the
//! candidate is wrong everywhere); a differing block is located with `dump` on both sides. `run.sh` drives the native, wasm and x86-64 runs, and
//! `fetch-worst-cases.sh` fetches CORE-MATH's worst-case lists at the vendored commit.

mod args;
mod candidates;

use args::Set;
use candidates::{functions, Function, Imp, C_REFERENCE};
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const BLOCK: u64 = 1 << 16;
const CANONICAL_NAN: u64 = 0x7ff8_0000_0000_0000;

struct Options {
    millions: f64,
    worst_cases: Option<PathBuf>,
    only: Option<Vec<String>>,
    candidates: Option<Vec<String>>,
    digests_out: Option<PathBuf>,
    digests_in: Option<PathBuf>,
    threads: usize,
    examples: usize,
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let (command, rest) = argv.split_first().map(|(c, r)| (c.as_str(), r)).unwrap_or(("", &[]));
    match command {
        "compare" => compare(&options(rest)),
        "time" => time(&options(rest)),
        "dump" => dump(rest),
        _ => {
            eprintln!("usage: libm-compare compare|time|dump ... (see src/main.rs)");
            std::process::exit(2);
        }
    }
}

fn options(rest: &[String]) -> Options {
    let mut o = Options {
        millions: 4.0,
        worst_cases: None,
        only: None,
        candidates: None,
        digests_out: None,
        digests_in: None,
        threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
        examples: 3,
    };
    let mut it = rest.iter();
    while let Some(flag) = it.next() {
        let mut value = || it.next().unwrap_or_else(|| panic!("{flag} needs a value")).clone();
        let list = |s: String| Some(s.split(',').map(str::to_owned).collect());
        match flag.as_str() {
            "--millions" => o.millions = value().parse().expect("--millions"),
            "--worst-cases" => o.worst_cases = Some(value().into()),
            "--only" => o.only = list(value()),
            "--candidates" => o.candidates = list(value()),
            "--digests-out" => o.digests_out = Some(value().into()),
            "--digests-in" => o.digests_in = Some(value().into()),
            "--threads" => o.threads = value().parse().expect("--threads"),
            "--examples" => o.examples = value().parse().expect("--examples"),
            other => panic!("unknown option {other}"),
        }
    }
    if cfg!(target_family = "wasm") {
        o.threads = 1;
    }
    o
}

fn selected(o: &Options) -> Vec<Function> {
    let mut fs = functions();
    if let Some(only) = &o.only {
        fs.retain(|f| only.iter().any(|n| n == f.name));
    }
    if let Some(cands) = &o.candidates {
        for f in &mut fs {
            f.imps
                .retain(|(n, _)| *n == C_REFERENCE || cands.iter().any(|c| c == n));
        }
    }
    fs
}

/// NaN-insensitive bits, so a digest is the same whatever payload a NaN carries
#[inline(always)]
fn canonical(x: f64) -> u64 {
    if x.is_nan() {
        CANONICAL_NAN
    } else {
        x.to_bits()
    }
}

fn is_snan(x: f64) -> bool {
    x.is_nan() && x.to_bits() & (1 << 51) == 0
}

#[inline(always)]
fn fnv(mut h: u64, v: u64) -> u64 {
    for b in v.to_le_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;

#[derive(Default, Clone)]
struct Tally {
    /// arguments where the result differs (not both NaN)
    mismatches: u64,
    /// arguments where both are NaN with different bits
    nan_bits: u64,
    /// mismatches with a signalling NaN among the arguments: no arithmetic yields one (IEEE 754
    /// operations return quiet NaNs), so yFoil can only meet them by building the bits itself
    snan_input: u64,
    examples: Vec<String>,
}

struct BlockResult {
    /// digest per implementation (index as in `Function::imps`)
    digests: Vec<u64>,
    /// tallies against the reference (index 0) per implementation
    tallies: Vec<Tally>,
}

fn run_block(f: &Function, set: &Set, block: u64, examples: usize, have_c: bool) -> BlockResult {
    let lo = block * BLOCK;
    let hi = (lo + BLOCK).min(set.len());
    let k = f.imps.len();
    let mut digests = vec![FNV_OFFSET; k];
    let mut tallies = vec![Tally::default(); k];
    for i in lo..hi {
        let a = set.arg(i);
        let reference = f.imps[0].1.eval(a);
        for (j, (_, imp)) in f.imps.iter().enumerate() {
            let got = if j == 0 { reference } else { imp.eval(a) };
            digests[j] = fnv(fnv(digests[j], canonical(got[0])), canonical(got[1]));
            if j == 0 || !have_c {
                continue;
            }
            let mut real = false;
            let mut nan = false;
            for c in 0..2 {
                if got[c].to_bits() != reference[c].to_bits() {
                    if got[c].is_nan() && reference[c].is_nan() {
                        nan = true;
                    } else {
                        real = true;
                    }
                }
            }
            let t = &mut tallies[j];
            let snan = a
                .iter()
                .take(if f.arity == candidates::Arity::Binary { 2 } else { 1 })
                .any(|x| is_snan(*x));
            if real && snan {
                t.snan_input += 1;
            } else if real {
                t.mismatches += 1;
                if t.examples.len() < examples {
                    t.examples.push(example(f, a, reference, got));
                }
            } else if nan {
                t.nan_bits += 1;
            }
        }
    }
    BlockResult { digests, tallies }
}

fn example(f: &Function, a: [f64; 2], want: [f64; 2], got: [f64; 2]) -> String {
    let args = if f.arity == candidates::Arity::Binary {
        format!(
            "({:e} [{:#018x}], {:e} [{:#018x}])",
            a[0],
            a[0].to_bits(),
            a[1],
            a[1].to_bits()
        )
    } else {
        format!("({:e} [{:#018x}])", a[0], a[0].to_bits())
    };
    let n = if f.name == "sin_cos" { 2 } else { 1 };
    let mut s = format!("{}{args}:", f.name);
    for c in 0..n {
        let ulps = (want[c].to_bits() as i64).wrapping_sub(got[c].to_bits() as i64);
        let _ = write!(
            s,
            " C {:e} [{:#018x}] got {:e} [{:#018x}] ({ulps:+} ulp)",
            want[c],
            want[c].to_bits(),
            got[c],
            got[c].to_bits()
        );
    }
    s
}

fn run_set(f: &Function, set: &Set, o: &Options, have_c: bool) -> Vec<BlockResult> {
    let blocks = set.len().div_ceil(BLOCK);
    let next = AtomicU64::new(0);
    let mut results: Vec<Option<BlockResult>> = (0..blocks).map(|_| None).collect();
    let slots: Vec<std::sync::Mutex<&mut Option<BlockResult>>> =
        results.iter_mut().map(std::sync::Mutex::new).collect();
    let work = || loop {
        let b = next.fetch_add(1, Ordering::Relaxed);
        if b >= blocks {
            break;
        }
        let r = run_block(f, set, b, o.examples, have_c);
        **slots[b as usize].lock().unwrap() = Some(r);
    };
    if o.threads <= 1 {
        work();
    } else {
        std::thread::scope(|s| {
            for _ in 0..o.threads {
                s.spawn(work);
            }
        });
    }
    drop(slots);
    results.into_iter().map(|r| r.expect("every block ran")).collect()
}

/// (function, set, block or `all`, implementation) -> digest
type Digests = BTreeMap<(String, String, String, String), u64>;

fn read_digests(path: &PathBuf) -> Digests {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let key = (f[0].to_owned(), f[1].to_owned(), f[2].to_owned(), f[3].to_owned());
            (key, u64::from_str_radix(f[4], 16).expect("digest"))
        })
        .collect()
}

/// The blocks whose digest differs from the recorded one of `recorded` (None: nothing recorded)
fn differing_blocks(
    rd: &Digests,
    f: &Function,
    set: &Set,
    blocks: &[BlockResult],
    j: usize,
    recorded: &str,
) -> Option<Vec<usize>> {
    let key = |b: &str| (f.name.to_owned(), set.name.clone(), b.to_owned(), recorded.to_owned());
    let all = blocks.iter().fold(FNV_OFFSET, |h, b| fnv(h, b.digests[j]));
    let d = *rd.get(&key("all"))?;
    if d == all {
        return Some(Vec::new());
    }
    Some(
        blocks
            .iter()
            .enumerate()
            .filter(|(b, r)| rd.get(&key(&b.to_string())) != Some(&r.digests[j]))
            .map(|(b, _)| b)
            .collect(),
    )
}

fn compare(o: &Options) {
    let have_c = cfg!(feature = "c-reference");
    let recorded = o.digests_in.as_ref().map(read_digests);
    if !have_c && recorded.is_none() {
        panic!("built without c-reference: pass --digests-in FILE (from a native run's --digests-out)");
    }
    let mut out_digests =
        String::from("# function\tset\tblock\timplementation\tdigest (FNV-1a over the result bits, NaNs canonical)\n");
    let n = (o.millions * 1e6) as u64;
    println!(
        "target {} | {} thread(s) | {} M random arguments per set | reference {}",
        target(),
        o.threads,
        o.millions,
        if have_c {
            "core-math (C) via cr-libm"
        } else {
            "recorded digests"
        }
    );
    // (function, candidate) -> [args, mismatches, NaN-only, sNaN input, blocks != C, blocks != recorded self]
    let mut summary: BTreeMap<(usize, String), [u64; 6]> = BTreeMap::new();
    let fs = selected(o);
    for (fi, f) in fs.iter().enumerate() {
        let sets = args::sets(f.name, f.arity, n, f.worst_cases, o.worst_cases.as_deref());
        for set in &sets {
            let t0 = Instant::now();
            let blocks = run_set(f, set, o, have_c);
            let mut line = format!("{:8} {:32} {:>10}", f.name, set.name, set.len());
            for (j, (cand, _)) in f.imps.iter().enumerate() {
                let all = blocks.iter().fold(FNV_OFFSET, |h, b| fnv(h, b.digests[j]));
                let _ = writeln!(out_digests, "{}\t{}\tall\t{cand}\t{all:016x}", f.name, set.name);
                for (b, r) in blocks.iter().enumerate() {
                    let _ = writeln!(
                        out_digests,
                        "{}\t{}\t{b}\t{cand}\t{:016x}",
                        f.name, set.name, r.digests[j]
                    );
                }
                if *cand == C_REFERENCE {
                    continue;
                }
                let s = summary.entry((fi, cand.to_string())).or_default();
                s[0] += set.len();
                let _ = write!(line, " | {cand}");
                if have_c {
                    let mism: u64 = blocks.iter().map(|b| b.tallies[j].mismatches).sum();
                    let nan: u64 = blocks.iter().map(|b| b.tallies[j].nan_bits).sum();
                    let snan: u64 = blocks.iter().map(|b| b.tallies[j].snan_input).sum();
                    s[1] += mism;
                    s[2] += nan;
                    s[3] += snan;
                    let _ = write!(line, " {mism}");
                    if nan > 0 {
                        let _ = write!(line, " (+{nan} NaN bits)");
                    }
                    if snan > 0 {
                        let _ = write!(line, " (+{snan} sNaN input)");
                    }
                    for e in blocks.iter().flat_map(|b| &b.tallies[j].examples).take(o.examples) {
                        eprintln!("  {cand} {}: {e}", set.name);
                    }
                }
                if let Some(rd) = &recorded {
                    for (k, who, label) in [(4, C_REFERENCE, "C"), (5, *cand, "recorded")] {
                        match differing_blocks(rd, f, set, &blocks, j, who) {
                            None => {
                                let _ = write!(line, " no-{label}-digest");
                            }
                            Some(bad) if bad.is_empty() => {
                                let _ = write!(line, " digest={label}");
                            }
                            Some(bad) => {
                                s[k] += bad.len() as u64;
                                let _ = write!(
                                    line,
                                    " DIGEST!={label} in {} block(s) {:?}",
                                    bad.len(),
                                    &bad[..bad.len().min(6)]
                                );
                            }
                        }
                    }
                }
            }
            println!("{line}  [{:.1}s]", t0.elapsed().as_secs_f64());
        }
    }
    // without the C the per-argument counts are not known; without recorded digests the block counts
    println!("\nsummary: arguments | mismatches | sNaN-input mismatches | NaN-vs-NaN bit differences | digest blocks != C | != recorded self");
    for ((fi, cand), s) in &summary {
        let count = |known: bool, v: u64| if known { v.to_string() } else { "-".into() };
        println!(
            "  {:8} {:12} {:>11} args {:>9} mismatches {:>5} sNaN-input {:>6} NaN-bits {:>4} blocks!=C {:>4} blocks!=self",
            fs[*fi].name,
            cand,
            s[0],
            count(have_c, s[1]),
            count(have_c, s[3]),
            count(have_c, s[2]),
            count(recorded.is_some(), s[4]),
            count(recorded.is_some(), s[5])
        );
    }
    if let Some(path) = &o.digests_out {
        std::fs::write(path, out_digests).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        println!("wrote {}", path.display());
    }
}

/// Rough ns/call (throughput) over 2^20 arguments from each function's first `value[...]` set, its
/// typical range
fn time(o: &Options) {
    println!(
        "target {} | ns/call, best of 5 passes over 2^20 arguments from the typical-range set",
        target()
    );
    for f in selected(o) {
        let sets = args::sets(f.name, f.arity, 1 << 20, &[], None);
        let set = sets
            .iter()
            .find(|s| s.name.starts_with("value") || s.name.contains(":value"))
            .unwrap_or(&sets[0]);
        let a: Vec<[f64; 2]> = (0..set.len()).map(|i| set.arg(i)).collect();
        let mut line = format!("{:8} {:30}", f.name, set.name);
        for (cand, imp) in &f.imps {
            let mut best = f64::INFINITY;
            for _ in 0..5 {
                let t0 = Instant::now();
                let mut acc = 0u64;
                for &x in &a {
                    let r = imp.eval(std::hint::black_box(x));
                    acc ^= r[0].to_bits() ^ r[1].to_bits();
                }
                std::hint::black_box(acc);
                best = best.min(t0.elapsed().as_secs_f64() * 1e9 / a.len() as f64);
            }
            let _ = write!(line, " | {cand} {best:.1}");
        }
        println!("{line}");
    }
}

/// Every argument and result bit pattern of one block, for differencing a target without the C
/// against a native run
fn dump(rest: &[String]) {
    let [function, set_name, block, candidate, rest @ ..] = rest else {
        panic!("dump FUNCTION SET BLOCK CANDIDATE [--millions N] [--worst-cases DIR]");
    };
    let o = options(rest);
    let block: u64 = block.parse().expect("block");
    let f = functions().into_iter().find(|f| f.name == function).expect("function");
    let imp: Imp = f.imps.iter().find(|(n, _)| n == candidate).expect("candidate").1;
    let sets = args::sets(
        f.name,
        f.arity,
        (o.millions * 1e6) as u64,
        f.worst_cases,
        o.worst_cases.as_deref(),
    );
    let set = sets.iter().find(|s| &s.name == set_name).expect("set");
    for i in block * BLOCK..((block + 1) * BLOCK).min(set.len()) {
        let a = set.arg(i);
        let r = imp.eval(a);
        println!(
            "{i} {:016x} {:016x} {:016x} {:016x}",
            a[0].to_bits(),
            a[1].to_bits(),
            canonical(r[0]),
            canonical(r[1])
        );
    }
}

fn target() -> String {
    // aarch64 has fused multiply-add in its base ISA; x86-64 has it only from -Ctarget-cpu=haswell
    // (x86-64-v3) up; wasm32 has none, so `f64::mul_add` is compiler-builtins' software fma
    let fma = if cfg!(any(target_arch = "aarch64", target_feature = "fma")) {
        "hardware fma"
    } else {
        "no hardware fma"
    };
    format!("{}-{} ({fma})", std::env::consts::ARCH, std::env::consts::OS)
}
