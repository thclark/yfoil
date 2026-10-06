//! (b) The route measurement behind `cargo xtask route`, compiled only with `--cfg yfoil_route`
//! and LLVM source coverage (`-C instrument-coverage`). Without the cfg this module is empty.
//!
//! XFOIL's route through a step is measured by `cargo xtask steps` as the difference of the gcov
//! counts of two prefix runs: the script through step (c, i) less the script through the step
//! before. Both prefixes end with VISCAL's end-of-call code, so the difference is the step's own
//! work plus the change in what the end of the call does — and, for the first iteration of a
//! call, the operating-point command and VISCAL's set-up before it. yFoil's route is measured as
//! the same quantity from windows of profile counts, each written by this test to
//! `$YFOIL_ROUTE_OUT/<case>_<call>_<iteration>_<part>.profraw`, with the counters reset immediately
//! before the window:
//!
//! - `A<j>`: the session seeded with XFOIL's state on entry to iteration j, then one iteration
//!   (`solve_viscous(1)`: SETBL … UPDATE, the tail and the end of the call);
//! - `B<j>`: the same seed, then no iteration (`solve_viscous(0)`: only the set-up and the end of
//!   the call, which `A<j>` also runs), so `A<j> − B<j>` is iteration j plus the change in the end
//!   of the call, as XFOIL's difference is;
//! - `W` (the first iteration of a call, and a whole call): the operating-point command and VISCAL
//!   with no iteration, from the end of the previous call — the previous call's last iteration
//!   replayed from XFOIL's state. A call that marches its BL afresh (the first call, from a session
//!   made from the panels, or the first after an `INIT`, applied to the previous call's end)
//!   instead runs its first iteration inside `W`, since SETBL's MRCHUE march is part of it. For an
//!   inviscid case `W` is the whole step: the call's command, after every earlier command from the
//!   panels.
//!
//! A step's route is then `W + Σ (A<j> − B<j>)` over its iterations (from the second, for a call
//! that marches). The session is made from the panels outside the window: XFOIL's `LOAD` is
//! likewise outside its step (`cargo xtask steps` measures the first step from the script before its first
//! operating point). `cargo xtask route` reads
//! the profiles and compares the call count of every translated subroutine with XFOIL's.
//!
//! The steps are listed in the file named by `YFOIL_ROUTE_STEPS`, one `case call iteration` per
//! line (iteration 0: the whole call). A step that cannot be replayed (a panic) is printed as
//! `ROUTE-UNMEASURED case call iteration` and the measurement goes on.

#![cfg(yfoil_route)]

use crate::fixtures::cases::{dump, load, prologue_jogged, seed, Case, Op};
use crate::utilities::tolerances::TOL_CROSS_HOST;
use std::ffi::CString;
use std::os::raw::c_char;
use yfoil::geometry::{panel_foil, read_geometry_from_file};
use yfoil::solver::analysis::Session;
use yfoil::solver::specal::{alpha_command, cl_command, sequence_command};
use yfoil::solver::viscal::solve_viscous;

extern "C" {
    fn __llvm_profile_reset_counters();
    fn __llvm_profile_set_filename(name: *const c_char);
    fn __llvm_profile_write_file() -> i32;
}

/// Run `work` with the counters reset before it and written to `path` after it. The runtime also
/// writes the counters at exit to the last file named, so the name is then pointed elsewhere.
fn window(path: &str, work: impl FnOnce()) {
    let p = CString::new(path).unwrap();
    let away = CString::new(format!(
        "{}/exit.profraw",
        path.rsplit_once('/').map_or(".", |(d, _)| d)
    ))
    .unwrap();
    // SAFETY: the LLVM profile runtime is linked whenever the crate is built with
    // `-C instrument-coverage`, which `cargo xtask route` always does with this cfg
    unsafe { __llvm_profile_reset_counters() };
    work();
    unsafe {
        __llvm_profile_set_filename(p.as_ptr());
        assert_eq!(__llvm_profile_write_file(), 0, "write {path}");
        __llvm_profile_set_filename(away.as_ptr());
    }
}

/// The operating-point command of call `call` (SPECAL or SPECCL).
fn command(c: &Case, session: &mut Session, call: usize) {
    let (st, sys) = session.parts_mut();
    match c.ops[call - 1] {
        Op::Alfa(a) => alpha_command(st, sys, a.to_radians()),
        Op::Aseq(a) => sequence_command(st, sys, a.to_radians()),
        Op::Cl(cl) => {
            cl_command(st, sys, cl);
        }
    }
}

/// The first SETBL call of VISCAL call `call`.
fn first_setbl(c: &Case, call: usize) -> usize {
    c.points[..call - 1]
        .iter()
        .map(|p| p["NITDONE"].parse::<usize>().unwrap())
        .sum::<usize>()
        + 1
}

/// A session seeded with XFOIL's state on entry to SETBL call `k`, of VISCAL call `call`.
fn seeded(c: &Case, call: usize, k: usize) -> Session {
    let d = dump(c, &format!("mrchdu_input_{k}.dat"));
    let alpha = (!matches!(c.ops[call - 1], Op::Cl(_))).then(|| d.real("ALFA"));
    let mut session = prologue_jogged(c, call, alpha, None);
    seed(&mut session, &d, TOL_CROSS_HOST);
    session
}

fn measure_step(c: &Case, call: usize, iteration: usize, out: &str) {
    let tag = format!("{out}/{}_{call}_{iteration}", c.name);
    let wake_length = c.spec.wake_length;
    if c.points.is_empty() {
        // inviscid: every earlier command from the panels, then this one in the window
        let g = read_geometry_from_file(c.dir.join("panels.json").to_str().unwrap()).expect("panels.json");
        let foil = panel_foil(&g);
        let mut session = Session::new(&foil, c.spec.clone());
        for k in 1..call {
            command(c, &mut session, k);
        }
        window(&format!("{tag}_W.profraw"), || command(c, &mut session, call));
        return;
    }
    let first = first_setbl(c, call);
    let nit: usize = c.points[call - 1]["NITDONE"].parse().unwrap();
    // a call that marches its BL afresh (the first, or the first after an INIT) is measured from
    // the start: SETBL's MRCHUE march is in its first iteration, which a seeded replay skips
    let marches = call == 1 || c.inits.contains(&(call - 1));
    let previous = (call > 1).then(|| {
        // the end of the previous call: its last iteration from XFOIL's state
        let mut s = seeded(c, call - 1, first - 1);
        let (st, sys) = s.parts_mut();
        solve_viscous(st, sys.as_mut(), 1, wake_length, None);
        if c.inits.contains(&(call - 1)) {
            s.init();
        }
        s
    });
    // the first call's session is made outside the window, like XFOIL's LOAD (the set-up
    // evaluates the leading edge, as GEOPAR does at LOAD)
    let mut previous = previous.or_else(|| {
        let g = read_geometry_from_file(c.dir.join("panels.json").to_str().unwrap()).expect("panels.json");
        Some(Session::new(&panel_foil(&g), c.spec.clone()))
    });
    if iteration <= 1 {
        let n = if marches { 1 } else { 0 };
        let mut s = previous.take().unwrap();
        window(&format!("{tag}_W.profraw"), || {
            command(c, &mut s, call);
            let (st, sys) = s.parts_mut();
            solve_viscous(st, sys.as_mut(), n, wake_length, None);
        });
    }
    let from = if marches { 2 } else { 1 };
    let iterations = if iteration == 0 {
        from..=nit
    } else {
        iteration.max(from)..=iteration
    };
    for j in iterations {
        let k = first + j - 1;
        for (part, n) in [("A", 1), ("B", 0)] {
            let mut s = seeded(c, call, k);
            window(&format!("{tag}_{part}{j}.profraw"), || {
                let (st, sys) = s.parts_mut();
                solve_viscous(st, sys.as_mut(), n, wake_length, None);
            });
        }
    }
}

#[test]
fn measure() {
    let list = std::env::var("YFOIL_ROUTE_STEPS").expect("YFOIL_ROUTE_STEPS");
    let out = std::env::var("YFOIL_ROUTE_OUT").expect("YFOIL_ROUTE_OUT");
    let mut current: Option<Case> = None;
    for line in std::fs::read_to_string(list)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty())
    {
        let t: Vec<&str> = line.split_whitespace().collect();
        if current.as_ref().is_none_or(|c| c.name != t[0]) {
            current = Some(load(t[0]));
        }
        let c = current.as_ref().unwrap();
        // a step that cannot be replayed is reported, not fatal: its route is then unmeasured
        let run = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            measure_step(c, t[1].parse().unwrap(), t[2].parse().unwrap(), &out)
        }));
        if run.is_err() {
            println!("ROUTE-UNMEASURED {line}");
        }
    }
}
