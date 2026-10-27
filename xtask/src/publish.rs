//! `cargo xtask publish <study> [--run <UTC>|latest] [--dry-run]`: copy one run of a study into
//! the docs.
//!
//! A study never writes under `docs/`. Each run leaves a bundle in its run folder,
//! `studies/<study>/runs/<UTC>/publish/{tables,figures}/` (`studies/study-support/publish.rs`), and
//! this command replaces `docs/validation/<study>/{tables,figures}/` with the chosen run's bundle.
//! It never writes a page: the pages of `docs/validation/<study>/` are written by hand, and only
//! include the tables (`--8<-- "validation/<study>/tables/<name>.md"`) and embed the figures
//! (`<ironlab-figure src="figures/<name>.fig">`).
//!
//! Before copying, every reference the study's pages make into its `tables/` or `figures/` is
//! checked against the bundle, and the publish is refused if one is missing (Zensical would fail on
//! a missing snippet, but render a missing figure as an empty box). Products of the run that no
//! page uses yet are listed, so a new table or figure is noticed rather than lost.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

/// The studies that publish into `docs/validation/<study>/`, each from `studies/<study>/runs/`.
const STUDIES: [&str; 6] = [
    "aerofoil-series",
    "series-cases",
    "known-issues",
    "branch-coverage",
    "tolerances",
    "conditioning",
];

pub fn run(args: &[String]) {
    let usage = || -> ! {
        eprintln!(
            "usage: cargo xtask publish <{}> [--run <UTC>|latest] [--dry-run]",
            STUDIES.join("|")
        );
        std::process::exit(2);
    };
    let Some(study) = args.first().filter(|s| STUDIES.contains(&s.as_str())) else {
        usage()
    };
    let mut run = "latest".to_string();
    let mut dry_run = false;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--run" => {
                run = args.get(i + 1).cloned().unwrap_or_else(|| usage());
                i += 1;
            }
            "--dry-run" => dry_run = true,
            _ => usage(),
        }
        i += 1;
    }

    let root = crate::root();
    let runs = root.join("studies").join(study).join("runs");
    let run_dir = if run == "latest" {
        latest_run(&runs).unwrap_or_else(|| {
            eprintln!("no run of {study} with a publish/ folder under {}", runs.display());
            std::process::exit(1);
        })
    } else {
        runs.join(&run)
    };
    let bundle = run_dir.join("publish");
    if !bundle.is_dir() {
        eprintln!(
            "{} has no publish/ folder: redraw it with the study's --render <run_dir>",
            run_dir.display()
        );
        std::process::exit(1);
    }
    let docs = root.join("docs/validation").join(study);
    println!(
        "publishing {} → {}",
        bundle.strip_prefix(&root).unwrap_or(&bundle).display(),
        docs.strip_prefix(&root).unwrap().display()
    );

    // what the run offers, as paths relative to the study's docs directory
    let products: BTreeSet<String> = ["tables", "figures"]
        .iter()
        .flat_map(|d| files(&bundle.join(d)).into_iter().map(move |f| format!("{d}/{f}")))
        .collect();

    // what the pages use
    let mut missing = vec![];
    let mut used = BTreeSet::new();
    for page in pages(&docs) {
        let text = fs::read_to_string(&page).unwrap();
        let rel_page = page.strip_prefix(&root).unwrap().display().to_string();
        for r in references(&root, &page, &text) {
            match r.strip_prefix(&docs) {
                Ok(inside) if inside.starts_with("tables") || inside.starts_with("figures") => {
                    let key = inside.display().to_string();
                    if products.contains(&key) {
                        used.insert(key);
                    } else {
                        missing.push(format!("{rel_page}: {key}"));
                    }
                }
                // a reference elsewhere in the docs must simply exist
                _ => {
                    if !r.exists() {
                        missing.push(format!("{rel_page}: {}", r.display()));
                    }
                }
            }
        }
    }
    if !missing.is_empty() {
        eprintln!("refused: the pages use files this run does not provide:");
        for m in &missing {
            eprintln!("  {m}");
        }
        std::process::exit(1);
    }

    // what changes
    let mut added = 0;
    let mut changed = 0;
    let mut removed = 0;
    for p in &products {
        match fs::read(docs.join(p)) {
            Ok(old) if old == fs::read(bundle.join(p)).unwrap() => {}
            Ok(_) => changed += 1,
            Err(_) => added += 1,
        }
    }
    for d in ["tables", "figures"] {
        for f in files(&docs.join(d)) {
            if !products.contains(&format!("{d}/{f}")) {
                removed += 1;
            }
        }
    }
    println!(
        "  {} files: {added} new, {changed} changed, {removed} removed",
        products.len()
    );
    let unused: Vec<&String> = products.iter().filter(|p| !used.contains(*p)).collect();
    if !unused.is_empty() {
        println!("  not used by any page yet:");
        for u in unused {
            println!("    {u}");
        }
    }
    if dry_run {
        println!("  dry run: nothing written");
        return;
    }
    // write in place: overwrite only what changed and remove only what the run no longer has. An
    // editor that updates links on file moves reads a deleted-then-recreated file as moved, and
    // repoints the pages at the run folder
    for d in ["tables", "figures"] {
        let target = docs.join(d);
        fs::create_dir_all(&target).unwrap();
        for f in files(&target) {
            if !products.contains(&format!("{d}/{f}")) {
                fs::remove_file(target.join(&f)).unwrap();
            }
        }
        for f in files(&bundle.join(d)) {
            let new = fs::read(bundle.join(d).join(&f)).unwrap();
            if fs::read(target.join(&f)).ok().as_ref() != Some(&new) {
                fs::write(target.join(&f), new).unwrap();
            }
        }
    }
    println!("  published");
}

/// The newest run folder (UTC stamps sort by time) that holds a `publish/` bundle.
fn latest_run(runs: &Path) -> Option<PathBuf> {
    let mut dirs: Vec<PathBuf> = fs::read_dir(runs)
        .ok()?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.join("publish").is_dir())
        .collect();
    dirs.sort();
    dirs.pop()
}

/// The file names in a directory (not recursive), sorted.
fn files(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = fs::read_dir(dir)
        .map(|rd| {
            rd.flatten()
                .filter(|e| e.path().is_file())
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default();
    v.sort();
    v
}

/// Every hand-written page of the study: the Markdown files of its docs directory, outside its
/// `tables/`.
fn pages(dir: &Path) -> Vec<PathBuf> {
    let mut out = vec![];
    let Ok(rd) = fs::read_dir(dir) else {
        return out;
    };
    for e in rd.flatten() {
        let p = e.path();
        if p.is_dir() {
            if p.file_name().is_some_and(|n| n != "tables" && n != "figures") {
                out.extend(pages(&p));
            }
        } else if p.extension().is_some_and(|x| x == "md") {
            out.push(p);
        }
    }
    out.sort();
    out
}

/// The files a page refers to: snippet includes (relative to `docs/`), ironlab figure embeds and
/// Markdown images (relative to the page). External links are not files and are skipped.
fn references(root: &Path, page: &Path, text: &str) -> Vec<PathBuf> {
    let dir = page.parent().unwrap();
    let mut out = vec![];
    for line in text.lines() {
        let l = line.trim();
        if let Some(rest) = l.strip_prefix("--8<--") {
            let target = rest.trim().trim_matches('"');
            if !target.is_empty() {
                out.push(root.join("docs").join(target));
            }
        }
        for (open, close) in [("<ironlab-figure src=\"", "\""), ("](", ")")] {
            let mut s = line;
            while let Some(at) = s.find(open) {
                s = &s[at + open.len()..];
                let Some(end) = s.find(close) else { break };
                let target = &s[..end];
                s = &s[end..];
                let is_file = open.starts_with('<')
                    || ((target.starts_with("figures/") || target.starts_with("tables/")) && !target.contains("://"));
                if is_file {
                    out.push(normalise(&dir.join(target.split('#').next().unwrap())));
                }
            }
        }
    }
    out
}

/// `a/b/../c` → `a/c`, without touching the file system.
fn normalise(p: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            std::path::Component::ParentDir => {
                out.pop();
            }
            std::path::Component::CurDir => {}
            other => out.push(other),
        }
    }
    out
}
