//! A study run's docs bundle: `runs/<UTC>/publish/`, the run's tables and figures laid out as
//! `docs/validation/<study>/` holds them.
//!
//! A study never writes under `docs/`. Every run writes its bundle, and `cargo xtask publish
//! <study> [--run <UTC>]` copies a chosen run's bundle into the docs, where hand-written pages
//! include the tables (`--8<-- "validation/<study>/tables/<name>.md"`) and embed the figures
//! (`figures/<name>.fig`). The file names are the contract between a study and its pages: rename
//! one and the pages that use it fail to build. A fragment is included into a page of
//! `docs/validation/<study>/`, so a link it carries is relative to that directory
//! (`../../user-guide/...`).
//!
//! ```text
//! publish/
//!   tables/run.md     which run this is (included at the top of the study's index page)
//!   tables/*.md       Markdown table fragments
//!   figures/*         ironlab `.fig`
//! ```

use std::fs;
use std::path::{Path, PathBuf};

pub struct Bundle {
    root: PathBuf,
}

impl Bundle {
    /// An empty bundle in `run_dir/publish/`; a bundle already there (an earlier render of the
    /// same run) is replaced.
    pub fn create(run_dir: &Path) -> Bundle {
        let root = run_dir.join("publish");
        if root.exists() {
            fs::remove_dir_all(&root).expect("remove the old publish/");
        }
        fs::create_dir_all(root.join("tables")).expect("create publish/tables");
        fs::create_dir_all(root.join("figures")).expect("create publish/figures");
        Bundle { root }
    }

    /// `tables/run.md`: which run the published tables and figures come from (`host` empty when
    /// the run did not record it).
    pub fn provenance(&self, command: &str, stamp: &str, git_describe: &str, git_dirty: bool, host: &str) {
        self.table(
            "run",
            &format!(
                "Results from the run of `{command}` at {stamp} (git {git_describe}{}{}).\n",
                if git_dirty { ", dirty" } else { "" },
                if host.is_empty() {
                    String::new()
                } else {
                    format!(", host {host}")
                }
            ),
        );
    }

    /// `tables/<name>.md`.
    pub fn table(&self, name: &str, markdown: &str) {
        fs::write(self.root.join("tables").join(format!("{name}.md")), markdown).expect("write a table fragment");
    }

    /// Copy a figure the run drew into `figures/`, under the same file name. A figure the run did
    /// not draw (one that did not validate) is left out, and the pages that embed it fail the
    /// publish check.
    pub fn figure(&self, path: &Path) {
        if path.exists() {
            fs::copy(path, self.root.join("figures").join(path.file_name().unwrap())).expect("copy a figure");
        }
    }
}
