---
icon: lucide/file-text
---

# Documentation

How the docs site is laid out, and how a study's results reach it.

## The site's structure

The site has three tabs: the **User guide** (using yFoil and background about XFOIL), **Development** (how it is built and
tested) and **Validation** (how it compares with XFOIL, for engineers). The navigation is written out in
`zensical.toml`, so the reading order is deliberate, and the file tree mirrors it:

- each tab and section is a directory with an `index.md` (or, in Development, a `README.md`): a
  tab's is its landing page, and a Validation study's is the page the section's heading opens;
- each page is the kebab-case of its title (`Validation › Series cases › Reynolds number` is
  `docs/validation/series-cases/reynolds-number.md`);
- moving or retitling a page moves its file, and there are no redirects;
- every subject has one home, and no two pages share a name.

## Icons

The sidebar follows the user guide's pattern, and `scripts/docs-check.py` enforces it before every
build and serve (`scripts/docs.sh`), so a page that breaks it fails the build:

- the entries visible when the sidebar opens have an `icon:` in their front matter: pages directly
  under a tab, pages directly under a group header, and the heading page of a collapsible section
  (a known-issues category);
- the heading page of a section directly under a tab, its `index.md` or `README.md`, is shown as the
  section's always-open group header (Validation › Series cases, Development › Testing) and has none;
- pages inside a collapsible section, collapsed by default (the issue pages, the run pages, the
  subroutine pages), have none;
- the home page and each tab's landing page are not in the sidebar and need none.

The check also fails on front matter that opens with `---` and never closes. Icons are Lucide names
(`icon: lucide/<name>`) from the set bundled with the pinned Zensical.

## Validation pages are written by hand

Each study is a section of the Validation tab, and every section has the same shape:

1. **The index page**: what the study is for and why, its main results, and a table of the section's
   pages (*In this section*).
2. **Chart pages**: the figures, each with the explanation of what it shows. When a study has many
   figures they are split over pages by subject (a regime, an issue), each page listed on the index.

The pages are prose written by people. The numbers are not: every table comes from a study run and is
**included** (`--8<-- "validation/<study>/tables/<name>.md"`), and every figure is **embedded** from the
run (`<ironlab-figure src="figures/<name>.fig" …>`). Anything that changes from run to run belongs in
an included table or a figure, never typed into the prose, so publishing a new run cannot leave a page
contradicting itself. Each index page includes `tables/run.md`, so the reader can see which run the
section shows.

## From a run to the docs

A study never writes under `docs/`. Every run writes a dated folder, `studies/<study>/runs/<UTC>/`
(gitignored), and in it a `publish/` bundle laid out as the docs hold it:

```text
studies/<study>/runs/<UTC>/publish/
  tables/run.md     which run this is
  tables/*.md       the table fragments
  figures/*         the .fig files (and the branch-coverage study's SVGs)
```

So any number of runs can be made and compared without touching the site. To show a run in the docs:

```sh
cargo xtask publish <study> --dry-run       # what would change; the newest run by default
cargo xtask publish <study> [--run <UTC>]   # copy it in
```

`publish` replaces `docs/validation/<study>/tables/` and `figures/` with the run's bundle and never
touches a page. Before copying it checks that every table and figure the section's pages use is in the
bundle, and refuses if one is missing. It lists anything in the bundle that no page uses yet, so a new
table or figure is noticed rather than lost. A run folder made before a change to the study's output can
be redrawn without rerunning the study: each study takes `--render <run_dir>`.

The table names are the contract between a study and its pages: renaming one in the study means
editing the pages that include it. A fragment's links are written relative to
`docs/validation/<study>/`, the directory of every page that includes it. The `tables/` directories are
fragments, not pages; `zensical.toml` excludes them from the build.

| Study | Run with | Publish with |
|---|---|---|
| Aerofoil series | `cargo run --release -p aerofoil-series` | `cargo xtask publish aerofoil-series` |
| Series cases | `cargo xtask fixtures --group series`, then `cargo run --release -p series-cases` | `cargo xtask publish series-cases` |
| Known issues | `cargo xtask fixtures --group known-issues`, then `cargo run --release -p known-issues` | `cargo xtask publish known-issues` |
| Branch coverage | `cargo xtask fixtures --group branch-coverage` and `--group non-finite`, then `cargo run --release -p branch-coverage` | `cargo xtask publish branch-coverage` |
| Tolerances | `studies/tolerances/tolerances.sh` | `cargo xtask publish tolerances` |

A series or known-issues case is placed on a page by its `topic` in `xtask/fixtures-config/cases.toml`:
the run tables it as `tables/cases-<topic>.md`, and the page of that name includes the table and embeds
the case's figure. Adding a case to a page means setting its `topic` and adding its figure to the page.

## Development pages are generated

The pages of Development › Testing (coverage, branch tests, subroutines) are written whole by the tools
that measure them (`cargo xtask coverage`, `cargo xtask steps`, `cargo run --bin
generate_subroutine_validation`) and are not edited by hand. Nothing under `docs/` is reference data:
no XFOIL dump, `.dat`, `.pol`, `DUMP` or `CPWR` output is ever committed there (`.gitignore` enforces it).
