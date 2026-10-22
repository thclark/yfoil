# Building the documentation

This site is built by [Zensical](https://zensical.org), the successor to Material for MkDocs,
from `docs/` and `zensical.toml`. Nothing is installed into the repository:
`scripts/docs.sh` has [`uv`](https://docs.astral.sh/uv/) fetch the pinned Zensical version into
its own cache, and fetches the pinned ironlab-web bundle that draws the embedded figures.

## Previewing and building

``` sh
scripts/docs.sh serve            # live preview on http://localhost:8000
scripts/docs.sh serve --open     # ... and open a browser
scripts/docs.sh build --clean    # the static site into site/
scripts/docs-linkcheck.py        # every internal link of the built site resolves
```

`build` also stamps every reference to the site's own stylesheets, scripts and images with a
content hash (`scripts/docs-cachebust.py`), since Zensical does not fingerprint them. A build
should report no warnings, and the link check no broken links.

## How the site is laid out

The navigation is written out in `zensical.toml`, so the order is a deliberate reading order.
The file tree mirrors it:

- each tab and section is a directory, kebab-case of its label, with an `index.md` that orients
  the reader: what it is for, who it is for, and its pages in reading order;
- each page is the kebab-case of its title (*Using yFoil › Commands › Sweeping a polar* is
  `docs/using-yfoil/commands/sweeping-a-polar.md`);
- a page's own figures go in a directory beside it named after it
  (`sweeping-a-polar.md` + `sweeping-a-polar/naca4412_polar.fig`; an `index.md`'s in `index/`).

Moving or retitling a page moves its file, and every link to it is updated; there are no
redirects. Every heading is unique across the site, and every subject has one home that the others
link to.

## Figures in a page

A figure is an `<ironlab-figure>` element pointing at a `.fig` file beside the page, written by
`yfoil::output::figures::embed_html` in the generators; the element is defined by the ironlab-web
bundle, which `overrides/main.html` loads ([figures](conventions/figures.md#embedding)).
