"""Check the docs pages against the site's sidebar rules before Zensical builds them.

Run by scripts/docs.sh (so by every `build` and `serve`, and by the docs workflow); exits non-zero,
naming each page and the rule it breaks.

1. Front matter is well formed: a page that opens with `---` closes it with a `---` line.
2. Icons (docs/development/conventions/documentation.md, "Icons"), by the page's place in the
   `nav` of zensical.toml:
   - the entries visible when the sidebar opens have an `icon:`: pages directly under a tab, pages
     directly under a group header, and the heading page of a collapsible section;
   - the heading page of a section directly under a tab (its `index.md` or `README.md`, which
     `navigation.indexes` makes the section's heading and `navigation.sections` renders as an
     always-open group header, like Validation › Series cases) has none;
   - pages inside a collapsible section (a section nested below a group header, collapsed by
     default, like a known-issues category) have none, and such a section starts with its heading
     page (`index.md` or `README.md`), the only place its icon can come from;
   - the site's home page and each tab's landing page are not in the sidebar and are exempt.
"""

import pathlib
import re
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"


def front_matter(path):
    """(well_formed, has_icon) for a page."""
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        return True, False
    end = re.search(r"^---\s*$", text[4:], re.M)
    if end is None:
        return False, False
    block = text[4 : 4 + end.start()]
    return True, re.search(r"^icon:\s*\S", block, re.M) is not None


def walk(entries, depth, rules):
    """Collect (page, rule) for each page of a nav list. depth 1 = a tab's entries."""
    for i, entry in enumerate(entries):
        if isinstance(entry, str):
            label, value = None, entry
        else:
            ((label, value),) = entry.items()
        if isinstance(value, list):
            # a collapsible section (below a group header) shows its heading page's icon, so it
            # must have a heading page
            first = value[0] if value else ""
            first = first if isinstance(first, str) else next(iter(first.values()))
            if depth >= 2 and not (isinstance(first, str) and first.endswith(("/index.md", "/README.md"))):
                rules.append((f"section '{label}'", "needs-heading"))
            walk(value, depth + 1, rules)
            continue
        # Zensical, like MkDocs, takes a directory's README.md as its index page
        is_heading = i == 0 and value.endswith(("/index.md", "/README.md"))
        if depth == 1:
            # directly under a tab: the tab's landing page is not in the sidebar
            rules.append((value, "exempt" if is_heading else "icon"))
        elif depth == 2:
            # a section directly under a tab is a group header, always open: its heading page has no
            # icon, the pages under it do
            rules.append((value, "no-icon" if is_heading else "icon"))
        else:
            # a deeper section is collapsible and collapsed by default: its heading has an icon, the
            # pages inside it do not
            rules.append((value, "icon" if is_heading else "no-icon"))


def main():
    nav = tomllib.loads((ROOT / "zensical.toml").read_text(encoding="utf-8"))["project"]["nav"]
    rules = []
    for tab in nav:
        ((_, value),) = tab.items()
        if isinstance(value, list):
            walk(value, 1, rules)
    errors = []
    for page in sorted(DOCS.rglob("*.md")):
        if "/tables/" in page.as_posix():
            continue
        ok, _ = front_matter(page)
        if not ok:
            errors.append(f"{page.relative_to(DOCS)}: front matter opened with '---' is never closed")
    for rel, rule in rules:
        if rule == "needs-heading":
            errors.append(f"{rel}: a collapsible section needs an index.md first, to carry its icon")
            continue
        path = DOCS / rel
        if not path.exists():
            errors.append(f"{rel}: in the nav but missing")
            continue
        ok, has_icon = front_matter(path)
        if rule == "icon" and not has_icon:
            errors.append(f"{rel}: a sidebar page needs an `icon:` in its front matter")
        if rule == "no-icon" and has_icon:
            errors.append(f"{rel}: no icon here (a group-header heading page, or a page inside a collapsible section)")
    if errors:
        print("docs check failed:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
