#!/usr/bin/env python3
"""Report internal links in the built documentation site that do not resolve.

Zensical does not (yet) validate cross-page links, and several of the docs were
written as repository files rather than site pages, so they link to paths like
``../../tests/fixtures/`` that exist in the checkout but not in the site.

Usage:

    scripts/docs.sh build --clean
    scripts/docs-linkcheck.py [SITE_DIR]

Exits non-zero if any link is broken, so it can be used as a gate later. Links
to external hosts, anchors and mailto: are not checked. Root-absolute links are
resolved against the site root, which is how they will behave once served.
"""

from __future__ import annotations

import os
import re
import sys
from collections import defaultdict
from urllib.parse import unquote, urldefrag

HREF = re.compile(r'(?:href|src)="([^"]+)"')
SKIP_PREFIXES = ("http://", "https://", "//", "#", "mailto:", "data:", "javascript:")


def resolve(site_dir: str, page: str, link: str) -> str:
    """Return the filesystem path a link points at, as a directory index if needed."""
    target = unquote(urldefrag(link)[0])
    if target.startswith("/"):
        full = os.path.join(site_dir, target.lstrip("/"))
    else:
        full = os.path.join(os.path.dirname(page), target)
    full = os.path.normpath(full)
    if os.path.isdir(full):
        full = os.path.join(full, "index.html")
    return full


def main(argv: list[str]) -> int:
    site_dir = argv[1] if len(argv) > 1 else "site"
    if not os.path.isdir(site_dir):
        print(f"error: {site_dir}/ does not exist — run scripts/docs.sh build first", file=sys.stderr)
        return 2

    pages = [
        os.path.join(root, name)
        for root, _, names in os.walk(site_dir)
        for name in names
        if name.endswith(".html")
    ]

    broken: dict[str, list[str]] = defaultdict(list)
    for page in sorted(pages):
        with open(page, encoding="utf-8") as handle:
            html = handle.read()
        for link in HREF.findall(html):
            if link.startswith(SKIP_PREFIXES) or not urldefrag(link)[0]:
                continue
            if not os.path.exists(resolve(site_dir, page, link)):
                broken[os.path.relpath(page, site_dir)].append(link)

    if not broken:
        print(f"{len(pages)} pages checked, no broken internal links")
        return 0

    total = sum(len(links) for links in broken.values())
    for page, links in sorted(broken.items()):
        print(page)
        for link in sorted(set(links)):
            print(f"    {link}")
    print(f"\n{total} broken internal link(s) across {len(broken)} of {len(pages)} pages")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
