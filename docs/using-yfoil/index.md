---
icon: lucide/book-open
---

# Using yFoil

This tab is for people who want to **run** yFoil: install it, analyse an aerofoil, and read
what comes out. It assumes you know what a polar is, not that you know XFOIL's internals.
Whether yFoil's numbers are fit for engineering work is the next tab,
[Engineering validation](../engineering-validation/index.md); whether it is a faithful
translation of XFOIL is [Fidelity to XFOIL](../fidelity-to-xfoil/index.md).

## Start here

[Install yFoil](installing-yfoil.md), then [run your first analysis](your-first-analysis.md).
Ten minutes, and you will have a converged viscous solution and a plot of it, like
this comparison of two polars drawn with [`yfoil plot`](commands/plotting-results.md) as an
[ironlab](https://github.com/thclark/ironlab) figure (pan, zoom, or click the key to hide a curve):

<ironlab-figure src="index/intro.fig" format="fig" name="intro" alt="NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha">NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha</ironlab-figure>

If you already know XFOIL, the mapping is direct — `yfoil analyse` is one OPER
point, `yfoil polar` is a `PACC`/`ASEQ` sweep — and each command page gives the
equivalent XFOIL session.

## Pages of Using yFoil, in reading order

| Page | Read it when |
| --- | --- |
| [Installing yFoil](installing-yfoil.md) | you have nothing installed yet |
| [Your first analysis](your-first-analysis.md) | you want one worked session, geometry to plotted polar |
| [Commands](commands/index.md) | you need an option: one page per command, each with an example and the XFOIL equivalent |
| [Aerofoil families](aerofoil-families/index.md) | you want to know exactly which section a designation produces, and how it is panelled |
| [Output data](output-data.md) | you are reading the JSON yFoil writes, in a script or another tool |
| [Whether to trust a result](whether-to-trust-a-result.md) | a point came back `Invalid`, or you are about to rely on one that did not |
| [Bibliography](bibliography.md) | you want the source of a definition |

The first two are a tutorial; the rest are reference pages, each complete on its subject, to be
read when the question comes up.

## Scope

yFoil reproduces XFOIL's OPER *analysis* functionality: fixed-alpha and fixed-CL
points, forced transition, polars and Kármán–Tsien compressibility. Inverse
design, interactive mode, multi-element sections and flap hinge moments are out
of scope.

## Open source

!!! note "Stub — to be written"

    yFoil is GPL-3.0-or-later, because XFOIL is vendored in-tree and a
    translation is a derivative work. This section will cover what that means
    for you: what you can do with yFoil, what you owe back, and how to
    contribute ([Contributing](../contributing/index.md) has the method).

## FAIR data and software

!!! note "Stub — to be written"

    yFoil's file formats are designed for Findable, Accessible, Interoperable
    and Reusable practice — one documented JSON schema in and out, full
    round-tripping precision, and no interactive state to reproduce. This
    section will set out how yFoil supports FAIR workflows.

## How to cite yFoil

TODO - add citations for our own papers and Mark Drela's.

Years work went into building and validating XFOIL. Authors, maintainers and
users of yFoil owe a greate debt of thanks to Prof. Drela, Harold Youngren and
any other contributors to the original XFOIL.
