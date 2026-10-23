---
icon: lucide/book-open
---

# User guide

## About yFoil

If you're an aerodynamicist, you've almost certainly used XFOIL for 2D aerofoil section analysis. At [OST](https://www.ost.ch/de/), we love XFOIL and use it for building digital twins of wind turbines. It's used by amateurs and professionals alike - for drones, gliders, planes, hydrofoils, turbines, you name it.

So we'd like to introduce [yFoil](https://github.com/thclark/yfoil/), our ground-up reimplementation of XFOIL, featuring:

|| Feature | Benefit |
|---|---|---|
| :white_check_mark: | Platform-independence | *Same results on different machines, always* |
| :white_check_mark: | Robust handling of ill-conditioned cases | *No random spikes post-stall* |
| :white_check_mark: | Error and warning management | *No silent errors or warnings lost in swamps of terminal output* | 
| :white_check_mark: | Memory-safe execution | *No crashes on some cases* | 
| :white_check_mark: | A complete CLI | *Straightforwardly automatable* |
| :white_check_mark: | Published data file schema | *Universally readable inputs-and-outputs in a standard file format* |
| :white_check_mark: | Integration with [ironlab](https://ironlab.org/) | *Interactive and editable output plots, as publication-quality vector graphics* |
| :white_check_mark: | Exact reproduction of XFOIL results | *If your processes already use XFOIL, you can drop yFoil right in* |

!!! note "What does 'Exact reproduction' mean?"

    Across every [well conditioned](../development/conventions/terminology.md#ill-conditioned-solution) and convergent case in the massive [validation data set](../validation/index.md) ( ≈ 2,050 operating points covering every line of code and logical branch across XFOIL's program), the **maximum absolute error is <= 1E-09** between yFoil and XFOIL results (values of CL, CD, CDp, CM and the transition locations)
    

To get started, [install yFoil](get-started/install.md), then [run your first analysis](get-started/your-first-analysis.md).
Ten minutes, and you will have a converged viscous solution and a plot of it, like
these polars drawn with [ironlab](https://ironlab.org) via [`yfoil plot`](commands/plotting.md):

<ironlab-figure src="figures/intro.fig" format="fig" name="intro" alt="NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha">NACA 0012 and NACA 4412 polars at Re 1e6, alpha -15 to 15 degrees: CL against alpha, CL against CD, CM against alpha and CD against alpha</ironlab-figure>


## Scope

yFoil reproduces XFOIL's OPER *analysis* functionality: fixed-alpha and fixed-CL
points, forced transition, polars and Kármán–Tsien compressibility. Inverse
design, interactive mode, multi-element sections and flap hinge moments are out
of scope.

!!! warning "Help us understand what you need!"
    If you need these functions, which we've deliberately placed out of scope to keep the initial project simple, please [raise an issue on GitHub](https://github.com/thclark/yfoil/issues) to let us know what you need!


## Open source license

!!! warning "Copyleft licensing requirements"

    yFoil is licensed as GPL-3.0-or-later for three reasons:

    - the original XFOIL source (licensed GPL-2.0-or-later) is modified and vendored in-tree
    - a translation is considered derivative work even if we hadn't vendored XFOIL
    - unlike non-copyleft licenses like MIT and Apache, GPL variants respect that development is often done for the good of the community with no commercial goal, and requires consuming projects to 'pay it forward'.

|  | In simple terms, this means: |
|---|---|
| :white_check_mark: | You can use the yFoil program executable with no restrictions at all, to make as many analyses as you wish. |
| :white_check_mark: | You can use the yfoil plotting tools and the [ironlab viewer](https://ironlab.org/) the same way. |
| :white_check_mark: | You can export figures and your results to PDF documents and use them in publications with no restrictions. |
| :x: | You must not modify yFoil's source code without *immediately* publishing it under a compatible license. |
| :x: | You must not link an executable to yFoil's underlying library without *immediately* also publishing that executable's source code under a compatible license. |

 If your project is incompatible with GPL license, we suggest either automating analyses via the CLI (which leaves no copyleft requirement), or getting in touch to discuss whether licensing under commercial terms (to support ongoing development) is possible.


## FAIR data and software

!!! note "Stub — to be written"

    yFoil's file formats are designed for Findable, Accessible, Interoperable
    and Reusable practice — one documented JSON schema in and out, full
    round-tripping precision, and no interactive state to reproduce. This
    section will set out how yFoil supports FAIR workflows.


## Related projects

**Fun fact:** there were *SO* many figures generated by successive rounds of validation studies, that the [ironlab](https://ironlab.org) application was
conceived and built JUST so we could manage the amount of figures!! This let us sort and search within studies using figure parameters, then quickly edit the view 
and export at publication quality for the paper. It was genuinely tons of fun to build a plotting library.
