---
icon: lucide/quote
---

# Citing yFoil

If you use yFoil in any research or work leading to a publication or technical report, please cite
both yFoil *and* XFOIL.

Years of work went into building and validating XFOIL. Authors, maintainers and users of yFoil owe a
great debt of thanks to Prof. Drela, Harold Youngren and any other contributors to the original XFOIL.

## yFoil

!!! note "Stub — to be written"

    yFoil has no paper or DOI yet. Until it does, cite the software and the version you used:

    ```bibtex
    @software{yfoil,
      author  = {Clark, Thomas and {the yFoil contributors}},
      title   = {{yFoil}: {XFOIL} in {Rust}},
      url     = {https://github.com/thclark/yfoil},
      version = {0.1.0},
      license = {GPL-3.0-or-later}
    }
    ```

    Replace this with the paper and archived-release DOI once they exist.

## XFOIL

Cite Drela's paper describing XFOIL ([drela1989](../reference/bibliography.md#drela1989)):

```bibtex
@incollection{drela1989,
  author    = {Drela, Mark},
  title     = {{XFOIL}: An Analysis and Design System for Low {Reynolds} Number Airfoils},
  booktitle = {Low Reynolds Number Aerodynamics},
  editor    = {Mueller, Thomas J.},
  series    = {Lecture Notes in Engineering},
  volume    = {54},
  publisher = {Springer},
  address   = {Berlin},
  pages     = {1--12},
  year      = {1989},
  doi       = {10.1007/978-3-642-84010-4_1}
}
```

## Licence

yFoil is GPL-3.0-or-later, because it is a translation of XFOIL (GPL-2.0-or-later) and carries XFOIL's
source in the repository; see the [licence](https://github.com/thclark/yfoil/blob/main/LICENSE).
