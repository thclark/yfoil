---
icon: lucide/plug
---

# Bindings

!!! warning "Planned, not yet available"

    yFoil is currently a CLI and a Rust library. This page describes where the
    other language bindings are going; nothing here works yet.

The solver core is a plain Rust library with no interactive state and no global
mutable data, which makes it straightforward to expose elsewhere.

## WebAssembly

A WASM target makes yFoil callable from any JavaScript runtime, and — more
interestingly — lets it run entirely in the browser, with no server. That is the
route to interactive documentation: a live polar on this site, computed on the
reader's machine.

## Python

The obvious host for aerofoil work. The intent is a thin, typed wrapper around
the same JSON data model, so that a result is a normal Python object and arrays
arrive as NumPy arrays.

## Julia

Same story, for the same reasons.

## Rust

Already available — add the crate as a dependency and call the solver directly.
The library API mirrors XFOIL's data model exactly inside `solver/` and `bl/`,
which is a deliberate choice: it keeps every translated statement diffable
against its Fortran original.
