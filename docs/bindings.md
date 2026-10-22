---
icon: lucide/plug
---

# Bindings

!!! warning "WebAssembly first, the rest planned"

    yFoil is a CLI and a Rust library, with a first WebAssembly binding. The
    Python and Julia bindings described below do not exist yet.

The solver core is a plain Rust library with no interactive state and no global
mutable data, which makes it straightforward to expose elsewhere.

## WebAssembly

A WASM target makes yFoil callable from any JavaScript runtime, and — more
interestingly — lets it run entirely in the browser, with no server. That is the
route to interactive documentation: a live polar on this site, computed on the
reader's machine.

The binding is the `yfoil-wasm` crate in `bindings/wasm/`. It has three calls,
each taking and returning JSON text in the shapes the CLI reads and writes, so a
result computed in the browser is the file `yfoil` would have written for the
same inputs, byte for byte:

| Call | Returns | CLI equivalent |
|---|---|---|
| `naca(designation, panelling, a)` | geometry JSON | `yfoil geometry naca <designation> -o x.json` |
| `analyse(geometry, options)` | analysis record | `yfoil analyse <geometry> -o x.json` |
| `polar(geometry, options)` | polar record | `yfoil polar <geometry> --json -o x.json` |

`options` is a JSON object whose keys are the command's long flags with `_` for
`-` (`{"alpha": 2, "reynolds": 3e5, "include_lagged_closures": true}`), plus
`name`, the foil name the CLI takes from the geometry file's stem. An absent key
takes the CLI's default and an unknown key is an error. `panelling` is a
panelling record, or `""` for the default; `a` is the 6-series `--a`, or `NaN`.
A fourth call, `evaluate(function, x, y, z)`, applies one of the maths library's
functions elementwise; it is there for the bit-identity check below.

```js
import init, { naca, analyse } from "./yfoil_wasm.js";
await init();
const geometry = naca("4412", "", NaN);
const record = JSON.parse(analyse(geometry, JSON.stringify({ alpha: 4, name: "naca4412" })));
```

### Building it

```bash
scripts/wasm.sh build            # target/wasm-pkg/web/ and target/wasm-pkg/nodejs/
scripts/wasm.sh check [--quick]  # build, then the native-against-wasm bit-identity check
```

This needs the `wasm32-unknown-unknown` Rust target, a clang with the
WebAssembly backend and its `llvm-ar` (Apple's clang has no such backend:
`brew install llvm`; on Linux the distribution's `clang` and `llvm` packages
have it; or set `WASM_CC`), and Node for the check. `wasm-bindgen-cli` must be
the exact version `bindings/wasm/Cargo.toml` pins; the script installs it into
`target/wasm-tools/` on first use. `wasm-opt` (binaryen) is used if it is
installed.

The C sources and their flags are the native build's (CLAUDE.md, one vendored
copy and one set of flags). WebAssembly has no C library headers, and wasi-libc's
`fenv.h` lacks the rounding-mode and flag macros CORE-MATH names, so for wasm32
`cr-libm/build.rs` adds one include path, `cr-libm/wasm-include/`: six small
headers stating what WebAssembly is — one rounding mode, round to nearest even,
and no floating-point exception flags. No returned value depends on them:
`fegetround()` is round-to-nearest natively too, and CORE-MATH uses the flag
functions only to set, save and restore the flags. CORE-MATH's `__builtin_fma`
becomes a call to `fma`, since WebAssembly has no fused multiply-add instruction,
resolved from Rust's compiler-builtins: a software implementation, correctly
rounded like the hardware's. The main crate needs no change: everything it
depends on compiles for wasm32, and the file I/O in `geometry::io` is simply not
called.

### Numerical identity

`scripts/wasm.sh check` compares bytes, never values within a tolerance. The
records are compared as text: serde_json writes each f64 as its shortest
round-trip text, so equal text is equal bits, and it writes a NaN as `null`. The
maths library is compared as raw bits, with one rule for NaNs: any NaN equals any
NaN, a NaN against a number is a difference, and every other value, signed zeros
and infinities included, must match bit for bit. IEEE 754 and the WebAssembly
specification leave a NaN's sign and payload unspecified, x86-64 and ARM native
builds already differ in them, and yFoil only ever tests whether a value is a NaN
and serialises it as `null`. The check counts the NaNs whose bits differ, as
information.

Measured 2026-10-22, 227 comparisons, every one identical:

- the geometry of eight sections, one per generator family (`0012`, `4412`,
  `23012`, `23112`, `0012-34`, `16-212`, `63-415`, `64A010`);
- five analyses of every tracked fixture's panels (41 panel sets, N = 60 to 160):
  α = 2° with the lagged closures, fixed CL 0.5 at Re 3e5 and M 0.3, α = 6° at
  Re 1e4, inviscid at M 0.6, and a −10° → 20° polar at ITER 50 embedding every
  point's full distributions — 1129 solved points, every per-station BL array of
  each, past stall (917 valid, 212 invalid and reported with `--allow-invalid`;
  142 more not attempted after a sweep halted);
- each `cr_libm` function and `fma` on 262 144 arguments (random bit patterns and
  their domains) plus 26 edge values (every pair of them for `pow`, `atan2`,
  `hypot` and `fma`).

| Host | Records | NaNs with other bits |
|---|---|---|
| aarch64 macOS (Node 22) | identical | 0 |
| aarch64 Linux (CI runner, gcc native) | identical | 0 |
| x86-64 Linux (CI runner; also Ubuntu 24.04 emulated) | identical | 97 990 of 3.7 million results |

On x86-64 the NaNs come from invalid operations (`ln` and `log10` of a negative
number, `asin(2)`, `cos(∞)`): the native build writes them as `0xfff8…`, the
module as `0x7ff8…`. The CLI itself, built for `wasm32-wasip1` and run by
wasmtime, also writes the reference case's record byte for byte as the native CLI
does.

Speed: the CI reference case's analysis (NACA 0012, N = 60, α = 2°) takes 11 ms in
wasm against 2.8 ms native; an N = 160 NACA 4412 at α = 4° 49 ms against 20 ms, and
its default 41-point polar 1.1 s against 0.5 s. The module is 824 kB (360 kB
gzipped, 267 kB brotli).

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
