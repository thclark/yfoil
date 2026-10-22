---
icon: lucide/download
---

# Installing yFoil

yFoil comes two ways: as the `yfoil` command-line tool, which is what the rest of this tab
uses, and as a Rust library. There are no binary releases and no bindings for other
languages yet ([calling yFoil from code](commands/calling-yfoil-from-code.md) says what
exists).

## The command-line tool

yFoil is a Rust crate. With a [recent stable
toolchain](https://rustup.rs/) and a C compiler installed (Xcode's command-line
tools on macOS, `build-essential` or equivalent on Linux; it compiles the
correctly rounded maths functions that make yFoil's results the same on every
host):

``` sh
git clone https://github.com/thclark/yfoil.git
cd yfoil
cargo build --release
```

The binary is at `target/release/yfoil`. Put it on your `PATH`, or install it
into `~/.cargo/bin`:

``` sh
cargo install --path .
```

Plotting is behind a cargo feature; build with it enabled if you want
`yfoil plot`, and install the ironlab viewer to open the figures it writes:

``` sh
cargo install --path . --features plotting
cargo install ironlab-viewer --version 0.12.0
```

Check it works:

``` sh
yfoil --version
yfoil --help
```

!!! note "Binary releases"

    Prebuilt, cross-platform binaries are part of the plan — one of the reasons
    yFoil exists is that keeping XFOIL building on current toolchains and
    architectures has become painful. Until they are published, build from
    source.

## As a Rust library

The same crate is a library. Add it to a Rust project from the repository (it is not on
crates.io yet):

``` sh
cargo add --git https://github.com/thclark/yfoil.git yfoil
```

The same C compiler is needed, for the same reason. What the library offers is on
[calling yFoil from code](commands/calling-yfoil-from-code.md).

Next: [your first analysis](your-first-analysis.md).
