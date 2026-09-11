---
icon: lucide/download
---

# Install

## From source

yFoil is a Rust crate. With a [recent stable
toolchain](https://rustup.rs/) installed:

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
`yfoil plot`:

``` sh
cargo install --path . --features plotting
```

## Check it works

``` sh
yfoil --version
yfoil --help
```

!!! note "Binary releases"

    Prebuilt, cross-platform binaries are part of the plan — one of the reasons
    yFoil exists is that keeping XFOIL building on current toolchains and
    architectures has become painful. Until they are published, build from
    source.
