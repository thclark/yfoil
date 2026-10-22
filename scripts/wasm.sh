#!/usr/bin/env bash
#
# Build yFoil's WebAssembly binding (bindings/wasm), check it bit for bit against native, and
# package what was checked.
#
#   scripts/wasm.sh build            # target/wasm-pkg/{web,nodejs}/: the module and its JS glue
#   scripts/wasm.sh check [--quick]  # build, then run bindings/wasm/tests/identity.mjs
#   scripts/wasm.sh package [--tag TAG] [--out DIR]
#                                    # the built packages as release assets in DIR (default
#                                    # target/wasm-dist/): yfoil-wasm-<version>-{web,nodejs}.tar.gz
#                                    # and SHA256SUMS. Does not build: run check first.
#
# package names the assets after the crate version, or after TAG (v1.2.3 or 1.2.3), which must
# then equal the crate version; without a tag the version is <crate version>-dev.<commit>, so a
# CI artifact is never mistaken for a release.
#
# Needs: the wasm32-unknown-unknown Rust target; a clang with the WebAssembly backend (Apple's
# has none: `brew install llvm`, or set WASM_CC) and its llvm-ar; Node for the check; wasm-opt
# (binaryen) if present, or $WASM_OPT. wasm-bindgen-cli is installed into target/wasm-tools/ on
# first use, or taken from $WASM_BINDGEN. CORE-MATH is compiled with the native flags;
# cr-libm/build.rs adds only the include path of cr-libm/wasm-include/ (docs/bindings.md).
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
target_dir="${CARGO_TARGET_DIR:-$repo_root/target}"
pkg="target/wasm-pkg"
mode="${1:-build}"
shift || true

find_cc() {
    for cc in "${WASM_CC:-}" /opt/homebrew/opt/llvm/bin/clang /usr/local/opt/llvm/bin/clang clang; do
        [ -n "$cc" ] && command -v "$cc" >/dev/null 2>&1 || continue
        if "$cc" --print-targets 2>/dev/null | grep -q wasm32; then echo "$cc"; return; fi
    done
    echo "error: no clang with the wasm32 backend (brew install llvm, or set WASM_CC)" >&2
    exit 1
}

build() {
    local cc ar pinned bindgen wasm_opt module out wasm
    cc="$(find_cc)"
    ar="$("$cc" -print-prog-name=llvm-ar)"
    [ -x "$ar" ] || ar="$(command -v llvm-ar)" || { echo "error: no llvm-ar for $cc" >&2; exit 1; }

    # wasm-bindgen-cli must be the version the crate pins: $WASM_BINDGEN, else one installed into
    # target/wasm-tools/ (done here when missing, so another version elsewhere is left alone)
    pinned="$(sed -n 's/^wasm-bindgen = "=\(.*\)"/\1/p' bindings/wasm/Cargo.toml)"
    bindgen="${WASM_BINDGEN:-$target_dir/wasm-tools/bin/wasm-bindgen}"
    if [ "$("$bindgen" --version 2>/dev/null | awk '{print $2}')" != "$pinned" ]; then
        if [ -n "${WASM_BINDGEN:-}" ]; then
            echo "error: $WASM_BINDGEN is not wasm-bindgen-cli $pinned" >&2
            exit 1
        fi
        cargo install -q wasm-bindgen-cli --version "$pinned" --locked --root "$target_dir/wasm-tools"
    fi
    wasm_opt="${WASM_OPT:-$(command -v wasm-opt || true)}"

    export CC_wasm32_unknown_unknown="$cc" AR_wasm32_unknown_unknown="$ar"
    cargo build -q -p yfoil-wasm --target wasm32-unknown-unknown --release
    module="$target_dir/wasm32-unknown-unknown/release/yfoil_wasm.wasm"
    for target in web nodejs; do
        out="$pkg/$target"
        rm -rf "$out"
        "$bindgen" --target "$target" --out-dir "$out" "$module"
        if [ -n "$wasm_opt" ]; then
            # -O3 keeps the IEEE semantics (binaryen never reassociates or fuses floating point)
            "$wasm_opt" -O3 --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext \
                --enable-mutable-globals --enable-reference-types --enable-multivalue \
                "$out/yfoil_wasm_bg.wasm" -o "$out/yfoil_wasm_bg.wasm"
        fi
    done
    wasm="$pkg/web/yfoil_wasm_bg.wasm"
    echo "built $pkg/{web,nodejs}: $(wc -c <"$wasm" | tr -d ' ') bytes," \
        "$(gzip -9c "$wasm" | wc -c | tr -d ' ') gzipped; wasm-opt: ${wasm_opt:-not used}"
}

crate_version() { sed -n 's/^version = "\(.*\)"/\1/p' "$1" | head -1; }

package() {
    local tag="" out="target/wasm-dist" version yfoil wasm_crate name target dir tar_cmd
    while [ $# -gt 0 ]; do
        case "$1" in
            --tag) tag="$2"; shift 2 ;;
            --out) out="$2"; shift 2 ;;
            *) echo "error: package takes --tag TAG and --out DIR" >&2; exit 1 ;;
        esac
    done
    for target in web nodejs; do
        [ -f "$pkg/$target/yfoil_wasm_bg.wasm" ] || { echo "error: no $pkg/$target: run check first" >&2; exit 1; }
    done
    yfoil="$(crate_version Cargo.toml)"
    wasm_crate="$(crate_version bindings/wasm/Cargo.toml)"
    if [ "$yfoil" != "$wasm_crate" ]; then
        echo "error: yfoil is $yfoil but yfoil-wasm is $wasm_crate" >&2
        exit 1
    fi
    if [ -n "$tag" ]; then
        version="${tag#v}"
        if [ "$version" != "$yfoil" ]; then
            echo "error: tag $tag is version $version, but the crates are $yfoil" >&2
            exit 1
        fi
    else
        version="$yfoil-dev.$(git rev-parse --short=12 HEAD)"
    fi

    # reproducible archives: sorted entries, fixed times and owners, no gzip timestamp (GNU tar)
    tar_cmd="$(command -v gtar || command -v tar)"
    "$tar_cmd" --version 2>/dev/null | grep -q 'GNU tar' || { echo "error: GNU tar needed (gtar)" >&2; exit 1; }
    rm -rf "$out"
    mkdir -p "$out/stage"
    for target in web nodejs; do
        name="yfoil-wasm-$version-$target"
        dir="$out/stage/$name"
        mkdir -p "$dir"
        cp "$pkg/$target"/* "$dir/"
        # the module is GPL (yFoil) and contains CORE-MATH (MIT)
        cp LICENSE "$dir/LICENSE"
        cp third-party/core-math/LICENSE "$dir/LICENSE-core-math"
        "$tar_cmd" --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -C "$out/stage" -cf - "$name" \
            | gzip -9n >"$out/$name.tar.gz"
    done
    rm -rf "$out/stage"
    (cd "$out" && sha256sum yfoil-wasm-*.tar.gz 2>/dev/null >SHA256SUMS || shasum -a 256 yfoil-wasm-*.tar.gz >SHA256SUMS)
    echo "packaged $version into $out:"
    (cd "$out" && ls -l && cat SHA256SUMS)
}

case "$mode" in
    build) build ;;
    check)
        build
        cargo build -q --release --bin yfoil
        cargo build -q --release -p yfoil-wasm --example native
        node bindings/wasm/tests/identity.mjs --pkg "$pkg/nodejs" --yfoil "$target_dir/release/yfoil" \
            --native "$target_dir/release/examples/native" --work .tmp/wasm-identity "$@"
        ;;
    package) package "$@" ;;
    *) echo "usage: scripts/wasm.sh build | check [--quick] | package [--tag TAG] [--out DIR]" >&2; exit 1 ;;
esac
