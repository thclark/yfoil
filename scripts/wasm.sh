#!/usr/bin/env bash
#
# Build yFoil's WebAssembly binding (bindings/wasm) and check it bit for bit against native.
#
#   scripts/wasm.sh build            # target/wasm-pkg/{web,nodejs}/: the module and its JS glue
#   scripts/wasm.sh check [--quick]  # build, then run bindings/wasm/tests/identity.mjs
#
# Needs: the wasm32-unknown-unknown Rust target; a clang with the WebAssembly backend (Apple's
# has none: `brew install llvm`, or set WASM_CC) and its llvm-ar; Node for the check; wasm-opt
# (binaryen) if present. wasm-bindgen-cli is installed into target/wasm-tools/ on first use.
# CORE-MATH is compiled with the native flags; cr-libm/build.rs adds only the include path of
# cr-libm/wasm-include/ (docs/bindings.md).
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
target_dir="${CARGO_TARGET_DIR:-$repo_root/target}"
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

export CC_wasm32_unknown_unknown="$cc" AR_wasm32_unknown_unknown="$ar"
cargo build -q -p yfoil-wasm --target wasm32-unknown-unknown --release
module="$target_dir/wasm32-unknown-unknown/release/yfoil_wasm.wasm"
for target in web nodejs; do
    out="target/wasm-pkg/$target"
    rm -rf "$out"
    "$bindgen" --target "$target" --out-dir "$out" "$module"
    if command -v wasm-opt >/dev/null 2>&1; then
        # -O3 keeps the IEEE semantics (binaryen never reassociates or fuses floating point)
        wasm-opt -O3 --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext \
            --enable-mutable-globals --enable-reference-types --enable-multivalue \
            "$out/yfoil_wasm_bg.wasm" -o "$out/yfoil_wasm_bg.wasm"
    fi
done
wasm="target/wasm-pkg/web/yfoil_wasm_bg.wasm"
echo "built target/wasm-pkg/{web,nodejs}: $(wc -c <"$wasm" | tr -d ' ') bytes, $(gzip -9c "$wasm" | wc -c | tr -d ' ') gzipped"

case "$mode" in
    build) ;;
    check)
        cargo build -q --release --bin yfoil
        cargo build -q --release -p yfoil-wasm --example native
        node bindings/wasm/tests/identity.mjs --pkg target/wasm-pkg/nodejs --yfoil "$target_dir/release/yfoil" \
            --native "$target_dir/release/examples/native" --work .tmp/wasm-identity "$@"
        ;;
    *) echo "usage: scripts/wasm.sh build | check [--quick]" >&2; exit 1 ;;
esac
