#!/usr/bin/env bash
#
# The gate for replacing cr-libm's C with a pure-Rust crate: compare every candidate with CORE-MATH
# bit for bit, natively, on WebAssembly and on x86-64 Linux.
#
#   scripts/libm-compare/run.sh native [MILLIONS]  # against the C; records every digest (default 16 M/set)
#   scripts/libm-compare/run.sh wasm               # wasm32-wasip1 under wasmtime, no C: against the
#                                                  # digests `native 2` recorded (C's and each candidate's own)
#   scripts/libm-compare/run.sh x86                # x86_64 Linux in Docker (gcc compiles the C), 2 M/set
#   scripts/libm-compare/run.sh time               # ns/call, natively and under wasmtime
#
# Results go to target/libm-compare/. Needs the worst-case lists (fetch-worst-cases.sh, run first
# if missing); wasm needs the wasm32-wasip1 target and wasmtime; x86 needs Docker and an amd64
# image with Rust and a C compiler ($X86_IMAGE, default rust:1-slim-bookworm, which has gcc).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$here/../.." && pwd)"
out="$repo_root/target/libm-compare"
wc="$out/worst-cases"
manifest="$here/Cargo.toml"
mode="${1:-native}"
mkdir -p "$out"
[ -f "$wc/COMMIT" ] || "$here/fetch-worst-cases.sh" "$wc"

native_bin() {
    cargo build -q --release --manifest-path "$manifest"
    echo "$here/target/release/libm-compare"
}

case "$mode" in
    native)
        millions="${2:-16}"
        bin="$(native_bin)"
        "$bin" compare --millions "$millions" --worst-cases "$wc" --digests-out "$out/digests-${millions}M.tsv" \
            >"$out/native-${millions}M.txt" 2>"$out/native-${millions}M.examples.txt"
        sed -n '/^summary/,$p' "$out/native-${millions}M.txt"
        ;;
    wasm)
        [ -f "$out/digests-2M.tsv" ] || "$0" native 2 >/dev/null
        cargo build -q --release --no-default-features --target wasm32-wasip1 --manifest-path "$manifest"
        wasm="$here/target/wasm32-wasip1/release/libm-compare.wasm"
        # no threads in wasip1: one process per candidate and function group
        pids=()
        for cand in metallic pxfm libm-0.2.16 libm-main; do
            for fns in exp,ln,log10,pow tanh,sin,cos,sin_cos atan,atan2,asin,hypot; do
                wasmtime run --dir "$out::/data" "$wasm" compare --millions 2 --worst-cases /data/worst-cases \
                    --digests-in /data/digests-2M.tsv --only "$fns" --candidates "$cand" \
                    >"$out/wasm-$cand-$fns.txt" 2>&1 &
                pids+=($!)
            done
        done
        for p in "${pids[@]}"; do wait "$p"; done
        cat "$out"/wasm-*.txt | grep -E '^  [a-z_0-9]+ +[a-z0-9.-]+ +[0-9]+ args' | sort
        ;;
    x86)
        [ -f "$out/digests-2M.tsv" ] || "$0" native 2 >/dev/null
        docker run --rm --platform linux/amd64 -v "$repo_root:/w" -w /w/scripts/libm-compare \
            "${X86_IMAGE:-rust:1-slim-bookworm}" sh -c '
                export CARGO_TARGET_DIR=/w/scripts/libm-compare/target/x86_64-docker
                cargo build -q --release
                $CARGO_TARGET_DIR/release/libm-compare compare --millions 2 --worst-cases /w/target/libm-compare/worst-cases \
                    --digests-in /w/target/libm-compare/digests-2M.tsv >/w/target/libm-compare/x86.txt 2>/w/target/libm-compare/x86.examples.txt'
        sed -n '/^summary/,$p' "$out/x86.txt"
        ;;
    time)
        "$(native_bin)" time | tee "$out/time-native.txt"
        cargo build -q --release --no-default-features --target wasm32-wasip1 --manifest-path "$manifest"
        wasmtime run "$here/target/wasm32-wasip1/release/libm-compare.wasm" time | tee "$out/time-wasm.txt"
        ;;
    *)
        echo "usage: scripts/libm-compare/run.sh native [MILLIONS] | wasm | x86 | time" >&2
        exit 1
        ;;
esac
