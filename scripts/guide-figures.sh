#!/usr/bin/env bash
#
# Regenerate the Using yFoil pages' figures by running the commands the pages show, so each
# embedded figure is the output of the command printed beside it. Each page's figures go in the
# directory beside it named after it (docs/using-yfoil/your-first-analysis/, ...).
#
#   scripts/guide-figures.sh
#
# The work directory is .tmp/guide-figures/; only the .fig files are copied into the docs.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

cargo build --release --quiet --features plotting
yfoil="$repo_root/target/release/yfoil"

work=".tmp/guide-figures"
out="docs/using-yfoil"
rm -rf "$work"
mkdir -p "$work"
cd "$work"

# Your first analysis (your-first-analysis.md)
"$yfoil" geometry naca 0012 -n 160 -o naca0012.json
"$yfoil" analyse naca0012.json --alpha 5 --reynolds 1e9 -o naca0012_a5.json
"$yfoil" polar naca0012.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --label "NACA 0012" -o naca0012_polar.json
"$yfoil" plot polar naca0012_polar.json -o naca0012_polar.fig
"$yfoil" plot analysis naca0012_a5.json -o naca0012_a5.fig

# Sweeping a polar (commands/sweeping-a-polar.md)
"$yfoil" geometry naca 4412 -n 160 -o naca4412.json
"$yfoil" polar naca4412.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --max-iterations 20 \
    --label "NACA 4412" --distributions -o naca4412_polar.json
"$yfoil" plot polar naca4412_polar.json -o naca4412_polar.fig

# Plotting results (commands/plotting-results.md)
"$yfoil" plot polar naca0012_polar.json naca4412_polar.json \
    --title "NACA 0012 vs 4412" -o compare.fig
"$yfoil" analyse naca4412.json --alpha 5 --reynolds 1e6 -o a5.json
"$yfoil" plot analysis a5.json -o a5.fig
"$yfoil" plot foil naca4412.json --panels notches -o naca4412_foil.fig
"$yfoil" plot foil naca4412_polar.json --alpha 0,5,10 \
    --quantity dstar,theta --wake -o foil.fig

# The Using yFoil index's introductory comparison: the same two sections over a symmetric range,
# so the camber's effect on stall shows at both ends
for foil in 0012 4412; do
    "$yfoil" polar "naca$foil.json" \
        --alpha-min -15 --alpha-max 15 --alpha-step 0.5 \
        --reynolds 1e6 --label "NACA $foil" -o "naca${foil}_polar_symmetric.json"
done
"$yfoil" plot polar naca0012_polar_symmetric.json naca4412_polar_symmetric.json \
    --title "NACA 0012 vs 4412" -o intro.fig

# each page's figures, into the directory beside the page (index.md's in index/)
place() {
    local dir="$repo_root/$out/$1"
    shift
    mkdir -p "$dir"
    rm -f "$dir"/*.fig
    cp "$@" "$dir/"
}
place index intro.fig
place your-first-analysis naca0012_polar.fig naca0012_a5.fig
place commands/sweeping-a-polar naca4412_polar.fig
place commands/plotting-results compare.fig a5.fig naca4412_foil.fig foil.fig
echo "wrote 8 figures under $out"
