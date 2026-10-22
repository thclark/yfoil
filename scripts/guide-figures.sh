#!/usr/bin/env bash
#
# Regenerate the user guide's figures (docs/user-guide/figures/) by running the commands the
# guide shows, so each embedded figure is the output of the command printed beside it.
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
out="docs/user-guide/figures"
rm -rf "$work"
mkdir -p "$work" "$out"
cd "$work"

# Your first analysis (get-started/your-first-analysis.md)
"$yfoil" geometry naca 0012 -n 160 -o naca0012.json
"$yfoil" analyse naca0012.json --alpha 5 --reynolds 1e9 -o naca0012_a5.json
"$yfoil" polar naca0012.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --label "NACA 0012" -o naca0012_polar.json
"$yfoil" plot polar naca0012_polar.json -o naca0012_polar.fig
"$yfoil" plot analysis naca0012_a5.json -o naca0012_a5.fig

# Polars (commands/polars.md)
"$yfoil" geometry naca 4412 -n 160 -o naca4412.json
"$yfoil" polar naca4412.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --max-iterations 20 \
    --label "NACA 4412" --distributions -o naca4412_polar.json
"$yfoil" plot polar naca4412_polar.json -o naca4412_polar.fig

# Plotting (commands/plotting.md); the guide index shows the comparison
"$yfoil" plot polar naca0012_polar.json naca4412_polar.json \
    --title "NACA 0012 vs 4412" -o compare.fig
"$yfoil" analyse naca4412.json --alpha 5 --reynolds 1e6 -o a5.json
"$yfoil" plot analysis a5.json -o a5.fig
"$yfoil" plot foil naca4412.json --panels notches -o naca4412_foil.fig
"$yfoil" plot foil naca4412_polar.json --alpha 0,5,10 \
    --quantity dstar,theta --wake -o foil.fig

# The guide index's introductory comparison: the same two sections over a symmetric range, so
# the camber's effect on stall shows at both ends
for foil in 0012 4412; do
    "$yfoil" polar "naca$foil.json" \
        --alpha-min -15 --alpha-max 15 --alpha-step 0.5 \
        --reynolds 1e6 --label "NACA $foil" -o "naca${foil}_polar_symmetric.json"
done
"$yfoil" plot polar naca0012_polar_symmetric.json naca4412_polar_symmetric.json \
    --title "NACA 0012 vs 4412" -o intro.fig

rm -f "$repo_root/$out"/*.fig
cp intro.fig naca0012_polar.fig naca0012_a5.fig naca4412_polar.fig compare.fig a5.fig naca4412_foil.fig foil.fig \
    "$repo_root/$out/"
echo "wrote $(ls "$repo_root/$out"/*.fig | wc -l | tr -d ' ') figures to $out"
