# Test Fixtures

This directory contains validation fixtures generated from XFOIL for testing yFoil numerical accuracy.

## Directory Structure

```
fixtures/
├── README.md
├── naca0012/
│   ├── alpha_0_re_1e6/
│   │   ├── input.json
│   │   ├── geometry.json
│   │   ├── inviscid.json
│   │   ├── bl_stations.json
│   │   ├── viscal_iters.json
│   │   └── final.json
│   ├── alpha_2_re_1e6/
│   │   └── ...
│   └── alpha_5_re_1e6/
│       └── ...
└── naca4412/
    ├── alpha_0_re_1e6/
    │   └── ...
    └── alpha_4_re_1e6/
        └── ...
```

## Test Cases

| Airfoil | Alpha (°) | Re | Purpose |
|---------|-----------|-----|---------|
| NACA 0012 | 0 | 1e6 | Symmetric baseline |
| NACA 0012 | 2 | 1e6 | Small angle, attached |
| NACA 0012 | 5 | 1e6 | Moderate angle |
| NACA 4412 | 0 | 1e6 | Cambered, zero lift angle |
| NACA 4412 | 4 | 1e6 | Near design point |

## Fixture Files

### input.json

Test case parameters:
```json
{
  "airfoil": "NACA 0012",
  "alpha_deg": 2.0,
  "reynolds": 1000000.0,
  "mach": 0.0,
  "n_panels": 160,
  "n_crit": 9.0
}
```

### geometry.json

Panel geometry:
```json
{
  "n": 160,
  "x": [...],
  "y": [...],
  "s": [...],
  "nx": [...],
  "ny": [...]
}
```

### inviscid.json

Inviscid solution:
```json
{
  "alpha_rad": 0.0349066,
  "gamma": [...],
  "qinv": [...],
  "cpi": [...],
  "cl_inv": 0.2199,
  "cm_inv": -0.0244
}
```

### bl_stations.json

Boundary layer at each station:
```json
{
  "upper": {
    "n_stations": 84,
    "stations": [
      {
        "ibl": 1,
        "x": 0.0,
        "s": 0.0,
        "ue": 0.0,
        "delta_star": 0.0,
        "theta": 0.0,
        "hk": 2.59,
        "cf": 0.0,
        "ctau": 0.0,
        "regime": "stagnation"
      },
      ...
    ]
  },
  "lower": { ... }
}
```

### viscal_iters.json

Iteration history:
```json
{
  "n_iterations": 8,
  "converged": true,
  "iterations": [
    {
      "iter": 1,
      "alpha_deg": 2.0,
      "cl": 0.2195,
      "cd": 0.00712,
      "cdf": 0.00523,
      "cdp": 0.00189,
      "cm": -0.0245,
      "rmsbl": 0.0234,
      "rmxbl": 0.0891,
      "rlx": 1.0
    },
    ...
  ]
}
```

### final.json

Converged results:
```json
{
  "alpha_deg": 2.0,
  "cl": 0.2199,
  "cd": 0.00689,
  "cdf": 0.00512,
  "cdp": 0.00177,
  "cm": -0.0244,
  "xtr_upper": 0.412,
  "xtr_lower": 0.623
}
```

## Generating Fixtures

Use the `generate_fixtures.rs` example:

```bash
cargo run --example generate_fixtures
```

This runs the instrumented XFOIL binary and parses the output into JSON fixtures.

## Numerical Precision

All floating-point values are stored with full double precision (16 significant figures).

When comparing yFoil results:
- Relative tolerance: 1e-10
- Absolute tolerance: 1e-14 (for values near zero)

## Regenerating Fixtures

If XFOIL source is modified:

1. Rebuild XFOIL: `cd xfoil/xfoil6.99/bin && make clean && make`
2. Regenerate fixtures: `cargo run --example generate_fixtures`
3. Run validation: `cargo test --test fixtures`
