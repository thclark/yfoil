# 5.2 LIMAGE is permanently false — Out of scope

The ground-image toggle is commented out (`xoper.f:405-410`), so PSILIN's image-airfoil block
(`xpanel.f:488` and 12 dependent lines) is dead. Not translated (`src/solver/psilin.rs:6-7`).
GEOLIN (geometric sensitivities) is inverse-design only and also not translated.
