# Fair 1-D kernel comparison (per-step prediction for every method, covering inducing grid, seed-level statistics)

Seeds: 3 (different reference / plant initial conditions). Steady state = last 40% of 30 s. Step time = median full online step (prediction + amortised training) after warm-up, same process type, 2 BLAS threads.

| method | mean SS ‖e‖ (mean ± 95% CI over seeds) | max SS ‖e‖ (mean) | median step (ms) | mean step (ms) | seeds where proposed has lower mean SS |
|---|---|---|---|---|---|
| Proposed (Proj) | 0.1978 ± 0.0503 | 1.147 | 0.060 | 0.069 | — |
| Exact GP | 0.0591 ± 0.0061 | 0.780 | 6.309 | 13.936 | 0/3 |
| Sparse GP (9×3 grid) | 0.1867 ± 0.0210 | 0.699 | 1.353 | 2.732 | 2/3 |
| Deep GP | 0.3087 ± 0.1121 | 1.262 | 4.461 | 5.495 | 2/3 |
| NTK | 0.4700 ± 0.0250 | 1.438 | 0.212 | 0.503 | 3/3 |
