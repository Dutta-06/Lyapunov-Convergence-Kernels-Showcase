# Proposition D‴.56: the transit under measurement noise (`exp_noisy_certificate.py`)

The sliding path on y = f* + ξ (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]) against its first-order response along the noiseless path (`path_response`); targets: three features at least 0.45 apart (D‴.55's recipe), certificate margin m > 0.8q.

## A. The response is the first order

Per target (one noise direction ξ, unit variance per sample) and level h: the noisy path's entry times and its duals at the K-th entry against the first-order prediction — the largest error, and the first-order term's own size.

| target | h | same entries | entry times: error (first-order term) | duals: error (first-order term) |
|---|---|---|---|---|
| 0 | 0.001 | True | 1.26e-05 (1.24e-02) | 4.81e-07 (6.57e-04) |
| 0 | 0.0001 | True | 1.26e-07 (1.24e-03) | 4.82e-09 (6.57e-05) |
| 0 | 1e-05 | True | 1.26e-09 (1.24e-04) | 4.82e-11 (6.57e-06) |
| 1 | 0.001 | True | 2.06e-05 (1.84e-02) | 5.03e-07 (1.13e-03) |
| 1 | 0.0001 | True | 2.06e-07 (1.84e-03) | 5.02e-09 (1.13e-04) |
| 1 | 1e-05 | True | 2.06e-09 (1.84e-04) | 5.03e-11 (1.13e-05) |
| 2 | 0.001 | True | 1.43e-05 (1.53e-02) | 1.26e-06 (1.33e-03) |
| 2 | 0.0001 | True | 1.43e-07 (1.53e-03) | 1.26e-08 (1.33e-04) |
| 2 | 1e-05 | True | 1.43e-09 (1.53e-04) | 1.26e-10 (1.33e-05) |
| 3 | 0.001 | True | 3.06e-06 (7.89e-03) | 8.54e-07 (1.78e-03) |
| 3 | 0.0001 | True | 3.07e-08 (7.89e-04) | 8.54e-09 (1.78e-04) |
| 3 | 1e-05 | True | 3.07e-10 (7.89e-05) | 8.54e-11 (1.78e-05) |
| 4 | 0.001 | True | 7.64e-06 (1.06e-02) | 3.31e-07 (4.79e-04) |
| 4 | 0.0001 | True | 7.64e-08 (1.06e-03) | 3.30e-09 (4.79e-05) |
| 4 | 1e-05 | True | 7.64e-10 (1.06e-04) | 3.30e-11 (4.79e-06) |
| 5 | 0.001 | True | 3.55e-06 (6.23e-03) | 1.42e-07 (3.89e-04) |
| 5 | 0.0001 | True | 3.56e-08 (6.23e-04) | 1.42e-09 (3.89e-05) |
| 5 | 1e-05 | True | 3.56e-10 (6.23e-05) | 1.42e-11 (3.89e-06) |
| 6 | 0.001 | True | 3.87e-06 (6.60e-03) | 9.74e-08 (1.35e-03) |
| 6 | 0.0001 | True | 3.87e-08 (6.60e-04) | 9.74e-10 (1.35e-04) |
| 6 | 1e-05 | True | 3.87e-10 (6.60e-05) | 9.74e-12 (1.35e-05) |
| 7 | 0.001 | True | 3.16e-05 (2.55e-02) | 1.83e-06 (1.66e-03) |
| 7 | 0.0001 | True | 3.16e-07 (2.55e-03) | 1.83e-08 (1.66e-04) |
| 7 | 1e-05 | True | 3.16e-09 (2.55e-04) | 1.83e-10 (1.66e-05) |
| 8 | 0.001 | True | 2.84e-05 (2.08e-02) | 1.74e-06 (9.28e-04) |
| 8 | 0.0001 | True | 2.84e-07 (2.08e-03) | 1.74e-08 (9.28e-05) |
| 8 | 1e-05 | True | 2.84e-09 (2.08e-04) | 1.74e-10 (9.28e-06) |
| 9 | 0.001 | True | 1.42e-05 (1.33e-02) | 6.08e-07 (1.55e-03) |
| 9 | 0.0001 | True | 1.43e-07 (1.33e-03) | 6.08e-09 (1.55e-04) |
| 9 | 1e-05 | True | 1.43e-09 (1.33e-04) | 6.08e-11 (1.55e-05) |

The duals' error falls by 99.8–100.1 per decade of h (second order: 100).

## B. The support

Per noise level s: 10 draws on each of the 40 targets (40 with the noiseless path exactly E). Observed: the path enters exactly E first; predicted per draw: the first-order criterion max|u⁰_j + U_jξ| < 1 over the entries and the kernels off E; the expected count from the linear model's probability (20000 Monte Carlo draws per target).

| s | exact (observed) | criterion < 1 | agree per draw | expected (linear model) | observed exact, criterion ≥ 1 | not exact, criterion < 1 |
|---|---|---|---|---|---|---|
| 0.003 | 400/400 | 400 | 400 | 400.0 | 0 | 0 |
| 0.01 | 400/400 | 400 | 400 | 400.0 | 0 | 0 |
| 0.03 | 396/400 | 395 | 399 | 398.0 | 1 | 0 |
| 0.1 | 374/400 | 370 | 388 | 359.9 | 8 | 4 |
| 0.3 | 191/400 | 155 | 338 | 159.9 | 49 | 13 |

## C. The window

On the exact draws: the first noise kernel's entry after the K-th, τ_{K+1} − τ_K, against the first-order W (the ratio's median and range), the kernel as predicted, and W·s (median, quartiles).

| s | exact draws | W_obs/W_pred: median [min, max] | next kernel as predicted | W·s: median [quartiles] |
|---|---|---|---|---|
| 0.003 | 400 | 1.0000 [0.975, 1.043] | 397/400 | 13.6 [10.1, 18.7] |
| 0.01 | 400 | 1.0003 [0.918, 1.193] | 386/400 | 13.2 [9.5, 18.5] |
| 0.03 | 396 | 1.0007 [0.888, 1.550] | 364/396 | 13.8 [10.1, 19.9] |
| 0.1 | 374 | 1.0080 [-2.314, 2.815] | 313/374 | 12.3 [8.08, 17.8] |
| 0.3 | 191 | 1.0004 [-16.305, 15.250] | 138/191 | 8.94 [5.35, 16.1] |

## D. Stopping and the oracle

On the exact draws: the discrepancy rule (the first entry after which the relaxed residual's rms ≤ 1.1s) stopping at the K-th; the largest relative gap (τ_{k+1} − τ_k)/τ_k over the K + 1 entries falling after the K-th; the fit stopped at the K-th entry: its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against s(2K/N)^{1/2}.

| s | exact draws | discrepancy rule at K | largest gap at K | error / rms(P_Tξ): median [range] | error / (s(2K/N)^{1/2}): median |
|---|---|---|---|---|---|
| 0.003 | 400 | 383/400 | 400/400 | 1.0000 [0.995, 1.005] | 0.973 |
| 0.01 | 400 | 389/400 | 400/400 | 1.0001 [0.980, 1.016] | 0.925 |
| 0.03 | 396 | 379/396 | 396/396 | 1.0003 [0.945, 1.044] | 0.952 |
| 0.1 | 374 | 367/374 | 374/374 | 0.9986 [0.863, 1.188] | 0.928 |
| 0.3 | 191 | 107/191 | 161/191 | 0.9957 [0.764, 1.754] | 0.831 |

Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.
