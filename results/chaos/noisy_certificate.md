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

On the exact draws: the first noise kernel's entry after the K-th, τ_{K+1} − τ_K, against the first-order W — the ratio's median, range and 5–95% quantiles, and their half-width over s (constant if the relative error is O(s)) —; the kernel as predicted, and on the draws where it is not, the relative gap between the two smallest predicted windows (W₍₂₎ − W₍₁₎)/W₍₁₎ (median; all draws' median in brackets); W·s (median, quartiles).

| s | exact draws | W_obs/W_pred: median [min, max] | 5–95% | half-width / s | next kernel as predicted | missed: (W₍₂₎ − W₍₁₎)/W₍₁₎ [all] | W·s: median [quartiles] |
|---|---|---|---|---|---|---|---|
| 0.003 | 400 | 1.0000 [0.975, 1.043] | [0.9934, 1.0112] | 2.97 | 397/400 | 0.0115 [0.124] | 13.6 [10.1, 18.7] |
| 0.01 | 400 | 1.0003 [0.918, 1.193] | [0.9799, 1.0337] | 2.69 | 386/400 | 0.0177 [0.126] | 13.2 [9.5, 18.5] |
| 0.03 | 396 | 1.0007 [0.888, 1.550] | [0.9537, 1.1234] | 2.83 | 364/396 | 0.0187 [0.128] | 13.8 [10.1, 19.9] |
| 0.1 | 374 | 1.0080 [-2.314, 2.815] | [0.8859, 1.4738] | 2.94 | 313/374 | 0.117 [0.181] | 12.3 [8.08, 17.8] |
| 0.3 | 191 | 1.0004 [-16.305, 15.250] | [-0.3378, 3.1301] | 5.78 | 138/191 | 0.348 [0.306] | 8.94 [5.35, 16.1] |

## D. Stopping and the oracle

The discrepancy rule: the first entry after which the relaxed residual's rms is ≤ 1.1s, on the exact draws: at the K-th, not by the K-th (late), earlier; its first-order prediction — at the K-th iff rms(R_Kξ) ≤ 1.1s, late otherwise, R_Kξ the first-order residual (|R_Kξ|²/s² ~ χ²_(N−2K)) — and its agreement per draw; the first-order probability of a late stop, P(χ²_(N−2K) > 1.1²N) = 0.0461, times the draws (on all draws, exact or not, the count with rms(R_Kξ) > (1 + δ)s in brackets); the relaxed residual after the K-th entry against rms(R_Kξ).

| s | exact draws | at K | late | early | first order: at K | agree per draw | late: expected [all draws] | rms_K / rms(R_Kξ): median [5–95%] |
|---|---|---|---|---|---|---|---|---|
| 0.003 | 400 | 383 | 17 | 0 | 383 | 400 | 18.4 [17/400] | 1.0000 [0.9999, 1.0001] |
| 0.01 | 400 | 389 | 11 | 0 | 389 | 400 | 18.4 [11/400] | 1.0000 [0.9996, 1.0004] |
| 0.03 | 396 | 379 | 17 | 0 | 378 | 395 | 18.2 [18/400] | 1.0000 [0.9991, 1.0010] |
| 0.1 | 374 | 367 | 7 | 0 | 367 | 374 | 17.2 [9/400] | 1.0000 [0.9970, 1.0030] |
| 0.3 | 191 | 107 | 9 | 75 | 182 | 116 | 8.8 [18/400] | 1.0002 [0.9934, 1.0071] |

The largest relative gap (τ_{k+1} − τ_k)/τ_k over the first K + 1 entries falling after the K-th; the fit stopped at the K-th entry: its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against s(2K/N)^{1/2}.

| s | exact draws | largest gap at K | error / rms(P_Tξ): median [range] | 5–95% | error / (s(2K/N)^{1/2}): median |
|---|---|---|---|---|---|
| 0.003 | 400 | 400/400 | 1.0000 [0.995, 1.005] | [0.9982, 1.0018] | 0.973 |
| 0.01 | 400 | 400/400 | 1.0001 [0.980, 1.016] | [0.9937, 1.0066] | 0.925 |
| 0.03 | 396 | 396/396 | 1.0003 [0.945, 1.044] | [0.9861, 1.0188] | 0.952 |
| 0.1 | 374 | 374/374 | 0.9986 [0.863, 1.188] | [0.9443, 1.0562] | 0.928 |
| 0.3 | 191 | 161/191 | 0.9957 [0.764, 1.754] | [0.8754, 1.2388] | 0.831 |

## E. Past the first noise kernel

The first 20 targets, 5 draws per level, the path to its (K + 4)-th entry (each path under 180 s): on the draws whose first K entries are E, the largest relative gap over the K + 4 entries falling after the K-th, the discrepancy rule over the relaxed fits after entries 1 … K + 3 (at the K-th, later, never), and the fits the two rules stop at against the oracle's (the K-th) error. Paths with a chunk of the fast flow redone by Radau (`relax_chunk`, LSODA's stall): 79 of 399 finished.

| s | paths (over the limit) | exact | largest gap at K | discrepancy: at K, later, never | discrepancy-stopped error / oracle: median [max] | gap-stopped error / oracle: median [max] |
|---|---|---|---|---|---|---|
| 0.01 | 100 (0) | 100 | 100/100 | 94, 3, 3 | 1.000 [1.77] | 1.000 [1] |
| 0.03 | 100 (0) | 100 | 100/100 | 99, 1, 0 | 1.000 [2.57] | 1.000 [1] |
| 0.1 | 100 (0) | 90 | 87/90 | 88, 1, 1 | 1.000 [1.82] | 1.000 [1.39] |
| 0.3 | 100 (1) | 52 | 34/52 | 29, 1, 2 | 1.000 [8.3] | 1.000 [3.13] |

Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.
