# Proposition D‴.58: the transit under width mismatch (`exp_width_mismatch.py`)

Kernels of width σ = 0.12; a feature Aφ_σ*(· − c*), σ* = σ(1 + (τ/σ)²)^{1/2}; t = (σ* − σ)/(σ* + σ), λ = τ²/(σ² + σ*²) (`src/core/wide_feature.py`).

## A. The rules on ℝ

Relative L² errors: the half-squeeze rule (nodes τ(σ/σ*)^{1/2}z_k, amplitudes A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*)))) by the trapezoid rule and by the Mehler series (their relative difference), the series' first term's share of the squared error; moment matching (nodes τz_k, amplitudes A(σ*/σ)w_k); the half-squeeze nodes with least-squares amplitudes; the optimal K-kernel fit (variable projection; not computed below 1e-10, where its solve stops short of the optimum — A3 marks the optimum unconverged where the one-parameter rule beats it). Rows stop where the series falls below 1e-13 (the trapezoid rule's floor). Ratios per kernel e_{K−1}/e_K; the optimum's nodes over the half-squeeze nodes (range over the nodes).

### τ/σ = 0.25 (1/t = 65.9848, 1/λ = 33.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.1431e-02 | 2.1431e-02 (1e-15) | 1.00000 | 2.645e-02 | 2.1430e-02 | 2.1430e-02 | — | — |
| 2 | 3.7504e-04 | 3.7504e-04 (5e-15) | 0.99999 | 7.939e-04 | 3.7504e-04 | 3.7502e-04 | 57.1443, 57.1432, 33.3161 | 0.9999 … 0.9999 |
| 3 | 6.2262e-06 | 6.2262e-06 (1e-12) | 0.99998 | 2.399e-05 | 6.2262e-06 | 6.2260e-06 | 60.2354, 60.2352, 33.0854 | 1.0000 … 1.0000 |
| 4 | 1.0087e-07 | 1.0087e-07 (2e-10) | 0.99998 | 7.263e-07 | 1.0087e-07 | 1.0087e-07 | 61.7230, 61.7229, 33.0349 | 1.0000 … 1.0000 |
| 5 | 1.6114e-09 | 1.6114e-09 (4e-09) | 0.99997 | 2.200e-08 | 1.6114e-09 | 1.6114e-09 | 62.5986, 62.5985, 33.0177 | 1.0000 … 1.0000 |
| 6 | 2.5507e-11 | 2.5507e-11 (4e-07) | 0.99997 | 6.664e-10 | 2.5507e-11 | — | 63.1756, —, 33.0102 | — |
| 7 | 4.0116e-13 | 4.0115e-13 (3e-05) | 0.99997 | 2.019e-11 | 4.0116e-13 | — | 63.5825, —, 33.0064 | — |

### τ/σ = 0.5 (1/t = 17.9443, 1/λ = 9.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 7.8750e-02 | 7.8750e-02 (1e-15) | 0.99999 | 9.924e-02 | 7.8689e-02 | 7.8689e-02 | — | — |
| 2 | 5.0678e-03 | 5.0678e-03 (7e-16) | 0.99989 | 1.105e-02 | 5.0678e-03 | 5.0651e-03 | 15.5394, 15.5354, 8.9774 | 0.9990 … 0.9990 |
| 3 | 3.0939e-04 | 3.0939e-04 (9e-16) | 0.99978 | 1.232e-03 | 3.0939e-04 | 3.0925e-04 | 16.3798, 16.3790, 8.9698 | 0.9994 … 0.9994 |
| 4 | 1.8433e-05 | 1.8433e-05 (5e-13) | 0.99969 | 1.373e-04 | 1.8433e-05 | 1.8425e-05 | 16.7846, 16.7842, 8.9786 | 0.9996 … 0.9996 |
| 5 | 1.0828e-06 | 1.0828e-06 (1e-11) | 0.99962 | 1.528e-05 | 1.0828e-06 | 1.0824e-06 | 17.0229, 17.0227, 8.9849 | 0.9997 … 0.9997 |
| 6 | 6.3030e-08 | 6.3030e-08 (3e-10) | 0.99957 | 1.700e-06 | 6.3030e-08 | 6.3003e-08 | 17.1799, 17.1798, 8.9889 | 0.9997 … 0.9997 |
| 7 | 3.6452e-09 | 3.6452e-09 (4e-09) | 0.99953 | 1.890e-07 | 3.6452e-09 | 3.6437e-09 | 17.2912, 17.2911, 8.9916 | 0.9998 … 0.9998 |
| 8 | 2.0980e-10 | 2.0980e-10 (8e-08) | 0.99950 | 2.102e-08 | 2.0980e-10 | 2.0972e-10 | 17.3742, 17.3738, 8.9934 | 0.9998 … 0.9998 |
| 9 | 1.2031e-11 | 1.2031e-11 (1e-06) | 0.99948 | 2.337e-09 | 1.2031e-11 | — | 17.4385, —, 8.9947 | — |
| 10 | 6.8790e-13 | 6.8790e-13 (1e-05) | 0.99946 | 2.597e-10 | 6.8792e-13 | — | 17.4895, —, 8.9957 | — |

### τ/σ = 1 (1/t = 5.8284, 1/λ = 3.0000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 2.4090e-01 | 2.4090e-01 (1e-16) | 0.99946 | 3.237e-01 | 2.3915e-01 | 2.3915e-01 | — | — |
| 2 | 4.7745e-02 | 4.7745e-02 (3e-16) | 0.99867 | 1.112e-01 | 4.7723e-02 | 4.7497e-02 | 5.0456, 5.0350, 2.9107 | 0.9901 … 0.9901 |
| 3 | 8.9780e-03 | 8.9780e-03 (1e-15) | 0.99768 | 3.757e-02 | 8.9771e-03 | 8.9365e-03 | 5.3180, 5.3149, 2.9608 | 0.9939 … 0.9939 |
| 4 | 1.6474e-03 | 1.6474e-03 (4e-15) | 0.99691 | 1.261e-02 | 1.6474e-03 | 1.6402e-03 | 5.4499, 5.4485, 2.9793 | 0.9956 … 0.9957 |
| 5 | 2.9802e-04 | 2.9802e-04 (5e-14) | 0.99634 | 4.221e-03 | 2.9802e-04 | 2.9676e-04 | 5.5277, 5.5269, 2.9873 | 0.9966 … 0.9967 |
| 6 | 5.3418e-05 | 5.3418e-05 (2e-13) | 0.99590 | 1.411e-03 | 5.3418e-05 | 5.3197e-05 | 5.5791, 5.5785, 2.9915 | 0.9972 … 0.9973 |
| 7 | 9.5126e-06 | 9.5126e-06 (1e-12) | 0.99555 | 4.713e-04 | 9.5126e-06 | 9.4739e-06 | 5.6154, 5.6151, 2.9939 | 0.9976 … 0.9977 |
| 8 | 1.6859e-06 | 1.6859e-06 (6e-12) | 0.99527 | 1.573e-04 | 1.6859e-06 | 1.6791e-06 | 5.6426, 5.6423, 2.9954 | 0.9979 … 0.9980 |
| 9 | 2.9767e-07 | 2.9767e-07 (5e-11) | 0.99504 | 5.251e-05 | 2.9767e-07 | 2.9649e-07 | 5.6636, 5.6633, 2.9964 | 0.9982 … 0.9983 |
| 10 | 5.2404e-08 | 5.2404e-08 (1e-10) | 0.99484 | 1.752e-05 | 5.2404e-08 | 5.2197e-08 | 5.6803, 5.6801, 2.9971 | 0.9984 … 0.9984 |

### τ/σ = 2 (1/t = 2.6180, 1/λ = 1.5000)

| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |
|---|---|---|---|---|---|---|---|---|
| 1 | 5.2280e-01 | 5.2280e-01 (0e+00) | 0.98665 | 8.088e-01 | 5.0462e-01 | 5.0462e-01 | — | — |
| 2 | 2.3059e-01 | 2.3059e-01 (2e-16) | 0.98657 | 5.371e-01 | 2.2815e-01 | 2.2364e-01 | 2.2672, 2.2564, 1.5059 | 0.9485 … 0.9485 |
| 3 | 9.6611e-02 | 9.6611e-02 (4e-16) | 0.98401 | 3.600e-01 | 9.6247e-02 | 9.4014e-02 | 2.3868, 2.3787, 1.4920 | 0.9657 … 0.9657 |
| 4 | 3.9497e-02 | 3.9497e-02 (7e-16) | 0.98168 | 2.408e-01 | 3.9438e-02 | 3.8509e-02 | 2.4460, 2.4414, 1.4946 | 0.9746 … 0.9770 |
| 5 | 1.5918e-02 | 1.5918e-02 (9e-16) | 0.97978 | 1.609e-01 | 1.5907e-02 | 1.5539e-02 | 2.4813, 2.4783, 1.4966 | 0.9799 … 0.9818 |
| 6 | 6.3555e-03 | 6.3555e-03 (2e-16) | 0.97824 | 1.074e-01 | 6.3532e-03 | 6.2092e-03 | 2.5046, 2.5025, 1.4977 | 0.9834 … 0.9854 |
| 7 | 2.5208e-03 | 2.5208e-03 (4e-15) | 0.97697 | 7.171e-02 | 2.5203e-03 | 2.4644e-03 | 2.5212, 2.5196, 1.4983 | 0.9858 … 0.9876 |
| 8 | 9.9499e-04 | 9.9499e-04 (5e-15) | 0.97591 | 4.785e-02 | 9.9488e-04 | 9.7318e-04 | 2.5335, 2.5323, 1.4987 | 0.9877 … 0.9894 |
| 9 | 3.9125e-04 | 3.9125e-04 (8e-16) | 0.97502 | 3.192e-02 | 3.9122e-04 | 3.8282e-04 | 2.5431, 2.5421, 1.4990 | 0.9891 … 0.9906 |
| 10 | 1.5339e-04 | 1.5339e-04 (1e-14) | 0.97426 | 2.129e-02 | 1.5338e-04 | 1.5013e-04 | 2.5507, 2.5499, 1.4992 | 0.9902 … 0.9917 |

### A2. The series' tail

D_{2K+2j,K}/D_{2K,K} (integers) for j = 1, …, 4 and K = 1, …, 30; their j-th differences in K, constant if they are polynomials of degree j, against j!·binom(1/2, j)(−2)^j = −1, −1, −3, −15.

| j | K = 1, 2, 3, 4, 5, …, 30 | j-th differences (all K) | j!·binom(1/2, j)(−2)^j |
|---|---|---|---|
| 1 | 0, -1, -2, -3, -4, …, -29 | -1 | -1 |
| 2 | 15, 26, 36, 45, 53, …, -72 | -1 | -1 |
| 3 | 0, 86, 270, 549, 920, …, 31320 | -3 | -3 |
| 4 | 945, 751, 513, 1206, 3790, …, 3183840 | -15 | -15 |

### A3. Scaled rules

The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`) at the half-squeeze variance v_hs = τ²σ/σ*, the natural Fock scaling v_hs/cosh r and the Mercer-zero scaling τ²(1 + 2τ²/σ²)^{−1/2}: their ratios per kernel at the last K against 1/t, (2 − λ)/λ = (1 − t + t²)/t and 1/q, q = (β² − 1)/(β² + 1), β² = (1 + 2τ²/σ²)^{1/2}; the variance of least error (`best_scale`): v/v_hs, its error and the optimal fit's over the rule's, at K = 1 and at the last K with an optimum (in brackets the share of the gap to the optimum that the variance recovers, (e_hs − e_best)/(e_hs − min(e_opt, e_best)); 'optimum unconverged' where e_best < e_opt by more than a relative 1e-09); over all K with an optimum the shares 93.0%–100.0% (1 optima unconverged):

| τ/σ | half-squeeze (1/t) | Fock ((2 − λ)/λ) | Mercer-zero (1/q) | K = 1: v/v_hs, e_best/e_hs, e_opt/e_hs | the last K: K, v/v_hs, e_best/e_hs, e_opt/e_hs |
|---|---|---|---|---|---|
| 0.25 | 63.583 (65.985) | 63.500 (65.000) | 34.018 (33.971) | 0.99954, 0.99994, 0.99994 (100.0%) | 5, 0.99995, 0.99997, 0.99997 (100.0%) |
| 0.5 | 17.489 (17.944) | 17.079 (17.000) | 9.916 (9.899) | 0.99381, 0.99922, 0.99922 (100.0%) | 8, 0.99959, 0.99959, 0.99960 (100.0%, optimum unconverged) |
| 1 | 5.680 (5.828) | 5.034 (5.000) | 3.751 (3.732) | 0.94281, 0.99272, 0.99272 (100.0%) | 10, 0.99678, 0.99608, 0.99605 (99.1%) |
| 2 | 2.551 (2.618) | 2.021 (2.000) | 2.021 (2.000) | 0.74536, 0.96523, 0.96523 (100.0%) | 10, 0.98131, 0.97968, 0.97877 (95.7%) |

### A4. Narrower features: the continued frame

σ* < σ: the nodes imaginary in conjugate pairs, the K-term sum e^{−x²/(2σ²)}Σ_k b_k cos(η_kx/σ²) at the Hermite-zero frequencies η_k = (2σ²|sinh r|)^{1/2}z_k (`narrow_rule`); its relative L² error (trapezoid) against the Mehler series at |t| = (σ − σ*)/(σ + σ*), and the ratios per kernel against 1/|t|:

| σ*/σ | 1/|t| | K: error (relative difference from the series) | ratio per kernel at the last K |
|---|---|---|---|
| 0.5 | 3.0000 | 1: 4.595e-01 (2e-16), 2: 1.769e-01 (0e+00), 3: 6.469e-02 (7e-16), 4: 2.308e-02 (0e+00), 5: 8.115e-03 (2e-16), 6: 2.827e-03 (2e-15), 7: 9.786e-04 (2e-15), 8: 3.370e-04 (2e-15), 9: 1.156e-04 (1e-14), 10: 3.956e-05 (3e-14) | 2.9231 |
| 0.8 | 9.0000 | 1: 1.567e-01 (0e+00), 2: 2.010e-02 (2e-16), 3: 2.447e-03 (1e-15), 4: 2.908e-04 (2e-15), 5: 3.406e-05 (9e-15), 6: 3.953e-06 (6e-13), 7: 4.558e-07 (6e-13), 8: 5.231e-08 (3e-11), 9: 5.981e-09 (7e-11), 10: 6.819e-10 (1e-10) | 8.7717 |

The tail's share Σ_{n>2K}/(first term) against ₂F₁(−½, −½; 1; t²) − 1 (the share's limit if the ratios' leading coefficients govern the tail):

| τ/σ | t | K = 5 | 10 | 20 | 30 | ₂F₁(−½, −½; 1; t²) − 1 |
|---|---|---|---|---|---|---|
| 0.25 | 0.01515 | 0.00003 | 0.00004 | 0.00005 | 0.00005 | 0.00006 |
| 0.5 | 0.05573 | 0.00038 | 0.00054 | 0.00065 | 0.00069 | 0.00078 |
| 1 | 0.17157 | 0.00367 | 0.00518 | 0.00617 | 0.00655 | 0.00737 |
| 2 | 0.38197 | 0.02064 | 0.02642 | 0.03074 | 0.03251 | 0.03682 |

## B. The transit

The sliding path (polished; a kernel leaves only if its correlation stays inward with its own contribution removed) on f* = 0.8φ_σ*(· − 0.03) sampled at 321 points of [-1.6, 1.6], 21 kernels started at the grid −1, −0.9, …, 1, to 13 entries, the fast flow's chunks up to 10⁵·2^24: after each entry the relaxed block's relative rms error against the optimal fit with as many kernels on the same samples (the better of variable projection from the block and from the half-squeeze nodes: the best local optimum found); a block is optimal within 0.001. The numerical variants: base (polish step 1σ, grid shifted by 0, the polish's budget 5000 evaluations, its gradient test 1e-15), polish 2σ (polish step 2σ, grid shifted by 0, the polish's budget 5000 evaluations, its gradient test 1e-15), grid +0.05 (polish step 1σ, grid shifted by 0.05, the polish's budget 5000 evaluations, its gradient test 1e-15), polish 2σ, 5·10⁴ evaluations (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, its gradient test 1e-15), polish 2σ, 5·10⁴, repeated (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, repeated up to 3 times, its gradient test 1e-15), polish 2σ, 5·10⁴, repeated, no gradient test (polish step 2σ, grid shifted by 0, the polish's budget 50000 evaluations, repeated up to 3 times, no gradient test). 'Relaxed': the block passed the stationarity test. The first non-optimal block's polishes as status/evaluations/accepted (status 0: the budget was spent first; 1: the gradient test). 'none': every block optimal, to the path's cap of 12 kernels (the range censored).

| τ/σ | variant | optimal blocks (kernels) | the first non-optimal: kernels, its optimum's error, block / optimum, relaxed | its polishes | relaxed blocks | chunks redone (stall, failure), failed | polish rejected | exits |
|---|---|---|---|---|---|---|---|---|
| 0.25 | base | 1, 2, 3, 4 | 5, 1.61e-09, 10.85, True | 0/5000/a | 12/12 | 0, 0, 0 | 13 | 1 |
| 0.25 | polish 2σ | 1, 2, 3, 4 | 5, 1.61e-09, 10.95, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 2 |
| 0.25 | grid +0.05 | 1, 2, 3, 4 | 5, 1.61e-09, 10.96, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 2 |
| 0.25 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4 | 5, 1.61e-09, 1.012, True | 1/30240/a | 12/12 | 0, 0, 0 | 2 | 1 |
| 0.25 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4 | 5, 1.61e-09, 1.012, True | 1/30240/a 1/1/r | 12/12 | 0, 0, 0 | 2 | 0 |
| 0.25 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5 | 6, 2.55e-11, 4.951, True | 0/50000/a 0/50000/a 0/50000/a 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.5 | base | 1, 2, 3, 4, 5 | 6, 6.30e-08, 14.34, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 3, 5, 0 | 182 | 0 |
| 0.5 | polish 2σ | 1, 2, 3, 4, 5 | 6, 6.30e-08, 1.46, True | 0/5000/a | 12/12 | 0, 0, 0 | 3 | 1 |
| 0.5 | grid +0.05 | 1, 2, 3, 4, 5 | 6, 6.30e-08, 14.19, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 2, 5, 0 | 143 | 0 |
| 0.5 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 4.605, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.5 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 2.677, True | 0/50000/a 0/50000/a 1/1/r | 12/12 | 0, 0, 0 | 1 | 0 |
| 0.5 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7 | 8, 2.10e-10, 1.233, True | 0/50000/a 0/50000/a 0/50000/a 0/50000/a | 12/12 | 1, 0, 0 | 0 | 0 |
| 0.75 | base | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 7.53, True | 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r 0/5000/r | 12/12 | 1, 7, 0 | 138 | 0 |
| 0.75 | polish 2σ | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 1.642, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7 | 8, 5.22e-08, 8.656, True | 0/5000/r | 12/12 | 1, 3, 0 | 88 | 0 |
| 0.75 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 6.81e-10, 2.305, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 6.81e-10, 1.121, True | 0/50000/a 0/50000/a 1/1/r | 12/12 | 0, 0, 0 | 0 | 0 |
| 0.75 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 | 11, 7.75e-11, 2.685, True | 0/50000/a 0/50000/a 0/50000/a 0/50000/a | 12/12 | 1, 0, 0 | 0 | 0 |
| 1 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 5.22e-08, 1.771, True | 0/5000/a | 12/12 | 0, 0, 0 | 12 | 0 |
| 1 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9 | 10, 5.22e-08, 1.771, True | 0/5000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 10 | 9, 2.96e-07, 4.241, True | 1/3635/r 1/3620/r 1/3720/r 1/3725/r 1/3618/r 1/3600/r 1/3648/r 1/3682/r 1/3722/r 1/3604/r 1/3685/r 1/3689/r 1/3684/r 1/3608/r 1/3624/r 1/3594/r 1/3881/r 1/3658/r 1/4129/r 1/3699/r 1/3648/r 1/3718/r 1/3616/r 1/3553/r | 12/12 | 2, 3, 0 | 110 | 0 |
| 1 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 | 12, 1.61e-09, 1.465, True | 0/50000/a | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 1, 0 | 71 | 0 |
| 1.5 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 1.5 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | base | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 13 | 0 |
| 2 | polish 2σ | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | grid +0.05 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴ evaluations | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴, repeated | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |
| 2 | polish 2σ, 5·10⁴, repeated, no gradient test | 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | none | — | 12/12 | 0, 0, 0 | 0 | 0 |

Over every width and variant: the optimal blocks' optimum errors down to 6.81e-10; the non-optimal blocks' from 2.96e-07 down (0 of 116 failed the stationarity test; 53 of those that passed it pass it at any state by Cauchy–Schwarz, the test being absolute — their errors up to 4.0e-10); non-optimal blocks with an optimum above 10⁻⁷: 1; the non-optimal blocks' errors 1.01e-12 … 1.26e-06.

Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in 33 of the 38 blocks whose optimum is at least 1e-10; 1 ends within 1% of it (1.001 times: in its basin, unconverged); 4 end 16.36–24.74 times above it, the block's nodes up to 2.14–3.99σ from the optimum's. The 78 blocks with an optimum below 1e-10, where the solve stops short of the optimum (A), are not counted (0 reached). From the first non-optimal block of each path: 20 of 20; on 2 paths the first non-optimal block's optimum is below 1e-10. There the Jacobian's condition number is 1.3e+06–5e+08 for the joint fit in (w, c), 5.4e+03–3.6e+05 for variable projection's in c alone.


The entries over the variants that record them, by the block they join: the entrant had the largest current rate |⟨rφ_j⟩|; its rate over the largest; its dual |u_j| at the block's start; whether it starts beyond the block's outermost node (how far) or between its nodes (their duals). The entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩| (`sliding_iss`'s rule), its dual charged during the earlier blocks.

| joining | entries | largest rate | rate / largest: median (range) | dual: median (least) | beyond: entrants, median (range) | between: entrants (duals) |
|---|---|---|---|---|---|---|
| none (the first entry) | 36 | 36 | 1.00 (1.000–1.00) | 0.000 (0.000) | — | — |
| an optimal block | 316 | 5 | 0.21 (0.001–1.00) | 0.959 (0.570) | 314, 1.22σ (0.11–2.75σ) | 2 (0.989–0.992) |
| a non-optimal block | 116 | 39 | 0.71 (0.004–1.00) | 0.786 (0.020) | 97, 2.50σ (0.17–4.38σ) | 19 (0.237–0.999) |
| all | 468 | 80 | 0.30 (0.001–1.00) | 0.947 (0.000) | 411, 1.42σ (0.11–4.38σ) | 21 (0.237–0.999) |

### τ/σ = 0.25, base (events: 10+, 11+, 8+, 12+, 13+, 7+, 16+, 14+, 4+, 5+, 4−, 17+, 3+, 9+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 18.87 | 2.1430e-02 | 2.1430e-02 | 1.000000 | True | 0.0e+00 | 0.000 | 2.1431e-02 | — |
| 2 | 529.4 | 3.7502e-04 | 3.7502e-04 | 1.000000 | True | 4.2e-11 | 0.000 | 3.7504e-04 | 57.143 |
| 3 | 5.093e+04 | 6.2260e-06 | 6.2260e-06 | 1.000000 | True | 2.2e-09 | 0.000 | 6.2262e-06 | 60.235 |
| 4 | 2.389e+05 | 1.0087e-07 | 1.0087e-07 | 1.000000 | True | 3.9e-07 | 0.000 | 1.0087e-07 | 61.723 |
| 5 | 6.067e+07 | 1.7479e-08 | 1.6114e-09 | 10.847501 | True | 3.1e-01 | 0.314 | 1.6114e-09 | 5.771 |
| 6 | 6.766e+08 | 2.4073e-09 | 2.5507e-11 | 94.379669 | True | 6.5e-01 | 0.649 | 2.5507e-11 | 7.261 |
| 7 | 2.499e+10 | 1.2145e-07 | 4.0115e-13 | 302749.882294 | True | 3.8e+00 | 3.826 | 4.0116e-13 | 0.020 |
| 8 | 2.521e+10 | 3.5753e-10 | 6.2753e-15 | 56974.096023 | True | 1.3e+00 | 1.321 | 6.2929e-15 | 339.685 |
| 9 | 1.021e+11 | 2.8072e-10 | 1.7429e-16 | 1610643.059713 | True | 4.1e+00 | 4.100 | 1.6856e-16 | 1.274 |
| 9 | 2.267e+11 | 2.3861e-10 | 1.7429e-16 | 1369005.248164 | True | 3.3e+00 | 3.290 | 1.6856e-16 | — |
| 10 | 4.159e+11 | 1.7099e-10 | 2.2206e-16 | 770019.312199 | True | 4.4e+00 | 4.355 | 2.9816e-16 | 1.395 |
| 11 | 7.292e+11 | 1.3883e-10 | 1.6488e-16 | 842053.741460 | True | 4.8e+00 | 4.789 | 2.9420e-16 | 1.232 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4): min C/max C within ±1e-07, its peaks beyond the outermost nodes at (left, right) (2.02, 2.02)σ, (2.61, 2.61)σ, (3.08, 3.08)σ, (3.48, 3.48)σ.

### τ/σ = 0.5, base (events: 10+, 11+, 8+, 13+, 7+, 14+, 6+, 17+, 3+, 2+, 19+, 4+, 18+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 18.15 | 7.8689e-02 | 7.8689e-02 | 1.000000 | True | 0.0e+00 | 0.000 | 7.8750e-02 | — |
| 2 | 141 | 5.0651e-03 | 5.0651e-03 | 1.000000 | True | 6.2e-10 | 0.000 | 5.0678e-03 | 15.535 |
| 3 | 3729 | 3.0925e-04 | 3.0925e-04 | 1.000000 | True | 2.7e-09 | 0.001 | 3.0939e-04 | 16.379 |
| 4 | 9700 | 1.8425e-05 | 1.8425e-05 | 1.000000 | True | 1.4e-08 | 0.000 | 1.8433e-05 | 16.784 |
| 5 | 6.597e+05 | 1.0824e-06 | 1.0824e-06 | 1.000000 | True | 2.9e-08 | 0.000 | 1.0828e-06 | 17.023 |
| 6 | 1.07e+07 | 9.0368e-07 | 6.3003e-08 | 14.343496 | True | 2.0e+00 | 1.976 | 6.3030e-08 | 1.198 |
| 7 | 1.518e+07 | 3.0205e-07 | 3.6437e-09 | 82.896419 | True | 1.2e+00 | 1.228 | 3.6452e-09 | 2.992 |
| 8 | 3.371e+08 | 2.9130e-07 | 2.0972e-10 | 1388.997007 | True | 3.6e+00 | 3.623 | 2.0980e-10 | 1.037 |
| 9 | 3.763e+08 | 2.8760e-07 | 1.2031e-11 | 23904.519332 | True | 3.9e+00 | 3.949 | 1.2031e-11 | 1.013 |
| 10 | 7.815e+08 | 2.5457e-07 | 6.8789e-13 | 370068.553314 | True | 4.6e+00 | 4.599 | 6.8790e-13 | 1.130 |
| 11 | 9.679e+08 | 2.4758e-07 | 3.9238e-14 | 6309816.493600 | True | 4.8e+00 | 4.789 | 3.9246e-14 | 1.028 |
| 12 | 1.101e+09 | 2.4274e-07 | 2.2524e-15 | 107768417.576135 | True | 4.6e+00 | 4.639 | 2.2392e-15 | 1.020 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5): min C/max C within ±9e-10, its peaks beyond the outermost nodes at (left, right) (2.06, 2.06)σ, (2.47, 2.47)σ, (2.81, 2.81)σ, (3.09, 3.09)σ, (3.34, 3.34)σ.

### τ/σ = 0.75, base (events: 10+, 11+, 13+, 8+, 7+, 14+, 15+, 5+, 16+, 19+, 4+, 2+, 18+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 17.29 | 1.5617e-01 | 1.5617e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 1.5666e-01 | — |
| 2 | 69.35 | 2.0060e-02 | 2.0060e-02 | 1.000000 | True | 0.0e+00 | 0.003 | 2.0103e-02 | 7.785 |
| 3 | 927.4 | 2.4427e-03 | 2.4427e-03 | 1.000000 | True | 1.2e-09 | 0.003 | 2.4474e-03 | 8.212 |
| 4 | 1281 | 2.9023e-04 | 2.9023e-04 | 1.000000 | True | 1.4e-08 | 0.003 | 2.9075e-04 | 8.416 |
| 5 | 3.709e+04 | 3.3999e-05 | 3.3999e-05 | 1.000000 | True | 0.0e+00 | 0.003 | 3.4058e-05 | 8.537 |
| 6 | 4.95e+05 | 3.9461e-06 | 3.9461e-06 | 1.000000 | True | 1.8e-08 | 0.003 | 3.9529e-06 | 8.616 |
| 7 | 2.59e+06 | 4.5505e-07 | 4.5505e-07 | 1.000000 | True | 2.1e-07 | 0.002 | 4.5582e-07 | 8.672 |
| 8 | 6.861e+06 | 3.9327e-07 | 5.2224e-08 | 7.530460 | True | 1.8e+00 | 1.838 | 5.2311e-08 | 1.157 |
| 9 | 1.494e+07 | 2.9214e-07 | 5.9713e-09 | 48.924738 | True | 1.9e+00 | 1.851 | 5.9811e-09 | 1.346 |
| 10 | 8.66e+08 | 2.8715e-07 | 6.8078e-10 | 421.802203 | True | 4.0e+00 | 3.988 | 6.8186e-10 | 1.017 |
| 11 | 9.907e+08 | 2.8108e-07 | 7.7548e-11 | 3624.569070 | True | 3.8e+00 | 3.765 | 7.7548e-11 | 1.022 |
| 12 | 1.16e+09 | 2.1206e-07 | 8.8020e-12 | 24091.841122 | True | 3.5e+00 | 3.475 | 8.8021e-12 | 1.325 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7): min C/max C within ±4e-08, its peaks beyond the outermost nodes at (left, right) (2.13, 2.12)σ, (2.40, 2.40)σ, (2.62, 2.62)σ, (2.83, 2.83)σ, (3.01, 3.01)σ, (3.18, 3.18)σ, (3.34, 3.34)σ.

### τ/σ = 1, base (events: 10+, 11+, 13+, 7+, 14+, 6+, 15+, 5+, 16+, 4+, 17+, 3+, 2+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 16.51 | 2.3915e-01 | 2.3915e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 2.4090e-01 | — |
| 2 | 44.49 | 4.7497e-02 | 4.7497e-02 | 1.000000 | True | 0.0e+00 | 0.008 | 4.7745e-02 | 5.035 |
| 3 | 374 | 8.9365e-03 | 8.9365e-03 | 1.000000 | True | 2.9e-09 | 0.009 | 8.9780e-03 | 5.315 |
| 4 | 643.1 | 1.6402e-03 | 1.6402e-03 | 1.000000 | True | 0.0e+00 | 0.009 | 1.6474e-03 | 5.449 |
| 5 | 5788 | 2.9676e-04 | 2.9676e-04 | 1.000000 | True | 5.1e-09 | 0.008 | 2.9802e-04 | 5.527 |
| 6 | 5.893e+04 | 5.3197e-05 | 5.3197e-05 | 1.000000 | True | 1.5e-07 | 0.008 | 5.3418e-05 | 5.579 |
| 7 | 1.212e+05 | 9.4739e-06 | 9.4739e-06 | 1.000000 | True | 8.4e-08 | 0.007 | 9.5126e-06 | 5.615 |
| 8 | 4.83e+05 | 1.6791e-06 | 1.6791e-06 | 1.000000 | True | 0.0e+00 | 0.007 | 1.6859e-06 | 5.642 |
| 9 | 5.348e+06 | 2.9649e-07 | 2.9649e-07 | 1.000000 | True | 0.0e+00 | 0.007 | 2.9767e-07 | 5.663 |
| 10 | 5.839e+07 | 9.2441e-08 | 5.2197e-08 | 1.771002 | True | 2.7e-01 | 0.260 | 5.2404e-08 | 3.207 |
| 11 | 8.35e+07 | 1.6136e-08 | 9.1674e-09 | 1.760191 | True | 2.6e-01 | 0.255 | 9.2035e-09 | 5.729 |
| 12 | 5.188e+08 | 6.0433e-08 | 1.6069e-09 | 37.609199 | True | 1.5e+00 | 1.458 | 1.6131e-09 | 0.267 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9): min C/max C within ±1e-08, its peaks beyond the outermost nodes at (left, right) (2.21, 2.21)σ, (2.38, 2.38)σ, (2.53, 2.53)σ, (2.67, 2.67)σ, (2.80, 2.80)σ, (2.92, 2.92)σ, (3.04, 3.04)σ, (3.14, 3.14)σ, (3.25, 3.25)σ.

### τ/σ = 1.5, base (events: 10+, 11+, 7+, 14+, 6+, 15+, 5+, 16+, 3+, 18+, 2+, 19+, 1+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 15.37 | 3.8940e-01 | 3.8940e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 3.9732e-01 | — |
| 2 | 26.96 | 1.2942e-01 | 1.2942e-01 | 1.000000 | True | 0.0e+00 | 0.032 | 1.3149e-01 | 3.009 |
| 3 | 124.3 | 4.0725e-02 | 4.0725e-02 | 1.000000 | True | 0.0e+00 | 0.035 | 4.1300e-02 | 3.178 |
| 4 | 172.1 | 1.2492e-02 | 1.2492e-02 | 1.000000 | True | 0.0e+00 | 0.034 | 1.2657e-02 | 3.260 |
| 5 | 1003 | 3.7762e-03 | 3.7762e-03 | 1.000000 | True | 2.1e-08 | 0.033 | 3.8240e-03 | 3.308 |
| 6 | 3159 | 1.1307e-03 | 1.1307e-03 | 1.000000 | True | 2.5e-08 | 0.031 | 1.1446e-03 | 3.340 |
| 7 | 7097 | 3.3631e-04 | 3.3631e-04 | 1.000000 | True | 1.8e-08 | 0.030 | 3.4036e-04 | 3.362 |
| 8 | 1.718e+04 | 9.9542e-05 | 9.9542e-05 | 1.000000 | True | 9.0e-08 | 0.029 | 1.0072e-04 | 3.379 |
| 9 | 1.939e+05 | 2.9351e-05 | 2.9351e-05 | 1.000000 | True | 1.2e-07 | 0.028 | 2.9693e-05 | 3.391 |
| 10 | 2.05e+05 | 8.6284e-06 | 8.6284e-06 | 1.000000 | True | 2.6e-07 | 0.027 | 8.7280e-06 | 3.402 |
| 11 | 7.793e+05 | 2.5303e-06 | 2.5303e-06 | 1.000000 | True | 1.1e-07 | 0.026 | 2.5593e-06 | 3.410 |
| 12 | 3.246e+06 | 7.4053e-07 | 7.4053e-07 | 1.000000 | True | 1.3e-07 | 0.025 | 7.4893e-07 | 3.417 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12): min C/max C within ±5e-07, its peaks beyond the outermost nodes at (left, right) (2.38, 2.38)σ, (2.41, 2.41)σ, (2.47, 2.47)σ, (2.53, 2.53)σ, (2.60, 2.60)σ, (2.66, 2.66)σ, (2.72, 2.72)σ, (2.77, 2.77)σ, (2.84, 2.84)σ, (2.89, 2.89)σ, (2.94, 2.94)σ, (2.99, 2.99)σ.

### τ/σ = 2, base (events: 10+, 11+, 7+, 14+, 5+, 16+, 4+, 17+, 18+, 3+, 20+, 1+, 0+)

| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | half-squeeze on ℝ | ratio per kernel |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 14.69 | 5.0462e-01 | 5.0462e-01 | 1.000000 | True | 0.0e+00 | 0.000 | 5.2280e-01 | — |
| 2 | 20.93 | 2.2364e-01 | 2.2364e-01 | 1.000000 | True | 3.9e-09 | 0.069 | 2.3059e-01 | 2.256 |
| 3 | 65.67 | 9.4014e-02 | 9.4014e-02 | 1.000000 | True | 0.0e+00 | 0.079 | 9.6611e-02 | 2.379 |
| 4 | 69.8 | 3.8509e-02 | 3.8509e-02 | 1.000000 | True | 1.9e-09 | 0.079 | 3.9497e-02 | 2.441 |
| 5 | 426 | 1.5539e-02 | 1.5539e-02 | 1.000000 | True | 8.8e-09 | 0.077 | 1.5918e-02 | 2.478 |
| 6 | 430.4 | 6.2092e-03 | 6.2092e-03 | 1.000000 | True | 5.4e-09 | 0.074 | 6.3555e-03 | 2.502 |
| 7 | 1616 | 2.4644e-03 | 2.4644e-03 | 1.000000 | True | 3.9e-09 | 0.071 | 2.5208e-03 | 2.520 |
| 8 | 5276 | 9.7318e-04 | 9.7318e-04 | 1.000000 | True | 4.0e-08 | 0.068 | 9.9499e-04 | 2.532 |
| 9 | 8671 | 3.8282e-04 | 3.8282e-04 | 1.000000 | True | 9.9e-08 | 0.066 | 3.9125e-04 | 2.542 |
| 10 | 1.156e+04 | 1.5013e-04 | 1.5013e-04 | 1.000000 | True | 1.7e-07 | 0.064 | 1.5339e-04 | 2.550 |
| 11 | 1.214e+05 | 5.8730e-05 | 5.8730e-05 | 1.000000 | True | 2.5e-07 | 0.062 | 5.9987e-05 | 2.556 |
| 12 | 1.231e+05 | 2.2927e-05 | 2.2927e-05 | 1.000000 | True | 1.7e-07 | 0.060 | 2.3413e-05 | 2.562 |

The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12): min C/max C within ±1e-06, its peaks beyond the outermost nodes at (left, right) (2.57, 2.57)σ, (2.51, 2.51)σ, (2.51, 2.51)σ, (2.53, 2.53)σ, (2.56, 2.56)σ, (2.58, 2.58)σ, (2.60, 2.60)σ, (2.63, 2.63)σ, (2.66, 2.66)σ, (2.70, 2.70)σ, (2.72, 2.72)σ, (2.75, 2.75)σ.

## C. Under noise

τ/σ = 1, 10 draws per level, the path to 11 entries: the discrepancy rule's stop (kernels), its prediction from (i) and the predicted error there, the stopped fit's relative error (median [range]), the kernels within 4σ* of the feature at the stop, the largest relative gap's kernel count (in how many draws it is the stop's).

| s | stop: kernels (median [range]) | predicted | error at the stop: median [range] | predicted | kernels at the feature | largest gap at |
|---|---|---|---|---|---|---|
| 0.0001 | 6 [6, 6] (10/10 stopped) | 6 | 8.712e-05 [7.949e-05, 1.053e-04] | 9.533e-05 | 6 [6, 6] | 5 [4, 9] (the stop in 0) |
| 0.001 | 4 [4, 5] (10/10 stopped) | 4 | 1.712e-03 [5.808e-04, 1.813e-03] | 1.769e-03 | 4 [4, 5] | 4 [2, 5] (the stop in 4) |
| 0.01 | 3 [3, 4] (10/10 stopped) | 3 | 1.079e-02 [3.814e-03, 1.274e-02] | 1.057e-02 | 3 [3, 4] | 2 [2, 4] (the stop in 1) |

At the stop (K kernels): the error over the oracle's s(2K/N)^{1/2} (the variance alone, the right support), over the bias–variance prediction (e_K²‖f*‖²_N + 2Ks²/N)^{1/2} — e_K the noiseless optimum's error with K kernels (B's base path at τ/σ = 1) — and over the discrepancy rule's bound s((1 + δ)² − 1 + 4K/N)^{1/2}; paths with a runaway kernel before the stop:

| s | error / oracle: median [range] | / bias–variance | / bound: max | runaway |
|---|---|---|---|---|
| 0.0001 | 1.10 [1.01, 1.33] | 0.915 [0.835, 1.106] | 0.483 | 0 |
| 0.001 | 2.66 [0.81, 2.81] | 0.990 [0.745, 1.348] | 0.871 | 0 |
| 0.01 | 1.93 [0.59, 2.28] | 1.024 [0.573, 1.209] | 0.627 | 0 |

The stop as a margin (`stop_law`): the first-term formula's root K*(s) of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/(s(2δ + δ²)^{1/2}))/log(1/t), whose ceiling is the predicted stop, and the probability of the stop at each K — the rms's expectation s(b_K² + 1 − 2K/N)^{1/2} and its scatter from the noise's chi-square and its cross term with the bias, b_K = e_K‖f*‖_N/s, e_K the half-squeeze rule's error, in the normal approximation (the rms falls along the path) —; the expected number of the draws stopping at the predicted K against the observed:

| s | K*(s) | P(stop at K) | at the predicted K | expected | observed |
|---|---|---|---|---|---|
| 0.0001 | 5.26 | 5: 0.005, 6: 0.993, 7: 0.001 | 0.993 | 9.9 of 10 | 10 |
| 0.001 | 3.91 | 4: 0.789, 5: 0.210, 6: 0.001 | 0.789 | 7.9 of 10 | 7 |
| 0.01 | 2.55 | 3: 0.986, 4: 0.012 | 0.986 | 9.9 of 10 | 9 |

## D. The tail of the exact series (the sharp constant, conjectural)

The half-squeeze rule's squared error over the series' first term, Σ_j t^{2j}(2K)!/(2K + 2j)!·(D_{2K+2j,K}/D_{2K,K})² (the defects exact, by the remainders modulo He_K), against ₂F₁(−½, −½; 1; t²), the sum of the terms' limits binom(1/2, j)²t^{2j}; K times the gap in parentheses:

| t | ₂F₁ | K = 10 | K = 20 | K = 40 | K = 80 | K = 160 | K = 320 |
|---|---|---|---|---|---|---|---|
| 0.1716 | 1.007373 | 1.005184 (-0.022) | 1.006173 (-0.024) | 1.006743 (-0.025) | 1.007050 (-0.026) | 1.007209 (-0.026) | 1.007291 (-0.026) |
| 0.382 | 1.036820 | 1.026417 (-0.104) | 1.030739 (-0.122) | 1.033484 (-0.133) | 1.035063 (-0.141) | 1.035916 (-0.145) | 1.036361 (-0.147) |
| 0.5 | 1.063544 | 1.047067 (-0.165) | 1.053239 (-0.206) | 1.057594 (-0.238) | 1.060295 (-0.260) | 1.061834 (-0.274) | 1.062665 (-0.281) |
| 0.8 | 1.167810 | 1.162342 (-0.055) | 1.156639 (-0.223) | 1.157221 (-0.424) | 1.159788 (-0.642) | 1.162436 (-0.860) | 1.164502 (-1.058) |

The ratios D_{2K+2j,K}/D_{2K,K} are polynomials in K of degree j (the (j + 1)-th divided difference over K = 1 … j + 2: 0, 0, 0, 0, 0, 0, 0) with the leading coefficients -1, -1/2, -1/2, -5/8, -7/8, -21/16, -33/16 at j = 1 … 7 — −C_{j−1}/2^{j−1} at every j (C the Catalan numbers), binom(1/2, j)(−2)^j.

Numbers are produced by `experiments/exp_width_mismatch.py`; nothing is edited by hand.
