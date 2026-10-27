# Proposition D‴.43: the colour of the measurement noise selects the kernels (`exp_colour_selection.py`)

w = u⊙u − v⊙v on a fixed dictionary Φ (N × d), fitted by gradient descent with label noise; along the exact fits the reduced flow dP/ds = −ω⊙P, dw/ds = −Π_D(ω⊙w), ω_i = φ_iᵀf(K)φ_i, the gradient flow of R = tr Ψ(K) (`colour_selection`). Success: within 0.0001·max|w*| of the truth. The start: the noiseless gradient flow's interpolant from u = v = 0.5 (the hyperbolic entropy's minimiser on M, P = α²); the end interpolants likewise from each path's end state.

## A. The potential's identities

| N × d | τ | det(I + τK) against its subset expansion (relative) | ∂R/∂D against (4/N)ω (relative) | R's second differences along 20 random segments (max, min) | white: |second difference| |
|---|---|---|---|---|---|
| 4 × 7 | 0.3 | 0.0e+00 | 5.6e-10 | -7.12e-03, -2.54e-01 | 1.8e-15 |
| 4 × 7 | 10 | 7.8e-16 | 1.6e-09 | -3.25e-03, -7.93e-02 | 3.6e-15 |
| 4 × 7 | 1000 | 2.6e-15 | 2.8e-09 | -8.03e-05, -2.18e-03 | 8.9e-16 |
| 5 × 9 | 0.3 | 4.4e-16 | 7.4e-10 | -2.29e-02, -3.87e-01 | 0.0e+00 |
| 5 × 9 | 10 | 8.9e-16 | 2.5e-09 | -2.16e-02, -1.29e-01 | 3.6e-15 |
| 5 × 9 | 1000 | 3.3e-16 | 3.9e-09 | -1.66e-04, -1.62e-03 | 0.0e+00 |

On a k-sparse exact fit (P = 0, D = |w|): τ_cR − k log(4τ_c/N) against log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) —

| N × d | k | the limit | τ_c = 10², 10⁴, 10⁶ |
|---|---|---|---|
| 10 × 40 | 3 | -0.065900 | +0.012561, -0.065104, -0.065892 |
| 10 × 40 | 6 | -1.168243 | -0.941086, -1.165909, -1.168220 |
| 20 × 80 | 5 | -1.141188 | -0.788862, -1.137498, -1.141151 |

## B. Label-noise descent against the reduced flow (10 × 40, a 4-sparse truth, step 0.5, σ = 0.1, 8 paths)

| seed | the noise | steps | s at the end | the paths' mean against the flow (max over w, at s = 0, S/8, …, S) | the flow's error at S | the paths' mean's error at S | the end interpolants' errors (per path) | their supports' sizes | basis pursuit's error |
|---|---|---|---|---|---|---|---|---|---|
| 0 | fresh | 3000000 | 300 | 0.000, 0.021, 0.018, 0.020, 0.028, 0.019, 0.019, 0.032, 0.033 | 2.4e+00 | 2.4e+00 | 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00 | 13, 13, 13, 12, 13, 13, 13, 13 | 2.98 |
| 0 | AR(1), τ_c = 10 | 1500000 | 150 | 0.000, 0.035, 0.053, 0.021, 0.078, 0.061, 0.041, 0.015, 0.049 | 2.9e-10 | 4.9e-02 | 4e-10, 3e-10, 4e-10, 4e-10, 4e-10, 4e-10, 4e-10, 4e-10 | 4, 4, 4, 4, 4, 4, 4, 4 | 2.98 |
| 0 | AR(1), τ_c = 100 | 4000000 | 400 | 0.000, 0.018, 0.051, 0.206, 0.046, 0.063, 0.051, 0.053, 0.040 | 5.1e-12 | 4.0e-02 | 5e-13, 7e-13, 4e-13, 1e-13, 6e-13, 5e-13, 2e-13, 4e-14 | 4, 4, 4, 4, 4, 4, 4, 4 | 2.98 |
| 1 | fresh | 3000000 | 300 | 0.000, 0.016, 0.007, 0.008, 0.007, 0.009, 0.021, 0.022, 0.022 | 1.5e+00 | 1.5e+00 | 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00 | 10, 10, 10, 10, 10, 10, 10, 10 | 1.51 |
| 1 | AR(1), τ_c = 10 | 1500000 | 150 | 0.000, 0.064, 0.130, 0.143, 0.248, 0.212, 0.065, 0.049, 0.041 | 1.5e+00 | 1.5e+00 | 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00, 2e+00 | 11, 11, 11, 11, 10, 10, 11, 11 | 1.51 |
| 1 | AR(1), τ_c = 100 | 4000000 | 400 | 0.000, 0.020, 0.115, 0.161, 0.320, 0.147, 0.153, 0.098, 0.075 | 6.4e-12 | 7.5e-02 | 2e-10, 5e-11, 1e-10, 1e-10, 2e-10, 4e-11, 1e-10, 4e-11 | 4, 4, 4, 4, 4, 4, 4, 4 | 1.51 |
| 2 | fresh | 3000000 | 300 | 0.000, 0.016, 0.023, 0.013, 0.016, 0.009, 0.009, 0.025, 0.014 | 1.3e+00 | 1.3e+00 | 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00 | 10, 10, 10, 10, 10, 10, 10, 10 | 1.32 |
| 2 | AR(1), τ_c = 10 | 1500000 | 150 | 0.000, 0.087, 0.096, 0.052, 0.050, 0.024, 0.034, 0.121, 0.056 | 9.3e-13 | 5.6e-02 | 2e-13, 1e-13, 5e-13, 8e-13, 5e-13, 9e-13, 1e-13, 1e-13 | 4, 4, 4, 4, 4, 4, 4, 4 | 1.32 |
| 2 | AR(1), τ_c = 100 | 4000000 | 400 | 0.000, 0.033, 0.031, 0.254, 0.062, 0.065, 0.060, 0.163, 0.078 | 1.7e-11 | 7.8e-02 | 2e-13, 3e-13, 9e-14, 1e-13, 1e-13, 3e-13, 1e-13, 3e-13 | 4, 4, 4, 4, 4, 4, 4, 4 | 1.32 |
| 3 | fresh | 3000000 | 300 | 0.000, 0.014, 0.014, 0.014, 0.016, 0.026, 0.016, 0.014, 0.016 | 1.3e+00 | 1.3e+00 | 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00 | 10, 10, 10, 10, 10, 10, 10, 10 | 1.33 |
| 3 | AR(1), τ_c = 10 | 1500000 | 150 | 0.000, 0.058, 0.057, 0.078, 0.068, 0.076, 0.064, 0.067, 0.056 | 1.1e+00 | 1.1e+00 | 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00, 1e+00 | 10, 10, 11, 10, 10, 10, 10, 11 | 1.33 |
| 3 | AR(1), τ_c = 100 | 4000000 | 400 | 0.000, 0.026, 0.100, 0.142, 0.143, 0.090, 0.058, 0.034, 0.084 | 1.7e-11 | 8.4e-02 | 3e-13, 5e-13, 1e-13, 3e-13, 2e-13, 2e-14, 1e-13, 1e-13 | 4, 4, 4, 4, 4, 4, 4, 4 | 1.33 |

## C. Recovery against sparsity (the reduced flows, continuum filters)

| N × d | k | instances | bp | white | slow 10 | slow 100 | slow 1000 | white, then slow 100 | the slow flow's (τ_c = 100) mean support size | basis pursuit's |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 × 40 | 1 | 40 | 0.97 | 0.97 | 1.00 | 1.00 | 1.00 | 0.97 | 1.0 | 1.2 |
| 10 × 40 | 2 | 40 | 0.90 | 0.88 | 0.95 | 0.97 | 1.00 | 0.90 | 2.2 | 2.8 |
| 10 × 40 | 3 | 40 | 0.38 | 0.35 | 0.62 | 0.75 | 0.78 | 0.40 | 4.8 | 7.3 |
| 10 × 40 | 4 | 40 | 0.17 | 0.17 | 0.42 | 0.60 | 0.62 | 0.20 | 6.4 | 8.9 |
| 10 × 40 | 5 | 40 | 0.00 | 0.00 | 0.07 | 0.15 | 0.17 | 0.00 | 9.2 | 10.0 |
| 10 × 40 | 6 | 40 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 10.1 | 10.0 |
| 20 × 80 | 2 | 20 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 2.0 | 2.0 |
| 20 × 80 | 4 | 20 | 0.90 | 0.90 | 0.95 | 1.00 | 1.00 | 0.95 | 4.0 | 5.6 |
| 20 × 80 | 6 | 20 | 0.35 | 0.35 | 0.65 | 0.75 | 0.75 | 0.40 | 9.4 | 15.1 |
| 20 × 80 | 8 | 20 | 0.10 | 0.10 | 0.20 | 0.35 | 0.30 | 0.10 | 15.6 | 18.6 |
| 20 × 80 | 10 | 20 | 0.05 | 0.05 | 0.00 | 0.10 | 0.10 | 0.05 | 18.9 | 19.4 |

## D. Gaussian kernel dictionaries on [0, 1] (40 centres, 20 equispaced samples), random k-atom truths (40 per row)

| width | the largest off-diagonal |φ_iᵀφ_j| | k | bp | white | slow 10 | slow 100 | mean support size: basis pursuit, slow 100 |
|---|---|---|---|---|---|---|---|
| 0.03 | 0.90 | 3 | 0.85 | 0.85 | 0.90 | 0.95 | 3.7, 3.2 |
| 0.03 | 0.90 | 4 | 0.80 | 0.80 | 0.80 | 0.95 | 5.0, 4.2 |
| 0.03 | 0.90 | 5 | 0.70 | 0.70 | 0.80 | 0.90 | 6.6, 5.3 |
| 0.03 | 0.90 | 6 | 0.45 | 0.45 | 0.55 | 0.68 | 8.8, 6.9 |
| 0.05 | 0.96 | 3 | 0.95 | 0.95 | 0.95 | 1.00 | 3.5, 3.0 |
| 0.05 | 0.96 | 4 | 0.85 | 0.85 | 0.93 | 1.00 | 5.3, 4.0 |
| 0.05 | 0.96 | 5 | 0.82 | 0.82 | 0.85 | 0.93 | 6.7, 5.6 |
| 0.05 | 0.96 | 6 | 0.65 | 0.65 | 0.60 | 0.78 | 9.2, 7.4 |

Label-noise descent on the first instances where basis pursuit fails (width 0.05, k = 4; 8 paths, σ = 0.1, the step min(0.5, 1/(2λ_max)) at the start):

| instance | step | the noise | steps | the end interpolants' errors (per path) | their supports' sizes | the ends' atoms alive (D > 10⁻⁸ max D) | the raw ends' distance to the truth, to basis pursuit's solution (max over paths) | the reduced flow's error, support, distance to basis pursuit's solution |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., 1e+00 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., 14 | 15, 14, 15, 15, 15, 14, 15, 14 | 1.20, 0.13 | 1.2e+00, 14, 2.2e-08 |
| 0 | 0.5 | AR(1), τ_c = 100 | 16000000 | 2e-12, 2e-12, 3e-12, 3e-13, 1e-12, 2e-12, 4e-13, 2e-12 | 4, 4, 4, 4, 4, 4, 4, 4 | 6, 9, 6, 6, 8, 8, 5, 5 | 0.78, 1.86 | 6.9e-12, 4, 1.2e+00 |
| 7 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | 14, 14, 14, 13, 15, 14, 14, 15 | 0.73, 0.10 | 7.3e-01, 18, 1.3e-09 |
| 7 | 0.5 | AR(1), τ_c = 100 | 16000000 | 5e-08, 5e-08, 6e-08, 4e-08, 4e-08, 5e-08, 6e-08, 4e-08 | 4, 4, 4, 4, 4, 4, 4, 4 | 7, 8, 9, 10, 9, 7, 10, 7 | 0.37, 0.79 | 1.8e-08, 4, 7.3e-01 |
| 13 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | 11, 11, 11, 11, 11, 11, 11, 11 | 1.07, 0.07 | 1.0e+00, 15, 3.8e-07 |
| 13 | 0.5 | AR(1), τ_c = 100 | 16000000 | 7e-01, 6e-01, 7e-01, 7e-01, 7e-01, 6e-01, 7e-01, 7e-01 | 12, 13, 13, 12, 13, 12, 13, 12 | 14, 14, 11, 16, 12, 13, 12, 14 | 1.24, 0.47 | 8.8e-03, 8, 1.0e+00 |

n.c.: the end's flow interpolant does not converge — the exact fit needs atoms whose masses have underflowed (fewer than N alive, a residual outside their span), its maximiser is at infinity; the raw end (residual O(σ)) is reported instead.

## E. Other colours at the step 0.5 (the reduced flow with each scheme's filter; 10 × 40, 20 instances)

| k | fresh | AR(1), τ_c = 10 | smooth (Gaussian), τ_c = 10 | kept 40 steps |
|---|---|---|---|---|
| 3 | 0.30 | 0.65 | 0.70 | 0.70 |
| 4 | 0.20 | 0.40 | 0.55 | 0.40 |

## F. The start: the interpolant of u = v = α (10 × 40)

| k | α | the interpolant itself | the white flow | the slow flow (τ_c = 100) |
|---|---|---|---|---|
| 3 | 0.1 | 0.00 | 0.40 | 0.70 |
| 3 | 0.5 | 0.00 | 0.40 | 0.75 |
| 3 | 2 | 0.00 | 0.40 | 0.00 |
| 4 | 0.1 | 0.00 | 0.25 | 0.35 |
| 4 | 0.5 | 0.00 | 0.25 | 0.50 |
| 4 | 2 | 0.00 | 0.25 | 0.00 |

## G. References: the same penalties minimised by linear programmes (10 × 40, 40 instances; 12 iterations; ε = 0.1·max|w*|; τ_c = 100)

| k | bp | reweighted ℓ¹ | log-det MM from basis pursuit | log-det MM from the start | their mean support sizes |
|---|---|---|---|---|---|
| 2 | 0.90 | 0.90 | 0.93 | 0.95 | 2.8, 2.8, 2.6, 2.4 |
| 3 | 0.38 | 0.38 | 0.53 | 0.70 | 7.3, 7.3, 6.3, 5.1 |
| 4 | 0.17 | 0.17 | 0.33 | 0.45 | 8.9, 8.9, 8.0, 7.2 |
| 5 | 0.00 | 0.00 | 0.05 | 0.10 | 10.0, 9.9, 9.6, 9.4 |
| 6 | 0.00 | 0.00 | 0.00 | 0.00 | 10.0, 9.8, 10.0, 10.0 |

## H. Which exact fits the selection holds: the dual certificate max_{i∉S}|c_i|/ω_i (below 1: the fit is stable)

| k | instances | at the truth: the fraction certified, white; τ_c = 10, 100, 1000 | the median ratio at the truth, white; τ_c = 10, 100, 1000 | where basis pursuit fails: its solution's support size, and the fraction of its solutions certified for τ_c = 10, 100, 1000 |
|---|---|---|---|---|
| 2 | 40 | 0.38; 0.97; 1.00; 1.00 | 1.2; 0.439; 0.0814; 0.0087 | 10.0; 0.00; 0.00; 0.00 |
| 3 | 40 | 0.03; 0.93; 1.00; 1.00 | 1.72; 0.585; 0.111; 0.0126 | 10.0; 0.24; 0.12; 0.24 |
| 4 | 40 | 0.03; 0.75; 1.00; 1.00 | 1.98; 0.788; 0.185; 0.0217 | 10.0; 0.21; 0.24; 0.30 |
| 5 | 40 | 0.00; 0.50; 1.00; 1.00 | 2.37; 0.997; 0.268; 0.0393 | 10.0; 0.25; 0.23; 0.23 |
| 6 | 40 | 0.00; 0.30; 0.93; 1.00 | 2.79; 1.19; 0.407; 0.0661 | 10.0; 0.12; 0.07; 0.10 |

Numbers are produced by `experiments/exp_colour_selection.py`; nothing is edited by hand.
