# Proposition D‴.28: the adaptive prior's floor removal is a gauge property (`exp_gauge_floor.py`)

The three-bump instance with a fourth bump ε·exp(−(x − 0.8)²/(2·0.3²)) that K₀ = 3 kernels cannot represent; masses q = 1/3, σ₀ = 0.15σ_crit = 0.00707; ϑ^b(q) the best weighted approximant (a critical point of the unregularised loss, continued in ε).

## A. The approximant, the gauge and the sensitivities (no dynamics)

Per ε: the approximant's relative fit error and gradient residual; the smallest eigenvalue of the whitened unregularised Hessian Ã = H_L − R̃ and ‖R̃‖/‖H_L‖; the gauge check at six weight vectors of W_g (the largest distance of Newton's critical point from G_wϑ^b(q), and of its f̂ from f̂_b(q)); the adaptive prior's function-level sensitivity relative to the fixed prior's (central differences); the fixed matching prior's sensitivity against σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and against D‴.24's form with H_L; the fixed prior's floor law per unit of 1/m − 1/M.

| ε | fit error | residual | min eig Ã | ‖R̃‖/‖H_L‖ | gauge: distance, Δf̂ | ‖F_ad‖/‖F_fixed‖ | F_fixed vs σ₀Ĵ(Ã + σ₀I)⁻¹g̃ | vs H_L form | law (fixed) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1.104e-16 | 2.3e-17 | 3.059e-04 | 1.95e-11 | 0.0e+00, 2.2e-16 | 0.0e+00 | 1.9e-09 | 1.9e-09 | 6.1662e-02 |
| 0.05 | 3.591e-03 | 9.0e-17 | 2.873e-04 | 5.29e-04 | 0.0e+00, 7.8e-16 | 3.3e-12 | 1.9e-09 | 7.0e-04 | 6.2119e-02 |
| 0.1 | 6.841e-03 | 1.7e-16 | 2.729e-04 | 9.32e-04 | 0.0e+00, 6.7e-16 | 1.4e-11 | 2.2e-09 | 1.3e-03 | 6.2718e-02 |
| 0.2 | 1.239e-02 | 1.6e-15 | 2.574e-04 | 1.43e-03 | 0.0e+00, 2.2e-16 | 6.2e-12 | 2.0e-09 | 2.3e-03 | 6.4655e-02 |
| 0.3 | 1.695e-02 | 4.6e-16 | 2.426e-04 | 1.72e-03 | 0.0e+00, 1.1e-16 | 2.3e-11 | 2.1e-09 | 3.0e-03 | 6.8178e-02 |
| 0.5 | 2.528e-02 | 1.4e-16 | 4.602e-05 | 2.83e-03 | 0.0e+00, 5.6e-16 | 0.0e+00 | 2.0e-09 | 5.6e-03 | 8.5230e-02 |

The residual curvature at ϑ^b(q): the largest entry of its amplitude row and column (zero at a critical point), and λ⁺_k = λ_max(G^{1/2}R_kG^{1/2}) per class; W_g's thresholds q_k|a^b_k|/a_max.

| ε | amplitude row | largest entry of R | λ⁺_k | W_g: w > |
|---|---|---|---|---|
| 0 | 6.1e-18 | 1.08e-16 | 1.71e-19, 4.89e-20, -8.52e-19 | (0.2500, 0.1750, 0.2000) |
| 0.05 | 2.1e-14 | 2.64e-04 | 2.49e-05, 6.34e-05, 6.44e-05 | (0.2501, 0.1763, 0.2006) |
| 0.1 | 3.9e-14 | 4.80e-04 | 4.47e-05, 1.16e-04, 1.08e-04 | (0.2502, 0.1782, 0.2016) |
| 0.2 | 6.5e-14 | 8.26e-04 | 7.41e-05, 2.00e-04, 1.37e-04 | (0.2503, 0.1852, 0.2045) |
| 0.3 | 7.4e-14 | 1.17e-03 | 1.00e-04, 2.86e-04, 1.08e-04 | (0.2505, 0.1994, 0.2097) |
| 0.5 | 6.5e-14 | 2.91e-03 | 2.10e-04, 7.20e-04, 4.98e-06 | (0.2497, 0.2680, 0.2335) |

The adaptive prior's splitting margins at G_wϑ^b(q) over every count vector of m in W_g (with no empty class): the largest difference from the prediction σ₀ − (q_k/w_k)λ⁺_k, the probability of a negative margin (some class can split), and the smallest margin.

| ε | m | count vectors | largest |computed − predicted| | P(negative margin) | smallest margin |
|---|---|---|---|---|---|
| 0.05 | 15 | 15 | 2.2e-15 | 0.0000 | +6.97e-03 |
| 0.05 | 60 | 231 | 2.4e-15 | 0.0000 | +6.96e-03 |
| 0.1 | 15 | 15 | 6.1e-15 | 0.0000 | +6.88e-03 |
| 0.1 | 60 | 231 | 6.7e-15 | 0.0000 | +6.86e-03 |
| 0.2 | 15 | 15 | 5.3e-15 | 0.0000 | +6.74e-03 |
| 0.2 | 60 | 210 | 5.3e-15 | 0.0000 | +6.74e-03 |
| 0.3 | 15 | 15 | 9.8e-15 | 0.0000 | +6.60e-03 |
| 0.3 | 60 | 210 | 9.8e-15 | 0.0000 | +6.60e-03 |
| 0.5 | 15 | 6 | 1.4e-14 | 0.0000 | +6.35e-03 |
| 0.5 | 60 | 105 | 2.1e-14 | 0.0000 | +6.23e-03 |

## B. The nearest-atom leak from warm clouds around ϑ^b(q): fixed (η = 0) against adaptive priors (6000 slow units)

Per ε, η and m: ‖f̂_m − f̂_M‖² at s = 600, 3000, 6000 (min–max over the seeds), the fit error at the end (min–max) against the approximant's, the class means' whitened weighted distance from G_wϑ^b(q) at the end (max), the largest class diameter at the end (max), and the counts conserved.

| ε | η | m | ‖f̂_m − f̂_M‖² at those s | fit error at the end (approximant) | distance from G_wϑ^b | class diameter | counts conserved |
|---|---|---|---|---|---|---|---|
| 0.1 | 0.00000 | 15 | 3.1e-04–5.8e-03, 3.1e-04–5.9e-03, 3.1e-04–5.9e-03 | 3.0137e-02–1.2184e-01 (6.8407e-03) | 1.1e+00 | 1.5e+00 | yes |
| 0.1 | 0.00000 | 60 | 1.1e-04–1.1e-03, 1.1e-04–1.0e-03, 1.1e-04–1.0e-03 | 1.5002e-02–4.3962e-02 (6.8407e-03) | 4.0e-01 | 9.4e-01 | yes |
| 0.1 | 0.00000 | 240 | 5.0e-06–8.2e-05, 5.0e-06–8.2e-05, 5.0e-06–8.2e-05 | 8.9652e-03–3.6328e-02 (6.8407e-03) | 3.4e-01 | 1.2e-07 | yes |
| 0.1 | 0.00707 | 15 | 3.0e-05–1.4e-03, 2.5e-07–3.4e-04, 1.8e-09–1.8e-04 | 6.8411e-03–2.1723e-02 (6.8407e-03) | 9.1e-02 | 1.2e-07 | yes |
| 0.1 | 0.00707 | 60 | 6.9e-06–1.2e-04, 9.9e-08–3.6e-07, 1.0e-09–1.4e-08 | 6.8412e-03–6.8418e-03 (6.8407e-03) | 4.0e-03 | 0.0e+00 | yes |
| 0.1 | 0.00707 | 240 | 4.6e-07–1.9e-05, 1.9e-09–4.1e-07, 5.6e-11–1.9e-08 | 6.8408e-03–6.8496e-03 (6.8407e-03) | 1.3e-02 | 1.2e-07 | yes |
| 0.1 | ∞ | 15 | 6.7e-06–6.7e-04, 1.9e-09–1.7e-04, 8.7e-11–9.6e-05 | 6.8407e-03–1.6494e-02 (6.8407e-03) | 1.6e-02 | 1.2e-07 | yes |
| 0.1 | ∞ | 60 | 2.1e-06–2.4e-05, 6.5e-10–2.6e-08, 3.2e-11–3.8e-09 | 6.8408e-03–6.8410e-03 (6.8407e-03) | 2.3e-03 | 1.2e-07 | yes |
| 0.1 | ∞ | 240 | 9.4e-08–5.9e-06, 4.7e-11–1.6e-08, 1.0e-11–2.6e-09 | 6.8407e-03–6.8419e-03 (6.8407e-03) | 4.9e-03 | 1.2e-07 | yes |
| 0.3 | 0.00000 | 15 | 3.7e-04–7.3e-03, 3.7e-04–7.4e-03, 3.7e-04–7.4e-03 | 3.5470e-02–1.3426e-01 (1.6953e-02) | 4.2e+00 | 1.6e+00 | yes |
| 0.3 | 0.00000 | 60 | 1.2e-04–1.2e-03, 1.2e-04–1.1e-03, 1.2e-04–1.1e-03 | 2.1689e-02–4.8465e-02 (1.6953e-02) | 6.1e-01 | 9.6e-01 | yes |
| 0.3 | 0.00000 | 240 | 6.0e-06–8.8e-05, 6.0e-06–8.8e-05, 6.0e-06–8.8e-05 | 1.7943e-02–3.9370e-02 (1.6953e-02) | 3.5e-01 | 1.2e-07 | yes |
| 0.3 | 0.00707 | 15 | 7.3e-05–2.4e-03, 3.9e-06–6.8e-04, 1.1e-07–4.8e-04 | 1.6956e-02–3.7051e-02 (1.6953e-02) | 3.1e+00 | 1.2e-07 | yes |
| 0.3 | 0.00707 | 60 | 1.6e-05–2.0e-04, 1.0e-06–1.1e-05, 1.0e-07–2.6e-06 | 1.6957e-02–1.7157e-02 (1.6953e-02) | 1.1e-01 | 1.2e-07 | yes |
| 0.3 | 0.00707 | 240 | 1.1e-06–1.9e-05, 5.8e-08–5.0e-07, 1.1e-08–1.0e-07 | 1.6954e-02–1.6957e-02 (1.6953e-02) | 9.1e-03 | 1.2e-07 | yes |
| 0.3 | ∞ | 15 | 3.2e-05–1.2e-03, 8.6e-08–4.7e-04, 1.2e-10–3.5e-04 | 1.6954e-02–3.3012e-02 (1.6953e-02) | 2.8e+00 | 1.2e-07 | yes |
| 0.3 | ∞ | 60 | 7.8e-06–6.3e-05, 9.9e-08–2.5e-06, 1.6e-08–4.9e-07 | 1.6954e-02–1.6991e-02 (1.6953e-02) | 4.9e-02 | 1.2e-07 | yes |
| 0.3 | ∞ | 240 | 4.2e-07–5.9e-06, 1.1e-08–9.6e-08, 2.4e-09–2.2e-08 | 1.6953e-02–1.6954e-02 (1.6953e-02) | 4.0e-03 | 8.4e-08 | yes |

The samples within 0.02 of the boundary of W_g or outside it (the margin min_k(w_k − q_k|a^b_k|/a_max); negative: outside, no gauge transform): the fit error at the end under the fixed prior and the adaptive priors (η = σ₀, ∞) against the approximant's, and the adaptive class means' distance from G_wϑ^b(q).

| ε | m | seed | w | margin | fit error: fixed, η = σ₀, η = ∞ (approximant) | distance from G_wϑ^b: η = σ₀, ∞ |
|---|---|---|---|---|---|---|
| 0.1 | 15 | 3 | (0.533, 0.267, 0.200) | -0.0016 | 1.2184e-01, 2.1723e-02, 1.6494e-02 (6.8407e-03) | —, — |
| 0.3 | 15 | 1 | (0.400, 0.200, 0.400) | +0.0006 | 7.9743e-02, 2.2916e-02, 2.0038e-02 (1.6953e-02) | 3.07e+00, 2.79e+00 |
| 0.3 | 15 | 3 | (0.533, 0.267, 0.200) | -0.0097 | 1.3426e-01, 3.7051e-02, 3.3012e-02 (1.6953e-02) | —, — |

At the other samples (margin ≥ 0.02) the adaptive priors' fit error at the end exceeds the approximant's by at most 1.1e-01 (m = 15, median 1.2e-04), 1.2e-02 (m = 60, median 1.2e-04), 1.3e-03 (m = 240, median 1.6e-05) (relative): still converging toward the gauge transform, slowest at small m, where the samples lie nearest the boundary.

Numbers are produced by `experiments/exp_gauge_floor.py`; nothing is edited by hand.
