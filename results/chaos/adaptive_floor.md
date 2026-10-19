# Proposition D‴.24: the finite-size floor is the leak's bias; an adaptive prior removes it without a mass constraint (`exp_adaptive_floor.py`)

The three-bump instance (σ_crit = 0.04716); the warm clouds of Proposition D‴.12 (width 0.02), m = 15, 60, 240 against M = 510, 4 seeds; midpoint rule dt = 0.25, the plan every 4 steps.

## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)

F: Proposition D‴.9's function-level sensitivity along e₁ − e₃, e₂ − e₃ (central differences of the fixed-prior nearest-atom stationary point); F_bias = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰. The law per unit of 1/m − 1/M, tr(FC̃Fᵀ)/N, from F and as the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j; the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ).

| σ₀/σ_crit | ‖F − F_bias‖/‖F‖ | law from F | sum of the mode terms | bound | law/bound | floor (1/m − 1/M)·law at m = 15, 60, 240 |
|---|---|---|---|---|---|---|
| 0.05 | 1.1e-06 | 1.3559e-02 | 1.3559e-02 | 2.7310e-02 | 0.496 | 8.77e-04, 1.99e-04, 2.99e-05 |
| 0.15 | 1.9e-07 | 6.1662e-02 | 6.1662e-02 | 8.1930e-02 | 0.753 | 3.99e-03, 9.07e-04, 1.36e-04 |
| 0.5 | 1.2e-08 | 2.2913e-01 | 2.2913e-01 | 2.7310e-01 | 0.839 | 1.48e-02, 3.37e-03, 5.05e-04 |

The parameter level: the fixed prior's weight sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e with S^e = ∂_wϑ^e (per mode the fraction h_j/(h_j + σ₀) of the gauge's displacement, the rest σ₀/(h_j + σ₀) being the bias), and Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M, Σ tr(S̃C̃S̃ᵀ), for the fixed prior and for the adaptive one (S^e: the whole displacement in the parameters).

| σ₀/σ_crit | ‖S − H_L(H_L + σ₀I)⁻¹S^e‖/‖S‖ | coupled-distance law, fixed prior | adaptive prior | ratio |
|---|---|---|---|---|
| 0.05 | 1.4e-07 | 3.2274e+01 | 4.6325e+01 | 1.435 |
| 0.15 | 9.2e-08 | 2.1321e+01 | 4.6325e+01 | 2.173 |
| 0.5 | 2.9e-08 | 8.7922e+00 | 4.6325e+01 | 5.269 |

Mode by mode (whitened weighted Gauss–Newton eigenvalues h_j; the share of each term in the law):

| σ₀/σ_crit | mode 1 | mode 2 | mode 3 | mode 4 | mode 5 | mode 6 | mode 7 | mode 8 | mode 9 |
|---|---|---|---|---|---|---|---|---|---|
| 0.05: h_j/σ₀ | 0.13 | 0.83 | 1.32 | 6.42 | 8.95 | 13.77 | 37.58 | 56.91 | 82.78 |
| 0.05: share | 0.000 | 0.109 | 0.313 | 0.261 | 0.220 | 0.091 | 0.005 | 0.000 | 0.000 |
| 0.15: h_j/σ₀ | 0.04 | 0.28 | 0.44 | 2.14 | 2.98 | 4.59 | 12.53 | 18.97 | 27.59 |
| 0.15: share | 0.000 | 0.050 | 0.179 | 0.321 | 0.302 | 0.140 | 0.009 | 0.000 | 0.000 |
| 0.5: h_j/σ₀ | 0.01 | 0.08 | 0.13 | 0.64 | 0.89 | 1.38 | 3.76 | 5.69 | 8.28 |
| 0.5: share | 0.000 | 0.019 | 0.078 | 0.316 | 0.359 | 0.208 | 0.020 | 0.001 | 0.000 |

Admissible weights: an exact representation at the weights w needs |A_k|/w_k < a_max = 4, i.e. w > (0.250, 0.175, 0.200). Per seed and size: the sampled counts, the margin min_k(w_k − |A_k|/a_max), and the exact representation's gradient residual and fit error.

| m | seed 0 | seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| 15 | (6, 5, 4): +0.0667, 1e-16, 2e-16 | (6, 3, 6): +0.0250, 3e-17, 2e-16 | (5, 4, 6): +0.0833, 2e-16, 2e-16 | (8, 4, 3): +0.0000 (none) |
| 60 | (19, 21, 20): +0.0667, 3e-17, 2e-16 | (24, 15, 21): +0.0750, 7e-14, 4e-13 | (19, 18, 23): +0.0667, 2e-17, 9e-17 | (25, 17, 18): +0.1000, 4e-17, 1e-16 |
| 240 | (71, 76, 93): +0.0458, 7e-17, 2e-16 | (77, 76, 87): +0.0708, 2e-16, 3e-16 | (76, 81, 83): +0.0667, 2e-17, 1e-16 | (82, 79, 79): +0.0917, 5e-14, 6e-13 |
| 510 | (158, 161, 191): +0.0598, 1e-14, 4e-14 | (163, 166, 181): +0.0696, 3e-17, 1e-16 | (170, 167, 173): +0.0833, 5e-16, 3e-15 | (168, 184, 158): +0.0794, 3e-16, 3e-16 |

## B. The nearest-atom leak at the matching prior: fixed (η = 0) against adaptive priors (6000 slow units)

‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = 600, 3000, 6000, per seed; W²(μ_m, μ_M) at the end against the W² between the exact representations of the two systems' sampled weights (the parameter-level floor); the class means' whitened weighted distance from the exact representation of the sampled weights at the end; the counts (all conserved). —: no exact representation at those weights (on the boundary of the admissible set).

| η | m | ‖f̂_m − f̂_M‖² at s = 600, 3000, 6000 (per seed) | fit error at 6000 (per seed) | W²(μ_m, μ_M) at the end / exact-representation W² (per seed) | distance from ϑ^e(w) at the end | counts conserved |
|---|---|---|---|---|---|---|
| 0.00000 | 15 | 1.8e-03, 1.7e-03, 1.7e-03; 2.5e-03, 2.4e-03, 2.4e-03; 3.0e-04, 3.0e-04, 3.0e-04; 5.5e-03, 5.5e-03, 5.5e-03 | 4.2e-02; 7.4e-02; 2.9e-02; 1.2e-01 | 4.48e+00 / 5.92e+00; 5.51e+00 / 7.27e+00; 2.22e+00 / 2.32e+00; 7.99e+00 / — | 3.7e-01; 9.5e-01; 2.0e-01; — | yes |
| 0.00000 | 60 | 1.1e-04, 1.1e-04, 1.1e-04; 1.1e-03, 9.9e-04, 9.9e-04; 1.6e-04, 1.6e-04, 1.6e-04; 7.3e-04, 7.1e-04, 7.1e-04 | 1.3e-02; 4.3e-02; 2.1e-02; 3.8e-02 | 1.57e+00 / 1.59e+00; 3.70e+00 / 4.16e+00; 1.80e+00 / 1.86e+00; 3.63e+00 / 3.95e+00 | 1.1e-01; 3.8e-01; 1.4e-01; 3.6e-01 | yes |
| 0.00000 | 240 | 8.2e-05, 8.2e-05, 8.2e-05; 4.9e-06, 4.8e-06, 4.8e-06; 7.3e-05, 7.3e-05, 7.3e-05; 7.5e-05, 7.5e-05, 7.5e-05 | 3.6e-02; 1.4e-02; 1.3e-02; 5.8e-03 | 6.81e-01 / 7.59e-01; 3.52e-01 / 3.55e-01; 8.34e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 3.3e-01; 9.0e-02; 1.1e-01; 4.1e-02 | yes |
| 0.00707 | 15 | 1.2e-04, 1.3e-06, 8.8e-09; 5.4e-04, 2.3e-05, 2.0e-06; 2.3e-05, 1.4e-07, 4.6e-09; 1.3e-03, 2.9e-04, 1.5e-04 | 8.7e-05; 2.2e-03; 9.5e-05; 1.9e-02 | 5.92e+00 / 5.92e+00; 7.15e+00 / 7.27e+00; 2.32e+00 / 2.32e+00; 1.34e+01 / — | 2.3e-03; 6.0e-02; 3.0e-03; — | yes |
| 0.00707 | 60 | 5.3e-06, 8.0e-08, 2.6e-09; 1.1e-04, 2.8e-07, 5.4e-10; 9.2e-06, 1.2e-07, 6.2e-09; 6.2e-05, 1.8e-07, 1.7e-09 | 8.0e-05; 5.7e-05; 1.1e-04; 4.4e-05 | 1.58e+00 / 1.59e+00; 4.16e+00 / 4.16e+00; 1.86e+00 / 1.86e+00; 3.95e+00 / 3.95e+00 | 3.0e-03; 1.8e-03; 3.7e-03; 8.3e-04 | yes |
| 0.00707 | 240 | 1.9e-05, 4.0e-07, 1.9e-08; 3.7e-07, 1.1e-09, 2.5e-12; 8.9e-06, 4.8e-08, 3.3e-09; 6.8e-06, 1.2e-08, 1.3e-09 | 3.7e-04; 7.9e-05; 7.7e-05; 1.9e-05 | 7.58e-01 / 7.59e-01; 3.55e-01 / 3.55e-01; 8.67e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 1.3e-02; 2.8e-03; 2.8e-03; 6.8e-04 | yes |
| ∞ | 15 | 3.4e-05, 7.2e-09, 6.5e-10; 1.8e-04, 1.5e-06, 5.4e-08; 4.4e-06, 2.0e-09, 1.9e-10; 6.3e-04, 1.4e-04, 7.0e-05 | 1.1e-05; 3.9e-04; 1.1e-05; 1.3e-02 | 5.92e+00 / 5.92e+00; 7.26e+00 / 7.27e+00; 2.32e+00 / 2.32e+00; 1.48e+01 / — | 4.2e-04; 1.2e-02; 3.9e-04; — | yes |
| ∞ | 60 | 1.5e-06, 2.1e-09, 2.4e-10; 2.1e-05, 4.5e-09, 6.2e-10; 2.5e-06, 4.8e-09, 6.2e-10; 1.3e-05, 3.5e-09, 4.7e-10 | 2.7e-05; 7.4e-06; 2.9e-05; 1.8e-05 | 1.58e+00 / 1.59e+00; 4.16e+00 / 4.16e+00; 1.86e+00 / 1.86e+00; 3.95e+00 / 3.95e+00 | 1.0e-03; 2.5e-04; 1.0e-03; 6.2e-04 | yes |
| ∞ | 240 | 5.8e-06, 1.6e-08, 2.2e-09; 7.2e-08, 1.5e-12, 7.1e-14; 1.8e-06, 3.0e-09, 4.8e-10; 1.1e-06, 9.9e-10, 1.7e-10 | 1.2e-04; 3.1e-05; 2.4e-05; 4.5e-06 | 7.59e-01 / 7.59e-01; 3.55e-01 / 3.55e-01; 8.67e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 4.5e-03; 1.1e-03; 8.9e-04; 1.7e-04 | yes |

## C. The design table at the misplaced prior (centres +0.15), 3000 slow units

The fit error (the bias; 0 at an exact representation) and ‖f̂_m − f̂_M‖² (the floor) at the end, means over the seeds (min–max); the OT leak toward the adaptive prior with the settled gate of Theorem D‴.23.

| leak | prior | m | fit error at the end | ‖f̂_m − f̂_M‖² at the end |
|---|---|---|---|---|
| nearest atom | fixed | 15 | 6.76e-02 (3.0e-02–1.2e-01) | 2.50e-03 (3.0e-04–5.6e-03) |
| nearest atom | fixed | 60 | 3.12e-02 (1.6e-02–4.5e-02) | 5.10e-04 (1.0e-04–1.0e-03) |
| nearest atom | fixed | 240 | 2.01e-02 (1.2e-02–3.7e-02) | 5.89e-05 (4.8e-06–8.2e-05) |
| nearest atom | adaptive (η = σ₀) | 15 | 8.93e-03 (5.8e-04–2.7e-02) | 7.96e-05 (1.4e-07–2.9e-04) |
| nearest atom | adaptive (η = σ₀) | 60 | 5.53e-04 (2.8e-04–7.3e-04) | 1.64e-07 (8.0e-08–2.8e-07) |
| nearest atom | adaptive (η = σ₀) | 240 | 5.97e-04 (1.4e-04–1.6e-03) | 1.14e-07 (1.1e-09–3.9e-07) |
| OT | fixed | 15 | 9.78e-03 (9.8e-03–9.8e-03) | 5.08e-22 (2.3e-22–8.0e-22) |
| OT | fixed | 60 | 9.78e-03 (9.8e-03–9.8e-03) | 2.47e-22 (2.4e-23–6.2e-22) |
| OT | fixed | 240 | 9.78e-03 (9.8e-03–9.8e-03) | 1.16e-21 (4.1e-23–2.4e-21) |
| OT | adaptive (η = σ₀, settled gate) | 15 | 1.59e-04 (1.5e-04–1.6e-04) | 1.05e-11 (2.6e-12–2.8e-11) |
| OT | adaptive (η = σ₀, settled gate) | 60 | 1.62e-04 (1.6e-04–1.6e-04) | 1.37e-11 (9.6e-15–3.2e-11) |
| OT | adaptive (η = σ₀, settled gate) | 240 | 1.65e-04 (1.6e-04–1.7e-04) | 1.44e-11 (2.3e-13–4.2e-11) |

## D. The real closed loop (k_e = 25, dt = 0.001, 3000 time units), matching prior, nearest-atom leak, seed 0's sampled weights

The rms tracking error over the five periods ending at t = 31, 1000, 2000, 3000, against (1/k_e)‖f − f̂‖_ρ of the averaged loop from the same start at s = t/k_e (the slaved error of Theorem J at leading order); from the exact representation of w with the adaptive prior at its atoms the tracking error stays at the integrator's level.

| m | counts | start, prior | e_rms at those times | (1/k_e)‖f − f̂‖_ρ averaged at s = t/k_e | ratio at the end |
|---|---|---|---|---|---|
| 15 | (6, 5, 4) | Θ*(w), fixed prior | 1.173e-03, 1.180e-03, 1.180e-03, 1.182e-03 | 1.181e-03, 1.181e-03, 1.181e-03, 1.181e-03 | 1.001 |
| 15 | (6, 5, 4) | Θ*(w), adaptive prior (η = ∞) | 1.157e-03, 6.102e-04, 4.062e-04, 3.273e-04 | 1.151e-03, 6.056e-04, 4.093e-04, 3.317e-04 | 0.987 |
| 15 | (6, 5, 4) | ϑ^e(w), adaptive prior at ϑ^e(w) | 9.0e-12, 3.6e-12, 2.1e-12, 1.4e-12 | 6.9e-18, 6.8e-18, 6.8e-18, 6.8e-18 | — |
| 60 | (19, 21, 20) | Θ*(w), fixed prior | 3.427e-04, 3.429e-04, 3.433e-04, 3.441e-04 | 3.456e-04, 3.456e-04, 3.456e-04, 3.456e-04 | 0.996 |
| 60 | (19, 21, 20) | Θ*(w), adaptive prior (η = ∞) | 3.397e-04, 2.294e-04, 1.794e-04, 1.522e-04 | 3.401e-04, 2.312e-04, 1.812e-04, 1.534e-04 | 0.992 |
| 60 | (19, 21, 20) | ϑ^e(w), adaptive prior at ϑ^e(w) | 9.0e-12, 3.7e-12, 2.2e-12, 1.5e-12 | 8.7e-18, 8.5e-18, 8.5e-18, 8.6e-18 | — |

Numbers are produced by `experiments/exp_adaptive_floor.py`; nothing is edited by hand.
