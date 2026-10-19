# Remark D‴.16: dither and the mass constraint — the dithered averaged slow loop under the nearest-atom and the OT leak (`exp_dither_ot.py`)

The three-bump instance with the matching prior (three atoms at the exact representation, masses 1/3), σ₀ = 0.15σ_crit = 0.00707, m = 15 particles started at the collapsed configuration with counts (6, 5, 4) (jitter 0.02), the Langevin dynamics dϑ_i = −G∇J_m ds + √(2τ₀)G^{1/2}dW_i in slow time (Euler–Maruyama, ds = 0.05; the leak target recomputed every 10 steps), cluster masses by nearest atom recorded every 10 slow units; horizons 50000, 200000, 50000 slow units at τ₀ = 0.005, 0.01, 0.02; seeds 0, 1. The Gibbs prediction for the nearest-atom leak is the Laplace approximation over the count vectors, P(n) ∝ multinomial(n)·exp(−m𝒢*_{n/m}/τ₀)·Π_k det(σ₀I − R^G_k)^{−(n_k−1)/2}·det(H_n)^{−1/2} (count vectors whose replicated configuration is a saddle carry no Laplace mode and are omitted); 'hop records' counts the records at which some particle's nearest atom changed.

## τ₀ = 0.005 (horizon 50000 slow units; the OT leak's barrier per misplaced particle σ₀Δ_ν²/2 = 0.139 = 28τ₀)

Gibbs prediction for the nearest-atom leak at leading order (entropy × energy, multinomial(n)·exp(−m𝒢*_n/τ₀)), P(n): (5, 5, 5): 0.862, (5, 4, 6): 0.075, (5, 6, 4): 0.043, (6, 4, 5): 0.013, (6, 5, 4): 0.004; with the Laplace determinant factors: (5, 5, 5): 0.720, (5, 6, 4): 0.243, (5, 4, 6): 0.035 — the determinant factors over-weight the count vectors whose over-weighted class has a nearly flat separation block (a soft mode, where the quadratic approximation fails; Remark D‴.11), and are not used in the reading.

| leak | seed | hop records / records | time fraction by counts (top four) | mean |w − q|₁ | final counts |
|---|---|---|---|---|---|
| nearest | 0 | 0/4999 | (6, 5, 4): 1.000 | 0.133 | (np.int64(6), np.int64(5), np.int64(4)) |
| nearest | 1 | 1/4999 | (6, 5, 4): 0.976, (5, 6, 4): 0.024 | 0.133 | (np.int64(5), np.int64(6), np.int64(4)) |
| ot | 0 | 1/4999 | (5, 5, 5): 0.996, (6, 5, 4): 0.004 | 0.001 | (np.int64(5), np.int64(5), np.int64(5)) |
| ot | 1 | 1/4999 | (5, 5, 5): 0.997, (6, 5, 4): 0.003 | 0.000 | (np.int64(5), np.int64(5), np.int64(5)) |

## τ₀ = 0.01 (horizon 200000 slow units; the OT leak's barrier per misplaced particle σ₀Δ_ν²/2 = 0.139 = 14τ₀)

Gibbs prediction for the nearest-atom leak at leading order (entropy × energy, multinomial(n)·exp(−m𝒢*_n/τ₀)), P(n): (5, 5, 5): 0.576, (5, 4, 6): 0.156, (5, 6, 4): 0.118, (6, 4, 5): 0.064, (6, 5, 4): 0.037; with the Laplace determinant factors: (5, 6, 4): 0.527, (5, 5, 5): 0.383, (5, 4, 6): 0.058 — the determinant factors over-weight the count vectors whose over-weighted class has a nearly flat separation block (a soft mode, where the quadratic approximation fails; Remark D‴.11), and are not used in the reading.

| leak | seed | hop records / records | time fraction by counts (top four) | mean |w − q|₁ | final counts |
|---|---|---|---|---|---|
| nearest | 0 | 111/19999 | (5, 5, 5): 0.562, (6, 4, 5): 0.208, (5, 4, 6): 0.097, (6, 5, 4): 0.069 | 0.059 | (np.int64(5), np.int64(4), np.int64(6)) |
| nearest | 1 | 93/19999 | (5, 5, 5): 0.511, (6, 4, 5): 0.203, (6, 5, 4): 0.127, (5, 4, 6): 0.088 | 0.065 | (np.int64(5), np.int64(5), np.int64(5)) |
| ot | 0 | 131/19999 | (5, 5, 5): 0.991, (5, 4, 6): 0.004, (5, 6, 4): 0.002, (6, 4, 5): 0.001 | 0.001 | (np.int64(5), np.int64(5), np.int64(5)) |
| ot | 1 | 115/19999 | (5, 5, 5): 0.994, (5, 4, 6): 0.003, (5, 6, 4): 0.001, (6, 5, 4): 0.001 | 0.001 | (np.int64(5), np.int64(5), np.int64(5)) |

## τ₀ = 0.02 (horizon 50000 slow units; the OT leak's barrier per misplaced particle σ₀Δ_ν²/2 = 0.139 = 7τ₀)

Gibbs prediction for the nearest-atom leak at leading order (entropy × energy, multinomial(n)·exp(−m𝒢*_n/τ₀)), P(n): (5, 5, 5): 0.332, (5, 4, 6): 0.157, (5, 6, 4): 0.137, (6, 4, 5): 0.101, (6, 5, 4): 0.076; with the Laplace determinant factors: (5, 6, 4): 0.520, (5, 5, 5): 0.187, (5, 7, 3): 0.098 — the determinant factors over-weight the count vectors whose over-weighted class has a nearly flat separation block (a soft mode, where the quadratic approximation fails; Remark D‴.11), and are not used in the reading.

| leak | seed | hop records / records | time fraction by counts (top four) | mean |w − q|₁ | final counts |
|---|---|---|---|---|---|
| nearest | 0 | 281/4999 | (5, 5, 5): 0.211, (6, 4, 5): 0.163, (5, 4, 6): 0.160, (5, 6, 4): 0.112 | 0.125 | (np.int64(5), np.int64(4), np.int64(6)) |
| nearest | 1 | 275/4999 | (5, 5, 5): 0.287, (6, 4, 5): 0.165, (6, 5, 4): 0.154, (5, 4, 6): 0.152 | 0.105 | (np.int64(5), np.int64(5), np.int64(5)) |
| ot | 0 | 364/4999 | (5, 5, 5): 0.861, (5, 4, 6): 0.057, (5, 6, 4): 0.026, (6, 5, 4): 0.026 | 0.019 | (np.int64(5), np.int64(5), np.int64(5)) |
| ot | 1 | 344/4999 | (5, 5, 5): 0.894, (5, 4, 6): 0.037, (6, 4, 5): 0.025, (5, 6, 4): 0.020 | 0.014 | (np.int64(5), np.int64(5), np.int64(5)) |

## Reading

Nearest-atom leak: hop rates 1.00e-05 per slow unit at τ₀ = 0.005, 5.10e-04 per slow unit at τ₀ = 0.01, 5.56e-03 per slow unit at τ₀ = 0.02 — Arrhenius barriers ΔV = 0.048 from τ₀ = 0.01, 0.02, 0.039 from τ₀ = 0.005, 0.01 against the leak's ridge height σ₀Δ_ν²/8 = 0.035. Time fraction at the balanced counts (5, 5, 5), nearest-atom leak against the OT leak against the leading-order Gibbs prediction: τ₀ = 0.005: 0.000 / 0.996 / 0.862; τ₀ = 0.01: 0.537 / 0.992 / 0.576; τ₀ = 0.02: 0.249 / 0.878 / 0.332. At the lowest temperature the nearest-atom leak's masses do not move within the horizon (the hop time exceeds it) while the OT leak holds the balanced counts; at the highest the nearest-atom leak samples its Gibbs law over the count vectors by rare hops and the OT leak still holds the balanced counts most of the time, the excursions being the transits that follow a swap of assignments: dither does on the Arrhenius time scale, and only in distribution, what the mass constraint does deterministically in ln 2/σ₀.

Numbers are produced by `experiments/exp_dither_ot.py`; nothing is edited by hand.
