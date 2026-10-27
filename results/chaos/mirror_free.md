# Proposition D‴.45: the free mass–amplitude law is a mirror flow (`exp_mirror_free.py`)

Free masses (λ ≡ 0): dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, a = a_max tanh(α/a_max), w = β⊙a/m; m = 12, G = 1, γ = 0.1, a_max = 1 (c = γa_max²/(4G) = 0.025) unless said otherwise; Q_i = β_i − c cosh(2α_i/a_max); on a dictionary the mirror flow dη = Φᵀ(r + ξ)dt/N of Ψ_q, ψ′ = η = ∫dα/(Gβa′) (`mirror_free`).

## A. Identities

η's closed form against quadrature (q = 100c, 3, 0.5, c, c/2, 0, −c/2 and the branch q = −2c, 17 and 5 points): largest relative difference 1.53e-14; (ψ*)′ = w against finite differences 7.67e-11; dw/dη = mD 4.28e-10.

Q along the loop (5 × 12, 2-sparse truths, β(0) = 0.3, α(0) ~ N(0, 0.05²)), largest |Q(t) − Q(0)| over 60 times; the conserved law's shift dQ_i/dt = −γβ_iλ at the end state (relative):

| instance | noiseless, t ≤ 10⁴ | a smooth noise (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes, amplitude 0.3), t ≤ 3·10⁴ | the conserved law's shift |
|---|---|---|---|
| 0 | 2.04e-12 | 2.89e-10 | 4.35e-14 |
| 1 | 1.40e-12 | 2.50e-10 | 9.09e-15 |
| 2 | 1.41e-11 | 3.61e-10 | 4.85e-14 |

Q's drift along Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update; at the lazy projection) against the Euler step's closed form (2|φ_i|²/(mN))γβ_i(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) per unit s: largest relative difference 9.23e-16, 9.11e-16, 1.44e-15.

## B. The implicit bias: the start's masses set the regime (5 × 12, 20 instances per k)

From α = 0 and equal masses β = s·m/d: the Bregman projection (Newton on the dual) — ‖w‖₁ over basis pursuit's, the atoms used (|w| > 10⁻³), the fraction recovering the truth (10⁻³·max|w*|), the median distances to basis pursuit and (relative) to the minimum-norm fit Φ⁺y, the largest |a|/a_max and β; for s ≥ 0.01 the noiseless loop's end at t = 1e+07 against the projection (its largest residual and distance).

| k | s | ‖w‖₁/‖w_bp‖₁ | used | the truth | median ‖w − w_bp‖∞ | median ‖w − Φ⁺y‖∞/‖Φ⁺y‖∞ | max |a|/a_max, max β | the loop's end: residual, distance to the projection |
|---|---|---|---|---|---|---|---|---|
| 2 | basis pursuit | 1 | | 0.80 | | | | |
| 2 | 10 | 1.3715 | 11.8 | 0.00 | 0.14 | 0.00499 | 0.241, 10 | 1.14e-16, 2.03e-12 |
| 2 | 1 | 1.5248 | 11.8 | 0.00 | 0.159 | 0.197 | 0.979, 2.13 | 8.44e-17, 1.90e-12 |
| 2 | 0.1 | 1.1662 | 11.4 | 0.00 | 0.0973 | 0.408 | 0.993, 3.76 | 5.69e-15, 8.51e-12 |
| 2 | 0.01 | 1.0274 | 8.1 | 0.00 | 0.0235 | 0.814 | 0.994, 3.86 | 4.79e-15, 4.74e-11 |
| 2 | 0.001 | 1.0022 | 3.6 | 0.15 | 9.04e-4 | 1.03 | 0.993, 3.58 | — |
| 2 | 0.0001 | 1.0004 | 3.0 | 0.75 | 1.38e-5 | 1.16 | 0.993, 3.59 | — |
| 3 | basis pursuit | 1 | | 0.30 | | | | |
| 3 | 10 | 1.2770 | 11.9 | 0.00 | 0.111 | 0.00548 | 0.331, 10 | 1.45e-16, 2.73e-12 |
| 3 | 1 | 1.4152 | 11.9 | 0.00 | 0.13 | 0.187 | 0.989, 3.16 | 8.36e-17, 3.82e-12 |
| 3 | 0.1 | 1.1022 | 11.6 | 0.00 | 0.0659 | 0.322 | 0.994, 4.18 | 8.78e-17, 5.37e-12 |
| 3 | 0.01 | 1.0175 | 8.0 | 0.00 | 0.0286 | 0.587 | 0.994, 3.92 | 7.34e-17, 4.10e-11 |
| 3 | 0.001 | 1.0019 | 5.6 | 0.05 | 0.0041 | 0.727 | 0.993, 3.56 | — |
| 3 | 0.0001 | 1.0001 | 4.5 | 0.30 | 1.30e-4 | 0.747 | 0.993, 3.67 | — |

Small random amplitudes (α(0) ~ N(0, 0.05²)) at tiny masses: the particles with q < −c sit on a branch (the sign of w frozen); the projection's atoms used, those on a branch, the fraction of instances whose branches kept their start's signs, the distance to basis pursuit:

| k | s | branches (mean) | used | used on a branch | signs kept | median ‖w − w_bp‖∞ |
|---|---|---|---|---|---|---|
| 2 | 0.01 | 0.0 | 7.8 | 0.00 | 1.00 | 0.0268 |
| 2 | 0.0001 | 4.2 | 3.2 | 1.15 | 1.00 | 6.58e-4 |
| 3 | 0.01 | 0.0 | 8.1 | 0.00 | 1.00 | 0.0236 |
| 3 | 0.0001 | 4.2 | 4.8 | 1.45 | 1.00 | 0.17 |

Unequal masses β⁰_i = s·r_i (r_i log-uniform in [1/4, 4]), α = 0: the projection against weighted basis pursuit (min Σ|w_i|/√r_i) as s → 0 — and against the unweighted one — and against the weighted minimum norm (min Σw_i²/r_i²) as s → ∞ (relative), medians over the instances:

| k | s = 0.001: weighted, unweighted ℓ¹ | s = 1e-05: weighted, unweighted ℓ¹ | s = 1e-07: weighted, unweighted ℓ¹ | s = 100: weighted ℓ² | s = 10000: weighted ℓ² |
|---|---|---|---|---|---|
| 2 | 0.00622, 0.0203 | 2.33e-6, 2.70e-5 | 2.09e-9, 0.0057 | 2.44e-5 | 2.44e-9 |
| 3 | 0.0128, 0.0833 | 3.04e-5, 0.084 | 2.06e-8, 0.091 | 2.52e-5 | 2.52e-9 |

## C. Noise

(1) Continuous time: the loop from α = 0, β = 0.1·m/d driven by a smooth noise (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes, amplitude 0.3) — |Q − Q(0)|, the residual, and the distance of the state's projection to the noiseless loop's limit:

| instance | t | max |Q − Q(0)| | |Φw − y| | the projection's distance |
|---|---|---|---|---|
| 0 | 1000 | 2.52e-12 | 0.826 | 2.73e-13 |
| 0 | 10000 | 2.75e-11 | 0.475 | 4.63e-12 |
| 0 | 30000 | 8.21e-11 | 1.246 | 2.98e-11 |
| 1 | 1000 | 2.55e-12 | 0.276 | 6.66e-14 |
| 1 | 10000 | 2.58e-11 | 0.497 | 5.96e-12 |
| 1 | 30000 | 7.10e-11 | 0.274 | 8.43e-12 |
| 2 | 1000 | 3.38e-12 | 1.053 | 1.09e-12 |
| 2 | 10000 | 3.21e-11 | 1.389 | 3.65e-12 |
| 2 | 30000 | 9.39e-11 | 1.931 | 5.68e-12 |
| 3 | 1000 | 1.81e-12 | 0.455 | 2.56e-13 |
| 3 | 10000 | 2.45e-11 | 0.448 | 2.35e-12 |
| 3 | 30000 | 8.02e-11 | 1.068 | 1.10e-11 |

(2) The Euler step on (α, log β) with label noise (σ = 0.3, 4·10⁵ draws at random states): the mean change of Q against the closed form (h²E[g²]γβ/2)(γa² − Gβa′(1 + a²/a_max²)) — largest relative difference (largest deviation in standard errors): h = 0.5: 0.0313 (14.7); h = 0.5: 0.0332 (16.6); h = 0.5: 0.0376 (24.8); h = 0.05: 0.00191 (2.1); h = 0.05: 0.00409 (2.5); h = 0.05: 0.0047 (3.1).

(3) The selection in the mirror coordinates: Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update) from the lazy projection (s = 1) — ‖w‖₁, the distance to basis pursuit, the atoms used, Q's mean and range:

| instance (basis pursuit = truth) | s | ‖w‖₁ | ‖w − w_bp‖∞ | used | Q: mean, min, max |
|---|---|---|---|---|---|
| 0 (yes) | 0e+00 | 0.7060 | 0.196 | 12 | 0.975, 0.975, 0.975 |
| 0 (yes) | 1e+05 | 0.6146 | 0.154 | 12 | -0.001371, -0.003571, -0.000425 |
| 0 (yes) | 1e+06 | 0.4563 | 1.58e-6 | 2 | -0.02056, -0.02466, -0.0001863 |
| 0 (yes) | 1e+07 | 0.4563 | 4.80e-8 | 2 | -0.02085, -0.02497, -0.0001863 |
| 0 (yes) | 1e+08 | 0.4563 | 6.54e-9 | 2 | -0.02087, -0.025, -0.0001863 |
| 1 (yes) | 0e+00 | 0.5690 | 0.123 | 12 | 0.975, 0.975, 0.975 |
| 1 (yes) | 1e+05 | 0.5053 | 0.0945 | 10 | -0.004291, -0.0195, -0.0004748 |
| 1 (yes) | 1e+06 | 0.3863 | 1.17e-6 | 2 | -0.02058, -0.02468, -0.0002576 |
| 1 (yes) | 1e+07 | 0.3863 | 4.27e-8 | 2 | -0.02085, -0.02497, -0.0002576 |
| 1 (yes) | 1e+08 | 0.3863 | 7.96e-9 | 2 | -0.02088, -0.025, -0.0002576 |
| 2 (yes) | 0e+00 | 0.7899 | 0.178 | 12 | 0.975, 0.975, 0.975 |
| 2 (yes) | 1e+05 | 0.6570 | 0.135 | 12 | -0.001433, -0.003571, -0.0003782 |
| 2 (yes) | 1e+06 | 0.4894 | 1.01e-6 | 2 | -0.02056, -0.02466, -0.0001892 |
| 2 (yes) | 1e+07 | 0.4894 | 3.55e-8 | 2 | -0.02084, -0.02497, -0.0001892 |
| 2 (yes) | 1e+08 | 0.4894 | 4.54e-9 | 2 | -0.02087, -0.025, -0.0001892 |
| 3 (yes) | 0e+00 | 0.7513 | 0.169 | 12 | 0.975, 0.975, 0.975 |
| 3 (yes) | 1e+05 | 0.6050 | 0.158 | 11 | -0.00311, -0.01572, -0.0003515 |
| 3 (yes) | 1e+06 | 0.4857 | 1.41e-6 | 2 | -0.02055, -0.02467, -0.0002009 |
| 3 (yes) | 1e+07 | 0.4857 | 4.65e-8 | 2 | -0.02084, -0.02497, -0.0002009 |
| 3 (yes) | 1e+08 | 0.4857 | 6.88e-9 | 2 | -0.02087, -0.025, -0.0002009 |
| 4 (yes) | 0e+00 | 0.7598 | 0.169 | 12 | 0.975, 0.975, 0.975 |
| 4 (yes) | 1e+05 | 0.5881 | 0.136 | 11 | -0.003355, -0.02053, -0.000337 |
| 4 (yes) | 1e+06 | 0.4273 | 9.99e-7 | 2 | -0.02058, -0.02469, -0.0001962 |
| 4 (yes) | 1e+07 | 0.4273 | 2.26e-8 | 2 | -0.02085, -0.02497, -0.0001962 |
| 4 (yes) | 1e+08 | 0.4273 | 2.67e-9 | 2 | -0.02087, -0.025, -0.0001962 |
| 5 (yes) | 0e+00 | 0.7133 | 0.187 | 12 | 0.975, 0.975, 0.975 |
| 5 (yes) | 1e+05 | 0.5995 | 0.123 | 12 | -0.001489, -0.004242, -0.0003432 |
| 5 (yes) | 1e+06 | 0.4299 | 1.17e-6 | 2 | -0.02057, -0.02466, -0.000187 |
| 5 (yes) | 1e+07 | 0.4299 | 5.98e-8 | 2 | -0.02085, -0.02497, -0.000187 |
| 5 (yes) | 1e+08 | 0.4299 | 1.16e-8 | 2 | -0.02087, -0.025, -0.000187 |

(4) The sampled-data loop at the step h = 0.5 with label noise drawn afresh (σ = 0.3, 16 paths, 20000 steps) from the lazy projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror update η⁺ = η + hg on each particle's curve — Q's change (mean over the particles and paths, largest), the projection of the end state: its distance to the start's (mean over the paths), to basis pursuit's, its ‖w‖₁; the end's residual:

| instance | the start: ‖w‖₁, ‖w − w_bp‖∞ | the update | Q − Q(0): mean, max | ‖P − P(0)‖∞ | ‖P − w_bp‖∞ | ‖P‖₁ | |Φw − y| |
|---|---|---|---|---|---|---|---|
| 0 | 0.7060, 0.196 | Euler | -3.379e-01, 0.411 | 0.03 | 0.214 | 0.7782 | 0.021 |
| 0 | 0.7060, 0.196 | mirror | +9.252e-18, 1.11e-16 | 2.95e-15 | 0.196 | 0.7060 | 0.046 |
| 1 | 0.5690, 0.123 | Euler | -3.559e-01, 0.477 | 0.0429 | 0.132 | 0.7139 | 0.022 |
| 1 | 0.5690, 0.123 | mirror | +5.204e-18, 1.11e-16 | 4.36e-15 | 0.123 | 0.5690 | 0.048 |

Proposition D‴.44's reduced flow (the step's filter) to the same s = (m v/4)t: its ‖w‖₁ and mean change of Q, and the Euler paths' mean projection against its state: instance 0 (s = 270): 0.7793, -3.381e-01, 7.42e-4; instance 1 (s = 270): 0.7232, -3.552e-01, 0.00568.

## D. The real loop (k_e = 5, γ = 0.5, m = 6, six Gaussian atoms of width 0.4 on [−1.2, 1.2] and a second particle on the third; x_m = sin t + 0.4 sin 2.7t; 2000 time units)

The pair's Bregman split of w*₃ = 0.05111: 0.03685 + 0.01426.

| the adaptation's error | V(0) | V(T) | V monotone | max |V̇ + k_ee²|/max k_ee² | the pair's mirror coordinates (max |Δη_i − Δη_j|, relative) | max |Q − Q(0)| | max |e| (last 10%) | the distinct atoms' |w − w*| | the pair's sum, split |
|---|---|---|---|---|---|---|---|---|---|
| noiseless | 0.175 | 0.022 | yes | 2.83e-9 | 4.88e-15 | 9.99e-15 | 0.00544 | 0.048 | 0.0272, 0.0181 |
| smooth sensor noise | 0.175 | 0.136 | no | — | 6.12e-13 | 6.78e-12 | 0.148 | 0.247 | 0.0234, 0.0157 |

## E. The conserved law (Σβ = m) from the same start (α = 0, β = m/d; t = 10⁶)

| instance | its residual | max |Q − Q(0)| | ‖w‖₁: conserved, free (the projection) | the ends' distance | Σβ |
|---|---|---|---|---|---|
| 0 | 4.16e-17 | 0.0817 | 0.7097, 0.7060 | 0.00158 | 12.000000 |
| 1 | 1.43e-17 | 0.0584 | 0.5721, 0.5690 | 0.00177 | 12.000000 |
| 2 | 6.80e-17 | 0.186 | 0.7886, 0.7899 | 0.00219 | 12.000000 |
| 3 | 4.13e-17 | 0.224 | 0.7454, 0.7513 | 0.00409 | 12.000000 |
| 4 | 4.86e-17 | 0.191 | 0.7529, 0.7598 | 0.00388 | 12.000000 |
| 5 | 7.47e-17 | 0.0984 | 0.7181, 0.7133 | 0.00209 | 12.000000 |

Numbers are produced by `experiments/exp_mirror_free.py`; nothing is edited by hand.
