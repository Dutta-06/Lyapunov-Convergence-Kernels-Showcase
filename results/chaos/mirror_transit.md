# Proposition D‴.48: the rich transit (`exp_mirror_transit.py`)

The free mass–amplitude law (Proposition D‴.45; m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025) from small masses, in the time rescaled by the transit's scale (`mirror_free.iss_path`, `rich_transit`).

## A. From α = 0: the inverse scale space flow on the scale Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}

Three 6 × 14 Gaussian designs, 3-sparse truths, masses εb_i (b_i log-uniform in [1/4, 4]); the path of Σ|w_i|/√b_i against the loop in τ = t/Λ, u = η/Λ, on 2001 times in [0, 1.5τ*] (τ* the last event): sup_τ|u − u_ISS|, the residual's rms difference over τ, the largest difference of the entry times (the first sample where |w_i| > 10⁻³ max|w_bp|, on the loop and on the path; in brackets the atoms crossing on one of them only; the grid's step given), the end's distance to weighted basis pursuit.

Design 3: the path's events — enter 8 at 0.7057; enter 13 at 1.356; enter 7 at 2.028; its end against weighted basis pursuit 2.22e-16; the grid's step 0.00152.

| ε | Λ | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 0.536 | 1.17 | 0.931 (4) | 0.368 |
| 0.001 | 222.1 | 0.203 | 0.728 | 0.389 (3) | 0.0171 |
| 0.0001 | 702.5 | 0.0673 | 0.42 | 0.14 (0) | 1.84e-4 |
| 1e-05 | 2221 | 0.0216 | 0.238 | 0.0441 (0) | 2.82e-6 |
| 1e-06 | 7025 | 0.00685 | 0.133 | 0.0137 (0) | 7.06e-8 |
| 1e-07 | 2.221e+04 | 0.00217 | 0.0726 | 0.00456 (0) | 2.08e-9 |
| 1e-08 | 7.025e+04 | 6.87e-4 | 0.0402 | 0.00152 (0) | 6.42e-11 |

sup|u − u_ISS| between ε = 10⁻⁶ and 10⁻⁸: the slope 0.500 in log ε (1/Λ ∝ ε^{1/2}: 0.5).

Design 4: the path's events — enter 0 at 1.226; enter 4 at 1.453; enter 2 at 19.16; enter 11 at 38.43; enter 10 at 39.41; enter 3 at 55.48; exit 10 at 55.48; enter 7 at 78.31; its end against weighted basis pursuit 4.54e-16; the grid's step 0.0587.

| ε | Λ | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 0.61 | 0.172 | 57.6 (3) | 0.02 |
| 0.001 | 222.1 | 0.235 | 0.106 | 42.9 (0) | 0.00345 |
| 0.0001 | 702.5 | 0.0865 | 0.065 | 37.5 (0) | 5.04e-4 |
| 1e-05 | 2221 | 0.0284 | 0.0498 | 11.6 (0) | 6.62e-5 |
| 1e-06 | 7025 | 0.00892 | 0.0154 | 3.76 (0) | 4.92e-6 |
| 1e-07 | 2.221e+04 | 0.00282 | 0.00125 | 1.17 (0) | 1.85e-7 |
| 1e-08 | 7.025e+04 | 8.91e-4 | 6.23e-4 | 0.352 (0) | 6.15e-9 |

sup|u − u_ISS| between ε = 10⁻⁶ and 10⁻⁸: the slope 0.500 in log ε (1/Λ ∝ ε^{1/2}: 0.5).

Design 5: the path's events — enter 0 at 2.413; enter 9 at 2.556; enter 1 at 2.988; enter 2 at 3.18; exit 0 at 3.18; enter 12 at 4.348; enter 3 at 4.615; exit 9 at 4.615; enter 8 at 12.15; enter 10 at 58.54; its end against weighted basis pursuit 3.33e-16; the grid's step 0.0439.

| ε | Λ | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 1.15 | 0.235 | 55.6 (5) | 0.178 |
| 0.001 | 222.1 | 0.5 | 0.147 | 53 (0) | 0.00341 |
| 0.0001 | 702.5 | 0.166 | 0.0847 | 26.7 (0) | 4.37e-6 |
| 1e-05 | 2221 | 0.0534 | 0.0453 | 8.43 (0) | 2.17e-7 |
| 1e-06 | 7025 | 0.017 | 0.0169 | 2.68 (0) | 8.04e-9 |
| 1e-07 | 2.221e+04 | 0.00537 | 0.0101 | 0.834 (0) | 2.95e-10 |
| 1e-08 | 7.025e+04 | 0.00194 | 0.00484 | 0.263 (0) | 4.00e-11 |

sup|u − u_ISS| between ε = 10⁻⁶ and 10⁻⁸: the slope 0.471 in log ε (1/Λ ∝ ε^{1/2}: 0.5).

## B. On branches: the logarithmic scale

Design 3 of A; amplitudes s_i·U(0.2, 0.6) (the truth's signs on its atoms, random elsewhere), masses ε; the one-sided path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in τ = t/log(1/ε), u = (η − η(0))/log(1/ε). The path's events: enter 8 at 21.77; enter 13 at 35.26; enter 7 at 92.15; its end against the sign-constrained weighted basis pursuit 4.44e-16.

| ε | all on branches | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_lim|∞ |
|---|---|---|---|---|---|
| 0.0001 | yes | 15.7 | 0.652 | 55 (0) | 1.26e-5 |
| 1e-08 | yes | 7.97 | 0.465 | 27.9 (0) | 3.59e-9 |
| 1e-16 | yes | 3.98 | 0.328 | 14 (0) | 4.44e-15 |
| 1e-32 | yes | 1.99 | 0.232 | 6.98 (0) | 6.66e-16 |
| 1e-64 | yes | 0.996 | 0.164 | 3.52 (0) | 7.77e-16 |

A ± pair per atom (amplitudes ±0.4, equal masses, w(0) = u): the atom's common increment over log(1/u) against the two-sided path with the weights 1/(γa₀):

| u | sup|u − u_ISS| | residual: rms difference | |W(1.5τ*) − w_end|∞ | entries on the path |
|---|---|---|---|---|
| 0.0001 | 20.8 | 1.49 | 0.0396 | 3 |
| 1e-08 | 7.86 | 0.998 | 0.00712 | 3 |
| 1e-16 | 3.96 | 0.705 | 2.32e-4 | 3 |
| 1e-32 | 1.98 | 0.493 | 1.94e-6 | 3 |
| 1e-64 | 0.991 | 0.347 | 2.78e-10 | 3 |

## C. Moving centres

Gaussian kernels of width 0.12 on 61 samples of [−1, 1]; 21 kernels started on a grid of [−1, 1] (spacing 0.1), α = 0, masses ε; the target 3 bumps of the same width at [-0.52, 0.07, 0.61] with amplitudes [0.9, -0.6, 0.75]. A kernel has formed when |w| > 0.01·max|A*|.

(1) The paradigm's metric (ċ = Gβa⟨r∂_cφ⟩). The sliding inverse scale space: enter kernel 5 (start -0.5) at τ = 10.72; enter kernel 16 (start +0.6) at τ = 12.85; enter kernel 11 (start +0.1) at τ = 16.39; it ends with 3 kernels at [-0.52, 0.07, 0.61] with coefficients [0.9, -0.6, 0.75] (rms residual 1.36e-16). The loop (the entries at the first sample — step 0.0484 in τ — where |w| > 0.01·max|A*|):

| ε | Λ | entries (kernel: τ) | formed | their centres | their coefficients | rms residual | inactive: max |Δc|, max |w| |
|---|---|---|---|---|---|---|---|
| 0.001 | 222.1 | 5: 10.7, 4: 11.23, 6: 12.25, 16: 12.83, 17: 13.7, 15: 14.52, 11: 16.22, 10: 16.8 | 8 | [-0.52, -0.52, -0.5202, 0.07, 0.07, 0.6103, 0.6101, 0.6101] | [0.2103, 0.6566, 0.0323, -0.1823, -0.3884, 0.0382, 0.5812, 0.1299] | 0.00702 | 0.181, 0.00714 |
| 0.0001 | 702.5 | 5: 10.75, 16: 12.88, 11: 16.36, 10: 17.53 | 4 | [-0.52, 0.07, 0.07, 0.61] | [0.8969, -0.0127, -0.5872, 0.7493] | 1.21e-6 | 0.0877, 0.00312 |
| 1e-05 | 2221 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 3.80e-7 | 0.0331, 2.91e-5 |
| 1e-06 | 7025 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 1.55e-8 | 0.0026, 4.26e-7 |
| 1e-07 | 2.221e+04 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 4.40e-10 | 2.34e-4, 1.12e-8 |

(2) The Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩): before the first entry (τ = 10) the inactive kernels' limit (du/dτ = ⟨f*φ⟩, dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩) gathers the 21 kernels at 3 positions [-0.521, 0.068, 0.613]; the loop's centres at τ = 10 against that limit, and its end:

| ε | max |c − c_lim| at τ = 10 | max |w| at τ = 10 | formed | their positions (clusters) | rms residual |
|---|---|---|---|---|---|
| 0.0001 | 2.55e-6 | 3.42e-4 | 15 | [-0.52, 0.07, 0.61] (3) | 5.32e-10 |
| 1e-05 | 1.05e-7 | 1.48e-5 | 12 | [-0.52, 0.07, 0.61] (3) | 1.25e-10 |
| 1e-06 | 3.44e-9 | 4.86e-7 | 8 | [-0.52, 0.07, 0.61] (3) | 1.05e-10 |

(3) The centres held on the grid (ε = 1e-06, to τ = 63.9): 5 kernels formed, rms residual 0.0192 (the grid's least-squares fit: 8.38e-5).

(4) 16 random targets (three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)): the sliding inverse scale space ends with exactly three kernels within 10⁻³ of the bumps' centres on 16 of 16; the loop with the paradigm's metric on 11 of 16 at ε = 1e-05 (formed kernels [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5], the inactive kernels' largest displacement 0.12) and on 13 of 16 at ε = 1e-06 (formed kernels [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4], the inactive kernels' largest displacement 0.0699; largest rms residual 6.69e-5) and on 16 of 16 at ε = 1e-08 (the entering kernels those of the sliding path on 14 of 16); with the Wasserstein–Fisher–Rao metric (ε = 1e-06) it forms [5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7] kernels in [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4] clusters (largest rms residual 1.56e-6); with the centres held on the grid [4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6] kernels (rms residual 0.0114–0.0257).

(5) Harder targets (ε = 10⁻⁵, the paradigm's metric):

| target | sliding ISS: kernels, their centres, rms residual (events) | loop: formed, centres, coefficients (clusters) | rms residual |
|---|---|---|---|
| two bumps 0.2 apart | 2, [-0.1, 0.1], 4.95e-10 (2) | 2, [-0.1, 0.1], [0.8, 0.6] (2) | 5.81e-8 |
| two bumps 0.3 apart | 2, [-0.15, 0.15], 7.00e-13 (2) | 2, [-0.15, 0.15], [0.8, 0.6] (2) | 6.57e-8 |
| three wider bumps (width 0.2) | 4, [-0.681, -0.453, 0.06, 0.632], 0.0831 (56, not terminated) | 6, [-0.675, -0.446, -0.097, 0.08, 0.561, 0.782], [0.563, 0.763, -0.081, -0.591, 0.645, 0.41] (6) | 0.0301 |


(6) A lost near-tie (target 10 of (4)): the sliding path enters kernels [2, 14, 12] at τ = [11.492, 17.494, 21.386]. After the second entry the inactive kernels left near their walls linger at their margins μ and creep at ≈ (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩| in τ ((iii)(a)'s bound):

| ε | the lingering kernels (margin < 0.01 after the second entry): kernel, μ, ⟨r∂_cφ⟩, the speed (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩|, the displacement before the third entry | the loop's entrants | the formed kernels' centres |
|---|---|---|---|
| 1e-08 | 3: 0.00434, -9.28e-6, 3.57e-5, 0.0038; 15: 0.00646, -0.033, 0.0385, 0.162 | [2, 14, 15] | [-0.7519, 0.1116, 0.4879] |
| 1e-10 | 3: 0.00484, -1.06e-5, 2.94e-7, 3.19e-5; 15: 0.00696, -0.0326, 3.04e-4, 8.02e-4 | [2, 14, 12] | [-0.7519, 0.1116, 0.4879] |

Numbers are produced by `experiments/exp_mirror_transit.py`; nothing is edited by hand.
