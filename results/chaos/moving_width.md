# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)

The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = 0.12 (`rich_transit_nd.BumpsND`); 21 kernels started at −1, −0.9, …, 1; 321 samples of [-1.6, 1.6]; exact: as many formed kernels (|w| above 1% of the smallest amplitude) as features, each feature's centre, width and amplitude within 1e-06, the residual's rms below 1e-08.

## A. One feature

A = 0.8 at c* = 0.03, the width σ* = ratio·σ₀; moving widths: the kernels, the largest error of (c, σ, w), the residual; fixed widths (up to 4 entries, 3e+05 evaluations — † the budget spent first): the kernels and the residual:

| σ*/σ₀ | moving: kernels | exact | max error (c, σ, w) | residual | fixed: kernels | residual |
|---|---|---|---|---|---|---|
| 0.25 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 6.8e-17 | 2† | 7.49e-02 |
| 0.5 | 1 | yes | 0.0e+00, 1.2e-16, 1.4e-15 | 2.3e-16 | 2† | 6.51e-02 |
| 0.75 | 1 | yes | 3.5e-18, 1.4e-17, 5.6e-16 | 1.1e-16 | 2† | 3.57e-02 |
| 1.031 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 1.5e-16 | 4 | 1.30e-06 |
| 1.118 | 1 | yes | 0.0e+00, 0.0e+00, 5.6e-16 | 1.5e-16 | 3† | 1.10e-03 |
| 1.414 | 1 | yes | 1.4e-17, 0.0e+00, 5.6e-16 | 1.7e-16 | 4 | 2.19e-03 |
| 2.236 | 1 | yes | 3.5e-18, 0.0e+00, 5.6e-16 | 2.1e-16 | 4 | 2.90e-02 |
| 3 | 1 | yes | 3.5e-17, 0.0e+00, 5.6e-16 | 2.4e-16 | 4 | 7.17e-02 |
| 4 | 1 | yes | 1.7e-17, 0.0e+00, 5.6e-16 | 2.8e-16 | 4 | 1.35e-01 |

## B. The one-kernel reduced energy (the continuum)

J(c, σ) = ½‖f*‖² − ½⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² on a 241 × 241 grid (c within ±6σ*, σ from σ*/8 to 8σ*, A = 1, c* = 0): its grid-local minima and the closed-form gradient's largest relative error against central differences (20 random points):

| σ*/σ₀ | local minima (c, σ/σ*) | min J | gradient error |
|---|---|---|---|
| 0.25 | (+0.0e+00, 1.000) | 0.0e+00 | 4.6e-10 |
| 0.5 | (+0.0e+00, 1.000) | 6.9e-18 | 2.2e-09 |
| 0.75 | (+0.0e+00, 1.000) | 0.0e+00 | 1.3e-09 |
| 1.031 | (+0.0e+00, 1.000) | 0.0e+00 | 6.9e-10 |
| 1.118 | (+0.0e+00, 1.000) | 0.0e+00 | 8.9e-10 |
| 1.414 | (+0.0e+00, 1.000) | -2.8e-17 | 1.3e-09 |
| 2.236 | (+0.0e+00, 1.000) | 5.6e-17 | 1.4e-09 |
| 3 | (+0.0e+00, 1.000) | 0.0e+00 | 1.0e-09 |
| 4 | (+0.0e+00, 1.000) | 0.0e+00 | 7.9e-10 |

## C. Random three-feature targets

40 targets per minimal separation (seed 590; widths log-uniform in [0.05, 0.3], amplitudes ±U(0.4, 1)):

| Δ ≥ | moving widths: exact | formed kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median [range]) | residual (median [max]) |
|---|---|---|---|---|---|
| 0.4 | 36 of 40 | 3 [1, 4] | 3.7e-16 [2.7e-01] | 4 [4, 4] (0† stopped by the budget) | 8.9e-02 [1.9e-01] |
| 0.6 | 35 of 40 | 3 [1, 4] | 3.8e-16 [2.2e-01] | 4 [4, 4] (0† stopped by the budget) | 1.1e-01 [2.1e-01] |

The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; a coinciding pair within 0.01 in |Δc|/σ + |Δlog σ|; a budget stop a dipole run-off if G's stalled block has an opposite-sign pair under 0.05 with |w| beyond 10 times the largest amplitude): Δ ≥ 0.4: 36 exact, 2 coinciding pair, 0 cancellation, 2 dipole run-off, 0 stopped by the budget, 0 not fitted; Δ ≥ 0.6: 35 exact, 1 coinciding pair, 0 cancellation, 4 dipole run-off, 0 stopped by the budget, 0 not fitted.

Not exact with moving widths (the state before a stalled block; G gives the block): Δ ≥ 0.4, target 1 (widths 0.064, 0.273, 0.079; centres -0.06, +0.35, -0.61; amplitudes -0.85, +0.78, +0.64): coinciding pair, 4 kernels (the pair: +0.334/+0.444, Δc 9.57e-04, the width split +7.9453e-08 (predicted +7.9449e-08)), residual 2.5e-13; Δ ≥ 0.4, target 7 (widths 0.064, 0.234, 0.077; centres +0.27, -0.68, -0.23; amplitudes +0.58, -0.91, -0.43): coinciding pair, 4 kernels (the pair: -0.786/+0.357, Δc 1.53e-04, the width split -1.3404e-07 (predicted -1.3405e-07)), residual 1.5e-13; Δ ≥ 0.4, target 9 (widths 0.141, 0.176, 0.144; centres -0.61, +0.18, +0.65; amplitudes +0.79, -0.85, +0.56): dipole run-off, 2 kernels, residual 2.7e-01; Δ ≥ 0.4, target 18 (widths 0.102, 0.280, 0.185; centres -0.58, -0.02, +0.54; amplitudes +0.72, -0.51, +0.60): dipole run-off, 2 kernels, residual 2.2e-01; Δ ≥ 0.6, target 2 (widths 0.062, 0.295, 0.071; centres +0.69, +0.07, -0.70; amplitudes +0.43, +0.82, +0.47): coinciding pair, 4 kernels (the pair: +0.231/+0.587, Δc 2.78e-08, the width split +3.2834e-08 (predicted +1.8958e-16)), residual 1.0e-15; Δ ≥ 0.6, target 12 (widths 0.219, 0.274, 0.133; centres +0.04, -0.67, +0.68; amplitudes +0.83, -0.49, -0.46): dipole run-off, 2 kernels, residual 2.2e-01; Δ ≥ 0.6, target 20 (widths 0.080, 0.123, 0.272; centres -0.63, +0.01, +0.66; amplitudes +0.72, +0.83, -0.51): dipole run-off, 3 kernels, residual 1.6e-01; Δ ≥ 0.6, target 27 (widths 0.096, 0.124, 0.142; centres -0.68, +0.63, -0.05; amplitudes +0.48, +0.48, -0.56): dipole run-off, 2 kernels, residual 1.7e-01; Δ ≥ 0.6, target 29 (widths 0.197, 0.054, 0.088; centres +0.04, +0.69, -0.64; amplitudes -0.59, +0.85, +0.72): dipole run-off, 2 kernels, residual 2.2e-01. The width split (σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), the pair's third moment cancelled) against the observed.

## D. A narrow feature on a wide pedestal

The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):

| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |
|---|---|---|---|---|---|---|
| 0 | 3 | coinciding pair (+0.250/+0.250, Δc 2.94e-04, the width split -1.6515e-15 (predicted -4.0779e-18)) | 1.5e-04, 1.8e-07, 2.5e-01 | 2.6e-13 | 5 | 3.53e-02 |
| 0.1 | 2 | exact | 7.8e-18, 6.9e-18, 5.6e-16 | 1.8e-16 | 4 | 4.35e-02 |
| 0.2 | 2 | exact | 1.5e-17, 6.9e-18, 5.6e-16 | 2.1e-16 | 4 | 4.52e-02 |
| 0.3 | 2 | exact | 1.0e-15, 3.3e-16, 7.2e-15 | 1.0e-15 | 4 | 3.68e-02 |
| 0.45 | 2 | exact | 5.3e-18, 0.0e+00, 5.6e-16 | 2.2e-16 | 4 | 5.05e-02 |
| 0.6 | 3 | coinciding pair (+0.606/+0.194, Δc 1.03e-03, the width split +3.0498e-07 (predicted +3.0516e-07)) | 2.5e-04, 2.5e-07, 1.9e-01 | 1.9e-13 | 4 | 7.25e-02 |

## E. Under noise: the discrepancy stop

C's first 20 targets at Δ ≥ 0.6, one draw each, the path with moving widths to 5 entries (past the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against the data) is within (1 + 0.1)s; its kernels and its error against the truth (rms on the samples), median [range] and the rms over the draws, against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3, the rms of the errors), and per draw against the oracle fit (the K-kernel least-squares fit of the data started at the truth):

| s | stopped | kernels | error: median [range]; rms | the oracle | the stopped error over the oracle fit's: median [range] (within 10⁻⁴ of it) | not stopped: the last block's rms |
|---|---|---|---|---|---|---|
| 0.001 | 19 of 20 | 3 [3, 4] | 1.42e-04 [7.68e-05, 2.48e-04]; 1.648e-04 | 1.674e-04 | 1.0000 [1.0000, 1.0486] (18 of 19) | 2.16e-01 (budget) |
| 0.01 | 19 of 20 | 3 [3, 4] | 1.77e-03 [6.58e-04, 2.17e-03]; 1.750e-03 | 1.674e-03 | 1.0000 [1.0000, 1.2670] (18 of 19) | 2.16e-01 (budget) |

## F. The real loop at a small mass

The free law with moving centres and widths from α = 0, masses ε = 1e-06 (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to 1.3, 6 times the path's last entry, on C's first 4 targets at Δ ≥ 0.6: its formed kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest |w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair and residual:

| target | the path's kernels (closest pair; residual) | times | formed by the loop | max \|Δθ\| | max \|Δw\| | the others' largest \|w\| | the loop's closest pair: distance (\|Δc\|, \|Δσ\|) | residual |
|---|---|---|---|---|---|---|---|---|
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1.3 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 6 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01) | 4.8e-09 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1.3 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 6 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03) | 9.0e-08 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1.3 | [3, 10, 11, 17] | 9.0e-04 | 9.8e-02 | 5.4e-06 | 7.0e-03 (1.5e-03, 5.5e-04) | 5.2e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 6 | [3, 10, 11, 17] | 1.3e-03 | 9.8e-02 | 5.4e-06 | 8.8e-03 (2.2e-03, 3.9e-04) | 5.0e-07 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1.3 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 6 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02) | 1.4e-08 |

## G. The blocks that never settle: a dipole

C's paths stopped by their budget. The stalled block (C's path again; its state when the budget ran out): its kernels, the largest |w|, the closest opposite-sign pair's distance (|Δc|/σ + |Δlog σ|), its residual against the optimal dipole fit's (below), its polishes (status/evaluations/accepted, the largest coefficient over the block's). The block from its entry, flowed on the critical curve to s = 1e+11: the largest |w| and the pair's distance at s = 1000, 100000, 1e+07, 1e+09, 3.2e+10, 1e+11; over the last factor in s the exponents of the largest |w|, the distance, the pair's |Δc| and |Δσ| and the residual's excess over the dipole infimum; the dipole law's constant K = v²(E − E_inf) at the last two outputs (v the pair's half-difference |w_a − w_b|/2, E_inf the optimal dipole fit's — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄) by least squares), the growth of v³ over the last interval against 3μ₀K (μ₀ = 4Gc/A: the coefficient's mobility μ₀|w| on the exponential branch), the moment's σ-share |p_σ|/|p|:

| Δ ≥ | target | entries | the stalled block: kernels, max\|w\|, distance, residual (infimum), polishes | the flow: largest \|w\| | distance | exponents: \|w\|, distance, \|Δc\|, \|Δσ\|, the residual's excess | K | Δv³/Δs, 3μ₀K | \|p_σ\|/\|p\| |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 2 | 2, 3910, 3.12e-04, 1.229967e-01 (1.229967e-01), 3/3358/a/846 | 0.895, 4.62, 20.2, 92, 291, 425 | 5.95e+00, 2.91e-01, 6.18e-02, 1.33e-02, 4.20e-03, 2.87e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 2.55977e-03, 2.55976e-03 | 7.67931e-04, 7.67927e-04 | 0.962 |
| 0.4 | 18 | 2 | 2, 3964, 2.30e-04, 9.689918e-02 (9.689918e-02), 3/2904/a/1.01e+03 | 0.572, 3.92, 17.8, 81.4, 258, 377 | 4.88e+00, 2.53e-01, 5.20e-02, 1.12e-02, 3.53e-03, 2.42e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 1.77965e-03, 1.77964e-03 | 5.33894e-04, 5.33893e-04 | 0.987 |
| 0.6 | 12 | 2 | 2, 1821, 6.13e-04, 5.423548e-02 (5.423547e-02), 3/4572/a/592 | 0.867, 3.07, 9.64, 40.4, 127, 185 | 5.55e+00, 4.29e-01, 1.21e-01, 2.79e-02, 8.82e-03, 6.03e-03 | +0.3322, -0.3331, -0.3331, -0.3331, -0.6663 | 2.10595e-04, 2.10560e-04 | 6.31839e-05, 6.31681e-05 | 0.956 |
| 0.6 | 20 | 3 | 3, 1431, 4.08e-04, 1.185147e-01 (1.185147e-01), 3/2968/a/1.04e+03 | 0.906, 1.37, 7.25, 32.9, 104, 152 | 1.80e+00, 5.19e-01, 8.31e-02, 1.78e-02, 5.62e-03, 3.84e-03 | +0.3328, -0.3333, -0.3333, -0.3333, -0.6667 | 1.16958e-04, 1.16958e-04 | 3.50872e-05, 3.50873e-05 | 0.734 |
| 0.6 | 27 | 2 | 2, 3002, 2.30e-04, 7.510333e-02 (7.510333e-02), 3/2981/a/1.09e+03 | 0.566, 2.75, 12.5, 57.1, 181, 264 | 5.16e+00, 2.78e-01, 5.63e-02, 1.21e-02, 3.82e-03, 2.61e-03 | +0.3329, -0.3333, -0.3333, -0.3333, -0.6667 | 6.10916e-04, 6.10914e-04 | 1.83275e-04, 1.83274e-04 | 0.994 |
| 0.6 | 29 | 2 | 2, 3956, 1.92e-04, 1.430698e-01 (1.430698e-01), 3/2876/a/1.12e+03 | 0.614, 3.52, 16.2, 74.2, 235, 343 | 4.63e+00, 2.38e-01, 4.78e-02, 1.03e-02, 3.24e-03, 2.22e-03 | +0.3330, -0.3333, -0.3333, -0.3333, -0.6667 | 1.34468e-03, 1.34468e-03 | 4.03404e-04, 4.03404e-04 | 0.987 |

The real loop at the masses ε from α = 0 at the times τ: its formed kernels, the closest pair among them (its distance, coefficients, and its midpoint's distance from the nearest feature in that feature's width), the coinciding pairs (within 0.01), the largest |w| over the largest amplitude, the residual:

| Δ ≥ | target | ε | τ | formed | the closest pair: distance, coefficients, offset (σ*) | coinciding pairs | max\|w\|/max\|A\| | residual |
|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 1e-06 | 400 | 4 | 1.12e-03, +4.526/-3.962, 0.008 | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 1600 | 4 | 1.12e-03, +4.526/-3.962, 0.008 | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 6400 | 4 | 1.12e-03, +4.526/-3.962, 0.008 | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-08 | 400 | 4 | 4.65e-04, +9.749/-9.185, 0.008 | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 1600 | 4 | 4.65e-04, +9.749/-9.185, 0.008 | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 6400 | 4 | 4.65e-04, +9.749/-9.185, 0.008 | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 18 | 1e-06 | 400 | 4 | 1.29e-01, -4.429/+5.218, 0.988 | none | 7.27 | 1.0e-03 |
| 0.4 | 18 | 1e-06 | 1600 | 4 | 2.29e-02, -4.522/+5.115, 0.175 | none | 7.13 | 2.3e-06 |
| 0.4 | 18 | 1e-06 | 6400 | 4 | 1.47e-02, -4.519/+5.118, 0.114 | none | 7.13 | 4.3e-07 |
| 0.4 | 18 | 1e-08 | 400 | 4 | 2.13e-02, -5.588/+5.089, 0.210 | none | 7.79 | 1.4e-06 |
| 0.4 | 18 | 1e-08 | 1600 | 4 | 1.61e-02, -5.591/+5.087, 0.161 | none | 7.79 | 4.5e-07 |
| 0.4 | 18 | 1e-08 | 6400 | 4 | 1.24e-02, -5.592/+5.086, 0.125 | none | 7.79 | 1.7e-07 |
| 0.6 | 12 | 1e-06 | 400 | 4 | 1.55e-01, -5.185/+4.604, 1.214 | none | 6.26 | 4.1e-04 |
| 0.6 | 12 | 1e-06 | 1600 | 5 | 1.31e-01, -5.164/+4.637, 1.200 | none | 6.23 | 2.2e-04 |
| 0.6 | 12 | 1e-06 | 6400 | 5 | 8.21e-02, -5.103/+4.687, 1.154 | none | 6.16 | 1.3e-04 |
| 0.6 | 12 | 1e-08 | 400 | 4 | 2.14e-02, -5.286/+6.107, 0.141 | none | 7.37 | 1.0e-06 |
| 0.6 | 12 | 1e-08 | 1600 | 4 | 1.51e-02, -5.285/+6.109, 0.101 | none | 7.37 | 2.6e-07 |
| 0.6 | 12 | 1e-08 | 6400 | 4 | 1.12e-02, -5.284/+6.110, 0.075 | none | 7.38 | 8.5e-08 |
| 0.6 | 20 | 1e-06 | 400 | 4 | 3.86e-01, -0.227/-0.289, 0.023 | none | 1.00 | 1.2e-05 |
| 0.6 | 20 | 1e-06 | 1600 | 4 | 2.76e-01, -0.225/-0.287, 0.017 | none | 1.00 | 3.2e-06 |
| 0.6 | 20 | 1e-06 | 6400 | 4 | 1.99e-01, -0.224/-0.286, 0.012 | none | 1.00 | 8.7e-07 |
| 0.6 | 20 | 1e-08 | 400 | 4 | 3.02e-02, +1.865/-2.369, 0.122 | none | 2.86 | 7.5e-07 |
| 0.6 | 20 | 1e-08 | 1600 | 4 | 2.29e-02, +1.865/-2.369, 0.093 | none | 2.86 | 2.5e-07 |
| 0.6 | 20 | 1e-08 | 6400 | 4 | 1.77e-02, +1.864/-2.370, 0.072 | none | 2.86 | 8.9e-08 |
| 0.6 | 27 | 1e-06 | 400 | 5 | 2.98e-01, -2.075/+3.872, 0.725 | none | 6.92 | 3.8e-05 |
| 0.6 | 27 | 1e-06 | 1600 | 5 | 2.23e-01, -2.061/+3.898, 0.545 | none | 6.96 | 8.9e-06 |
| 0.6 | 27 | 1e-06 | 6400 | 5 | 1.79e-01, -2.055/+3.911, 0.440 | none | 6.99 | 2.7e-06 |
| 0.6 | 27 | 1e-08 | 400 | 4 | 7.52e-03, -5.562/+6.043, 0.088 | -5.562/+6.043 | 10.79 | 2.2e-07 |
| 0.6 | 27 | 1e-08 | 1600 | 4 | 5.54e-03, -5.561/+6.044, 0.065 | -5.561/+6.044 | 10.79 | 6.8e-08 |
| 0.6 | 27 | 1e-08 | 6400 | 4 | 4.21e-03, -5.561/+6.044, 0.050 | -5.561/+6.044 | 10.80 | 2.3e-08 |
| 0.6 | 29 | 1e-06 | 400 | 4 | 8.57e-02, +0.993/-1.578, 0.177 | none | 1.85 | 5.8e-06 |
| 0.6 | 29 | 1e-06 | 1600 | 4 | 5.98e-02, +0.991/-1.580, 0.125 | none | 1.85 | 1.5e-06 |
| 0.6 | 29 | 1e-06 | 6400 | 4 | 4.33e-02, +0.990/-1.581, 0.091 | none | 1.85 | 4.1e-07 |
| 0.6 | 29 | 1e-08 | 400 | 4 | 1.89e-02, +2.798/-3.389, 0.096 | none | 3.96 | 4.7e-07 |
| 0.6 | 29 | 1e-08 | 1600 | 4 | 1.40e-02, +2.798/-3.390, 0.072 | none | 3.96 | 1.5e-07 |
| 0.6 | 29 | 1e-08 | 6400 | 4 | 1.07e-02, +2.797/-3.390, 0.055 | none | 3.96 | 5.0e-08 |

## H. Two features of one sign, close

Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 unless stated — closer than the kernels' starting width σ₀ = 0.12 —, the moving-width path (as C's, its budget 1e+06): the outcome, the most kernels of the other sign in a block along the path and the most negative coefficient, the end's formed kernels of the other sign, the smallest formed width, the largest |w|, the residual:

| target | outcome | negative kernels along the path (most negative w) | at the end | smallest width | max\|w\| | residual |
|---|---|---|---|---|---|---|
| σ* 0.04, equal | coinciding pair | 10 (-0.427) | 2 | 4.00e-02 | 1 | 1.3e-13 |
| σ* 0.06, equal | cancellation | 14 (-0.610) | 7 | 4.96e-02 | 3.24 | 4.3e-10 |
| σ* 0.06, A₂ 0.99 | cancellation | 13 (-0.573) | 6 | 5.19e-02 | 2.93 | 9.6e-10 |
| σ* 0.06, A₂ 0.95 | cancellation | 14 (-0.587) | 7 | 4.92e-02 | 3.14 | 6.0e-10 |
| σ* 0.06, A₂ 0.9 | cancellation | 12 (-0.555) | 7 | 5.12e-02 | 2.84 | 8.1e-10 |
| σ* 0.06, c₂ +0.01 | cancellation | 11 (-0.561) | 7 | 1.75e-04 | 2.85 | 4.4e-10 |
| σ* 0.06, σ₂ 0.065 | cancellation | 14 (-0.603) | 7 | 4.94e-02 | 3.16 | 6.1e-10 |
| σ* 0.04, A₂ 0.6 | exact | 0 (+0.600) | 0 | 4.00e-02 | 1 | 5.0e-17 |
| σ* 0.06, A₂ 0.6 | exact | 0 (+0.600) | 0 | 6.00e-02 | 1 | 2.5e-16 |

Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.
