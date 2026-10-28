# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)

The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = 0.12 (`rich_transit_nd.BumpsND`); 21 kernels started at −1, −0.9, …, 1; 321 samples of [-1.6, 1.6]; exact: as many kernels as features, each feature's centre, width and amplitude within 1e-06, the residual's rms below 1e-08.

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

The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude): Δ ≥ 0.4: 36 exact, 2 coinciding pair, 2 stopped by the budget, 0 not fitted; Δ ≥ 0.6: 35 exact, 1 coinciding pair, 4 stopped by the budget, 0 not fitted.

Not exact with moving widths: Δ ≥ 0.4, target 1 (widths 0.064, 0.273, 0.079; centres -0.06, +0.35, -0.61; amplitudes -0.85, +0.78, +0.64): coinciding pair, 4 kernels (the pair's coefficients +0.334/+0.444), residual 2.5e-13; Δ ≥ 0.4, target 7 (widths 0.064, 0.234, 0.077; centres +0.27, -0.68, -0.23; amplitudes +0.58, -0.91, -0.43): coinciding pair, 4 kernels (the pair's coefficients -0.786/+0.357), residual 1.5e-13; Δ ≥ 0.4, target 9 (widths 0.141, 0.176, 0.144; centres -0.61, +0.18, +0.65; amplitudes +0.79, -0.85, +0.56): stopped by the budget, 2 kernels, residual 2.7e-01; Δ ≥ 0.4, target 18 (widths 0.102, 0.280, 0.185; centres -0.58, -0.02, +0.54; amplitudes +0.72, -0.51, +0.60): stopped by the budget, 2 kernels, residual 2.2e-01; Δ ≥ 0.6, target 2 (widths 0.062, 0.295, 0.071; centres +0.69, +0.07, -0.70; amplitudes +0.43, +0.82, +0.47): coinciding pair, 4 kernels (the pair's coefficients +0.231/+0.587), residual 1.0e-15; Δ ≥ 0.6, target 12 (widths 0.219, 0.274, 0.133; centres +0.04, -0.67, +0.68; amplitudes +0.83, -0.49, -0.46): stopped by the budget, 2 kernels, residual 2.2e-01; Δ ≥ 0.6, target 20 (widths 0.080, 0.123, 0.272; centres -0.63, +0.01, +0.66; amplitudes +0.72, +0.83, -0.51): stopped by the budget, 3 kernels, residual 1.6e-01; Δ ≥ 0.6, target 27 (widths 0.096, 0.124, 0.142; centres -0.68, +0.63, -0.05; amplitudes +0.48, +0.48, -0.56): stopped by the budget, 2 kernels, residual 1.7e-01; Δ ≥ 0.6, target 29 (widths 0.197, 0.054, 0.088; centres +0.04, +0.69, -0.64; amplitudes -0.59, +0.85, +0.72): stopped by the budget, 2 kernels, residual 2.2e-01.

## D. A narrow feature on a wide pedestal

The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):

| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |
|---|---|---|---|---|---|---|
| 0 | 3 | coinciding pair (+0.250/+0.250) | 1.5e-04, 1.8e-07, 2.5e-01 | 2.6e-13 | 5 | 3.53e-02 |
| 0.1 | 2 | exact | 7.8e-18, 6.9e-18, 5.6e-16 | 1.8e-16 | 4 | 4.35e-02 |
| 0.2 | 2 | exact | 1.5e-17, 6.9e-18, 5.6e-16 | 2.1e-16 | 4 | 4.52e-02 |
| 0.3 | 2 | exact | 1.0e-15, 3.3e-16, 7.2e-15 | 1.0e-15 | 4 | 3.68e-02 |
| 0.45 | 2 | exact | 5.3e-18, 0.0e+00, 5.6e-16 | 2.2e-16 | 4 | 5.05e-02 |
| 0.6 | 3 | coinciding pair (+0.606/+0.194) | 2.5e-04, 2.5e-07, 1.9e-01 | 1.9e-13 | 4 | 7.25e-02 |

## E. Under noise: the discrepancy stop

C's first 20 targets at Δ ≥ 0.6, one draw each, the path with moving widths to 5 entries (past the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against the data) is within (1 + 0.1)s; its kernels and its error against the truth (rms on the samples), median [range], against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3):

| s | stopped | kernels | error | the oracle |
|---|---|---|---|---|
| 0.001 | 19 of 20 | 3 [3, 4] | 1.42e-04 [7.68e-05, 2.48e-04] | 1.67e-04 |
| 0.01 | 19 of 20 | 3 [3, 4] | 1.77e-03 [6.58e-04, 2.17e-03] | 1.67e-03 |

## F. The real loop at a small mass

The free law with moving centres and widths from α = 0, masses ε = 1e-06 (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to 1.3 times the path's last entry, on C's first 4 targets at Δ ≥ 0.6: its formed kernels (|w| above 1% of the smallest amplitude) against the path's, and the path's kernels' largest deviation:

| target | the path's kernels | formed by the loop | max \|Δθ\| | max \|Δw\| | the others' largest \|w\| |
|---|---|---|---|---|---|
| 0 | [10, 11, 17] | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 |
| 1 | [3, 9, 10] | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 |
| 2 | [3, 10, 11, 17] | [3, 10, 11, 17] | 9.0e-04 | 9.8e-02 | 5.4e-06 |
| 3 | [3, 4, 16] | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 |

## G. The blocks that never settle: a dipole

C's paths stopped by their budget: the stalled block from its entry (the path to that many entries), its flow on the critical curve to s = 1e+11 — the largest |w| and the closest opposite-sign pair (|Δc|/σ + |Δlog σ|) at s = 1000, 100000, 1e+07, 1e+09, 1e+11 —; and the real loop at ε = 1e-06, 1e-08 to τ = 400: its formed kernels, the coinciding pairs' coefficients, the largest |w| over the largest amplitude and the residual:

| Δ ≥ | target | entries | the flow: largest \|w\| | the opposite-sign pair's distance | the last factor 100 in s: the exponents of \|w\| and the distance | residual | the loop: formed, pairs, max\|w\|/max\|A\|, residual |
|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 2 | 0.895, 4.62, 20.2, 92, 425 | 5.95e+00, 2.91e-01, 6.18e-02, 1.33e-02, 2.87e-03 | +0.333, -0.333 | 1.23e-01 | ε = 1e-06: 4, +4.526/-3.962, 5.33, 2.0e-10; ε = 1e-08: 4, +9.749/-9.185, 11.47, 1.1e-11 |
| 0.4 | 18 | 2 | 0.572, 3.92, 17.8, 81.4, 377 | 4.88e+00, 2.53e-01, 5.20e-02, 1.12e-02, 2.42e-03 | +0.333, -0.333 | 9.69e-02 | ε = 1e-06: 4, none, 7.27, 1.0e-03; ε = 1e-08: 4, -5.588/+5.089, 7.79, 1.4e-06 |
| 0.6 | 12 | 2 | 0.867, 3.07, 9.64, 40.4, 185 | 5.55e+00, 4.29e-01, 1.21e-01, 2.79e-02, 6.03e-03 | +0.331, -0.332 | 5.42e-02 | ε = 1e-06: 4, none, 6.26, 4.1e-04; ε = 1e-08: 4, -5.286/+6.107, 7.37, 1.0e-06 |
| 0.6 | 20 | 3 | 0.906, 1.37, 7.25, 32.9, 152 | 1.80e+00, 5.19e-01, 8.31e-02, 1.78e-02, 3.84e-03 | +0.332, -0.333 | 1.19e-01 | ε = 1e-06: 4, none, 1.00, 1.2e-05; ε = 1e-08: 4, +1.865/-2.369, 2.86, 7.5e-07 |
| 0.6 | 27 | 2 | 0.566, 2.75, 12.5, 57.1, 264 | 5.16e+00, 2.78e-01, 5.63e-02, 1.21e-02, 2.61e-03 | +0.333, -0.333 | 7.51e-02 | ε = 1e-06: 5, none, 6.92, 3.8e-05; ε = 1e-08: 4, -5.562/+6.043, 10.79, 2.2e-07 |
| 0.6 | 29 | 2 | 0.614, 3.52, 16.2, 74.2, 343 | 4.63e+00, 2.38e-01, 4.78e-02, 1.03e-02, 2.22e-03 | +0.333, -0.333 | 1.43e-01 | ε = 1e-06: 4, none, 1.85, 5.8e-06; ε = 1e-08: 4, +2.798/-3.389, 3.96, 4.7e-07 |

Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.
