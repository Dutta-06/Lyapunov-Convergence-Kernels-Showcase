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

The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; exact with a null cluster where each feature's nearest kernel is exact and the others sum to nothing; a coinciding pair within 0.01 in |Δc|/σ + |Δlog σ| whose net coefficient and centre of mass fit a feature within 0.01, a null cluster where it fits none; a budget stop a dipole run-off if the stalled block has an opposite-sign pair under 0.05 with |w| beyond 10 times the largest amplitude): Δ ≥ 0.4: 36 exact, 2 coinciding pair, 2 dipole run-off; Δ ≥ 0.6: 35 exact, 1 coinciding pair, 4 dipole run-off.

Not exact with moving widths (the end's state; for a budget stop the state before the stalled block, and the stalled block): Δ ≥ 0.4, target 1 (widths 0.064, 0.273, 0.079; centres -0.06, +0.35, -0.61; amplitudes -0.85, +0.78, +0.64): coinciding pair, 4 kernels (the pair: +0.334/+0.444, Δc 9.57e-04, the width split +7.9453e-08 (predicted +7.9449e-08)), residual 2.5e-13; Δ ≥ 0.4, target 7 (widths 0.064, 0.234, 0.077; centres +0.27, -0.68, -0.23; amplitudes +0.58, -0.91, -0.43): coinciding pair, 4 kernels (the pair: -0.786/+0.357, Δc 1.53e-04, the width split -1.3404e-07 (predicted -1.3405e-07)), residual 1.5e-13; Δ ≥ 0.4, target 9 (widths 0.141, 0.176, 0.144; centres -0.61, +0.18, +0.65; amplitudes +0.79, -0.85, +0.56): dipole run-off, 2 kernels, residual 2.7e-01; the stalled block: 2 kernels, max|w| 3910 (the accepted polish's, solver-dependent), its closest opposite pair 3.12e-04 apart; Δ ≥ 0.4, target 18 (widths 0.102, 0.280, 0.185; centres -0.58, -0.02, +0.54; amplitudes +0.72, -0.51, +0.60): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 3964 (the accepted polish's, solver-dependent), its closest opposite pair 2.30e-04 apart; Δ ≥ 0.6, target 2 (widths 0.062, 0.295, 0.071; centres +0.69, +0.07, -0.70; amplitudes +0.43, +0.82, +0.47): coinciding pair, 4 kernels (the pair: +0.231/+0.587, Δc 2.78e-08, the width split +3.2834e-08 (predicted +1.8958e-16)), residual 1.0e-15; Δ ≥ 0.6, target 12 (widths 0.219, 0.274, 0.133; centres +0.04, -0.67, +0.68; amplitudes +0.83, -0.49, -0.46): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 1821 (the accepted polish's, solver-dependent), its closest opposite pair 6.13e-04 apart; Δ ≥ 0.6, target 20 (widths 0.080, 0.123, 0.272; centres -0.63, +0.01, +0.66; amplitudes +0.72, +0.83, -0.51): dipole run-off, 3 kernels, residual 1.6e-01; the stalled block: 3 kernels, max|w| 1431 (the accepted polish's, solver-dependent), its closest opposite pair 4.08e-04 apart; Δ ≥ 0.6, target 27 (widths 0.096, 0.124, 0.142; centres -0.68, +0.63, -0.05; amplitudes +0.48, +0.48, -0.56): dipole run-off, 2 kernels, residual 1.7e-01; the stalled block: 2 kernels, max|w| 3002 (the accepted polish's, solver-dependent), its closest opposite pair 2.30e-04 apart; Δ ≥ 0.6, target 29 (widths 0.197, 0.054, 0.088; centres +0.04, +0.69, -0.64; amplitudes -0.59, +0.85, +0.72): dipole run-off, 2 kernels, residual 2.2e-01; the stalled block: 2 kernels, max|w| 3956 (the accepted polish's, solver-dependent), its closest opposite pair 1.92e-04 apart. The width split (σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), the pair's third moment cancelled) against the observed.

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

| s | stopped | kernels | error: median [range]; rms | the oracle | the stopped error over the oracle fit's: median [range] (within 10⁻⁴ of it: how many, the largest deviation among them) | not stopped: the last block's rms |
|---|---|---|---|---|---|---|
| 0.001 | 19 of 20 | 3 [3, 4] | 1.42e-04 [7.68e-05, 2.48e-04]; 1.648e-04 | 1.674e-04 | 1.0000 [1.0000, 1.0486] (18 of 19, 7e-08) | 2.16e-01 (budget) |
| 0.01 | 19 of 20 | 3 [3, 4] | 1.77e-03 [6.58e-04, 2.17e-03]; 1.750e-03 | 1.674e-03 | 1.0000 [1.0000, 1.2670] (18 of 19, 7e-08) | 2.16e-01 (budget) |

## F. The real loop at a small mass

The free law with moving centres and widths from α = 0, masses ε = 1e-06, 1e-08 (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to 1.3, 6, 30, 150 times the path's last entry, on C's first 4 targets at Δ ≥ 0.6: its formed kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest |w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair and residual:

| target | the path's kernels (closest pair; residual) | ε | times | formed by the loop | max \|Δθ\| | max \|Δw\| | the others' largest \|w\| | the loop's closest pair: distance (\|Δc\|, \|Δσ\|; coefficients) | residual |
|---|---|---|---|---|---|---|---|---|---|
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 1.3 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 6 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 30 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-06 | 150 | [10, 11, 17] | 6.5e-09 | 7.2e-08 | 8.3e-08 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 4.8e-09 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 1.3 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 6 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 30 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 0 | [10, 11, 17] (4.8e+00; 7.9e-16) | 1e-08 | 150 | [10, 11, 17] | 8.5e-12 | 8.4e-11 | 9.9e-11 | 4.8e+00 (6.9e-01, 1.8e-01; -0.859/-0.519) | 5.9e-12 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 1.3 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 6 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 30 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 9.0e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-06 | 150 | [3, 9, 10] | 7.0e-08 | 9.2e-07 | 1.0e-06 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 8.9e-08 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 1.3 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 6 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 30 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 1 | [3, 9, 10] (9.6e+00; 6.7e-15) | 1e-08 | 150 | [3, 9, 10] | 3.1e-10 | 3.7e-09 | 4.3e-09 | 9.6e+00 (6.7e-01, 2.4e-03; +0.655/+0.773) | 3.5e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 1.3 | [3, 10, 11, 17] | 9.0e-04 | 9.8e-02 | 5.4e-06 | 7.0e-03 (1.5e-03, 5.5e-04; +0.329/+0.489) | 5.2e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 6 | [3, 10, 11, 17] | 1.3e-03 | 9.8e-02 | 5.4e-06 | 8.8e-03 (2.2e-03, 3.9e-04; +0.329/+0.489) | 5.0e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 30 | [3, 10, 11, 17] | 1.9e-03 | 9.8e-02 | 5.3e-06 | 1.2e-02 (3.2e-03, 2.8e-04; +0.329/+0.489) | 4.8e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-06 | 150 | [3, 10, 11, 17] | 2.7e-03 | 9.8e-02 | 5.0e-06 | 1.6e-02 (4.5e-03, 1.9e-04; +0.329/+0.489) | 4.5e-07 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 1.3 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 1.2e-06; +0.329/+0.489) | 2.6e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 6 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 30 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 2 | [3, 10, 11, 17] (2.1e-07; 1.0e-15) | 1e-08 | 150 | [3, 10, 11, 17] | 4.5e-04 | 9.8e-02 | 3.4e-09 | 2.6e-03 (7.6e-04, 8.1e-07; +0.329/+0.489) | 2.1e-10 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 1.3 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 6 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 30 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-06 | 150 | [3, 4, 16] | 3.0e-08 | 1.6e-07 | 2.1e-07 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.4e-08 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 1.3 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 6 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 30 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |
| 3 | [3, 4, 16] (4.0e+00; 1.6e-16) | 1e-08 | 150 | [3, 4, 16] | 3.3e-11 | 1.7e-10 | 2.3e-10 | 4.0e+00 (6.7e-01, 5.8e-02; -0.500/-0.843) | 1.6e-11 |

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

The real loop at the masses ε from α = 0 (`loop_nd`, rtol 10⁻⁹) at τ = 400, 1600, 6400: the formed kernels, the dipole pair (the largest |w| and the closest kernel of the other sign: its distance |Δc|/σ + |Δlog σ|, coefficients, its centre of mass's offset from the nearest feature in that feature's width and its net coefficient W against the feature's amplitude, (W − A)/|A|; landed (L) if both within 0.01), the coinciding pairs (within 0.01), the largest |w| over the largest amplitude, the residual:

| Δ ≥ | target | ε | τ | formed | the dipole pair: distance, coefficients, offset (σ*), (W − A)/\|A\| | coinciding pairs | max\|w\|/max\|A\| | residual |
|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 1e-10 | 400 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-10 | 1600 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-10 | 6400 | 4 | 2.09e-04, +20.712/-20.148, 1.4e-07, -2.8e-05 L | +20.712/-20.148 | 24.37 | 9.4e-12 |
| 0.4 | 9 | 1e-09 | 400 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-09 | 1600 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-09 | 6400 | 4 | 3.10e-04, +14.197/-13.633, 1.5e-07, -2.9e-05 L | +14.197/-13.633 | 16.70 | 9.8e-12 |
| 0.4 | 9 | 1e-08 | 400 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 1600 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-08 | 6400 | 4 | 4.65e-04, +9.749/-9.185, 1.6e-07, -3.0e-05 L | +9.749/-9.185 | 11.47 | 1.1e-11 |
| 0.4 | 9 | 1e-07 | 400 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-07 | 1600 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-07 | 6400 | 4 | 7.05e-04, +6.696/-6.132, 1.7e-07, -3.2e-05 L | +6.696/-6.132 | 7.88 | 1.5e-11 |
| 0.4 | 9 | 1e-06 | 400 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 1600 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-06 | 6400 | 4 | 1.12e-03, +4.526/-3.962, 2.1e-07, -3.5e-05 L | +4.526/-3.962 | 5.33 | 2.0e-10 |
| 0.4 | 9 | 1e-05 | 400 | 4 | 3.12e-01, +3.637/-1.563, 1.3e+00, +2.7e+00 | none | 4.28 | 2.1e-03 |
| 0.4 | 9 | 1e-05 | 1600 | 4 | 2.88e-01, +4.106/-1.818, 1.3e+00, +3.1e+00 | none | 4.83 | 1.9e-03 |
| 0.4 | 9 | 1e-05 | 6400 | 5 | 2.27e-01, +4.290/-1.976, 1.2e+00, +3.1e+00 | none | 5.05 | 6.2e-04 |
| 0.4 | 18 | 1e-10 | 400 | 4 | 6.42e-03, -11.635/+11.130, 7.9e-04, +9.4e-03 L | -11.635/+11.130 | 16.21 | 2.4e-07 |
| 0.4 | 18 | 1e-10 | 1600 | 4 | 4.81e-03, -11.636/+11.129, 3.5e-04, +5.4e-03 L | -11.636/+11.129 | 16.21 | 8.2e-08 |
| 0.4 | 18 | 1e-10 | 6400 | 4 | 3.70e-03, -11.637/+11.128, 1.7e-04, +3.2e-03 L | -11.637/+11.128 | 16.21 | 3.1e-08 |
| 0.4 | 18 | 1e-09 | 400 | 4 | 1.17e-02, -8.025/+7.522, 1.4e-03, +1.4e-02 | none | 11.18 | 5.6e-07 |
| 0.4 | 18 | 1e-09 | 1600 | 4 | 8.79e-03, -8.026/+7.521, 6.5e-04, +8.3e-03 L | -8.026/+7.521 | 11.18 | 1.9e-07 |
| 0.4 | 18 | 1e-09 | 6400 | 4 | 6.75e-03, -8.027/+7.520, 3.1e-04, +5.0e-03 L | -8.027/+7.520 | 11.18 | 7.0e-08 |
| 0.4 | 18 | 1e-08 | 400 | 4 | 2.13e-02, -5.588/+5.089, 2.6e-03, +2.2e-02 | none | 7.79 | 1.4e-06 |
| 0.4 | 18 | 1e-08 | 1600 | 4 | 1.61e-02, -5.591/+5.087, 1.2e-03, +1.3e-02 | none | 7.79 | 4.5e-07 |
| 0.4 | 18 | 1e-08 | 6400 | 4 | 1.24e-02, -5.592/+5.086, 5.9e-04, +7.7e-03 L | none | 7.79 | 1.7e-07 |
| 0.4 | 18 | 1e-07 | 400 | 4 | 1.91e-02, +8.468/-7.882, 3.8e-03, -2.8e-02 | none | 11.80 | 8.1e-06 |
| 0.4 | 18 | 1e-07 | 1600 | 4 | 8.87e-03, +8.474/-7.876, 4.8e-04, -6.7e-03 L | -7.876/+8.474 | 11.81 | 4.7e-07 |
| 0.4 | 18 | 1e-07 | 6400 | 4 | 6.44e-03, +8.475/-7.875, 1.9e-04, -3.6e-03 L | -7.875/+8.475 | 11.81 | 1.4e-07 |
| 0.4 | 18 | 1e-06 | 400 | 4 | 1.29e-01, +5.218/-4.429, 4.7e-01, +3.1e-01 | none | 7.27 | 1.0e-03 |
| 0.4 | 18 | 1e-06 | 1600 | 4 | 2.29e-02, +5.115/-4.522, 1.5e-03, -1.5e-02 | none | 7.13 | 2.3e-06 |
| 0.4 | 18 | 1e-06 | 6400 | 4 | 1.47e-02, +5.118/-4.519, 4.5e-04, -6.4e-03 L | none | 7.13 | 4.3e-07 |
| 0.4 | 18 | 1e-05 | 400 | 4 | 2.23e-01, +3.233/-2.424, 4.8e-01, +3.4e-01 | none | 4.50 | 1.0e-03 |
| 0.4 | 18 | 1e-05 | 1600 | 5 | 1.91e-01, +3.128/-2.515, 2.7e-01, +1.8e-02 | none | 4.36 | 4.2e-04 |
| 0.4 | 18 | 1e-05 | 6400 | 5 | 8.67e-02, +3.091/-2.548, 2.8e-03, -9.9e-02 | none | 4.31 | 6.0e-06 |
| 0.6 | 12 | 1e-10 | 400 | 4 | 6.35e-03, +11.468/-10.642, 1.9e-04, -3.4e-03 L | -10.642/+11.468 | 13.84 | 1.2e-07 |
| 0.6 | 12 | 1e-10 | 1600 | 4 | 4.74e-03, +11.468/-10.642, 8.0e-05, -1.9e-03 L | -10.642/+11.468 | 13.85 | 3.9e-08 |
| 0.6 | 12 | 1e-10 | 6400 | 4 | 3.66e-03, +11.469/-10.641, 3.7e-05, -1.2e-03 L | -10.641/+11.469 | 13.85 | 1.4e-08 |
| 0.6 | 12 | 1e-09 | 400 | 4 | 1.15e-02, +8.265/-7.441, 4.0e-04, -5.6e-03 L | none | 9.98 | 3.3e-07 |
| 0.6 | 12 | 1e-09 | 1600 | 4 | 8.42e-03, +8.266/-7.440, 1.6e-04, -3.0e-03 L | -7.440/+8.266 | 9.98 | 9.6e-08 |
| 0.6 | 12 | 1e-09 | 6400 | 4 | 6.42e-03, +8.266/-7.439, 7.1e-05, -1.8e-03 L | -7.439/+8.266 | 9.98 | 3.3e-08 |
| 0.6 | 12 | 1e-08 | 400 | 4 | 2.14e-02, +6.107/-5.286, 9.5e-04, -9.9e-03 L | none | 7.37 | 1.0e-06 |
| 0.6 | 12 | 1e-08 | 1600 | 4 | 1.51e-02, +6.109/-5.285, 3.4e-04, -5.0e-03 L | none | 7.37 | 2.6e-07 |
| 0.6 | 12 | 1e-08 | 6400 | 4 | 1.12e-02, +6.110/-5.284, 1.4e-04, -2.8e-03 L | none | 7.38 | 8.5e-08 |
| 0.6 | 12 | 1e-07 | 400 | 4 | 4.04e-02, +4.648/-3.835, 2.5e-03, -1.9e-02 | none | 5.61 | 3.7e-06 |
| 0.6 | 12 | 1e-07 | 1600 | 4 | 2.73e-02, +4.652/-3.831, 8.1e-04, -8.9e-03 L | none | 5.62 | 8.3e-07 |
| 0.6 | 12 | 1e-07 | 6400 | 4 | 1.97e-02, +4.654/-3.830, 3.1e-04, -4.7e-03 L | none | 5.62 | 2.4e-07 |
| 0.6 | 12 | 1e-06 | 400 | 4 | 1.55e-01, -5.185/+4.604, 1.2e+00, -1.7e+00 | none | 6.26 | 4.1e-04 |
| 0.6 | 12 | 1e-06 | 1600 | 5 | 1.31e-01, -5.164/+4.637, 1.3e+00, -1.6e+00 | none | 6.23 | 2.2e-04 |
| 0.6 | 12 | 1e-06 | 6400 | 5 | 8.21e-02, -5.103/+4.687, 1.1e+00, +1.5e-01 | none | 6.16 | 1.3e-04 |
| 0.6 | 12 | 1e-05 | 400 | 5 | 4.93e-01, -3.409/+0.988, 1.5e+00, -3.9e+00 | none | 4.12 | 1.0e-03 |
| 0.6 | 12 | 1e-05 | 1600 | 5 | 4.44e-01, -3.347/+1.018, 1.1e+00, -3.8e+00 | none | 4.04 | 2.9e-04 |
| 0.6 | 12 | 1e-05 | 6400 | 5 | 1.53e-01, -3.285/+1.038, 6.1e-01, -3.7e+00 | none | 3.97 | 7.2e-06 |
| 0.6 | 20 | 1e-10 | 400 | 4 | 8.80e-03, -4.838/+4.332, 1.6e-04, +3.0e-03 L | +4.332/-4.838 | 5.84 | 1.2e-07 |
| 0.6 | 20 | 1e-10 | 1600 | 4 | 6.84e-03, -4.838/+4.332, 7.5e-05, +1.8e-03 L | +4.332/-4.838 | 5.84 | 4.3e-08 |
| 0.6 | 20 | 1e-10 | 6400 | 4 | 5.38e-03, -4.838/+4.332, 3.7e-05, +1.1e-03 L | +4.332/-4.838 | 5.84 | 1.7e-08 |
| 0.6 | 20 | 1e-09 | 400 | 4 | 1.61e-02, -3.369/+2.864, 3.0e-04, +4.6e-03 L | none | 4.07 | 2.8e-07 |
| 0.6 | 20 | 1e-09 | 1600 | 4 | 1.24e-02, -3.369/+2.864, 1.4e-04, +2.8e-03 L | none | 4.07 | 1.0e-07 |
| 0.6 | 20 | 1e-09 | 6400 | 4 | 9.70e-03, -3.369/+2.863, 6.7e-05, +1.7e-03 L | +2.863/-3.369 | 4.07 | 3.8e-08 |
| 0.6 | 20 | 1e-08 | 400 | 4 | 3.02e-02, -2.369/+1.865, 6.1e-04, +7.3e-03 L | none | 2.86 | 7.5e-07 |
| 0.6 | 20 | 1e-08 | 1600 | 4 | 2.29e-02, -2.369/+1.865, 2.7e-04, +4.3e-03 L | none | 2.86 | 2.5e-07 |
| 0.6 | 20 | 1e-08 | 6400 | 4 | 1.77e-02, -2.370/+1.864, 1.3e-04, +2.6e-03 L | none | 2.86 | 8.9e-08 |
| 0.6 | 20 | 1e-07 | 400 | 4 | 5.76e-02, -1.694/+1.193, 1.3e-03, +1.2e-02 | none | 2.05 | 2.2e-06 |
| 0.6 | 20 | 1e-07 | 1600 | 4 | 4.30e-02, -1.695/+1.192, 5.5e-04, +6.8e-03 L | none | 2.05 | 6.8e-07 |
| 0.6 | 20 | 1e-07 | 6400 | 4 | 3.27e-02, -1.696/+1.191, 2.4e-04, +4.0e-03 L | none | 2.05 | 2.3e-07 |
| 0.6 | 20 | 1e-06 | 400 | 4 | 3.81e+00, +0.828/-0.227, 1.8e+00, -2.7e-01 | none | 1.00 | 1.2e-05 |
| 0.6 | 20 | 1e-06 | 1600 | 4 | 3.88e+00, +0.828/-0.225, 1.8e+00, -2.7e-01 | none | 1.00 | 3.2e-06 |
| 0.6 | 20 | 1e-06 | 6400 | 4 | 3.93e+00, +0.828/-0.224, 1.9e+00, -2.7e-01 | none | 1.00 | 8.7e-07 |
| 0.6 | 20 | 1e-05 | 400 | 4 | 3.80e+00, +0.826/-0.285, 2.5e+00, -3.5e-01 | none | 1.00 | 4.0e-05 |
| 0.6 | 20 | 1e-05 | 1600 | 4 | 3.86e+00, +0.825/-0.280, 2.5e+00, -3.4e-01 | none | 1.00 | 1.1e-05 |
| 0.6 | 20 | 1e-05 | 6400 | 4 | 3.91e+00, +0.825/-0.278, 2.5e+00, -3.4e-01 | none | 1.00 | 3.2e-06 |
| 0.6 | 27 | 1e-10 | 400 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-10 | 1600 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-10 | 6400 | 4 | 1.91e-04, -13.138/+12.578, 2.7e-08, +9.5e-06 L | -13.138/+12.578 | 23.47 | 3.4e-12 |
| 0.6 | 27 | 1e-09 | 400 | 4 | 3.88e-03, +8.812/-8.330, 1.0e-04, -2.3e-03 L | -8.330/+8.812 | 15.74 | 7.7e-08 |
| 0.6 | 27 | 1e-09 | 1600 | 4 | 2.93e-03, +8.813/-8.330, 4.5e-05, -1.3e-03 L | -8.330/+8.813 | 15.74 | 2.6e-08 |
| 0.6 | 27 | 1e-09 | 6400 | 4 | 2.27e-03, +8.813/-8.330, 2.1e-05, -7.9e-04 L | -8.330/+8.813 | 15.74 | 9.3e-09 |
| 0.6 | 27 | 1e-08 | 400 | 4 | 7.52e-03, +6.043/-5.562, 2.3e-04, -3.9e-03 L | -5.562/+6.043 | 10.79 | 2.2e-07 |
| 0.6 | 27 | 1e-08 | 1600 | 4 | 5.54e-03, +6.044/-5.561, 9.2e-05, -2.1e-03 L | -5.561/+6.044 | 10.79 | 6.8e-08 |
| 0.6 | 27 | 1e-08 | 6400 | 4 | 4.21e-03, +6.044/-5.561, 4.1e-05, -1.2e-03 L | -5.561/+6.044 | 10.80 | 2.3e-08 |
| 0.6 | 27 | 1e-07 | 400 | 4 | 1.53e-02, +4.152/-3.672, 5.7e-04, -7.1e-03 L | none | 7.42 | 7.8e-07 |
| 0.6 | 27 | 1e-07 | 1600 | 4 | 1.09e-02, +4.153/-3.671, 2.1e-04, -3.7e-03 L | none | 7.42 | 2.1e-07 |
| 0.6 | 27 | 1e-07 | 6400 | 4 | 8.07e-03, +4.153/-3.671, 8.7e-05, -2.0e-03 L | -3.671/+4.153 | 7.42 | 6.5e-08 |
| 0.6 | 27 | 1e-06 | 400 | 5 | 2.98e-01, +3.872/-2.075, 2.5e-01, +4.2e+00 | none | 6.92 | 3.8e-05 |
| 0.6 | 27 | 1e-06 | 1600 | 5 | 2.23e-01, +3.898/-2.061, 2.0e-01, +4.3e+00 | none | 6.96 | 8.9e-06 |
| 0.6 | 27 | 1e-06 | 6400 | 5 | 1.79e-01, +3.911/-2.055, 1.6e-01, +4.3e+00 | none | 6.99 | 2.7e-06 |
| 0.6 | 27 | 1e-05 | 400 | 4 | 1.01e-01, +1.745/-1.282, 7.7e-03, -4.3e-02 | none | 3.12 | 1.3e-05 |
| 0.6 | 27 | 1e-05 | 1600 | 4 | 8.25e-02, +1.748/-1.279, 3.8e-03, -3.0e-02 | none | 3.12 | 4.1e-06 |
| 0.6 | 27 | 1e-05 | 6400 | 4 | 6.91e-02, +1.751/-1.278, 1.9e-03, -2.2e-02 | none | 3.13 | 1.5e-06 |
| 0.6 | 29 | 1e-10 | 400 | 4 | 5.07e-03, -7.381/+6.789, 7.0e-05, +1.8e-03 L | +6.789/-7.381 | 8.63 | 6.8e-08 |
| 0.6 | 29 | 1e-10 | 1600 | 4 | 3.90e-03, -7.381/+6.788, 3.2e-05, +1.0e-03 L | +6.788/-7.381 | 8.63 | 2.4e-08 |
| 0.6 | 29 | 1e-10 | 6400 | 4 | 3.05e-03, -7.381/+6.788, 1.6e-05, +6.5e-04 L | +6.788/-7.381 | 8.63 | 9.2e-09 |
| 0.6 | 29 | 1e-09 | 400 | 4 | 9.53e-03, -5.034/+4.442, 1.4e-04, +2.7e-03 L | +4.442/-5.034 | 5.89 | 1.7e-07 |
| 0.6 | 29 | 1e-09 | 1600 | 4 | 7.23e-03, -5.034/+4.442, 6.1e-05, +1.6e-03 L | +4.442/-5.034 | 5.89 | 5.7e-08 |
| 0.6 | 29 | 1e-09 | 6400 | 4 | 5.61e-03, -5.035/+4.442, 2.9e-05, +9.7e-04 L | +4.442/-5.035 | 5.89 | 2.1e-08 |
| 0.6 | 29 | 1e-08 | 400 | 4 | 1.89e-02, -3.389/+2.798, 2.9e-04, +4.6e-03 L | none | 3.96 | 4.7e-07 |
| 0.6 | 29 | 1e-08 | 1600 | 4 | 1.40e-02, -3.390/+2.798, 1.2e-04, +2.5e-03 L | none | 3.96 | 1.5e-07 |
| 0.6 | 29 | 1e-08 | 6400 | 4 | 1.07e-02, -3.390/+2.797, 5.5e-05, +1.5e-03 L | none | 3.96 | 5.0e-08 |
| 0.6 | 29 | 1e-07 | 400 | 4 | 4.07e-02, -2.248/+1.660, 7.0e-04, +8.1e-03 L | none | 2.63 | 1.6e-06 |
| 0.6 | 29 | 1e-07 | 1600 | 4 | 2.92e-02, -2.250/+1.659, 2.7e-04, +4.3e-03 L | none | 2.63 | 4.3e-07 |
| 0.6 | 29 | 1e-07 | 6400 | 4 | 2.17e-02, -2.250/+1.658, 1.1e-04, +2.4e-03 L | none | 2.63 | 1.4e-07 |
| 0.6 | 29 | 1e-06 | 400 | 4 | 8.57e-02, -1.578/+0.993, 1.8e-03, +1.5e-02 | none | 1.85 | 5.8e-06 |
| 0.6 | 29 | 1e-06 | 1600 | 4 | 5.98e-02, -1.580/+0.991, 6.2e-04, +7.4e-03 L | none | 1.85 | 1.5e-06 |
| 0.6 | 29 | 1e-06 | 6400 | 4 | 4.33e-02, -1.581/+0.990, 2.4e-04, +3.9e-03 L | none | 1.85 | 4.1e-07 |
| 0.6 | 29 | 1e-05 | 400 | 4 | 6.32e+00, +0.854/-0.399, 6.1e+00, +1.8e+00 | none | 1.00 | 2.3e-05 |
| 0.6 | 29 | 1e-05 | 1600 | 4 | 6.36e+00, +0.854/-0.395, 6.1e+00, +1.8e+00 | none | 1.00 | 6.0e-06 |
| 0.6 | 29 | 1e-05 | 6400 | 4 | 6.39e+00, +0.854/-0.394, 6.1e+00, +1.8e+00 | none | 1.00 | 1.7e-06 |

(ii)(a)'s width split on the landed pairs (L, 72 rows): the observed σ_b − σ_a over (w_b − w_a)Δc²/(6Wσ) 0.9828–1.0458 (median 0.9993), at |w_a − w_b||Δc|/(|W|σ) = 8.74e-03–2.84e-01; its largest deviation 4.6e-02 at the parameter 2.76e-01.

The three phases (`phases`): the dipole block's entry τ_d (the path's) and the next kernel's formation τ₃ in the loop; over the grid times between them the slopes of v³ (v the dipole pair's half-difference) against 3μ₀KΛ(ε) — the block's law in its fast time Λ(ε)τ, K and μ₀ from the flow's last output — (median [range]); the landed v (at τ = 400) against the law's (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}; landed at τ = 400, 1600, 6400; over 1600 → 6400 the exponents in τ of the pair's distance and of the residual:

| Δ ≥ | target | ε | τ_d | τ₃ | d(v³)/dτ / 3μ₀KΛ | v at 400 | the law's v | landed | exponents: distance, residual |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 9 | 1e-05 | 17.75 | 22 | 1.360 [1.360, 1.360] | 2.6 | 1.936 | no, no, no | -0.173, -0.794 |
| 0.4 | 9 | 1e-06 | 17.75 | 28 | 1.108 [1.077, 1.223] | 4.244 | 3.81 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-07 | 17.75 | 30 | 1.046 [1.032, 1.123] | 6.414 | 5.935 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-08 | 17.75 | 30 | 1.022 [1.015, 1.063] | 9.467 | 8.711 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-09 | 17.75 | 30 | 1.011 [1.007, 1.031] | 13.92 | 12.79 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 9 | 1e-10 | 17.75 | 30 | 1.005 [1.003, 1.015] | 20.43 | 18.77 | yes, yes, yes | -0.000, -0.000 |
| 0.4 | 18 | 1e-05 | 25 | 26 | — | 2.829 | 1.057 | no, no, no | -0.571, -3.061 |
| 0.4 | 18 | 1e-06 | 25 | 26 | — | 4.823 | 1.552 | no, no, yes | -0.322, -1.211 |
| 0.4 | 18 | 1e-07 | 25 | 28 | — | 8.175 | 3.288 | no, yes, yes | -0.230, -0.880 |
| 0.4 | 18 | 1e-08 | 25 | 34 | 1.013 [1.010, 1.021] | 5.339 | 6.962 | no, no, yes | -0.190, -0.723 |
| 0.4 | 18 | 1e-09 | 25 | 36 | 1.005 [1.004, 1.010] | 7.773 | 10.93 | no, yes, yes | -0.191, -0.718 |
| 0.4 | 18 | 1e-10 | 25 | 36 | 1.003 [1.002, 1.005] | 11.38 | 16.04 | yes, yes, yes | -0.189, -0.713 |
| 0.6 | 12 | 1e-05 | 24.35 | 26 | — | 2.199 | 0.6147 | no, no, no | -0.771, -2.659 |
| 0.6 | 12 | 1e-06 | 24.35 | 26 | — | 4.894 | 0.9022 | no, no, no | -0.336, -0.374 |
| 0.6 | 12 | 1e-07 | 24.35 | 57.4349 | 1.768 [1.569, 2.688] | 4.241 | 3.594 | no, yes, yes | -0.235, -0.901 |
| 0.6 | 12 | 1e-08 | 24.35 | 57.4349 | 1.427 [1.310, 1.988] | 5.697 | 5.276 | yes, yes, yes | -0.212, -0.821 |
| 0.6 | 12 | 1e-09 | 24.35 | 57.4349 | 1.228 [1.162, 1.558] | 7.853 | 7.744 | yes, yes, yes | -0.196, -0.762 |
| 0.6 | 12 | 1e-10 | 24.35 | 65.9754 | 1.112 [1.073, 1.304] | 11.06 | 12.27 | yes, yes, yes | -0.185, -0.726 |
| 0.6 | 20 | 1e-05 | 32.22 | 34 | — | 0.5552 | 0.5173 | no, no, no | +0.009, -0.886 |
| 0.6 | 20 | 1e-06 | 32.22 | 34 | — | 0.5273 | 0.7593 | no, no, no | +0.010, -0.937 |
| 0.6 | 20 | 1e-07 | 32.22 | 36 | — | 1.443 | 1.433 | no, yes, yes | -0.198, -0.791 |
| 0.6 | 20 | 1e-08 | 32.22 | 38 | 0.984 [0.984, 0.984] | 2.117 | 2.424 | yes, yes, yes | -0.186, -0.741 |
| 0.6 | 20 | 1e-09 | 32.22 | 38 | 0.993 [0.993, 0.993] | 3.116 | 3.557 | yes, yes, yes | -0.178, -0.708 |
| 0.6 | 20 | 1e-10 | 32.22 | 38 | 0.997 [0.997, 0.997] | 4.585 | 5.222 | yes, yes, yes | -0.173, -0.689 |
| 0.6 | 27 | 1e-05 | 31.54 | 32 | — | 1.513 | 0.5708 | no, no, no | -0.128, -0.740 |
| 0.6 | 27 | 1e-06 | 31.54 | 32 | — | 2.974 | 0.8378 | no, no, no | -0.159, -0.846 |
| 0.6 | 27 | 1e-07 | 31.54 | 36 | 1.064 [1.064, 1.064] | 3.912 | 2.628 | yes, yes, yes | -0.219, -0.843 |
| 0.6 | 27 | 1e-08 | 31.54 | 38 | 1.026 [1.018, 1.033] | 5.802 | 4.364 | yes, yes, yes | -0.198, -0.777 |
| 0.6 | 27 | 1e-09 | 31.54 | 46 | 1.006 [1.004, 1.016] | 8.571 | 8.38 | yes, yes, yes | -0.185, -0.730 |
| 0.6 | 27 | 1e-10 | 31.54 | 48 | 1.002 [1.002, 1.008] | 12.86 | 12.84 | yes, yes, yes | -0.000, -0.000 |
| 0.6 | 29 | 1e-05 | 26.75 | 28 | — | 0.6267 | 1.038 | no, no, no | +0.003, -0.925 |
| 0.6 | 29 | 1e-06 | 26.75 | 30 | — | 1.285 | 2.096 | no, yes, yes | -0.233, -0.906 |
| 0.6 | 29 | 1e-07 | 26.75 | 30 | — | 1.954 | 3.076 | yes, yes, yes | -0.213, -0.833 |
| 0.6 | 29 | 1e-08 | 26.75 | 30 | — | 3.094 | 4.515 | yes, yes, yes | -0.194, -0.763 |
| 0.6 | 29 | 1e-09 | 26.75 | 30 | — | 4.738 | 6.628 | yes, yes, yes | -0.183, -0.720 |
| 0.6 | 29 | 1e-10 | 26.75 | 30 | — | 7.085 | 9.728 | yes, yes, yes | -0.177, -0.698 |

The landed v against the mass: v at τ = 400 and its ratio per decade of ε (10^{1/6} = 1.4678 if v ∝ ε^{−1/6}):

| Δ ≥ | target | v at ε = 1e-05, 1e-06, 1e-07, 1e-08, 1e-09, 1e-10 | ratios per decade |
|---|---|---|---|
| 0.4 | 9 | 2.6, 4.244, 6.414, 9.467, 13.92, 20.43 | 1.6324, 1.5112, 1.4760, 1.4699, 1.4682 |
| 0.4 | 18 | 2.829, 4.823, 8.175, 5.339, 7.773, 11.38 | 1.7052, 1.6948, 0.6531, 1.4560, 1.4643 |
| 0.6 | 12 | 2.199, 4.894, 4.241, 5.697, 7.853, 11.06 | 2.2262, 0.8666, 1.3431, 1.3785, 1.4078 |
| 0.6 | 20 | 0.5552, 0.5273, 1.443, 2.117, 3.116, 4.585 | 0.9496, 2.7370, 1.4669, 1.4721, 1.4713 |
| 0.6 | 27 | 1.513, 2.974, 3.912, 5.802, 8.571, 12.86 | 1.9648, 1.3156, 1.4832, 1.4772, 1.5002 |
| 0.6 | 29 | 0.6267, 1.285, 1.954, 3.094, 4.738, 7.085 | 2.0507, 1.5203, 1.5832, 1.5316, 1.4953 |

Over the six targets: landed at τ = 6400 in 0 of 6 at ε = 1e-05, 3 of 6 at ε = 1e-06, 6 of 6 at ε = 1e-07, 6 of 6 at ε = 1e-08, 6 of 6 at ε = 1e-09, 6 of 6 at ε = 1e-10; at ε = 1e-10 the dipole phase's slope over the law 0.997–1.112 (the cells' medians, 5 cells with a phase), the landed v over the law's 0.71–1.09, v's ratio over the last decade 1.4078–1.5002; on the 25 cells landed at both 1600 and 6400, 6 frozen (9 at 1e-06, 9 at 1e-07, 9 at 1e-08, 9 at 1e-09, 9 at 1e-10, 27 at 1e-10), the others collapsing with the exponents -0.235 … -0.173 (distance) and -0.906 … -0.689 (residual).

## H. Two features of one sign, close

Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 — 0.16 apart, 1.33 times the kernels' starting width σ₀ = 0.12 — unless stated, the moving-width path (as C's: its tolerance 1e-09, its budget 1e+06): the outcome (as C's), the most kernels of the other sign in a block along the path and the most negative coefficient, the end's formed kernels of the other sign, the smallest formed width, over all the end's active kernels those narrower than the samples' spacing and those of negative width (the widths are not constrained), the largest |w|, the residual:

| target | outcome | negative kernels along the path (most negative w) | at the end | smallest width | sub-sample, negative widths (all active) | max\|w\| | residual |
|---|---|---|---|---|---|---|---|
| σ* 0.04, equal | exact, null cluster | 10 (-0.427) | 2 | 4.00e-02 | 0, 0 | 1 | 1.3e-13 |
| σ* 0.06, equal | cancellation | 14 (-0.610) | 7 | 4.96e-02 | 6, 0 | 3.24 | 4.3e-10 |
| σ* 0.06, A₂ 0.99 | cancellation | 13 (-0.573) | 6 | 5.19e-02 | 5, 1 | 2.93 | 9.6e-10 |
| σ* 0.06, A₂ 0.95 | cancellation | 14 (-0.587) | 7 | 4.92e-02 | 4, 1 | 3.14 | 6.0e-10 |
| σ* 0.06, A₂ 0.9 | cancellation | 12 (-0.555) | 7 | 5.12e-02 | 3, 1 | 2.84 | 8.1e-10 |
| σ* 0.06, c₂ +0.01 | cancellation | 11 (-0.561) | 7 | 1.75e-04 | 5, 1 | 2.85 | 4.4e-10 |
| σ* 0.06, σ₂ 0.065 | cancellation | 14 (-0.603) | 7 | 4.94e-02 | 5, 2 | 3.16 | 6.1e-10 |
| σ* 0.04, A₂ 0.6 | exact | 0 (+0.600) | 0 | 4.00e-02 | 0, 0 | 1 | 5.0e-17 |
| σ* 0.06, A₂ 0.6 | exact | 0 (+0.600) | 0 | 6.00e-02 | 0, 0 | 1 | 2.5e-16 |

### H2. The σ* = 0.06 paths to a tolerance of 1e-13

The same paths with the tolerance lowered (the 1e-09 ends were the tolerance's): the outcome, the entries and the last one's τ, the residuals of the last blocks, the formed kernels by the feature nearest each (their count and mass over the feature's, Σw_i|σ_i|/(Aσ*)), and H's columns:

| target | outcome | entries (last τ) | the last blocks' residuals | formed per feature (count: mass) | negative along the path (most negative w) | at the end | smallest width | sub-sample, negative widths | max\|w\| |
|---|---|---|---|---|---|---|---|---|---|
| σ* 0.06, equal | cancellation | 24 (1.10e+15) | 2.83e-10, 2.32e-10, 2.58e-10, 1.66e-10, 1.12e-10, 8.39e-11 | 1: 1.0000, 8: 1.0000 | 14 (-0.610) | 7 | 5.08e-02 | 6, 0 | 3.26 |
| σ* 0.06, A₂ 0.99 | stopped by the budget | 25 (7.36e+13) | 1.15e-10, 7.94e-11, 5.88e-11, 4.66e-11, 3.73e-11, 2.89e-11 | 7: 1.0000, 1: 1.0000 | 13 (-0.573) | 6 | 5.43e-02 | 5, 1 | 2.96 |
| σ* 0.06, A₂ 0.95 | cancellation | 27 (2.29e+16) | 1.01e-10, 7.49e-11, 6.48e-11, 5.37e-11, 4.97e-11, 4.13e-11 | 8: 1.0000, 1: 1.0000 | 14 (-0.587) | 7 | 5.14e-02 | 4, 1 | 3.18 |
| σ* 0.06, A₂ 0.9 | stopped by the budget | 20 (1.16e+12) | 5.67e-09, 8.13e-10, 5.16e-10, 4.00e-10, 2.92e-10, 2.02e-10 | 7: 1.0002, 2: 0.9998 | 12 (-0.555) | 7 | 5.22e-02 | 3, 1 | 2.85 |
| σ* 0.06, c₂ +0.01 | cancellation | 32 (2.34e+67) | 1.71e-11, 1.70e-11, 1.77e-11, 4.21e-11, 2.60e-11, 1.49e-11 | 7: 1.0000, 2: 0.9999 | 11 (-0.561) | 7 | 1.75e-04 | 6, 2 | 2.88 |
| σ* 0.06, σ₂ 0.065 | cancellation | 24 (2.83e+14) | 1.92e-10, 1.31e-10, 9.62e-11, 7.51e-11, 6.68e-11, 4.99e-11 | 8: 1.0000, 1: 1.0000 | 14 (-0.603) | 7 | 5.15e-02 | 5, 2 | 3.2 |
| σ* 0.06, A₂ 0.6 | exact | 2 (3.09e+02) | 2.28e-01, 4.88e-02, 2.48e-16 | 1: 1.0000, 1: 1.0000 | 0 (+0.600) | 0 | 6.00e-02 | 0, 0 | 1 |

### H3. The separation scan

σ* = 0.06, equal amplitudes 1, at ±c (2c in units of σ₀ in brackets; H's ±0.08 is 1.33σ₀), H's columns:

| separation | outcome | negative kernels along the path (most negative w) | at the end | smallest width | sub-sample, negative widths (all active) | max\|w\| | residual |
|---|---|---|---|---|---|---|---|
| ±0.05 (0.83σ₀) | stopped by the budget | 10 (-0.240) | 10 | 3.23e-02 | 0, 0 | 1.42 | 6.5e-06 |
| ±0.08 (1.33σ₀) | cancellation | 14 (-0.610) | 7 | 4.96e-02 | 6, 0 | 3.24 | 4.3e-10 |
| ±0.1 (1.67σ₀) | exact | 0 (+0.000) | 0 | 6.00e-02 | 0, 0 | 1 | 3.3e-10 |
| ±0.12 (2.00σ₀) | exact | 0 (-0.000) | 0 | 6.00e-02 | 0, 0 | 1 | 2.5e-12 |
| ±0.15 (2.50σ₀) | coinciding pair | 0 (+0.305) | 0 | 6.00e-02 | 0, 0 | 1 | 5.8e-16 |
| ±0.2 (3.33σ₀) | coinciding pair | 0 (+0.231) | 0 | 6.00e-02 | 0, 0 | 0.631 | 2.3e-15 |

### HL. The real loop on H's equal targets

The free law from α = 0 at the masses ε (`loop_nd`, rtol 10⁻⁹) at the times τ: the formed kernels (|w| above 1% of the smallest amplitude), of the other sign, the smallest formed width, the formed kernels centred outside the samples' window, the largest |w|, the residual:

| target | ε | τ | formed | of the other sign (most negative w) | smallest width | outside the window | max\|w\| | residual |
|---|---|---|---|---|---|---|---|---|
| σ* 0.04, equal | 1e-06 | 1600 | 9 | 6 (-0.027) | 0.036 | 0 | 0.951 | 7.8e-04 |
| σ* 0.04, equal | 1e-06 | 6400 | 2 | 0 | 0.040 | 0 | 0.995 | 6.7e-05 |
| σ* 0.04, equal | 1e-08 | 1600 | 9 | 8 (-0.408) | 0.024 | 0 | 0.716 | 7.2e-02 |
| σ* 0.04, equal | 1e-08 | 6400 | 2 | 0 | 0.040 | 0 | 0.998 | 2.0e-05 |
| σ* 0.06, equal | 1e-06 | 1600 | 7 | 6 (-0.389) | 0.037 | 0 | 0.979 | 2.7e-02 |
| σ* 0.06, equal | 1e-06 | 6400 | 9 | 6 (-0.112) | 0.047 | 0 | 0.572 | 1.4e-04 |
| σ* 0.06, equal | 1e-08 | 1600 | 7 | 6 (-0.389) | 0.037 | 0 | 0.979 | 2.7e-02 |
| σ* 0.06, equal | 1e-08 | 6400 | 13 | 12 (-0.456) | 0.024 | 0 | 0.958 | 2.1e-02 |

Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.
