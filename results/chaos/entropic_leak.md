# Proposition D‴.20: the entropic OT leak — a soft mass constraint with a critical excess (`exp_entropic_leak.py`)

The three-bump instance with the matching prior (masses q = 1/3) at σ₀ = 0.15σ_crit = 0.00707; the entropic plan by log-domain Sinkhorn (tolerance 10⁻¹³).

## 1. The critical excess: the separation margin of the over-weighted class of Θ*(w_δ), w_δ = q + δ(1.5e₁ − 0.5), against ε

| δ | w_δ | smallest eigenvalue of the separation block of class 1 at ε = 1, 4, 16, 64, 256 | ε_c (bisection) | excess fraction p sent away at ε_c | rank-one prediction 2σ₀p(1 − p)Δ̂ᵀ(σ₀I − R^G)⁻¹Δ̂ at ε_c |
|---|---|---|---|---|---|
| 0.02 | (0.3533, 0.3233, 0.3233) | -0.02337, -0.00154, +0.00386, +0.00729, +0.00860 | 5.107 | 0.0566 | 6.218 |
| 0.05 | (0.3833, 0.3083, 0.3083) | -0.05746, -0.01103, +0.00094, +0.00627, +0.00798 | 13.508 | 0.1515 | 16.584 |
| 0.1 | (0.4333, 0.2833, 0.2833) | -0.09392, -0.02117, -0.00243, +0.00417, +0.00580 | 23.105 | 0.3345 | 30.275 |
| 0.2 | (0.5333, 0.2333, 0.2333) | -0.17299, -0.03114, -0.04808, -0.36601, +0.00123 | 136.138 | 0.6233 | 88.091 |

Below ε_c the over-weighted collapsed class is a saddle — it splits, and the split-off particles carry the excess to the under-weighted atom (the hard OT leak's transit); above ε_c it is stable and the excess is absorbed by a bias of the class atom toward its target T₁ = Σ_j(π_{1j}/w₁)ϑ⁰_j.

## 2. Dynamics: m = 15, 60, 240 against M = 510 from the warm clouds of `exp_ot_leak.py`, 600 slow units, 4 seeds

| ε | m | cluster weights at the end: max |w − q|_∞ over seeds (start) | fit error at the end (mean over seeds) | W²(μ_m, ν) at the end (mean) | ‖f̂_m − f̂_M‖² on the last 150 slow units (mean) |
|---|---|---|---|---|---|
| 1 | 15 | 0.0000 (0.2000) | 3.036e-04 | 2.694e-03 | 5.573e-07 |
| 1 | 60 | 0.0000 (0.0833) | 1.662e-04 | 9.936e-04 | 1.648e-07 |
| 1 | 240 | 0.0000 (0.0542) | 1.908e-04 | 7.602e-04 | 9.973e-08 |
| 1 | 510 (reference) | 0.0000 (0.0412) | 2.517e-04 | 3.163e-03 | — |
| 4 | 15 | 0.0000 (0.2000) | 1.469e-03 | 1.989e-02 | 7.054e-05 |
| 4 | 60 | 0.0167 (0.0833) | 8.799e-03 | 4.202e-01 | 7.681e-05 |
| 4 | 240 | 0.0292 (0.0542) | 9.179e-03 | 6.741e-01 | 4.546e-05 |
| 4 | 510 (reference) | 0.0275 (0.0412) | 9.603e-03 | 7.383e-01 | — |
| 16 | 15 | 0.1333 (0.2000) | 6.936e-02 | 3.979e+00 | 1.546e-03 |
| 16 | 60 | 0.0833 (0.0833) | 4.066e-02 | 2.506e+00 | 2.219e-04 |
| 16 | 240 | 0.0542 (0.0542) | 3.635e-02 | 1.211e+00 | 3.531e-05 |
| 16 | 510 (reference) | 0.0412 (0.0412) | 3.310e-02 | 1.007e+00 | — |
| 64 | 15 | 0.2000 (0.2000) | 1.412e-01 | 5.175e+00 | 1.929e-03 |
| 64 | 60 | 0.0833 (0.0833) | 1.245e-01 | 2.876e+00 | 3.928e-04 |
| 64 | 240 | 0.0542 (0.0542) | 1.256e-01 | 1.621e+00 | 3.995e-05 |
| 64 | 510 (reference) | 0.0412 (0.0412) | 1.231e-01 | 1.414e+00 | — |

Numbers are produced by `experiments/exp_entropic_leak.py`; nothing is edited by hand.
