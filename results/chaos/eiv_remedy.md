# Proposition D‴.31: the adaptive prior is a washout leak; the desired regressor (`exp_eiv_remedy.py`)

The instance of Proposition D‴.30 (σ₀ = 0.00707, x_m = 2 sin t, the lift of the exact representation of q = 1/3 with the prior at its atoms) under the measurement y = x + ν, with an Ornstein–Uhlenbeck noise of rate B. MR: the measured regressor and compensation (Proposition D‴.30); DR: the regressor at the reference, the compensation f̂(y); DR*: both at the reference (the desired compensation adaptation law).

## A. The second-order force and its stationary responses (no dynamics)

Per unit (Eν²)², function level: the stationary bias under the fixed prior and under the adaptive prior (every η > 0; = FᵀH⁻¹F, the squared norm of the smallest target perturbation exerting F), the least-squares bound ‖g‖²_ρ (g = k_eE[e]; DR, DR*), and the projected mean tracking error at the adaptive prior's stationary point, ‖P(g − Ĵȳ_∞)‖²_ρ = F₁ᵀH⁻¹F₁ (F₁ the regressor's errors in variables, zero for DR and DR*). The last column is the noise variance at which the adaptive prior's bias equals its floor of the parameters (N₀/2k_e)tr Q, N₀ = 2Eν²/B; the bias dominates above it.

| k_e | B | law | fixed prior | adaptive prior | ‖g‖²_ρ | projected mean error (adaptive) | Eν² at bias = floor (std) |
|---|---|---|---|---|---|---|---|
| 25 | 100 | MR | 4.9106e+02 | 2.3969e+03 | — | 2.394e+03 | 1.47e-06 (0.0012) |
| 25 | 100 | DR | 2.1060e-01 | 4.7053e-01 | 4.7057e-01 | 0.000e+00 | 7.47e-03 (0.0864) |
| 25 | 100 | DR* | 2.3601e-02 | 5.2621e-02 | 5.2663e-02 | 0.000e+00 | 6.68e-02 (0.2584) |
| 25 | 1000 | MR | 7.2955e+02 | 3.5637e+03 | — | 3.557e+03 | 9.86e-08 (0.0003) |
| 25 | 1000 | DR | 5.2899e-01 | 1.1821e+00 | 1.1821e+00 | 0.000e+00 | 2.97e-04 (0.0172) |
| 25 | 1000 | DR* | 3.5056e-04 | 7.8187e-04 | 7.8240e-04 | 0.000e+00 | 4.49e-01 (0.6703) |
| 50 | 100 | MR | 1.3658e+03 | 6.6485e+03 | — | 6.646e+03 | 2.64e-07 (0.0005) |
| 50 | 100 | DR | 6.5041e-02 | 1.4529e-01 | 1.4531e-01 | 0.000e+00 | 1.21e-02 (0.1100) |
| 50 | 100 | DR* | 6.5109e-02 | 1.4542e-01 | 1.4545e-01 | 0.000e+00 | 1.21e-02 (0.1099) |
| 50 | 1000 | MR | 2.7847e+03 | 1.3569e+04 | — | 1.356e+04 | 1.29e-08 (0.0001) |
| 50 | 1000 | DR | 4.7858e-01 | 1.0695e+00 | 1.0695e+00 | 0.000e+00 | 1.64e-04 (0.0128) |
| 50 | 1000 | DR* | 1.3281e-03 | 2.9666e-03 | 2.9671e-03 | 0.000e+00 | 5.92e-02 (0.2433) |
| 100 | 100 | MR | 3.0756e+03 | 1.4949e+04 | — | 1.495e+04 | 5.88e-08 (0.0002) |
| 100 | 100 | DR | 3.1075e-05 | 4.5092e-05 | 5.8406e-05 | 0.000e+00 | 1.95e+01 (4.4135) |
| 100 | 100 | DR* | 1.4623e-01 | 3.2675e-01 | 3.2677e-01 | 0.000e+00 | 2.69e-03 (0.0518) |
| 100 | 1000 | MR | 1.0158e+04 | 4.9429e+04 | — | 4.941e+04 | 1.78e-09 (0.0000) |
| 100 | 1000 | DR | 3.9137e-01 | 8.7459e-01 | 8.7459e-01 | 0.000e+00 | 1.00e-04 (0.0100) |
| 100 | 1000 | DR* | 4.8334e-03 | 1.0800e-02 | 1.0801e-02 | 0.000e+00 | 8.13e-03 (0.0902) |
| 200 | 100 | MR | 5.4701e+03 | 2.6568e+04 | — | 2.657e+04 | 1.65e-08 (0.0001) |
| 200 | 100 | DR | 6.4958e-02 | 1.4516e-01 | 1.4517e-01 | 0.000e+00 | 3.03e-03 (0.0550) |
| 200 | 100 | DR* | 2.5986e-01 | 5.8069e-01 | 5.8070e-01 | 0.000e+00 | 7.56e-04 (0.0275) |
| 200 | 1000 | MR | 3.4156e+04 | 1.6610e+05 | — | 1.661e+05 | 2.64e-10 (0.0000) |
| 200 | 1000 | DR | 2.5984e-01 | 5.8066e-01 | 5.8066e-01 | 0.000e+00 | 7.56e-05 (0.0087) |
| 200 | 1000 | DR* | 1.6240e-02 | 3.6292e-02 | 3.6292e-02 | 0.000e+00 | 1.21e-03 (0.0348) |

The adaptive prior's bias, MR over DR: 5.09e+03 (k_e = 25, B = 100), 3.01e+03 (k_e = 25, B = 1000), 4.58e+04 (k_e = 50, B = 100), 1.27e+04 (k_e = 50, B = 1000), 3.32e+08 (k_e = 100, B = 100), 5.65e+04 (k_e = 100, B = 1000), 1.83e+05 (k_e = 200, B = 100), 2.86e+05 (k_e = 200, B = 1000). MR's grows like k_e² where B ≫ k_e; DR's and DR*'s stay below ¼‖∂²_xf‖²_ρ(1 + O(1/k_e)) = 1.3065e+00 for every k_e and B.

## B. The stationary identity HR + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a)

η = σ₀; R and B from the averaged collective block at zero frequency (k_e = 25, B = 100). h₁ = 3.0587e-04, h_max = 1.9521e-01. The errors-in-variables biases are function-level per unit (Eν²)²; the anchor's biases are function-level for a unit misplacement along the weakest and the strongest mode.

| prior | S₀ (zero-frequency leak) | weakest mode: response r₁ | bias factor b₁ | h₁r₁ + b₁ | EIV bias, MR | EIV bias, DR | anchor bias, weakest | anchor bias, strongest | max|HR + B − I| |
|---|---|---|---|---|---|---|---|---|---|
| adaptive, κ = 0 | 0.0000e+00 | 3.2694e+03 | 0.0000 | 1.000000 | 2.3969e+03 | 4.7053e-01 | 0.000e+00 | 0.000e+00 | 2.4e-15 |
| κ = 0.1η | 6.4313e-04 | 1.0537e+03 | 0.6777 | 1.000000 | 1.5667e+03 | 4.0005e-01 | 1.405e-04 | 2.105e-06 | 1.0e-15 |
| κ = 1η | 3.5372e-03 | 2.6021e+02 | 0.9204 | 1.000000 | 6.8172e+02 | 2.8736e-01 | 2.591e-04 | 6.183e-05 | 5.6e-16 |
| κ = 10η | 6.4313e-03 | 1.4843e+02 | 0.9546 | 1.000000 | 5.1129e+02 | 2.2184e-01 | 2.787e-04 | 1.986e-04 | 6.7e-16 |
| fixed (η = 0) | 7.0744e-03 | 1.3550e+02 | 0.9586 | 1.000000 | 4.9106e+02 | 2.1060e-01 | 2.810e-04 | 2.388e-04 | 4.3e-16 |

Along the family the zero-frequency leak S₀ = σ₀κ/(η + κ) trades the anchor's bias for the response to persistent forces, mode by mode, at h_jr_j + b_j = 1: no prior removes the bias and bounds the response at once. The desired regressor changes the force instead: its bias is three to four orders of magnitude below MR's under every prior of the family.

## C. The real loop: MR, DR and DR* driven by the same noise paths (k_e = 25, B = 100, std 0.01, η = ∞, RK4 dt = 0.001)

E‖f̂ − f‖²_ρ at the records (mean ± standard error over 12 runs) against the prediction (covariance part + mean part); MR and DR share the covariance part (their linearisations coincide).

| t | MR: Monte Carlo | MR: predicted | DR: Monte Carlo | DR: predicted | DR*: Monte Carlo | DR*: predicted | covariance part (MR, DR) | mean part: MR | DR | DR* |
|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 1.2049e-07 ± 3.2e-08 | 1.4519e-07 | 8.2906e-08 ± 2.0e-08 | 7.5793e-08 | 8.5242e-08 ± 2.1e-08 | 7.5144e-08 | 7.5790e-08 | 6.9404e-08 | 2.8142e-12 | 3.5996e-13 |
| 50 | 3.6510e-07 ± 8.9e-08 | 3.6891e-07 | 1.3140e-07 ± 3.7e-08 | 1.3184e-07 | 1.3303e-07 ± 3.6e-08 | 1.3086e-07 | 1.3183e-07 | 2.3708e-07 | 1.0366e-11 | 1.3097e-12 |
| 75 | 5.7835e-07 ± 7.4e-08 | 6.3186e-07 | 1.5430e-07 ± 3.0e-08 | 1.7319e-07 | 1.5968e-07 ± 3.3e-08 | 1.7203e-07 | 1.7317e-07 | 4.5869e-07 | 2.1615e-11 | 2.7001e-12 |
| 100 | 7.1395e-07 ± 1.3e-07 | 9.1005e-07 | 1.8124e-07 ± 3.0e-08 | 2.0428e-07 | 1.8119e-07 ± 3.1e-08 | 2.0303e-07 | 2.0424e-07 | 7.0581e-07 | 3.5819e-11 | 4.4282e-12 |
| 125 | 1.0566e-06 ± 1.5e-07 | 1.1868e-06 | 2.0191e-07 ± 3.2e-08 | 2.2626e-07 | 1.9823e-07 ± 3.4e-08 | 2.2504e-07 | 2.2620e-07 | 9.6060e-07 | 5.2439e-11 | 6.4221e-12 |
| 150 | 1.1344e-06 ± 1.6e-07 | 1.4552e-06 | 2.8096e-07 ± 8.1e-08 | 2.4312e-07 | 2.7769e-07 ± 8.2e-08 | 2.4184e-07 | 2.4305e-07 | 1.2121e-06 | 7.1081e-11 | 8.6315e-12 |
| 175 | 1.5203e-06 ± 1.8e-07 | 1.7118e-06 | 3.3507e-07 ± 9.4e-08 | 2.5790e-07 | 3.3152e-07 ± 9.4e-08 | 2.5653e-07 | 2.5780e-07 | 1.4540e-06 | 9.1452e-11 | 1.1021e-11 |
| 200 | 1.6688e-06 ± 2.1e-07 | 1.9532e-06 | 3.0488e-07 ± 5.3e-08 | 2.7057e-07 | 2.9782e-07 ± 5.2e-08 | 2.6917e-07 | 2.7046e-07 | 1.6828e-06 | 1.1333e-10 | 1.3564e-11 |
| 225 | 1.9434e-06 ± 2.1e-07 | 2.1774e-06 | 2.4968e-07 ± 4.6e-08 | 2.8069e-07 | 2.4401e-07 ± 4.4e-08 | 2.7928e-07 | 2.8056e-07 | 1.8968e-06 | 1.3653e-10 | 1.6243e-11 |
| 250 | 2.1220e-06 ± 2.0e-07 | 2.3841e-06 | 2.1120e-07 ± 4.5e-08 | 2.8862e-07 | 2.0796e-07 ± 4.4e-08 | 2.8717e-07 | 2.8846e-07 | 2.0956e-06 | 1.6093e-10 | 1.9042e-11 |
| 275 | 2.3755e-06 ± 2.2e-07 | 2.5743e-06 | 2.8778e-07 ± 6.1e-08 | 2.9497e-07 | 2.8120e-07 ± 5.7e-08 | 2.9348e-07 | 2.9479e-07 | 2.2795e-06 | 1.8640e-10 | 2.1951e-11 |
| 300 | 2.4310e-06 ± 2.4e-07 | 2.7491e-06 | 2.7989e-07 ± 5.5e-08 | 3.0023e-07 | 2.7468e-07 ± 5.3e-08 | 2.9869e-07 | 3.0001e-07 | 2.4490e-06 | 2.1285e-10 | 2.4960e-11 |

The run mean of f̂_MR − f̂_DR on the orbit's grid (common random numbers) against Ĵ(ȳ_MR − ȳ_DR): its mean square (with the run-to-run standard error of the mean function, in the same units), the prediction's, and the relative L²_ρ error of the prediction.

| t | ⟨(run mean)²⟩ | predicted | relative error ‖run mean − Ĵ(ȳ_MR − ȳ_DR)‖/‖Ĵ(ȳ_MR − ȳ_DR)‖ | sampling error of the run mean (rms) |
|---|---|---|---|---|
| 25 | 7.0616e-08 | 6.9660e-08 | 0.0527 | 1.02e-05 |
| 50 | 2.3844e-07 | 2.3790e-07 | 0.0283 | 1.15e-05 |
| 75 | 4.7282e-07 | 4.6016e-07 | 0.0337 | 1.36e-05 |
| 100 | 7.3366e-07 | 7.0792e-07 | 0.0297 | 1.38e-05 |
| 125 | 9.9049e-07 | 9.6327e-07 | 0.0215 | 1.47e-05 |
| 150 | 1.2483e-06 | 1.2153e-06 | 0.0229 | 1.38e-05 |
| 175 | 1.4628e-06 | 1.4575e-06 | 0.0152 | 1.41e-05 |
| 200 | 1.6906e-06 | 1.6865e-06 | 0.0146 | 1.33e-05 |
| 225 | 1.9022e-06 | 1.9007e-06 | 0.0123 | 1.42e-05 |
| 250 | 2.1028e-06 | 2.0997e-06 | 0.0052 | 1.78e-05 |
| 275 | 2.2983e-06 | 2.2836e-06 | 0.0070 | 1.42e-05 |
| 300 | 2.4552e-06 | 2.4532e-06 | 0.0079 | 1.36e-05 |

Numbers are produced by `experiments/exp_eiv_remedy.py`; nothing is edited by hand.
