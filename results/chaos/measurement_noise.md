# Proposition D‴.30: measurement noise at perfect tracking (`exp_measurement_noise.py`)

The three-bump instance at σ₀ = 0.15σ_crit = 0.00707, x_m = 2 sin t, the nearest-atom leak toward a fixed (η = 0) or an adaptive prior, at the lift of the exact representation of q = 1/3 with the prior at its atoms; the measurement y = x + ν in the controller and the adaptive law.

## A. The linear theory and the second-order force (no dynamics)

k_e = 25, N₀ = 1. Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩, a = ∂_xf̂ on the orbit: tr Q = 8.7836 = 8.1049 (⟨Ĵ̇ᵀĴ̇⟩) +0.0000 (cross) +0.6786 (⟨a²ĴᵀĴ⟩). Per mode (the eigenbasis of H): h_j, Q_jj, Q_jj/h_j, and the function-level shares of the floor, (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀), for each η.

| j | h_j | Q_jj | Q_jj/h_j | share at η = 0 (fixed) | share at η = 0.00707 | share at η = 0.07074 | share at η = ∞ |
|---|---|---|---|---|---|---|---|
| 1 | 3.0587e-04 | 3.1925e-02 | 104.37 | 2.6462e-05 | 3.2601e-04 | 5.8068e-04 | 6.3850e-04 |
| 2 | 1.9644e-03 | 1.0304e-01 | 52.45 | 4.4787e-04 | 1.1560e-03 | 1.8780e-03 | 2.0607e-03 |
| 3 | 3.1144e-03 | 1.3748e-01 | 44.14 | 8.4049e-04 | 1.6229e-03 | 2.5093e-03 | 2.7497e-03 |
| 4 | 1.5149e-02 | 2.7295e-01 | 18.02 | 3.7213e-03 | 4.1409e-03 | 5.0437e-03 | 5.4591e-03 |
| 5 | 2.1104e-02 | 3.1305e-01 | 14.83 | 4.6891e-03 | 5.0046e-03 | 5.8132e-03 | 6.2610e-03 |
| 6 | 3.2479e-02 | 1.7803e-01 | 5.48 | 2.9237e-03 | 3.0203e-03 | 3.3322e-03 | 3.5605e-03 |
| 7 | 8.8615e-02 | 2.2602e+00 | 25.51 | 4.1863e-02 | 4.2093e-02 | 4.3283e-02 | 4.5205e-02 |
| 8 | 1.3420e-01 | 2.4514e+00 | 18.27 | 4.6574e-02 | 4.6691e-02 | 4.7393e-02 | 4.9029e-02 |
| 9 | 1.9521e-01 | 3.0354e+00 | 15.55 | 5.8585e-02 | 5.8657e-02 | 5.9135e-02 | 6.0708e-02 |

| η | floor (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) | relative to η = ∞ | the weakest mode's parameter variance | its bias-removal rate λ_1^− |
|---|---|---|---|---|
| 0 (fixed) | 1.59670e-01 | 0.9089 | 8.6515e-02 | 0.000e+00 |
| 0.00707 | 1.62711e-01 | 0.9262 | 1.0658e+00 | 1.513e-04 |
| 0.07074 | 1.68968e-01 | 0.9618 | 1.8985e+00 | 2.780e-04 |
| ∞ | 1.75671e-01 | 1.0000 | 2.0875e+00 | 3.059e-04 |

At η = ∞ the floor is (N₀/2k_e)tr Q = 1.75671e-01 whatever H. The second-order force (errors in variables; OU noise of rate B = 100), slow time, per unit of the noise's variance, in the eigenbasis of H: +2.133e-02, +1.880e+00, +6.659e-01, +1.583e-01, +4.332e-01, +9.929e-01, -4.724e+00, -4.150e-02, +5.619e+00. The stationary mean deviation's function-level size per unit of the variance squared: fixed prior Σ_jF_j²h_j/(h_j + σ₀)² = 4.9106e+02, adaptive prior (every η > 0) Σ_jF_j²/h_j = 2.3969e+03 (ratio 4.88).

The fast parts: tr H² = 0.0657 (the function-level floor) and ⟨|Ĵ|⁴⟩ = 0.2991 (the compensation error at the current point), each times N₀/(2k_e). With the OU noise (N₀ = 2Eν²/B) the errors-in-variables bias equals the floor of the parameters at Eν² = 6.50e-06 (std 0.0026) under the fixed prior and 1.47e-06 (std 0.0012) under the adaptive one, and dominates above.

## B. The exactly linearised loop with white measurement noise against the averaged prediction

The periodic stationary covariance P(t) (the monodromy and a discrete Lyapunov equation; N₀ = 1), averaged over a period: the function-level floor E‖f̂ − f‖²_ρ = ⟨tr(HP_yy)⟩ and the compensation error at the plant's current point ⟨ĴP_yyĴᵀ⟩, against the prediction — the slow (parameter) part (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) plus the fast part (N₀/2k_e)tr H² and (N₀/2k_e)⟨|Ĵ|⁴⟩ respectively.

| k_e | η | E‖f̂ − f‖²_ρ: exact | predicted | ratio | at the current point: exact | predicted | ratio | the slow part alone |
|---|---|---|---|---|---|---|---|---|
| 25 | 0 (fixed) | 1.57161e-01 | 1.60985e-01 | 0.9762 | 1.61703e-01 | 1.65653e-01 | 0.9762 | 1.59670e-01 |
| 50 | 0 (fixed) | 7.99990e-02 | 8.04923e-02 | 0.9939 | 8.23123e-02 | 8.28263e-02 | 0.9938 | 7.98351e-02 |
| 100 | 0 (fixed) | 4.01892e-02 | 4.02461e-02 | 0.9986 | 4.13526e-02 | 4.14131e-02 | 0.9985 | 3.99176e-02 |
| 200 | 0 (fixed) | 2.01175e-02 | 2.01231e-02 | 0.9997 | 2.07002e-02 | 2.07066e-02 | 0.9997 | 1.99588e-02 |
| 25 | 0.00707 | 1.60323e-01 | 1.64025e-01 | 0.9774 | 1.64871e-01 | 1.68693e-01 | 0.9773 | 1.62711e-01 |
| 50 | 0.00707 | 8.15347e-02 | 8.20126e-02 | 0.9942 | 8.38487e-02 | 8.43466e-02 | 0.9941 | 8.13554e-02 |
| 100 | 0.00707 | 4.09513e-02 | 4.10063e-02 | 0.9987 | 4.21148e-02 | 4.21733e-02 | 0.9986 | 4.06777e-02 |
| 200 | 0.00707 | 2.04978e-02 | 2.05031e-02 | 0.9997 | 2.10806e-02 | 2.10866e-02 | 0.9997 | 2.03389e-02 |
| 25 | ∞ | 1.73501e-01 | 1.76985e-01 | 0.9803 | 1.78066e-01 | 1.81653e-01 | 0.9803 | 1.75671e-01 |
| 50 | ∞ | 8.80421e-02 | 8.84927e-02 | 0.9949 | 9.03585e-02 | 9.08267e-02 | 0.9948 | 8.78356e-02 |
| 100 | ∞ | 4.41949e-02 | 4.42464e-02 | 0.9988 | 4.53586e-02 | 4.54134e-02 | 0.9988 | 4.39178e-02 |
| 200 | ∞ | 2.21183e-02 | 2.21232e-02 | 0.9998 | 2.27011e-02 | 2.27067e-02 | 0.9998 | 2.19589e-02 |

## C. The real loop with an Ornstein–Uhlenbeck measurement noise (k_e = 25, B = 100, RK4 dt = 0.001, 300 time units)

E‖f̂ − f‖²_ρ at the records (mean ± standard error over the runs) against the linear prediction: the transient covariance of the exactly linearised loop with the OU noise plus the second-order mean deviation (their sum, and the two parts).

Noise std 0.003 (variance 9.0e-06), η = 0 (fixed), 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 6.6096e-09 ± 1.9e-09 | 7.3368e-09 | 6.7785e-09 | 5.5832e-10 |
| 50 | 1.1271e-08 ± 2.4e-09 | 1.3615e-08 | 1.1720e-08 | 1.8949e-09 |
| 75 | 2.5858e-08 ± 9.1e-09 | 1.8956e-08 | 1.5312e-08 | 3.6437e-09 |
| 100 | 3.2351e-08 ± 8.9e-09 | 2.3547e-08 | 1.7972e-08 | 5.5746e-09 |
| 125 | 3.0619e-08 ± 7.0e-09 | 2.7362e-08 | 1.9816e-08 | 7.5457e-09 |
| 150 | 3.0917e-08 ± 7.9e-09 | 3.0681e-08 | 2.1208e-08 | 9.4728e-09 |
| 175 | 3.3660e-08 ± 6.1e-09 | 3.3732e-08 | 2.2424e-08 | 1.1308e-08 |
| 200 | 3.5042e-08 ± 5.8e-09 | 3.6488e-08 | 2.3461e-08 | 1.3027e-08 |
| 225 | 3.7653e-08 ± 7.7e-09 | 3.8899e-08 | 2.4278e-08 | 1.4621e-08 |
| 250 | 3.5510e-08 ± 7.4e-09 | 4.0995e-08 | 2.4908e-08 | 1.6087e-08 |
| 275 | 4.3328e-08 ± 8.5e-09 | 4.2835e-08 | 2.5405e-08 | 1.7430e-08 |
| 300 | 3.7097e-08 ± 7.7e-09 | 4.4468e-08 | 2.5810e-08 | 1.8657e-08 |

Noise std 0.003 (variance 9.0e-06), η = ∞, 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 6.6504e-09 ± 1.9e-09 | 7.3833e-09 | 6.8211e-09 | 5.6218e-10 |
| 50 | 1.1399e-08 ± 2.4e-09 | 1.3785e-08 | 1.1864e-08 | 1.9204e-09 |
| 75 | 2.6229e-08 ± 9.2e-09 | 1.9301e-08 | 1.5585e-08 | 3.7154e-09 |
| 100 | 3.3026e-08 ± 9.1e-09 | 2.4099e-08 | 1.8382e-08 | 5.7171e-09 |
| 125 | 3.1524e-08 ± 7.2e-09 | 2.8139e-08 | 2.0358e-08 | 7.7809e-09 |
| 150 | 3.1970e-08 ± 8.1e-09 | 3.1692e-08 | 2.1874e-08 | 9.8183e-09 |
| 175 | 3.4940e-08 ± 6.4e-09 | 3.4980e-08 | 2.3202e-08 | 1.1778e-08 |
| 200 | 3.6446e-08 ± 6.0e-09 | 3.7971e-08 | 2.4341e-08 | 1.3630e-08 |
| 225 | 3.9364e-08 ± 8.0e-09 | 4.0614e-08 | 2.5250e-08 | 1.5364e-08 |
| 250 | 3.7337e-08 ± 7.8e-09 | 4.2936e-08 | 2.5961e-08 | 1.6974e-08 |
| 275 | 4.5534e-08 ± 8.9e-09 | 4.4995e-08 | 2.6531e-08 | 1.8464e-08 |
| 300 | 3.9309e-08 ± 8.1e-09 | 4.6839e-08 | 2.7001e-08 | 1.9837e-08 |

Noise std 0.01 (variance 1.0e-04), η = 0 (fixed), 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 1.5887e-07 ± 3.2e-08 | 1.4425e-07 | 7.5317e-08 | 6.8929e-08 |
| 50 | 4.0936e-07 ± 5.1e-08 | 3.6416e-07 | 1.3022e-07 | 2.3394e-07 |
| 75 | 8.1396e-07 ± 1.7e-07 | 6.1997e-07 | 1.7013e-07 | 4.4984e-07 |
| 100 | 1.1820e-06 ± 1.7e-07 | 8.8791e-07 | 1.9969e-07 | 6.8822e-07 |
| 125 | 1.3691e-06 ± 1.7e-07 | 1.1518e-06 | 2.2018e-07 | 9.3157e-07 |
| 150 | 1.5286e-06 ± 1.8e-07 | 1.4051e-06 | 2.3565e-07 | 1.1695e-06 |
| 175 | 1.7691e-06 ± 1.6e-07 | 1.6452e-06 | 2.4915e-07 | 1.3961e-06 |
| 200 | 1.9403e-06 ± 1.8e-07 | 1.8690e-06 | 2.6067e-07 | 1.6083e-06 |
| 225 | 2.1198e-06 ± 2.0e-07 | 2.0748e-06 | 2.6976e-07 | 1.8050e-06 |
| 250 | 2.2769e-06 ± 2.1e-07 | 2.2628e-06 | 2.7675e-07 | 1.9860e-06 |
| 275 | 2.5985e-06 ± 2.4e-07 | 2.4341e-06 | 2.8227e-07 | 2.1519e-06 |
| 300 | 2.5140e-06 ± 2.6e-07 | 2.5902e-06 | 2.8678e-07 | 2.3034e-06 |

Noise std 0.01 (variance 1.0e-04), η = ∞, 12 runs:

| t | Monte Carlo | prediction | covariance part | mean part |
|---|---|---|---|---|
| 25 | 1.5990e-07 ± 3.2e-08 | 1.4519e-07 | 7.5790e-08 | 6.9404e-08 |
| 50 | 4.1443e-07 ± 5.2e-08 | 3.6891e-07 | 1.3183e-07 | 2.3708e-07 |
| 75 | 8.2804e-07 ± 1.7e-07 | 6.3186e-07 | 1.7317e-07 | 4.5869e-07 |
| 100 | 1.2094e-06 ± 1.7e-07 | 9.1005e-07 | 2.0424e-07 | 7.0581e-07 |
| 125 | 1.4107e-06 ± 1.8e-07 | 1.1868e-06 | 2.2620e-07 | 9.6060e-07 |
| 150 | 1.5839e-06 ± 1.9e-07 | 1.4552e-06 | 2.4305e-07 | 1.2121e-06 |
| 175 | 1.8411e-06 ± 1.7e-07 | 1.7118e-06 | 2.5780e-07 | 1.4540e-06 |
| 200 | 2.0287e-06 ± 1.9e-07 | 1.9532e-06 | 2.7046e-07 | 1.6828e-06 |
| 225 | 2.2266e-06 ± 2.1e-07 | 2.1774e-06 | 2.8056e-07 | 1.8968e-06 |
| 250 | 2.4013e-06 ± 2.2e-07 | 2.3841e-06 | 2.8846e-07 | 2.0956e-06 |
| 275 | 2.7458e-06 ± 2.5e-07 | 2.5743e-06 | 2.9479e-07 | 2.2795e-06 |
| 300 | 2.6739e-06 ± 2.7e-07 | 2.7491e-06 | 3.0001e-07 | 2.4490e-06 |

The lift m = 15 (counts 5, 5, 5) driven by the same noise path as the K₀-atom run (seed 0, std 0.01, η = ∞): the largest difference of f̂ on the grid over the records 0.0e+00, of E‖f̂ − f‖² 6.0e-19 — the noise is coherent: a replicated start stays replicated.

Numbers are produced by `experiments/exp_measurement_noise.py`; nothing is edited by hand.
