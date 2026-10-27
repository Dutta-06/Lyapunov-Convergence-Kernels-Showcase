# Proposition D‴.49: the rich transit in the closed loop (`exp_closed_transit.py`)

ẋ = f(x) + u, u = ẋ_m − k_ee − f̂(x), k_e = 10, x_m = x₀ + a sin(ωt), ω = 0.2; f three bumps of width 0.12 at [-0.52, 0.07, 0.61] with amplitudes [0.9, -0.6, 0.75]; 21 kernels of width 0.12 started on a grid of [−1, 1] at α = 0 with the masses ε; the free law (m = 12, G = 1, γ = 0.1, a_max = 1) with moving centres; τ = t/Λ(ε), Λ = (πa_max/(4G))(2/(cε))^{1/2}. A kernel has formed when |w| > 0.01·max|A*|; entries listed as kernel: τ (kernel i started at −1 + 0.1i).

## A. The averaged loop (the reference's orbit sampled over a period, the error filter 1/(s + k_e))

| ε | entries (kernel: τ) | the last entry's time t | formed | centres | coefficients | rms residual | inactive: max |Δc| |
|---|---|---|---|---|---|---|---|
| 0.0001 | 5: 142, 4: 149.5, 16: 156.5, 17: 168, 11: 254, 10: 273.5 | 1.92e+05 | 6 | [-0.5201, -0.52, 0.07, 0.07, 0.61, 0.6104] | [0.0316, 0.8684, -0.0115, -0.5885, 0.737, 0.013] | 1.31e-6 | 0.0591 |
| 1e-05 | 5: 142.5, 16: 156.5, 11: 254.5 | 5.65e+05 | 3 | [-0.52, 0.07, 0.61] | [0.8999, -0.6, 0.75] | 3.87e-7 | 0.062 |
| 1e-06 | 5: 142.5, 16: 157, 11: 254.5 | 1.79e+06 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 3.38e-8 | 0.00442 |
| 1e-07 | 5: 142.5, 16: 157, 11: 254.5 | 5.65e+06 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 9.45e-10 | 3.76e-4 |

## B. The real loop against the averaged loop (x_m = sin ωt)

| ε | | entries (kernel: τ) | formed | centres | coefficients | inactive: max |Δc| | max |e|: first, last period |
|---|---|---|---|---|---|---|---|
| 0.0001 | averaged | 5: 142, 4: 149.5, 16: 156.5, 17: 168, 11: 254, 10: 273.5 | 6 | [-0.5201, -0.52, 0.07, 0.07, 0.61, 0.6104] | [0.0316, 0.8684, -0.0115, -0.5885, 0.737, 0.013] | 0.0591 | — |
| 0.0001 | real (T = 2.5e+05, 1227 s) | 5: 132.9, 17: 157.3, 16: 164, 10: 253.3, 11: 260 | 5 | [-0.52, 0.07, 0.07, 0.6398, 0.6083] | [0.8997, -0.5351, -0.0649, 0.0306, 0.7203] | 0.119 | 0.0893, 1.63e-5 |
| 1e-05 | averaged | 5: 142.5, 16: 156.5, 11: 254.5 | 3 | [-0.52, 0.07, 0.61] | [0.8999, -0.6, 0.75] | 0.062 | — |
| 1e-05 | real (T = 7.35e+05, 4155 s) | 5: 133.2, 17: 158, 10: 254.3 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.5993, 0.7499] | 0.0346 | 0.0893, 3.58e-8 |

## C. A reference that does not visit the bump at −0.52 (x_m = 0.45 + 0.5 sin ωt, the range [-0.05, 0.95])

| | ε | entries (kernel: τ) | formed | centres | coefficients |
|---|---|---|---|---|---|
| averaged | 1e-06 | 10: 74.5, 16: 92.5 | 2 | [0.07, 0.61] | [-0.5999, 0.75] |
| averaged | 0.0001 | 10: 74, 11: 76, 16: 92.5 | 3 | [0.07, 0.07, 0.61] | [-0.5276, -0.0723, 0.742] |
| real | 0.0001 | 10: 76.36, 17: 92.29, 16: 95.15 | 3 | [0.07, 0.61, 0.6099] | [-0.5998, 0.0501, 0.6992] |

## D. The positions' metric (the averaged loop, ε = 1e-06)

- Wasserstein–Fisher–Rao (ċ = Ga e∂_cφ): 6 kernels formed in 3 clusters at [-0.52, 0.07, 0.61], rms residual 8.31e-11.
- fixed centres: 6 kernels formed in 6 clusters at [-0.6, -0.5, 0.0, 0.1, 0.6, 0.7], rms residual 0.0127.

Numbers are produced by `experiments/exp_closed_transit.py`; nothing is edited by hand.
