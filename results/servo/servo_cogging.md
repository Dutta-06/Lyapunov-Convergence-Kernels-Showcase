# An application: cogging-torque compensation in a position servo (`exp_servo_cogging.py`)

The servo ẍ = −τ_cog(x) + u, τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x, tracking x_m = sin t under the sliding controller (s = ė + λe); K₀ = 4 saturated Gaussian kernels on the position, masses 1/4, σ₀ = 0.01; the averaged slow loop along the position orbit.

## A. The local minimum ϑ^b(q), the rates, and the global picture (no dynamics)

ϑ^b(q), the nondegenerate local minimum of the unregularised weighted loss that BFGS reaches from START: a = (0.982, 1.362, -1.896, -0.125), c = (-0.905, -0.335, 0.297, 0.782), s = (0.256, 0.247, 0.266, 0.097); gradient residual 5.4e-17, relative fit error 1.4620e-02. The whitened unregularised Hessian Ã's spectrum (units of σ₀): 2.89e-04, 5.40e-04, 1.13e-02, 9.22e-02, 2.58e-01, 4.43e-01, 1.11e+00, 1.70e+00, 3.50e+00, 5.57e+00, 1.51e+01, 3.01e+01 — 4 modes below 0.1σ₀.

The collective rates (units of σ₀; the slowest three): the fixed prior's ã_j + σ₀: 1.0003, 1.0005, 1.0113; the adaptive prior's λ_j^− at η = σ₀: 1.44e-04, 2.70e-04, 5.65e-03, at η = ∞: 2.89e-04, 5.40e-04, 1.13e-02.

The residual curvature: λ⁺_k/σ₀ = 0.0004, 0.0016, 0.0098, 0.0075 (Proposition D‴.28(iii): margins σ₀ − (q_k/w_k)λ⁺_k); the gauge set W_g: w > (0.0614, 0.0851, 0.1185, 0.0078). The fixed matching prior's floor law (D‴.28(ii)) per unit of 1/m − 1/M: 1.2940e-02.

The misplaced prior (the approximant of the cogging model with the phase 0.5 in place of 0.3, the engineer's model): its own fit error on the true target 2.0005e-01; the fixed prior's replicated stationary point there (residual 6.5e-17; the same from the prior's atoms and from ϑ^b(q), to 1.3e-13) has the fit error 5.4520e-02 — the leak's bias — and Proposition D‴.5's splitting margin -3.97e-03 (κ_R/σ₀ = 0.474, -0.003, 0.263, 1.397): a saddle, which the m-systems leave by splitting their classes (Part B). At the matching prior the replicated point ϑ^b(q) has the margin +9.90e-03 at w = q.

The global picture: of 120 BFGS runs from random starts (physical parameters uniform in their boxes), 1 end at ϑ^b(q) (up to relabelling), 33 below its fit error 1.4620e-02, and those all with a saturated amplitude (the largest |α| of each at least 44.1; a = a_max·tanh(α/a_max)); the end fit errors' quartiles 0.0064, 0.0116, 0.0193, 0.0264. The best end point (fit error 6.3901e-03): a = (-3.498, 0.465, -4.000, 4.000), c = (0.781, -0.955, 0.224, 0.498), s = (0.422, 0.166, 0.298, 0.735). With its 2 saturated amplitudes held at |α| = 6, 8, 10, 14 (|a| = 3.621, 3.856, 3.946, 3.993) and the rest re-optimised, the fit error is 7.1372e-03, 6.6537e-03, 6.4855e-03, 6.4029e-03: the infimum of the loss is approached at the saturation, not attained — globally the servo is Proposition D‴.27's case (critical points at infinity), and ϑ^b(q) is a local minimum, where Proposition D‴.28 applies.

## B. The averaged loop (2000 slow units, 3 seeds)

The fit error at s = 600, 1500 and the end (mean over the seeds), ‖f̂_m − f̂_M‖² at the end (mean, min–max), and the largest whitened class diameter at the end over the seeds (the classes by the nearest atom of the design's prior at the end; —: σ-modification and no leak). Start: warm clouds around ϑ^b(q) (oracle — computed from the true torque) or around the engineer's model.

| start | prior | design | m | fit error at 600, 1500, end | ‖f̂_m − f̂_M‖² at the end | class diameter |
|---|---|---|---|---|---|---|
| oracle | matching | nearest atom, fixed prior | 16 | 1.740e-01, 1.735e-01, 1.735e-01 | 1.85e-03 (6.5e-04–2.6e-03) | 1.20 |
| oracle | matching | nearest atom, fixed prior | 64 | 6.644e-02, 6.639e-02, 6.639e-02 | 2.24e-04 (2.4e-05–5.0e-04) | 0.45 |
| oracle | matching | nearest atom, fixed prior | 128 | 3.820e-02, 3.786e-02, 3.781e-02 | 4.49e-05 (6.1e-06–8.6e-05) | 0.39 |
| oracle | matching | nearest atom, adaptive prior (η = ∞) | 16 | 4.442e-02, 3.174e-02, 2.978e-02 | 4.85e-05 (9.1e-07–9.5e-05) | 0.78 |
| oracle | matching | nearest atom, adaptive prior (η = ∞) | 64 | 1.795e-02, 1.503e-02, 1.481e-02 | 2.88e-07 (1.4e-08–6.3e-07) | 0.00 |
| oracle | matching | nearest atom, adaptive prior (η = ∞) | 128 | 1.487e-02, 1.468e-02, 1.466e-02 | 6.36e-08 (2.6e-09–1.4e-07) | 0.00 |
| oracle | matching | nearest atom, adaptive prior (η = σ₀) | 16 | 6.865e-02, 4.042e-02, 3.589e-02 | 7.81e-05 (5.8e-06–1.3e-04) | 0.82 |
| oracle | matching | nearest atom, adaptive prior (η = σ₀) | 64 | 2.611e-02, 1.700e-02, 1.587e-02 | 2.10e-06 (1.8e-08–5.6e-06) | 0.00 |
| oracle | matching | nearest atom, adaptive prior (η = σ₀) | 128 | 1.621e-02, 1.479e-02, 1.473e-02 | 1.17e-07 (6.8e-09–2.9e-07) | 0.00 |
| oracle | matching | OT, fixed prior | 16 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| oracle | matching | OT, fixed prior | 64 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| oracle | matching | OT, fixed prior | 128 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| oracle | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 16 | 3.279e-01, 3.285e-01, 3.285e-01 | 1.46e-05 (9.0e-06–2.0e-05) | — |
| oracle | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 64 | 3.275e-01, 3.281e-01, 3.281e-01 | 4.37e-06 (3.1e-08–9.0e-06) | — |
| oracle | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 128 | 3.278e-01, 3.284e-01, 3.284e-01 | 5.32e-07 (2.8e-07–8.1e-07) | — |
| oracle | — | no leak (σ₀ = 0) | 16 | 2.948e-02, 1.910e-02, 1.801e-02 | 7.71e-06 (3.6e-07–1.9e-05) | — |
| oracle | — | no leak (σ₀ = 0) | 64 | 1.760e-02, 1.486e-02, 1.466e-02 | 2.43e-07 (3.4e-08–5.4e-07) | — |
| oracle | — | no leak (σ₀ = 0) | 128 | 1.473e-02, 1.459e-02, 1.457e-02 | 6.62e-08 (2.5e-09–1.0e-07) | — |
| oracle | misplaced (phase 0.5) | nearest atom, fixed prior | 16 | 1.025e-01, 1.017e-01, 1.017e-01 | 4.94e-04 (3.0e-04–7.0e-04) | 0.95 |
| oracle | misplaced (phase 0.5) | nearest atom, fixed prior | 64 | 6.236e-02, 5.974e-02, 5.969e-02 | 8.65e-05 (1.8e-05–1.5e-04) | 1.02 |
| oracle | misplaced (phase 0.5) | nearest atom, fixed prior | 128 | 5.014e-02, 4.850e-02, 4.800e-02 | 3.00e-05 (6.2e-06–5.3e-05) | 1.07 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 16 | 4.442e-02, 3.174e-02, 2.978e-02 | 4.85e-05 (9.1e-07–9.5e-05) | 0.78 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 64 | 1.795e-02, 1.503e-02, 1.481e-02 | 2.88e-07 (1.4e-08–6.3e-07) | 0.00 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 128 | 1.487e-02, 1.468e-02, 1.466e-02 | 6.36e-08 (2.6e-09–1.4e-07) | 0.00 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 16 | 4.908e-02, 3.638e-02, 3.517e-02 | 6.28e-05 (5.3e-06–1.2e-04) | 0.68 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 64 | 3.542e-02, 2.402e-02, 1.992e-02 | 5.99e-06 (4.0e-08–1.7e-05) | 0.00 |
| oracle | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 128 | 2.612e-02, 1.994e-02, 1.854e-02 | 2.15e-07 (1.2e-08–5.4e-07) | 0.00 |
| oracle | misplaced (phase 0.5) | OT, fixed prior | 16 | 5.009e-02, 4.849e-02, 4.847e-02 | 1.59e-08 (7.9e-09–2.4e-08) | 0.32 |
| oracle | misplaced (phase 0.5) | OT, fixed prior | 64 | 5.161e-02, 4.850e-02, 4.832e-02 | 9.91e-09 (2.6e-09–2.4e-08) | 0.34 |
| oracle | misplaced (phase 0.5) | OT, fixed prior | 128 | 5.194e-02, 4.854e-02, 4.839e-02 | 2.06e-09 (8.0e-11–5.4e-09) | 0.34 |
| engineer | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 16 | 3.280e-01, 3.285e-01, 3.285e-01 | 1.46e-05 (9.0e-06–2.0e-05) | — |
| engineer | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 64 | 3.275e-01, 3.281e-01, 3.281e-01 | 4.37e-06 (3.1e-08–9.0e-06) | — |
| engineer | — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 128 | 3.278e-01, 3.284e-01, 3.284e-01 | 5.32e-07 (2.8e-07–8.1e-07) | — |
| engineer | — | no leak (σ₀ = 0) | 16 | 1.873e-02, 1.566e-02, 1.519e-02 | 8.49e-07 (6.7e-07–1.1e-06) | — |
| engineer | — | no leak (σ₀ = 0) | 64 | 1.756e-02, 1.520e-02, 1.466e-02 | 3.09e-07 (9.2e-08–5.7e-07) | — |
| engineer | — | no leak (σ₀ = 0) | 128 | 1.741e-02, 1.524e-02, 1.467e-02 | 1.00e-07 (1.5e-08–1.4e-07) | — |
| engineer | misplaced (phase 0.5) | nearest atom, fixed prior | 16 | 1.026e-01, 1.017e-01, 1.017e-01 | 4.94e-04 (3.0e-04–6.9e-04) | 0.96 |
| engineer | misplaced (phase 0.5) | nearest atom, fixed prior | 64 | 6.143e-02, 5.980e-02, 5.974e-02 | 8.62e-05 (1.8e-05–1.5e-04) | 1.01 |
| engineer | misplaced (phase 0.5) | nearest atom, fixed prior | 128 | 4.946e-02, 4.810e-02, 4.764e-02 | 2.98e-05 (6.2e-06–5.3e-05) | 1.05 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 16 | 3.509e-02, 2.778e-02, 2.660e-02 | 2.32e-05 (1.6e-06–5.1e-05) | 0.50 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 64 | 3.084e-02, 2.567e-02, 2.510e-02 | 1.48e-05 (4.4e-09–4.4e-05) | 0.32 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 128 | 3.170e-02, 2.393e-02, 2.189e-02 | 4.99e-08 (7.7e-09–9.6e-08) | 0.00 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 16 | 3.989e-02, 3.424e-02, 3.189e-02 | 5.39e-05 (7.7e-06–8.6e-05) | 0.57 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 64 | 3.470e-02, 2.959e-02, 2.750e-02 | 1.89e-05 (3.2e-08–5.6e-05) | 0.49 |
| engineer | misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 128 | 3.433e-02, 3.086e-02, 2.806e-02 | 2.79e-07 (8.5e-08–5.4e-07) | 0.02 |
| engineer | misplaced (phase 0.5) | OT, fixed prior | 16 | 4.900e-02, 4.848e-02, 4.847e-02 | 2.83e-08 (9.9e-09–6.2e-08) | 0.32 |
| engineer | misplaced (phase 0.5) | OT, fixed prior | 64 | 4.965e-02, 4.839e-02, 4.830e-02 | 1.51e-08 (8.8e-10–3.8e-08) | 0.32 |
| engineer | misplaced (phase 0.5) | OT, fixed prior | 128 | 4.991e-02, 4.843e-02, 4.831e-02 | 1.10e-09 (1.2e-10–2.2e-09) | 0.33 |

The collapsed K₀-atom flows toward the misplaced prior (one particle per class), from ϑ^b(q) and from the engineer's model: the fit error and the parameter distance from ϑ^b(q) (largest coordinate difference, up to relabelling) at s = 600, 2000, 10000, 100000.

| design | start | fit error | distance from ϑ^b(q) | largest |α| at the end |
|---|---|---|---|---|
| nearest atom, fixed prior | oracle | 5.4540e-02, 5.4520e-02, 5.4520e-02, 5.4520e-02 | 0.635, 0.637, 0.637, 0.637 | 2.01 |
| nearest atom, fixed prior | engineer | 5.4522e-02, 5.4520e-02, 5.4520e-02, 5.4520e-02 | 0.637, 0.637, 0.637, 0.637 | 2.01 |
| nearest atom, adaptive prior (η = ∞) | oracle | 1.4620e-02, 1.4620e-02, 1.4620e-02, 1.4620e-02 | 0.000, 0.000, 0.000, 0.000 | 2.06 |
| nearest atom, adaptive prior (η = ∞) | engineer | 2.9742e-02, 2.1329e-02, 1.6890e-02, 1.5542e-02 | 0.694, 0.755, 0.792, 0.553 | 2.03 |
| nearest atom, adaptive prior (η = σ₀) | oracle | 2.4043e-02, 1.8269e-02, 1.5445e-02, 1.4912e-02 | 0.324, 0.338, 0.351, 0.309 | 2.04 |
| nearest atom, adaptive prior (η = σ₀) | engineer | 3.4536e-02, 2.5922e-02, 1.8118e-02, 1.5982e-02 | 0.669, 0.721, 0.790, 0.683 | 2.03 |
| no leak (σ₀ = 0) | oracle | 1.4620e-02, 1.4620e-02, 1.4620e-02, 1.4620e-02 | 0.000, 0.000, 0.000, 0.000 | 2.06 |
| no leak (σ₀ = 0) | engineer | 2.9742e-02, 2.1329e-02, 1.6890e-02, 1.5542e-02 | 0.694, 0.755, 0.792, 0.553 | 2.03 |

## C. The real servo (k_s = 25, λ = 5, dt = 0.002, 1000 time units), seed 0, from the oracle-start end states

The position-tracking error's rms over the five periods ending at t = 31, 500 and the end (records every 0.02 time units), against the averaged predictions at the parameters reached then: the residual r(x_m) filtered by (d/dt + λ)⁻¹(d/dt + k_s)⁻¹ (exact) and by (d/dt + λ)⁻¹ after the slaved s = r/k_s; their ratios.

| prior | design | m | e_rms at t = 31, 500, end | exact-filter prediction | ratios | slaved prediction | ratios |
|---|---|---|---|---|---|---|---|
| matching | nearest atom, fixed prior | 16 | 1.7986e-04, 1.8120e-04, 1.8156e-04 | 1.8145e-04, 1.8117e-04, 1.8143e-04 | 0.9912, 1.0002, 1.0007 | 1.8299e-04, 1.8273e-04, 1.8299e-04 | 0.9829, 0.9916, 0.9922 |
| matching | nearest atom, fixed prior | 64 | 9.1443e-05, 9.2229e-05, 9.2374e-05 | 9.2156e-05, 9.2145e-05, 9.2294e-05 | 0.9923, 1.0009, 1.0009 | 9.2634e-05, 9.2633e-05, 9.2787e-05 | 0.9871, 0.9956, 0.9956 |
| matching | nearest atom, adaptive prior (η = ∞) | 16 | 1.7042e-05, 1.7034e-05, 1.7034e-05 | 1.7006e-05, 1.7007e-05, 1.7008e-05 | 1.0021, 1.0016, 1.0015 | 1.7710e-05, 1.7710e-05, 1.7710e-05 | 0.9623, 0.9618, 0.9618 |
| matching | nearest atom, adaptive prior (η = ∞) | 64 | 1.8687e-05, 1.8683e-05, 1.8679e-05 | 1.8651e-05, 1.8650e-05, 1.8646e-05 | 1.0019, 1.0018, 1.0018 | 1.9198e-05, 1.9196e-05, 1.9192e-05 | 0.9734, 0.9733, 0.9733 |
| matching | nearest atom, adaptive prior (η = σ₀) | 16 | 1.8834e-05, 1.8702e-05, 1.8607e-05 | 1.8782e-05, 1.8679e-05, 1.8580e-05 | 1.0028, 1.0012, 1.0015 | 1.9677e-05, 1.9573e-05, 1.9472e-05 | 0.9571, 0.9555, 0.9556 |
| matching | nearest atom, adaptive prior (η = σ₀) | 64 | 1.8255e-05, 1.8251e-05, 1.8247e-05 | 1.8222e-05, 1.8220e-05, 1.8216e-05 | 1.0018, 1.0017, 1.0017 | 1.8840e-05, 1.8837e-05, 1.8834e-05 | 0.9690, 0.9689, 0.9689 |
| matching | OT, fixed prior | 16 | 1.8886e-05, 1.8877e-05, 1.8869e-05 | 1.8849e-05, 1.8842e-05, 1.8834e-05 | 1.0020, 1.0019, 1.0019 | 1.9368e-05, 1.9361e-05, 1.9353e-05 | 0.9751, 0.9750, 0.9750 |
| matching | OT, fixed prior | 64 | 1.8886e-05, 1.8877e-05, 1.8869e-05 | 1.8849e-05, 1.8842e-05, 1.8834e-05 | 1.0020, 1.0019, 1.0019 | 1.9368e-05, 1.9361e-05, 1.9353e-05 | 0.9751, 0.9750, 0.9750 |
| — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 16 | 5.5336e-04, 5.5230e-04, 5.5292e-04 | 5.5068e-04, 5.5141e-04, 5.5131e-04 | 1.0049, 1.0016, 1.0029 | 5.5714e-04, 5.5787e-04, 5.5777e-04 | 0.9932, 0.9900, 0.9913 |
| — | σ-modification (toward ϑ⁰ = 0, the same σ₀) | 64 | 5.4569e-04, 5.4465e-04, 5.4527e-04 | 5.4302e-04, 5.4375e-04, 5.4368e-04 | 1.0049, 1.0017, 1.0029 | 5.4940e-04, 5.5013e-04, 5.5006e-04 | 0.9932, 0.9900, 0.9913 |
| — | no leak (σ₀ = 0) | 16 | 1.8960e-05, 1.8931e-05, 1.8903e-05 | 1.8923e-05, 1.8896e-05, 1.8869e-05 | 1.0020, 1.0018, 1.0018 | 1.9463e-05, 1.9438e-05, 1.9411e-05 | 0.9741, 0.9740, 0.9738 |
| — | no leak (σ₀ = 0) | 64 | 1.8819e-05, 1.8798e-05, 1.8775e-05 | 1.8783e-05, 1.8763e-05, 1.8741e-05 | 1.0019, 1.0019, 1.0019 | 1.9284e-05, 1.9265e-05, 1.9244e-05 | 0.9759, 0.9758, 0.9757 |
| misplaced | nearest atom, fixed prior | 16 | 1.3818e-04, 1.3986e-04, 1.4021e-04 | 1.3962e-04, 1.3978e-04, 1.4013e-04 | 0.9896, 1.0006, 1.0005 | 1.4098e-04, 1.4115e-04, 1.4151e-04 | 0.9801, 0.9909, 0.9908 |
| misplaced | nearest atom, fixed prior | 64 | 9.4523e-05, 9.5872e-05, 9.6122e-05 | 9.5509e-05, 9.5756e-05, 9.6109e-05 | 0.9897, 1.0012, 1.0001 | 9.6446e-05, 9.6698e-05, 9.7055e-05 | 0.9801, 0.9915, 0.9904 |
| misplaced | nearest atom, adaptive prior (η = ∞) | 16 | 1.7042e-05, 1.7034e-05, 1.7034e-05 | 1.7006e-05, 1.7007e-05, 1.7008e-05 | 1.0021, 1.0016, 1.0015 | 1.7710e-05, 1.7710e-05, 1.7710e-05 | 0.9623, 0.9618, 0.9618 |
| misplaced | nearest atom, adaptive prior (η = ∞) | 64 | 1.8687e-05, 1.8683e-05, 1.8679e-05 | 1.8651e-05, 1.8650e-05, 1.8646e-05 | 1.0019, 1.0018, 1.0018 | 1.9198e-05, 1.9196e-05, 1.9192e-05 | 0.9734, 0.9733, 0.9733 |
| misplaced | nearest atom, adaptive prior (η = σ₀) | 16 | 2.0821e-05, 2.0642e-05, 2.0512e-05 | 2.0752e-05, 2.0614e-05, 2.0498e-05 | 1.0033, 1.0014, 1.0007 | 2.2223e-05, 2.2086e-05, 2.1969e-05 | 0.9369, 0.9346, 0.9337 |
| misplaced | nearest atom, adaptive prior (η = σ₀) | 64 | 1.7806e-05, 1.7727e-05, 1.7681e-05 | 1.7761e-05, 1.7703e-05, 1.7670e-05 | 1.0026, 1.0013, 1.0006 | 1.9086e-05, 1.9028e-05, 1.8990e-05 | 0.9329, 0.9316, 0.9310 |
| misplaced | OT, fixed prior | 16 | 7.7738e-05, 7.7884e-05, 7.8048e-05 | 7.7684e-05, 7.7822e-05, 7.8039e-05 | 1.0007, 1.0008, 1.0001 | 7.8619e-05, 7.8767e-05, 7.8989e-05 | 0.9888, 0.9888, 0.9881 |
| misplaced | OT, fixed prior | 64 | 7.7484e-05, 7.7639e-05, 7.7807e-05 | 7.7431e-05, 7.7578e-05, 7.7799e-05 | 1.0007, 1.0008, 1.0001 | 7.8359e-05, 7.8515e-05, 7.8742e-05 | 0.9888, 0.9888, 0.9881 |

## D. Phase-plane kernels: the whitened Gauss–Newton spectrum along the orbit (Proposition D‴.29(iv))

K = 5 isotropic kernels (tests/test_sliding.py's atoms), weights 1/K, σ₀ = 0.01: the number of eigenvalues below 10⁻⁶σ₀ and below 10⁻¹⁵ of the largest, the smallest (units of σ₀), and the (K + 1)-th relative to the largest.

| orbit | dimension | below 10⁻⁶σ₀ | below 10⁻¹⁵ of the largest | smallest/σ₀ | (K + 1)-th/largest |
|---|---|---|---|---|---|
| circle | 20 | 8 | 5 | -3.49e-16 | 2.40e-10 |
| ellipse, velocity axis ×1.1 | 20 | 8 | 0 | 8.80e-15 | 2.09e-09 |
| ellipse, velocity axis ×2 | 20 | 6 | 0 | 2.73e-11 | 5.61e-07 |
| x_m = 1.2 sin t + 0.3 sin 2t | 20 | 3 | 0 | 5.41e-09 | 1.98e-06 |
| circle, the centres at the origin | 20 | 17 | 16 | -4.42e-15 | -4.57e-17 |

Numbers are produced by `experiments/exp_servo_cogging.py`; nothing is edited by hand.
