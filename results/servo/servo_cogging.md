# An application: cogging-torque compensation in a position servo (`exp_servo_cogging.py`)

The servo ẍ = −τ_cog(x) + u, τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x, tracking x_m = sin t under the sliding controller (s = ė + λe); K₀ = 4 saturated Gaussian kernels on the position, masses 1/4, σ₀ = 0.01; the averaged slow loop along the position orbit.

## A. The approximant and the rates (no dynamics)

The best weighted approximant ϑ^b(q): a = (0.982, 1.362, -1.896, -0.125), c = (-0.905, -0.335, 0.297, 0.782), s = (0.256, 0.247, 0.266, 0.097); gradient residual 5.4e-17, relative fit error 1.4620e-02 (no exact representation). The whitened unregularised Hessian Ã's spectrum (units of σ₀): 2.89e-04, 5.40e-04, 1.13e-02, 9.22e-02, 2.58e-01, 4.43e-01, 1.11e+00, 1.70e+00, 3.50e+00, 5.57e+00, 1.51e+01, 3.01e+01.

The collective rates (units of σ₀; the slowest three): the fixed prior's ã_j + σ₀: 1.0003, 1.0005, 1.0113; the adaptive prior's λ_j^− at η = σ₀: 1.44e-04, 2.70e-04, 5.65e-03, at η = ∞: 2.89e-04, 5.40e-04, 1.13e-02.

The residual curvature: λ⁺_k/σ₀ = 0.0004, 0.0016, 0.0098, 0.0075 (Proposition D‴.28(iii): margins σ₀ − (q_k/w_k)λ⁺_k); the gauge set W_g: w > (0.0614, 0.0851, 0.1185, 0.0078). The fixed matching prior's floor law (D‴.28(ii)) per unit of 1/m − 1/M: 1.2940e-02.

The misplaced prior (the approximant of the cogging model with the phase 0.5 in place of 0.3): its own fit error on the true target 2.0005e-01; the fixed prior's replicated stationary point there (residual 6.5e-17; the same from the prior's atoms and from ϑ^b(q), to 1.3e-13) has the fit error 5.4520e-02 — the leak's bias, which the adaptive prior removes — and Proposition D‴.5's splitting margin -3.97e-03 (κ_R/σ₀ = 0.474, -0.003, 0.263, 1.397): a saddle, which the m-systems leave by splitting their classes (Part B). At the matching prior the replicated point ϑ^b(q) has the margin +9.90e-03 at w = q.

## B. The averaged loop from warm clouds around ϑ^b(q) (2000 slow units, 3 seeds)

The fit error at s = 600, 1500 and the end (mean over the seeds), ‖f̂_m − f̂_M‖² at the end (mean, min–max), and the largest whitened class diameter at the end over the seeds (the classes by the nearest atom of the design's prior at the end; —: σ-modification, a single atom).

| prior | design | m | fit error at 600, 1500, end | ‖f̂_m − f̂_M‖² at the end | class diameter |
|---|---|---|---|---|---|
| matching | nearest atom, fixed prior | 16 | 1.740e-01, 1.735e-01, 1.735e-01 | 1.85e-03 (6.5e-04–2.6e-03) | 1.20 |
| matching | nearest atom, fixed prior | 64 | 6.644e-02, 6.639e-02, 6.639e-02 | 2.24e-04 (2.4e-05–5.0e-04) | 0.45 |
| matching | nearest atom, fixed prior | 128 | 3.820e-02, 3.786e-02, 3.781e-02 | 4.49e-05 (6.1e-06–8.6e-05) | 0.39 |
| matching | nearest atom, adaptive prior (η = ∞) | 16 | 4.442e-02, 3.174e-02, 2.978e-02 | 4.85e-05 (9.1e-07–9.5e-05) | 0.78 |
| matching | nearest atom, adaptive prior (η = ∞) | 64 | 1.795e-02, 1.503e-02, 1.481e-02 | 2.88e-07 (1.4e-08–6.3e-07) | 0.00 |
| matching | nearest atom, adaptive prior (η = ∞) | 128 | 1.487e-02, 1.468e-02, 1.466e-02 | 6.36e-08 (2.6e-09–1.4e-07) | 0.00 |
| matching | nearest atom, adaptive prior (η = σ₀) | 16 | 6.865e-02, 4.042e-02, 3.589e-02 | 7.81e-05 (5.8e-06–1.3e-04) | 0.82 |
| matching | nearest atom, adaptive prior (η = σ₀) | 64 | 2.611e-02, 1.700e-02, 1.587e-02 | 2.10e-06 (1.8e-08–5.6e-06) | 0.00 |
| matching | nearest atom, adaptive prior (η = σ₀) | 128 | 1.621e-02, 1.479e-02, 1.473e-02 | 1.17e-07 (6.8e-09–2.9e-07) | 0.00 |
| matching | OT, fixed prior | 16 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| matching | OT, fixed prior | 64 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| matching | OT, fixed prior | 128 | 1.462e-02, 1.462e-02, 1.462e-02 | 0.00e+00 (0.0e+00–0.0e+00) | 0.00 |
| σ-mod | σ-modification (ϑ⁰ = 0) | 16 | 3.279e-01, 3.285e-01, 3.285e-01 | 1.46e-05 (9.0e-06–2.0e-05) | — |
| σ-mod | σ-modification (ϑ⁰ = 0) | 64 | 3.275e-01, 3.281e-01, 3.281e-01 | 4.37e-06 (3.1e-08–9.0e-06) | — |
| σ-mod | σ-modification (ϑ⁰ = 0) | 128 | 3.278e-01, 3.284e-01, 3.284e-01 | 5.32e-07 (2.8e-07–8.1e-07) | — |
| misplaced (phase 0.5) | nearest atom, fixed prior | 16 | 1.025e-01, 1.017e-01, 1.017e-01 | 4.94e-04 (3.0e-04–7.0e-04) | 0.95 |
| misplaced (phase 0.5) | nearest atom, fixed prior | 64 | 6.236e-02, 5.974e-02, 5.969e-02 | 8.65e-05 (1.8e-05–1.5e-04) | 1.02 |
| misplaced (phase 0.5) | nearest atom, fixed prior | 128 | 5.014e-02, 4.850e-02, 4.800e-02 | 3.00e-05 (6.2e-06–5.3e-05) | 1.07 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 16 | 4.442e-02, 3.174e-02, 2.978e-02 | 4.85e-05 (9.1e-07–9.5e-05) | 0.78 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 64 | 1.795e-02, 1.503e-02, 1.481e-02 | 2.88e-07 (1.4e-08–6.3e-07) | 0.00 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = ∞) | 128 | 1.487e-02, 1.468e-02, 1.466e-02 | 6.36e-08 (2.6e-09–1.4e-07) | 0.00 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 16 | 4.908e-02, 3.638e-02, 3.517e-02 | 6.28e-05 (5.3e-06–1.2e-04) | 0.68 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 64 | 3.542e-02, 2.402e-02, 1.992e-02 | 5.99e-06 (4.0e-08–1.7e-05) | 0.00 |
| misplaced (phase 0.5) | nearest atom, adaptive prior (η = σ₀) | 128 | 2.612e-02, 1.994e-02, 1.854e-02 | 2.15e-07 (1.2e-08–5.4e-07) | 0.00 |
| misplaced (phase 0.5) | OT, fixed prior | 16 | 5.009e-02, 4.849e-02, 4.847e-02 | 1.59e-08 (7.9e-09–2.4e-08) | 0.32 |
| misplaced (phase 0.5) | OT, fixed prior | 64 | 5.161e-02, 4.850e-02, 4.832e-02 | 9.91e-09 (2.6e-09–2.4e-08) | 0.34 |
| misplaced (phase 0.5) | OT, fixed prior | 128 | 5.194e-02, 4.854e-02, 4.839e-02 | 2.06e-09 (8.0e-11–5.4e-09) | 0.34 |

## C. The real servo (k_s = 25, λ = 5, dt = 0.002, 1000 time units), seed 0

From the end state of each averaged design: the position-tracking error's rms over the five periods ending at the times below, against the averaged prediction at the parameters reached then (the slaved s = r(x_m)/k_s filtered by (d/dt + λ)⁻¹ along the orbit), and at the end against the prediction with the error filter kept exactly, s = (d/dt + k_s)⁻¹r(x_m).

| prior | design | m | e_rms at t = 31, 500, end | predicted | ratios | exact filter at the end: predicted, ratio |
|---|---|---|---|---|---|---|
| matching | nearest atom, fixed prior | 16 | 1.804e-04, 1.813e-04, 1.815e-04 | 1.830e-04, 1.827e-04, 1.830e-04 | 0.986, 0.992, 0.992 | 1.814e-04, 1.0007 |
| matching | nearest atom, fixed prior | 64 | 9.155e-05, 9.224e-05, 9.237e-05 | 9.263e-05, 9.263e-05, 9.279e-05 | 0.988, 0.996, 0.995 | 9.229e-05, 1.0008 |
| matching | nearest atom, adaptive prior (η = ∞) | 16 | 1.705e-05, 1.704e-05, 1.704e-05 | 1.771e-05, 1.771e-05, 1.771e-05 | 0.963, 0.962, 0.962 | 1.701e-05, 1.0021 |
| matching | nearest atom, adaptive prior (η = ∞) | 64 | 1.869e-05, 1.869e-05, 1.869e-05 | 1.920e-05, 1.920e-05, 1.919e-05 | 0.973, 0.974, 0.974 | 1.865e-05, 1.0022 |
| matching | nearest atom, adaptive prior (η = σ₀) | 16 | 1.885e-05, 1.872e-05, 1.866e-05 | 1.968e-05, 1.958e-05, 1.951e-05 | 0.958, 0.956, 0.956 | 1.862e-05, 1.0022 |
| matching | nearest atom, adaptive prior (η = σ₀) | 64 | 1.826e-05, 1.826e-05, 1.825e-05 | 1.884e-05, 1.884e-05, 1.883e-05 | 0.969, 0.969, 0.969 | 1.821e-05, 1.0022 |
| matching | OT, fixed prior | 16 | 1.889e-05, 1.888e-05, 1.887e-05 | 1.937e-05, 1.936e-05, 1.935e-05 | 0.975, 0.975, 0.975 | 1.883e-05, 1.0022 |
| matching | OT, fixed prior | 64 | 1.889e-05, 1.888e-05, 1.887e-05 | 1.937e-05, 1.936e-05, 1.935e-05 | 0.975, 0.975, 0.975 | 1.883e-05, 1.0022 |
| σ-mod | σ-modification (ϑ⁰ = 0) | 16 | 5.545e-04, 5.524e-04, 5.529e-04 | 5.571e-04, 5.579e-04, 5.578e-04 | 0.995, 0.990, 0.991 | 5.513e-04, 1.0029 |
| σ-mod | σ-modification (ϑ⁰ = 0) | 64 | 5.467e-04, 5.447e-04, 5.453e-04 | 5.494e-04, 5.501e-04, 5.501e-04 | 0.995, 0.990, 0.991 | 5.437e-04, 1.0029 |
| misplaced | nearest atom, fixed prior | 16 | 1.385e-04, 1.399e-04, 1.402e-04 | 1.410e-04, 1.412e-04, 1.415e-04 | 0.982, 0.991, 0.991 | 1.401e-04, 1.0004 |
| misplaced | nearest atom, fixed prior | 64 | 9.472e-05, 9.588e-05, 9.611e-05 | 9.645e-05, 9.670e-05, 9.705e-05 | 0.982, 0.992, 0.990 | 9.611e-05, 1.0000 |
| misplaced | nearest atom, adaptive prior (η = ∞) | 16 | 1.705e-05, 1.704e-05, 1.704e-05 | 1.771e-05, 1.771e-05, 1.771e-05 | 0.963, 0.962, 0.962 | 1.701e-05, 1.0021 |
| misplaced | nearest atom, adaptive prior (η = ∞) | 64 | 1.869e-05, 1.869e-05, 1.869e-05 | 1.920e-05, 1.920e-05, 1.919e-05 | 0.973, 0.974, 0.974 | 1.865e-05, 1.0022 |
| misplaced | nearest atom, adaptive prior (η = σ₀) | 16 | 2.084e-05, 2.065e-05, 2.058e-05 | 2.222e-05, 2.210e-05, 2.202e-05 | 0.938, 0.935, 0.935 | 2.055e-05, 1.0016 |
| misplaced | nearest atom, adaptive prior (η = σ₀) | 64 | 1.782e-05, 1.773e-05, 1.772e-05 | 1.909e-05, 1.903e-05, 1.901e-05 | 0.934, 0.932, 0.932 | 1.768e-05, 1.0019 |
| misplaced | OT, fixed prior | 16 | 7.782e-05, 7.789e-05, 7.806e-05 | 7.862e-05, 7.877e-05, 7.899e-05 | 0.990, 0.989, 0.988 | 7.804e-05, 1.0003 |
| misplaced | OT, fixed prior | 64 | 7.757e-05, 7.765e-05, 7.782e-05 | 7.836e-05, 7.852e-05, 7.874e-05 | 0.990, 0.989, 0.988 | 7.780e-05, 1.0003 |

Numbers are produced by `experiments/exp_servo_cogging.py`; nothing is edited by hand.
