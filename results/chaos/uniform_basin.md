# Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate (`exp_uniform_basin.py`)

The real closed loop of Theorem J on the three-bump instance: k_e = 25, σ = σ₀/k_e with σ₀ = 0.15σ_crit = 0.00707, RK4 plant with the control and the reference at the stages and the parameters frozen over the step, Euler parameters, dt = 0.005, horizon 8000 time units (1273 periods, 320 slow units); the leak target refreshed every 200 steps; records every 5 time units. Fitted rates: least-squares slope of log d_∞ over the second half of the run.

## A. The matching prior at the matching masses (the lift is the exact representation, a constant orbit)

Prediction: λ_sep = 0.00707 (= σ₀ at zero residual), λ_coll = 0.00738, rate min{λ_coll, λ_sep}/k_e = 2.830e-04 per time unit. Every particle is displaced from its class atom by ρ in a random direction of the G⁻¹-unit sphere (the same displacements under both leaks); the Voronoi cells extend at least Δ_ν/2 = 3.13 from the atoms.

| leak | ρ | m | final counts by nearest atom | d_∞(0) → d_∞(end) | fitted rate of d_∞ (second half) | C = max_t d_∞(t)e^{μt}/d_∞(0) (μ the predicted rate) | W²(μ_Θ, ν): start → end | rms tracking error (last 10 periods) |
|---|---|---|---|---|---|---|---|---|
| nearest | 0.3 | 3 | (1, 1, 1) | 0.300 → 1.04e-02 | 3.484e-04 | 1.00 | 9.000e-02 → 8.69e-05 | 1.4e-05 |
| nearest | 0.3 | 15 | (5, 5, 5) | 0.300 → 3.63e-02 | 2.694e-04 | 1.16 | 9.000e-02 → 8.29e-04 | 1.1e-05 |
| nearest | 0.3 | 60 | (20, 20, 20) | 0.300 → 4.01e-02 | 2.802e-04 | 1.29 | 9.000e-02 → 8.71e-04 | 4.5e-06 |
| nearest | 0.3 | 240 | (80, 80, 80) | 0.300 → 3.42e-02 | 2.888e-04 | 1.12 | 9.000e-02 → 9.24e-04 | 1.7e-06 |
| nearest | 1 | 3 | (1, 1, 1) | 1.000 → 6.25e-02 | 2.948e-04 | 1.00 | 1.000e+00 → 2.29e-03 | 5.9e-05 |
| nearest | 1 | 15 | (5, 5, 5) | 1.000 → 1.32e-01 | 2.756e-04 | 1.27 | 1.000e+00 → 9.50e-03 | 2.4e-05 |
| nearest | 1 | 60 | (20, 20, 20) | 1.000 → 1.07e-01 | 2.967e-04 | 1.11 | 1.000e+00 → 9.18e-03 | 1.5e-05 |
| nearest | 1 | 240 | (80, 80, 80) | 1.000 → 1.11e-01 | 2.934e-04 | 1.14 | 1.000e+00 → 8.62e-03 | 9.8e-06 |
| nearest | 2.5 | 3 | (1, 1, 1) | 2.500 → 1.13e-01 | 4.062e-04 | 1.00 | 6.250e+00 → 8.83e-03 | 2.0e-04 |
| nearest | 2.5 | 15 | (5, 5, 5) | 2.500 → 2.51e-01 | 2.900e-04 | 1.03 | 6.250e+00 → 3.89e-02 | 6.2e-05 |
| nearest | 2.5 | 60 | (20, 20, 20) | 2.500 → 2.87e-01 | 2.942e-04 | 1.16 | 6.250e+00 → 3.61e-02 | 2.4e-05 |
| nearest | 2.5 | 240 | (80, 80, 80) | 2.500 → 2.79e-01 | 2.971e-04 | 1.14 | 6.250e+00 → 4.64e-02 | 3.5e-05 |
| nearest | 4.5 | 3 | (2, 0, 1) | 4.500 → 2.50e+00 | 1.508e-05 | 5.33 | 2.025e+01 → 1.23e+01 | 9.0e-03 |
| nearest | 4.5 | 15 | (6, 5, 4) | 4.500 → 1.16e+00 | 1.392e-04 | 2.47 | 1.866e+01 → 2.81e+00 | 1.2e-03 |
| nearest | 4.5 | 60 | (23, 19, 18) | 4.500 → 1.27e+00 | 1.341e-04 | 2.71 | 1.836e+01 → 2.30e+00 | 7.2e-04 |
| nearest | 4.5 | 240 | (87, 70, 83) | 4.500 → 8.93e-01 | 2.350e-04 | 1.91 | 1.832e+01 → 1.92e+00 | 5.4e-04 |
| ot | 0.3 | 3 | (1, 1, 1) | 0.300 → 1.04e-02 | 3.484e-04 | 1.00 | 9.000e-02 → 8.69e-05 | 1.4e-05 |
| ot | 0.3 | 15 | (5, 5, 5) | 0.300 → 3.63e-02 | 2.694e-04 | 1.16 | 9.000e-02 → 8.29e-04 | 1.1e-05 |
| ot | 0.3 | 60 | (20, 20, 20) | 0.300 → 4.01e-02 | 2.802e-04 | 1.29 | 9.000e-02 → 8.71e-04 | 4.5e-06 |
| ot | 0.3 | 240 | (80, 80, 80) | 0.300 → 3.42e-02 | 2.888e-04 | 1.12 | 9.000e-02 → 9.24e-04 | 1.7e-06 |
| ot | 1 | 3 | (1, 1, 1) | 1.000 → 6.25e-02 | 2.948e-04 | 1.00 | 1.000e+00 → 2.29e-03 | 5.9e-05 |
| ot | 1 | 15 | (5, 5, 5) | 1.000 → 1.32e-01 | 2.756e-04 | 1.27 | 1.000e+00 → 9.50e-03 | 2.4e-05 |
| ot | 1 | 60 | (20, 20, 20) | 1.000 → 1.07e-01 | 2.967e-04 | 1.11 | 1.000e+00 → 9.18e-03 | 1.5e-05 |
| ot | 1 | 240 | (80, 80, 80) | 1.000 → 1.11e-01 | 2.934e-04 | 1.14 | 1.000e+00 → 8.62e-03 | 9.8e-06 |
| ot | 2.5 | 3 | (1, 1, 1) | 2.500 → 1.13e-01 | 4.062e-04 | 1.00 | 6.250e+00 → 8.83e-03 | 2.0e-04 |
| ot | 2.5 | 15 | (5, 5, 5) | 2.500 → 2.51e-01 | 2.900e-04 | 1.03 | 6.250e+00 → 3.89e-02 | 6.2e-05 |
| ot | 2.5 | 60 | (20, 20, 20) | 2.500 → 2.87e-01 | 2.942e-04 | 1.16 | 6.250e+00 → 3.61e-02 | 2.4e-05 |
| ot | 2.5 | 240 | (80, 80, 80) | 2.500 → 2.79e-01 | 2.971e-04 | 1.14 | 6.250e+00 → 4.64e-02 | 3.5e-05 |
| ot | 4.5 | 3 | (1, 1, 1) | 4.500 → 1.57e-01 | 7.377e-04 | 1.88 | 2.025e+01 → 1.01e-02 | 2.1e-04 |
| ot | 4.5 | 15 | (5, 5, 5) | 4.500 → 9.20e-01 | 3.017e-04 | 2.11 | 1.866e+01 → 2.26e-01 | 2.1e-04 |
| ot | 4.5 | 60 | (20, 20, 20) | 4.500 → 7.55e-01 | 2.690e-04 | 1.61 | 1.836e+01 → 1.65e-01 | 1.2e-04 |
| ot | 4.5 | 240 | (80, 80, 80) | 4.500 → 8.77e-01 | 2.474e-04 | 1.87 | 1.832e+01 → 1.56e-01 | 1.5e-04 |

Inside the cells (ρ = 0.3, 1.0, 2.5; both leaks; all m): fitted rates 2.694e-04–4.062e-04 against the predicted 2.830e-04 (ratios 0.952–1.435); constants C = 1.00–1.29. At ρ = 4.5: the nearest-atom leak ends with counts (2, 0, 1) (m = 3), (6, 5, 4) (m = 15), (23, 19, 18) (m = 60), (87, 70, 83) (m = 240) and the OT leak with (1, 1, 1), (5, 5, 5), (20, 20, 20), (80, 80, 80); W²(μ_Θ, ν) at the end 1.23e+01 / 1.01e-02, 2.81e+00 / 2.26e-01, 2.30e+00 / 1.65e-01, 1.92e+00 / 1.56e-01 (nearest / OT).

## B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (non-zero residual), masses q

Θ*(q) by Newton (residual 4.3e-17); prediction λ_sep = 0.00567, λ_coll = 0.00740, rate 2.270e-04 per time unit. The reference is the 3-atom loop from Θ*(q); d_∞ is measured against the lift of the reference trajectory.

| ρ | m | d_∞(0) → d_∞(end) | fitted rate | C | rms tracking error (m-system / reference) |
|---|---|---|---|---|---|
| 0.3 | 15 | 0.300 → 6.02e-02 | 2.184e-04 | 1.23 | 2.55e-04 / 2.55e-04 |
| 1 | 15 | 1.000 → 1.36e-01 | 2.632e-04 | 1.08 | 2.55e-04 / 2.55e-04 |
| 0.3 | 60 | 0.300 → 5.08e-02 | 2.339e-04 | 1.14 | 2.55e-04 / 2.55e-04 |
| 1 | 60 | 1.000 → 1.48e-01 | 2.290e-04 | 1.08 | 2.56e-04 / 2.55e-04 |
| 0.3 | 240 | 0.300 → 4.82e-02 | 2.242e-04 | 1.04 | 2.55e-04 / 2.55e-04 |
| 1 | 240 | 1.000 → 1.53e-01 | 2.336e-04 | 1.12 | 2.56e-04 / 2.55e-04 |

Fitted rates 2.184e-04–2.632e-04 against the predicted 2.270e-04 (ratios 0.962–1.159).

Numbers are produced by `experiments/exp_uniform_basin.py`; nothing is edited by hand.
