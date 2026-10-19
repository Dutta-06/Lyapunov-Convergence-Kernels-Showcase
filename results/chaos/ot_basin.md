# Theorem D‴.21: a measure-level basin for the OT leak — the Polyak–Łojasiewicz inequality along the runs (`exp_ot_basin.py`)

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = 0.00707 (2σ₀ = 0.01415); sizes m = 15, 60, 240 and M = 510 from the warm clouds of Proposition D‴.12 (width 0.02 around the exact representation; the m-system the first m rows); 4 seeds; 600 slow units, midpoint rule dt = 0.02, leak target every 10 steps. E = 𝒥(μ_Θ) − 𝒥* with 𝒥 = L̄ + (σ₀/2)W²(·, ν) the OT functional, evaluated along both leaks' trajectories; D_a the class-respecting distance to μ*; the observed PL ratio |∇J_a|²_G/(mE), which along the OT flow is −d log E/ds; λ the theorem's rate from the plan-segment curvatures (9 points per segment, every 25 slow units). Entries: mean over the seeds (min–max).

## Prior: matching (centres shifted by 0)

μ*: Newton residual 2.3e-17, stationarity of the atoms 2.7e-17, 𝒥* = 5.0447e-33; diagonal plan: ε = 0.000 < Δ_ν/2 = 3.131. Identities (E), (P) along every run: largest relative residual 5.4e-15 (E) and 4.1e-15 (P).

| leak | m | E at s = 0, 25, 100, 300, 600 | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |
|---|---|---|---|---|---|
| ot | 15 | 3.88e-02, 1.51e-02, 8.00e-03, 6.29e-04, 8.08e-06 | 5.24e+00, 1.78e+00, 2.28e-03 | 5.24e+00, 1.78e+00, 2.28e-03 | 0–0 |
| ot | 60 | 1.44e-02, 7.11e-03, 2.84e-03, 1.59e-04, 2.25e-06 | 2.75e+00, 7.19e-01, 6.34e-04 | 2.75e+00, 7.19e-01, 6.34e-04 | 0–0 |
| ot | 240 | 5.34e-03, 3.17e-03, 1.41e-03, 7.59e-05, 1.06e-06 | 1.17e+00, 3.45e-01, 2.99e-04 | 1.17e+00, 3.45e-01, 2.99e-04 | 0–0 |
| ot | 510 | 4.26e-03, 2.60e-03, 1.05e-03, 5.85e-05, 8.29e-07 | 9.94e-01, 2.72e-01, 2.34e-04 | 9.94e-01, 2.72e-01, 2.34e-04 | 0–0 |
| nearest | 15 | 3.88e-02, 2.22e-02, 1.84e-02, 1.81e-02, 1.79e-02 | 5.24e+00, 4.80e+00, 4.71e+00 | 5.24e+00, 4.80e+00, 4.71e+00 | 0–2 |
| nearest | 60 | 1.44e-02, 1.06e-02, 9.39e-03, 9.20e-03, 9.04e-03 | 2.75e+00, 2.57e+00, 2.49e+00 | 2.75e+00, 2.57e+00, 2.49e+00 | 0–8 |
| nearest | 240 | 5.34e-03, 4.47e-03, 4.33e-03, 4.35e-03, 4.37e-03 | 1.17e+00, 1.19e+00, 1.21e+00 | 1.17e+00, 1.19e+00, 1.21e+00 | 0–10 |
| nearest | 510 | 4.26e-03, 3.69e-03, 3.59e-03, 3.61e-03, 3.62e-03 | 9.94e-01, 1.00e+00, 1.01e+00 | 9.94e-01, 1.00e+00, 1.01e+00 | 2–17 |

The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).

| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |
|---|---|---|---|---|---|---|---|---|
| 15 | 150, 25, 25, 300 | 1.3e-01, 4.9e-01, 4.8e-01, 2.0e-02 | 1.00, 1.00 | 0.01417 | 0.01396 | 0.00703 | -3.5e-19, 9.4e-05 | 3.1e-04 ≥ 3.7e-06; 1.2e-04 ≥ 4.1e-06; 1.7e-05 ≥ 1.9e-06; 4.4e-03 ≥ 2.3e-05 |
| 60 | 0, 25, 125, 150 | 1.0e+00, 5.0e-01, 1.7e-01, 1.0e-01 | 1.00, 1.00 | 0.01417 | 0.01404 | 0.00705 | -3.5e-19, 5.2e-05 | 1.6e-06 ≥ 6.2e-07; 4.3e-05 ≥ 2.9e-06; 3.0e-05 ≥ 1.8e-06; 9.3e-05 ≥ 3.7e-06 |
| 240 | 150, 0, 0, 0 | 1.5e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.00, 1.00 | 0.01417 | 0.01401 | 0.00704 | -3.5e-19, 7.0e-05 | 9.4e-05 ≥ 2.4e-06; 9.8e-06 ≥ 9.8e-07; 5.7e-06 ≥ 6.1e-07; 6.4e-07 ≥ 2.9e-07 |
| 510 | 125, 0, 0, 0 | 1.9e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.00, 1.00 | 0.01416 | 0.01406 | 0.00705 | -3.5e-19, 4.6e-05 | 3.4e-05 ≥ 1.6e-06; 7.2e-06 ≥ 7.5e-07; 2.1e-07 ≥ 1.7e-07; 2.0e-06 ≥ 8.4e-07 |

The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the nearest-atom energy; its drift under the nearest-atom leak and under the OT leak).

| m | sampled counts (seeds) | E(600) | E of the replicated configuration of the sampled weights | its drift: nearest-atom / OT leak |
|---|---|---|---|---|
| 15 | (6, 5, 4); (6, 3, 6); (5, 4, 6); (8, 4, 3) | 1.08e-02; 2.10e-02; 8.74e-03; 3.09e-02 | 1.18e-02; 2.11e-02; 8.74e-03; 3.67e-02 | 5e-17 / 2.3e-02; 1e-17 / 4.4e-02; 5e-17 / 3.6e-02; 6e-14 / 4.4e-02 |
| 60 | (19, 21, 20); (24, 15, 21); (19, 18, 23); (25, 17, 18) | 3.10e-03; 1.35e-02; 7.21e-03; 1.24e-02 | 3.10e-03; 1.38e-02; 7.22e-03; 1.39e-02 | 3e-15 / 4.4e-02; 6e-17 / 4.4e-02; 3e-17 / 3.6e-02; 6e-17 / 4.4e-02 |
| 240 | (71, 76, 93); (77, 76, 87); (76, 81, 83); (82, 79, 79) | 8.87e-03; 4.31e-03; 2.85e-03; 1.42e-03 | 8.89e-03; 4.31e-03; 2.85e-03; 1.42e-03 | 1e-17 / 3.6e-02; 2e-14 / 3.6e-02; 2e-15 / 4.4e-02; 4e-17 / 4.4e-02 |
| 510 | (158, 161, 191); (163, 166, 181); (170, 167, 173); (168, 184, 158) | 6.37e-03; 3.33e-03; 8.10e-04; 3.95e-03 | 6.38e-03; 3.33e-03; 8.11e-04; 3.95e-03 | 7e-17 / 3.6e-02; 2e-15 / 3.6e-02; 5e-17 / 3.6e-02; 1e-14 / 4.4e-02 |

The a-priori basin of Theorem D‴.21(iv): sup λ_max R^G[μ*] over the box of the runs = 5.8500e-17, so κ = σ₀ − sup = 0.00707 and E₀ = κ²/(2(L₂^G)²) = 2.03e-07 (L₂^G = 11.111).

## Prior: shifted (centres shifted by 0.15)

μ*: Newton residual 4.3e-17, stationarity of the atoms 5.5e-17, 𝒥* = 4.6447e-04; diagonal plan: ε = 0.389 < Δ_ν/2 = 3.162. Identities (E), (P) along every run: largest relative residual 5.6e-13 (E) and 4.2e-15 (P).

| leak | m | E at s = 0, 25, 100, 300, 600 | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |
|---|---|---|---|---|---|
| ot | 15 | 3.98e-02, 1.58e-02, 8.36e-03, 4.08e-04, 5.21e-06 | 5.25e+00, 1.71e+00, 1.51e-03 | 5.25e+00, 1.71e+00, 1.51e-03 | 0–0 |
| ot | 60 | 1.46e-02, 7.26e-03, 2.86e-03, 1.52e-04, 2.34e-06 | 2.76e+00, 7.04e-01, 6.84e-04 | 2.76e+00, 7.04e-01, 6.84e-04 | 0–0 |
| ot | 240 | 5.03e-03, 2.88e-03, 1.32e-03, 1.26e-04, 2.10e-06 | 1.16e+00, 3.76e-01, 5.96e-04 | 1.16e+00, 3.76e-01, 5.96e-04 | 0–0 |
| ot | 510 | 4.11e-03, 2.44e-03, 9.94e-04, 8.64e-05, 1.46e-06 | 9.93e-01, 2.89e-01, 4.17e-04 | 9.93e-01, 2.89e-01, 4.17e-04 | 0–0 |
| nearest | 15 | 3.98e-02, 2.31e-02, 1.93e-02, 1.89e-02, 1.86e-02 | 5.25e+00, 4.80e+00, 4.69e+00 | 5.25e+00, 4.80e+00, 4.69e+00 | 0–2 |
| nearest | 60 | 1.46e-02, 1.07e-02, 9.56e-03, 9.38e-03, 9.18e-03 | 2.76e+00, 2.57e+00, 2.49e+00 | 2.76e+00, 2.57e+00, 2.49e+00 | 0–8 |
| nearest | 240 | 5.03e-03, 4.18e-03, 4.04e-03, 4.06e-03, 4.07e-03 | 1.16e+00, 1.19e+00, 1.21e+00 | 1.16e+00, 1.19e+00, 1.21e+00 | 0–14 |
| nearest | 510 | 4.11e-03, 3.54e-03, 3.44e-03, 3.45e-03, 3.46e-03 | 9.93e-01, 1.00e+00, 1.01e+00 | 9.93e-01, 1.00e+00, 1.01e+00 | 2–23 |

The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).

| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |
|---|---|---|---|---|---|---|---|---|
| 15 | 150, 125, 25, 275 | 9.5e-02, 1.5e-01, 4.6e-01, 1.6e-02 | 1.67, 1.64 | 0.01367 | 0.00766 | 0.00561 | 1.4e-03, 1.5e-03 | 1.7e-03 ≥ 2.3e-06; 1.6e-03 ≥ 4.7e-06; 3.8e-04 ≥ 3.4e-06; 1.3e-02 ≥ 1.0e-05 |
| 60 | 0, 25, 175, 125 | 1.0e+00, 5.1e-01, 1.0e-01, 1.3e-01 | 1.72, 1.69 | 0.01389 | 0.00778 | 0.00566 | 1.4e-03, 1.4e-03 | 5.0e-05 ≥ 8.5e-07; 6.6e-04 ≥ 2.6e-06; 7.6e-04 ≥ 3.6e-06; 9.2e-04 ≥ 2.3e-06 |
| 240 | 200, 25, 0, 0 | 1.0e-01, 5.8e-01, 1.0e+00, 1.0e+00 | 1.68, 1.64 | 0.01384 | 0.00778 | 0.00566 | 1.4e-03, 1.4e-03 | 1.2e-03 ≥ 4.9e-06; 3.4e-04 ≥ 2.1e-06; 2.0e-04 ≥ 1.2e-06; 3.1e-05 ≥ 2.2e-07 |
| 510 | 175, 0, 0, 0 | 1.2e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.61, 1.58 | 0.01381 | 0.00779 | 0.00567 | 1.4e-03, 1.4e-03 | 6.7e-04 ≥ 3.3e-06; 2.5e-04 ≥ 1.6e-06; 8.4e-06 ≥ 3.3e-07; 6.8e-05 ≥ 6.4e-07 |

The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the nearest-atom energy; its drift under the nearest-atom leak and under the OT leak).

| m | sampled counts (seeds) | E(600) | E of the replicated configuration of the sampled weights | its drift: nearest-atom / OT leak |
|---|---|---|---|---|
| 15 | (6, 5, 4); (6, 3, 6); (5, 4, 6); (8, 4, 3) | 1.22e-02; 2.11e-02; 8.30e-03; 3.30e-02 | 1.30e-02; 2.12e-02; 8.31e-03; 3.95e-02 | 4e-17 / 2.3e-02; 2e-17 / 4.4e-02; 6e-17 / 3.6e-02; 5e-14 / 4.4e-02 |
| 60 | (19, 21, 20); (24, 15, 21); (19, 18, 23); (25, 17, 18) | 2.96e-03; 1.38e-02; 6.73e-03; 1.32e-02 | 2.96e-03; 1.41e-02; 6.74e-03; 1.49e-02 | 2e-15 / 4.4e-02; 3e-17 / 4.4e-02; 3e-17 / 3.6e-02; 4e-17 / 4.4e-02 |
| 240 | (71, 76, 93); (77, 76, 87); (76, 81, 83); (82, 79, 79) | 8.16e-03; 4.00e-03; 2.62e-03; 1.52e-03 | 8.18e-03; 4.01e-03; 2.62e-03; 1.53e-03 | 4e-17 / 3.6e-02; 2e-14 / 3.6e-02; 2e-15 / 4.4e-02; 7e-17 / 4.4e-02 |
| 510 | (158, 161, 191); (163, 166, 181); (170, 167, 173); (168, 184, 158) | 5.87e-03; 3.06e-03; 7.73e-04; 4.13e-03 | 5.89e-03; 3.06e-03; 7.73e-04; 4.13e-03 | 1e-16 / 3.6e-02; 3e-15 / 3.6e-02; 1e-16 / 3.6e-02; 1e-14 / 4.4e-02 |

The a-priori basin of Theorem D‴.21(iv): sup λ_max R^G[μ*] over the box of the runs = 4.3487e-03, so κ = σ₀ − sup = 0.00273 and E₀ = κ²/(2(L₂^G)²) = 3.01e-08 (L₂^G = 11.111).

## The measure-level basin: particles moved into another cell (matching prior)

The lift of q with k particles of class 0 moved 0.6 of the way to atom 1 (into cell 1); the same fractions p = k/m at m = 30 and 300. At s = 0: whether the plan keeps the moved particles in class 0, E, D_π, the current-residual curvature κ^Θ along the plan segments, the margin α and the rate λ (κ* = 0 at the matching prior).

| m | k | p | kept in class 0 | E | D_π | κ^Θ | α | λ | observed ratio |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 1 | 0.0333 | True | 3.346e-03 | 6.229e-01 | 9.339e-03 | +0.00240 | 0.00163 | 0.03141 |
| 30 | 2 | 0.0667 | True | 8.976e-03 | 1.246e+00 | 1.868e-02 | -0.00227 | 0.00000 | 0.03763 |
| 30 | 3 | 0.1000 | True | 1.689e-02 | 1.869e+00 | 2.802e-02 | -0.00693 | 0.00000 | 0.03952 |
| 30 | 5 | 0.1667 | True | 3.957e-02 | 3.114e+00 | 4.670e-02 | -0.01627 | 0.00000 | 0.03781 |
| 300 | 10 | 0.0333 | True | 3.346e-03 | 6.229e-01 | 9.339e-03 | +0.00240 | 0.00163 | 0.03141 |
| 300 | 20 | 0.0667 | True | 8.976e-03 | 1.246e+00 | 1.868e-02 | -0.00227 | 0.00000 | 0.03763 |
| 300 | 30 | 0.1000 | True | 1.689e-02 | 1.869e+00 | 2.802e-02 | -0.00693 | 0.00000 | 0.03952 |
| 300 | 50 | 0.1667 | True | 3.957e-02 | 3.114e+00 | 4.670e-02 | -0.01627 | 0.00000 | 0.03781 |

One particle moved, 300 slow units (midpoint rule, ds = 0.5, plan every step): 

- m = 30, ot: E 3.346e-03 → 2.546e-05 (a-posteriori bound 1.017e-04); observed ratio at the end 0.01415; the moved particle labelled 0 by the OT plan, in cell 0 at the end.
- m = 30, nearest: E 3.346e-03 → 5.719e-03 (E is the OT functional's excess, not this leak's Lyapunov function); observed ratio at the end 0.01379; the moved particle labelled 0 by the OT plan, in cell 1 at the end.
- m = 300, ot: E 2.317e-04 → 3.093e-06 (a-posteriori bound 3.705e-06); observed ratio at the end 0.01415; the moved particle labelled 0 by the OT plan, in cell 0 at the end.
- m = 300, nearest: E 2.317e-04 → 5.558e-04 (E is the OT functional's excess, not this leak's Lyapunov function); observed ratio at the end 0.01413; the moved particle labelled 0 by the OT plan, in cell 1 at the end.

Numbers are produced by `experiments/exp_ot_basin.py`; nothing is edited by hand.
