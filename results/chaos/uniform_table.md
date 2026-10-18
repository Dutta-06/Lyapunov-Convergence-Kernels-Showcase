# Theorem D″: uniform-in-time propagation of chaos under a dominant leak (`exp_chaos_uniform.py`)

n = 1, k_e = 25, G = diag(0.2, 0.05, 0.05), saturated features with a_max = 4.0, c_max = 3.0, s ∈ [0.3, 2.0], prior ϑ⁰ = 0, τ = 0; the m-system shares its particles with an M = 2048 reference; T = 80 s, dt = 0.001; 8 seeds. Lipschitz constants (numerical suprema over a sampled parameter box and |x| ≤ 2.5; the runs stay within |x| ≤ 1.75): L_f = 2.3, L_x = 7.98, L_ϑ = 7.70, L₂ = 31.58, L_xϑ = 41.08; measured E_X = sup‖ē‖ = 0.0522 on the reference runs. κ_e = k_e − L_f − L_x = 14.72; L_b′ = λ_G(L_ϑ + E_X L_xϑ) = 1.969; the theorem applies for σ > σ_cert = λ_G E_X L₂ + 2L_b′L_ϑ/κ_e = 2.391 (the runs below used a provisional E_X = 0.1 to set the leaks: σ = 4.656, 0.310, 0.020).

## σ = 4.656 (certified: det A = 33.340); particle spread of the reference at T (max_i ‖ϑ_i − ϑ̄‖, seed mean): 4.27e-17

| m | max D_m on (10, 20] | on (20, 40] | on (40, 80] | time average on (40, 80] | bound (K_Φκ_ϑ + 2K_ΦL_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 8.078e-26 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 2.479e-01 |
| 64 | 1.155e-25 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 1.240e-01 |
| 256 | 6.668e-26 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 6.198e-02 |

Slope of the (40, 80] window maximum against m: nan (Theorem D″: −0.5). Ratio of the (40, 80] to the (10, 20] window maxima: 0.00, 0.00, 0.00 (uniform in time means no growth).

## σ = 0.310 (below the threshold); particle spread of the reference at T (max_i ‖ϑ_i − ϑ̄‖, seed mean): 2.09e-10

| m | max D_m on (10, 20] | on (20, 40] | on (40, 80] | time average on (40, 80] | bound (K_Φκ_ϑ + 2K_ΦL_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 4.987e-04 | 2.429e-05 | 1.028e-07 | 8.864e-09 | nan |
| 64 | 4.040e-04 | 1.975e-05 | 8.267e-08 | 7.121e-09 | nan |
| 256 | 2.073e-04 | 1.017e-05 | 4.262e-08 | 3.669e-09 | nan |

Slope of the (40, 80] window maximum against m: -0.32 (Theorem D″: −0.5). Ratio of the (40, 80] to the (10, 20] window maxima: 0.00, 0.00, 0.00 (uniform in time means no growth).

## σ = 0.020 (below the threshold); particle spread of the reference at T (max_i ‖ϑ_i − ϑ̄‖, seed mean): 2.57e+00

| m | max D_m on (10, 20] | on (20, 40] | on (40, 80] | time average on (40, 80] | bound (K_Φκ_ϑ + 2K_ΦL_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 9.928e-03 | 1.053e-02 | 9.874e-03 | 7.770e-03 | nan |
| 64 | 5.579e-03 | 6.408e-03 | 6.513e-03 | 5.467e-03 | nan |
| 256 | 2.785e-03 | 3.198e-03 | 3.210e-03 | 2.728e-03 | nan |

Slope of the (40, 80] window maximum against m: -0.41 (Theorem D″: −0.5). Ratio of the (40, 80] to the (10, 20] window maxima: 0.99, 1.17, 1.15 (uniform in time means no growth).

Reading: the condition κ_ϑ > 0 makes the particle map a contraction uniformly in time, so it also synchronises the particles of one system: at the certified leak the reference collapses to a single kernel (particle spread at the rounding level) and the coupled distance vanishes identically — propagation of chaos is uniform in time but trivial. At σ_cert/8 the particles still synchronise, more slowly, and the coupled distance decays with them. At σ = 0.02 the network stays spread, the coupled distance settles at the m^{−1/2} level and does not grow over the horizon: uniform in time empirically, beyond what the theorem certifies. The learning regime is the open part of §8, branch 3.

Numbers are produced by `experiments/exp_chaos_uniform.py`; nothing is edited by hand.
