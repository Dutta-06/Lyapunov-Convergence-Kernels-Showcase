# Corollary D‴.14: the real closed loop with the OT leak — the orbit is the lift of the weighted K₀-atom loop's, the same for every commensurate size (`exp_ot_real_loop.py`)

The real closed loop of Theorem J on the three-bump instance (k_e = 25, σ = σ₀/k_e with σ₀ = 0.15σ_crit = 0.00707, RK4 plant with the control held over the step, Euler parameters, dt = 0.001), 100 time units = 15.9 periods = 4.0 slow units; statistics over the last 5 periods on the 400-point grid of the reference orbit. Starts: the contracted clouds at the end of `exp_ot_leak.py`'s averaged runs (commensurate set, m = 15, 60, M = 510; the OT runs' clouds carry the prior's masses q = (1/3, 1/3, 1/3), the nearest-atom runs' clouds their conserved sampled weights), and the collapsed lifts: the three atoms of the exact representation for the OT leak (the weighted 3-atom loop with weights q, itself the m = 3 loop), and the replicated weighted stationary configurations Θ*(w^m) for the nearest-atom leak. Seeds 0, 1; the leak target (the nearest atom, or the barycentric projection of the vertex plan) recomputed every 100 steps.

## 1. The OT leak: every size on the same orbit

| seed | system | cluster weights at the start | ‖f̂_m − f̂_M‖²_{L²(orbit)} (M = 510) | ‖f̂_m − f̂_lift‖² | fit error along the orbit | rms tracking error | W²_{G⁻¹}(μ_m, μ_M) at the end | W² to the lift |
|---|---|---|---|---|---|---|---|---|
| 0 | m15 | (0.3333, 0.3333, 0.3333) | 2.743e-08 | 3.557e-08 | 3.748e-04 | 3.940e-05 | 1.390e-03 | 9.875e-04 |
| 0 | m60 | (0.3333, 0.3333, 0.3333) | 7.452e-09 | 5.198e-10 | 2.124e-04 | 3.927e-05 | 3.346e-04 | 1.647e-04 |
| 0 | m510 | (0.3333, 0.3333, 0.3333) | 0.000e+00 | 1.167e-08 | 2.634e-04 | 4.002e-05 | 1.226e-15 | 4.160e-04 |
| 0 | lift (m = 3) | (0.3333, 0.3333, 0.3333) | 1.167e-08 | 0.000e+00 | 2.111e-04 | 3.912e-05 | 4.160e-04 | 4.737e-15 |
| 1 | m15 | (0.3333, 0.3333, 0.3333) | 8.577e-09 | 8.591e-09 | 2.568e-04 | 3.949e-05 | 1.188e-03 | 1.088e-03 |
| 1 | m60 | (0.3333, 0.3333, 0.3333) | 8.525e-09 | 4.938e-09 | 2.396e-04 | 3.941e-05 | 8.627e-04 | 7.775e-04 |
| 1 | m510 | (0.3333, 0.3333, 0.3333) | 0.000e+00 | 4.249e-09 | 2.302e-04 | 3.958e-05 | 1.811e-15 | 1.989e-04 |
| 1 | lift (m = 3) | (0.3333, 0.3333, 0.3333) | 4.249e-09 | 0.000e+00 | 2.111e-04 | 3.912e-05 | 1.989e-04 | 4.737e-15 |

## 2. The nearest-atom leak: every size on the orbit of its own weights

| seed | m | conserved weights w^m | ‖f̂_m − f̂_M‖² | ‖f̂_m − f̂_{lift(w^m)}‖² | ‖f̂_{lift(w^m)} − f̂_{lift(w^M)}‖² | fit error along the orbit (m-system / its lift) | rms tracking error (m-system / its lift) | W²(μ_m, μ_M) | W² to its lift |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 15 | (0.4000, 0.3333, 0.2667) | 1.823e-03 | 3.321e-06 | 1.899e-03 | 4.428e-02 / 4.549e-02 | 1.162e-03 / 1.193e-03 | 4.883e+00 | 1.252e-02 |
| 0 | 60 | (0.3167, 0.3500, 0.3333) | 1.064e-04 | 1.468e-11 | 1.064e-04 | 1.337e-02 / 1.337e-02 | 3.311e-04 / 3.310e-04 | 1.569e+00 | 3.067e-07 |
| 0 | 510 | (0.3098, 0.3157, 0.3745) | 0.000e+00 | 2.522e-11 | 0.000e+00 | 2.274e-02 / 2.274e-02 | 5.656e-04 / 5.656e-04 | 1.533e-15 | 1.812e-06 |
| 1 | 15 | (0.4000, 0.2000, 0.4000) | 2.467e-03 | 4.154e-08 | 2.473e-03 | 7.470e-02 / 7.477e-02 | 1.914e-03 / 1.915e-03 | 5.547e+00 | 1.800e-03 |
| 1 | 60 | (0.4000, 0.2500, 0.3500) | 1.074e-03 | 1.838e-06 | 1.118e-03 | 4.451e-02 / 4.537e-02 | 1.150e-03 / 1.171e-03 | 3.734e+00 | 1.074e-02 |
| 1 | 510 | (0.3196, 0.3255, 0.3549) | 0.000e+00 | 9.607e-12 | 0.000e+00 | 1.252e-02 / 1.252e-02 | 3.033e-04 / 3.032e-04 | 1.811e-15 | 5.168e-07 |

Reading. Under the OT leak the learned functions of m = 15, 60 and M = 510 agree along the orbit to ‖f̂_m − f̂_M‖² = 7.5e-09–2.7e-08 and with the 3-atom lift to 5.2e-10–3.6e-08 (the within-class spread of the clouds, second order); under the nearest-atom leak they disagree by 1.1e-04–2.5e-03, and this disagreement is the distance between the lifted orbits of their conserved weights to a ratio of 0.960–1.000 (each m-system within 3.3e-06 of its own lift). The rms tracking error is 3.9e-05–4.0e-05 under the OT leak for every size against 3.0e-04–1.9e-03 under the nearest-atom leak — the weight bias of the fit is a tracking error, and the mass constraint removes it for every size at once.

Numbers are produced by `experiments/exp_ot_real_loop.py`; nothing is edited by hand.
