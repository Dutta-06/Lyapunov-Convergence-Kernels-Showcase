# Proposition D‴.47: the free law's branches (`exp_mirror_branch.py`)

The free mass–amplitude law (Proposition D‴.45) with particles on branches (Q < −c), represented without cancellation near their branch points (`mirror_free.ParticleSet`).

## A. Identities

dη/dlog|w| against a/(γa² + Gβa′²) at 40 random states (free and on branches; finite differences): largest relative difference 2.33e-9.

The temperature at the branch point against 1/(γa₀):

| α₀ | a/(γa² + Gβa′²) | 1/(γa₀) |
|---|---|---|
| 0.05 | 200.1666389 | 200.1666389 |
| 0.3 | 34.3273843 | 34.3273843 |
| 1 | 13.13035285 | 13.13035285 |
| 2 | 10.37314721 | 10.37314721 |

B(w̄, w) on a branch as w̄ → 0 (δ̄ the distance to the branch point) against β/(γm):

| δ̄ | w̄ | B(w̄, w) | β/(γm) |
|---|---|---|---|
| 1e-06 | 1.24e-9 | 0.00185989952352 | 0.00186030524931 |
| 1e-12 | 1.24e-15 | 0.00186030524842 | 0.00186030524931 |
| 1e-30 | 1.24e-33 | 0.00186030524931 | 0.00186030524931 |
| 1e-100 | 1.24e-103 | 0.00186030524931 | 0.00186030524931 |

The fraction of starts on branches, α(0) uniform on [−0.3, 0.3]: 0.5248 at β(0) = 0.001, 0.9512 at β(0) = 1e-05, 0.9951 at β(0) = 1e-07, 0.9996 at β(0) = 1e-09.

A ± pair (α(0) = ±0.3, u = w₊(0) = 2.43e-9, a₀ = 0.29131): the atom's coefficient W(Δ) against 2u sinh(γa₀Δ):

| Δ | W/u | relative difference |
|---|---|---|
| 1 | 0.0583 | 1.98e-5 |
| 10 | 0.591 | 2.14e-5 |
| 30 | 1.98 | 3.58e-5 |
| 60 | 5.57 | 9.70e-5 |
| 90 | 13.7 | 2.53e-4 |
| 120 | 33 | 6.32e-4 |

The critical curve Q = −c: η·w^{1/3} (η by quadrature, up to its constant) against −(a_max²/(2Gc))(2c/(ma_max²))^{1/3} = -3.2183: -3.20801 at α = 0.01, -3.21729 at α = 0.001, -3.2182 at α = 0.0001.

## B. The limit on random Gaussian designs (6 × 12, the target on four atoms with signs +, −, +, −; masses 10⁻³, amplitudes 0.2–0.4)

The flow (LSODA) to t = 3·10⁷ against the Bregman projection of its start onto S (computed independently); |y − ŷ| the sign-constrained residual, |Z| the coordinates vanishing on S; J − min_K J at the end; the largest mass on Z at the end; B(w̄, w(t))'s largest relative increase; the extinction rates (fitted where 10⁻²⁵⁰ < β < 10⁻¹²) against γa₀|g(∞)|:

| start | seed | branches | |y − ŷ| | |Z| | projection (it, residual) | max |w(T) − w̄| | J − J* | max β on Z | B's increase | rates: largest relative error (n) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive | 0 | 12 | 0.703 | 8 | 9, 3.46e-14 | 1.18e-14 | 5.00e-16 | 0 | 5.78e-17 | 0.0112 (8) |
| positive | 1 | 12 | 0.0443 | 7 | 11, 6.01e-14 | 5.17e-14 | 6.07e-18 | 7.22e-104 | 2.07e-17 | 0.101 (7) |
| positive | 2 | 12 | 4.58e-16 | 0 | 12, 5.30e-14 | 1.00e-13 | 1.55e-29 | 0 | 7.71e-17 | — (0) |
| positive | 3 | 12 | 0.304 | 7 | 11, 3.95e-14 | 2.11e-14 | -9.71e-17 | 0 | 3.61e-17 | 0.0111 (7) |
| positive | 4 | 12 | 3.72e-16 | 0 | 12, 9.88e-15 | 8.44e-15 | 1.86e-29 | 0 | 4.23e-17 | — (0) |
| positive | 5 | 12 | 0.806 | 8 | 11, 3.41e-15 | 4.81e-15 | -2.22e-16 | 0 | 6.81e-17 | 0.00952 (8) |
| positive | 6 | 12 | 4.08e-16 | 0 | 12, 5.18e-15 | 1.22e-10 | 3.91e-30 | 0 | 7.19e-17 | — (0) |
| positive | 7 | 12 | 1.7 | 8 | 11, 8.78e-13 | 4.11e-13 | 0 | 0 | 5.13e-17 | 0.00453 (8) |
| random signs | 0 | 12 | 2.25 | 8 | 12, 1.57e-14 | 8.88e-15 | 1.78e-15 | 0 | 7.30e-17 | 0.0834 (8) |
| random signs | 1 | 12 | 7.13e-16 | 0 | 13, 4.98e-15 | 4.73e-14 | 3.08e-31 | 0 | 3.99e-17 | — (0) |
| random signs | 2 | 12 | 1.12e-15 | 0 | 15, 1.04e-12 | 3.30e-13 | 4.28e-29 | 0 | 2.08e-17 | — (0) |
| random signs | 3 | 12 | 0.556 | 7 | 11, 2.53e-14 | 5.22e-15 | -1.67e-16 | 0 | 1.29e-16 | 0.0311 (7) |
| random signs | 4 | 12 | 2.99e-16 | 0 | 11, 1.07e-14 | 1.11e-10 | 8.11e-30 | 0 | 6.50e-17 | — (0) |
| random signs | 5 | 12 | 4.78e-16 | 0 | 10, 3.14e-15 | 3.55e-12 | 1.85e-30 | 0 | 1.03e-16 | — (0) |
| random signs | 6 | 12 | 3.47e-16 | 0 | 13, 1.10e-14 | 1.73e-13 | 1.76e-30 | 0 | 4.74e-17 | — (0) |
| random signs | 7 | 12 | 3.51e-16 | 0 | 12, 6.49e-15 | 3.28e-12 | 2.62e-30 | 0 | 5.43e-17 | — (0) |
| mixed | 0 | 6 | 9.61e-16 | 0 | 14, 1.29e-13 | 3.11e-14 | 4.91e-30 | 0 | 1.09e-16 | — (0) |
| mixed | 1 | 6 | 1.57e-16 | 0 | 13, 8.55e-14 | 1.48e-10 | 3.93e-30 | 0 | 1.53e-16 | — (0) |
| mixed | 2 | 6 | 1.02e-15 | 0 | 16, 4.96e-14 | 3.49e-11 | 4.49e-30 | 0 | 6.35e-17 | — (0) |
| mixed | 3 | 6 | 3.25e-16 | 0 | 12, 1.22e-14 | 1.08e-11 | 2.78e-30 | 0 | 1.00e-16 | — (0) |
| mixed | 4 | 6 | 3.51e-16 | 0 | 13, 1.07e-14 | 5.00e-10 | 1.33e-30 | 0 | 3.76e-17 | — (0) |
| mixed | 5 | 6 | 2.24e-16 | 0 | 11, 4.84e-15 | 1.82e-11 | 1.06e-31 | 0 | 2.07e-17 | — (0) |
| mixed | 6 | 6 | 7.63e-16 | 0 | 14, 2.61e-14 | 7.56e-14 | 2.00e-30 | 0 | 4.79e-17 | — (0) |
| mixed | 7 | 6 | 4.61e-16 | 0 | 14, 6.00e-15 | 2.06e-12 | 1.23e-32 | 0 | 2.19e-16 | — (0) |

## C. The rich limits on one design (8 × 20; the target on atoms [1, 5, 9])

Sign-constrained basis pursuit weighted by 1/|a(0)| (the frozen signs: the target's on its atoms, random elsewhere): support [0, 3, 4, 5, 6, 8, 9, 11] (the target: no); basis pursuit: support [2, 4, 5, 7, 9, 12, 14, 19].

The vertices' certificates (|S| = N, ξ from Φ_Sᵀξ = ρ_S·sign): the smallest relative slacks off the support μ_j = 1 − s_j(Φᵀξ)_j/ρ_j — sign-constrained, weights 1/(γ|a(0)|): 0.1039, 0.1301, 0.1873; basis pursuit (equal weights): 0.0157, 0.0197, 0.0465.

The projection's distance to its limit against the leading order computed from the certificate (the support's mirror coordinates at the vertex, the off-support coefficients through each particle's mirror map: ≈ e^{O(1)}ε^{μ_j} on branches, a_maxε^{3/2}cos(πμ_j/2)/(m(2c)^{1/2}sin³(πμ_j/2)) from α(0) = 0 — the closed form's largest relative error off the support in the last column):

| start | ε | all on branches | projection (it, residual) | distance to its limit | the leading order's distance | |w − w₁|∞ | the closed form |
|---|---|---|---|---|---|---|---|
| branches (amplitudes 0.2–0.6) | 0.0001 | yes | 14, 8.08e-15 | 0.232 | 493 | 493 | — |
| branches (amplitudes 0.2–0.6) | 1e-08 | yes | 15, 1.75e-13 | 0.0298 | 6.12 | 6.09 | — |
| branches (amplitudes 0.2–0.6) | 1e-16 | yes | 18, 4.78e-14 | 0.00133 | 0.00174 | 4.14e-4 | — |
| branches (amplitudes 0.2–0.6) | 1e-32 | yes | 55, 4.98e-14 | 5.07e-6 | 5.07e-6 | 2.41e-9 | — |
| branches (amplitudes 0.2–0.6) | 1e-64 | yes | 37, 1.41e-13 | 1.05e-9 | 1.05e-9 | 3.43e-14 | — |
| branches (amplitudes 0.2–0.6) | 1e-128 | yes | 38, 2.12e-12 | 6.17e-13 | 1.22e-15 | 6.17e-13 | — |
| α(0) = 0 | 0.01 | no | 9, 5.16e-15 | 0.291 | 3.46e+05 | 3.46e+05 | 20.6 |
| α(0) = 0 | 0.0001 | no | 15, 1.84e-12 | 0.138 | 2.37 | 2.23 | 432 |
| α(0) = 0 | 1e-06 | no | 22, 1.79e-13 | 4.03e-4 | 4.07e-4 | 3.48e-6 | 7.49 |
| α(0) = 0 | 1e-08 | no | 28, 3.02e-12 | 5.19e-8 | 5.19e-8 | 8.36e-13 | 0.186 |
| α(0) = 0 | 1e-10 | no | 34, 6.08e-12 | 4.26e-11 | 4.37e-11 | 3.51e-12 | 0.0169 |

A mixed start (the even atoms at α = 0, the odd ones on branches): the free particles' total |w| and the distance to the branch atoms' own sign-constrained weighted basis pursuit:

| ε | projection (it, residual) | Σ|w| on free particles | distance |
|---|---|---|---|
| 0.0001 | 10, 1.67e-12 | 2.51e-6 | 1.96e-6 |
| 1e-08 | 15, 1.81e-12 | 7.84e-14 | 9.71e-13 |
| 1e-12 | 12, 3.43e-14 | 7.82e-22 | 2.23e-14 |

## D. A pair per atom (± 0.4, equal masses; the design of C)

The atom's coefficient W = w₊ + w₋ of the projection against basis pursuit weighted by 1/a₀; the atoms where both members carry more than 10⁻³ of the largest |W| (basis pursuit's support [2, 4, 5, 7, 9, 12, 14, 19]; its certificate's smallest relative slacks 0.0157, 0.0197, 0.0465); the leading order from the certificate as in C:

| u = w(0) | projection (it, residual) | distance | the leading order's distance | |W − W₁|∞ | atoms with both members active |
|---|---|---|---|---|---|
| 0.0001 | 8, 3.02e-12 | 0.291 | 3.61e+05 | 3.61e+05 | 0 |
| 1e-08 | 12, 8.34e-14 | 0.213 | 53 | 52.8 | 0 |
| 1e-16 | 15, 1.60e-14 | 0.123 | 2.05 | 1.93 | 0 |
| 1e-32 | 21, 3.60e-14 | 0.0203 | 0.436 | 0.416 | 0 |
| 1e-64 | 50, 1.26e-13 | 0.0111 | 0.0211 | 0.00995 | 0 |
| 1e-128 | 68, 1.74e-13 | 2.59e-4 | 2.60e-4 | 1.41e-6 | 0 |
| 1e-256 | 62, 7.47e-13 | 1.87e-6 | 1.87e-6 | 5.69e-11 | 0 |

## E. Positive kernel dictionaries (Gaussian atoms of one width on a grid of [−1, 1], samples evenly on [−1.1, 1.1], the target 1.0, 0.7, 0.5 on two adjacent atoms and one more)

### 21 atoms, width 0.15, 7 samples (the target on [5, 6, 15])

The nonnegative certificate's margin μ = 0.0146; basis pursuit's max |w − w*|: 1 with equal weights, 1.18e-12 with the columns' norms.

| start | α(0) | β(0) | on branches | max |w − w*| at t = 10³, 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ (the certificate's bound) | |Φw − y| | support at 10⁹ |
|---|---|---|---|---|---|---|---|
| branch | 1 | 1e-06 | yes | 0.854, 0.285, 0.00579, 5.35e-5 | 5.07e-5 (0.0922) | 6.82e-7 | [5, 6, 15] |
| branch | 1 | 0.01 | yes | 0.703, 0.222, 0.00584, 5.35e-5 | 5.07e-5 (0.0571) | 6.83e-7 | [5, 6, 15] |
| branch | 1 | 1 | yes | 0.589, 0.2, 0.00629, 5.95e-5 | 5.64e-5 (0.0325) | 7.59e-7 | [5, 6, 15] |
| branch | 2 | 10 | yes | 0.512, 0.185, 0.00489, 4.32e-5 | 4.09e-5 (0.0545) | 5.51e-7 | [5, 6, 15] |
| alpha=0 | 0 | 1e-06 | no | 1, 1, 1, 1 | 0.745 | 4.71e-16 | [0, 3, 6, 10, 14, 17, 20] |
| alpha=0 | 0 | 100 | no | 0.632, 0.632, 0.632, 0.632 | 1.59 | 6.29e-17 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20] |

### 41 atoms, width 0.15, 7 samples (the target on [10, 11, 30])

The nonnegative certificate's margin μ = 0.003; basis pursuit's max |w − w*|: 1.28 with equal weights, 1.41e-13 with the columns' norms.

| start | α(0) | β(0) | on branches | max |w − w*| at t = 10³, 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ (the certificate's bound) | |Φw − y| | support at 10⁹ |
|---|---|---|---|---|---|---|---|
| branch | 1 | 1e-06 | yes | 0.911, 0.647, 0.078, 0.0012 | 0.00118 (0.448) | 3.28e-6 | [10, 11, 29, 30, 31] |
| branch | 1 | 0.01 | yes | 0.834, 0.539, 0.0779, 0.0012 | 0.00118 (0.278) | 3.28e-6 | [10, 11, 29, 30, 31] |
| branch | 1 | 1 | yes | 0.774, 0.487, 0.0778, 0.00132 | 0.0013 (0.181) | 3.61e-6 | [10, 11, 29, 30, 31] |
| branch | 2 | 10 | yes | 0.732, 0.455, 0.0735, 9.82e-4 | 9.63e-4 (0.385) | 2.67e-6 | [10, 11, 29, 30, 31] |
| alpha=0 | 0 | 1e-06 | no | 1, 1.33, 1.28, 1.28 | 2.11 | 6.24e-17 | [0, 6, 12, 20, 28, 34, 40] |
| alpha=0 | 0 | 100 | no | 0.824, 0.824, 0.824, 0.824 | 2.14 | 2.22e-16 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40] |

### 41 atoms, width 0.1, 11 samples (the target on [10, 11, 30])

The nonnegative certificate's margin μ = 0.0272; basis pursuit's max |w − w*|: 2.22e-16 with equal weights, 2.22e-16 with the columns' norms.

| start | α(0) | β(0) | on branches | max |w − w*| at t = 10³, 10⁵, 10⁷, 10⁹ | the other atoms' mass at 10⁹ (the certificate's bound) | |Φw − y| | support at 10⁹ |
|---|---|---|---|---|---|---|---|
| branch | 1 | 1e-06 | yes | 0.482, 0.163, 0.00238, 2.10e-5 | 2.44e-5 (0.0952) | 5.94e-7 | [10, 11, 30] |
| branch | 1 | 0.01 | yes | 0.496, 0.162, 0.00238, 2.10e-5 | 2.44e-5 (0.0591) | 5.94e-7 | [10, 11, 30] |
| branch | 1 | 1 | yes | 0.476, 0.159, 0.00261, 2.33e-5 | 2.72e-5 (0.0384) | 6.60e-7 | [10, 11, 30] |
| branch | 2 | 10 | yes | 0.448, 0.155, 0.00195, 1.69e-5 | 1.97e-5 (0.0817) | 4.79e-7 | [10, 11, 30] |
| alpha=0 | 0 | 1e-06 | no | 1, 1, 1.69e-4, 3.04e-5 | 3.13e-5 | 1.08e-6 | [10, 11, 30] |
| alpha=0 | 0 | 100 | no | 0.609, 0.609, 0.609, 0.609 | 2.09 | 5.66e-18 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40] |

## F. The real loop on positive branches (k_e = 5, G = γ = 1, twelve atoms of width 0.25, the target on atoms [2, 6, 9], α(0) = 0.5, β(0) = 10⁻³, 4000 time units)

V̇ against −k_ee² at 12 random states (w* vanishing on nine atoms, their terms β/(γm)): largest relative difference 1.69e-8.

By quarters: max V, the rms error, the largest mass on the unused atoms; V's largest increase between samples 0 (V(0) = 3.64); max |w − w*| at the end 0.0797:

| | quarter 1 | quarter 2 | quarter 3 | quarter 4 |
|---|---|---|---|---|
| max V | 3.64 | 2.67 | 1.72 | 0.871 |
| rms e | 0.0139 | 0.0138 | 0.013 | 0.00898 |
| max β, unused atoms | 0.00552 | 0.0315 | 0.204 | 0.781 |

Numbers are produced by `experiments/exp_mirror_branch.py`; nothing is edited by hand.
