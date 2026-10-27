# Proposition D‴.46: a leak in the mirror coordinates (`exp_mirror_leak.py`)

The free law with the leak −σ(η_i − η_p,i) on each particle's scalar drive (Q kept); η, B the mirror coordinate and Bregman divergence of Proposition D‴.45 (`mirror_free`).

## A. Identities

V̇ against −k_ee² + ed − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop (k_e = 5, a disturbance, σ = 0 and 0.3, 12 random states): largest relative difference 1.18e-10.

The averaged equilibrium against its limits on random positive definite H (10 × 10, 4 draws): the ridge (1/k_e)(w − w*)ᵀH(w − w*) + σ(m/(Gβ⁰²))|w|² at large starting masses (σ = 10⁻²), the LASSO with penalty σλ₁(β⁰) = 0.02 at small ones (relative distance, largest over the draws; the supports equal):

| the limit | β⁰ | distance | supports equal |
|---|---|---|---|
| ridge | 10000 | 2.26e-13 | — |
| ridge | 100 | 2.24e-8 | — |
| lasso | 0.0001 | 0.00721 | 4 of 4 |
| lasso | 1e-06 | 7.13e-4 | 4 of 4 |
| lasso | 1e-08 | 7.13e-5 | 4 of 4 |

## B. The real loop under an approximation error (0.1 sin 4x) and a disturbance (0.2 sin 3.1t + 0.1), 3000 time units (k_e = 5, m = 6, G = 1, γ = 0.5; six atoms of width 0.4)

By quarters: the largest |w|, the largest V (V = ½e² + ΣB(w*, w), w* the representable part), the rms error; the largest |Q − Q(0)|; and the ultimate bound max(V(0), (d̄²/(2k_e) + σB(w*, w_p))/min(k_e, σ)) with d̄ = 0.4 (the disturbance and the approximation error):

| σ | max |w| by quarters | max V by quarters | rms e by quarters | max |Q − Q(0)| | the bound |
|---|---|---|---|---|---|
| 0 | 0.198, 0.246, 0.273, 0.288 | 0.182, 0.281, 0.344, 0.37 | 0.0299, 0.0250, 0.0245, 0.0244 | 1.29e-11 | — |
| 0.01 | 0.081, 0.081, 0.080, 0.081 | 0.182, 0.0847, 0.0847, 0.0847 | 0.0332, 0.0309, 0.0309, 0.0309 | 5.03e-12 | 1.74 |
| 0.1 | 0.040, 0.040, 0.040, 0.040 | 0.182, 0.13, 0.13, 0.13 | 0.0382, 0.0366, 0.0366, 0.0366 | 1.63e-11 | 0.297 |

## C. Slow adaptation (k_e = 20, G = 0.5, γ = 0.05, 40000 time units, from α = 0 and equal masses β⁰, the leak toward the start): the late mean of w (the last quarter) against the averaged loop's equilibrium w_σ

| start | β⁰ | σ | max |w̄ − w_σ| | max |w_σ − w*| | max |w̄ − w*| | the last quarter's drift |
|---|---|---|---|---|---|---|
| lazy | 3 | 0.0005 | 6.94e-4 | 0.0374 | 0.038 | 1.32e-4 |
| intermediate | 0.3 | 0.0001 | 7.06e-4 | 0.123 | 0.124 | 4.36e-4 |

## D. The regularisation paths on the orbit's Gram (12 atoms of width 0.25 on [−1.4, 1.4]; its eigenvalues 1.42e-4 to 0.553; a 3-sparse truth on atoms [2, 6, 9])

The equilibrium against the ridge (lazy start, β⁰ = 30) and the LASSO of penalty μ = σλ₁(β⁰) (rich start, β⁰ = 10⁻⁶); its support (|w| > 10⁻³ max|w|), whether it keeps the truth's atoms, the atoms beyond them, max |w − w*|:

| start | μ (rich) or σ (lazy) | distance to the reference | support | the truth kept | extra atoms | max |w − w*| |
|---|---|---|---|---|---|---|
| lazy | 0.0001 | 3.86e-6 | [2, 6, 7, 8, 9, 10] | yes | 3 | 4.84e-4 |
| lazy | 0.0003 | 1.13e-5 | [2, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 6 | 0.00143 |
| lazy | 0.001 | 3.50e-5 | [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 8 | 0.00458 |
| lazy | 0.003 | 8.70e-5 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 9 | 0.0124 |
| lazy | 0.01 | 1.80e-4 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 9 | 0.0312 |
| lazy | 0.03 | 2.49e-4 | [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11] | yes | 8 | 0.0571 |
| rich | 0.0001 | 3.71e-6 | [2, 6, 9] | yes | 0 | 7.52e-4 |
| rich | 0.0003 | 1.11e-5 | [2, 6, 9] | yes | 0 | 0.00226 |
| rich | 0.001 | 3.68e-5 | [2, 6, 9] | yes | 0 | 0.00752 |
| rich | 0.003 | 1.09e-4 | [2, 6, 9] | yes | 0 | 0.0226 |
| rich | 0.01 | 3.40e-4 | [2, 6, 9] | yes | 0 | 0.0752 |
| rich | 0.03 | 6.53e-4 | [2, 9] | no | 0 | 0.214 |

Numbers are produced by `experiments/exp_mirror_leak.py`; nothing is edited by hand.
