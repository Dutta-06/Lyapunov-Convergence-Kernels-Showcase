# Lemma E.0: the explicit rate against the exact Floquet rate (`exp_explicit_rate.py`)

The explicit rate is a rigorous lower bound of the exact rate in every row.

## 1. Remark E.1's instance, x_m = 1.5 sin ωt (m = n = 1, G = I, k_e = 3)

| ω | window T | α(T) | Ψ_M | Ψ_D | exact Floquet rate | explicit rate, closed-form ε | explicit rate, optimised ε | exact / explicit |
|---|---|---|---|---|---|---|---|---|
| 1 | 6.28 | 1.13 | 1.78 | 3.69 | 4.660e-02 | 3.688e-05 | 4.886e-05 | 954 |
| 2 | 3.14 | 0.566 | 1.78 | 7.38 | 2.158e-02 | 7.858e-05 | 1.040e-04 | 208 |
| 4 | 1.57 | 0.283 | 1.78 | 14.8 | 6.834e-03 | 8.696e-05 | 1.152e-04 | 59.3 |
| 8 | 0.785 | 0.142 | 1.78 | 29.5 | 1.831e-03 | 4.736e-05 | 6.301e-05 | 29.1 |
| 16 | 0.393 | 0.0708 | 1.78 | 59 | 4.662e-04 | 1.650e-05 | 2.199e-05 | 21.2 |
| 32 | 0.196 | 0.0354 | 1.78 | 118 | 1.171e-04 | 4.720e-06 | 6.286e-06 | 18.6 |
| 64 | 0.0982 | 0.0177 | 1.78 | 236 | 2.931e-05 | 1.249e-06 | 1.665e-06 | 17.6 |

Log-log slopes over ω ≥ 4: exact -1.97, explicit -1.56 (Remark E.1 reports the exact ω⁻² law).

## 2. Slow-adaptation limit: the Theorem E ideal (m = 3, n = 2, k_e = 6) with G = g·diag(2, 2, 1)

| g | Ψ_M | Ψ_D | α | Ψ_M²/k_e | κ = Ψ_M²/λ_min⟨Ψ_GᵀΨ_G⟩ | exact Floquet rate | explicit rate (optimised ε) | exact / explicit | 3κ(1 + Ψ_D/(k_eΨ_M))² | averaging value λ_min⟨Ψ_GᵀΨ_G⟩/k_e |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 9.08 | 52.3 | 0.0322 | 13.7 | 1.27e+05 | 7.465e-05 | 8.882e-16 | 8.4e+10 | 1.47e+06 | 1.079e-04 |
| 10 | 5.24 | 30.2 | 0.0107 | 4.58 | 1.27e+05 | 2.358e-05 | 2.665e-15 | 8.85e+09 | 1.47e+06 | 3.597e-05 |
| 3 | 2.87 | 16.5 | 0.00322 | 1.37 | 1.27e+05 | 7.130e-06 | 6.661e-15 | 1.07e+09 | 1.47e+06 | 1.079e-05 |
| 1 | 1.66 | 9.55 | 0.00107 | 0.458 | 1.27e+05 | 2.381e-06 | 1.865e-14 | 1.28e+08 | 1.47e+06 | 3.597e-06 |
| 0.3 | 0.908 | 5.23 | 0.000322 | 0.137 | 1.27e+05 | 7.132e-07 | 4.752e-14 | 1.5e+07 | 1.47e+06 | 1.079e-06 |
| 0.1 | 0.524 | 3.02 | 0.000107 | 0.0458 | 1.27e+05 | 2.375e-07 | 7.550e-14 | 3.15e+06 | 1.47e+06 | 3.597e-07 |
| 0.03 | 0.287 | 1.65 | 3.22e-05 | 0.0137 | 1.27e+05 | 7.124e-08 | 6.661e-14 | 1.07e+06 | 1.47e+06 | 1.079e-07 |
| 0.01 | 0.166 | 0.955 | 1.07e-05 | 0.00458 | 1.27e+05 | 2.374e-08 | 3.508e-14 | 6.77e+05 | 1.47e+06 | 3.597e-08 |

As g → 0 (slow adaptation) the ratio exact/explicit saturates near 3κ(1 + Ψ_D/(k_eΨ_M))², where κ is the condition number of the averaged Gramian (Remark E.0a): the exact rate pays for ill-conditioning once, through λ_min, and the bound pays twice.  The ratio grows further with the excitation strength Ψ_M²/k_e (fast adaptation).

## 3. The PE window (same ideal, G = diag(3, 3, 1.5))

| T / period | α(T) | explicit rate | ε |
|---|---|---|---|
| 0.25 | 0.000115 | 1.377e-14 | 2.124e-10 |
| 0.5 | 0.000658 | 3.242e-14 | 9.388e-11 |
| 1 | 0.00161 | 1.288e-14 | 1.359e-11 |
| 2 | 0.00301 | 3.109e-15 | 1.385e-12 |

Numbers are produced by `experiments/exp_explicit_rate.py`; nothing is edited by hand.
