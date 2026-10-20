# Corollary D‴.26: perfect tracking with an adaptive prior in the real loop — the Floquet spectrum (`exp_adaptive_floquet.py`)

The three-bump instance at σ₀ = 0.15σ_crit = 0.00707, x_m = 2 sin t, the nearest-atom leak toward an adaptive prior in the real closed loop; the sampled weights w = (6, 5, 4)/15 (the amplitudes A_k/w_k = 2.500, -2.100, 3.000), the lift of their exact representation with the prior at its atoms and the plant on the reference. Floquet exponents: rates −log|λ|/T₀ of the monodromy matrix (central differences of the period map, step 10⁻⁵ relative, DOP853 at rtol 10⁻¹³); the error channel's multiplier (≈ e^{−k_eT₀}) is below the resolution of the differences and is left out.

The whitened weighted Gauss–Newton spectrum at ϑ^e(w): h_j(w) = 3.0798e-04, 2.0149e-03, 3.0884e-03, 1.2820e-02, 2.0512e-02, 3.5144e-02, 9.4012e-02, 1.3952e-01, 1.8336e-01.

## The weighted K₀-atom loop against the averaged prediction

| η | k_e | return residual | collective exponents × k_e (actual / predicted (1/k_e)λ_j^±(w) × k_e) | largest relative deviation |
|---|---|---|---|---|
| 0.00707 | 25 | 1.7e-15 | 1.3166e-04/1.5231e-04, 8.7380e-04/9.3607e-04, 1.3006e-03/1.3776e-03, 3.8882e-03/3.9378e-03, 4.8329e-03/4.8712e-03, 5.6903e-03/5.7038e-03, 6.5225e-03/6.5450e-03, 6.7116e-03/6.7166e-03, 6.7950e-03/6.8018e-03, 1.4283e-02/1.4304e-02, 1.5146e-02/1.5228e-02, 1.5742e-02/1.5860e-02, 2.2782e-02/2.3030e-02, 2.9403e-02/2.9790e-02, 4.3236e-02/4.3589e-02, 9.7778e-02/1.0162e-01, 1.4508e-01/1.4696e-01, 1.8627e-01/1.9071e-01 | 1.36e-01 |
| 0.00707 | 50 | 5.8e-16 | 1.4647e-04/1.5231e-04, 9.1948e-04/9.3607e-04, 1.3573e-03/1.3776e-03, 3.9252e-03/3.9378e-03, 4.8615e-03/4.8712e-03, 5.7004e-03/5.7038e-03, 6.5392e-03/6.5450e-03, 6.7153e-03/6.7166e-03, 6.8001e-03/6.8018e-03, 1.4298e-02/1.4304e-02, 1.5206e-02/1.5228e-02, 1.5828e-02/1.5860e-02, 2.2967e-02/2.3030e-02, 2.9691e-02/2.9790e-02, 4.3499e-02/4.3589e-02, 1.0059e-01/1.0162e-01, 1.4647e-01/1.4696e-01, 1.8954e-01/1.9071e-01 | 3.84e-02 |
| ∞ | 25 | 1.4e-16 | 2.6581e-04/3.0798e-04, 1.8708e-03/2.0149e-03, 2.8942e-03/3.0884e-03, 1.2521e-02/1.2820e-02, 2.0087e-02/2.0512e-02, 3.4778e-02/3.5144e-02, 9.0151e-02/9.4012e-02, 1.3764e-01/1.3952e-01, 1.7891e-01/1.8336e-01 | 1.37e-01 |
| ∞ | 50 | 7.2e-17 | 2.9603e-04/3.0798e-04, 1.9763e-03/2.0149e-03, 3.0370e-03/3.0884e-03, 1.2743e-02/1.2820e-02, 2.0403e-02/2.0512e-02, 3.5050e-02/3.5144e-02, 9.2983e-02/9.4012e-02, 1.3904e-01/1.3952e-01, 1.8219e-01/1.8336e-01 | 3.88e-02 |

η = 0.00707: the largest deviation of the slow-time rates (exponent × k_e) from λ_j^±(w) is 4.44e-03 at k_e = 25 and 1.17e-03 at k_e = 50 (ratio 3.80, i.e. k_e^{−1.93}; relative 1.36e-01 and 3.84e-02) — the averaging correction is of second order in 1/k_e in slow time, so the exponents are (1/k_e)λ_j^±(w) + O(k_e⁻³) (Corollary D‴.26(c)).
η = ∞: the largest deviation of the slow-time rates (exponent × k_e) from λ_j^±(w) is 4.45e-03 at k_e = 25 and 1.17e-03 at k_e = 50 (ratio 3.80, i.e. k_e^{−1.93}; relative 1.37e-01 and 3.88e-02) — the averaging correction is of second order in 1/k_e in slow time, so the exponents are (1/k_e)λ_j^±(w) + O(k_e⁻³) (Corollary D‴.26(c)).

## The lifts: separations at σ₀/k_e, the collective exponents the K₀-atom loop's

| m | η | k_e | return residual | separations: count, largest relative distance from σ₀/k_e | collective exponents: count, largest relative difference from the K₀-atom loop's |
|---|---|---|---|---|---|
| 15 | 0.00707 | 25 | 1.3e-15 | 36, 1.7e-07 | 18, 5.5e-05 |
| 30 | 0.00707 | 25 | 1.4e-15 | 81, 2.5e-07 | 18, 3.9e-05 |
| 15 | 0.00707 | 50 | 1.2e-15 | 36, 1.0e-06 | 18, 2.7e-04 |
| 30 | 0.00707 | 50 | 1.7e-15 | 81, 1.7e-06 | 18, 6.4e-04 |
| 15 | ∞ | 25 | 1.0e-16 | 36, 1.6e-07 | 9, 1.4e-05 |
| 30 | ∞ | 25 | 4.4e-16 | 81, 2.7e-07 | 9, 4.0e-05 |
| 15 | ∞ | 50 | 4.3e-16 | 36, 8.8e-07 | 9, 2.4e-05 |
| 30 | ∞ | 50 | 5.2e-17 | 81, 1.6e-06 | 9, 9.4e-05 |

## C. The variational equation: the remainder of (c), to k_e = 400

The variational equation along the lift on the collapsed subspace, in the whitened weighted coordinates (δe, y, ζ): δė = −k_eδe − Ĵ(t)y, ẏ = Ĵ(t)ᵀδe − (σ₀/k_e)(y − ζ), ζ̇ = −(η/k_e)(ζ − y) (η = ∞: ζ = y); its monodromy by integrating it (DOP853, rtol 10⁻¹³; no differences), the error channel's exponent by Liouville's formula (the exponents sum to k_e + (σ₀ + η)·9/k_e; η = ∞: k_e). Two references: x_m = 2 sin t, where the first-order averaging correction Ā₁ vanishes identically (in one dimension ⟨ĴᵀĴ̇⟩ is the integral of an exact differential, and the commutator term vanishes by the symmetry t → π − t), and the asymmetric x_m = 2 sin t + 0.8 sin 2t, where it does not.

- 2 sin t, η = 0.00707: ‖Ā₁‖_max = 4.62e-17, its self-adjoint part 5.5e-17 (Ā₁ is skew); the generator against the central differences of the real loop's field at six times: 5.8e-11 (relative); h_j = 3.0798e-04, 2.0149e-03, 3.0884e-03, 1.2820e-02, 2.0512e-02, 3.5144e-02, 9.4012e-02, 1.3952e-01, 1.8336e-01.
- 2 sin t, η = ∞: ‖Ā₁‖_max = 4.62e-17, its self-adjoint part 5.7e-17 (Ā₁ is skew); the generator against the central differences of the real loop's field at six times: 5.8e-11 (relative); h_j = 3.0798e-04, 2.0149e-03, 3.0884e-03, 1.2820e-02, 2.0512e-02, 3.5144e-02, 9.4012e-02, 1.3952e-01, 1.8336e-01.
- 2 sin t + 0.8 sin 2t, η = 0.00707: ‖Ā₁‖_max = 5.63e-03, its self-adjoint part 1.4e-16 (Ā₁ is skew); the generator against the central differences of the real loop's field at six times: 1.6e-11 (relative); h_j = 4.3720e-04, 2.6482e-03, 4.0371e-03, 1.0272e-02, 2.0935e-02, 3.7038e-02, 8.4042e-02, 1.1098e-01, 1.8694e-01.
- 2 sin t + 0.8 sin 2t, η = ∞: ‖Ā₁‖_max = 5.63e-03, its self-adjoint part 1.3e-16 (Ā₁ is skew); the generator against the central differences of the real loop's field at six times: 1.6e-11 (relative); h_j = 4.3720e-04, 2.6482e-03, 4.0371e-03, 1.0272e-02, 2.0935e-02, 3.7038e-02, 8.4042e-02, 1.1098e-01, 1.8694e-01.

Against the finite differences of the first table (x_m = 2 sin t, k_e = 25 and 50, both η) the variational exponents differ by at most 4.1e-05 (relative): the resolution of the differences.

| reference | η | k_e | slowest mode: (actual − first order)/first order | largest over the modes | largest \|actual − second order\|/first order | error channel (Liouville) | k_e − tr H/k_e |
|---|---|---|---|---|---|---|---|
| 2 sin t | 0.00707 | 12.5 | -3.7534e-01 | 3.753e-01 | 2.95e-01 | 12.463921 | 12.460737 |
| 2 sin t | 0.00707 | 25 | -1.3562e-01 | 1.356e-01 | 2.67e-02 | 24.980835 | 24.980369 |
| 2 sin t | 0.00707 | 50 | -3.8400e-02 | 3.840e-02 | 1.91e-03 | 49.990246 | 49.990184 |
| 2 sin t | 0.00707 | 100 | -9.9368e-03 | 9.937e-03 | 1.24e-04 | 99.995100 | 99.995092 |
| 2 sin t | 0.00707 | 200 | -2.5063e-03 | 2.506e-03 | 7.82e-06 | 199.997547 | 199.997546 |
| 2 sin t | 0.00707 | 400 | -6.2799e-04 | 6.280e-04 | 4.90e-07 | 399.998773 | 399.998773 |
| 2 sin t | ∞ | 12.5 | -3.7793e-01 | 3.779e-01 | 2.95e-01 | 12.463923 | 12.460737 |
| 2 sin t | ∞ | 25 | -1.3692e-01 | 1.369e-01 | 2.69e-02 | 24.980835 | 24.980369 |
| 2 sin t | ∞ | 50 | -3.8810e-02 | 3.881e-02 | 1.93e-03 | 49.990246 | 49.990184 |
| 2 sin t | ∞ | 100 | -1.0046e-02 | 1.005e-02 | 1.25e-04 | 99.995100 | 99.995092 |
| 2 sin t | ∞ | 200 | -2.5341e-03 | 2.534e-03 | 7.90e-06 | 199.997547 | 199.997546 |
| 2 sin t | ∞ | 400 | -6.3497e-04 | 6.350e-04 | 4.95e-07 | 399.998773 | 399.998773 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 12.5 | -2.9312e-01 | 2.931e-01 | 4.28e-01 | 12.467004 | 12.463413 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 25 | -1.2305e-01 | 1.231e-01 | 5.17e-02 | 24.982298 | 24.981707 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 50 | -3.8958e-02 | 3.896e-02 | 4.42e-03 | 49.990936 | 49.990853 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 100 | -1.0520e-02 | 1.052e-02 | 3.07e-04 | 99.995437 | 99.995427 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 200 | -2.6859e-03 | 2.686e-03 | 1.98e-05 | 199.997715 | 199.997713 |
| 2 sin t + 0.8 sin 2t | 0.00707 | 400 | -6.7510e-04 | 6.751e-04 | 1.24e-06 | 399.998857 | 399.998857 |
| 2 sin t + 0.8 sin 2t | ∞ | 12.5 | -2.9639e-01 | 2.964e-01 | 4.28e-01 | 12.467005 | 12.463413 |
| 2 sin t + 0.8 sin 2t | ∞ | 25 | -1.2477e-01 | 1.248e-01 | 5.23e-02 | 24.982298 | 24.981707 |
| 2 sin t + 0.8 sin 2t | ∞ | 50 | -3.9554e-02 | 3.955e-02 | 4.49e-03 | 49.990936 | 49.990853 |
| 2 sin t + 0.8 sin 2t | ∞ | 100 | -1.0686e-02 | 1.069e-02 | 3.12e-04 | 99.995437 | 99.995427 |
| 2 sin t + 0.8 sin 2t | ∞ | 200 | -2.7286e-03 | 2.729e-03 | 2.01e-05 | 199.997715 | 199.997713 |
| 2 sin t + 0.8 sin 2t | ∞ | 400 | -6.8584e-04 | 6.858e-04 | 1.26e-06 | 399.998857 | 399.998857 |

- 2 sin t, η = 0.00707: per doubling of k_e the largest first-order deviation falls by 2.77, 3.53, 3.86, 3.96, 3.99 (→ 4: O(k_e⁻²) relative) and the second-order residual by 11.05, 14.00, 15.41, 15.84, 15.97 (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order prediction at every k_e: True (on this instance: every c_j < 0 here; see the gain-scaled case below). Ā₀ is self-adjoint in the metric diag(I, (σ₀/η)I) to 0.0e+00 (in I: 0.0e+00; |Ā₀| = 1.5e-01).
- 2 sin t, η = ∞: per doubling of k_e the largest first-order deviation falls by 2.76, 3.53, 3.86, 3.96, 3.99 (→ 4: O(k_e⁻²) relative) and the second-order residual by 10.96, 13.96, 15.40, 15.84, 15.96 (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order prediction at every k_e: True (on this instance: every c_j < 0 here; see the gain-scaled case below). Ā₀ is self-adjoint in the metric diag(I, (σ₀/η)I) to 0.0e+00 (in I: 0.0e+00; |Ā₀| = 1.4e-01).
- 2 sin t + 0.8 sin 2t, η = 0.00707: per doubling of k_e the largest first-order deviation falls by 2.38, 3.16, 3.70, 3.92, 3.98 (→ 4: O(k_e⁻²) relative) and the second-order residual by 8.27, 11.68, 14.41, 15.54, 15.87 (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order prediction at every k_e: True (on this instance: every c_j < 0 here; see the gain-scaled case below). Ā₀ is self-adjoint in the metric diag(I, (σ₀/η)I) to 0.0e+00 (in I: 0.0e+00; |Ā₀| = 1.7e-01).
- 2 sin t + 0.8 sin 2t, η = ∞: per doubling of k_e the largest first-order deviation falls by 2.38, 3.15, 3.70, 3.92, 3.98 (→ 4: O(k_e⁻²) relative) and the second-order residual by 8.18, 11.64, 14.39, 15.54, 15.90 (→ 16: O(k_e⁻⁴), the odd orders vanishing); every exponent below its first-order prediction at every k_e: True (on this instance: every c_j < 0 here; see the gain-scaled case below). Ā₀ is self-adjoint in the metric diag(I, (σ₀/η)I) to 0.0e+00 (in I: 0.0e+00; |Ā₀| = 1.6e-01).

The second-order coefficients (the slow-time rates are λ_j + c_j/k_e² + O(k_e⁻⁴)), c_j/λ_j in the order of λ_j; the error filter's part of the slowest mode's, −κ₁ω₁² with ω₁² = v₁ᵀ⟨Ĵ̇ᵀĴ̇⟩v₁/h₁ = Σ_n n²|Ĵ_nv₁|²/Σ_n|Ĵ_nv₁|² — the mode's mean squared frequency in the regressor, independent of η — and κ₁ = h₁(dλ₁/dh)/λ₁ (1 at η = ∞); and the measured (actual·k_e − λ₁)k_e²/λ₁ of the slowest mode at k_e = 100, 200, 400:

- 2 sin t, η = 0.00707: c_j/λ_j = -100.55, -45.38, -37.53, -8.05, -4.98, -1.50, -2.25, -0.47, -0.66, -1.12, -3.70, -5.03, -6.99, -8.38, -5.21, -25.76, -8.31, -15.58; the slowest c₁/λ₁ = -100.55, of which the filter's -102.14 = −0.989ω₁², ω₁ = 10.1625; measured -99.37, -100.25, -100.48.
- 2 sin t, η = ∞: c_j/λ_j = -101.67, -49.08, -42.55, -15.04, -13.37, -6.71, -28.00, -8.78, -16.23; the slowest c₁/λ₁ = -101.67, of which the filter's -103.28 = −1.000ω₁², ω₁ = 10.1625; measured -100.46, -101.37, -101.59.
- 2 sin t + 0.8 sin 2t, η = 0.00707: c_j/λ_j = -108.21, -53.12, -47.97, -12.35, -6.18, -1.59, -3.08, -1.32, -0.94, -1.72, -5.93, -8.89, -7.86, -10.70, -5.90, -31.07, -18.22, -22.78; the slowest c₁/λ₁ = -108.21, of which the filter's -109.72 = −0.984ω₁², ω₁ = 10.5578; measured -105.20, -107.44, -108.02.
- 2 sin t + 0.8 sin 2t, η = ∞: c_j/λ_j = -109.93, -59.05, -56.86, -20.21, -16.87, -7.49, -34.15, -19.54, -23.72; the slowest c₁/λ₁ = -109.93, of which the filter's -111.47 = −1.000ω₁², ω₁ = 10.5578; measured -106.86, -109.14, -109.73.

The adaptation gain scaled by 30 (2 sin t, η = ∞): c_j/λ_j = -62.59, -23.76, -12.10, -12.05, -16.14, -9.78, -13.06, +4.17, +12.65 — 2 positive (the slaving's −(J̄ᵀJ̄)² and the level repulsion outweigh the filter where |Ĵ|² is large against the frequencies). The relative deviations (actual − first order)/first order of those modes at k_e = 100, 200, 400: +4.171e-04, +1.270e-03; +1.043e-04, +3.166e-04; +2.608e-05, +7.911e-05, against c_j/(λ_jk_e²) = +4.173e-04, +1.265e-03; +1.043e-04, +3.164e-04; +2.608e-05, +7.909e-05: those exponents lie above their first-order prediction, and the slowest below it (-6.207e-03, -1.562e-03, -3.910e-04).

At η = 10σ₀ (2 sin t), where the metric W = diag(I, (σ₀/η)I) is not the identity: |WĀ₀ − (WĀ₀)ᵀ| = 8.7e-18 against |Ā₀ − Ā₀ᵀ| = 6.4e-02 (|Ā₀| = 1.5e-01): Ā₀ is self-adjoint in W only; |WĀ₁ + (WĀ₁)ᵀ| = 5.5e-17 (|Ā₁| = 4.6e-17).

At η = 10σ₀ (2 sin t + 0.8 sin 2t), where the metric W = diag(I, (σ₀/η)I) is not the identity: |WĀ₀ − (WĀ₀)ᵀ| = 8.7e-18 against |Ā₀ − Ā₀ᵀ| = 6.4e-02 (|Ā₀| = 1.7e-01): Ā₀ is self-adjoint in W only; |WĀ₁ + (WĀ₁)ᵀ| = 1.4e-16 (|Ā₁| = 5.6e-03).

## D. Lemma D‴.26a: the exact slow subspace (x_m = 2 sin t + 0.8 sin 2t)

Chang's P(t; ε) by its contraction on a 512-point Fourier grid: the iterations to 10⁻¹⁶ and the invariance residual; the truncation |P − P₂|_∞/ε³ (P₂ = −K + εK̇ − ε²(K̈ + KKᵀK + KΛ)); the slow block's slow-time rates (from log Φ(ε)/(εT₀)) against the full variational equation's; the rates at ±ε (the slow block at −ε defined by the same equation) — max|ν(ε) − ν(−ε)|/max|ν| — against their first-order deviation max|ν(ε) − ν₀|/max|ν₀|; the same at ±ε with the plant's slope as a mismatch in the error channel; the Magnus parameter εT₀max‖A‖₂ (the series converges below π).

η = 0.00707:

| k_e | iterations | invariance residual | |P − P₂|/ε³ | slow block vs full | rates at ±ε | first-order deviation | ±ε with a mismatch | Magnus |
|---|---|---|---|---|---|---|---|---|
| 12.5 | 9 | 3.2e-15 | 1849.7804 | 1.5e-11 | 4.3e-14 | 9.45e-02 | 1.08e-02 | 0.480 |
| 25 | 7 | 1.9e-15 | 2621.0622 | 3.3e-11 | 7.6e-14 | 3.10e-02 | 5.18e-03 | 0.247 |
| 50 | 6 | 1.0e-15 | 3074.0166 | 1.7e-11 | 1.4e-13 | 8.68e-03 | 2.52e-03 | 0.125 |
| 100 | 5 | 6.7e-16 | 3268.2155 | 4.7e-11 | 3.2e-13 | 2.25e-03 | 1.25e-03 | 0.062 |
| 200 | 4 | 4.8e-16 | 3305.7286 | 1.2e-10 | 9.0e-13 | 5.68e-04 | 6.24e-04 | 0.031 |
| 400 | 4 | 3.3e-16 | 3332.2482 | 8.2e-10 | 1.6e-12 | 1.42e-04 | 3.12e-04 | 0.016 |

Under white measurement noise (N₀ = 1): the second transformation's decoupling residual, the slow noise b against b₀ = −(K̇ + aK)ᵀ (max over the period), and the exact periodic covariance at t = 0 in the coordinates (v, Y′) — Cov(v)·2ε³ − 1, ‖E[vY′ᵀ] + bᵀ‖/‖b‖, ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ — and the function-level floor at t = 0, tr(H Cov(y)), against tr(H(Σ̄ + (ε/2)KᵀK)) (the fast jitter's share).

| k_e | decoupling residual | ‖b − b₀‖/‖b₀‖ | Cov(v)·2ε³ − 1 | ‖E[vY′ᵀ] + bᵀ‖/‖b‖ | ‖Cov(Y′) − Σ̄‖/‖Σ̄‖ | floor at t = 0: exact | predicted | ratio | fast share |
|---|---|---|---|---|---|---|---|---|---|
| 12.5 | 3.0e-16 | 6.271e-01 | -9.600e-02 | 1.543e-01 | 3.581e-01 | 3.23656e-01 | 4.28924e-01 | 0.7546 | 0.001 |
| 25 | 8.3e-17 | 4.023e-01 | -2.840e-02 | 7.719e-02 | 1.801e-01 | 1.93447e-01 | 2.14462e-01 | 0.9020 | 0.001 |
| 50 | 2.0e-17 | 2.244e-01 | -7.708e-03 | 3.863e-02 | 6.180e-02 | 1.04025e-01 | 1.07231e-01 | 0.9701 | 0.001 |
| 100 | 6.9e-18 | 1.168e-01 | -2.083e-03 | 1.931e-02 | 1.714e-02 | 5.31876e-02 | 5.36155e-02 | 0.9920 | 0.001 |
| 200 | 2.9e-18 | 5.874e-02 | -5.888e-04 | 9.650e-03 | 4.405e-03 | 2.67531e-02 | 2.68078e-02 | 0.9980 | 0.001 |
| 400 | 1.1e-18 | 2.950e-02 | -1.806e-04 | 4.825e-03 | 1.108e-03 | 1.33969e-02 | 1.34039e-02 | 0.9995 | 0.001 |

η = ∞:

| k_e | iterations | invariance residual | |P − P₂|/ε³ | slow block vs full | rates at ±ε | first-order deviation | ±ε with a mismatch | Magnus |
|---|---|---|---|---|---|---|---|---|
| 12.5 | 9 | 3.8e-15 | 1849.8753 | 5.7e-12 | 2.8e-14 | 9.84e-02 | 1.13e-02 | 0.476 |
| 25 | 8 | 2.2e-15 | 2621.2101 | 1.2e-12 | 3.9e-14 | 3.23e-02 | 5.42e-03 | 0.245 |
| 50 | 6 | 8.9e-16 | 3074.1950 | 2.3e-11 | 6.7e-14 | 9.04e-03 | 2.63e-03 | 0.124 |
| 100 | 5 | 5.6e-16 | 3268.4044 | 5.7e-11 | 3.0e-13 | 2.34e-03 | 1.30e-03 | 0.062 |
| 200 | 5 | 4.4e-16 | 3305.9193 | 4.5e-10 | 4.5e-13 | 5.91e-04 | 6.50e-04 | 0.031 |
| 400 | 4 | 3.3e-16 | 3332.4406 | 1.5e-09 | 9.3e-13 | 1.48e-04 | 3.25e-04 | 0.016 |


Numbers are produced by `experiments/exp_adaptive_floquet.py`; nothing is edited by hand.
