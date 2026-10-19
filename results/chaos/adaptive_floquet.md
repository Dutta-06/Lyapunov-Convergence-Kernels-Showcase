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

Numbers are produced by `experiments/exp_adaptive_floquet.py`; nothing is edited by hand.
